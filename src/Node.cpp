// Defining the Node class
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <map>


#include "Node.h"
#include "Grid.h"
#include "Common.h"
#include "Config.h"

Node::Node(std::string ID, Point location, int numCattle, Grid *gridPtr, const Config *config, std::mt19937 *rng)
        : m_epiID{std::move(ID)}, m_location{location}, m_config(config), m_parentGrid(gridPtr), m_pRNG(rng),
          m_numCattle(numCattle) {
    // assume always have serotypes, no serotypes caught by input checking
    for (const auto &serotype : m_config->m_serotypes) {
        m_serotypeCompartments.emplace_back(m_config, serotype, numCattle);
    }
}

std::vector<SerotypeCompartmentalModel>::iterator Node::findSerotype(std::string_view serotype) {
    auto it{std::find_if(m_serotypeCompartments.begin(), m_serotypeCompartments.end(),
                         [&](SerotypeCompartmentalModel &scm) { return scm.getSerotype() == serotype; })};
    if (it != m_serotypeCompartments.end()) {
        return it;
    } else {
        throw std::runtime_error("Node::findSerotype() somehow cannot find specified serotype");
    }
}

std::vector<SerotypeCompartmentalModel>::const_iterator Node::findSerotype(std::string_view serotype) const {
    auto it{std::find_if(m_serotypeCompartments.begin(), m_serotypeCompartments.end(),
                         [&](const SerotypeCompartmentalModel &scm) { return scm.getSerotype() == serotype; })};
    if (it != m_serotypeCompartments.end()) {
        return it;
    } else {
        throw std::runtime_error("Node::findSerotype() somehow cannot find specified serotype");
    }
}


std::string Node::getID() const { return m_epiID; }

Point Node::getLocation() const { return m_location; }

bool Node::isInfected() const { return m_isInfected; }

void Node::update() {
    // also advance daysSinceInfected, MovementBan etc.
    checkIsInfected();
    updateDSI();
    updateDaysSinceMovementBan();
    updateDaysSinceVaccinated();
    if (m_isSentinelNode) {
        if (m_parentCell->getMaxNodeSize() != m_numCattle) {
            m_parentCell->updateMaxNodeSize(m_numCattle);
        }
    }
}


std::pair<bool, std::string> Node::attemptToInfect(std::map<std::string, double> &serotypeInfectionProbabilities) {
    // trim serotypes cannot be infected with
    for (auto it = serotypeInfectionProbabilities.begin(); it != serotypeInfectionProbabilities.end();) {
        if (canBeInfectedWithSerotype(it->first)) {
            ++it;
        } else {
            it = serotypeInfectionProbabilities.erase(it);
        }
    }
    if (!serotypeInfectionProbabilities.empty()) {
        // select between serotype based on relative probabilities of infection.
        std::string selectedSerotype{"none"};
        int numInfected;
        if (serotypeInfectionProbabilities.size() == 1) {
            selectedSerotype = serotypeInfectionProbabilities.begin()->first;
        } else {
            std::vector<double> serotypeNumbers(serotypeInfectionProbabilities.size());
            std::vector<std::string> serotypeNames(serotypeInfectionProbabilities.size());
            for (const auto&[serotype, P] : serotypeInfectionProbabilities) {
                serotypeNames.push_back(serotype);
                serotypeNumbers.push_back(P);
            }
            std::discrete_distribution<size_t> discreteDistribution(serotypeNumbers.begin(), serotypeNumbers.end());
            selectedSerotype = serotypeNames.at(discreteDistribution(*(m_pRNG)));
        }
        auto it{findSerotype(selectedSerotype)};
        /* multiply by the proportion of the population that is actually susceptible, because assuming a random mixing of the population
         * the probability that the cattle is actually susceptible needs to be taken into account. */
        double trueProb{serotypeInfectionProbabilities.at(selectedSerotype) *
                        (double((*it).getSumSusceptible()) / double((*it).getTotalPopulation()))};
        std::binomial_distribution<int> binomialDistribution((*it).getSumSusceptible(), trueProb);
        numInfected = binomialDistribution(*(m_pRNG));

        if (numInfected > 0) {
            // infect available susceptible cattle
            bool successfullyInfected = (*it).infect(numInfected);
            // only update properties if successful and Node not already infected
            if (successfullyInfected && !m_isInfected && m_serotypeInfectedWith == "none") {
                m_isInfected = true;
                m_serotypeInfectedWith = selectedSerotype;
                // decide based on detection rate whether this is detectable
                std::uniform_real_distribution<double> uniformDist(0.0, 100.0);
                double detectable = uniformDist((*(m_pRNG)));
                detectable <= m_config->m_detectionRate ? (m_infectionDetectable = true)
                                                        : (m_infectionDetectable = false);
            }
            return std::make_pair(true, selectedSerotype);
        } else {
            return std::make_pair(false, "none");
        }
    } else {
        return std::make_pair(false, "none");
    }
}


void Node::setMovement(bool newMovementPermission) {
    m_allowedMovements = newMovementPermission;
    m_daysSinceMovementBan = 0;
}

std::pair<int, int> Node::calculateBirthsAndDeaths() const {
    // Separate the events that are not associated with specific serotypes so don't do them multiple times
    // e.g. number of births has nothing to do with serotypes
    std::vector<double> rates{m_config->m_dailyBirthRate * double(m_numCattle),
                              m_config->m_dailyNaturalMortality * double(m_numCattle)};
    std::pair<int, int> events{0, 0};

    // births
    if (rates[0] > 0.0) {
        std::poisson_distribution<int> birthDist(rates[0]);
        events.first = birthDist(*(m_pRNG));
    }

    // deaths
    if (rates[1] > 0.0) {
        std::poisson_distribution<int> deathsDist(rates[1]);
        events.second = deathsDist(*(m_pRNG));
    }

    return events;
}

void Node::tauLeap() {

    // Calculate births and non-serotype deaths before resolve serotype-specific stuff, then resolve births etc.
    auto birthsAndDeaths{calculateBirthsAndDeaths()};

    // can only be 1 serotype at a time, so only 1 SCM will have infection mortality
    int serotypeDeaths{0};
    for (auto &scm : m_serotypeCompartments) {
        scm.tauLeap(*(m_pRNG));
        if (scm.getSerotype() == m_serotypeInfectedWith) {
            serotypeDeaths += scm.infectionMortality(*(m_pRNG));
        }
    }

    // harmonise serotype deaths (each SCM should have the same population at the end)
    for (auto &scm : m_serotypeCompartments) {
        if (scm.getSerotype() != m_serotypeInfectedWith) {
            if (birthsAndDeaths.second + serotypeDeaths > 0) {
                scm.imposeDeaths(*(m_pRNG), birthsAndDeaths.second + serotypeDeaths);
            }
        } else {
            if (birthsAndDeaths.second > 0) {
                scm.imposeDeaths(*(m_pRNG), birthsAndDeaths.second);
            }
        }
        if (birthsAndDeaths.first > 0) {
            scm.births(*(m_pRNG), birthsAndDeaths.first);
        }
    }
    // double check all populations are in sync
    // TODO: make this a test
    int newPop{m_serotypeCompartments.front().getTotalPopulation()};
    for (const auto &scm : m_serotypeCompartments) {
        if (scm.getTotalPopulation() != newPop) {
            throw std::runtime_error(
                    "Error, multiple SerotypeCompartmentalModels disagree on Node total population. Called from Node::tauLeap().");
        }
    }
    m_numCattle = newPop;
    updateTransmission();
}


Shipment Node::generateShipment(int numShipped) {
    // TODO: Refactor shipping to be single function interaction between Nodes
    for (const auto &scm : m_serotypeCompartments) {
        if (scm.getTotalPopulation() != m_numCattle) {
            throw std::runtime_error(
                    "Error, multiple SerotypeCompartmentalModels disagree on Node total population. Called from Node::tauLeap().");
        }
    }
    // generate the actual shipment of cattle
    Shipment shipment;
    shipment.sourceIsInfected = m_isInfected;
    shipment.infectionSerotype = m_serotypeInfectedWith;
    // ASSUMPTION: all cattle in any class are equally likely to be shipped

    for (auto &scm : m_serotypeCompartments) {
        SerotypeCompartmentalModel shipmentSCM{scm.sampleForShipment(*(m_pRNG), numShipped)};
        shipment.serotypeCompartments.push_back(shipmentSCM);
        shipment.shipmentIsInfected = shipmentSCM.getIsInfected();
    }
    // if > population, all SCM should truncate to same number and so have same population
    shipment.size = shipment.serotypeCompartments.front().getTotalPopulation();
    m_numCattle = m_serotypeCompartments.front().getTotalPopulation();

    return shipment;
}

ShipmentResult Node::acceptShipment(Shipment &incomingShipment) {
    ShipmentResult result;
    result.targetInfectedPrior = m_isInfected;

    std::uniform_real_distribution<double> unifDist(0.0, 100.0);
    result.fomiteTransmissionOccured = (incomingShipment.sourceIsInfected &&
                                        unifDist(*m_pRNG) <= m_config->m_shipmentFomiteSpreadProbability);

    if (m_isInfected) {
        if ((incomingShipment.shipmentIsInfected || result.fomiteTransmissionOccured)) {
            if (incomingShipment.infectionSerotype != m_serotypeInfectedWith) {
                for (auto &scm : incomingShipment.serotypeCompartments) {
                    if (scm.getSerotype() == incomingShipment.infectionSerotype) {
                        scm.uninfect();
                    }
                }
                incomingShipment.shipmentIsInfected = false;
                incomingShipment.sourceIsInfected = false;
                incomingShipment.infectionSerotype = "none";
            } else {
                result.targetNowInfected = true;
                // otherwise do nothing as already infected with same serotype and shipment will be added
            }
        }
    } else {
        // can now proceed with assumption that infection is valid
        if (incomingShipment.shipmentIsInfected || result.fomiteTransmissionOccured) {
            m_isInfected = true;
            result.targetNowInfected = true;
            m_serotypeInfectedWith = incomingShipment.infectionSerotype;
            double detectable = unifDist((*(m_pRNG)));
            detectable <= m_config->m_detectionRate ? (m_infectionDetectable = true)
                                                    : (m_infectionDetectable = false);
            if (result.fomiteTransmissionOccured && not incomingShipment.shipmentIsInfected) {
                // need to add infected animals as none in shipment or farm
                for (auto &scm : m_serotypeCompartments) {
                    if (scm.getSerotype() == incomingShipment.infectionSerotype && scm.getSumSusceptible() > 0) {
                        // TODO: decide on better way of deciding number infected, should be low, e^-x?
                        scm.infect(1);
                    }
                }
            }
        }
    }
    for (auto &shipmentSCM : incomingShipment.serotypeCompartments) {
        auto it{findSerotype(shipmentSCM.getSerotype())};
        (*it) += shipmentSCM;
    }
    m_numCattle += incomingShipment.size;
    return result;
}

int Node::getDaysSinceInfected() const {
    return m_daysSinceInfected;
}

void Node::updateDSI() {
    if (m_isInfected) {
        ++m_daysSinceInfected;
    } else {
        m_daysSinceInfected = 0;
    }
}

void Node::vaccinate() {
    /* ASSUMPTION: Calling this function takes place AFTER logic has been run to decide whether vaccination should take place*/
    /* ASSUMPTION: Vaccinating cattle who have maternal immunity, are already vaccinated, or are infected doesn't do anything*/

    auto vaccine = m_config->m_vaccine;
    for (auto &scm : m_serotypeCompartments) {
        // std::normal_distribution has to return double/float/long double
        std::normal_distribution<double> normDist(vaccine.serotypeEfficacies.at(scm.getSerotype()).meanEfficacy,
                                                  vaccine.serotypeEfficacies.at(scm.getSerotype()).stdevEfficacy);
        double serotypeEfficacy{normDist((*(m_pRNG))) / 100.0};
        // bounds check
        serotypeEfficacy = std::clamp(serotypeEfficacy, 0.0, 1.0);

        // Calculate proportion that would be successfully vaccinated and update numbers
        scm.vaccinate(*(m_pRNG), serotypeEfficacy);
    }
    m_isQueuedForVaccination = false;
    m_hasBeenVaccinated = true;
    m_daysSinceVaccination = 0;
}

void Node::updateDaysSinceMovementBan() {
    // update days since movement ban
    // to be used in tau-leap
    if (!m_allowedMovements) {
        // if movement banned
        ++m_daysSinceMovementBan;
    }
}

bool Node::getMoveStatus() const {
    return m_allowedMovements;
}

int Node::getDaysSinceMovementBan() const {
    return m_daysSinceMovementBan;
}

void Node::checkIsInfected() {
    /* Double check the infection status of the Node, and adjust as necessary.
     * Three scenarios:
     * - Not infected
     * - Infected with 1 serotype
     * - Infected with 2+ serotypes -> prefer one already infected with, otherwise choose randomly */

    std::vector<size_t> infectedIndices;
    for (size_t i = 0; i < m_serotypeCompartments.size(); ++i) {
        if (m_serotypeCompartments[i].getIsInfected()) {
            infectedIndices.push_back(i);
        }
    }
    auto numInfected{infectedIndices.size()};
    switch (numInfected) {
        case 0:
            m_isInfected = false;
            m_serotypeInfectedWith = "none";
            m_infectionDetectable = true;
            break;
        case 1:
            m_isInfected = true;
            m_serotypeInfectedWith = m_serotypeCompartments[infectedIndices[0]].getSerotype();
            break;
        default:
            // hopefully never gets here. Should lof if it does, once I implement logging.
            if (m_isInfected) {
                for (auto &scm : m_serotypeCompartments) {
                    if (scm.getSerotype() != m_serotypeInfectedWith) {
                        scm.uninfect();
                    }
                }
            } else {
                m_isInfected = true;
                std::uniform_int_distribution<size_t> dist(0, numInfected - 1);
                size_t selected{dist(*(m_pRNG))};
                for (size_t i = 0; i < numInfected; ++i) {
                    size_t ind{infectedIndices[i]};
                    if (i == selected) {
                        m_serotypeInfectedWith = m_serotypeCompartments[ind].getSerotype();
                    } else {
                        m_serotypeCompartments[ind].uninfect();
                    }
                }
            }
            break;
    }
}

bool Node::detectInfected() {
    m_hasBeenDetected = (m_isInfected && m_infectionDetectable &&
                         (m_daysSinceInfected > int(m_config->m_detectionDelay)));
    return m_hasBeenDetected;
}

void Node::uninfect() {
    if (m_isInfected) {
        auto it{findSerotype(m_serotypeInfectedWith)};
        (*it).uninfect();
    }
}

const std::string &Node::getSerotypeInfectedWith() const {
    return m_serotypeInfectedWith;
}

bool Node::canBeInfectedWithSerotype(const std::string &infectingSerotype) const {
    // ASSUMPTION: Already infected premises can still be infected with the same serotype, so long as S > 0
    // it just converts one S -> E
    // ASSUMPTIOM: Node infected with another serotype cannot be infected
    auto it{findSerotype(infectingSerotype)};
    bool susceptibleToSerotype = (*it).getCompartmentPopulation(Status::SUS) > 0;
    bool noConflictingSerotype = (m_serotypeInfectedWith == infectingSerotype || m_serotypeInfectedWith == "none");
    return ((susceptibleToSerotype && noConflictingSerotype));
}

Cell *Node::getCellIn() const {
    return m_parentCell;
}

void Node::setCellIn(Cell *cell) {
    m_parentCell = cell;
}

void Node::setGridIn(Grid *gridPtr) {
    m_parentGrid = gridPtr;
}

int Node::getNumCattle() const {
    return m_numCattle;
}

double Node::getCattleInfectionProb(const Node *otherNode) {
    // This is intended to calculate the probability of this node infecting the other node.
    if (m_isInfected) {
        // distance
        double distance = common::euclideanDistance(m_location, otherNode->getLocation());
        double kernel = m_parentGrid->powerLawKernel(distance);
        // transmission_pars take into account number of infectious cattle
        return common::infectionProbabilty(m_transmissionPars.first, m_transmissionPars.second, kernel);
    } else {
        return 0.0;
    }
}


void Node::setIsSentinel(bool newVal) {
    m_isSentinelNode = newVal;
}

void Node::updateDaysSinceVaccinated() {
    if (m_hasBeenVaccinated) {
        ++m_daysSinceVaccination;
    }
}

int Node::getDaysSinceVaccinated() const {
    return m_daysSinceVaccination;
}

bool Node::hasBeenVaccinated() const {
    return m_hasBeenVaccinated;
}

void Node::updateTransmission() {

    if (m_isInfected) {
        auto it{findSerotype(m_serotypeInfectedWith)};
        int totalInfectious = (*it).getSumInfectious();

        m_transmissionPars.first = totalInfectious *
                m_config->m_serotypeInfectionParameters.at(
                        m_serotypeInfectedWith).transmission;
        // don't multiply susceptibility, is multiplied by num susceptible in recipient node.
        m_transmissionPars.second = m_config->m_serotypeInfectionParameters.at(
                m_serotypeInfectedWith).susceptibility;
    } else {
        m_transmissionPars.first = 0.0;
        m_transmissionPars.second = 0.0;
    }

}

bool Node::isInfectiousWith(const std::string &serotype) const {
    auto it{findSerotype(serotype)};
    return (m_isInfected && (*it).getSumInfectious() > 0);
}

std::vector<int> Node::getSumCompartments(std::string_view serotype) const {
    auto it{findSerotype(serotype)};
    std::vector<int> compartmentSums{(*it).getCompartmentSums()};
    return compartmentSums;
}

bool Node::getIsQueuedForVaccination() const {
    return m_isQueuedForVaccination;
}

void Node::setIsQueuedForVaccination(bool mIsQueuedForVaccination) {
    m_isQueuedForVaccination = mIsQueuedForVaccination;
}

bool Node::getHasBeenDetected() const {
    return m_hasBeenDetected;
}

