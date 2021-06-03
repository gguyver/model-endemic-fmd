//
// Created by glen on 13/02/2020.
//

#include "Reporting.h"

#include "Config.h"
#include "Node.h"

std::ostream &operator<<(std::ostream &os, DayContext const &context) {
    os << context.run << ", " << context.day << ", " << context.serotype;
    return os;
}

std::ostream &operator<<(std::ostream &os, const InfectionEvent &event) {
    os << event.context << ", " << event.type << ", " << event.sourceEpi
       << ", " << event.targetEpi << ", " << (event.targetInfectedPrior ? "true\n" : "false\n");
    return os;
}

std::ostream &operator<<(std::ostream &os, const InfectedNodeSummary &ins) {
    os << ins.context << ", " << ins.nodeID << ", " << ins.daysInfected << '\n';
    return os;
}

std::ostream &operator<<(std::ostream &os, DailySerotypeSummary const &dss) {
    os << dss.context << ", true_incidence, " << dss.true_incidence << '\n';
    os << dss.context << ", true_prevalence, " << dss.true_prevalence << '\n';
    os << dss.context << ", detected_incidence, " << dss.detected_incidence << '\n';
    os << dss.context << ", detected_prevalence, " << dss.detected_prevalence << '\n';
    return os;
}

std::ostream &operator<<(std::ostream &os, DailyCompartmentSums const &dcs) {
    os << dcs.context << ", MAT, " << dcs.sums.at(0) << '\n';
    os << dcs.context << ", SUS, " << dcs.sums.at(1) << '\n';
    os << dcs.context << ", EXP, " << dcs.sums.at(2) << '\n';
    os << dcs.context << ", INF, " << dcs.sums.at(3) << '\n';
    os << dcs.context << ", REC, " << dcs.sums.at(4) << '\n';
    os << dcs.context << ", VAC, " << dcs.sums.at(5) << '\n';
    os << dcs.context << ", CAR, " << dcs.sums.at(6) << '\n';
    return os;
}

Reporting::Reporting(Config *configPtr) : m_config(configPtr), m_outputDirectory(configPtr->m_outputDirectory),
                                          m_uniqueID(configPtr->m_uniqueID) {
    if (m_outputDirectory.back() != '/') {
        m_outputDirectory += "/";
    }
    initModelReport();
    initCompartmentSumsReport();
}

void Reporting::initModelReport() {
    // calculate size so can reserve
    size_t entriesPerDay{(m_config->m_serotypes.size())};
    size_t numEntries{((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)) * entriesPerDay)};
    if (m_config->m_burnIn) {
        numEntries += ((static_cast<unsigned long long int>(m_config->m_burnInDuration)) * entriesPerDay);
    }
    m_nodeInfectionsReport.reserve(numEntries);
}

void Reporting::initCompartmentSumsReport() {
    size_t itemsPerDay{m_config->m_serotypes.size()};
    size_t numItems{((m_config->m_numModelRuns) * (static_cast<unsigned int>(m_config->m_numDays)) * itemsPerDay) + 1};
    if (m_config->m_burnIn) {
        numItems += ((static_cast<unsigned long long int>(m_config->m_burnInDuration)) * itemsPerDay);
    }
    m_compartmentSums.reserve(numItems);
}

DailySerotypeSummary
Reporting::sumSerotypeStatistics(std::vector<Node> &nodes, DayContext const &context) {

    DailySerotypeSummary summary(context);
    for (auto &nodePtr : nodes) {
        if (nodePtr.isInfected() && nodePtr.getSerotypeInfectedWith() == context.serotype) {
            ++summary.true_prevalence;
            int daysSinceInf{nodePtr.getDaysSinceInfected()};
            summary.true_incidence += (daysSinceInf == 1) ? 1 : 0;
            bool detected{nodePtr.detectInfected()};
            summary.detected_prevalence += (detected) ? 1 : 0;
            summary.detected_incidence += (detected && daysSinceInf == int(m_config->m_detectionDelay) + 1) ? 1 : 0;
        }
    }
    return summary;
}

DailyCompartmentSums
Reporting::sumGlobalCompartments(const std::vector<Node> &nodes, DayContext const &context) const {
    DailyCompartmentSums compartmentSums(context);
    for (const auto &node : nodes) {
        // relying on being in same order
        auto nodeSums{node.getSumCompartments(context.serotype)};
        for (size_t i = 0; i < nodeSums.size(); ++i) {
            compartmentSums.sums.at(i) += nodeSums[i];
        }
    }
    return compartmentSums;
}

void Reporting::updateModelReport(std::vector<Node> &nodes, const std::string &run,
                                  const std::string &day) {
    // columns: run, day, serotype, value_type, value_sum
    for (const auto &serotype : m_config->m_serotypes) {
        DailySerotypeSummary summary{sumSerotypeStatistics(nodes, DayContext(run, day, serotype))};
        m_nodeInfectionsReport.push_back(summary);
    }
}



void Reporting::updateCompartmentSums(const std::vector<Node> &nodes, const std::string &runName,
                                      const std::string &dayName) {
    for (const auto &s : m_config->m_serotypes) {
        const auto globalSerotypeSums{sumGlobalCompartments(nodes, DayContext(runName, dayName, s))};
        m_compartmentSums.push_back(globalSerotypeSums);
    }
}

void Reporting::addInfectionEvent(const std::string &run, const std::string &day, const std::string &type,
                                  const std::string &serotype,
                                  const std::string &sourceID, const std::string &destinationID, bool prior) {
    m_infectionEvents.emplace_back(run, day, type, serotype, sourceID, destinationID, prior);
}

void Reporting::updateDetailedReport(const std::string &run, const std::string &day, const std::vector<Node> &nodes) {
    for (const auto &node : nodes) {
        if (node.isInfected()) {
            m_dailyInfectedNodes.emplace_back(run, day, node.getID(), node.getSerotypeInfectedWith(),
                                              node.getDaysSinceInfected());
        }
    }

}

void Reporting::writeReport(std::string_view report) const {
    if (report == "model-report") {
        writeTidyReport<DailySerotypeSummary>("-node-infection-report.csv",
                                              "run, day, serotype, value_type, value_sum\n",
                                              m_nodeInfectionsReport);
    }
    if (report == "compartment-sums") {
        writeTidyReport<DailyCompartmentSums>("-global-compartment-sums.csv",
                                              "run, day, serotype, compartment, value_sum\n", m_compartmentSums);
    }
    if (report == "infection-events") {
        writeTidyReport<InfectionEvent>("-infection-events.csv",
                                        "run, day, serotype, context, source_node, target_node, target_infected_prior\n",
                                        m_infectionEvents);
    }
    if (report == "detailed-report") {
        writeTidyReport<InfectedNodeSummary>("-detailed-report.csv",
                                             "run, day, serotype_infected_with, node_id, days_infected\n",
                                             m_dailyInfectedNodes);
    }
}


