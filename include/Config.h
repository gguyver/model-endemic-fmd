//
// Created by glen on 03/02/2020.
//

#ifndef MODEL_SINGLESEROTYPE_CONFIG_H
#define MODEL_SINGLESEROTYPE_CONFIG_H

// FORWARD DECLARATIONS

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <iostream>
#include "yaml.h"
#include "Reporting.h"


struct VaccineEfficacy {
    double meanEfficacy{0};
    double stdevEfficacy{0};
};

struct Vaccine {
    std::map<std::string, VaccineEfficacy> serotypeEfficacies;
};

struct ConfigSerotypeParameters {
    double beta; // infectiousness
    double symptomaticRate; // rate progress to infected
    double recoveryRate; // rate recover
    double waningNatural; // natural immunity duration
    double maternalImmunityDecayRate;
    double waningVaccine; // vaccine immunity duration
    double associatedMortalityRate;
    double transmission;
    double susceptibility;
    double proportionCarriers;
    double waningCarriers;
    double carrierTransmission;
};

struct ShipmentRecord {
    int day;
    std::string sourceNodeID;
    std::string targetNodeID;
    int numMoved;

    ShipmentRecord(int d, std::string sid, std::string tid, int size) : day(d), sourceNodeID(std::move(sid)),
                                                                        targetNodeID(std::move(tid)), numMoved(size) {}
};

class Config {
public:
    explicit Config(const std::string &configFile);

    std::unique_ptr<Reporting> m_reporting;
    std::unique_ptr<std::mt19937> m_pRNG;
    std::unique_ptr<std::vector<ShipmentRecord>> m_shipmentRecords;

    // PARAMETERS
    std::string m_nodeFile, m_shipmentFile, m_outputDirectory, m_uniqueID;
    bool m_outputDetailedReports, m_outputInfectionEvents, m_outputCompartmentSums;
    unsigned int m_verbosity{1};
    bool m_seedsAreFixed; // seed specified or random?
    int m_numberOfRandomSeeds{0};
    std::vector<std::string> m_seedIDs; // IDs of nodes that are to be seeded
    std::string m_seededSerotype;

    unsigned int m_numModelRuns;
    int m_startDay, m_endDay, m_numDays;
    unsigned int m_maxNodeInfectionDuration;
    bool m_burnIn, m_runShipments;
    int m_burnInDuration;

    // Infection
    std::vector<std::string> m_serotypes;
    double m_dailyBirthRate;
    double m_dailyNaturalMortality;
    bool m_maternalImmunityIsActive;
    std::map<std::string, ConfigSerotypeParameters> m_serotypeInfectionParameters;
    double m_maxTransmission{0.0};
    double m_maxSusceptibility{0.0};
    // Kernel
    double m_a, m_scale, m_shape;
    bool m_forceInfection;
    double m_forcingRate;

    double m_shipmentFomiteSpreadProbability;
    // Control
    double m_detectionRate;
    unsigned int m_detectionDelay; // delay between Node becoming infected and being able to be detected

    Vaccine m_vaccine;
    double m_averageVaccineDuration;
    bool m_reactiveVaccinationImplemented;
    double m_reactiveVaccinationRadius;
    double m_reactiveVaccinationCoverage;
    int m_dailyVaccinationCapacity;

    bool m_movementBansImplemented;
    double m_movementBansRadius;
    double m_movementBansCompliance;
    unsigned int m_movementBanDuration;

    bool m_massVaccinationImplemented;
    double m_massVaccinationCoverage;
    unsigned int m_massVaccinationInterval;

    void writeReports() const;

private:

    YAML::Node m_configYAML;

    template<typename T>
    T GetYAMLParam(const std::vector<std::string> &keywords) const {
        std::string path{"config.yaml"};
        try {
            // .reset() does what copy assignment should do.
            YAML::Node current;
            current.reset(m_configYAML);
            for (const auto &key : keywords) {
                path += (" -> " + key);
                if (current[key]) {
                    current.reset(current[key]);
                } else {
                    throw std::exception();
                }
            }
            return current.as<T>();
        } catch (...) {
            throw std::runtime_error(std::string("Cannot find/convert parameter: ") + path +
                                     std::string("\n Are you sure it is present/correctly formatted?"));
        }
    }

    void translateYAML();

    void validateParameters();

    void initReporting();

    void readShipmentFile();


};


#endif //MODEL_SINGLESEROTYPE_CONFIG_H
