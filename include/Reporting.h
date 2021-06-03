//
// Created by glen on 13/02/2020.
//

#ifndef MODEL_FMD_REPORTING_H
#define MODEL_FMD_REPORTING_H

// Forward declarations
class Config;

class Node;

#include <unordered_map>
#include <memory>
#include <iostream>
#include "Matrix.h"
#include "DiseaseCompartment.h"

struct DayContext {
    std::string const run;
    std::string const day;
    std::string const serotype;

    DayContext(std::string const &run, std::string const &day, std::string const &serotype) : run(run), day(day),
                                                                                              serotype(serotype) {}

    friend std::ostream &operator<<(std::ostream &os, DayContext const &context);
};

struct InfectionEvent {
    DayContext const context;
    std::string const type;
    std::string const sourceEpi;
    std::string const targetEpi;
    bool const targetInfectedPrior;

    InfectionEvent(const std::string &run, const std::string &day,
                   const std::string &eventType,
                   const std::string &serotype,
                   const std::string &sourceID, const std::string &targetID,
                   bool prior)
            : context(run, day, serotype), type(eventType), sourceEpi(sourceID), targetEpi(targetID),
              targetInfectedPrior(prior) {}

    friend std::ostream &operator<<(std::ostream &os, const InfectionEvent &event);
};

struct InfectedNodeSummary {
    DayContext const context;
    std::string const nodeID;
    int const daysInfected;

    InfectedNodeSummary(const std::string &run, const std::string &day, const std::string &ID,
                        const std::string &serotype, int daysInf)
            : context(run, day, serotype), nodeID(ID), daysInfected(daysInf) {}

    friend std::ostream &operator<<(std::ostream &os, const InfectedNodeSummary &ins);
};

struct DailySerotypeSummary {
    DayContext const context;
    int true_incidence{0};
    int true_prevalence{0};
    int detected_incidence{0};
    int detected_prevalence{0};

    explicit DailySerotypeSummary(DayContext const &context) : context(context) {}

    friend std::ostream &operator<<(std::ostream &os, DailySerotypeSummary const &dss);
};

struct DailyCompartmentSums {
    DayContext const context;
    std::array<int, 7> sums{};

    explicit DailyCompartmentSums(DayContext const &context) : context(context) {}

    friend std::ostream &operator<<(std::ostream &os, DailyCompartmentSums const &dcs);
};


class Reporting {
private:
    Config const *const m_config;
    std::vector<DailySerotypeSummary> m_nodeInfectionsReport;
    std::vector<DailyCompartmentSums> m_compartmentSums;
    std::vector<InfectionEvent> m_infectionEvents;
    std::vector<InfectedNodeSummary> m_dailyInfectedNodes;

    std::string m_outputDirectory;
    std::string const m_uniqueID;

    void initModelReport();

    void initCompartmentSumsReport();

    DailySerotypeSummary
    sumSerotypeStatistics(std::vector<Node> &nodes, DayContext const &context);

    [[nodiscard]] DailyCompartmentSums
    sumGlobalCompartments(const std::vector<Node> &nodes, DayContext const &context) const;

    template<class Event>
    void writeTidyReport(std::string const &fileSuffix, std::string const &columnHeaders,
                         std::vector<Event> const &eventsToWrite) const {
        std::string outputFile{m_outputDirectory + m_uniqueID + fileSuffix};
        std::ofstream outfile(outputFile.c_str());
        if (outfile.is_open()) {
            outfile << columnHeaders;
            for (const auto &event : eventsToWrite) {
                outfile << event;
            }
            outfile.close();
        } else {
            std::cerr << "\nError: Cannot write events to: " << outputFile;
        }
    }

public:
    explicit Reporting(Config *configPtr);


    void updateModelReport(std::vector<Node> &nodes, const std::string &run,
                           const std::string &day);

    void updateCompartmentSums(const std::vector<Node> &nodes, const std::string &runName,
                               const std::string &dayName);

    void addInfectionEvent(const std::string &run, const std::string &day, const std::string &type,
                           const std::string &serotype,
                           const std::string &sourceID, const std::string &destinationID, bool prior);

    void updateDetailedReport(const std::string &run, const std::string &day, const std::vector<Node> &nodes);

    void writeReport(std::string_view report) const;

};


#endif //MODEL_FMD_REPORTING_H
