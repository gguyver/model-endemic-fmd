#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <Common.h>
#include "Point.h"


std::vector<std::string> common::split(std::stringstream &inputString, char delim) {
    std::vector<std::string> tokens;
    std::string token;
    while (getline(inputString, token, delim)) {
        tokens.emplace_back(token);
    }
    return (tokens);
}


/// Based on algorithm described at http://www.johndcook.com/blog/cpp_expm1/
/// \param[in]	x	Exponent values
double common::oneMinusExp(double x) {
    if (x == 0.0) {
        return 0.0;
    } else if (std::abs(x) < 1e-5) {
        return -(x + 0.5 * x * x);
    } else {
        return -(exp(x) - 1.0);
    }
}

double common::euclideanDistance(const Point &pointA, const Point &pointB) {
    double vertical{pointA.getXCoord() - pointB.getXCoord()};
    double horizontal{pointA.getYCoord() - pointB.getYCoord()};
    return std::sqrt(std::pow(vertical, 2) + std::pow(horizontal, 2));
}

bool common::checkFilePathWorks(std::string &filePath) {
    std::fstream file{filePath.c_str()};
    return file.good();
}

double common::infectionProbabilty(double transmission, double susceptibility, double kernel) {
    return (1.0 - std::exp((-transmission * susceptibility * kernel)));
}

