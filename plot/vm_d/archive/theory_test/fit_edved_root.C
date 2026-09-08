#include <TMinuit.h>
#include <TSystem.h>
#include <TString.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
void edved_c(const int* in, const int* ivm, const double* ei,
            const double* q2, const double* q0, const double* epsl,
            const double* t, double* crs, const double* params);
}

namespace {
struct Point {
    double photonEnergy;
    double minusT;
    double crossSection;
    double statisticalError;
    double p2pError;
    double normalizationError;
    double totalSystematicError;
};

std::vector<Point> points;
int ivm = 3;

double chiSquare(double* parameters) {
    const int initialize = 1;
    const double zero = 0.0;
    double unusedCrossSection = 0.0;
    double fortranParameters[8] = {
        4.5,
        11.0,
        0.0,
        parameters[0],
        parameters[1],
        0.0,
        0.0,
        1.0
    };
    edved_c(&initialize, &ivm, &zero, &zero, &zero, &zero, &zero,
            &unusedCrossSection, fortranParameters);

    const int calculate = 0;
    constexpr double normShift = 0.0;
    double chi2 = 0.0;
    for (const Point& point : points) {
        const double q2 = 0.0;
        const double epsl = 0.0;
        const double ei = 0.0;
        const double t = -point.minusT;
        double theory = 0.0;
        edved_c(&calculate, &ivm, &ei, &q2, &point.photonEnergy,
            &epsl, &t, &theory, fortranParameters);
        const double pointwiseError = std::sqrt(
            point.statisticalError * point.statisticalError +
            point.p2pError * point.p2pError);
        const double theoryNormalizationError = 0.05 * point.crossSection;
        const double combinedNormalizationError = std::sqrt(
            point.normalizationError * point.normalizationError +
            theoryNormalizationError * theoryNormalizationError);
        const double shiftedData = point.crossSection +
                       normShift * combinedNormalizationError;
        const double residual = (shiftedData - theory) / pointwiseError;
        chi2 += residual * residual;
    }
    return chi2;
}

void fcn(int&, double*, double& value, double* parameters, int) {
    value = chiSquare(parameters);
}

bool loadData(const char* filename) {
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "Could not open data file: " << filename << '\n';
        return false;
    }

    points.clear();
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream stream(line);
        Point point{};
        if (stream >> point.photonEnergy >> point.minusT >> point.crossSection
                   >> point.statisticalError >> point.p2pError
                   >> point.normalizationError >> point.totalSystematicError) {
            if (point.statisticalError > 0.0 && point.p2pError >= 0.0 &&
                point.normalizationError >= 0.0 && point.totalSystematicError >= 0.0) {
                points.push_back(point);
            }
        }
    }
    std::cout << "Loaded " << points.size() << " data points\n";
    return !points.empty();
}
}

void fit_edved_root(const char* library = "./libedved.so",
                    const char* dataFile = "output/table_vm_d_dsdt.txt") {
    if (gSystem->Load(library) < 0)
        throw std::runtime_error(std::string("Could not load ") + library);
    if (!loadData(dataFile)) return;

    TMinuit minuit(2);
    minuit.SetFCN(fcn);
    double command[2] = {1.0, 0.0};
    int status = 0;
    minuit.mnexcm("SET ERR", command, 1, status);

    double strategy = 2.0;
    minuit.mnexcm("SET STRATEGY", &strategy, 1, status);

    // Fortran receives: bgamman, sgamman, bphin, sphin.
    minuit.mnparm(0, "bphin", 11.0, 0.3, 0.01, 30.0, status);
    minuit.mnparm(1, "sphin", 30.0, 1.0, 0.01, 200.0, status);

    command[0] = 5000.0;
    command[1] = 1.0;
    int simplexStatus = 0;
    minuit.mnexcm("SIMPLEX", command, 2, simplexStatus);
    std::cout << "SIMPLEX status = " << simplexStatus << '\n';

    minuit.mnexcm("MIGRAD", command, 2, status);

    std::cout << "MIGRAD status = " << status << '\n';

    double minosCommand[2] = {0.0, 0.0};
    int minosStatus = 0;
    minuit.mnexcm("MINOS", minosCommand, 0, minosStatus);
    std::cout << "MINOS status = " << minosStatus << '\n';

    if (status != 0 && minosStatus == 0) {
        std::cout << "Using the SIMPLEX/MINOS result; MIGRAD covariance was not valid.\n";
    }
    if (minosStatus != 0) {
        std::cerr << "MINOS failed; do not use the reported symmetric errors.\n";
    }

    for (int i = 0; i < 2; ++i) {
        TString name;
        double value = 0.0;
        double error = 0.0;
        double lowerLimit = 0.0;
        double upperLimit = 0.0;

        minuit.mnpout(i, name, value, error,
                      lowerLimit, upperLimit, status);
        std::cout << name << ": " << value << " +/- " << error << '\n';

        if (status == 0) {
            double errorLow = 0.0;
            double errorHigh = 0.0;
            double errorParabolic = 0.0;
            double gcc = 0.0;
            minuit.mnerrs(i, errorLow, errorHigh, errorParabolic, gcc);
            std::cout << "  MINOS: -" << std::abs(errorLow)
                      << " / +" << errorHigh << '\n';
        }
    }

    minuit.mnmatu(1);
}

void scan_edved_step_sizes(const char* library = "./libedved.so",
                           double bphin = 11.0,
                           double sphin = 30.0) {
    if (gSystem->Load(library) < 0)
        throw std::runtime_error(std::string("Could not load ") + library);
    if (points.empty() && !loadData("output/table_vm_d_dsdt.txt")) return;

    double base[2] = {bphin, sphin};
    const double baseValue = chiSquare(base);
    const double repeatedValue = chiSquare(base);
    std::cout << std::setprecision(12)
              << "chi2(base) = " << baseValue
              << ", repeat difference = " << repeatedValue - baseValue << '\n';
    std::cout << "parameter step, bphin derivative, bphin curvature, "
                 "sphin derivative, sphin curvature\n";

    for (double step : {1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001}) {
        double bphinPlus[2] = {bphin + step, sphin};
        double bphinMinus[2] = {bphin - step, sphin};
        double sphinPlus[2] = {bphin, sphin + step};
        double sphinMinus[2] = {bphin, sphin - step};

        const double bphinDerivative =
            (chiSquare(bphinPlus) - chiSquare(bphinMinus)) / (2.0 * step);
        const double bphinCurvature =
            (chiSquare(bphinPlus) - 2.0 * baseValue + chiSquare(bphinMinus)) /
            (step * step);
        const double sphinDerivative =
            (chiSquare(sphinPlus) - chiSquare(sphinMinus)) / (2.0 * step);
        const double sphinCurvature =
            (chiSquare(sphinPlus) - 2.0 * baseValue + chiSquare(sphinMinus)) /
            (step * step);

        std::cout << step << ' '
                  << bphinDerivative << ' ' << bphinCurvature << ' '
                  << sphinDerivative << ' ' << sphinCurvature << '\n';
    }
}
