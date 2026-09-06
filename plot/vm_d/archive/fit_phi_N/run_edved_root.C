#include <TSystem.h>
#include <iostream>
#include <stdexcept>

extern "C" {
void edved_c(const int* in, const int* ivm, const double* ei,
            const double* q2, const double* q0, const double* epsl,
            const double* t, double* crs, const double* params);
}

void run_edved_root(const char* library = "./libedved.so",
                    double photonEnergy = 8.3,
                    double minusT = 0.574,
                    double sgamman = 11.0,
                    double bgamman = 4.0,
                    double sphin = 30.0,
                    double bphin = 11.0) {
    if (gSystem->Load(library) < 0)
        throw std::runtime_error(std::string("Could not load ") + library);
    if (photonEnergy <= 0.0 || minusT < 0.0)
        throw std::invalid_argument("photonEnergy must be positive and minusT non-negative");

    const int ivm = 3;
    const double params[4] = {bgamman, sgamman, bphin, sphin};
    const double zero = 0.0;
    double crossSection = 0.0;
    int in = 1;
    edved_c(&in, &ivm, &zero, &zero, &zero, &zero, &zero, &crossSection, params);

    const double t = -minusT;
    in = 0;
    edved_c(&in, &ivm, &zero, &zero, &photonEnergy, &zero, &t, &crossSection, params);
    std::cout << "d sigma/dt = " << crossSection << " nb/GeV^2\n";
}
