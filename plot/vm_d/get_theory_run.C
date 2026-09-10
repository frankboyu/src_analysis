#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <string>
using namespace std;

R__LOAD_LIBRARY(/work/halld2/home/boyu/src_analysis/plot/vm_d/lib_edved.so)

extern "C" {
void edved_c(const int* in, const int* ivm, const double* ei,
            const double* q2, const double* q0, const double* epsl,
            const double* t, double* crs, const double* params);
}

void get_theory_run(string tag)
{
    const int ivm = 3;

    // Define the parameters for the edved_c function
    double photonEnergy, minusT;
    double sphin, bphin;
    double sphin_range[3] {20.0, 40.0, 0.2}; // start, end, step
    double bphin_range[3] {5.0, 15.0, 0.1};  // start, end, step
    if (tag.find("nominal") == string::npos)  // if systematic study
    {
        // only need to search for minimum instead of running on the full range
        sphin_range[0] = 25.0;
        sphin_range[1] = 35.0;
        bphin_range[0] = 7.0;
        bphin_range[1] = 13.0;
    }

    double sgamman  = 11.0;
    double bgamman  = 4.4;
    double agamman  = 0.0;
    double aphin    = 0.0;
    double t2term   = 0.0;
    double wfflag   = 1.0;

    // Parse both integer and decimal-point values, such as sgamman_11.5.
    auto readTagValue = [&tag](const string& key, double& value) {
        const size_t pos = tag.find(key);
        if (pos != string::npos)
            value = stod(tag.substr(pos + key.size()));
    };

    readTagValue("sgamman_", sgamman);
    readTagValue("bgamman_", bgamman);
    readTagValue("agamman_", agamman);
    readTagValue("aphin_",   aphin);
    readTagValue("t2term_",  t2term);
    readTagValue("wfflag_",  wfflag);

    // Prepare output file for results
    std::ofstream results("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/theory/theory_" + tag + ".txt");
    if (!results)
        throw std::runtime_error("Could not open /work/halld2/home/boyu/src_analysis/plot/vm_d/output/theory/theory_" + tag + ".txt");

    // Read kinematic points from input file
    std::ifstream kinematics_list("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/table_vm_d_dsdt.txt");
    std::string line;
    std::vector<double> energy_list;
    std::vector<double> minust_list;

    while (std::getline(kinematics_list, line)) {
        std::istringstream values(line);
        double energy;
        double minust;
        double unused[5];
        if (values >> energy >> minust >> unused[0] >> unused[1] >> unused[2] >> unused[3] >> unused[4]) {
            if (energy > 5.8 && energy < 7.8)
                energy_list.push_back(6.9);
            else if (energy > 7.8 && energy < 8.8)
                energy_list.push_back(8.3);
            else if (energy > 8.8 && energy < 10.8)
                energy_list.push_back(9.7);
            else
                std::cerr << "Warning: Energy value " << energy << " is out of expected ranges." << std::endl;
            minust_list.push_back(minust);
        }
    }
    kinematics_list.close();

    results << "Differential cross sections calculated using the edved_c function in units of nb/GeV^2.\n"
            << "s_gamma_n = " << sgamman << " mb\n"
            << "b_gamma_n = " << bgamman << " GeV^-2\n"
            << "a_gamma_n = " << agamman << "\n"
            << "a_phi_n = "   << aphin   << "\n"
            << "t2term = "    << t2term  << "\n"
            << "wfflag = "    << wfflag  << "\n";

    results << std::setw(20) << " ";
    results << std::setw(20) << " ";
    results << std::setw(20) << "E_gamma (GeV)";
    for (size_t i = 0; i < energy_list.size(); ++i)
        results << std::setw(20) << std::fixed << std::setprecision(1) << energy_list[i];
    results << '\n';

    results << std::setw(20) << " ";
    results << std::setw(20) << " ";
    results << std::setw(20) << "-t (GeV^2)";
    for (size_t i = 0; i < minust_list.size(); ++i)
        results << std::setw(20) << std::setprecision(4) << minust_list[i];
    results << '\n';

    results << std::setw(20) << "s_phi_n (mb)" << std::setw(20) << "b_phi_n (GeV^-2)" << '\n';

    // Loop over sphin and bphin values
    for (sphin = sphin_range[0]; sphin < sphin_range[1]; sphin += sphin_range[2])
    {
        cout << "sphin = " << sphin << std::endl;
        for (bphin = bphin_range[0]; bphin < bphin_range[1]; bphin += bphin_range[2])
        {
            // Initialize parameters for calculations
            const double params[8] = {
                bgamman, sgamman, agamman, bphin,
                sphin, aphin, t2term, wfflag
            };
            const double zero = 0.0;
            double crossSection = 0.0;
            int in = 1;
            edved_c(&in, &ivm, &zero, &zero, &zero, &zero, &zero, &crossSection, params);

            // Loop over kinematic points and calculate cross sections
            in = 0;
            std::vector<double> crossSection_list;
            for (size_t i = 0; i < energy_list.size(); ++i) {
                photonEnergy = energy_list[i];
                minusT = minust_list[i];
                const double t = -minusT;
                edved_c(&in, &ivm, &zero, &zero, &photonEnergy, &zero, &t, &crossSection, params);
                crossSection_list.push_back(crossSection);
            }

            // Write results to output file
            results << std::setw(20) << std::setprecision(1) << sphin << std::setw(20) << std::setprecision(1) << bphin << std::setw(20) << " ";
            for (size_t i = 0; i < crossSection_list.size(); ++i)
                results << std::setw(20) << std::setprecision(8) << crossSection_list[i];
            results << '\n';
        }
    }

    results.close();
}