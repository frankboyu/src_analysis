#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TMinuit.h>

// 1. Declare the Fortran bridge wrapper
extern "C" {
    void edved_c(const int* in, const int* ivm, const double* ei, 
                 const double* q2, const double* q0, const double* epsl, 
                 const double* t, double* crs, const double* params);
}

// 2. Global structures to hold your experimental data
std::vector<double> data_ei, data_q2, data_q0, data_epsl, data_t;
std::vector<double> exp_crs, exp_err;

// 3. The MINUIT FCN (Objective Function)
void fcn(int &npar, double *gin, double &f, double *par, int iflag) {
    double chisq = 0.0;
    int ivm = 3;          // 3 = phi meson
    double dummy_val = 0.0;
    double dummy_crs = 0.0;

    // STEP A: Initialize Fortran Common Blocks with current MINUIT parameters
    int in = 1; 
    edved_c(&in, &ivm, &dummy_val, &dummy_val, &dummy_val, &dummy_val, 
            &dummy_val, &dummy_crs, par);

    // STEP B: Calculate Chi-Square across all data points
    in = 0; 
    for (size_t i = 0; i < exp_crs.size(); ++i) {
        double theory_crs = 0.0;
        
        // Compute theoretical cross-section for this kinematic data point
        edved_c(&in, &ivm, &data_ei[i], &data_q2[i], &data_q0[i], 
                &data_epsl[i], &data_t[i], &theory_crs, par);
        
        // Calculate standard residual
        double delta = (exp_crs[i] - theory_crs) / exp_err[i];
        chisq += delta * delta;
    }
    
    // Return Chi-Square to MINUIT
    f = chisq;
}

// 4. Function to load data from text file
bool LoadData(const std::string& filename) {
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    int point_count = 0;
    
    // Read file line by line
    while (std::getline(file, line)) {
        // Skip empty lines or comment lines starting with #
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double photon_energy, minus_t, cross_section, error, error2, error3, error4;

        // Parse the 4 columns: photon energy, -t, cross section, error
        if (iss >> photon_energy >> minus_t >> cross_section >> error >> error2 >> error3 >> error4) {
            
            // Photoproduction kinematics
            data_q0.push_back(photon_energy);
            data_t.push_back(-minus_t);       // Fortran expects t to be negative!
            data_q2.push_back(0.0);           // q2 = 0 for photoproduction
            data_epsl.push_back(0.0);         // epsl = 0 for photoproduction
            data_ei.push_back(0.0);           // ei is ignored, pass 0.0
            
            // Observables
            exp_crs.push_back(cross_section);
            exp_err.push_back(error4);
            
            point_count++;
        }
    }
    
    std::cout << "Successfully loaded " << point_count << " data points." << std::endl;
    return point_count > 0;
}

int main(int argc, char** argv) {
    // Check if user provided a data file argument
    std::string filename = "output/table_vm_d_dsdt.txt"; // Default filename
    if (argc > 1) {
        filename = argv[1];
    }

    // Load the data
    if (!LoadData(filename)) {
        return 1; // Exit if data loading fails
    }

    // Initialize TMinuit with 4 parameters
    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;
    
    // Set error definition to 1.0 for Chi-square fits
    arglist[0] = 1.0; 
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Define Parameters: (Index, Name, Initial Guess, Step Size, Min, Max)
    // Note: Min=0 and Max=0 tells MINUIT not to put boundaries on the parameter
    // The order exactly matches the Fortran wrapper `params_c` array
    gMinuit->mnparm(0, "bgamman",  4.5,  0.1, 0, 0, ierflg); // par[0]
    gMinuit->mnparm(1, "sgamman", 11.0,  0.1, 0, 0, ierflg); // par[1]
    gMinuit->mnparm(2, "bphin",   10.0,  0.1, 0, 0, ierflg); // par[2]
    gMinuit->mnparm(3, "sphin",   30.0,  0.1, 0, 0, ierflg); // par[3]

    // Run MIGRAD (the minimizer)
    arglist[0] = 5000; // Maximum number of function calls
    arglist[1] = 1.0;  // Tolerance
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Run MINOS (to get better asymmetric error bars, optional but recommended)
    // gMinuit->mnexcm("MINOS", arglist, 0, ierflg);

    return 0;
}