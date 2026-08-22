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

    // Build the 4-parameter array expected by the Fortran wrapper
    double full_params[4];
    full_params[0] = 4.5;      // bgamman (fixed)
    full_params[1] = 11.0;     // sgamman (fixed)
    full_params[2] = par[0];   // bphin   (MINUIT par 0)
    full_params[3] = par[1];   // sphin   (MINUIT par 1)

    // STEP A: Initialize Fortran Common Blocks with current parameters
    int in = 1; 
    edved_c(&in, &ivm, &dummy_val, &dummy_val, &dummy_val, &dummy_val, 
            &dummy_val, &dummy_crs, full_params);

    // STEP B: Calculate Chi-Square across all data points
    in = 0; 
    for (size_t i = 0; i < exp_crs.size(); ++i) {
        double theory_crs = 0.0;
        
        // Compute theoretical cross-section using full_params
        edved_c(&in, &ivm, &data_ei[i], &data_q2[i], &data_q0[i], 
                &data_epsl[i], &data_t[i], &theory_crs, full_params);
        
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
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        double photon_energy, minus_t, cross_section, error;

        if (iss >> photon_energy >> minus_t >> cross_section >> error) {
            
            // Photoproduction kinematics
            data_q0.push_back(photon_energy);
            data_t.push_back(-minus_t);       // Fortran expects t to be negative
            data_q2.push_back(0.0);           
            data_epsl.push_back(0.0);         
            data_ei.push_back(0.0);           
            
            // Observables
            exp_crs.push_back(cross_section);
            exp_err.push_back(error);
            
            point_count++;
        }
    }
    
    std::cout << "Successfully loaded " << point_count << " data points." << std::endl;
    return point_count > 0;
}

int main(int argc, char** argv) {
    std::string filename = "data.txt";
    if (argc > 1) {
        filename = argv[1];
    }

    if (!LoadData(filename)) {
        return 1;
    }

    // Initialize TMinuit with 2 parameters instead of 4
    TMinuit *gMinuit = new TMinuit(2);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;
    
    // Set error definition to 1.0 for Chi-square fits
    arglist[0] = 1.0; 
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Define the 2 MINUIT parameters
    gMinuit->mnparm(0, "bphin",   10.0,  0.1, 0, 0, ierflg); 
    gMinuit->mnparm(1, "sphin",   30.0,  0.1, 0, 0, ierflg); 

    // Run MIGRAD (the minimizer)
    arglist[0] = 5000; 
    arglist[1] = 1.0;  
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    return 0;
}