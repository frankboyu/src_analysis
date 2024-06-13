using namespace std;
using namespace RooFit;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMath.h>
#include <TFitResult.h>

double get_yield(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    
    TH1F *hist = new TH1F("hist", "hist", 100, 0.4, 1.4);
    input_tree->Draw("MissingPMinus>>hist", Form("WeightFactor*(MissingPMinus>0.5 && MissingPMinus<1.3 && BeamEnergy>%f && BeamEnergy<%f && thetaCM>%f && thetaCM<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    double yield = hist->Integral();
    return yield;
}

int yield_thrown()
{
    TFile *input_file = new TFile("thrown_4He.root", "read");
	TTree *input_tree = (TTree*) input_file->Get("piminus_p_4He_thrown");
    TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);
    FILE *yieldfile = fopen("output/yield_piminus_p_4He_thrown.txt","w");
    // FILE *errorfile = fopen("output/error_piminus_p_4He_thrown.txt","w");
    // FILE *parasfile = fopen("output/paras_piminus_p_4He_thrown.txt","w");

    double energy_low  = 6.0;
    double energy_high = 10.5;
    double energy_step = 0.5;
    int energy_bins    = int((energy_high - energy_low) / energy_step);
    double theta_low   = 20;
    double theta_high  = 40;
    double theta_step  = 5;
    int theta_bins     = int((theta_high - theta_low) / theta_step);

    double this_energy, this_theta, this_yield;

    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            this_energy = energy_low + i * energy_step;
            this_theta = theta_low + j * theta_step;
            this_yield = get_yield(this_energy, this_energy + energy_step, this_theta, this_theta + theta_step, input_tree);
            fprintf(yieldfile, "%f %f %f %f %f\n", this_energy, this_energy + energy_step, this_theta, this_theta + theta_step, this_yield);
            canvas->Update();
            canvas->Print("output/fit_piminus_p_4He_thrown.pdf(");
            canvas->Clear();
        }
    }
    canvas->Print("output/fit_piminus_p_4He_thrown.pdf)");

    return 0;

}