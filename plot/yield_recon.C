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

    RooRealVar MissingPMinus("MissingPMinus", "MissingPMinus", 0.7, 1.2);

    RooRealVar mu("mu", "mu", 0.9, 0.8, 1.0);
    RooRealVar sigma("sigma", "sigma", 0.1, 0.01, 0.4);
    RooGaussian gauss("gauss", "gauss", MissingPMinus, mu, sigma);

    RooRealVar c0("c0", "c0", 0.5, 0, 2);
    RooRealVar c1("c1", "c1", -1, -3, 3);
    RooRealVar c2("c2", "c2", 1, 0, 15);
    RooPolynomial poly("poly", "poly", MissingPMinus, RooArgList(c0, c1, c2));

    RooRealVar fsig("fsig", "signal fraction", 0.9, 0.01, 1.0);
    RooAddPdf model("model", "model", RooArgList(gauss, poly), fsig);

    RooDataHist data("data", "data", MissingPMinus, hist);
    model.fitTo(data);

    RooPlot *frame = MissingPMinus.frame(Title("MissingPMinus"));
    data.plotOn(frame);
    model.plotOn(frame);
    frame->Draw();

    double yield = fsig.getVal() * hist->Integral();
    return yield;
}

int yield_recon()
{
    TFile *input_file = new TFile("/work/halld2/home/boyu/src_analysis/selection/output/flattree_piminus_p_2H_sim_samecut.root", "read");
    FILE *yieldfile = fopen("output/yield_piminus_p_2H_sim_samecut.txt","w");
    TTree *input_tree = (TTree*) input_file->Get("flattree_piminus_p_2H_sim");
    TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);

    double energy_low  = 6.0;
    double energy_high = 10.5;
    double energy_step = 0.5;
    int energy_bins    = int((energy_high - energy_low) / energy_step);
    double theta_low   = 20;
    double theta_high  = 50;
    double theta_step  =  5;
    int theta_bins     = int((theta_high - theta_low) / theta_step);

    double this_energy, this_theta, this_yield;

    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            this_energy = energy_low + i * energy_step;
            this_theta = theta_low + j * theta_step;
            cout << "Energy: " << this_energy << " GeV" << endl;
            cout << "Theta: " << this_theta << " deg" << endl;
            this_yield = get_yield(this_energy, this_energy + energy_step, this_theta, this_theta + theta_step, input_tree);
            fprintf(yieldfile, "%f %f %f %f %f\n", this_energy, this_energy + energy_step, this_theta, this_theta + theta_step, this_yield);
            canvas->Update();
            canvas->Print("output/fit_piminus_p_2H_sim_samecut.pdf(");
            canvas->Clear();
        }
    }
    canvas->Print("output/fit_piminus_p_2H_sim_samecut.pdf)");

    return 0;

}
