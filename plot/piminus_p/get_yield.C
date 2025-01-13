using namespace std;
using namespace RooFit;

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
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

double yield_fit(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);

int get_yield(string Reaction)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_piminus_p_recon_%s.root",Reaction.c_str());
    string txt_name = Form("output/yield_piminus_p_recon_%s.txt",Reaction.c_str());
    string pdf_name = Form("output/plot_piminus_p_recon_%s.pdf",Reaction.c_str());

    TFile *input_file = new TFile(input_name.c_str(), "read");
    TTree *input_tree = (TTree*) input_file->Get("filteredtree_piminus_p_recon");
    FILE *txt_file = fopen(txt_name.c_str(),"w");
    TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);

    // double energy_edges[] = {6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5};
    double energy_edges[] = {6.0, 6.5};
    double theta_edges[] = {20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
    int energy_bins = sizeof(energy_edges) / sizeof(energy_edges[0]) - 1;
    int theta_bins = sizeof(theta_edges) / sizeof(theta_edges[0]) - 1;

    double this_yield;

    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            cout << "Energy: " << energy_edges[i] << " GeV" << endl;
            cout << "Theta: " << theta_edges[j] << " deg" << endl;
            this_yield = yield_count(energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], input_tree);
            fprintf(txt_file, "%3.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], this_yield);
            canvas->Update();
            canvas->Print((pdf_name+"(").c_str());
            canvas->Clear();
        }
    }
    canvas->Print((pdf_name+")").c_str());

    return 0;

}

double yield_fit_roofit(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{

    TH1F *hist = new TH1F("hist", "hist", 100, 0.4, 1.4);
    input_tree->Draw("miss_pminus_kin>>hist", Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));

    RooRealVar miss_pminus_kin("miss_pminus_kin", "miss_pminus_kin", 0.7, 1.2);

    RooRealVar mu("mu", "mu", 0.9, 0.8, 1.0);
    RooRealVar sigma("sigma", "sigma", 0.1, 0.01, 0.4);
    RooGaussian gauss("gauss", "gauss", miss_pminus_kin, mu, sigma);

    RooRealVar c0("c0", "c0", 0.5, 0, 2);
    RooRealVar c1("c1", "c1", -1, -5, 5);
    RooPolynomial poly("poly", "poly", miss_pminus_kin, RooArgList(c0, c1));

    RooRealVar fsig("fsig", "signal fraction", 0.9, 0.01, 1.0);
    RooAddPdf model("model", "model", RooArgList(gauss, poly), fsig);

    RooDataHist data("data", "data", miss_pminus_kin, hist);
    model.fitTo(data);

    string hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high);
    RooPlot *frame = miss_pminus_kin.frame(Title(hist_name.c_str()));
    data.plotOn(frame);
    model.plotOn(frame);
    frame->SetTitle(hist_name.c_str());
    frame->Draw();

    double yield = fsig.getVal() * hist->Integral();
    delete hist;
    return yield;
}

double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    input_tree->Draw(Form("miss_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(miss_pminus_kin>0.7 && miss_pminus_kin<1.1 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    double yield = hist->Integral();
    return yield;
}