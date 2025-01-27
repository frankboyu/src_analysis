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

double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_roofit_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_roofit_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);

TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);

int get_yield(string Reaction)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_piminus_p_recon_%s.root",Reaction.c_str());
    string txt_name = Form("output/yield_piminus_p_recon_%s.txt",Reaction.c_str());
    string pdf_name = Form("output/plot_piminus_p_recon_%s.pdf",Reaction.c_str());

    TFile *input_file = new TFile(input_name.c_str(), "read");
    TTree *input_tree = (TTree*) input_file->Get("filteredtree_piminus_p_recon");
    FILE *txt_file = fopen(txt_name.c_str(),"w");
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    vector<double> energy_edges;
    vector<double> theta_edges;
    ifstream energy_file("input/bins_energy.txt");
    ifstream theta_file("input/bins_theta.txt");
    double value;
    while (energy_file >> value)
    {
        energy_edges.push_back(value);
    }
    while (theta_file >> value)
    {
        theta_edges.push_back(value);
    }
    energy_file.close();
    theta_file.close();
    int energy_bins = energy_edges.size() - 1;
    int theta_bins = theta_edges.size() - 1;

    double this_yield;

    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            cout << "Energy: " << energy_edges[i] << " GeV" << endl;
            cout << "Theta: " << theta_edges[j] << " deg" << endl;
            this_yield = yield_root_quadratic_bkg(energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], input_tree);
            fprintf(txt_file, "%3.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], this_yield);
            canvas->Update();
            canvas->Print((pdf_name+"(").c_str());
            canvas->Clear();
        }
    }
    canvas->Print((pdf_name+")").c_str());

    return 0;

}

double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    input_tree->Draw(Form("miss_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    double yield = hist->Integral();
    return yield;
}

double yield_root_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    input_tree->Draw(Form("miss_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->GetEntries() < 200)
        return 0;

    TF1 *model = new TF1("model", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9, 0.9, 0.05, 20, -20);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("same");

    TF1 *gauss = new TF1("gauss", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0.5, 1.3);
    gauss->SetParameters(model->GetParameter(0), model->GetParameter(1), model->GetParameter(2));

    double yield = gauss->Integral(0.5, 1.3) / hist->GetBinWidth(1);
    return yield;
}

double yield_root_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    input_tree->Draw(Form("miss_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->GetEntries() < 200)
        return 0;

    TF1 *model = new TF1("model", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9, 0.9, 0.05, 20, -20);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    // model->SetParLimits(3, 0, 100);
    // model->SetParLimits(4, -100, 100);
    hist->Fit(model, "R");
    model->Draw("same");

    TF1 *gauss = new TF1("gauss", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0.5, 1.3);
    gauss->SetParameters(model->GetParameter(0), model->GetParameter(1), model->GetParameter(2));

    double yield = gauss->Integral(0.5, 1.3) / hist->GetBinWidth(1);
    return yield;
}

double yield_root_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    input_tree->Draw(Form("miss_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->GetEntries() < 200)
        return 0;

    TF1 *model = new TF1("model", "[0]*exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x + [5]*x*x", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9, 0.9, 0.05, 0, 0, 0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    // model->SetParLimits(3, 0, 100);
    // model->SetParLimits(4, -100, 100);
    hist->Fit(model, "R");
    model->Draw("same");

    TF1 *gauss = new TF1("gauss", "[0]*exp(-0.5*((x-[1])/[2])**2)", 0.5, 1.3);
    gauss->SetParameters(model->GetParameter(0), model->GetParameter(1), model->GetParameter(2));

    double yield = gauss->Integral(0.5, 1.3) / hist->GetBinWidth(1);
    return yield;
}

double yield_roofit_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
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

double yield_roofit_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{

    TH1F *hist = new TH1F("hist", "hist", 100, 0.4, 1.4);
    input_tree->Draw("miss_pminus_kin>>hist", Form("accidweight*(miss_pminus_kin>0.5 && miss_pminus_kin<1.3 && beam_E_kin>%f && beam_E_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));

    RooRealVar miss_pminus_kin("miss_pminus_kin", "miss_pminus_kin", 0.7, 1.2);

    RooRealVar mu("mu", "mu", 0.9, 0.8, 1.0);
    RooRealVar sigma("sigma", "sigma", 0.1, 0.01, 0.4);
    RooGaussian gauss("gauss", "gauss", miss_pminus_kin, mu, sigma);

    RooRealVar c0("c0", "c0", 0.5, 0, 2);
    RooRealVar c1("c1", "c1", -1, -5, 5);
    RooRealVar c2("c2", "c2", 0.5, 0, 2);
    RooPolynomial poly("poly", "poly", miss_pminus_kin, RooArgList(c0, c1, c2));

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