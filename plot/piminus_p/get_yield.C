using namespace std;
using namespace RooFit;

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>

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

double pi = TMath::Pi();

double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_root_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);

TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);

int get_yield(string Reaction)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_piminus_p_%s_all.root",Reaction.c_str());
    string txt_name = Form("output/yield_piminus_p_%s.txt",Reaction.c_str());
    string pdf_name = Form("output/plot_piminus_p_%s.pdf",Reaction.c_str());

    TFile *input_file = new TFile(input_name.c_str(), "read");
    string input_tree_name;
    if (Reaction.find("recon") != string::npos)
        input_tree_name = "filteredtree_piminus_p_recon";
    else if (Reaction.find("thrown") != string::npos)
        input_tree_name = "filteredtree_piminus_p_thrown";
    TTree *input_tree = (TTree*) input_file->Get(input_tree_name.c_str());
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
            if (Reaction.find("sim") != string::npos || Reaction.find("thrown") != string::npos)
            {
                this_yield = yield_root_no_bkg(energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], input_tree);
            }
            else if (Reaction.find("data") != string::npos)
            {
                if (Reaction.find("miss") != string::npos)
                    this_yield = yield_root_no_bkg(energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], input_tree);
                else if (Reaction.find("inc") != string::npos)
                    this_yield = yield_root_quadratic_bkg(energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], input_tree);
            }
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
    if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_recon",   strlen("filteredtree_piminus_p_recon")))
        input_tree->Draw(Form("n_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(n_pminus_kin>0.5 && n_pminus_kin<1.3 && beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_thrown",   strlen("filteredtree_piminus_p_thrown")))
        input_tree->Draw(Form("n_pminus_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("n_pminus_truth>0.5 && n_pminus_truth<1.3 && beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    double yield = hist->Integral();
    return yield;
}

double yield_root_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_recon",   strlen("filteredtree_piminus_p_recon")))
        input_tree->Draw(Form("n_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(n_pminus_kin>0.5 && n_pminus_kin<1.3 && beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_thrown",   strlen("filteredtree_piminus_p_thrown")))
        input_tree->Draw(Form("n_pminus_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("n_pminus_truth>0.5 && n_pminus_truth<1.3 && beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2)", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("same");

    return model->GetParameter(0)/hist->GetBinWidth(1) > 10 ? model->GetParameter(0)/hist->GetBinWidth(1) : 0;
}

double yield_root_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_recon",   strlen("filteredtree_piminus_p_recon")))
        input_tree->Draw(Form("n_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(n_pminus_kin>0.5 && n_pminus_kin<1.3 && beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_thrown",   strlen("filteredtree_piminus_p_thrown")))
        input_tree->Draw(Form("n_pminus_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("n_pminus_truth>0.5 && n_pminus_truth<1.3 && beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, 0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("same");

    return model->GetParameter(0)/hist->GetBinWidth(1) > 10 ? model->GetParameter(0)/hist->GetBinWidth(1) : 0;
}

double yield_root_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_recon",   strlen("filteredtree_piminus_p_recon")))
        input_tree->Draw(Form("n_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(n_pminus_kin>0.5 && n_pminus_kin<1.3 && beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_thrown",   strlen("filteredtree_piminus_p_thrown")))
        input_tree->Draw(Form("n_pminus_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("n_pminus_truth>0.5 && n_pminus_truth<1.3 && beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]*x + [4]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, 0, 0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("same");

    return model->GetParameter(0)/hist->GetBinWidth(1) > 10 ? model->GetParameter(0)/hist->GetBinWidth(1) : 0;
}

double yield_root_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 100, 0.4, 1.4);
    if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_recon",   strlen("filteredtree_piminus_p_recon")))
        input_tree->Draw(Form("n_pminus_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(n_pminus_kin>0.5 && n_pminus_kin<1.3 && beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_piminus_p_thrown",   strlen("filteredtree_piminus_p_thrown")))
        input_tree->Draw(Form("n_pminus_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("n_pminus_truth>0.5 && n_pminus_truth<1.3 && beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]*x*x + [4]*x + [5]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, 0, 0, 0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("same");

    return model->GetParameter(0)/hist->GetBinWidth(1) > 10 ? model->GetParameter(0)/hist->GetBinWidth(1) : 0;
}