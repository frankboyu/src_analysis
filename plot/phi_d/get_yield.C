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
double yield_single_gaus_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_single_gaus_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_single_gaus_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_single_gaus_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_double_gaus_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_double_gaus_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_double_gaus_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);
double yield_double_gaus_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree);

TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);

int get_yield(string Reaction)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_%s_all.root",Reaction.c_str());
    string txt_name = Form("output/yield_phi_d_%s.txt",Reaction.c_str());
    string pdf_name = Form("output/plot_phi_d_%s.pdf",Reaction.c_str());

    TFile *input_file = new TFile(input_name.c_str(), "read");
    string input_tree_name;
    if (Reaction.find("recon") != string::npos)
        input_tree_name = "filteredtree_phi_d_recon";
    else if (Reaction.find("thrown") != string::npos)
        input_tree_name = "filteredtree_phi_d_thrown";
    TTree *input_tree = (TTree*) input_file->Get(input_tree_name.c_str());
    FILE *txt_file = fopen(txt_name.c_str(),"w");
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    ifstream matrix_file("output/bin_edges.txt");
    vector<vector<double>> matrix;
    string line;
    while (getline(matrix_file, line))
    {
        istringstream iss(line);
        vector<double> row;
        double value;
        while (iss >> value)
        {
            row.push_back(value);
        }
        matrix.push_back(row);
    }
    matrix_file.close();

    double this_yield;
    for (int i = 0; i < matrix.size(); i++)
    {
        cout << "Theta: " << matrix[i][0] << " deg" << endl;
        cout << "Energy: " << matrix[i][2] << " GeV" << endl;
        if (Reaction.find("sim") != string::npos || Reaction.find("thrown") != string::npos)
        {
            this_yield = yield_single_gaus_no_bkg(matrix[i][2], matrix[i][3], matrix[i][0], matrix[i][1], input_tree);
        }
        else if (Reaction.find("data") != string::npos)
        {
            // this_yield = yield_single_gaus_linear_bkg(matrix[i][2], matrix[i][3], matrix[i][0], matrix[i][1], input_tree);
            this_yield = yield_count(matrix[i][2], matrix[i][3], matrix[i][0], matrix[i][1], input_tree);
        }
        fprintf(txt_file, "%3.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", matrix[i][2], matrix[i][3], matrix[i][0], matrix[i][1], this_yield);
        canvas->Update();
        canvas->Print((pdf_name+"(").c_str());
        canvas->Clear();
    }
    canvas->Print((pdf_name+")").c_str());

    return 0;

}

double yield_count(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    return hist->Integral();
}

double yield_single_gaus_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2)", 0.9, 1.2);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 1.0, 0.01);
    model->SetParLimits(0, 0.1*hist->Integral()*hist->GetBinWidth(1), hist->Integral()*hist->GetBinWidth(1));
    model->SetParLimits(1, 1.01, 1.03);
    model->SetParLimits(2, 0.001, 0.01);
    hist->Fit(model, "R");
    model->Draw("csame");

    return model->GetParameter(0)/hist->GetBinWidth(1);
}

double yield_single_gaus_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]", 0.9, 1.2);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 1.0, 0.01, 0);
    model->SetParLimits(0, 0.1*hist->Integral()*hist->GetBinWidth(1), hist->Integral()*hist->GetBinWidth(1));
    model->SetParLimits(1, 1.01, 1.03);
    model->SetParLimits(2, 0.001, 0.01);
    hist->Fit(model, "R");
    model->SetNpx(1000);
    model->Draw("csame");

    return model->GetParameter(0)/hist->GetBinWidth(1);
}

double yield_single_gaus_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]*x + [4]", 0.9, 1.2);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 1.0, 0.01, 0, 0);
    model->SetParLimits(0, 0.1*hist->Integral()*hist->GetBinWidth(1), hist->Integral()*hist->GetBinWidth(1));
    model->SetParLimits(1, 1.01, 1.03);
    model->SetParLimits(2, 0.001, 0.01);
    hist->Fit(model, "R");
    model->Draw("csame");

    return model->GetParameter(0)/hist->GetBinWidth(1);
}

double yield_single_gaus_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]*x*x + [4]*x + [5]", 0.9, 1.2);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 1.0, 0.01, 0, 0, 0);
    model->SetParLimits(0, 0.1*hist->Integral()*hist->GetBinWidth(1), hist->Integral()*hist->GetBinWidth(1));
    model->SetParLimits(1, 1.01, 1.03);
    model->SetParLimits(2, 0.001, 0.01);
    hist->Fit(model, "R");
    model->Draw("csame");

    return model->GetParameter(0)/hist->GetBinWidth(1);
}

double yield_double_gaus_no_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]*[4])*exp(-0.5*((x-[1])/[4])**2)", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.05);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    model->SetParLimits(3, 0, hist->GetMaximum());
    model->SetParLimits(4, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("csame");

    return (model->GetParameter(0)+model->GetParameter(3))/hist->GetBinWidth(1);
}

double yield_double_gaus_const_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]*[4])*exp(-0.5*((x-[1])/[4])**2) + [5]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.05, 0.0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    model->SetParLimits(3, 0, hist->GetMaximum());
    model->SetParLimits(4, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("csame");

    return (model->GetParameter(0)+model->GetParameter(3))/hist->GetBinWidth(1);
}

double yield_double_gaus_linear_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]*[4])*exp(-0.5*((x-[1])/[4])**2) + [5]*x + [6]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.05, 0.0, 0.0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    model->SetParLimits(3, 0, hist->GetMaximum());
    model->SetParLimits(4, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("csame");

    return (model->GetParameter(0)+model->GetParameter(3))/hist->GetBinWidth(1);
}

double yield_double_gaus_quadratic_bkg(double energy_low, double energy_high, double theta_low, double theta_high, TTree *input_tree)
{
    TH1F *hist = new TH1F(Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), 60, 0.9, 1.2);
    if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_recon",   strlen("filteredtree_phi_d_recon")))
        input_tree->Draw(Form("phi_mass_kin>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && thetaCM_kin>%f && thetaCM_kin<%f)", energy_low, energy_high, theta_low, theta_high));
    else if (!strncmp(input_tree->GetName(), "filteredtree_phi_d_thrown",   strlen("filteredtree_phi_d_thrown")))
        input_tree->Draw(Form("phi_mass_truth>>hist_%.1f_%.1f_%.1f_%.1f", energy_low, energy_high, theta_low, theta_high), Form("beam_energy_truth>%f && beam_energy_truth<%f && thetaCM_truth>%f && thetaCM_truth<%f", energy_low, energy_high, theta_low, theta_high));
    hist->Draw();

    if (hist->Integral() < 50)
        return 0;

    TF1 *model = new TF1("model", "[0]/sqrt(2*TMath::Pi()*[2]*[2])*exp(-0.5*((x-[1])/[2])**2) + [3]/sqrt(2*TMath::Pi()*[4]*[4])*exp(-0.5*((x-[1])/[4])**2) + [5]*x*x + [6]*x + [7]", 0.5, 1.3);
    model->SetParameters(hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.9, 0.05, hist->GetMaximum()*0.9*hist->GetBinWidth(1), 0.05, 0.0, 0.0, 0.0);
    model->SetParLimits(0, 0, hist->GetMaximum());
    model->SetParLimits(1, 0.7, 1.1);
    model->SetParLimits(2, 0.01, 0.20);
    model->SetParLimits(3, 0, hist->GetMaximum());
    model->SetParLimits(4, 0.01, 0.20);
    hist->Fit(model, "R");
    model->Draw("csame");

    return (model->GetParameter(0)+model->GetParameter(3))/hist->GetBinWidth(1);
}