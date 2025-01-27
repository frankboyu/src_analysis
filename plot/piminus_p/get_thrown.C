using namespace std;

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


int get_gen(string Reaction)
{
    string input_name;
    if (Reaction == "2H_flat")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/generator/genOut_gen_MF_*.root";
    else if (Reaction == "2H_model")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/generator/genOut_gen_MF_*.root";
    else if (Reaction == "4He_flat")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/generator/genOut_gen_MF_*.root";
    else if (Reaction == "4He_model")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/generator/genOut_gen_MF_*.root";
    else if (Reaction == "12C_flat")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/generator/genOut_gen_MF_*.root";
    else if (Reaction == "12C_model")
        input_name  = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/generator/genOut_gen_MF_*.root";

    string txt_name = Form("output/yield_piminus_p_gen_%s.txt",Reaction.c_str());
    string pdf_name = Form("output/plot_piminus_p_gen_%s.txt",Reaction.c_str());

    TChain *input_chain = new TChain("genT");
    input_chain->Add(input_name.c_str());
    FILE *txt_file = fopen(txt_name.c_str(),"w");
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);

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
    double yield[energy_bins][theta_bins];
    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            yield[i][j] = 0;
        }
    }
    TH1D* histograms[energy_bins][theta_bins];
    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            histograms[i][j] = new TH1D(Form("hist_%d_%d", i, j), Form("Energy bin %f, Theta bin %f", energy_edges[i], theta_edges[j]), 80, 0.5, 1.3);
        }
    }
    TLorentzVector *pBeam = new TLorentzVector();
    TLorentzVector *pMeson = new TLorentzVector();
    TLorentzVector *pBaryon = new TLorentzVector();
    TLorentzVector pBeam_boost;
    TLorentzVector pMeson_boost;
    double energy, theta, pminus;

    input_chain->SetBranchAddress("pBeam", &pBeam);
    input_chain->SetBranchAddress("pMeson", &pMeson);
    input_chain->SetBranchAddress("pBaryon", &pBaryon);

    for(int i = 0; i < input_chain->GetEntries(); i++)
    {
        input_chain->GetEntry(i);
        double energy = pBeam->E();
        pBeam_boost = *pBeam;
        pMeson_boost = *pMeson;
        TVector3 boost = -(*pMeson+*pBaryon).BoostVector();
        pBeam_boost.Boost(boost);
        pMeson_boost.Boost(boost);
        theta = pBeam_boost.Angle(pMeson_boost.Vect()) * 180 / TMath::Pi();
        pminus = (*pMeson + *pBaryon - *pBeam).Minus();

        for (int j = 0; j < energy_bins; j++)
        {
            if (energy >= energy_edges[j] && energy < energy_edges[j+1])
            {
                for (int k = 0; k < theta_bins; k++)
                {
                    if (theta >= theta_edges[k] && theta < theta_edges[k+1])
                    {
                        if ((*pMeson + *pBaryon - *pBeam).P() < 0.2 && (*pBaryon).P() > 0.2)
                            {
                            yield[j][k] += 1;
                            histograms[j][k]->Fill(pminus);
                            }
                        break;
                    }
                }
                break;
            }
        }
    }


    for (int i = 0; i < energy_bins; i++)
    {
        for (int j = 0; j < theta_bins; j++)
        {
            fprintf(txt_file, "%3.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", energy_edges[i], energy_edges[i+1], theta_edges[j], theta_edges[j+1], yield[i][j]);
            histograms[i][j]->Draw();
            canvas->Update();
            canvas->Print((pdf_name+"(").c_str());
            canvas->Clear();
        }
    }
    canvas->Print((pdf_name+")").c_str());

    return 0;

}