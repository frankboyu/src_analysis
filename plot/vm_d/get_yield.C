#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;

int get_yield(string channel, string reaction, string observable)
{
    // Read tree from input root file
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_%s_%s.root", channel.c_str(), reaction.c_str());
    cout << "Input tree file: " << input_treefile_name << endl;
    TFile *input_treefile = new TFile(input_treefile_name.c_str(), "read");
    string input_tree_name;
    if (reaction.find("recon") != string::npos)
        input_tree_name = Form("filteredtree_%s_recon", channel.c_str());
    else if (reaction.find("thrown") != string::npos)
        input_tree_name = Form("filteredtree_%s_thrown", channel.c_str());
    cout << "Input tree: " << input_tree_name << endl;
    TTree *input_tree = (TTree*) input_treefile->Get(input_tree_name.c_str());

    // Read bin edges from input text file
    string target;
    if (reaction.find("2H") != string::npos)
        target = "2H";
    else if (reaction.find("4He") != string::npos)
        target = "4He";
    else if (reaction.find("12C") != string::npos)
        target = "12C";
    string input_txt_name = Form("output/bins_%s_%s_%s.txt", channel.c_str(), target.c_str(), observable.c_str());
    cout << "Input text file: " << input_txt_name << endl;
    ifstream input_txt(input_txt_name.c_str());
    vector<vector<double>> bins;
    string this_line;
    while (getline(input_txt, this_line))
    {
        istringstream iss(this_line);
        vector<double> this_row;
        double this_value;
        while (iss >> this_value)
            this_row.push_back(this_value);
        bins.push_back(this_row);
    }
    input_txt.close();

    // Prepare output files
    string output_textfile_name = Form("output/yield_%s_%s_%s.txt", channel.c_str(), reaction.c_str(), observable.c_str());
    string output_pdffile_name = Form("output/plot_%s_%s_%s.pdf", channel.c_str(), reaction.c_str(), observable.c_str());
    cout << "Output text file: " << output_textfile_name << endl;
    cout << "Output pdf file: " << output_pdffile_name << endl;
    FILE *output_textfile = fopen(output_textfile_name.c_str(),"w");

    // Declare variables used for the yield calculation
    string yield_observable;
    int hist_bins;
    double hist_min, hist_max;
    if (channel.find("phi") != string::npos)
    {
        yield_observable = "phi_mass";
        hist_bins = 60;
        hist_min = 0.9;
        hist_max = 1.2;
    }
    else if (channel.find("rho") != string::npos)
    {
        yield_observable = "rho_mass";
        hist_bins = 60;
        hist_min = 0.2;
        hist_max = 2.0;
    }
    else if (channel.find("omega") != string::npos)
    {
        yield_observable = "omega_mass";
        hist_bins = 60;
        hist_min = 0.6;
        hist_max = 1.2;
    }
    if (reaction.find("recon") != string::npos)
        yield_observable += "_kin";
    else if (reaction.find("thrown") != string::npos)
        yield_observable += "_truth";
    TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);

    // Loop over bins and calculate yield
    double this_yield;
    TH1F *this_hist;
    string this_hist_name, this_hist_cut;
    for (int i = 0; i < bins.size(); i++)
    {
        if (observable == "dsdt")
        {
            cout << "Photon energy: " << bins[i][0] << "-" << bins[i][1] << " GeV" << ", ";
            cout << "-t: " << bins[i][2] << "-" << bins[i][3] << " GeV^2" << endl;
            this_hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]);
            if (reaction.find("recon") != string::npos)
                this_hist_cut = Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && minust_kin>%f && minust_kin<%f)", bins[i][0], bins[i][1], bins[i][2], bins[i][3]);
            else if (reaction.find("thrown") != string::npos)
                this_hist_cut = Form("beam_energy_truth>%f && beam_energy_truth<%f && minust_truth>%f && minust_truth<%f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]);
        }
        else if (observable == "Wcostheta")
        {
            cout << "Photon energy: " << bins[i][0] << "-" << bins[i][1] << " GeV" << ", ";
            cout << "-t: " << bins[i][2] << "-" << bins[i][3] << " GeV^2" << ", ";
            cout << "cos(theta_H): " << bins[i][4] << "-" << bins[i][5] << endl;
            this_hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            if (reaction.find("recon") != string::npos)
                this_hist_cut = Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && minust_kin>%f && minust_kin<%f && decay_costheta_helicity_kin>%f && decay_costheta_helicity_kin<%f)", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            else if (reaction.find("thrown") != string::npos)
                this_hist_cut = Form("beam_energy_truth>%f && beam_energy_truth<%f && minust_truth>%f && minust_truth<%f && decay_costheta_helicity_truth>%f && decay_costheta_helicity_truth<%f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
        }
        else if (observable == "Wphi")
        {
            cout << "Photon energy: " << bins[i][0] << "-" << bins[i][1] << " GeV" << ", ";
            cout << "-t: " << bins[i][2] << "-" << bins[i][3] << " GeV^2" << ", ";
            cout << "phi: " << bins[i][4] << "-" << bins[i][5] << " deg" << endl;
            this_hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            if (reaction.find("recon") != string::npos)
                this_hist_cut = Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && minust_kin>%f && minust_kin<%f && decay_phi_helicity_kin>%f && decay_phi_helicity_kin<%f)", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            else if (reaction.find("thrown") != string::npos)
                this_hist_cut = Form("beam_energy_truth>%f && beam_energy_truth<%f && minust_truth>%f && minust_truth<%f && decay_phi_helicity_truth>%f && decay_phi_helicity_truth<%f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
        }
        else if (observable == "WPhi")
        {
            cout << "Photon energy: " << bins[i][0] << "-" << bins[i][1] << " GeV" << ", ";
            cout << "-t: " << bins[i][2] << "-" << bins[i][3] << " GeV^2" << ", ";
            cout << "Phi: " << bins[i][4] << "-" << bins[i][5] << " deg" << endl;
            this_hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            if (reaction.find("recon") != string::npos)
                this_hist_cut = Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && minust_kin>%f && minust_kin<%f && polarization_phi_com_kin>%f && polarization_phi_com_kin<%f)", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            else if (reaction.find("thrown") != string::npos)
                this_hist_cut = Form("beam_energy_truth>%f && beam_energy_truth<%f && minust_truth>%f && minust_truth<%f && polarization_phi_com_truth>%f && polarization_phi_com_truth<%f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
        }
        else if (observable == "Wpsi")
        {
            cout << "Photon energy: " << bins[i][0] << "-" << bins[i][1] << " GeV" << ", ";
            cout << "-t: " << bins[i][2] << "-" << bins[i][3] << " GeV^2" << ", ";
            cout << "psi: " << bins[i][4] << "-" << bins[i][5] << " deg" << endl;
            this_hist_name = Form("hist_%.1f_%.1f_%.1f_%.1f_%.1f_%.1f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            if (reaction.find("recon") != string::npos)
                this_hist_cut = Form("accidweight*(beam_energy_kin>%f && beam_energy_kin<%f && minust_kin>%f && minust_kin<%f && psi_helicity_kin>%f && psi_helicity_kin<%f)", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
            else if (reaction.find("thrown") != string::npos)
                this_hist_cut = Form("beam_energy_truth>%f && beam_energy_truth<%f && minust_truth>%f && minust_truth<%f && psi_helicity_truth>%f && psi_helicity_truth<%f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
        }

        this_hist = new TH1F(this_hist_name.c_str(), this_hist_name.c_str(), hist_bins, hist_min, hist_max);
        input_tree->Draw(Form("%s>>%s", yield_observable.c_str(), this_hist_name.c_str()), this_hist_cut.c_str());  // draw the yield observable into the histogram with the cut
        this_yield = this_hist->Integral();
        if (observable == "dsdt")
            fprintf(output_textfile, "%6.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3], this_yield);
        else if (observable == "Wcostheta" || observable == "Wphi" || observable == "WPhi" || observable == "Wpsi")
            fprintf(output_textfile, "%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%6.1f\t%f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5], this_yield);
        this_hist->Draw();
        canvas->Update();
        canvas->Print((output_pdffile_name+"(").c_str());
        canvas->Clear();
        delete this_hist;
    }
    canvas->Print((output_pdffile_name+")").c_str());  // close the pdf file

    return 0;

}