#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <cstring>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_deuteron    = 1.875612;
double mass_proton      = 0.938272;
double mass_pion        = 0.139570;
double RadToDeg         = 180.0 / 3.14159265;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_piminus_p_recon()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/piminus_p_12C/flattree_piminus_p_recon/*090262.root";
    string InputTree  = "flattree_piminus_p_recon";
    string OutputFile = "output/filtered_piminus_p_12C_data.root";
    string OutputTree = "filtered_piminus_p_recon";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("miss_p4_meas","pim_p4_meas + p_p4_meas - beam_p4_meas")
    .Define("miss_p4_kin","pim_p4_kin + p_p4_kin - beam_p4_kin")
    .Define("p_as_pip_p4_meas","TLorentzVector(p_p4_meas.Vect(), TMath::Sqrt(p_p4_meas.P()*p_p4_meas.P() + mass_pion*mass_pion))")
    .Define("p_as_pip_p4_kin","TLorentzVector(p_p4_kin.Vect(), TMath::Sqrt(p_p4_kin.P()*p_p4_kin.P() + mass_pion*mass_pion))")
    .Define("beam_p4com_meas","boost_lorentz_vector(beam_p4_meas, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("beam_p4com_kin","boost_lorentz_vector(beam_p4_kin, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("pim_p4com_meas","boost_lorentz_vector(pim_p4_meas, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("pim_p4com_kin","boost_lorentz_vector(pim_p4_kin, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("pim_p_meas","pim_p4_meas.P()")
    .Define("pim_p_kin","pim_p4_kin.P()")
    .Define("p_p_meas","p_p4_meas.P()")
    .Define("p_p_kin","p_p4_kin.P()")
    .Define("miss_mass_meas","miss_p4_meas.M()")
    .Define("miss_mass_kin","miss_p4_kin.M()")
    .Define("miss_p_meas","miss_p4_meas.P()")
    .Define("miss_p_kin","miss_p4_kin.P()")
    .Define("miss_pminus_meas","miss_p4_meas.Minus()")
    .Define("miss_pminus_kin","miss_p4_kin.Minus()")
    .Define("sqrt_s_meas", "(pim_p4_meas + p_p4_meas).Mag()")
    .Define("sqrt_s_kin", "(pim_p4_meas + p_p4_meas).Mag()")
    .Define("minus_t_meas", "-(pim_p4_meas - beam_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(pim_p4_kin - beam_p4_kin).Mag2()")
    .Define("minus_u_meas", "-(p_p4_meas - beam_p4_meas).Mag2()")
    .Define("minus_u_kin", "-(p_p4_kin - beam_p4_kin).Mag2()")
    .Define("rho_mass_meas","(pim_p4_meas + p_as_pip_p4_meas).M()")
    .Define("rho_mass_kin","(pim_p4_kin + p_as_pip_p4_kin).M()")
    .Define("coplanarity_meas","abs(pim_p4_meas.Phi() - p_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin","abs(pim_p4_kin.Phi() - p_p4_kin.Phi())*RadToDeg")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_no_filter          = rdf_def;
    auto rdf_cl_filtered        = rdf_no_filter.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered    = rdf_cl_filtered.Filter([](double pim_pidfom, double p_pidfom) {return (pim_pidfom > 0.01) && (p_pidfom > 0.01);}, {"pim_pidfom","p_pidfom"});
    auto rdf_miss_p_filtered   = rdf_pidfom_filtered.Filter([](double miss_p_kin) {return miss_p_kin < 2.0;}, {"miss_p_kin"});
    auto rdf_miss_pminus_filtered  = rdf_miss_p_filtered.Filter([](double miss_pminus_kin) {return miss_pminus_kin > 0.5 && miss_pminus_kin < 1.3;}, {"miss_pminus_kin"});
    auto rdf_all_filters = rdf_miss_pminus_filtered;

    // rdf_all_filters.Snapshot("filtered_piminus_p_recon",OutputFile);

    // Plot histograms
    cout << "Plotting histograms...\n";
    TFile * histFile = new TFile(OutputFile.c_str(), "recreate");
    histFile->cd();
    vector<TH1*> hist_list;

    int N_filters = 5;
    RNode rdfs [] = {rdf_no_filter, rdf_cl_filtered, rdf_pidfom_filtered, rdf_miss_p_filtered, rdf_miss_pminus_filtered};
    string labels [] = {"no_cut", "cl_cut", "pidfom_cut", "miss_p_cut", "miss_pminus_cut"};


    for (int i = 0; i < N_filters; i++)
    {
        auto rdf = rdfs[i];
        string label = labels[i];

        TDirectory * dir = histFile->mkdir(label.c_str());
        dir->cd();

        TH1D hist_miss_pminus = *rdf.Histo1D({("miss_pminus_"+ label).c_str(), ";P_{miss}^- (GeV/c);Counts", 100, 0.4, 1.4},"miss_pminus_kin","accidweight");
        hist_miss_pminus.Write();
        TH1D hist_miss_mass = *rdf.Histo1D({("miss_mass_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 100, 0.0, 4.0},"miss_mass_kin","accidweight");
        hist_miss_mass.Write();
        TH2D hist_miss_pminus_miss_p = *rdf.Histo2D({("miss_pminus_miss_p_"+ label).c_str(), ";P_{miss}^- (GeV/c);P_{miss} (GeV/c)", 100, 0.4, 1.4, 100, 0.0, 4.0},"miss_pminus_kin","miss_p_kin","accidweight");
        hist_miss_pminus_miss_p.Write();
        TH1D hist_miss_p = *rdf.Histo1D({("miss_p_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 100, 0.0, 4.0},"miss_p_kin","accidweight");
        hist_miss_p.Write();
        TH1D hist_minus_t = *rdf.Histo1D({("minus_t_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minus_t_kin","accidweight");
        hist_minus_t.Write();
        TH1D hist_sqrt_s = *rdf.Histo1D({("sqrt_s_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrt_s_kin","accidweight");
        hist_sqrt_s.Write();
        TH1D hist_rho_mass = *rdf.Histo1D({("rho_mass_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.4, 1.4},"rho_mass_kin","accidweight");
        hist_rho_mass.Write();
        TH1D hist_coplanarity = *rdf.Histo1D({("coplanarity_"+ label).c_str(), ";Coplanarity (deg);Counts", 100, -180.0, 180.0},"coplanarity_kin","accidweight");
        hist_coplanarity.Write();
        TH1D hist_accidweight = *rdf.Histo1D({("accidweight_"+ label).c_str(), ";accidweight;Counts", 20, -0.5, 1.5},"accidweight");
        hist_accidweight.Write();
        TH2D hist_pim_dEdx_cdc = *rdf.Histo2D({("pim_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_cdc","accidweight");
        hist_pim_dEdx_cdc.Write();
        TH2D hist_pim_dEdx_fdc = *rdf.Histo2D({("pim_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_fdc","accidweight");
        hist_pim_dEdx_fdc.Write();
        TH2D hist_p_dEdx_cdc = *rdf.Histo2D({("p_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_cdc","accidweight");
        hist_p_dEdx_cdc.Write();
        TH2D hist_p_dEdx_fdc = *rdf.Histo2D({("p_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_fdc","accidweight");
        hist_p_dEdx_fdc.Write();
    }

    histFile->Close();
}