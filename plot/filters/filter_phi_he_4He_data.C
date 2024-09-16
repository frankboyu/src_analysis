#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstring>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_proton      = 0.938272;
double mass_pion        = 0.13957018;
double mass_phi         = 1.019461;
double mass_unit     = 0.931494102;
double mass_electron = 0.0005109989461;
double mass_helium4  = 4.00260325415*mass_unit - 2*mass_electron;
double RadToDeg         = 180.0 / 3.14159265;

void filter_phi_he_4He_data()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/test/flattree_phi_he_4He_data_*.root";
    string InputTree  = "flattree_phi_he_4He_data";
    string OutputFile = "output/filteredtree_phi_he_4He_data.root";
    string OutputTree = "filteredtree_phi_he_4He_data";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4","TLorentzVector(0, 0, 0, mass_helium4)")
    .Define("phi_p4_meas","kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin","kp_p4_kin + km_p4_kin")
    // .Define("missd_p4_meas","beam_p4_meas + target_p4 - phi_p4_meas")  // missd_p4_kin already defined in default branches
    // .Define("pip_p4_meas","TLorentzVector(kp_p4_meas.X(), kp_p4_meas.Y(), kp_p4_meas.Z(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pip_p4_kin","TLorentzVector(kp_p4_kin.X(), kp_p4_kin.Y(), kp_p4_kin.Z(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_meas","TLorentzVector(km_p4_meas.X(), km_p4_meas.Y(), km_p4_meas.Z(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_kin","TLorentzVector(km_p4_kin.X(), km_p4_kin.Y(), km_p4_kin.Z(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("beam_p4com_meas","beam_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("beam_p4com_kin","beam_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    // .Define("phi_p4com_meas","phi_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("phi_p4com_kin","phi_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    .Define("kp_p_meas","kp_p4_meas.P()")
    .Define("kp_p_kin","kp_p4_kin.P()")
    .Define("km_p_meas","km_p4_meas.P()")
    .Define("km_p_kin","km_p4_kin.P()")
    .Define("he_p_meas","he_p4_meas.P()")
    .Define("he_p_kin","he_p4_kin.P()")
    .Define("phi_mass_meas","phi_p4_meas.M()")
    .Define("phi_mass_kin","phi_p4_kin.M()")
    .Define("phi_theta_meas","phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin","phi_p4_kin.Theta()*RadToDeg")
    // .Define("missd_mass_meas","missd_p4_meas.M()")
    // .Define("missd_mass_kin","missd_p4_kin.M()")
    .Define("sqrt_s_meas", "(beam_p4_meas + target_p4).Mag()")
    .Define("sqrt_s_kin", "(beam_p4_kin + target_p4).Mag()")
    .Define("minus_t_meas", "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(beam_p4_kin - phi_p4_kin).Mag2()")
    // .Define("minus_u_meas", "-(beam_p4_meas - missd_p4_meas).Mag2()")
    // .Define("minus_u_kin", "-(beam_p4_kin - missd_p4_kin).Mag2()")
    // .Define("rho_mass_meas","pip_p4_meas + pim_p4_meas")
    // .Define("rho_mass_kin","pip_p4_kin + pim_p4_kin")
    .Define("y_phi_meas","minus_t_meas/(2*mass_helium4*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin","minus_t_kin/(2*mass_helium4*(beam_p4_kin.E()-phi_p4_kin.E()))")
    // .Define("DeltaE_meas","(s_meas - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_meas)) - phi_p4com_meas.E()")
    // .Define("DeltaE_kin","(s_kin - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_kin)) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_no_filter = rdf_def;
    auto rdf_cl_filtered = rdf_no_filter.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered = rdf_cl_filtered.Filter([](double kp_pidfom, double km_pidfom) {return (kp_pidfom > 0.01) && (km_pidfom > 0.01);}, {"kp_pidfom","km_pidfom"});
    auto rdf_y_phi_filtered = rdf_pidfom_filtered.Filter([](double y_phi_meas) {return y_phi_meas > 0.4;}, {"y_phi_meas"});
    auto rdf_phi_mass_filtered = rdf_y_phi_filtered.Filter([](double phi_mass_meas) {return phi_mass_meas > 1.01 && phi_mass_meas < 1.03;}, {"phi_mass_meas"});

    rdf_phi_mass_filtered.Snapshot("filteredtree_phi_he_4He_data",OutputFile);

    // Plot histograms
    cout << "Plotting histograms...\n";
    TFile * histFile = new TFile(OutputFile.c_str(), "update");
    histFile->cd();
    vector<TH1*> hist_list;

    int N_filters = 5;
    RNode rdfs [] = {rdf_no_filter, rdf_cl_filtered, rdf_pidfom_filtered, rdf_y_phi_filtered, rdf_phi_mass_filtered};
    string labels [] = {"no_cut", "cl_cut", "pidfom_cut", "y_phi_cut", "phi_mass_cut"};


    for (int i = 0; i < N_filters; i++)
    {
        auto rdf = rdfs[i];
        string label = labels[i];

        TDirectory * dir = histFile->mkdir(label.c_str());
        dir->cd();

        TH2D hist_dEdx_kp = *rdf.Histo2D({("dEdx_kp_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"kp_p_meas","kp_dedx_fdc","accidweight");
        hist_dEdx_kp.Write();

        TH2D hist_dEdx_km = *rdf.Histo2D({("dEdx_km_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"km_p_meas","km_dedx_fdc","accidweight");
        hist_dEdx_km.Write();

        TH2D hist_dEdx_he = *rdf.Histo2D({("dEdx_he_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 3.0, 200, 0.0, 2e-5},"he_p_meas","he_dedx_cdc","accidweight");
        hist_dEdx_he.Write();

        TH1D hist_massKK = *rdf.Histo1D({("massKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","accidweight");
        hist_massKK.Write();

        // TH1D hist_massmiss = *rdf.Histo1D({("massmiss_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts",200, 1.0, 3.0},"missd_mass_meas","accidweight");
        // hist_massmiss.Write();

        // TH2D hist_massKK_massmiss = *rdf.Histo2D({("massKK_massmiss_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{miss} (GeV/c^{2})", 1000, 0.0, 2.0, 200, 1.0, 3.0},"phi_mass_meas","missd_mass_meas","accidweight");
        // hist_massKK_massmiss.Write();

        TH2D hist_massKK_yphi = *rdf.Histo2D({("massKK_yphi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);y_{#phi}", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","y_phi_meas","accidweight");
        hist_massKK_yphi.Write();

        TH2D hist_massKK_thetaKK = *rdf.Histo2D({("massKK_thetaKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#theta_{K^{+}K^{-}} (deg)", 400, 0.9, 1.3, 200, 0.0, 20.0},"phi_mass_meas","phi_theta_meas","accidweight");
        hist_massKK_thetaKK.Write();

        TH2D hist_massKK_minust = *rdf.Histo2D({("massKK_minust_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);-t (GeV^{2}/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","minus_t_meas","accidweight");
        hist_massKK_minust.Write();
    }

    histFile->Close();
}