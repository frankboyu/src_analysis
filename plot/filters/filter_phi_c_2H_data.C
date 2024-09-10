#include <iostream>
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
double mass_phi         = 1.019461;
double RadToDeg         = 180.0 / 3.14159265;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_phi_c_2H_data()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_c_2H_data_ver02/*090213.root";
    string InputTree  = "flattree_phi_c_2H_data";
    string OutputFile = "output/filteredtree_phi_c_2H_data.root";
    string OutputTree = "filteredtree_phi_c_2H_data";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4","TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("missd_p4_meas","beam_p4_meas + target_p4 - kp_p4_meas - km_p4_meas")  // missd_p4_kin already defined in default branches
    .Define("phi_p4_meas","kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin","kp_p4_kin + km_p4_kin")
    .Define("kp_as_pip_p4_meas","TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    .Define("kp_as_pip_p4_kin","TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    .Define("km_as_pim_p4_meas","TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    .Define("km_as_pim_p4_kin","TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    .Define("beam_p4com_meas","boost_lorentz_vector(beam_p4_meas, -(beam_p4_meas + target_p4).BoostVector())")
    .Define("beam_p4com_kin","boost_lorentz_vector(beam_p4_kin, -(beam_p4_kin + target_p4).BoostVector())")
    .Define("phi_p4com_meas","boost_lorentz_vector(phi_p4_meas, -(beam_p4_meas + target_p4).BoostVector())")
    .Define("phi_p4com_kin","boost_lorentz_vector(phi_p4_kin, -(beam_p4_kin + target_p4).BoostVector())")
    .Define("kp_p_meas","kp_p4_meas.P()")
    .Define("kp_p_kin","kp_p4_kin.P()")
    .Define("km_p_meas","km_p4_meas.P()")
    .Define("km_p_kin","km_p4_kin.P()")
    .Define("phi_mass_meas","phi_p4_meas.M()")
    .Define("phi_mass_kin","phi_p4_kin.M()")
    .Define("phi_theta_meas","phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin","phi_p4_kin.Theta()*RadToDeg")
    .Define("missd_mass_meas","missd_p4_meas.M()")
    .Define("missd_mass_kin","missd_p4_kin.M()")
    .Define("sqrt_s_meas", "(beam_p4_meas + target_p4).Mag()")
    .Define("sqrt_s_kin", "(beam_p4_kin + target_p4).Mag()")
    .Define("minus_t_meas", "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minus_u_meas", "-(beam_p4_meas - missd_p4_meas).Mag2()")
    .Define("minus_u_kin", "-(beam_p4_kin - missd_p4_kin).Mag2()")
    .Define("rho_mass_meas","(kp_as_pip_p4_meas + km_as_pim_p4_meas).M()")
    .Define("rho_mass_kin","(kp_as_pip_p4_kin + km_as_pim_p4_kin).M()")
    .Define("y_phi_meas","minus_t_meas/(2*mass_deuteron*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin","minus_t_kin/(2*mass_deuteron*(beam_p4_kin.E()-phi_p4_kin.E()))")
    .Define("deltaE_meas","(pow(sqrt_s_meas, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt_s_meas) - phi_p4com_meas.E()")
    .Define("deltaE_kin","(pow(sqrt_s_kin, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt_s_kin) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_no_filter = rdf_def;
    auto rdf_cl_filtered = rdf_no_filter.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered = rdf_cl_filtered.Filter([](double kp_pidfom, double km_pidfom) {return (kp_pidfom > 0.01) && (km_pidfom > 0.01);}, {"kp_pidfom","km_pidfom"});
    auto rdf_y_phi_filtered = rdf_pidfom_filtered.Filter([](double y_phi_meas) {return y_phi_meas > 0.4;}, {"y_phi_meas"});
    auto rdf_phi_mass_filtered = rdf_y_phi_filtered.Filter([](double phi_mass_meas) {return phi_mass_meas > 1.01 && phi_mass_meas < 1.03;}, {"phi_mass_meas"});

    rdf_phi_mass_filtered.Snapshot("filteredtree_phi_c_2H_data",OutputFile);

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
        TH1D hist_massKK = *rdf.Histo1D({("massKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","accidweight");
        hist_massKK.Write();
        TH1D hist_massmiss = *rdf.Histo1D({("massmiss_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts",200, 1.0, 3.0},"missd_mass_meas","accidweight");
        hist_massmiss.Write();
        TH2D hist_massKK_massmiss = *rdf.Histo2D({("massKK_massmiss_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{miss} (GeV/c^{2})", 400, 0.9, 1.3, 200, 1.0, 3.0},"phi_mass_meas","missd_mass_meas","accidweight");
        hist_massKK_massmiss.Write();
        TH2D hist_massKK_yphi = *rdf.Histo2D({("massKK_yphi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);y_{#phi}", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","y_phi_meas","accidweight");
        hist_massKK_yphi.Write();
        TH2D hist_massKK_thetaKK = *rdf.Histo2D({("massKK_thetaKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#theta_{K^{+}K^{-}} (deg)", 400, 0.9, 1.3, 200, 0.0, 20.0},"phi_mass_meas","phi_theta_meas","accidweight");
        hist_massKK_thetaKK.Write();
        TH2D hist_massKK_minust = *rdf.Histo2D({("massKK_minust_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);-t (GeV^{2}/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","minus_t_meas","accidweight");
        hist_massKK_minust.Write();
        TH2D hist_massKK_deltaE = *rdf.Histo2D({("massKK_deltaE_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#DeltaE (GeV)", 400, 0.9, 1.3, 200, -2.0, 2.0},"phi_mass_meas","deltaE_meas","accidweight");
        hist_massKK_deltaE.Write();
        TH2D hist_massKK_massPiPi = *rdf.Histo2D({("massKK_massPiPi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{#pi^{+}#pi^{-}} (GeV/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","rho_mass_meas","accidweight");
        hist_massKK_massPiPi.Write();
        TH2D hist_thetaKK_yphi = *rdf.Histo2D({("thetaKK_yphi_"+ label).c_str(), ";#theta_{K^{+}K^{-}} (deg);y_{#phi}", 200, 0.0, 20.0, 200, 0.0, 2.0},"phi_theta_meas","y_phi_meas","accidweight");
        hist_thetaKK_yphi.Write();
        TH2D hist_massmiss_yphi = *rdf.Histo2D({("massmiss_yphi_"+ label).c_str(), ";m_{miss} (GeV/c^{2});y_{#phi}", 200, 1.0, 3.0, 200, 0.0, 2.0},"missd_mass_meas","y_phi_meas","accidweight");
        hist_massmiss_yphi.Write();
        TH2D hist_massmiss_deltaE = *rdf.Histo2D({("massmiss_deltaE_"+ label).c_str(), ";m_{miss} (GeV/c^{2});#DeltaE (GeV)", 200, 1.0, 3.0, 200, -2.0, 2.0},"missd_mass_meas","deltaE_meas","accidweight");
        hist_massmiss_deltaE.Write();
    }

    histFile->Close();
}