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

void filter_phi_c_2H_kpkmd()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_c_2H_kpkmd_data_ver01.root";
    string InputTree  = "flattree_phi_c_2H_kpkmd";
    string OutputFile = "output/filtered_phi_c_2H_kpkmd_data.root";
    string OutputTree = "filtered_phi_c_2H_kpkmd";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("deuteron_p4",          "TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("nucleon_p4",           "TLorentzVector(0, 0, 0, mass_proton)")
    .Define("phi_p4_meas",          "kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin",           "kp_p4_kin + km_p4_kin")
    .Define("miss_p4_meas",         "beam_p4_meas + deuteron_p4 - phi_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",          "beam_p4_kin + deuteron_p4 - phi_p4_kin - d_p4_kin")
    .Define("missd_p4_meas",        "beam_p4_meas + deuteron_p4 - phi_p4_meas")
    .Define("missd_p4_kin",         "beam_p4_kin + deuteron_p4 - phi_p4_kin")
    .Define("kp_as_pip_p4_meas",    "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    .Define("kp_as_pip_p4_kin",     "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    .Define("km_as_pim_p4_meas",    "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    .Define("km_as_pim_p4_kin",     "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    .Define("beam_p4com_meas",      "boost_lorentz_vector(beam_p4_meas, -(beam_p4_meas + deuteron_p4).BoostVector())")
    .Define("beam_p4com_kin",       "boost_lorentz_vector(beam_p4_kin, -(beam_p4_kin + deuteron_p4).BoostVector())")
    .Define("phi_p4com_meas",       "boost_lorentz_vector(phi_p4_meas, -(beam_p4_meas + deuteron_p4).BoostVector())")
    .Define("phi_p4com_kin",        "boost_lorentz_vector(phi_p4_kin, -(beam_p4_kin + deuteron_p4).BoostVector())")
    .Define("kin_cl",               "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("kp_p_meas",            "kp_p4_meas.P()")
    .Define("kp_p_kin",             "kp_p4_kin.P()")
    .Define("kp_theta_meas",        "kp_p4_meas.Theta()*RadToDeg")
    .Define("kp_theta_kin",         "kp_p4_kin.Theta()*RadToDeg")
    .Define("kp_DeltaT_meas",       "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",        "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    .Define("km_p_meas",            "km_p4_meas.P()")
    .Define("km_p_kin",             "km_p4_kin.P()")
    .Define("km_theta_meas",        "km_p4_meas.Theta()*RadToDeg")
    .Define("km_theta_kin",         "km_p4_kin.Theta()*RadToDeg")
    .Define("km_DeltaT_meas",       "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",        "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    .Define("d_p_meas",             "d_p4_meas.P()")
    .Define("d_p_kin",              "d_p4_kin.P()")
    .Define("d_theta_meas",         "d_p4_meas.Theta()*RadToDeg")
    .Define("d_theta_kin",          "d_p4_kin.Theta()*RadToDeg")
    .Define("d_DeltaT_meas",        "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_meas.T()")
    .Define("d_DeltaT_kin",         "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_kin.T()")
    .Define("missd_m_meas",         "missd_p4_meas.M()")
    .Define("missd_m_kin",          "missd_p4_kin.M()")
    .Define("missd_pMinus_meas",    "missd_p4_meas.Minus()")
    .Define("missd_pMinus_kin",     "missd_p4_kin.Minus()")
    .Define("phi_mass_meas",        "phi_p4_meas.M()")
    .Define("phi_mass_kin",         "phi_p4_kin.M()")
    .Define("phi_theta_meas",       "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",        "phi_p4_kin.Theta()*RadToDeg")
    .Define("miss_energy_meas",     "miss_p4_meas.E()")
    .Define("miss_energy_kin",      "miss_p4_kin.E()")
    .Define("sqrts_meas",           "(beam_p4_meas + deuteron_p4).Mag()")
    .Define("sqrts_kin",            "(beam_p4_kin + deuteron_p4).Mag()")
    .Define("minust_meas",          "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minust_kin",           "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minusu_meas",          "-(deuteron_p4 - phi_p4_meas).Mag2()")
    .Define("minusu_kin",           "-(deuteron_p4 - phi_p4_kin).Mag2()")
    .Define("rho_mass_meas",        "(kp_as_pip_p4_meas + km_as_pim_p4_meas).M()")
    .Define("rho_mass_kin",         "(kp_as_pip_p4_kin + km_as_pim_p4_kin).M()")
    .Define("yphi_meas",            "minust_meas/(2*mass_deuteron*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("yphi_kin",             "minust_kin/(2*mass_deuteron*(beam_p4_kin.E()-phi_p4_kin.E()))")
    .Define("deltaE_meas",          "(pow(sqrts_meas, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrts_meas) - phi_p4com_meas.E()")
    .Define("deltaE_kin",           "(pow(sqrts_kin, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrts_kin) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_no_filtered        = rdf_def;
    auto rdf_dEdx_filtered      = rdf_no_filtered.Filter([](float d_dedx_cdc, double d_p_meas) {return 1000000*d_dedx_cdc > (TMath::Exp(-4.5*d_p_meas+5.0)+2) ;}, {"d_dedx_cdc", "d_p_meas"});
    auto rdf_cl_filtered        = rdf_dEdx_filtered.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered    = rdf_cl_filtered.Filter([](double kp_pidfom, double km_pidfom) {return (kp_pidfom > 0.01) && (km_pidfom > 0.01);}, {"kp_pidfom","km_pidfom"});
    auto rdf_yphi_filtered      = rdf_pidfom_filtered.Filter([](double yphi_meas) {return yphi_meas > 0.4;}, {"yphi_meas"});
    auto rdf_phi_mass_filtered  = rdf_yphi_filtered.Filter([](double phi_mass_meas) {return phi_mass_meas > 1.01 && phi_mass_meas < 1.03;}, {"phi_mass_meas"});

    rdf_phi_mass_filtered.Snapshot(OutputTree, OutputFile);

    // Plot histograms
    cout << "Plotting histograms...\n";
    TFile * histFile = new TFile(OutputFile.c_str(), "update");
    histFile->cd();
    vector<TH1*> hist_list;

    int N_filters = 6;
    RNode rdfs []       = {rdf_no_filtered, rdf_dEdx_filtered,  rdf_cl_filtered,    rdf_pidfom_filtered,    rdf_yphi_filtered,  rdf_phi_mass_filtered};
    string labels []    = {"no_cut",        "dEdx_cut",         "cl_cut",           "pidfom_cut",           "yphi_cut",         "phi_mass_cut"};

    for (int i = 0; i < N_filters; i++)
    {
        auto rdf = rdfs[i];
        string label = labels[i];

        TDirectory * dir = histFile->mkdir(label.c_str());
        dir->cd();

        TH1D hist_DeltaT_kp         = *rdf.Histo1D({("DeltaT_kp_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"kp_DeltaT_meas","accidweight");
        hist_DeltaT_kp.Write();
        TH1D hist_DeltaT_km         = *rdf.Histo1D({("DeltaT_km_"+ label).c_str(), ";#Delta t_{K^{-}} (ns);Counts", 1000, -10.0, 10.0},"km_DeltaT_meas","accidweight");
        hist_DeltaT_km.Write();
        TH1D hist_DeltaT_d          = *rdf.Histo1D({("DeltaT_d_"+ label).c_str(), ";#Delta t_{d} (ns);Counts", 1000, -10.0, 10.0},"d_DeltaT_meas","accidweight");
        hist_DeltaT_d.Write();
        TH1D hist_massKK            = *rdf.Histo1D({("massKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","accidweight");
        hist_massKK.Write();
        TH1D hist_yphi              = *rdf.Histo1D({("yphi_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"yphi_meas","accidweight");
        hist_yphi.Write();
        TH1D hist_minust            = *rdf.Histo1D({("minust_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 2.0},"minust_meas","accidweight");
        hist_minust.Write();
        TH1D hist_Emiss             = *rdf.Histo1D({("Emiss_"+ label).c_str(), ";E_{miss} (GeV);Counts", 400, -2.0, 2.0},"miss_energy_meas","accidweight");
        hist_Emiss.Write();
        TH1D hist_massMissd         = *rdf.Histo1D({("massMissd_"+ label).c_str(), ";m_{missd} (GeV/c^{2});Counts", 200, 1.0, 3.0},"missd_m_meas","accidweight");
        hist_massMissd.Write();
        TH1D hist_pMinusMissd       = *rdf.Histo1D({("pMinusMissd_"+ label).c_str(), ";p_{missd}^{-} (GeV/c);Counts", 100, 1.0, 3.0},"missd_pMinus_meas","accidweight");
        hist_pMinusMissd.Write();
        TH2D hist_kinematics_kp     = *rdf.Histo2D({("kinematics_kp_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 100, 0.0, 180.0},"kp_p_meas","kp_theta_meas","accidweight");
        hist_kinematics_kp.Write();
        TH2D hist_kinematics_km     = *rdf.Histo2D({("kinematics_km_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 100, 0.0, 180.0},"km_p_meas","km_theta_meas","accidweight");
        hist_kinematics_km.Write();
        TH2D hist_kinematics_d      = *rdf.Histo2D({("kinematics_d_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 3.0, 100, 0.0, 180.0},"d_p_meas","d_theta_meas","accidweight");
        hist_kinematics_d.Write();
        TH2D hist_dEdx_kp           = *rdf.Histo2D({("dEdx_kp_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"kp_p_meas","kp_dedx_fdc","accidweight");
        hist_dEdx_kp.Write();
        TH2D hist_dEdx_km           = *rdf.Histo2D({("dEdx_km_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"km_p_meas","km_dedx_fdc","accidweight");
        hist_dEdx_km.Write();
        TH2D hist_dEdx_d            = *rdf.Histo2D({("dEdx_d_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 3.0, 200, 0.0, 2e-5},"d_p_meas","d_dedx_cdc","accidweight");
        hist_dEdx_d.Write();
        TH2D hist_massKK_yphi       = *rdf.Histo2D({("massKK_yphi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);y_{#phi}", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","yphi_meas","accidweight");
        hist_massKK_yphi.Write();
        TH2D hist_massKK_thetaKK    = *rdf.Histo2D({("massKK_thetaKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#theta_{K^{+}K^{-}} (deg)", 400, 0.9, 1.3, 100, 0.0, 10.0},"phi_mass_meas","phi_theta_meas","accidweight");
        hist_massKK_thetaKK.Write();
        TH2D hist_massKK_minust     = *rdf.Histo2D({("massKK_minust_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);-t (GeV^{2}/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","minust_meas","accidweight");
        hist_massKK_minust.Write();
        TH2D hist_massKK_deltaE     = *rdf.Histo2D({("massKK_deltaE_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#DeltaE (GeV)", 400, 0.9, 1.3, 100, -0.5, 0.5},"phi_mass_meas","deltaE_meas","accidweight");
        hist_massKK_deltaE.Write();
        TH2D hist_massKK_massPiPi   = *rdf.Histo2D({("massKK_massPiPi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{#pi^{+}#pi^{-}} (GeV/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","rho_mass_meas","accidweight");
        hist_massKK_massPiPi.Write();
        TH2D hist_thetaKK_yphi      = *rdf.Histo2D({("thetaKK_yphi_"+ label).c_str(), ";#theta_{K^{+}K^{-}} (deg);y_{#phi}", 100, 0.0, 10.0, 200, 0.0, 2.0},"phi_theta_meas","yphi_meas","accidweight");
        hist_thetaKK_yphi.Write();
    }

    histFile->Close();
}