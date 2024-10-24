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

void filter_phi_c_2H_gen()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/sim/output/test/phi_d_2H_1_initial/root/generator/*.root";
    string InputTree  = "genT";
    string OutputFile = "output/filtered_phi_c_2H_gen.root";
    string OutputTree = "filtered_phi_c_2H_gen";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("beam_p4",   "pBeam")
    .Define("target_p4", "TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("phi_p4",    "pMeson")
    .Define("d_p4",      "pBaryon")
    .Define("gamma_E",   "beam_p4.E()")
    .Define("phi_E",     "phi_p4.E()")
    .Define("d_E",       "d_p4.E()")
    .Define("sqrts",     "(beam_p4 + target_p4).Mag()")
    .Define("minust",    "-(beam_p4 - phi_p4).Mag2()")
    .Define("minusu",    "-(beam_p4 - d_p4).Mag2()")
    .Define("yphi",      "minust/(2*mass_deuteron*(gamma_E-phi_E))")









    // .Define("missd_p4_meas",        "beam_p4_meas + deuteron_p4 - phi_p4_meas")
    // // .Define("missd_p4_kin",         "beam_p4_kin + deuteron_p4 - phi_p4_kin")
    // .Define("kp_as_pip_p4_meas",    "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("kp_as_pip_p4_kin",     "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("km_as_pim_p4_meas",    "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("km_as_pim_p4_kin",     "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("beam_p4com_meas",      "boost_lorentz_vector(beam_p4_meas, -(beam_p4_meas + deuteron_p4).BoostVector())")
    // .Define("beam_p4com_kin",       "boost_lorentz_vector(beam_p4_kin, -(beam_p4_kin + deuteron_p4).BoostVector())")
    // .Define("phi_p4com_meas",       "boost_lorentz_vector(phi_p4_meas, -(beam_p4_meas + deuteron_p4).BoostVector())")
    // .Define("phi_p4com_kin",        "boost_lorentz_vector(phi_p4_kin, -(beam_p4_kin + deuteron_p4).BoostVector())")
    // .Define("kin_cl",               "TMath::Prob(kin_chisq,kin_ndf)")
    // .Define("kp_p_meas",            "kp_p4_meas.P()")
    // .Define("kp_p_kin",             "kp_p4_kin.P()")
    // .Define("kp_theta_meas",        "kp_p4_meas.Theta()*RadToDeg")
    // .Define("kp_theta_kin",         "kp_p4_kin.Theta()*RadToDeg")
    // .Define("kp_DeltaT_meas",       "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    // .Define("kp_DeltaT_kin",        "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    // .Define("km_p_meas",            "km_p4_meas.P()")
    // .Define("km_p_kin",             "km_p4_kin.P()")
    // .Define("km_theta_meas",        "km_p4_meas.Theta()*RadToDeg")
    // .Define("km_theta_kin",         "km_p4_kin.Theta()*RadToDeg")
    // .Define("km_DeltaT_meas",       "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    // .Define("km_DeltaT_kin",        "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    // .Define("missd_m_meas",         "missd_p4_meas.M()")
    // .Define("missd_m_kin",          "missd_p4_kin.M()")
    // .Define("missd_pMinus_meas",    "missd_p4_meas.Minus()")
    // .Define("missd_pMinus_kin",     "missd_p4_kin.Minus()")
    // .Define("phi_mass_meas",        "phi_p4_meas.M()")
    // .Define("phi_mass_kin",         "phi_p4_kin.M()")
    // .Define("phi_theta_meas",       "phi_p4_meas.Theta()*RadToDeg")
    // .Define("phi_theta_kin",        "phi_p4_kin.Theta()*RadToDeg")
    // .Define("rho_mass_meas",        "(kp_as_pip_p4_meas + km_as_pim_p4_meas).M()")
    // .Define("rho_mass_kin",         "(kp_as_pip_p4_kin + km_as_pim_p4_kin).M()")
    // .Define("deltaE_meas",          "(pow(sqrts_meas, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrts_meas) - phi_p4com_meas.E()")
    // .Define("deltaE_kin",           "(pow(sqrts_kin, 2) - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrts_kin) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    rdf_def.Snapshot(OutputTree, OutputFile);

    // Plot histograms
    cout << "Plotting histograms...\n";
    TFile * histFile = new TFile(OutputFile.c_str(), "update");
    histFile->cd();
    vector<TH1*> hist_list;
    string label="gen";

    TH1D hist_yphi                 = *rdf_def.Histo1D({("yphi_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"yphi");
    hist_yphi.Write();
    // TH1D hist_minust            = *rdf.Histo1D({("minust_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 2.0},"minust_meas","accidweight");
    // hist_minust.Write();
    // TH1D hist_massMissd         = *rdf.Histo1D({("massMissd_"+ label).c_str(), ";m_{missd} (GeV/c^{2});Counts", 200, 1.0, 3.0},"missd_m_meas","accidweight");
    // hist_massMissd.Write();
    // TH1D hist_pMinusMissd       = *rdf.Histo1D({("pMinusMissd_"+ label).c_str(), ";p_{missd}^{-} (GeV/c);Counts", 100, 1.0, 3.0},"missd_pMinus_meas","accidweight");
    // hist_pMinusMissd.Write();
    // TH2D hist_massKK_yphi       = *rdf.Histo2D({("massKK_yphi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);y_{#phi}", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","yphi_meas","accidweight");
    // hist_massKK_yphi.Write();
    // TH2D hist_massKK_minust     = *rdf.Histo2D({("massKK_minust_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);-t (GeV^{2}/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","minust_meas","accidweight");
    // hist_massKK_minust.Write();
    // TH2D hist_massKK_deltaE     = *rdf.Histo2D({("massKK_deltaE_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#DeltaE (GeV)", 400, 0.9, 1.3, 100, -0.5, 0.5},"phi_mass_meas","deltaE_meas","accidweight");
    // hist_massKK_deltaE.Write();
    // TH2D hist_massKK_massPiPi   = *rdf.Histo2D({("massKK_massPiPi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{#pi^{+}#pi^{-}} (GeV/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","rho_mass_meas","accidweight");
    // hist_massKK_massPiPi.Write();
    // TH2D hist_thetaKK_yphi      = *rdf.Histo2D({("thetaKK_yphi_"+ label).c_str(), ";#theta_{K^{+}K^{-}} (deg);y_{#phi}", 100, 0.0, 10.0, 200, 0.0, 2.0},"phi_theta_meas","yphi_meas","accidweight");
    // hist_thetaKK_yphi.Write();

    histFile->Close();
}