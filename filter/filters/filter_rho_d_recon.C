#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include </work/halld2/home/boyu/src_analysis/filter/filters/const.h>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_rho_d_recon(string Reaction, string input_mode, string output_mode)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/test/selectedtree_rho_d_recon_%s.root",Reaction.c_str());
    string hist_name   = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredhist_rho_d_recon_%s.root",Reaction.c_str());
    string tree_name   = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredtree_rho_d_recon_%s.root",Reaction.c_str());

    // Read input files
    cout << "Reading input files...\n";
    TChain chain("selectedtree_rho_d_recon");
    chain.Add(input_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = RNode(rdf_raw);
    if (input_mode == "one" && Reaction.find("2H") != string::npos)
        rdf_def = rdf_def.Filter("run == 90213");
    else if (input_mode == "one" && Reaction.find("4He") != string::npos)
        rdf_def = rdf_def.Filter("run == 90061");
    else if (input_mode == "one" && Reaction.find("12C") != string::npos)
        rdf_def = rdf_def.Filter("run == 90291");

    auto rdf_input = rdf_def

    .Define("beam_p4com_meas",          "boost_lorentz_vector(beam_p4_meas, -(pip_p4_meas + pim_p4_meas + d_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",           "boost_lorentz_vector(beam_p4_kin, -(pip_p4_kin + pim_p4_kin + d_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",         "boost_lorentz_vector(beam_p4_truth, -(pip_p4_truth + pim_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_meas",         "beam_p4_meas.E()")
    .Define("beam_energy_kin",          "beam_p4_kin.E()")
    .Define("beam_energy_truth",        "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",         "rftime - beam_x4_meas.T()")
    .Define("beam_DeltaT_kin",          "rftime - beam_x4_kin.T()")
    .Define("beam_DeltaT_truth",        "rftime - beam_x4_truth.T()")

    .Define("target_p4_meas",           "pip_p4_meas + pim_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("target_p4_kin",            "pip_p4_kin + pim_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("target_p4_truth",          "pip_p4_truth + pim_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("target_p4rest",            "TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("target_energy_meas",       "target_p4_meas.E()")
    .Define("target_energy_kin",        "target_p4_kin.E()")
    .Define("target_energy_truth",      "target_p4_truth.E()")
    .Define("target_mass_meas",         "target_p4_meas.M()")
    .Define("target_mass_kin",          "target_p4_kin.M()")
    .Define("target_mass_truth",        "target_p4_truth.M()")
    .Define("target_momentum_meas",     "target_p4_meas.P()")
    .Define("target_momentum_kin",      "target_p4_kin.P()")
    .Define("target_momentum_truth",    "target_p4_truth.P()")
    .Define("target_pminus_meas",       "target_p4_meas.Minus()")
    .Define("target_pminus_kin",        "target_p4_kin.Minus()")
    .Define("target_pminus_truth",      "target_p4_truth.Minus()")
    .Define("target_theta_meas",        "target_p4_meas.Theta()*RadToDeg")
    .Define("target_theta_kin",         "target_p4_kin.Theta()*RadToDeg")
    .Define("target_theta_truth",       "target_p4_truth.Theta()*RadToDeg")

    .Define("pip_p4kaon_meas",          "TLorentzVector(pip_p4_meas.Vect(), TMath::Sqrt(pip_p4_meas.P()*pip_p4_meas.P() + mass_kplus*mass_kplus))")
    .Define("pip_p4kaon_kin",           "TLorentzVector(pip_p4_kin.Vect(), TMath::Sqrt(pip_p4_kin.P()*pip_p4_kin.P() + mass_kplus*mass_kplus))")
    .Define("pip_p4kaon_truth",         "TLorentzVector(pip_p4_truth.Vect(), TMath::Sqrt(pip_p4_truth.P()*pip_p4_truth.P() + mass_kplus*mass_kplus))")
    .Define("pip_energy_meas",          "pip_p4_meas.E()")
    .Define("pip_energy_kin",           "pip_p4_kin.E()")
    .Define("pip_energy_truth",         "pip_p4_truth.E()")
    .Define("pip_momentum_meas",        "pip_p4_meas.P()")
    .Define("pip_momentum_kin",         "pip_p4_kin.P()")
    .Define("pip_momentum_truth",       "pip_p4_truth.P()")
    .Define("pip_theta_meas",           "pip_p4_meas.Theta()*RadToDeg")
    .Define("pip_theta_kin",            "pip_p4_kin.Theta()*RadToDeg")
    .Define("pip_theta_truth",          "pip_p4_truth.Theta()*RadToDeg")
    .Define("pip_DeltaT_meas",          "rftime + (pip_x4_meas.Z()-65.0)/29.9792458 - pip_x4_meas.T()")
    .Define("pip_DeltaT_kin",           "rftime + (pip_x4_kin.Z()-65.0)/29.9792458 - pip_x4_kin.T()")
    .Define("pip_DeltaT_truth",         "rftime + (pip_x4_truth.Z()-65.0)/29.9792458 - pip_x4_truth.T()")
    .Define("pip_dedx_fdc_kev_per_cm",  "1000000*pip_dedx_fdc")
    .Define("pip_dedx_cdc_kev_per_cm",  "1000000*pip_dedx_cdc")
    .Define("pip_dedx_st_kev_per_cm",   "1000000*pip_dedx_st")

    .Define("pim_p4kaon_meas",          "TLorentzVector(pim_p4_meas.Vect(), TMath::Sqrt(pim_p4_meas.P()*pim_p4_meas.P() + mass_kminus*mass_kminus))")
    .Define("pim_p4kaon_kin",           "TLorentzVector(pim_p4_kin.Vect(), TMath::Sqrt(pim_p4_kin.P()*pim_p4_kin.P() + mass_kminus*mass_kminus))")
    .Define("pim_p4kaon_truth",         "TLorentzVector(pim_p4_truth.Vect(), TMath::Sqrt(pim_p4_truth.P()*pim_p4_truth.P() + mass_kminus*mass_kminus))")
    .Define("pim_energy_meas",          "pim_p4_meas.E()")
    .Define("pim_energy_kin",           "pim_p4_kin.E()")
    .Define("pim_energy_truth",         "pim_p4_truth.E()")
    .Define("pim_momentum_meas",        "pim_p4_meas.P()")
    .Define("pim_momentum_kin",         "pim_p4_kin.P()")
    .Define("pim_momentum_truth",       "pim_p4_truth.P()")
    .Define("pim_theta_meas",           "pim_p4_meas.Theta()*RadToDeg")
    .Define("pim_theta_kin",            "pim_p4_kin.Theta()*RadToDeg")
    .Define("pim_theta_truth",          "pim_p4_truth.Theta()*RadToDeg")
    .Define("pim_DeltaT_meas",          "rftime + (pim_x4_meas.Z()-65.0)/29.9792458 - pim_x4_meas.T()")
    .Define("pim_DeltaT_kin",           "rftime + (pim_x4_kin.Z()-65.0)/29.9792458 - pim_x4_kin.T()")
    .Define("pim_DeltaT_truth",         "rftime + (pim_x4_truth.Z()-65.0)/29.9792458 - pim_x4_truth.T()")
    .Define("pim_dedx_fdc_kev_per_cm",  "1000000*pim_dedx_fdc")
    .Define("pim_dedx_cdc_kev_per_cm",  "1000000*pim_dedx_cdc")
    .Define("pim_dedx_st_kev_per_cm",   "1000000*pim_dedx_st")

    .Define("rho_p4_meas",              "pip_p4_meas + pim_p4_meas")
    .Define("rho_p4_kin",               "pip_p4_kin + pim_p4_kin")
    .Define("rho_p4_truth",             "pip_p4_truth + pim_p4_truth")
    .Define("rho_p4com_meas",           "boost_lorentz_vector(rho_p4_meas, -(pip_p4_meas + pim_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_p4com_kin",            "boost_lorentz_vector(rho_p4_kin, -(pip_p4_meas + pim_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_p4com_truth",          "boost_lorentz_vector(rho_p4_truth, -(pip_p4_meas + pim_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_energy_meas",          "rho_p4_meas.E()")
    .Define("rho_energy_kin",           "rho_p4_kin.E()")
    .Define("rho_energy_truth",         "rho_p4_truth.E()")
    .Define("rho_momentum_meas",        "rho_p4_meas.P()")
    .Define("rho_momentum_kin",         "rho_p4_kin.P()")
    .Define("rho_momentum_truth",       "rho_p4_truth.P()")
    .Define("rho_mass_meas",            "rho_p4_meas.M()")
    .Define("rho_mass_kin",             "rho_p4_kin.M()")
    .Define("rho_mass_truth",           "rho_p4_truth.M()")
    .Define("rho_theta_meas",           "rho_p4_meas.Theta()*RadToDeg")
    .Define("rho_theta_kin",            "rho_p4_kin.Theta()*RadToDeg")
    .Define("rho_theta_truth",          "rho_p4_truth.Theta()*RadToDeg")

    .Define("d_energy_meas",            "d_p4_meas.E()")
    .Define("d_energy_kin",             "d_p4_kin.E()")
    .Define("d_energy_truth",           "d_p4_truth.E()")
    .Define("d_momentum_meas",          "d_p4_meas.P()")
    .Define("d_momentum_kin",           "d_p4_kin.P()")
    .Define("d_momentum_truth",         "d_p4_truth.P()")
    .Define("d_theta_meas",             "d_p4_meas.Theta()*RadToDeg")
    .Define("d_theta_kin",              "d_p4_kin.Theta()*RadToDeg")
    .Define("d_theta_truth",            "d_p4_truth.Theta()*RadToDeg")
    .Define("d_DeltaT_meas",            "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_meas.T()")
    .Define("d_DeltaT_kin",             "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_kin.T()")
    .Define("d_DeltaT_truth",           "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_truth.T()")
    .Define("d_dedx_fdc_kev_per_cm",    "1000000*d_dedx_fdc")
    .Define("d_dedx_cdc_kev_per_cm",    "1000000*d_dedx_cdc")
    .Define("d_dedx_st_kev_per_cm",     "1000000*d_dedx_st")

    .Define("sqrts_meas",               "(pip_p4_meas + pim_p4_meas + d_p4_meas).Mag()")
    .Define("sqrts_kin",                "(pip_p4_kin + pim_p4_kin + d_p4_kin).Mag()")
    .Define("sqrts_truth",              "(pip_p4_truth + pim_p4_truth + d_p4_truth).Mag()")
    .Define("minust_meas",              "-(beam_p4_meas - rho_p4_meas).Mag2()")
    .Define("minust_kin",               "-(beam_p4_kin - rho_p4_kin).Mag2()")
    .Define("minust_truth",             "-(beam_p4_truth - rho_p4_truth).Mag2()")
    .Define("minusu_meas",              "-(beam_p4_meas - d_p4_meas).Mag2()")
    .Define("minusu_kin",               "-(beam_p4_kin - d_p4_kin).Mag2()")
    .Define("minusu_truth",             "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("y_rho_meas",               "minust_meas/(2*mass_deuteron*(beam_p4_meas.E()-rho_p4_meas.E()))")
    .Define("y_rho_kin",                "minust_kin/(2*mass_deuteron*(beam_p4_kin.E()-rho_p4_kin.E()))")
    .Define("y_rho_truth",              "minust_truth/(2*mass_deuteron*(beam_p4_truth.E()-rho_p4_truth.E()))")
    .Define("phi_mass_meas",            "(pip_p4kaon_meas + pim_p4kaon_meas).M()")
    .Define("phi_mass_kin",             "(pip_p4kaon_kin + pim_p4kaon_kin).M()")
    .Define("phi_mass_truth",           "(pip_p4kaon_truth + pim_p4kaon_truth).M()")
    .Define("energy_balance_meas",      "target_energy_meas - mass_deuteron")
    .Define("energy_balance_kin",       "target_energy_kin - mass_deuteron")
    .Define("energy_balance_truth",     "target_energy_truth - mass_deuteron")

    .Define("kin_cl",                   "TMath::Prob(kin_chisq,kin_ndf)")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut              = rdf_input;
    auto rdf_dEdxCut            = rdf_NoCut.Filter("(d_dedx_cdc_kev_per_cm > (TMath::Exp(-4.5*d_momentum_meas+5.0)+2)) && (d_dedx_cdc_kev_per_cm > 5.0) && (d_dedx_cdc_kev_per_cm < 20.0)");
    auto rdf_KinFitFOMCut       = rdf_dEdxCut.Filter("kin_cl > 0.01");
    auto rdf_PIDFOMCut          = rdf_KinFitFOMCut.Filter("(pip_pidfom > 0.01) && (pim_pidfom > 0.01)");
    auto rdf_RhoMassCut         = rdf_PIDFOMCut.Filter("(rho_mass_kin > 0.5) && (rho_mass_kin < 0.9)");
    auto rdf_output             = rdf_RhoMassCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_RhoMassCut};
    string labels []    = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "RhoMassCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "all")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_rho_d_recon",tree_name);
    }

    // Save histograms
    if (output_mode == "hist" || output_mode == "all")
    {
        cout << "Plotting histograms...\n";
        TFile * hist_file = new TFile(hist_name.c_str(), "RECREATE");
        hist_file->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = hist_file->mkdir(label.c_str());
            dir->cd();

            TH1D hist_beam_DeltaT       = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 400, -20.0, 20.0},"beam_DeltaT_meas");
            hist_beam_DeltaT.Write();

            TH1D hist_target_mass       = *rdf.Histo1D({("target_mass_"+ label).c_str(), ";m_{d} (GeV/c);Counts", 500, 0.0, 5.0},"target_mass_kin","accidweight");
            hist_target_mass.Write();
            TH1D hist_target_momentum   = *rdf.Histo1D({("target_momentum_"+ label).c_str(), ";p_{d} (GeV/c);Counts", 100, 0.0, 5.0},"target_momentum_kin","accidweight");
            hist_target_momentum.Write();
            TH1D hist_target_pminus     = *rdf.Histo1D({("target_pminus_"+ label).c_str(), ";p_{d}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"target_pminus_kin","accidweight");
            hist_target_pminus.Write();
            TH2D hist_target_kinematics = *rdf.Histo2D({("target_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 1.0, 180, 0.0, 180.0},"target_momentum_kin","target_theta_kin","accidweight");
            hist_target_kinematics.Write();

            TH1D hist_pip_DeltaT         = *rdf.Histo1D({("pip_DeltaT_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"pip_DeltaT_meas","accidweight");
            hist_pip_DeltaT.Write();
            TH2D hist_pip_dEdx_fdc       = *rdf.Histo2D({("pip_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pip_momentum_meas","pip_dedx_fdc_kev_per_cm","accidweight");
            hist_pip_dEdx_fdc.Write();
            TH2D hist_pip_dEdx_cdc       = *rdf.Histo2D({("pip_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pip_momentum_meas","pip_dedx_cdc_kev_per_cm","accidweight");
            hist_pip_dEdx_cdc.Write();
            TH2D hist_pip_dEdx_st        = *rdf.Histo2D({("pip_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pip_momentum_meas","pip_dedx_st_kev_per_cm","accidweight");
            hist_pip_dEdx_st.Write();
            TH2D hist_pip_kinematics     = *rdf.Histo2D({("pip_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_meas","pip_theta_meas","accidweight");
            hist_pip_kinematics.Write();

            TH1D hist_pim_DeltaT         = *rdf.Histo1D({("pim_DeltaT_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"pim_DeltaT_meas","accidweight");
            hist_pim_DeltaT.Write();
            TH2D hist_pim_dEdx_fdc       = *rdf.Histo2D({("pim_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pim_momentum_meas","pim_dedx_fdc_kev_per_cm","accidweight");
            hist_pim_dEdx_fdc.Write();
            TH2D hist_pim_dEdx_cdc       = *rdf.Histo2D({("pim_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pim_momentum_meas","pim_dedx_cdc_kev_per_cm","accidweight");
            hist_pim_dEdx_cdc.Write();
            TH2D hist_pim_dEdx_st        = *rdf.Histo2D({("pim_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"pim_momentum_meas","pim_dedx_st_kev_per_cm","accidweight");
            hist_pim_dEdx_st.Write();
            TH2D hist_pim_kinematics     = *rdf.Histo2D({("pim_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_meas","pim_theta_meas","accidweight");
            hist_pim_kinematics.Write();

            TH1D hist_rho_mass          = *rdf.Histo1D({("rho_mass_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.0, 2.0},"rho_mass_kin","accidweight");
            hist_rho_mass.Write();
            TH2D hist_rho_kinematics    = *rdf.Histo2D({("rho_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"rho_momentum_meas","rho_theta_meas","accidweight");
            hist_rho_kinematics.Write();

            TH1D hist_d_DeltaT          = *rdf.Histo1D({("d_DeltaT_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"d_DeltaT_meas","accidweight");
            hist_d_DeltaT.Write();
            TH2D hist_d_dEdx_fdc        = *rdf.Histo2D({("d_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"d_momentum_meas","d_dedx_fdc_kev_per_cm","accidweight");
            hist_d_dEdx_fdc.Write();
            TH2D hist_d_dEdx_cdc        = *rdf.Histo2D({("d_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"d_momentum_meas","d_dedx_cdc_kev_per_cm","accidweight");
            hist_d_dEdx_cdc.Write();
            TH2D hist_d_dEdx_st         = *rdf.Histo2D({("d_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"d_momentum_meas","d_dedx_st_kev_per_cm","accidweight");
            hist_d_dEdx_st.Write();
            TH2D hist_d_kinematics      = *rdf.Histo2D({("d_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_meas","d_theta_meas","accidweight");
            hist_d_kinematics.Write();

            TH1D hist_sqrts             = *rdf.Histo1D({("sqrts_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 500, 2.0, 7.0},"sqrts_kin","accidweight");
            hist_sqrts.Write();
            TH1D hist_minust            = *rdf.Histo1D({("minust_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 400, 0.0, 4.0},"minust_kin","accidweight");
            hist_minust.Write();
            TH1D hist_minusu            = *rdf.Histo1D({("minusu_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 400, 8.0, 12.0},"minusu_kin","accidweight");
            hist_minusu.Write();
            TH1D hist_yrho              = *rdf.Histo1D({("yrho_"+ label).c_str(), ";y_{#rho};Counts", 200, 0.0, 2.0},"y_rho_kin","accidweight");
            hist_yrho.Write();
            TH1D hist_phi_mass          = *rdf.Histo1D({("phi_mass_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c^{2});Counts", 200, 0.0, 2.0},"phi_mass_kin","accidweight");
            hist_rho_mass.Write();
            TH1D hist_energy_balance    = *rdf.Histo1D({("energy_balance_"+ label).c_str(), ";E_{target} - M_{d} (GeV);Counts", 400, -2.0, 2.0},"energy_balance_kin","accidweight");
            hist_energy_balance.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}