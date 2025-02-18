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

void filter_phi_d_recon(string reaction_name, string input_mode, string output_mode)
{
    string input_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_recon_%s.root",reaction_name.c_str());
    string hist_name   = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_%s.root",reaction_name.c_str());
    string tree_name   = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_%s.root",reaction_name.c_str());

    // Determine reaction specific parameters
    double mass_target;
    if (reaction_name.find("2H") != string::npos)
        mass_target = mass_2H;
    else if (reaction_name.find("4He") != string::npos)
        mass_target = mass_4He;
    else if (reaction_name.find("12C") != string::npos)
        mass_target = mass_12C;

    // Read input files
    cout << "Reading input files...\n";
    TChain chain("selectedtree_phi_d_recon");
    chain.Add(input_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = RNode(rdf_raw);
    if (input_mode == "one" && reaction_name.find("2H") != string::npos)
        rdf_def = rdf_def.Filter("run == 90213");
    else if (input_mode == "one" && reaction_name.find("4He") != string::npos)
        rdf_def = rdf_def.Filter("run == 90061");
    else if (input_mode == "one" && reaction_name.find("12C") != string::npos)
        rdf_def = rdf_def.Filter("run == 90291");

    auto rdf_input = rdf_def

    .Define("kinfit_fom",               "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4",                "TLorentzVector(0, 0, 0, mass_target)")

    .Define("beam_p4com_meas",          "boost_lorentz_vector(beam_p4_meas, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",           "boost_lorentz_vector(beam_p4_kin, -(kp_p4_kin + km_p4_kin + d_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",         "boost_lorentz_vector(beam_p4_truth, -(kp_p4_truth + km_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_meas",         "beam_p4_meas.E()")
    .Define("beam_energy_kin",          "beam_p4_kin.E()")
    .Define("beam_energy_truth",        "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",         "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",          "beam_x4_kin.T() - rftime")
    .Define("beam_DeltaT_truth",        "beam_x4_truth.T() - rftime")

    .Define("kp_p4pion_meas",           "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_kin",            "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_truth",          "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_energy_meas",           "kp_p4_meas.E()")
    .Define("kp_energy_kin",            "kp_p4_kin.E()")
    .Define("kp_energy_truth",          "kp_p4_truth.E()")
    .Define("kp_momentum_meas",         "kp_p4_meas.P()")
    .Define("kp_momentum_kin",          "kp_p4_kin.P()")
    .Define("kp_momentum_truth",        "kp_p4_truth.P()")
    .Define("kp_theta_meas",            "kp_p4_meas.Theta()*RadToDeg")
    .Define("kp_theta_kin",             "kp_p4_kin.Theta()*RadToDeg")
    .Define("kp_theta_truth",           "kp_p4_truth.Theta()*RadToDeg")
    .Define("kp_DeltaT_meas",           "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",            "rftime + (kp_x4_kin.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    .Define("kp_DeltaT_truth",          "rftime + (kp_x4_truth.Z()-65.0)/29.9792458 - kp_x4_truth.T()")
    .Define("kp_dedx_fdc_kev_per_cm",   "1000000*kp_dedx_fdc")
    .Define("kp_dedx_cdc_kev_per_cm",   "1000000*kp_dedx_cdc")
    .Define("kp_dedx_st_kev_per_cm",    "1000000*kp_dedx_st")

    .Define("km_p4pion_meas",           "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_kin",            "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_truth",          "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_energy_meas",           "km_p4_meas.E()")
    .Define("km_energy_kin",            "km_p4_kin.E()")
    .Define("km_energy_truth",          "km_p4_truth.E()")
    .Define("km_momentum_meas",         "km_p4_meas.P()")
    .Define("km_momentum_kin",          "km_p4_kin.P()")
    .Define("km_momentum_truth",        "km_p4_truth.P()")
    .Define("km_theta_meas",            "km_p4_meas.Theta()*RadToDeg")
    .Define("km_theta_kin",             "km_p4_kin.Theta()*RadToDeg")
    .Define("km_theta_truth",           "km_p4_truth.Theta()*RadToDeg")
    .Define("km_DeltaT_meas",           "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",            "rftime + (km_x4_kin.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    .Define("km_DeltaT_truth",          "rftime + (km_x4_truth.Z()-65.0)/29.9792458 - km_x4_truth.T()")
    .Define("km_dedx_fdc_kev_per_cm",   "1000000*km_dedx_fdc")
    .Define("km_dedx_cdc_kev_per_cm",   "1000000*km_dedx_cdc")
    .Define("km_dedx_st_kev_per_cm",    "1000000*km_dedx_st")

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

    .Define("phi_p4_meas",              "kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin",               "kp_p4_kin + km_p4_kin")
    .Define("phi_p4_truth",             "kp_p4_truth + km_p4_truth")
    .Define("phi_p4com_meas",           "boost_lorentz_vector(phi_p4_meas, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("phi_p4com_kin",            "boost_lorentz_vector(phi_p4_kin, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())") // inital or final state?
    .Define("phi_p4com_truth",          "boost_lorentz_vector(phi_p4_truth, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("phi_energy_meas",          "phi_p4_meas.E()")
    .Define("phi_energy_kin",           "phi_p4_kin.E()")
    .Define("phi_energy_truth",         "phi_p4_truth.E()")
    .Define("phi_momentum_meas",        "phi_p4_meas.P()")
    .Define("phi_momentum_kin",         "phi_p4_kin.P()")
    .Define("phi_momentum_truth",       "phi_p4_truth.P()")
    .Define("phi_mass_meas",            "phi_p4_meas.M()")
    .Define("phi_mass_kin",             "phi_p4_kin.M()")
    .Define("phi_mass_truth",           "phi_p4_truth.M()")
    .Define("phi_theta_meas",           "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",            "phi_p4_kin.Theta()*RadToDeg")
    .Define("phi_theta_truth",          "phi_p4_truth.Theta()*RadToDeg")

    .Define("target_p4_meas",           "kp_p4_meas + km_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("target_p4_kin",            "kp_p4_kin + km_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("target_p4_truth",          "kp_p4_truth + km_p4_truth + d_p4_truth - beam_p4_truth")
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
    .Define("energy_balance_meas",      "target_energy_meas - mass_deuteron")
    .Define("energy_balance_kin",       "target_energy_kin - mass_deuteron")
    .Define("energy_balance_truth",     "target_energy_truth - mass_deuteron")

    .Define("sqrts_meas",               "(kp_p4_meas + km_p4_meas + d_p4_meas).Mag()")
    .Define("sqrts_kin",                "(kp_p4_kin + km_p4_kin + d_p4_kin).Mag()")
    .Define("sqrts_truth",              "(kp_p4_truth + km_p4_truth + d_p4_truth).Mag()")
    .Define("minust_meas",              "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minust_kin",               "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minust_truth",             "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minusu_meas",              "-(beam_p4_meas - d_p4_meas).Mag2()")
    .Define("minusu_kin",               "-(beam_p4_kin - d_p4_kin).Mag2()")
    .Define("minusu_truth",             "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("y_phi_meas",               "minust_meas/(2*mass_deuteron*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin",                "minust_kin/(2*mass_deuteron*(beam_p4_kin.E()-phi_p4_kin.E()))")
    .Define("y_phi_truth",              "minust_truth/(2*mass_deuteron*(beam_p4_truth.E()-phi_p4_truth.E()))")
    .Define("rho_mass_meas",            "(kp_p4pion_meas + km_p4pion_meas).M()")
    .Define("rho_mass_kin",             "(kp_p4pion_kin + km_p4pion_kin).M()")
    .Define("rho_mass_truth",           "(kp_p4pion_truth + km_p4pion_truth).M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut              = rdf_input;
    auto rdf_dEdxCut            = rdf_NoCut.Filter("(d_dedx_cdc_kev_per_cm > (TMath::Exp(-4.5*d_momentum_meas+5.0)+2)) && (d_dedx_cdc_kev_per_cm > 5.0) && (d_dedx_cdc_kev_per_cm < 20.0)");
    auto rdf_KinFitFOMCut       = rdf_dEdxCut.Filter("kinfit_fom > 0.01");
    auto rdf_PIDFOMCut          = rdf_KinFitFOMCut.Filter("(kp_pidfom > 0.01) && (km_pidfom > 0.01)");
    auto rdf_PhiMassCut         = rdf_PIDFOMCut.Filter("(phi_mass_kin > 1.01) && (phi_mass_kin < 1.03)");
    auto rdf_output             = rdf_PhiMassCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_PhiMassCut};
    string labels []    = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "PhiMassCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "all")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_phi_d_recon",tree_name);
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

            TH1D hist_kp_DeltaT         = *rdf.Histo1D({("kp_DeltaT_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"kp_DeltaT_meas","accidweight");
            hist_kp_DeltaT.Write();
            TH2D hist_kp_dEdx_fdc       = *rdf.Histo2D({("kp_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"kp_momentum_meas","kp_dedx_fdc_kev_per_cm","accidweight");
            hist_kp_dEdx_fdc.Write();
            TH2D hist_kp_dEdx_cdc       = *rdf.Histo2D({("kp_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"kp_momentum_meas","kp_dedx_cdc_kev_per_cm","accidweight");
            hist_kp_dEdx_cdc.Write();
            TH2D hist_kp_dEdx_st        = *rdf.Histo2D({("kp_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"kp_momentum_meas","kp_dedx_st_kev_per_cm","accidweight");
            hist_kp_dEdx_st.Write();
            TH2D hist_kp_kinematics     = *rdf.Histo2D({("kp_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","accidweight");
            hist_kp_kinematics.Write();

            TH1D hist_km_DeltaT         = *rdf.Histo1D({("km_DeltaT_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 1000, -10.0, 10.0},"km_DeltaT_meas","accidweight");
            hist_km_DeltaT.Write();
            TH2D hist_km_dEdx_fdc       = *rdf.Histo2D({("km_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"km_momentum_meas","km_dedx_fdc_kev_per_cm","accidweight");
            hist_km_dEdx_fdc.Write();
            TH2D hist_km_dEdx_cdc       = *rdf.Histo2D({("km_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"km_momentum_meas","km_dedx_cdc_kev_per_cm","accidweight");
            hist_km_dEdx_cdc.Write();
            TH2D hist_km_dEdx_st        = *rdf.Histo2D({("km_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 1000, 0.0, 10.0, 200, 0.0, 40.0},"km_momentum_meas","km_dedx_st_kev_per_cm","accidweight");
            hist_km_dEdx_st.Write();
            TH2D hist_km_kinematics     = *rdf.Histo2D({("km_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","accidweight");
            hist_km_kinematics.Write();

            TH1D hist_phi_mass          = *rdf.Histo1D({("phi_mass_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_kin","accidweight");
            hist_phi_mass.Write();
            TH2D hist_phi_kinematics    = *rdf.Histo2D({("phi_kinematics_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"phi_momentum_meas","phi_theta_meas","accidweight");
            hist_phi_kinematics.Write();

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
            TH1D hist_yphi              = *rdf.Histo1D({("yphi_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_kin","accidweight");
            hist_yphi.Write();
            TH1D hist_rho_mass          = *rdf.Histo1D({("rho_mass_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 200, 0.0, 2.0},"rho_mass_kin","accidweight");
            hist_rho_mass.Write();
            TH1D hist_energy_balance    = *rdf.Histo1D({("energy_balance_"+ label).c_str(), ";E_{target} - M_{d} (GeV);Counts", 400, -2.0, 2.0},"energy_balance_kin","accidweight");
            hist_energy_balance.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}