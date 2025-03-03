#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_target = 0.0;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_phi_d_recon(string reaction_name, string input_mode, string output_mode)
{
    string input_name   = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_recon_%s_%s.root",reaction_name.c_str(), input_mode.c_str());
    string hist_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_%s_%s.root",reaction_name.c_str(), input_mode.c_str());
    string tree_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_%s_%s.root",reaction_name.c_str(), input_mode.c_str());

    // Determine reaction specific parameters
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

    auto rdf_input = rdf_raw
    .Define("kinfit_fom",                   "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4",                    "TLorentzVector(0, 0, 0, mass_target)")
    .Define("N2_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_p4com_meas",              "boost_lorentz_vector(beam_p4_meas, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",               "boost_lorentz_vector(beam_p4_kin, -(kp_p4_kin + km_p4_kin + d_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(kp_p4_truth + km_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_meas",             "beam_p4_meas.E()")
    .Define("beam_energy_kin",              "beam_p4_kin.E()")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",             "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",              "beam_x4_kin.T() - rftime")
    .Define("beam_DeltaT_truth",            "beam_x4_truth.T() - rftime")

    .Define("kp_p4pion_meas",               "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_kin",                "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_energy_meas",               "kp_p4_meas.E()")
    .Define("kp_energy_kin",                "kp_p4_kin.E()")
    .Define("kp_energy_truth",              "kp_p4_truth.E()")
    .Define("kp_momentum_meas",             "kp_p4_meas.P()")
    .Define("kp_momentum_kin",              "kp_p4_kin.P()")
    .Define("kp_momentum_truth",            "kp_p4_truth.P()")
    .Define("kp_theta_meas",                "kp_p4_meas.Theta()*RadToDeg")
    .Define("kp_theta_kin",                 "kp_p4_kin.Theta()*RadToDeg")
    .Define("kp_theta_truth",               "kp_p4_truth.Theta()*RadToDeg")
    .Define("kp_in_fdc",                    "accidweight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_in_cdc",                    "accidweight*(kp_dedx_cdc > 0.0 && kp_dedx_fdc == 0.0)")
    .Define("kp_in_fdc_cdc",                "accidweight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc > 0.0)")
    .Define("kp_in_neither",                "accidweight*(kp_dedx_fdc == 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_DeltaT_meas",               "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",                "rftime + (kp_x4_kin.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    .Define("kp_DeltaT_truth",              "rftime + (kp_x4_truth.Z()-65.0)/29.9792458 - kp_x4_truth.T()")
    .Define("kp_dedx_fdc_keV_per_cm",       "kp_dedx_fdc*1e6")
    .Define("kp_dedx_cdc_keV_per_cm",       "kp_dedx_cdc*1e6")
    .Define("kp_dedx_st_keV_per_cm",        "kp_dedx_st*1e6")
    .Define("kp_dedx_tof_keV_per_cm",       "kp_dedx_tof*1e6")

    .Define("km_p4pion_meas",               "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_kin",                "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_energy_meas",               "km_p4_meas.E()")
    .Define("km_energy_kin",                "km_p4_kin.E()")
    .Define("km_energy_truth",              "km_p4_truth.E()")
    .Define("km_momentum_meas",             "km_p4_meas.P()")
    .Define("km_momentum_kin",              "km_p4_kin.P()")
    .Define("km_momentum_truth",            "km_p4_truth.P()")
    .Define("km_theta_meas",                "km_p4_meas.Theta()*RadToDeg")
    .Define("km_theta_kin",                 "km_p4_kin.Theta()*RadToDeg")
    .Define("km_theta_truth",               "km_p4_truth.Theta()*RadToDeg")
    .Define("km_in_fdc",                    "accidweight*(km_dedx_fdc > 0.0 && km_dedx_cdc == 0.0)")
    .Define("km_in_cdc",                    "accidweight*(km_dedx_cdc > 0.0 && km_dedx_fdc == 0.0)")
    .Define("km_in_fdc_cdc",                "accidweight*(km_dedx_fdc > 0.0 && km_dedx_cdc > 0.0)")
    .Define("km_in_neither",                "accidweight*(km_dedx_fdc == 0.0 && km_dedx_cdc == 0.0)")
    .Define("km_DeltaT_meas",               "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",                "rftime + (km_x4_kin.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    .Define("km_DeltaT_truth",              "rftime + (km_x4_truth.Z()-65.0)/29.9792458 - km_x4_truth.T()")
    .Define("km_dedx_fdc_keV_per_cm",       "km_dedx_fdc*1e6")
    .Define("km_dedx_cdc_keV_per_cm",       "km_dedx_cdc*1e6")
    .Define("km_dedx_st_keV_per_cm",        "km_dedx_st*1e6")
    .Define("km_dedx_tof_keV_per_cm",       "km_dedx_tof*1e6")

    .Define("d_energy_meas",                "d_p4_meas.E()")
    .Define("d_energy_kin",                 "d_p4_kin.E()")
    .Define("d_energy_truth",               "d_p4_truth.E()")
    .Define("d_momentum_meas",              "d_p4_meas.P()")
    .Define("d_momentum_kin",               "d_p4_kin.P()")
    .Define("d_momentum_truth",             "d_p4_truth.P()")
    .Define("d_theta_meas",                 "d_p4_meas.Theta()*RadToDeg")
    .Define("d_theta_kin",                  "d_p4_kin.Theta()*RadToDeg")
    .Define("d_theta_truth",                "d_p4_truth.Theta()*RadToDeg")
    .Define("d_DeltaT_meas",                "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_meas.T()")
    .Define("d_DeltaT_kin",                 "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_kin.T()")
    .Define("d_DeltaT_truth",               "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_truth.T()")
    .Define("d_dedx_fdc_keV_per_cm",        "d_dedx_fdc*1e6")
    .Define("d_dedx_cdc_keV_per_cm",        "d_dedx_cdc*1e6")
    .Define("d_dedx_st_keV_per_cm",         "d_dedx_st*1e6")
    .Define("d_dedx_tof_keV_per_cm",        "d_dedx_tof*1e6")

    .Define("phi_p4_meas",                  "kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin",                   "kp_p4_kin + km_p4_kin")
    .Define("phi_p4_truth",                 "kp_p4_truth + km_p4_truth")
    .Define("phi_p4com_meas",               "boost_lorentz_vector(phi_p4_meas, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("phi_p4com_kin",                "boost_lorentz_vector(phi_p4_kin, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())") // inital or final state?
    .Define("phi_p4com_truth",              "boost_lorentz_vector(phi_p4_truth, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("phi_energy_meas",              "phi_p4_meas.E()")
    .Define("phi_energy_kin",               "phi_p4_kin.E()")
    .Define("phi_energy_truth",             "phi_p4_truth.E()")
    .Define("phi_momentum_meas",            "phi_p4_meas.P()")
    .Define("phi_momentum_kin",             "phi_p4_kin.P()")
    .Define("phi_momentum_truth",           "phi_p4_truth.P()")
    .Define("phi_mass_meas",                "phi_p4_meas.M()")
    .Define("phi_mass_kin",                 "phi_p4_kin.M()")
    .Define("phi_mass_truth",               "phi_p4_truth.M()")
    .Define("phi_theta_meas",               "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",                "phi_p4_kin.Theta()*RadToDeg")
    .Define("phi_theta_truth",              "phi_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_meas",               "kp_p4_meas + km_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin",                "kp_p4_kin + km_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("struck_p4_truth",              "kp_p4_truth + km_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("struck_energy_meas",           "struck_p4_meas.E()")
    .Define("struck_energy_kin",            "struck_p4_kin.E()")
    .Define("struck_energy_truth",          "struck_p4_truth.E()")
    .Define("struck_mass_meas",             "struck_p4_meas.M()")
    .Define("struck_mass_kin",              "struck_p4_kin.M()")
    .Define("struck_mass_truth",            "struck_p4_truth.M()")
    .Define("struck_masssquared_meas",      "struck_p4_meas.M2()")
    .Define("struck_masssquared_kin",       "struck_p4_kin.M2()")
    .Define("struck_masssquared_truth",     "struck_p4_truth.M2()")
    .Define("struck_momentum_meas",         "struck_p4_meas.P()")
    .Define("struck_momentum_kin",          "struck_p4_kin.P()")
    .Define("struck_momentum_truth",        "struck_p4_truth.P()")
    .Define("struck_pminus_meas",           "struck_p4_meas.Minus()")
    .Define("struck_pminus_kin",            "struck_p4_kin.Minus()")
    .Define("struck_pminus_truth",          "struck_p4_truth.Minus()")
    .Define("struck_theta_meas",            "struck_p4_meas.Theta()*RadToDeg")
    .Define("struck_theta_kin",             "struck_p4_kin.Theta()*RadToDeg")
    .Define("struck_theta_truth",           "struck_p4_truth.Theta()*RadToDeg")
    .Define("struck_energy_balance_meas",   "struck_energy_meas - mass_2H")
    .Define("struck_energy_balance_kin",    "struck_energy_kin - mass_2H")
    .Define("struck_energy_balance_truth",  "struck_energy_truth - mass_2H")

    .Define("miss_p4_meas",                 "beam_p4_meas + target_p4 - kp_p4_meas - km_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                  "beam_p4_kin + target_p4 - kp_p4_kin - km_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                "beam_p4_truth + target_p4 - kp_p4_truth - km_p4_truth - d_p4_truth")
    .Define("miss_energy_meas",             "miss_p4_meas.E()")
    .Define("miss_energy_kin",              "miss_p4_kin.E()")
    .Define("miss_energy_truth",            "miss_p4_truth.E()")
    .Define("miss_mass_meas",               "miss_p4_meas.M()")
    .Define("miss_mass_kin",                "miss_p4_kin.M()")
    .Define("miss_mass_truth",              "miss_p4_truth.M()")
    .Define("miss_masssquared_meas",        "miss_p4_meas.M2()")
    .Define("miss_masssquared_kin",         "miss_p4_kin.M2()")
    .Define("miss_masssquared_truth",       "miss_p4_truth.M2()")
    .Define("miss_momentum_meas",           "miss_p4_meas.P()")
    .Define("miss_momentum_kin",            "miss_p4_kin.P()")
    .Define("miss_momentum_truth",          "miss_p4_truth.P()")
    .Define("miss_pminus_meas",             "miss_p4_meas.Minus()")
    .Define("miss_pminus_kin",              "miss_p4_kin.Minus()")
    .Define("miss_pminus_truth",            "miss_p4_truth.Minus()")
    .Define("miss_theta_meas",              "miss_p4_meas.Theta()*RadToDeg")
    .Define("miss_theta_kin",               "miss_p4_kin.Theta()*RadToDeg")
    .Define("miss_theta_truth",             "miss_p4_truth.Theta()*RadToDeg")
    .Define("miss_energy_balance_meas",     "miss_energy_meas - (mass_target - mass_2H)")
    .Define("miss_energy_balance_kin",      "miss_energy_kin - (mass_target - mass_2H)")
    .Define("miss_energy_balance_truth",    "miss_energy_truth - (mass_target - mass_2H)")

    .Define("sqrts_meas",                   "(kp_p4_meas + km_p4_meas + d_p4_meas).Mag()")
    .Define("sqrts_kin",                    "(kp_p4_kin + km_p4_kin + d_p4_kin).Mag()")
    .Define("sqrts_truth",                  "(kp_p4_truth + km_p4_truth + d_p4_truth).Mag()")
    .Define("minust_meas",                  "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minust_kin",                   "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minust_truth",                 "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minusu_meas",                  "-(beam_p4_meas - d_p4_meas).Mag2()")
    .Define("minusu_kin",                   "-(beam_p4_kin - d_p4_kin).Mag2()")
    .Define("minusu_truth",                 "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("coplanarity_meas",             "abs(phi_p4_meas.Phi() - d_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin",              "abs(phi_p4_kin.Phi() - d_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth",            "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_meas",                 "beam_p4com_meas.Vect().Angle(phi_p4com_meas.Vect())*RadToDeg")
    .Define("thetaCM_kin",                  "beam_p4com_kin.Vect().Angle(phi_p4com_kin.Vect())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(phi_p4com_truth.Vect())*RadToDeg")
    .Define("y_phi_meas",                   "minust_meas/(2*mass_2H*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin",                    "minust_kin/(2*mass_2H*(beam_p4_kin.E()-phi_p4_kin.E()))")
    .Define("y_phi_truth",                  "minust_truth/(2*mass_2H*(beam_p4_truth.E()-phi_p4_truth.E()))")
    .Define("rho_mass_meas",                "(kp_p4pion_meas + km_p4pion_meas).M()")
    .Define("rho_mass_kin",                 "(kp_p4pion_kin + km_p4pion_kin).M()")
    .Define("rho_mass_truth",               "(kp_p4pion_truth + km_p4pion_truth).M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut              = rdf_input;
    auto rdf_dEdxCut            = rdf_NoCut.Filter("(d_momentum_meas > 0.35) && (d_momentum_meas < 1.30) && (d_dedx_cdc_keV_per_cm > (TMath::Exp(-29.68353898*d_momentum_meas+13.50623694)+17.88279645*d_momentum_meas*d_momentum_meas-42.15473796*d_momentum_meas+28.83200736))");
    auto rdf_KinFitFOMCut       = rdf_dEdxCut.Filter("kinfit_fom > 0.01");
    auto rdf_PIDFOMCut          = rdf_KinFitFOMCut.Filter("(kp_pidfom > 0.01) && (km_pidfom > 0.01)");
    auto rdf_PhiMassCut         = rdf_PIDFOMCut.Filter("(phi_mass_kin > 1.01) && (phi_mass_kin < 1.03)");
    auto rdf_output             = rdf_PhiMassCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_PhiMassCut};
    string labels []    = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "PhiMassCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_phi_d_recon",tree_name);
    }

    // Save histograms
    if (output_mode == "hist" || output_mode == "both")
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

            TH1D hist_kinfit_fom                    = *rdf.Histo1D({("kinfit_fom_"+ label).c_str(), ";KinFit FOM;Counts", 100, 0.0, 1.0},"kinfit_fom","accidweight");
            hist_kinfit_fom.Write();

            TH1D hist_beam_DeltaT_kin               = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 400, -20.0, 20.0},"beam_DeltaT_kin");
            hist_beam_DeltaT_kin.Write();

            TH1D hist_kp_pidfom                     = *rdf.Histo1D({("kp_pidfom_"+ label).c_str(), ";kp_pidfom;Counts", 100, 0.0, 1.0},"kp_pidfom","accidweight");
            hist_kp_pidfom.Write();
            TH2D hist_kp_dEdx_cdc_kin               = *rdf.Histo2D({("kp_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"kp_momentum_kin","kp_dedx_cdc_keV_per_cm","accidweight");
            hist_kp_dEdx_cdc_kin.Write();
            TH2D hist_kp_dEdx_fdc_kin               = *rdf.Histo2D({("kp_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"kp_momentum_kin","kp_dedx_fdc_keV_per_cm","accidweight");
            hist_kp_dEdx_fdc_kin.Write();
            TH2D hist_kp_dEdx_tof_kin               = *rdf.Histo2D({("kp_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"kp_momentum_kin","kp_dedx_tof_keV_per_cm","accidweight");
            hist_kp_dEdx_tof_kin.Write();
            TH2D hist_kp_dEdx_st_kin                = *rdf.Histo2D({("kp_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"kp_momentum_kin","kp_dedx_st_keV_per_cm","accidweight");
            hist_kp_dEdx_st_kin.Write();
            TH1D hist_kp_DeltaT_kin                 = *rdf.Histo1D({("kp_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{-}} (ns);Counts", 400, -20.0, 20.0},"kp_DeltaT_kin","accidweight");
            hist_kp_DeltaT_kin.Write();
            TH2D hist_kp_kinematics_kin             = *rdf.Histo2D({("kp_kinematics_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","accidweight");
            hist_kp_kinematics_kin.Write();
            TH2D hist_kp_kinematics_fdc_kin         = *rdf.Histo2D({("kp_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc");
            hist_kp_kinematics_fdc_kin.Write();
            TH2D hist_kp_kinematics_fdc_cdc_kin     = *rdf.Histo2D({("kp_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc_cdc");
            hist_kp_kinematics_fdc_cdc_kin.Write();
            TH2D hist_kp_kinematics_cdc_kin         = *rdf.Histo2D({("kp_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_cdc");
            hist_kp_kinematics_cdc_kin.Write();
            TH2D hist_kp_kinematics_neither_kin     = *rdf.Histo2D({("kp_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_neither");
            hist_kp_kinematics_neither_kin.Write();

            TH1D hist_km_pidfom                     = *rdf.Histo1D({("km_pidfom_"+ label).c_str(), ";km_pidfom;Counts", 100, 0.0, 1.0},"km_pidfom","accidweight");
            hist_km_pidfom.Write();
            TH2D hist_km_dEdx_cdc_kin               = *rdf.Histo2D({("km_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"km_momentum_kin","km_dedx_cdc_keV_per_cm","accidweight");
            hist_km_dEdx_cdc_kin.Write();
            TH2D hist_km_dEdx_fdc_kin               = *rdf.Histo2D({("km_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"km_momentum_kin","km_dedx_fdc_keV_per_cm","accidweight");
            hist_km_dEdx_fdc_kin.Write();
            TH2D hist_km_dEdx_tof_kin               = *rdf.Histo2D({("km_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"km_momentum_kin","km_dedx_tof_keV_per_cm","accidweight");
            hist_km_dEdx_tof_kin.Write();
            TH2D hist_km_dEdx_st_kin                = *rdf.Histo2D({("km_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"km_momentum_kin","km_dedx_st_keV_per_cm","accidweight");
            hist_km_dEdx_st_kin.Write();
            TH1D hist_km_DeltaT_kin                 = *rdf.Histo1D({("km_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{-}} (ns);Counts", 400, -20.0, 20.0},"km_DeltaT_kin","accidweight");
            hist_km_DeltaT_kin.Write();
            TH2D hist_km_kinematics_kin             = *rdf.Histo2D({("km_kinematics_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","accidweight");
            hist_km_kinematics_kin.Write();
            TH2D hist_km_kinematics_fdc_kin         = *rdf.Histo2D({("km_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc");
            hist_km_kinematics_fdc_kin.Write();
            TH2D hist_km_kinematics_fdc_cdc_kin     = *rdf.Histo2D({("km_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc_cdc");
            hist_km_kinematics_fdc_cdc_kin.Write();
            TH2D hist_km_kinematics_cdc_kin         = *rdf.Histo2D({("km_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_cdc");
            hist_km_kinematics_cdc_kin.Write();
            TH2D hist_km_kinematics_neither_kin     = *rdf.Histo2D({("km_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_neither");
            hist_km_kinematics_neither_kin.Write();

            TH2D hist_d_dEdx_cdc_meas               = *rdf.Histo2D({("d_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_meas","d_dedx_cdc_keV_per_cm","accidweight");
            hist_d_dEdx_cdc_meas.Write();
            TH2D hist_d_dEdx_cdc_kin                = *rdf.Histo2D({("d_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_cdc_keV_per_cm","accidweight");
            hist_d_dEdx_cdc_kin.Write();
            TH2D hist_d_dEdx_tof_kin                = *rdf.Histo2D({("d_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_tof_keV_per_cm","accidweight");
            hist_d_dEdx_tof_kin.Write();
            TH2D hist_d_dEdx_st_kin                 = *rdf.Histo2D({("d_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_st_keV_per_cm","accidweight");
            hist_d_dEdx_st_kin.Write();
            TH1D hist_d_DeltaT_kin                  = *rdf.Histo1D({("d_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{-}} (ns);Counts", 400, -20.0, 20.0},"d_DeltaT_kin","accidweight");
            hist_d_DeltaT_kin.Write();
            TH2D hist_d_kinematics_kin              = *rdf.Histo2D({("d_kinematics_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_kin","d_theta_kin","accidweight");
            hist_d_kinematics_kin.Write();

            TH1D hist_phi_mass_kin                  = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_kin","accidweight");
            hist_phi_mass_kin.Write();
            TH2D hist_phi_kinematics_kin            = *rdf.Histo2D({("phi_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"phi_momentum_kin","phi_theta_kin","accidweight");
            hist_phi_kinematics_kin.Write();

            TH1D hist_sqrts_kin                     = *rdf.Histo1D({("sqrts_kin_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_kin","accidweight");
            hist_sqrts_kin.Write();
            TH1D hist_minust_kin                    = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minust_kin","accidweight");
            hist_minust_kin.Write();
            TH1D hist_minusu_kin                    = *rdf.Histo1D({("minusu_kin_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 300, 0.0, 30.0},"minusu_kin","accidweight");
            hist_minusu_kin.Write();
            TH1D hist_coplanarity_kin               = *rdf.Histo1D({("coplanarity_kin_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_kin","accidweight");
            hist_coplanarity_kin.Write();
            TH1D hist_thetaCM_kin                   = *rdf.Histo1D({("thetaCM_kin_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_kin","accidweight");
            hist_thetaCM_kin.Write();
            TH2D hist_minust_thetaCM_kin            = *rdf.Histo2D({("minust_thetaCM_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_kin","thetaCM_kin","accidweight");
            hist_minust_thetaCM_kin.Write();
            TH1D hist_rho_mass_kin                  = *rdf.Histo1D({("rho_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_mass_kin","accidweight");
            hist_rho_mass_kin.Write();
            TH1D hist_yphi_kin                      = *rdf.Histo1D({("yphi_kin_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_kin","accidweight");
            hist_yphi_kin.Write();

            TH1D hist_struck_mass_kin               = *rdf.Histo1D({("struck_mass_kin_"+ label).c_str(), ";m_{struck} (GeV/c^{2});Counts", 400, -4.0, 4.0},"struck_mass_kin","accidweight");
            hist_struck_mass_kin.Write();
            TH1D hist_struck_masssquared_kin        = *rdf.Histo1D({("struck_masssquared_kin_"+ label).c_str(), ";m_{struck}^{2} (GeV^{2}/c^{4});Counts", 400, 0.0, 4.0},"struck_masssquared_kin","accidweight");
            hist_struck_masssquared_kin.Write();
            TH1D hist_struck_momentum_kin           = *rdf.Histo1D({("struck_momentum_kin_"+ label).c_str(), ";P_{struck} (GeV/c);Counts", 100, 0.0, 1.0},"struck_momentum_kin","accidweight");
            hist_struck_momentum_kin.Write();
            TH1D hist_struck_pminus_kin             = *rdf.Histo1D({("struck_pminus_kin_"+ label).c_str(), ";P_{struck}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"struck_pminus_kin","accidweight");
            hist_struck_pminus_kin.Write();
            TH1D hist_struck_energy_balance_kin     = *rdf.Histo1D({("struck_energy_balance_kin_"+ label).c_str(), ";E_{struck} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"struck_energy_balance_kin","accidweight");
            hist_struck_energy_balance_kin.Write();
            TH2D hist_struck_momentum_pminus_kin    = *rdf.Histo2D({("struck_momentum_pminus_kin_"+ label).c_str(), ";P_{struck} (GeV/c);P_{struck}^{-} (GeV/c)", 100, 0.0, 1.0, 200, 1.0, 3.0},"struck_momentum_kin","struck_pminus_kin","accidweight");
            hist_struck_momentum_pminus_kin.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}