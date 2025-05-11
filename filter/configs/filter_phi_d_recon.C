#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;

void filter_phi_d_recon(string reaction, string output_mode)
{
    // Read input files
    cout << "Reading input files...\n";
    string input_root_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_recon_%s.root",reaction.c_str());
    TChain chain("selectedtree_phi_d_recon");
    chain.Add(input_root_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    if      (reaction.find("2H")   != string::npos)
        mass_target = mass_2H;
    else if (reaction.find("4He")  != string::npos)
        mass_target = mass_4He;
    else if (reaction.find("12C")  != string::npos)
        mass_target = mass_12C;

    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    auto rdf_input = rdf_def
    .Define("beam_p4com_meas",              "boost_lorentz_vector(beam_p4_meas, -(kp_p4_meas + km_p4_meas + d_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",               "boost_lorentz_vector(beam_p4_kin, -(kp_p4_kin + km_p4_kin + d_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(kp_p4_truth + km_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_meas",             "beam_p4_meas.E()")
    .Define("beam_energy_kin",              "beam_p4_kin.E()")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",             "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",              "beam_x4_kin.T() - rftime")

    .Define("kp_p4pion_meas",               "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_kin",                "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4pion_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4helicity_meas",           "boost_lorentz_vector(kp_p4_meas, -(kp_p4_meas + km_p4_meas).BoostVector())")
    .Define("kp_p4helicity_kin",            "boost_lorentz_vector(kp_p4_kin, -(kp_p4_kin + km_p4_kin).BoostVector())")
    .Define("kp_p4helicity_truth",          "boost_lorentz_vector(kp_p4_truth, -(kp_p4_truth + km_p4_truth).BoostVector())")
    .Define("kp_energy_meas",               "kp_p4_meas.E()")
    .Define("kp_energy_kin",                "kp_p4_kin.E()")
    .Define("kp_energy_truth",              "kp_p4_truth.E()")
    .Define("kp_momentum_meas",             "kp_p4_meas.P()")
    .Define("kp_momentum_kin",              "kp_p4_kin.P()")
    .Define("kp_momentum_truth",            "kp_p4_truth.P()")
    .Define("kp_theta_meas",                "kp_p4_meas.Theta()*RadToDeg")
    .Define("kp_theta_kin",                 "kp_p4_kin.Theta()*RadToDeg")
    .Define("kp_theta_truth",               "kp_p4_truth.Theta()*RadToDeg")
    .Define("kp_DeltaT_meas",               "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",                "rftime + (kp_x4_kin.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    .Define("kp_in_fdc",                    "accidweight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_in_cdc",                    "accidweight*(kp_dedx_cdc > 0.0 && kp_dedx_fdc == 0.0)")
    .Define("kp_in_fdc_cdc",                "accidweight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc > 0.0)")
    .Define("kp_in_neither",                "accidweight*(kp_dedx_fdc == 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_dedx_fdc_keV_per_cm",       "kp_dedx_fdc*1e6")
    .Define("kp_dedx_cdc_keV_per_cm",       "kp_dedx_cdc*1e6")
    .Define("kp_dedx_st_keV_per_cm",        "kp_dedx_st*1e6")
    .Define("kp_dedx_tof_keV_per_cm",       "kp_dedx_tof*1e6")

    .Define("km_p4pion_meas",               "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_kin",                "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_piminus*mass_piminus))")
    .Define("km_p4pion_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_p4helicity_meas",           "boost_lorentz_vector(km_p4_meas, -(kp_p4_meas + km_p4_meas).BoostVector())")
    .Define("km_p4helicity_kin",            "boost_lorentz_vector(km_p4_kin, -(kp_p4_kin + km_p4_kin).BoostVector())")
    .Define("km_p4helicity_truth",          "boost_lorentz_vector(km_p4_truth, -(kp_p4_truth + km_p4_truth).BoostVector())")
    .Define("km_energy_meas",               "km_p4_meas.E()")
    .Define("km_energy_kin",                "km_p4_kin.E()")
    .Define("km_energy_truth",              "km_p4_truth.E()")
    .Define("km_momentum_meas",             "km_p4_meas.P()")
    .Define("km_momentum_kin",              "km_p4_kin.P()")
    .Define("km_momentum_truth",            "km_p4_truth.P()")
    .Define("km_theta_meas",                "km_p4_meas.Theta()*RadToDeg")
    .Define("km_theta_kin",                 "km_p4_kin.Theta()*RadToDeg")
    .Define("km_theta_truth",               "km_p4_truth.Theta()*RadToDeg")
    .Define("km_DeltaT_meas",               "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",                "rftime + (km_x4_kin.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    .Define("km_in_fdc",                    "accidweight*(km_dedx_fdc > 0.0 && km_dedx_cdc == 0.0)")
    .Define("km_in_cdc",                    "accidweight*(km_dedx_cdc > 0.0 && km_dedx_fdc == 0.0)")
    .Define("km_in_fdc_cdc",                "accidweight*(km_dedx_fdc > 0.0 && km_dedx_cdc > 0.0)")
    .Define("km_in_neither",                "accidweight*(km_dedx_fdc == 0.0 && km_dedx_cdc == 0.0)")
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
    .Define("d_dedx_fdc_keV_per_cm",        "d_dedx_fdc*1e6")
    .Define("d_dedx_cdc_keV_per_cm",        "d_dedx_cdc*1e6")
    .Define("d_dedx_st_keV_per_cm",         "d_dedx_st*1e6")
    .Define("d_dedx_tof_keV_per_cm",        "d_dedx_tof*1e6")

    .Define("phi_p4_meas",                  "kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin",                   "kp_p4_kin + km_p4_kin")
    .Define("phi_p4_truth",                 "kp_p4_truth + km_p4_truth")
    .Define("phi_p4com_meas",               "boost_lorentz_vector(phi_p4_meas, -(phi_p4_meas + d_p4_meas).BoostVector())")
    .Define("phi_p4com_kin",                "boost_lorentz_vector(phi_p4_kin, -(phi_p4_kin + d_p4_kin).BoostVector())")
    .Define("phi_p4com_truth",              "boost_lorentz_vector(phi_p4_truth, -(phi_p4_truth + d_p4_truth).BoostVector())")
    .Define("phi_energy_meas",              "phi_p4_meas.E()")
    .Define("phi_energy_kin",               "phi_p4_kin.E()")
    .Define("phi_energy_truth",             "phi_p4_truth.E()")
    .Define("phi_momentum_meas",            "phi_p4_meas.P()")
    .Define("phi_momentum_kin",             "phi_p4_kin.P()")
    .Define("phi_momentum_truth",           "phi_p4_truth.P()")
    .Define("phi_mass_meas",                "phi_p4_meas.M()")
    .Define("phi_mass_kin",                 "phi_p4_kin.M()")
    .Define("phi_mass_truth",               "phi_p4_truth.M()")
    .Define("phi_proxymass_meas",           "TMath::Sqrt((kp_p4_meas.Minus() + km_p4_meas.Minus())*(2*beam_p4_meas.E() + mass_4He - d_p4_meas.Plus() - (mass_2H*mass_2H + (phi_p4_meas + d_p4_meas).Perp2())/(mass_4He - (phi_p4_meas + d_p4_meas).Minus())) - (kp_p4_meas.Px() + km_p4_meas.Px())*(kp_p4_meas.Px() + km_p4_meas.Px()) - (kp_p4_meas.Py() + km_p4_meas.Py())*(kp_p4_meas.Py() + km_p4_meas.Py()))")
    .Define("phi_proxymass_kin",            "TMath::Sqrt((kp_p4_kin.Minus() + km_p4_kin.Minus())*(2*beam_p4_kin.E() + mass_4He - d_p4_kin.Plus() - (mass_2H*mass_2H + (phi_p4_kin + d_p4_kin).Perp2())/(mass_4He - (phi_p4_kin + d_p4_kin).Minus())) - (kp_p4_kin.Px() + km_p4_kin.Px())*(kp_p4_kin.Px() + km_p4_kin.Px()) - (kp_p4_kin.Py() + km_p4_kin.Py())*(kp_p4_kin.Py() + km_p4_kin.Py()))")
    .Define("phi_proxymass_truth",          "TMath::Sqrt((kp_p4_truth.Minus() + km_p4_truth.Minus())*(2*beam_p4_truth.E() + mass_4He - d_p4_truth.Plus() - (mass_2H*mass_2H + (phi_p4_truth + d_p4_truth).Perp2())/(mass_4He - (phi_p4_truth + d_p4_truth).Minus())) - (kp_p4_truth.Px() + km_p4_truth.Px())*(kp_p4_truth.Px() + km_p4_truth.Px()) - (kp_p4_truth.Py() + km_p4_truth.Py())*(kp_p4_truth.Py() + km_p4_truth.Py()))")
    .Define("phi_theta_meas",               "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",                "phi_p4_kin.Theta()*RadToDeg")
    .Define("phi_theta_truth",              "phi_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_meas",               "phi_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin",                "phi_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("struck_p4_truth",              "phi_p4_truth + d_p4_truth - beam_p4_truth")
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

    .Define("miss_p4_meas",                 "beam_p4_meas + TLorentzVector(0, 0, 0, mass_target) - phi_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                  "beam_p4_kin + TLorentzVector(0, 0, 0, mass_target) - phi_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target) - phi_p4_truth - d_p4_truth")
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

    .Define("sqrts_meas",                   "(phi_p4_meas + d_p4_meas).Mag()")
    .Define("sqrts_kin",                    "(phi_p4_kin + d_p4_kin).Mag()")
    .Define("sqrts_truth",                  "(phi_p4_truth + d_p4_truth).Mag()")
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
    .Define("pi_helicity_meas",             "kp_p4helicity_meas.Vect().Unit()")
    .Define("pi_helicity_kin",              "kp_p4helicity_kin.Vect().Unit()")
    .Define("pi_helicity_truth",            "kp_p4helicity_truth.Vect().Unit()")
    .Define("z_helicity_meas",              "phi_p4com_meas.Vect().Unit()")
    .Define("z_helicity_kin",               "phi_p4com_kin.Vect().Unit()")
    .Define("z_helicity_truth",             "phi_p4com_truth.Vect().Unit()")
    .Define("y_helicity_meas",              "beam_p4com_meas.Vect().Cross(phi_p4com_meas.Vect()).Unit()")
    .Define("y_helicity_kin",               "beam_p4com_kin.Vect().Cross(phi_p4com_kin.Vect()).Unit()")
    .Define("y_helicity_truth",             "beam_p4com_truth.Vect().Cross(phi_p4com_truth.Vect()).Unit()")
    .Define("x_helicity_meas",              "y_helicity_meas.Cross(z_helicity_meas).Unit()")
    .Define("x_helicity_kin",               "y_helicity_kin.Cross(z_helicity_kin).Unit()")
    .Define("x_helicity_truth",             "y_helicity_truth.Cross(z_helicity_truth).Unit()")
    .Define("costheta_helicity_meas",       "pi_helicity_meas.Dot(z_helicity_meas)")
    .Define("costheta_helicity_kin",        "pi_helicity_kin.Dot(z_helicity_kin)")
    .Define("costheta_helicity_truth",      "pi_helicity_truth.Dot(z_helicity_truth)")
    .Define("phi_helicity_meas",            "TMath::ATan2(-x_helicity_meas.Dot(pi_helicity_meas.Cross(z_helicity_meas)), y_helicity_meas.Dot(pi_helicity_meas.Cross(z_helicity_meas)))*RadToDeg")
    .Define("phi_helicity_kin",             "TMath::ATan2(-x_helicity_kin.Dot(pi_helicity_kin.Cross(z_helicity_kin)), y_helicity_kin.Dot(pi_helicity_kin.Cross(z_helicity_kin)))*RadToDeg")
    .Define("phi_helicity_truth",           "TMath::ATan2(-x_helicity_truth.Dot(pi_helicity_truth.Cross(z_helicity_truth)), y_helicity_truth.Dot(pi_helicity_truth.Cross(z_helicity_truth)))*RadToDeg")
    .Define("kinfit_fom_kin",               "TMath::Prob(kin_chisq,kin_ndf)")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_dEdxCut        = rdf_NoCut.Filter("(d_dedx_cdc_keV_per_cm > (TMath::Exp(-29.68353898*d_momentum_meas+13.50623694)+17.88279645*d_momentum_meas*d_momentum_meas-42.15473796*d_momentum_meas+28.83200736)) && (d_dedx_cdc_keV_per_cm < (TMath::Exp(-26.69276323*d_momentum_meas+15.92466317)+17.1164272*d_momentum_meas*d_momentum_meas-48.7542903*d_momentum_meas+40.25692313))");
    auto rdf_KinFitFOMCut   = rdf_dEdxCut.Filter("kinfit_fom_kin > 0.01");
    auto rdf_PIDFOMCut      = rdf_KinFitFOMCut.Filter("(kp_pidfom > 0.01) && (km_pidfom > 0.01)");
    auto rdf_MissPCut       = rdf_PIDFOMCut.Filter("(abs(struck_energy_balance_kin) < 1.0)");
    auto rdf_output         = rdf_MissPCut;
    RNode rdfs []           = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_MissPCut};
    string labels []        = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "MissPCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_tree_file = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_%s.root",reaction.c_str());
        rdf_output.Snapshot("filteredtree_phi_d_recon",output_tree_file.c_str());
    }

    // Save histograms
    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_hist_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_%s.root",reaction.c_str());
        TFile * output_hist_file = new TFile(output_hist_name.c_str(), "RECREATE");
        output_hist_file->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = output_hist_file->mkdir(label.c_str());
            dir->cd();

            TH1D hist_beam_energy_kin               = *rdf.Histo1D({("beam_energy_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_kin");
            hist_beam_energy_kin.Write();
            TH1D hist_beam_DeltaT_kin               = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 400, -20.0, 20.0},"beam_DeltaT_kin");
            hist_beam_DeltaT_kin.Write();

            TH2D hist_kp_kinematics_kin             = *rdf.Histo2D({("kp_kinematics_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","accidweight");
            hist_kp_kinematics_kin.Write();
            TH2D hist_kp_kinematics_fdc_kin         = *rdf.Histo2D({("kp_kinematics_fdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc");
            hist_kp_kinematics_fdc_kin.Write();
            TH2D hist_kp_kinematics_fdc_cdc_kin     = *rdf.Histo2D({("kp_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc_cdc");
            hist_kp_kinematics_fdc_cdc_kin.Write();
            TH2D hist_kp_kinematics_cdc_kin         = *rdf.Histo2D({("kp_kinematics_cdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_cdc");
            hist_kp_kinematics_cdc_kin.Write();
            TH2D hist_kp_kinematics_neither_kin     = *rdf.Histo2D({("kp_kinematics_neither_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_neither");
            hist_kp_kinematics_neither_kin.Write();
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
            TH1D hist_kp_DeltaT_kin                 = *rdf.Histo1D({("kp_DeltaT_kin_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 400, -20.0, 20.0},"kp_DeltaT_kin","accidweight");
            hist_kp_DeltaT_kin.Write();

            TH2D hist_km_kinematics_kin             = *rdf.Histo2D({("km_kinematics_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","accidweight");
            hist_km_kinematics_kin.Write();
            TH2D hist_km_kinematics_fdc_kin         = *rdf.Histo2D({("km_kinematics_fdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc");
            hist_km_kinematics_fdc_kin.Write();
            TH2D hist_km_kinematics_fdc_cdc_kin     = *rdf.Histo2D({("km_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc_cdc");
            hist_km_kinematics_fdc_cdc_kin.Write();
            TH2D hist_km_kinematics_cdc_kin         = *rdf.Histo2D({("km_kinematics_cdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_cdc");
            hist_km_kinematics_cdc_kin.Write();
            TH2D hist_km_kinematics_neither_kin     = *rdf.Histo2D({("km_kinematics_neither_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_neither");
            hist_km_kinematics_neither_kin.Write();
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
            TH1D hist_km_DeltaT_kin                 = *rdf.Histo1D({("km_DeltaT_kin_"+ label).c_str(), ";#Delta t_{K^{-}} (ns);Counts", 400, -20.0, 20.0},"km_DeltaT_kin","accidweight");
            hist_km_DeltaT_kin.Write();

            TH2D hist_d_kinematics_kin              = *rdf.Histo2D({("d_kinematics_kin_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_kin","d_theta_kin","accidweight");
            hist_d_kinematics_kin.Write();
            TH2D hist_d_dEdx_cdc_meas               = *rdf.Histo2D({("d_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_meas","d_dedx_cdc_keV_per_cm","accidweight");
            hist_d_dEdx_cdc_meas.Write();
            TH2D hist_d_dEdx_cdc_kin                = *rdf.Histo2D({("d_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_cdc_keV_per_cm","accidweight");
            hist_d_dEdx_cdc_kin.Write();
            TH2D hist_d_dEdx_tof_kin                = *rdf.Histo2D({("d_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_tof_keV_per_cm","accidweight");
            hist_d_dEdx_tof_kin.Write();
            TH2D hist_d_dEdx_st_kin                 = *rdf.Histo2D({("d_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"d_momentum_kin","d_dedx_st_keV_per_cm","accidweight");
            hist_d_dEdx_st_kin.Write();
            TH1D hist_d_DeltaT_kin                  = *rdf.Histo1D({("d_DeltaT_kin_"+ label).c_str(), ";#Delta t_{d} (ns);Counts", 400, -20.0, 20.0},"d_DeltaT_kin","accidweight");
            hist_d_DeltaT_kin.Write();

            TH1D hist_phi_mass_kin                  = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_kin","accidweight");
            hist_phi_mass_kin.Write();
            TH1D hist_phi_mass_meas                 = *rdf.Histo1D({("phi_mass_meas_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","accidweight");
            hist_phi_mass_meas.Write();
            TH1D hist_phi_proxymass_kin             = *rdf.Histo1D({("phi_proxymass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}X} (GeV/c^{2});Counts", 400, 0.9, 1.3},"phi_proxymass_kin","accidweight");
            hist_phi_proxymass_kin.Write();
            TH1D hist_phi_proxymass_meas            = *rdf.Histo1D({("phi_proxymass_meas_"+ label).c_str(), ";m_{K^{+}K^{-}X} (GeV/c^{2});Counts", 400, 0.9, 1.3},"phi_proxymass_meas","accidweight");
            hist_phi_proxymass_meas.Write();
            TH2D hist_phi_kinematics_kin            = *rdf.Histo2D({("phi_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"phi_momentum_kin","phi_theta_kin","accidweight");
            hist_phi_kinematics_kin.Write();

            TH1D hist_miss_mass_kin                 = *rdf.Histo1D({("miss_mass_kin_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 400, -4.0, 4.0},"miss_mass_kin","accidweight");
            hist_miss_mass_kin.Write();
            TH1D hist_miss_mass_meas                = *rdf.Histo1D({("miss_mass_meas_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 400, -4.0, 4.0},"miss_mass_meas","accidweight");
            hist_miss_mass_meas.Write();
            TH1D hist_miss_masssquared_kin          = *rdf.Histo1D({("miss_masssquared_kin_"+ label).c_str(), ";m_{miss}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"miss_masssquared_kin","accidweight");
            hist_miss_masssquared_kin.Write();
            TH1D hist_miss_masssquared_meas         = *rdf.Histo1D({("miss_masssquared_meas_"+ label).c_str(), ";m_{miss}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"miss_masssquared_meas","accidweight");
            hist_miss_masssquared_meas.Write();
            TH1D hist_miss_momentum_kin             = *rdf.Histo1D({("miss_momentum_kin_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_kin","accidweight");
            hist_miss_momentum_kin.Write();
            TH1D hist_miss_momentum_meas            = *rdf.Histo1D({("miss_momentum_meas_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_meas","accidweight");
            hist_miss_momentum_meas.Write();
            TH1D hist_miss_pminus_kin               = *rdf.Histo1D({("miss_pminus_kin_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"miss_pminus_kin","accidweight");
            hist_miss_pminus_kin.Write();
            TH1D hist_miss_pminus_meas              = *rdf.Histo1D({("miss_pminus_meas_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"miss_pminus_meas","accidweight");
            hist_miss_pminus_meas.Write();
            TH1D hist_miss_energy_balance_kin       = *rdf.Histo1D({("miss_energy_balance_kin_"+ label).c_str(), ";E_{miss} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"miss_energy_balance_kin","accidweight");
            hist_miss_energy_balance_kin.Write();
            TH1D hist_miss_energy_balance_meas      = *rdf.Histo1D({("miss_energy_balance_meas_"+ label).c_str(), ";E_{miss} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"miss_energy_balance_meas","accidweight");
            hist_miss_energy_balance_meas.Write();
            TH2D hist_miss_momentum_pminus_kin      = *rdf.Histo2D({("miss_momentum_pminus_kin_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"miss_momentum_kin","miss_pminus_kin","accidweight");
            hist_miss_momentum_pminus_kin.Write();
            TH2D hist_miss_momentum_pminus_meas     = *rdf.Histo2D({("miss_momentum_pminus_meas_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"miss_momentum_meas","miss_pminus_meas","accidweight");
            hist_miss_momentum_pminus_meas.Write();
            TH2D hist_miss_momentum_energy_kin      = *rdf.Histo2D({("miss_momentum_energy_kin_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"miss_momentum_kin","miss_energy_kin","accidweight");
            hist_miss_momentum_energy_kin.Write();
            TH2D hist_miss_momentum_energy_meas     = *rdf.Histo2D({("miss_momentum_energy_meas_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"miss_momentum_meas","miss_energy_meas","accidweight");
            hist_miss_momentum_energy_meas.Write();

            TH1D hist_sqrts_kin                     = *rdf.Histo1D({("sqrts_kin_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_kin","accidweight");
            hist_sqrts_kin.Write();
            TH1D hist_minust_kin                    = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 30, 0.0, 3.0},"minust_kin","accidweight");
            hist_minust_kin.Write();
            TH1D hist_minusu_kin                    = *rdf.Histo1D({("minusu_kin_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 300, 15.0, 45.0},"minusu_kin","accidweight");
            hist_minusu_kin.Write();
            TH1D hist_coplanarity_kin               = *rdf.Histo1D({("coplanarity_kin_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_kin","accidweight");
            hist_coplanarity_kin.Write();
            TH1D hist_thetaCM_kin                   = *rdf.Histo1D({("thetaCM_kin_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_kin","accidweight");
            hist_thetaCM_kin.Write();
            TH2D hist_minust_thetaCM_kin            = *rdf.Histo2D({("minust_thetaCM_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_kin","thetaCM_kin","accidweight");
            hist_minust_thetaCM_kin.Write();
            TH2D hist_beam_energy_minust_kin        = *rdf.Histo2D({("beam_energy_minust_kin_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_kin","minust_kin","accidweight");
            hist_beam_energy_minust_kin.Write();
            TH1D hist_rho_mass_kin                  = *rdf.Histo1D({("rho_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_mass_kin","accidweight");
            hist_rho_mass_kin.Write();
            TH1D hist_yphi_kin                      = *rdf.Histo1D({("yphi_kin_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_kin","accidweight");
            hist_yphi_kin.Write();
            TH1D hist_costheta_helicity_kin         = *rdf.Histo1D({("costheta_helicity_kin_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"costheta_helicity_kin","accidweight");
            hist_costheta_helicity_kin.Write();
            TH1D hist_phi_helicity_kin              = *rdf.Histo1D({("phi_helicity_kin_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"phi_helicity_kin","accidweight");
            hist_phi_helicity_kin.Write();
            TH1D hist_kinfit_fom_kin                = *rdf.Histo1D({("kinfit_fom_kin_"+ label).c_str(), ";KinFit FOM;Counts", 100, 0.0, 1.0},"kinfit_fom_kin","accidweight");
            hist_kinfit_fom_kin.Write();

            if (reaction.find("sim") != string::npos)
            {
                TH1D hist_beam_energy_truth                 = *rdf.Histo1D({("beam_energy_truth_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_truth");
                hist_beam_energy_truth.Write();
                TH2D hist_kp_kinematics_truth               = *rdf.Histo2D({("kp_kinematics_truth_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_truth","kp_theta_truth");
                hist_kp_kinematics_truth.Write();
                TH2D hist_km_kinematics_truth               = *rdf.Histo2D({("km_kinematics_truth_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_truth","km_theta_truth");
                hist_km_kinematics_truth.Write();
                TH2D hist_d_kinematics_truth                = *rdf.Histo2D({("d_kinematics_truth_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_truth","d_theta_truth");
                hist_d_kinematics_truth.Write();
                TH1D hist_phi_mass_truth                    = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_truth");
                hist_phi_mass_truth.Write();
                TH2D hist_phi_kinematics_truth              = *rdf.Histo2D({("phi_kinematics_truth_"+ label).c_str(), ";P_{#phi} (GeV/c);#theta_{#phi} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"phi_momentum_truth","phi_theta_truth");
                hist_phi_kinematics_truth.Write();
                TH2D hist_phi_d_theta_truth                 = *rdf.Histo2D({("phi_d_theta_truth_"+ label).c_str(), ";#theta_{d} (deg);#theta_{#phi} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"d_theta_truth","phi_theta_truth");
                hist_phi_d_theta_truth.Write();
                TH2D hist_phi_d_momentum_truth              = *rdf.Histo2D({("phi_d_momentum_truth_"+ label).c_str(), ";P_{d} (GeV/c);P_{#phi} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"d_momentum_truth","phi_momentum_truth");
                hist_phi_d_momentum_truth.Write();

                TH1D hist_miss_mass_truth                   = *rdf.Histo1D({("miss_mass_truth_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 400, -4.0, 4.0},"miss_mass_truth");
                hist_miss_mass_truth.Write();
                TH1D hist_miss_masssquared_truth            = *rdf.Histo1D({("miss_masssquared_truth_"+ label).c_str(), ";m_{miss}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"miss_masssquared_truth");
                hist_miss_masssquared_truth.Write();
                TH1D hist_miss_momentum_truth               = *rdf.Histo1D({("miss_momentum_truth_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_truth");
                hist_miss_momentum_truth.Write();
                TH1D hist_miss_pminus_truth                 = *rdf.Histo1D({("miss_pminus_truth_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"miss_pminus_truth");
                hist_miss_pminus_truth.Write();
                TH1D hist_miss_energy_balance_truth         = *rdf.Histo1D({("miss_energy_balance_truth_"+ label).c_str(), ";E_{miss} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"miss_energy_balance_truth");
                hist_miss_energy_balance_truth.Write();
                TH2D hist_miss_momentum_pminus_truth        = *rdf.Histo2D({("miss_momentum_pminus_truth_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"miss_momentum_truth","miss_pminus_truth");
                hist_miss_momentum_pminus_truth.Write();
                TH2D hist_miss_momentum_energy_truth        = *rdf.Histo2D({("miss_momentum_energy_truth_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"miss_momentum_truth","miss_energy_truth");
                hist_miss_momentum_energy_truth.Write();

                TH1D hist_sqrts_truth                       = *rdf.Histo1D({("sqrts_truth_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_truth");
                hist_sqrts_truth.Write();
                TH1D hist_minust_truth                      = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minust_truth");
                hist_minust_truth.Write();
                TH1D hist_minusu_truth                      = *rdf.Histo1D({("minusu_truth_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 300, 15.0, 45.0},"minusu_truth");
                hist_minusu_truth.Write();
                TH1D hist_coplanarity_truth                 = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_truth");
                hist_coplanarity_truth.Write();
                TH1D hist_thetaCM_truth                     = *rdf.Histo1D({("thetaCM_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_truth");
                hist_thetaCM_truth.Write();
                TH2D hist_minust_thetaCM_truth              = *rdf.Histo2D({("minust_thetaCM_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_truth","thetaCM_truth");
                hist_minust_thetaCM_truth.Write();
                TH2D hist_beam_energy_minust_truth          = *rdf.Histo2D({("beam_energy_minust_truth_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_truth","minust_truth");
                hist_beam_energy_minust_truth.Write();
                TH1D hist_rho_mass_truth                    = *rdf.Histo1D({("rho_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_mass_truth");
                hist_rho_mass_truth.Write();
                TH1D hist_yphi_truth                        = *rdf.Histo1D({("yphi_truth_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_truth");
                hist_yphi_truth.Write();
                TH1D hist_costheta_helicity_truth           = *rdf.Histo1D({("costheta_helicity_truth_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"costheta_helicity_truth");
                hist_costheta_helicity_truth.Write();
                TH1D hist_phi_helicity_truth                = *rdf.Histo1D({("phi_helicity_truth_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"phi_helicity_truth");
                hist_phi_helicity_truth.Write();

                TH2D hist_beam_energy_kin_truth             = *rdf.Histo2D({("beam_energy_kin_truth_"+ label).c_str(), ";E_{beam} (GeV);E_{beam} (GeV)", 60, 5.0, 11.0, 60, 5.0, 11.0},"beam_energy_kin","beam_energy_truth");
                hist_beam_energy_kin_truth.Write();
                TH2D hist_minust_kin_truth                  = *rdf.Histo2D({("minust_kin_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});-t (GeV^{2}/c^{2})", 30, 0.0, 3.0, 30, 0.0, 3.0},"minust_kin","minust_truth");
                hist_minust_kin_truth.Write();
            }
        }
        output_hist_file->Close();
    }
    cout << "Done!\n";
}