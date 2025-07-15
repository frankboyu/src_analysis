#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;
double mass_missing = 0.0;

double sim_weight_func(double beam_energy_truth, double minust_truth)
{
    double a1 = 2813.72894997;
    double b1 = 15.13997936;
    double a2 = 17.88792021;
    double b2 = 2.98839991;
    double normalization = 10;
    if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
        return 1.0;
    else                            // simulation, weighted by the measured cross section
        return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}

void filter_phi_d_recon_inc(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_recon_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_phi_d_recon";
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());

    cout << "Defining data frame...\n";
    if      (reaction.find("2H")   != string::npos)
    {
        mass_target = mass_2H;
        mass_missing = 0.0;
    }
    else if (reaction.find("4He")  != string::npos)
    {
        mass_target = mass_4He;
        mass_missing = mass_2H;
    }
    else if (reaction.find("12C")  != string::npos)
    {
        mass_target = mass_12C;
        mass_missing = mass_10B;
    }

    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    auto rdf_input = rdf_def
    .Define("sim_weight",                       "sim_weight_func(beam_p4_truth.E(), -(beam_p4_truth - phi_p4_truth).Mag2())")
    .Define("event_weight",                     "accidental_weight*sim_weight")
    .Define("target_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_energy_meas",                 "beam_p4_meas.E()")
    .Define("beam_energy_kin",                  "beam_p4_kin.E()")
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",                 "rftime + (beam_x4_meas.Z()-65.0)/29.9792458 - beam_x4_meas.T()")
    .Define("beam_DeltaT_kin",                  "rftime + (beam_x4_kin.Z()-65.0)/29.9792458 - beam_x4_kin.T()")

    .Define("kp_as_pion_p4_meas",               "TLorentzVector(kp_p4_meas.Vect(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_piplus*mass_piplus))")
    .Define("kp_as_pion_p4_kin",                "TLorentzVector(kp_p4_kin.Vect(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_piplus*mass_piplus))")
    .Define("kp_as_pion_p4_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_energy_meas",                   "kp_p4_meas.E()")
    .Define("kp_energy_kin",                    "kp_p4_kin.E()")
    .Define("kp_energy_truth",                  "kp_p4_truth.E()")
    .Define("kp_momentum_meas",                 "kp_p4_meas.P()")
    .Define("kp_momentum_kin",                  "kp_p4_kin.P()")
    .Define("kp_momentum_truth",                "kp_p4_truth.P()")
    .Define("kp_theta_meas",                    "kp_p4_meas.Theta()*RadToDeg")
    .Define("kp_theta_kin",                     "kp_p4_kin.Theta()*RadToDeg")
    .Define("kp_theta_truth",                   "kp_p4_truth.Theta()*RadToDeg")
    .Define("kp_DeltaT_meas",                   "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",                    "rftime + (kp_x4_kin.Z()-65.0)/29.9792458 - kp_x4_kin.T()")
    .Define("kp_in_fdc_meas",                   "event_weight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_in_cdc_meas",                   "event_weight*(kp_dedx_cdc > 0.0 && kp_dedx_fdc == 0.0)")
    .Define("kp_in_fdc_cdc_meas",               "event_weight*(kp_dedx_fdc > 0.0 && kp_dedx_cdc > 0.0)")
    .Define("kp_in_neither_meas",               "event_weight*(kp_dedx_fdc == 0.0 && kp_dedx_cdc == 0.0)")
    .Define("kp_dedx_fdc_keV_per_cm_meas",      "kp_dedx_fdc*1e6")
    .Define("kp_dedx_cdc_keV_per_cm_meas",      "kp_dedx_cdc*1e6")
    .Define("kp_dedx_st_keV_per_cm_meas",       "kp_dedx_st*1e3")
    .Define("kp_dedx_tof_keV_per_cm_meas",      "kp_dedx_tof*1e3")

    .Define("km_as_pion_p4_meas",               "TLorentzVector(km_p4_meas.Vect(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_piminus*mass_piminus))")
    .Define("km_as_pion_p4_kin",                "TLorentzVector(km_p4_kin.Vect(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_piminus*mass_piminus))")
    .Define("km_as_pion_p4_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_energy_meas",                   "km_p4_meas.E()")
    .Define("km_energy_kin",                    "km_p4_kin.E()")
    .Define("km_energy_truth",                  "km_p4_truth.E()")
    .Define("km_momentum_meas",                 "km_p4_meas.P()")
    .Define("km_momentum_kin",                  "km_p4_kin.P()")
    .Define("km_momentum_truth",                "km_p4_truth.P()")
    .Define("km_theta_meas",                    "km_p4_meas.Theta()*RadToDeg")
    .Define("km_theta_kin",                     "km_p4_kin.Theta()*RadToDeg")
    .Define("km_theta_truth",                   "km_p4_truth.Theta()*RadToDeg")
    .Define("km_DeltaT_meas",                   "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",                    "rftime + (km_x4_kin.Z()-65.0)/29.9792458 - km_x4_kin.T()")
    .Define("km_in_fdc_meas",                   "event_weight*(km_dedx_fdc > 0.0 && km_dedx_cdc == 0.0)")
    .Define("km_in_cdc_meas",                   "event_weight*(km_dedx_cdc > 0.0 && km_dedx_fdc == 0.0)")
    .Define("km_in_fdc_cdc_meas",               "event_weight*(km_dedx_fdc > 0.0 && km_dedx_cdc > 0.0)")
    .Define("km_in_neither_meas",               "event_weight*(km_dedx_fdc == 0.0 && km_dedx_cdc == 0.0)")
    .Define("km_dedx_fdc_keV_per_cm_meas",      "km_dedx_fdc*1e6")
    .Define("km_dedx_cdc_keV_per_cm_meas",      "km_dedx_cdc*1e6")
    .Define("km_dedx_st_keV_per_cm_meas",       "km_dedx_st*1e3")
    .Define("km_dedx_tof_keV_per_cm_meas",      "km_dedx_tof*1e3")

    .Define("d_energy_meas",                    "d_p4_meas.E()")
    .Define("d_energy_kin",                     "d_p4_kin.E()")
    .Define("d_energy_truth",                   "d_p4_truth.E()")
    .Define("d_momentum_meas",                  "d_p4_meas.P()")
    .Define("d_momentum_kin",                   "d_p4_kin.P()")
    .Define("d_momentum_truth",                 "d_p4_truth.P()")
    .Define("d_theta_meas",                     "d_p4_meas.Theta()*RadToDeg")
    .Define("d_theta_kin",                      "d_p4_kin.Theta()*RadToDeg")
    .Define("d_theta_truth",                    "d_p4_truth.Theta()*RadToDeg")
    .Define("d_DeltaT_meas",                    "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_meas.T()")
    .Define("d_DeltaT_kin",                     "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_kin.T()")
    .Define("d_dedx_fdc_keV_per_cm_meas",       "d_dedx_fdc*1e6")
    .Define("d_dedx_cdc_keV_per_cm_meas",       "d_dedx_cdc*1e6")
    .Define("d_dedx_st_keV_per_cm_meas",        "d_dedx_st*1e3")
    .Define("d_dedx_tof_keV_per_cm_meas",       "d_dedx_tof*1e3")

    .Define("phi_p4_meas",                      "kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin",                       "kp_p4_kin + km_p4_kin")
    .Define("phi_p4_truth",                     "kp_p4_truth + km_p4_truth")
    .Define("phi_energy_meas",                  "phi_p4_meas.E()")
    .Define("phi_energy_kin",                   "phi_p4_kin.E()")
    .Define("phi_energy_truth",                 "phi_p4_truth.E()")
    .Define("phi_momentum_meas",                "phi_p4_meas.P()")
    .Define("phi_momentum_kin",                 "phi_p4_kin.P()")
    .Define("phi_momentum_truth",               "phi_p4_truth.P()")
    .Define("phi_mass_meas",                    "phi_p4_meas.M()")
    .Define("phi_mass_kin",                     "phi_p4_kin.M()")
    .Define("phi_mass_truth",                   "phi_p4_truth.M()")
    .Define("phi_proxymass_meas",               "TMath::Sqrt(phi_p4_meas.Minus()*(2*beam_energy_meas + mass_target - d_p4_meas.Plus() - (TMath::Sq(mass_missing) + (phi_p4_meas + d_p4_meas).Perp2())/(mass_target - (phi_p4_meas + d_p4_meas).Minus())) - phi_p4_meas.Perp2())")
    .Define("phi_proxymass_kin",                "TMath::Sqrt(phi_p4_kin.Minus()*(2*beam_energy_kin + mass_target - d_p4_kin.Plus() - (TMath::Sq(mass_missing) + (phi_p4_kin + d_p4_kin).Perp2())/(mass_target - (phi_p4_kin + d_p4_kin).Minus())) - phi_p4_kin.Perp2())")
    .Define("phi_proxymass_truth",              "TMath::Sqrt(phi_p4_truth.Minus()*(2*beam_energy_truth + mass_target - d_p4_truth.Plus() - (TMath::Sq(mass_missing) + (phi_p4_truth + d_p4_truth).Perp2())/(mass_target - (phi_p4_truth + d_p4_truth).Minus())) - phi_p4_truth.Perp2())")
    .Define("phi_theta_meas",                   "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",                    "phi_p4_kin.Theta()*RadToDeg")
    .Define("phi_theta_truth",                  "phi_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_meas",                   "phi_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin",                    "phi_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("struck_p4_truth",                  "phi_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("struck_energy_meas",               "struck_p4_meas.E()")
    .Define("struck_energy_kin",                "struck_p4_kin.E()")
    .Define("struck_energy_truth",              "struck_p4_truth.E()")
    .Define("struck_energy_balance_meas",       "struck_energy_meas - mass_2H")
    .Define("struck_energy_balance_kin",        "struck_energy_kin - mass_2H")
    .Define("struck_energy_balance_truth",      "struck_energy_truth - mass_2H")
    .Define("struck_mass_meas",                 "struck_p4_meas.M()")
    .Define("struck_mass_kin",                  "struck_p4_kin.M()")
    .Define("struck_mass_truth",                "struck_p4_truth.M()")
    .Define("struck_masssquared_meas",          "struck_p4_meas.M2()")
    .Define("struck_masssquared_kin",           "struck_p4_kin.M2()")
    .Define("struck_masssquared_truth",         "struck_p4_truth.M2()")
    .Define("struck_pminus_meas",               "struck_p4_meas.Minus()")
    .Define("struck_pminus_kin",                "struck_p4_kin.Minus()")
    .Define("struck_pminus_truth",              "struck_p4_truth.Minus()")
    .Define("struck_momentum_meas",             "struck_p4_meas.P()")
    .Define("struck_momentum_kin",              "struck_p4_kin.P()")
    .Define("struck_momentum_truth",            "struck_p4_truth.P()")
    .Define("struck_theta_meas",                "struck_p4_meas.Theta()*RadToDeg")
    .Define("struck_theta_kin",                 "struck_p4_kin.Theta()*RadToDeg")
    .Define("struck_theta_truth",               "struck_p4_truth.Theta()*RadToDeg")

    .Define("miss_p4_meas",                     "beam_p4_meas + TLorentzVector(0, 0, 0, mass_target) - phi_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                      "beam_p4_kin + TLorentzVector(0, 0, 0, mass_target) - phi_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                    "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target) - phi_p4_truth - d_p4_truth")
    .Define("miss_energy_meas",                 "miss_p4_meas.E()")
    .Define("miss_energy_kin",                  "miss_p4_kin.E()")
    .Define("miss_energy_truth",                "miss_p4_truth.E()")
    .Define("miss_energy_balance_meas",         "miss_energy_meas - mass_missing")
    .Define("miss_energy_balance_kin",          "miss_energy_kin - mass_missing")
    .Define("miss_energy_balance_truth",        "miss_energy_truth - mass_missing")
    .Define("miss_mass_meas",                   "miss_p4_meas.M()")
    .Define("miss_mass_kin",                    "miss_p4_kin.M()")
    .Define("miss_mass_truth",                  "miss_p4_truth.M()")
    .Define("miss_mass_balance_meas",           "miss_mass_meas - mass_missing")
    .Define("miss_mass_balance_kin",            "miss_mass_kin - mass_missing")
    .Define("miss_mass_balance_truth",          "miss_mass_truth - mass_missing")
    .Define("miss_masssquared_meas",            "miss_p4_meas.M2()")
    .Define("miss_masssquared_kin",             "miss_p4_kin.M2()")
    .Define("miss_masssquared_truth",           "miss_p4_truth.M2()")
    .Define("miss_masssquared_balance_meas",    "miss_masssquared_meas - mass_missing*mass_missing")
    .Define("miss_masssquared_balance_kin",     "miss_masssquared_kin - mass_missing*mass_missing")
    .Define("miss_masssquared_balance_truth",   "miss_masssquared_truth - mass_missing*mass_missing")
    .Define("miss_pminus_meas",                 "miss_p4_meas.Minus()")
    .Define("miss_pminus_kin",                  "miss_p4_kin.Minus()")
    .Define("miss_pminus_truth",                "miss_p4_truth.Minus()")
    .Define("miss_pminus_balance_meas",         "miss_pminus_meas - mass_missing")
    .Define("miss_pminus_balance_kin",          "miss_pminus_kin - mass_missing")
    .Define("miss_pminus_balance_truth",        "miss_pminus_truth - mass_missing")
    .Define("miss_momentum_meas",               "miss_p4_meas.P()")
    .Define("miss_momentum_kin",                "miss_p4_kin.P()")
    .Define("miss_momentum_truth",              "miss_p4_truth.P()")
    .Define("miss_theta_meas",                  "miss_p4_meas.Theta()*RadToDeg")
    .Define("miss_theta_kin",                   "miss_p4_kin.Theta()*RadToDeg")
    .Define("miss_theta_truth",                 "miss_p4_truth.Theta()*RadToDeg")

    .Define("total_initial_p4_meas",            "beam_p4_meas + TLorentzVector(0, 0, 0, mass_target)")
    .Define("total_initial_p4_kin",             "beam_p4_kin + TLorentzVector(0, 0, 0, mass_target)")
    .Define("total_initial_p4_truth",           "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target)")
    .Define("total_final_p4_meas",              "phi_p4_meas + d_p4_meas")
    .Define("total_final_p4_kin",               "phi_p4_kin + d_p4_kin")
    .Define("total_final_p4_truth",             "phi_p4_truth + d_p4_truth")
    .Define("minust_beam_meas",                 "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minust_beam_kin",                  "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minust_beam_truth",                "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minust_target_meas",               "-(target_p4 - d_p4_meas).Mag2()")
    .Define("minust_target_kin",                "-(target_p4 - d_p4_kin).Mag2()")
    .Define("minust_target_truth",              "-(target_p4 - d_p4_truth).Mag2()")
    .Define("y_beam_meas",                       "minust_beam_meas/(2*mass_2H*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_beam_kin",                        "minust_beam_kin/(2*mass_2H*(beam_p4_kin.E()-phi_p4_kin.E()))")
    .Define("y_beam_truth",                      "minust_beam_truth/(2*mass_2H*(beam_p4_truth.E()-phi_p4_truth.E()))")
    .Define("y_target_meas",                     "minust_target_meas/(2*mass_2H*(target_p4.E()-d_p4_meas.E()))")
    .Define("y_target_kin",                      "minust_target_kin/(2*mass_2H*(target_p4.E()-d_p4_kin.E()))")
    .Define("y_target_truth",                    "minust_target_truth/(2*mass_2H*(target_p4.E()-d_p4_truth.E()))")
    .Define("coplanarity_meas",                 "abs(phi_p4_meas.Phi() - d_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin",                  "abs(phi_p4_kin.Phi() - d_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth",                "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("rho_mass_meas",                    "(kp_as_pion_p4_meas + km_as_pion_p4_meas).M()")
    .Define("rho_mass_kin",                     "(kp_as_pion_p4_kin + km_as_pion_p4_kin).M()")
    .Define("rho_mass_truth",                   "(kp_as_pion_p4_truth + km_as_pion_p4_truth).M()")
    .Define("chisq_per_ndf_kin",                "kin_chisq/kin_ndf")
    .Define("kinfit_fom_kin",                   "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("log10_kinfit_fom_kin",             "TMath::Log10(kinfit_fom_kin)")

    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_meas",                 "boost_lorentz_vector(beam_p4_meas, -total_p4_meas.BoostVector())")
    .Define("beam_p4_com_kin",                  "boost_lorentz_vector(beam_p4_kin, -total_p4_kin.BoostVector())")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_p4_truth.BoostVector())")
    .Define("phi_p4_com_meas",                  "boost_lorentz_vector(phi_p4_meas, -total_p4_meas.BoostVector())")
    .Define("phi_p4_com_kin",                   "boost_lorentz_vector(phi_p4_kin, -total_p4_kin.BoostVector())")
    .Define("phi_p4_com_truth",                 "boost_lorentz_vector(phi_p4_truth, -total_p4_truth.BoostVector())")
    .Define("kp_p4_com_meas",                   "boost_lorentz_vector(kp_p4_meas, -total_p4_meas.BoostVector())")
    .Define("kp_p4_com_kin",                    "boost_lorentz_vector(kp_p4_kin, -total_p4_kin.BoostVector())")
    .Define("kp_p4_com_truth",                  "boost_lorentz_vector(kp_p4_truth, -total_p4_truth.BoostVector())")
    .Define("km_p4_com_meas",                   "boost_lorentz_vector(km_p4_meas, -total_p4_meas.BoostVector())")
    .Define("km_p4_com_kin",                    "boost_lorentz_vector(km_p4_kin, -total_p4_kin.BoostVector())")
    .Define("km_p4_com_truth",                  "boost_lorentz_vector(km_p4_truth, -total_p4_truth.BoostVector())")
    .Define("z_x3_com_meas",                    "beam_p4_com_meas.Vect().Unit()")
    .Define("z_x3_com_kin",                     "beam_p4_com_kin.Vect().Unit()")
    .Define("z_x3_com_truth",                   "beam_p4_com_truth.Vect().Unit()")
    .Define("y_x3_com_meas",                    "beam_p4_com_meas.Vect().Cross(phi_p4_com_meas.Vect()).Unit()")
    .Define("y_x3_com_kin",                     "beam_p4_com_kin.Vect().Cross(phi_p4_com_kin.Vect()).Unit()")
    .Define("y_x3_com_truth",                   "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_com_meas",                    "y_x3_com_meas.Cross(z_x3_com_meas).Unit()")
    .Define("x_x3_com_kin",                     "y_x3_com_kin.Cross(z_x3_com_kin).Unit()")
    .Define("x_x3_com_truth",                   "y_x3_com_truth.Cross(z_x3_com_truth).Unit()")
    .Define("polarization_phi_com_meas",        "TMath::ATan2(-x_x3_com_meas.Dot(epsilon_x3_com.Cross(z_x3_com_meas)), y_x3_com_meas.Dot(epsilon_x3_com.Cross(z_x3_com_meas)))*RadToDeg")
    .Define("polarization_phi_com_kin",         "TMath::ATan2(-x_x3_com_kin.Dot(epsilon_x3_com.Cross(z_x3_com_kin)), y_x3_com_kin.Dot(epsilon_x3_com.Cross(z_x3_com_kin)))*RadToDeg")
    .Define("polarization_phi_com_truth",       "TMath::ATan2(-x_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)), y_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)))*RadToDeg")
    .Define("scatter_theta_com_meas",           "phi_p4_com_meas.Vect().Angle(z_x3_com_meas)*RadToDeg")
    .Define("scatter_theta_com_kin",            "phi_p4_com_kin.Vect().Angle(z_x3_com_kin)*RadToDeg")
    .Define("scatter_theta_com_truth",          "phi_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")

    .Define("kp_p4_helicity_meas",              "boost_lorentz_vector(kp_p4_com_meas, -phi_p4_com_meas.BoostVector())")
    .Define("kp_p4_helicity_kin",               "boost_lorentz_vector(kp_p4_com_kin, -phi_p4_com_kin.BoostVector())")
    .Define("kp_p4_helicity_truth",             "boost_lorentz_vector(kp_p4_com_truth, -phi_p4_com_truth.BoostVector())")
    .Define("km_p4_helicity_meas",              "boost_lorentz_vector(km_p4_com_meas, -phi_p4_com_meas.BoostVector())")
    .Define("km_p4_helicity_kin",               "boost_lorentz_vector(km_p4_com_kin, -phi_p4_com_kin.BoostVector())")
    .Define("km_p4_helicity_truth",             "boost_lorentz_vector(km_p4_com_truth, -phi_p4_com_truth.BoostVector())")
    .Define("z_x3_helicity_meas",               "phi_p4_com_meas.Vect().Unit()")
    .Define("z_x3_helicity_kin",                "phi_p4_com_kin.Vect().Unit()")
    .Define("z_x3_helicity_truth",              "phi_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_meas",               "beam_p4_com_meas.Vect().Cross(phi_p4_com_meas.Vect()).Unit()")
    .Define("y_x3_helicity_kin",                "beam_p4_com_kin.Vect().Cross(phi_p4_com_kin.Vect()).Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_meas",               "y_x3_helicity_meas.Cross(z_x3_helicity_meas).Unit()")
    .Define("x_x3_helicity_kin",                "y_x3_helicity_kin.Cross(z_x3_helicity_kin).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_meas",              "kp_p4_helicity_meas.Vect().Unit()")
    .Define("pi_x3_helicity_kin",               "kp_p4_helicity_kin.Vect().Unit()")
    .Define("pi_x3_helicity_truth",             "kp_p4_helicity_truth.Vect().Unit()")
    .Define("decay_costheta_helicity_meas",     "pi_x3_helicity_meas.Dot(z_x3_helicity_meas)")
    .Define("decay_costheta_helicity_kin",      "pi_x3_helicity_kin.Dot(z_x3_helicity_kin)")
    .Define("decay_costheta_helicity_truth",    "pi_x3_helicity_truth.Dot(z_x3_helicity_truth)")
    .Define("decay_phi_helicity_meas",          "TMath::ATan2(-x_x3_helicity_meas.Dot(pi_x3_helicity_meas.Cross(z_x3_helicity_meas)), y_x3_helicity_meas.Dot(pi_x3_helicity_meas.Cross(z_x3_helicity_meas)))*RadToDeg")
    .Define("decay_phi_helicity_kin",           "TMath::ATan2(-x_x3_helicity_kin.Dot(pi_x3_helicity_kin.Cross(z_x3_helicity_kin)), y_x3_helicity_kin.Dot(pi_x3_helicity_kin.Cross(z_x3_helicity_kin)))*RadToDeg")
    .Define("decay_phi_helicity_truth",         "TMath::ATan2(-x_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)), y_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)))*RadToDeg")
    .Define("psi_helicity_meas",                "fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0) >= 180 ? fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0) - 360 : fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0)")
    .Define("psi_helicity_kin",                 "fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0) >= 180 ? fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0) - 360 : fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0)")
    .Define("psi_helicity_truth",               "fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0)")
    ;

    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_dEdxCut        = rdf_NoCut.Filter("(d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-29.68353898*d_momentum_meas+13.50623694)+17.88279645*d_momentum_meas*d_momentum_meas-42.15473796*d_momentum_meas+28.83200736)) && (d_dedx_cdc_keV_per_cm_meas < (TMath::Exp(-26.69276323*d_momentum_meas+15.92466317)+17.1164272*d_momentum_meas*d_momentum_meas-48.7542903*d_momentum_meas+40.25692313))");
    auto rdf_KinFitFOMCut   = rdf_dEdxCut.Filter("kinfit_fom_kin > 1e-5");
    auto rdf_PIDFOMCut      = rdf_KinFitFOMCut.Filter("(kp_pidfom > 0.01) && (km_pidfom > 0.01)");
    auto rdf_MissPCut       = rdf_PIDFOMCut.Filter("(abs(struck_energy_balance_kin) < 1.0)");
    auto rdf_PhiMassCut     = rdf_MissPCut.Filter("phi_mass_kin > 1.005 && phi_mass_kin < 1.04");
    auto rdf_output         = rdf_PhiMassCut;
    RNode rdfs []           = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_MissPCut,   rdf_PhiMassCut};
    string labels []        = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "MissPCut",     "PhiMassCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_phi_d_recon";
        rdf_output.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

            TH1D hist_num_unused_tracks_meas                = *rdf.Histo1D({("num_unused_tracks_meas_"+ label).c_str(), ";Number of unused tracks;Counts", 10, 0.0, 10.0},"num_unused_tracks","event_weight");
            hist_num_unused_tracks_meas.Write();
            TH1D hist_num_unused_showers_meas               = *rdf.Histo1D({("num_unused_showers_meas_"+ label).c_str(), ";Number of unused showers;Counts", 10, 0.0, 10.0},"num_unused_showers","event_weight");
            hist_num_unused_showers_meas.Write();

            TH1D hist_beam_energy_kin                       = *rdf.Histo1D({("beam_energy_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_kin","event_weight");
            hist_beam_energy_kin.Write();
            TH1D hist_beam_DeltaT_kin                       = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 100, -20.0, 20.0},"beam_DeltaT_kin","event_weight");
            hist_beam_DeltaT_kin.Write();

            TH1D hist_kp_pidfom                             = *rdf.Histo1D({("kp_pidfom_"+ label).c_str(), ";kp_pidfom;Counts", 100, 0.0, 1.0},"kp_pidfom","event_weight");
            hist_kp_pidfom.Write();
            TH1D hist_kp_DeltaT_kin                         = *rdf.Histo1D({("kp_DeltaT_kin_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 100, -2.0, 2.0},"kp_DeltaT_kin","event_weight");
            hist_kp_DeltaT_kin.Write();
            TH2D hist_kp_DeltaT_momentum_kin                = *rdf.Histo2D({("kp_DeltaT_momentum_kin_"+ label).c_str(), ";p (GeV/c);#Delta t_{K^{+}} (ns)", 100, 0.0, 10.0, 100, -2.0, 2.0},"kp_momentum_kin","kp_DeltaT_kin","event_weight");
            hist_kp_DeltaT_momentum_kin.Write();
            TH2D hist_kp_dEdx_cdc_meas                      = *rdf.Histo2D({("kp_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_cdc_meas.Write();
            TH2D hist_kp_dEdx_fdc_meas                      = *rdf.Histo2D({("kp_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_fdc_meas.Write();
            TH2D hist_kp_dEdx_tof_meas                      = *rdf.Histo2D({("kp_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_tof_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_tof_meas.Write();
            TH2D hist_kp_dEdx_st_meas                       = *rdf.Histo2D({("kp_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_st_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_st_meas.Write();
            TH2D hist_kp_kinematics_kin                     = *rdf.Histo2D({("kp_kinematics_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","event_weight");
            hist_kp_kinematics_kin.Write();
            TH2D hist_kp_kinematics_fdc_kin                 = *rdf.Histo2D({("kp_kinematics_fdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc_meas");
            hist_kp_kinematics_fdc_kin.Write();
            TH2D hist_kp_kinematics_fdc_cdc_kin             = *rdf.Histo2D({("kp_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_fdc_cdc_meas");
            hist_kp_kinematics_fdc_cdc_kin.Write();
            TH2D hist_kp_kinematics_cdc_kin                 = *rdf.Histo2D({("kp_kinematics_cdc_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_cdc_meas");
            hist_kp_kinematics_cdc_kin.Write();
            TH2D hist_kp_kinematics_neither_kin             = *rdf.Histo2D({("kp_kinematics_neither_kin_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_kin","kp_theta_kin","kp_in_neither_meas");
            hist_kp_kinematics_neither_kin.Write();

            TH1D hist_km_pidfom                             = *rdf.Histo1D({("km_pidfom_"+ label).c_str(), ";km_pidfom;Counts", 100, 0.0, 1.0},"km_pidfom","event_weight");
            hist_km_pidfom.Write();
            TH1D hist_km_DeltaT_kin                         = *rdf.Histo1D({("km_DeltaT_kin_"+ label).c_str(), ";#Delta t_{K^{-}} (ns);Counts", 100, -2.0, 2.0},"km_DeltaT_kin","event_weight");
            hist_km_DeltaT_kin.Write();
            TH2D hist_km_DeltaT_momentum_kin                = *rdf.Histo2D({("km_DeltaT_momentum_kin_"+ label).c_str(), ";p (GeV/c);#Delta t_{K^{-}} (ns)", 100, 0.0, 10.0, 100, -2.0, 2.0},"km_momentum_kin","km_DeltaT_kin","event_weight");
            hist_km_DeltaT_momentum_kin.Write();
            TH2D hist_km_dEdx_cdc_meas                      = *rdf.Histo2D({("km_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_km_dEdx_cdc_meas.Write();
            TH2D hist_km_dEdx_fdc_meas                      = *rdf.Histo2D({("km_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_km_dEdx_fdc_meas.Write();
            TH2D hist_km_dEdx_tof_meas                      = *rdf.Histo2D({("km_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_tof_keV_per_cm_meas","event_weight");
            hist_km_dEdx_tof_meas.Write();
            TH2D hist_km_dEdx_st_meas                       = *rdf.Histo2D({("km_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_st_keV_per_cm_meas","event_weight");
            hist_km_dEdx_st_meas.Write();
            TH2D hist_km_kinematics_kin                     = *rdf.Histo2D({("km_kinematics_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","event_weight");
            hist_km_kinematics_kin.Write();
            TH2D hist_km_kinematics_fdc_kin                 = *rdf.Histo2D({("km_kinematics_fdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc_meas");
            hist_km_kinematics_fdc_kin.Write();
            TH2D hist_km_kinematics_fdc_cdc_kin             = *rdf.Histo2D({("km_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_fdc_cdc_meas");
            hist_km_kinematics_fdc_cdc_kin.Write();
            TH2D hist_km_kinematics_cdc_kin                 = *rdf.Histo2D({("km_kinematics_cdc_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_cdc_meas");
            hist_km_kinematics_cdc_kin.Write();
            TH2D hist_km_kinematics_neither_kin             = *rdf.Histo2D({("km_kinematics_neither_kin_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_kin","km_theta_kin","km_in_neither_meas");
            hist_km_kinematics_neither_kin.Write();

            TH1D hist_d_DeltaT_kin                          = *rdf.Histo1D({("d_DeltaT_kin_"+ label).c_str(), ";#Delta t_{d} (ns);Counts", 100, -5.0, 5.0},"d_DeltaT_kin","event_weight");
            hist_d_DeltaT_kin.Write();
            TH2D hist_d_DeltaT_momentum_kin                 = *rdf.Histo2D({("d_DeltaT_momentum_kin_"+ label).c_str(), ";p (GeV/c);#Delta t_{d} (ns)", 100, 0.0, 10.0, 100, -5.0, 5.0},"d_momentum_kin","d_DeltaT_kin","event_weight");
            hist_d_DeltaT_momentum_kin.Write();
            TH2D hist_d_dEdx_cdc_meas                       = *rdf.Histo2D({("d_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_d_dEdx_cdc_meas.Write();
            TH2D hist_d_dEdx_fdc_meas                       = *rdf.Histo2D({("d_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_d_dEdx_fdc_meas.Write();
            TH2D hist_d_dEdx_tof_meas                       = *rdf.Histo2D({("d_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_tof_keV_per_cm_meas","event_weight");
            hist_d_dEdx_tof_meas.Write();
            TH2D hist_d_dEdx_st_meas                        = *rdf.Histo2D({("d_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_st_keV_per_cm_meas","event_weight");
            hist_d_dEdx_st_meas.Write();
            TH2D hist_d_kinematics_kin                      = *rdf.Histo2D({("d_kinematics_kin_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_kin","d_theta_kin","event_weight");
            hist_d_kinematics_kin.Write();

            TH1D hist_phi_mass_kin                          = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_kin","event_weight");
            hist_phi_mass_kin.Write();
            TH1D hist_phi_mass_meas                         = *rdf.Histo1D({("phi_mass_meas_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","event_weight");
            hist_phi_mass_meas.Write();
            TH1D hist_phi_proxymass_kin                     = *rdf.Histo1D({("phi_proxymass_kin_"+ label).c_str(), ";m^{proxy}_{K^{+}K^{-}} (GeV/c^{2});Counts", 400, 0.9, 1.3},"phi_proxymass_kin","event_weight");
            hist_phi_proxymass_kin.Write();
            TH1D hist_phi_proxymass_meas                    = *rdf.Histo1D({("phi_proxymass_meas_"+ label).c_str(), ";m^{proxy}_{K^{+}K^{-}} (GeV/c^{2});Counts", 400, 0.9, 1.3},"phi_proxymass_meas","event_weight");
            hist_phi_proxymass_meas.Write();
            TH2D hist_phi_mass_kin_meas                     = *rdf.Histo2D({("phi_mass_kin_meas_"+ label).c_str(), ";m^{kin}_{K^{+}K^{-}} (GeV/c);m^{meas}_{K^{+}K^{-}} (GeV/c^{2})", 400, 0.9, 1.3, 400, 0.9, 1.3},"phi_mass_kin","phi_mass_meas","event_weight");
            hist_phi_mass_kin_meas.Write();
            TH2D hist_phi_mass_kinfit_fom_kin               = *rdf.Histo2D({("phi_mass_kinfit_fom_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);KinFit FOM", 400, 0.9, 1.3, 100, 0.0, 1.0},"phi_mass_kin","kinfit_fom_kin","event_weight");
            hist_phi_mass_kinfit_fom_kin.Write();
            TH2D hist_phi_kinematics_kin                    = *rdf.Histo2D({("phi_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 110, 0.0, 11.0, 180, 0.0, 180.0},"phi_momentum_kin","phi_theta_kin","event_weight");
            hist_phi_kinematics_kin.Write();

            TH1D hist_miss_energy_kin                       = *rdf.Histo1D({("miss_energy_kin_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_kin","event_weight");
            hist_miss_energy_kin.Write();
            TH1D hist_miss_energy_meas                      = *rdf.Histo1D({("miss_energy_meas_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_meas","event_weight");
            hist_miss_energy_meas.Write();
            TH1D hist_miss_mass_kin                         = *rdf.Histo1D({("miss_mass_kin_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_kin","event_weight");
            hist_miss_mass_kin.Write();
            TH1D hist_miss_mass_meas                        = *rdf.Histo1D({("miss_mass_meas_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_meas","event_weight");
            hist_miss_mass_meas.Write();
            TH1D hist_miss_masssquared_kin                  = *rdf.Histo1D({("miss_masssquared_kin_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, -0.1, 0.1},"miss_masssquared_balance_kin","event_weight");
            hist_miss_masssquared_kin.Write();
            TH1D hist_miss_masssquared_meas                 = *rdf.Histo1D({("miss_masssquared_meas_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, -0.1, 0.1},"miss_masssquared_balance_meas","event_weight");
            hist_miss_masssquared_meas.Write();
            TH1D hist_miss_pminus_kin                       = *rdf.Histo1D({("miss_pminus_kin_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_balance_kin","event_weight");
            hist_miss_pminus_kin.Write();
            TH1D hist_miss_pminus_meas                      = *rdf.Histo1D({("miss_pminus_meas_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_balance_meas","event_weight");
            hist_miss_pminus_meas.Write();
            TH1D hist_miss_momentum_kin                     = *rdf.Histo1D({("miss_momentum_kin_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_kin","event_weight");
            hist_miss_momentum_kin.Write();
            TH1D hist_miss_momentum_meas                    = *rdf.Histo1D({("miss_momentum_meas_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_meas","event_weight");
            hist_miss_momentum_meas.Write();
            TH2D hist_miss_momentum_pminus_kin              = *rdf.Histo2D({("miss_momentum_pminus_kin_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_kin","miss_pminus_balance_kin","event_weight");
            hist_miss_momentum_pminus_kin.Write();
            TH2D hist_miss_momentum_pminus_meas             = *rdf.Histo2D({("miss_momentum_pminus_meas_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_meas","miss_pminus_balance_meas","event_weight");
            hist_miss_momentum_pminus_meas.Write();
            TH2D hist_miss_momentum_energy_kin              = *rdf.Histo2D({("miss_momentum_energy_kin_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2}  (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_kin","miss_energy_balance_kin","event_weight");
            hist_miss_momentum_energy_kin.Write();
            TH2D hist_miss_momentum_energy_meas             = *rdf.Histo2D({("miss_momentum_energy_meas_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2} (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_meas","miss_energy_balance_meas","event_weight");
            hist_miss_momentum_energy_meas.Write();
            TH2D hist_miss_energy_phi_mass_meas             = *rdf.Histo2D({("miss_energy_phi_mass_meas_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);m_{K^{+}K^{-}} (GeV/c)", 300, -3.0, 3.0, 400, 0.9, 1.3},"miss_energy_balance_meas","phi_mass_meas","event_weight");
            hist_miss_energy_phi_mass_meas.Write();

            TH1D hist_minust_kin                            = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 30, 0.0, 3.0},"minust_kin","event_weight");
            hist_minust_kin.Write();
            TH1D hist_coplanarity_kin                       = *rdf.Histo1D({("coplanarity_kin_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_kin","event_weight");
            hist_coplanarity_kin.Write();
            TH1D hist_coplanarity_meas                      = *rdf.Histo1D({("coplanarity_meas_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_meas","event_weight");
            hist_coplanarity_meas.Write();
            TH1D hist_yphi_kin                              = *rdf.Histo1D({("yphi_kin_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_kin","event_weight");
            hist_yphi_kin.Write();
            TH1D hist_yphi_meas                             = *rdf.Histo1D({("yphi_meas_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_meas","event_weight");
            hist_yphi_meas.Write();
            TH1D hist_rho_mass_kin                          = *rdf.Histo1D({("rho_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 100, 0.0, 1.0},"rho_mass_kin","event_weight");
            hist_rho_mass_kin.Write();
            TH1D hist_kinfit_fom_kin                        = *rdf.Histo1D({("kinfit_fom_kin_"+ label).c_str(), ";log(KinFit FOM);Counts", 30, -15.0, 0},"log10_kinfit_fom_kin","event_weight");
            hist_kinfit_fom_kin.Write();
            TH1D hist_chisq_per_ndf_kin                     = *rdf.Histo1D({("chisq_per_ndf_kin_"+ label).c_str(), ";#chi^{2}/NDF;Counts", 20, 0.0, 10.0},"chisq_per_ndf_kin","event_weight");
            hist_chisq_per_ndf_kin.Write();
            TH2D hist_kinfit_fom_chisq_per_ndf_kin          = *rdf.Histo2D({("kinfit_fom_chisq_per_ndf_kin_"+ label).c_str(), ";log(KinFit FOM);#chi^{2}/NDF", 30, -15.0, 0.0, 20, 0.0, 10.0},"log10_kinfit_fom_kin","chisq_per_ndf_kin","event_weight");
            hist_kinfit_fom_chisq_per_ndf_kin.Write();
            TH2D hist_beam_energy_minust_kin                = *rdf.Histo2D({("beam_energy_minust_kin_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_kin","minust_kin","event_weight");
            hist_beam_energy_minust_kin.Write();

            TH1D hist_scatter_theta_com_kin                 = *rdf.Histo1D({("scatter_theta_com_kin_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"scatter_theta_com_kin","event_weight");
            hist_scatter_theta_com_kin.Write();
            TH2D hist_minust_scatter_theta_com_kin          = *rdf.Histo2D({("minust_scatter_theta_com_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 20, 0.0, 2.0, 40, 0.0, 40.0},"minust_kin","scatter_theta_com_kin","event_weight");
            hist_minust_scatter_theta_com_kin.Write();
            TH1D hist_polarization_phi_com_kin              = *rdf.Histo1D({("polarization_phi_com_kin_"+ label).c_str(), ";#phi_{com} (deg);Counts", 9, -180, 180.0},"polarization_phi_com_kin","event_weight");
            hist_polarization_phi_com_kin.Write();

            TH1D hist_decay_costheta_helicity_kin           = *rdf.Histo1D({("decay_costheta_helicity_kin_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"decay_costheta_helicity_kin","event_weight");
            hist_decay_costheta_helicity_kin.Write();
            TH1D hist_decay_phi_helicity_kin                = *rdf.Histo1D({("decay_phi_helicity_kin_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"decay_phi_helicity_kin","event_weight");
            hist_decay_phi_helicity_kin.Write();
            TH1D hist_psi_helicity_kin                      = *rdf.Histo1D({("psi_helicity_kin_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 18, -180.0, 180.0},"psi_helicity_kin","event_weight");
            hist_psi_helicity_kin.Write();

            if (reaction.find("sim") != string::npos)
            {
                TH1D hist_beam_energy_truth                 = *rdf.Histo1D({("beam_energy_truth_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_truth","event_weight");
                hist_beam_energy_truth.Write();
                TH2D hist_kp_kinematics_truth               = *rdf.Histo2D({("kp_kinematics_truth_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_truth","kp_theta_truth","event_weight");
                hist_kp_kinematics_truth.Write();
                TH2D hist_km_kinematics_truth               = *rdf.Histo2D({("km_kinematics_truth_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_truth","km_theta_truth","event_weight");
                hist_km_kinematics_truth.Write();
                TH2D hist_d_kinematics_truth                = *rdf.Histo2D({("d_kinematics_truth_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_truth","d_theta_truth","event_weight");
                hist_d_kinematics_truth.Write();
                TH1D hist_phi_mass_truth                    = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_truth","event_weight");
                hist_phi_mass_truth.Write();
                TH2D hist_phi_kinematics_truth              = *rdf.Histo2D({("phi_kinematics_truth_"+ label).c_str(), ";P_{#phi} (GeV/c);#theta_{#phi} (deg)", 100, 0.0, 11.0, 180, 0.0, 180.0},"phi_momentum_truth","phi_theta_truth","event_weight");
                hist_phi_kinematics_truth.Write();
                TH2D hist_phi_d_theta_truth                 = *rdf.Histo2D({("phi_d_theta_truth_"+ label).c_str(), ";#theta_{d} (deg);#theta_{#phi} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"d_theta_truth","phi_theta_truth","event_weight");
                hist_phi_d_theta_truth.Write();
                TH2D hist_phi_d_momentum_truth              = *rdf.Histo2D({("phi_d_momentum_truth_"+ label).c_str(), ";P_{d} (GeV/c);P_{#phi} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"d_momentum_truth","phi_momentum_truth","event_weight");
                hist_phi_d_momentum_truth.Write();

                TH1D hist_miss_energy_truth                 = *rdf.Histo1D({("miss_energy_truth_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_truth","event_weight");
                hist_miss_energy_truth.Write();
                TH1D hist_miss_mass_truth                   = *rdf.Histo1D({("miss_mass_truth_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_truth","event_weight");
                hist_miss_mass_truth.Write();
                TH1D hist_miss_masssquared_truth            = *rdf.Histo1D({("miss_masssquared_truth_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 0.1},"miss_masssquared_balance_truth","event_weight");
                hist_miss_masssquared_truth.Write();
                TH1D hist_miss_pminus_truth                 = *rdf.Histo1D({("miss_pminus_truth_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_truth","event_weight");
                hist_miss_pminus_truth.Write();
                TH1D hist_miss_momentum_truth               = *rdf.Histo1D({("miss_momentum_truth_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_truth","event_weight");
                hist_miss_momentum_truth.Write();
                TH2D hist_miss_momentum_pminus_truth        = *rdf.Histo2D({("miss_momentum_pminus_truth_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_truth","miss_pminus_truth","event_weight");
                hist_miss_momentum_pminus_truth.Write();
                TH2D hist_miss_momentum_energy_truth        = *rdf.Histo2D({("miss_momentum_energy_truth_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2} (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_truth","miss_energy_truth","event_weight");
                hist_miss_momentum_energy_truth.Write();

                TH1D hist_minust_truth                      = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 30, 0.0, 3.0},"minust_truth","event_weight");
                hist_minust_truth.Write();
                TH1D hist_coplanarity_truth                 = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_truth","event_weight");
                hist_coplanarity_truth.Write();
                TH1D hist_yphi_truth                        = *rdf.Histo1D({("yphi_truth_"+ label).c_str(), ";y_{#phi};Counts", 200, 0.0, 2.0},"y_phi_truth","event_weight");
                hist_yphi_truth.Write();
                TH1D hist_rho_mass_truth                    = *rdf.Histo1D({("rho_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 100, 0.0, 1.0},"rho_mass_truth","event_weight");
                hist_rho_mass_truth.Write();
                TH2D hist_beam_energy_minust_truth          = *rdf.Histo2D({("beam_energy_minust_truth_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_truth","minust_truth","event_weight");
                hist_beam_energy_minust_truth.Write();

                TH1D hist_scatter_theta_com_truth           = *rdf.Histo1D({("scatter_theta_com_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"scatter_theta_com_truth","event_weight");
                hist_scatter_theta_com_truth.Write();
                TH2D hist_minust_scatter_theta_com_truth    = *rdf.Histo2D({("minust_scatter_theta_com_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 30, 0.0, 3.0, 180, 0.0, 180.0},"minust_truth","scatter_theta_com_truth","event_weight");
                hist_minust_scatter_theta_com_truth.Write();
                TH1D hist_polarization_phi_com_truth        = *rdf.Histo1D({("polarization_phi_com_truth_"+ label).c_str(), ";#phi_{com} (deg);Counts", 9, -180, 180.0},"polarization_phi_com_truth","event_weight");
                hist_polarization_phi_com_truth.Write();

                TH1D hist_decay_costheta_helicity_truth     = *rdf.Histo1D({("decay_costheta_helicity_truth_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"decay_costheta_helicity_truth","event_weight");
                hist_decay_costheta_helicity_truth.Write();
                TH1D hist_decay_phi_helicity_truth          = *rdf.Histo1D({("decay_phi_helicity_truth_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"decay_phi_helicity_truth","event_weight");
                hist_decay_phi_helicity_truth.Write();
                TH1D hist_psi_helicity_truth                = *rdf.Histo1D({("psi_helicity_truth_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 9, -270.0, 270.0},"psi_helicity_truth","event_weight");
                hist_psi_helicity_truth.Write();

                TH2D hist_beam_energy_kin_truth             = *rdf.Histo2D({("beam_energy_kin_truth_"+ label).c_str(), ";E_{beam} (GeV);E_{beam} (GeV)", 60, 5.0, 11.0, 60, 5.0, 11.0},"beam_energy_kin","beam_energy_truth","event_weight");
                hist_beam_energy_kin_truth.Write();
                TH2D hist_minust_kin_truth                  = *rdf.Histo2D({("minust_kin_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});-t (GeV^{2}/c^{2})", 30, 0.0, 3.0, 30, 0.0, 3.0},"minust_kin","minust_truth","event_weight");
                hist_minust_kin_truth.Write();
            }
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}