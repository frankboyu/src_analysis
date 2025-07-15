#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;
double mass_missing = 0.0;

void filter_omega_d_recon_exc(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/test/selectedtree_omega_d_recon_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_omega_d_recon";
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
    .Define("beam_energy_meas",                 "beam_p4_meas.E()")
    .Define("beam_energy_kin",                  "beam_p4_kin.E()")
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",                 "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",                  "beam_x4_kin.T() - rftime")

    .Define("g1_energy_meas",                   "g1_p4_meas.E()")
    .Define("g1_energy_kin",                    "g1_p4_kin.E()")
    .Define("g1_energy_truth",                  "g1_p4_truth.E()")
    .Define("g1_theta_meas",                    "g1_p4_meas.Theta()*RadToDeg")
    .Define("g1_theta_kin",                     "g1_p4_kin.Theta()*RadToDeg")
    .Define("g1_theta_truth",                   "g1_p4_truth.Theta()*RadToDeg")
    .Define("g1_DeltaT_meas",                   "rftime + (g1_x4_meas.Z()-65.0)/29.9792458 - g1_x4_meas.T()")
    .Define("g1_DeltaT_kin",                    "rftime + (g1_x4_kin.Z()-65.0)/29.9792458 - g1_x4_kin.T()")

    .Define("g2_energy_meas",                   "g2_p4_meas.E()")
    .Define("g2_energy_kin",                    "g2_p4_kin.E()")
    .Define("g2_energy_truth",                  "g2_p4_truth.E()")
    .Define("g2_theta_meas",                    "g2_p4_meas.Theta()*RadToDeg")
    .Define("g2_theta_kin",                     "g2_p4_kin.Theta()*RadToDeg")
    .Define("g2_theta_truth",                   "g2_p4_truth.Theta()*RadToDeg")
    .Define("g2_DeltaT_meas",                   "rftime + (g2_x4_meas.Z()-65.0)/29.9792458 - g2_x4_meas.T()")
    .Define("g2_DeltaT_kin",                    "rftime + (g2_x4_kin.Z()-65.0)/29.9792458 - g2_x4_kin.T()")

    .Define("pip_as_kaon_p4_meas",              "TLorentzVector(pip_p4_meas.Vect(), TMath::Sqrt(pip_p4_meas.P()*pip_p4_meas.P() + mass_kplus*mass_kplus))")
    .Define("pip_as_kaon_p4_kin",               "TLorentzVector(pip_p4_kin.Vect(), TMath::Sqrt(pip_p4_kin.P()*pip_p4_kin.P() + mass_kplus*mass_kplus))")
    .Define("pip_as_kaon_p4_truth",             "TLorentzVector(pip_p4_truth.Vect(), TMath::Sqrt(pip_p4_truth.P()*pip_p4_truth.P() + mass_kplus*mass_kplus))")
    .Define("pip_energy_meas",                  "pip_p4_meas.E()")
    .Define("pip_energy_kin",                   "pip_p4_kin.E()")
    .Define("pip_energy_truth",                 "pip_p4_truth.E()")
    .Define("pip_momentum_meas",                "pip_p4_meas.P()")
    .Define("pip_momentum_kin",                 "pip_p4_kin.P()")
    .Define("pip_momentum_truth",               "pip_p4_truth.P()")
    .Define("pip_theta_meas",                   "pip_p4_meas.Theta()*RadToDeg")
    .Define("pip_theta_kin",                    "pip_p4_kin.Theta()*RadToDeg")
    .Define("pip_theta_truth",                  "pip_p4_truth.Theta()*RadToDeg")
    .Define("pip_DeltaT_meas",                  "rftime + (pip_x4_meas.Z()-65.0)/29.9792458 - pip_x4_meas.T()")
    .Define("pip_DeltaT_kin",                   "rftime + (pip_x4_kin.Z()-65.0)/29.9792458 - pip_x4_kin.T()")
    .Define("pip_in_fdc_meas",                  "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_in_cdc_meas",                  "accidweight*(pip_dedx_cdc > 0.0 && pip_dedx_fdc == 0.0)")
    .Define("pip_in_fdc_cdc_meas",              "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc > 0.0)")
    .Define("pip_in_neither_meas",              "accidweight*(pip_dedx_fdc == 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_dedx_fdc_keV_per_cm_meas",     "pip_dedx_fdc*1e6")
    .Define("pip_dedx_cdc_keV_per_cm_meas",     "pip_dedx_cdc*1e6")
    .Define("pip_dedx_st_keV_per_cm_meas",      "pip_dedx_st*1e6")
    .Define("pip_dedx_tof_keV_per_cm_meas",     "pip_dedx_tof*1e6")

    .Define("pim_as_kaon_p4_meas",              "TLorentzVector(pim_p4_meas.Vect(), TMath::Sqrt(pim_p4_meas.P()*pim_p4_meas.P() + mass_kminus*mass_kminus))")
    .Define("pim_as_kaon_p4_kin",               "TLorentzVector(pim_p4_kin.Vect(), TMath::Sqrt(pim_p4_kin.P()*pim_p4_kin.P() + mass_kminus*mass_kminus))")
    .Define("pim_as_kaon_p4_truth",             "TLorentzVector(pim_p4_truth.Vect(), TMath::Sqrt(pim_p4_truth.P()*pim_p4_truth.P() + mass_kminus*mass_kminus))")
    .Define("pim_energy_meas",                  "pim_p4_meas.E()")
    .Define("pim_energy_kin",                   "pim_p4_kin.E()")
    .Define("pim_energy_truth",                 "pim_p4_truth.E()")
    .Define("pim_momentum_meas",                "pim_p4_meas.P()")
    .Define("pim_momentum_kin",                 "pim_p4_kin.P()")
    .Define("pim_momentum_truth",               "pim_p4_truth.P()")
    .Define("pim_theta_meas",                   "pim_p4_meas.Theta()*RadToDeg")
    .Define("pim_theta_kin",                    "pim_p4_kin.Theta()*RadToDeg")
    .Define("pim_theta_truth",                  "pim_p4_truth.Theta()*RadToDeg")
    .Define("pim_DeltaT_meas",                  "rftime + (pim_x4_meas.Z()-65.0)/29.9792458 - pim_x4_meas.T()")
    .Define("pim_DeltaT_kin",                   "rftime + (pim_x4_kin.Z()-65.0)/29.9792458 - pim_x4_kin.T()")
    .Define("pim_in_fdc_meas",                  "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc == 0.0)")
    .Define("pim_in_cdc_meas",                  "accidweight*(pim_dedx_cdc > 0.0 && pim_dedx_fdc == 0.0)")
    .Define("pim_in_fdc_cdc_meas",              "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc > 0.0)")
    .Define("pim_in_neither_meas",              "accidweight*(pim_dedx_fdc == 0.0 && pim_dedx_cdc == 0.0)")
    .Define("pim_dedx_fdc_keV_per_cm_meas",     "pim_dedx_fdc*1e6")
    .Define("pim_dedx_cdc_keV_per_cm_meas",     "pim_dedx_cdc*1e6")
    .Define("pim_dedx_st_keV_per_cm_meas",      "pim_dedx_st*1e6")
    .Define("pim_dedx_tof_keV_per_cm_meas",     "pim_dedx_tof*1e6")

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
    .Define("d_dedx_st_keV_per_cm_meas",        "d_dedx_st*1e6")
    .Define("d_dedx_tof_keV_per_cm_meas",       "d_dedx_tof*1e6")

    .Define("decaypi0_p4_meas",                 "g1_p4_meas + g2_p4_meas")
    .Define("decaypi0_energy_meas",             "decaypi0_p4_meas.E()")
    .Define("decaypi0_energy_kin",              "decaypi0_p4_kin.E()")
    .Define("decaypi0_energy_truth",            "decaypi0_p4_truth.E()")
    .Define("decaypi0_momentum_meas",           "decaypi0_p4_meas.P()")
    .Define("decaypi0_momentum_kin",            "decaypi0_p4_kin.P()")
    .Define("decaypi0_momentum_truth",          "decaypi0_p4_truth.P()")
    .Define("decaypi0_theta_meas",              "decaypi0_p4_meas.Theta()*RadToDeg")
    .Define("decaypi0_theta_kin",               "decaypi0_p4_kin.Theta()*RadToDeg")
    .Define("decaypi0_theta_truth",             "decaypi0_p4_truth.Theta()*RadToDeg")

    .Define("omega_p4_meas",                    "decaypi0_p4_meas + pip_p4_meas + pim_p4_meas")
    .Define("omega_p4_kin",                     "decaypi0_p4_kin + pip_p4_kin + pim_p4_kin")
    .Define("omega_p4_truth",                   "decaypi0_p4_truth + pip_p4_truth + pim_p4_truth")
    .Define("omega_energy_meas",                "omega_p4_meas.E()")
    .Define("omega_energy_kin",                 "omega_p4_kin.E()")
    .Define("omega_energy_truth",               "omega_p4_truth.E()")
    .Define("omega_momentum_meas",              "omega_p4_meas.P()")
    .Define("omega_momentum_kin",               "omega_p4_kin.P()")
    .Define("omega_momentum_truth",             "omega_p4_truth.P()")
    .Define("omega_mass_meas",                  "omega_p4_meas.M()")
    .Define("omega_mass_kin",                   "omega_p4_kin.M()")
    .Define("omega_mass_truth",                 "omega_p4_truth.M()")
    .Define("omega_proxymass_meas",             "TMath::Sqrt(omega_p4_meas.Minus()*(2*beam_energy_meas + mass_target - d_p4_meas.Plus() - (TMath::Sq(mass_missing) + (omega_p4_meas + d_p4_meas).Perp2())/(mass_target - (omega_p4_meas + d_p4_meas).Minus())) - omega_p4_meas.Perp2())")
    .Define("omega_proxymass_kin",              "TMath::Sqrt(omega_p4_kin.Minus()*(2*beam_energy_kin + mass_target - d_p4_kin.Plus() - (TMath::Sq(mass_missing) + (omega_p4_kin + d_p4_kin).Perp2())/(mass_target - (omega_p4_kin + d_p4_kin).Minus())) - omega_p4_kin.Perp2())")
    .Define("omega_proxymass_truth",            "TMath::Sqrt(omega_p4_truth.Minus()*(2*beam_energy_truth + mass_target - d_p4_truth.Plus() - (TMath::Sq(mass_missing) + (omega_p4_truth + d_p4_truth).Perp2())/(mass_target - (omega_p4_truth + d_p4_truth).Minus())) - omega_p4_truth.Perp2())")
    .Define("omega_theta_meas",                 "omega_p4_meas.Theta()*RadToDeg")
    .Define("omega_theta_kin",                  "omega_p4_kin.Theta()*RadToDeg")
    .Define("omega_theta_truth",                "omega_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_meas",                   "omega_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin",                    "omega_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("struck_p4_truth",                  "omega_p4_truth + d_p4_truth - beam_p4_truth")
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

    .Define("miss_p4_meas",                     "beam_p4_meas + TLorentzVector(0, 0, 0, mass_target) - omega_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                      "beam_p4_kin + TLorentzVector(0, 0, 0, mass_target) - omega_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                    "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target) - omega_p4_truth - d_p4_truth")
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
    .Define("total_final_p4_meas",              "omega_p4_meas + d_p4_meas")
    .Define("total_final_p4_kin",               "omega_p4_kin + d_p4_kin")
    .Define("total_final_p4_truth",             "omega_p4_truth + d_p4_truth")
    .Define("sqrts_meas",                       "total_final_p4_meas.Mag()")
    .Define("sqrts_kin",                        "total_final_p4_kin.Mag()")
    .Define("sqrts_truth",                      "total_final_p4_truth.Mag()")
    .Define("minust_meas",                      "-(beam_p4_meas - omega_p4_meas).Mag2()")
    .Define("minust_kin",                       "-(beam_p4_kin - omega_p4_kin).Mag2()")
    .Define("minust_truth",                     "-(beam_p4_truth - omega_p4_truth).Mag2()")
    .Define("minusu_meas",                      "-(beam_p4_meas - d_p4_meas).Mag2()")
    .Define("minusu_kin",                       "-(beam_p4_kin - d_p4_kin).Mag2()")
    .Define("minusu_truth",                     "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("coplanarity_meas",                 "abs(omega_p4_meas.Phi() - d_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin",                  "abs(omega_p4_kin.Phi() - d_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth",                "abs(omega_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("y_omega_meas",                     "minust_meas/(2*mass_2H*(beam_p4_meas.E()-omega_p4_meas.E()))")
    .Define("y_omega_kin",                      "minust_kin/(2*mass_2H*(beam_p4_kin.E()-omega_p4_kin.E()))")
    .Define("y_omega_truth",                    "minust_truth/(2*mass_2H*(beam_p4_truth.E()-omega_p4_truth.E()))")
    .Define("phi_mass_meas",                    "(pip_as_kaon_p4_meas + pim_as_kaon_p4_meas).M()")
    .Define("phi_mass_kin",                     "(pip_as_kaon_p4_kin + pim_as_kaon_p4_kin).M()")
    .Define("phi_mass_truth",                   "(pip_as_kaon_p4_truth + pim_as_kaon_p4_truth).M()")
    .Define("kinfit_fom_kin",                   "TMath::Prob(kin_chisq,kin_ndf)")

    .Define("pippim_p4_meas",                   "pip_p4_meas + pim_p4_meas")
    .Define("pippim_p4_kin",                    "pip_p4_kin + pim_p4_kin")
    .Define("pippim_p4_truth",                  "pip_p4_truth + pim_p4_truth")
    .Define("pimpi0_p4_meas",                   "pim_p4_meas + decaypi0_p4_meas")
    .Define("pimpi0_p4_kin",                    "pim_p4_kin + decaypi0_p4_kin")
    .Define("pimpi0_p4_truth",                  "pim_p4_truth + decaypi0_p4_truth")
    .Define("dalitz_s12_meas",                  "pippim_p4_meas.M2()")
    .Define("dalitz_s12_kin",                   "pippim_p4_kin.M2()")
    .Define("dalitz_s12_truth",                 "pippim_p4_truth.M2()")
    .Define("dalitz_s23_meas",                  "pimpi0_p4_meas.M2()")
    .Define("dalitz_s23_kin",                   "pimpi0_p4_kin.M2()")
    .Define("dalitz_s23_truth",                 "pimpi0_p4_truth.M2()")

    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_meas",                 "boost_lorentz_vector(beam_p4_meas, -total_final_p4_meas.BoostVector())")
    .Define("beam_p4_com_kin",                  "boost_lorentz_vector(beam_p4_kin, -total_final_p4_kin.BoostVector())")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("omega_p4_com_meas",                "boost_lorentz_vector(omega_p4_meas, -total_final_p4_meas.BoostVector())")
    .Define("omega_p4_com_kin",                 "boost_lorentz_vector(omega_p4_kin, -total_final_p4_kin.BoostVector())")
    .Define("omega_p4_com_truth",               "boost_lorentz_vector(omega_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("pip_p4_com_meas",                  "boost_lorentz_vector(pip_p4_meas, -total_final_p4_meas.BoostVector())")
    .Define("pip_p4_com_kin",                   "boost_lorentz_vector(pip_p4_kin, -total_final_p4_kin.BoostVector())")
    .Define("pip_p4_com_truth",                 "boost_lorentz_vector(pip_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("pim_p4_com_meas",                  "boost_lorentz_vector(pim_p4_meas, -total_final_p4_meas.BoostVector())")
    .Define("pim_p4_com_kin",                   "boost_lorentz_vector(pim_p4_kin, -total_final_p4_kin.BoostVector())")
    .Define("pim_p4_com_truth",                 "boost_lorentz_vector(pim_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("z_x3_com_meas",                    "beam_p4_com_meas.Vect().Unit()")
    .Define("z_x3_com_kin",                     "beam_p4_com_kin.Vect().Unit()")
    .Define("z_x3_com_truth",                   "beam_p4_com_truth.Vect().Unit()")
    .Define("y_x3_com_meas",                    "beam_p4_com_meas.Vect().Cross(omega_p4_com_meas.Vect()).Unit()")
    .Define("y_x3_com_kin",                     "beam_p4_com_kin.Vect().Cross(omega_p4_com_kin.Vect()).Unit()")
    .Define("y_x3_com_truth",                   "beam_p4_com_truth.Vect().Cross(omega_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_com_meas",                    "y_x3_com_meas.Cross(z_x3_com_meas).Unit()")
    .Define("x_x3_com_kin",                     "y_x3_com_kin.Cross(z_x3_com_kin).Unit()")
    .Define("x_x3_com_truth",                   "y_x3_com_truth.Cross(z_x3_com_truth).Unit()")
    .Define("polarization_phi_com_meas",        "TMath::ATan2(-x_x3_com_meas.Dot(epsilon_x3_com.Cross(z_x3_com_meas)), y_x3_com_meas.Dot(epsilon_x3_com.Cross(z_x3_com_meas)))*RadToDeg")
    .Define("polarization_phi_com_kin",         "TMath::ATan2(-x_x3_com_kin.Dot(epsilon_x3_com.Cross(z_x3_com_kin)), y_x3_com_kin.Dot(epsilon_x3_com.Cross(z_x3_com_kin)))*RadToDeg")
    .Define("polarization_phi_com_truth",       "TMath::ATan2(-x_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)), y_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)))*RadToDeg")
    .Define("scatter_theta_com_meas",           "omega_p4_com_meas.Vect().Angle(z_x3_com_meas)*RadToDeg")
    .Define("scatter_theta_com_kin",            "omega_p4_com_kin.Vect().Angle(z_x3_com_kin)*RadToDeg")
    .Define("scatter_theta_com_truth",          "omega_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")

    .Define("pip_p4_helicity_meas",             "boost_lorentz_vector(pip_p4_com_meas, -omega_p4_com_meas.BoostVector())")
    .Define("pip_p4_helicity_kin",              "boost_lorentz_vector(pip_p4_com_kin, -omega_p4_com_kin.BoostVector())")
    .Define("pip_p4_helicity_truth",            "boost_lorentz_vector(pip_p4_com_truth, -omega_p4_com_truth.BoostVector())")
    .Define("pim_p4_helicity_meas",             "boost_lorentz_vector(pim_p4_com_meas, -omega_p4_com_meas.BoostVector())")
    .Define("pim_p4_helicity_kin",              "boost_lorentz_vector(pim_p4_com_kin, -omega_p4_com_kin.BoostVector())")
    .Define("pim_p4_helicity_truth",            "boost_lorentz_vector(pim_p4_com_truth, -omega_p4_com_truth.BoostVector())")
    .Define("z_x3_helicity_meas",               "omega_p4_com_meas.Vect().Unit()")
    .Define("z_x3_helicity_kin",                "omega_p4_com_kin.Vect().Unit()")
    .Define("z_x3_helicity_truth",              "omega_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_meas",               "beam_p4_com_meas.Vect().Cross(omega_p4_com_meas.Vect()).Unit()")
    .Define("y_x3_helicity_kin",                "beam_p4_com_kin.Vect().Cross(omega_p4_com_kin.Vect()).Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(omega_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_meas",               "y_x3_helicity_meas.Cross(z_x3_helicity_meas).Unit()")
    .Define("x_x3_helicity_kin",                "y_x3_helicity_kin.Cross(z_x3_helicity_kin).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_meas",              "pip_p4_helicity_meas.Vect().Cross(pim_p4_helicity_meas.Vect()).Unit()")
    .Define("pi_x3_helicity_kin",               "pip_p4_helicity_kin.Vect().Cross(pim_p4_helicity_kin.Vect()).Unit()")
    .Define("pi_x3_helicity_truth",             "pip_p4_helicity_truth.Vect().Cross(pim_p4_helicity_truth.Vect()).Unit()")
    .Define("pi_x3_helicity_truth",             "boost_lorentz_vector(pip_p4_truth, -omega_p4_truth.BoostVector()).Vect().Unit()")
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
    auto rdf_KinFitFOMCut   = rdf_dEdxCut.Filter("kinfit_fom_kin > 0.01");
    auto rdf_PIDFOMCut      = rdf_KinFitFOMCut.Filter("(pip_pidfom > 0.01) && (pim_pidfom > 0.01)");
    auto rdf_MissPCut       = rdf_PIDFOMCut.Filter("(abs(struck_energy_balance_kin) < 1.0)");
    auto rdf_OmegaMassCut   = rdf_MissPCut.Filter("omega_mass_kin > 0.0");
    auto rdf_output         = rdf_OmegaMassCut;
    RNode rdfs []           = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_MissPCut,   rdf_OmegaMassCut};
    string labels []        = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "MissPCut",     "OmegaMassCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredtree_omega_d_recon_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_omega_d_recon";
        rdf_output.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredhist_omega_d_recon_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

            TH1D hist_beam_energy_kin                       = *rdf.Histo1D({("beam_energy_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_kin","accidweight");
            hist_beam_energy_kin.Write();
            TH1D hist_beam_DeltaT_kin                       = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 100, -20.0, 20.0},"beam_DeltaT_kin","accidweight");
            hist_beam_DeltaT_kin.Write();

            TH1D hist_pip_pidfom                            = *rdf.Histo1D({("pip_pidfom_"+ label).c_str(), ";pip_pidfom;Counts", 100, 0.0, 1.0},"pip_pidfom","accidweight");
            hist_pip_pidfom.Write();
            TH1D hist_pip_DeltaT_kin                        = *rdf.Histo1D({("pip_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{+}} (ns);Counts", 100, -2.0, 2.0},"pip_DeltaT_kin","accidweight");
            hist_pip_DeltaT_kin.Write();
            TH2D hist_pip_dEdx_cdc_kin                      = *rdf.Histo2D({("pip_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pip_momentum_kin","pip_dedx_cdc_keV_per_cm_meas","accidweight");
            hist_pip_dEdx_cdc_kin.Write();
            TH2D hist_pip_dEdx_fdc_kin                      = *rdf.Histo2D({("pip_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pip_momentum_kin","pip_dedx_fdc_keV_per_cm_meas","accidweight");
            hist_pip_dEdx_fdc_kin.Write();
            TH2D hist_pip_dEdx_tof_kin                      = *rdf.Histo2D({("pip_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pip_momentum_kin","pip_dedx_tof_keV_per_cm_meas","accidweight");
            hist_pip_dEdx_tof_kin.Write();
            TH2D hist_pip_dEdx_st_kin                       = *rdf.Histo2D({("pip_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pip_momentum_kin","pip_dedx_st_keV_per_cm_meas","accidweight");
            hist_pip_dEdx_st_kin.Write();
            TH2D hist_pip_kinematics_kin                    = *rdf.Histo2D({("pip_kinematics_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","accidweight");
            hist_pip_kinematics_kin.Write();
            TH2D hist_pip_kinematics_fdc_kin                = *rdf.Histo2D({("pip_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc_meas");
            hist_pip_kinematics_fdc_kin.Write();
            TH2D hist_pip_kinematics_fdc_cdc_kin            = *rdf.Histo2D({("pip_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc_cdc_meas");
            hist_pip_kinematics_fdc_cdc_kin.Write();
            TH2D hist_pip_kinematics_cdc_kin                = *rdf.Histo2D({("pip_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_cdc_meas");
            hist_pip_kinematics_cdc_kin.Write();
            TH2D hist_pip_kinematics_neither_kin            = *rdf.Histo2D({("pip_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_neither_meas");
            hist_pip_kinematics_neither_kin.Write();

            TH1D hist_pim_pidfom                            = *rdf.Histo1D({("pim_pidfom_"+ label).c_str(), ";pim_pidfom;Counts", 100, 0.0, 1.0},"pim_pidfom","accidweight");
            hist_pim_pidfom.Write();
            TH1D hist_pim_DeltaT_kin                        = *rdf.Histo1D({("pim_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{-}} (ns);Counts", 100, -2.0, 2.0},"pim_DeltaT_kin","accidweight");
            hist_pim_DeltaT_kin.Write();
            TH2D hist_pim_dEdx_cdc_kin                      = *rdf.Histo2D({("pim_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pim_momentum_kin","pim_dedx_cdc_keV_per_cm_meas","accidweight");
            hist_pim_dEdx_cdc_kin.Write();
            TH2D hist_pim_dEdx_fdc_kin                      = *rdf.Histo2D({("pim_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pim_momentum_kin","pim_dedx_fdc_keV_per_cm_meas","accidweight");
            hist_pim_dEdx_fdc_kin.Write();
            TH2D hist_pim_dEdx_tof_kin                      = *rdf.Histo2D({("pim_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pim_momentum_kin","pim_dedx_tof_keV_per_cm_meas","accidweight");
            hist_pim_dEdx_tof_kin.Write();
            TH2D hist_pim_dEdx_st_kin                       = *rdf.Histo2D({("pim_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"pim_momentum_kin","pim_dedx_st_keV_per_cm_meas","accidweight");
            hist_pim_dEdx_st_kin.Write();
            TH2D hist_pim_kinematics_kin                    = *rdf.Histo2D({("pim_kinematics_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","accidweight");
            hist_pim_kinematics_kin.Write();
            TH2D hist_pim_kinematics_fdc_kin                = *rdf.Histo2D({("pim_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_fdc_meas");
            hist_pim_kinematics_fdc_kin.Write();
            TH2D hist_pim_kinematics_fdc_cdc_kin            = *rdf.Histo2D({("pim_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_fdc_cdc_meas");
            hist_pim_kinematics_fdc_cdc_kin.Write();
            TH2D hist_pim_kinematics_cdc_kin                = *rdf.Histo2D({("pim_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_cdc_meas");
            hist_pim_kinematics_cdc_kin.Write();
            TH2D hist_pim_kinematics_neither_kin            = *rdf.Histo2D({("pim_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_neither_meas");
            hist_pim_kinematics_neither_kin.Write();

            TH1D hist_d_DeltaT_kin                          = *rdf.Histo1D({("d_DeltaT_kin_"+ label).c_str(), ";#Delta t_{d} (ns);Counts", 100, -5.0, 5.0},"d_DeltaT_kin","accidweight");
            hist_d_DeltaT_kin.Write();
            TH2D hist_d_dEdx_cdc_meas                       = *rdf.Histo2D({("d_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_cdc_keV_per_cm_meas","accidweight");
            hist_d_dEdx_cdc_meas.Write();
            TH2D hist_d_dEdx_cdc_kin                        = *rdf.Histo2D({("d_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_kin","d_dedx_cdc_keV_per_cm_meas","accidweight");
            hist_d_dEdx_cdc_kin.Write();
            TH2D hist_d_dEdx_tof_kin                        = *rdf.Histo2D({("d_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_kin","d_dedx_tof_keV_per_cm_meas","accidweight");
            hist_d_dEdx_tof_kin.Write();
            TH2D hist_d_dEdx_st_kin                         = *rdf.Histo2D({("d_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"d_momentum_kin","d_dedx_st_keV_per_cm_meas","accidweight");
            hist_d_dEdx_st_kin.Write();
            TH2D hist_d_kinematics_kin                      = *rdf.Histo2D({("d_kinematics_kin_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_kin","d_theta_kin","accidweight");
            hist_d_kinematics_kin.Write();

            TH1D hist_omega_mass_kin                        = *rdf.Histo1D({("omega_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c);Counts", 600, 0.0, 3.0},"omega_mass_kin","accidweight");
            hist_omega_mass_kin.Write();
            TH1D hist_omega_mass_meas                       = *rdf.Histo1D({("omega_mass_meas_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c);Counts", 600, 0.0, 3.0},"omega_mass_meas","accidweight");
            hist_omega_mass_meas.Write();
            TH1D hist_omega_proxymass_kin                   = *rdf.Histo1D({("omega_proxymass_kin_"+ label).c_str(), ";m^{proxy}_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 600, 0.0, 3.0},"omega_proxymass_kin","accidweight");
            hist_omega_proxymass_kin.Write();
            TH1D hist_omega_proxymass_meas                  = *rdf.Histo1D({("omega_proxymass_meas_"+ label).c_str(), ";m^{proxy}_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 600, 0.0, 3.0},"omega_proxymass_meas","accidweight");
            hist_omega_proxymass_meas.Write();
            TH2D hist_omega_mass_kin_meas                   = *rdf.Histo2D({("omega_mass_kin_meas_"+ label).c_str(), ";m^{kin}_{#pi^{+}#pi^{-}} (GeV/c);m^{meas}_{#pi^{+}#pi^{-}} (GeV/c^{2})", 600, 0.0, 3.0, 600, 0.0, 3.0},"omega_mass_kin","omega_mass_meas","accidweight");
            hist_omega_mass_kin_meas.Write();
            TH2D hist_omega_kinematics_kin                  = *rdf.Histo2D({("omega_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 110, 0.0, 11.0, 180, 0.0, 180.0},"omega_momentum_kin","omega_theta_kin","accidweight");
            hist_omega_kinematics_kin.Write();

            TH1D hist_miss_energy_kin                       = *rdf.Histo1D({("miss_energy_kin_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_kin","accidweight");
            hist_miss_energy_kin.Write();
            TH1D hist_miss_energy_meas                      = *rdf.Histo1D({("miss_energy_meas_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_meas","accidweight");
            hist_miss_energy_meas.Write();
            TH1D hist_miss_mass_kin                         = *rdf.Histo1D({("miss_mass_kin_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_kin","accidweight");
            hist_miss_mass_kin.Write();
            TH1D hist_miss_mass_meas                        = *rdf.Histo1D({("miss_mass_meas_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_meas","accidweight");
            hist_miss_mass_meas.Write();
            TH1D hist_miss_masssquared_kin                  = *rdf.Histo1D({("miss_masssquared_kin_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 0.1},"miss_masssquared_balance_kin","accidweight");
            hist_miss_masssquared_kin.Write();
            TH1D hist_miss_masssquared_meas                 = *rdf.Histo1D({("miss_masssquared_meas_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 0.1},"miss_masssquared_balance_meas","accidweight");
            hist_miss_masssquared_meas.Write();
            TH1D hist_miss_pminus_kin                       = *rdf.Histo1D({("miss_pminus_kin_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_balance_kin","accidweight");
            hist_miss_pminus_kin.Write();
            TH1D hist_miss_pminus_meas                      = *rdf.Histo1D({("miss_pminus_meas_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_balance_meas","accidweight");
            hist_miss_pminus_meas.Write();
            TH1D hist_miss_momentum_kin                     = *rdf.Histo1D({("miss_momentum_kin_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_kin","accidweight");
            hist_miss_momentum_kin.Write();
            TH1D hist_miss_momentum_meas                    = *rdf.Histo1D({("miss_momentum_meas_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_meas","accidweight");
            hist_miss_momentum_meas.Write();
            TH2D hist_miss_momentum_pminus_kin              = *rdf.Histo2D({("miss_momentum_pminus_kin_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_kin","miss_pminus_balance_kin","accidweight");
            hist_miss_momentum_pminus_kin.Write();
            TH2D hist_miss_momentum_pminus_meas             = *rdf.Histo2D({("miss_momentum_pminus_meas_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_meas","miss_pminus_balance_meas","accidweight");
            hist_miss_momentum_pminus_meas.Write();
            TH2D hist_miss_momentum_energy_kin              = *rdf.Histo2D({("miss_momentum_energy_kin_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2}  (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_kin","miss_energy_balance_kin","accidweight");
            hist_miss_momentum_energy_kin.Write();
            TH2D hist_miss_momentum_energy_meas             = *rdf.Histo2D({("miss_momentum_energy_meas_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2} (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_meas","miss_energy_balance_meas","accidweight");
            hist_miss_momentum_energy_meas.Write();
            TH2D hist_miss_energy_omega_mass_meas           = *rdf.Histo2D({("miss_energy_omega_mass_meas_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);m_{#pi^{+}#pi^{-}} (GeV/c)", 300, -3.0, 3.0, 600, 0.0, 3.0},"miss_energy_balance_meas","omega_mass_meas","accidweight");
            hist_miss_energy_omega_mass_meas.Write();

            TH1D hist_sqrts_kin                             = *rdf.Histo1D({("sqrts_kin_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 40, 4.0, 8.0},"sqrts_kin","accidweight");
            hist_sqrts_kin.Write();
            TH1D hist_minust_kin                            = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 30, 0.0, 3.0},"minust_kin","accidweight");
            hist_minust_kin.Write();
            TH1D hist_minusu_kin                            = *rdf.Histo1D({("minusu_kin_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 300, 15.0, 45.0},"minusu_kin","accidweight");
            hist_minusu_kin.Write();
            TH1D hist_coplanarity_kin                       = *rdf.Histo1D({("coplanarity_kin_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_kin","accidweight");
            hist_coplanarity_kin.Write();
            TH1D hist_coplanarity_meas                      = *rdf.Histo1D({("coplanarity_meas_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_meas","accidweight");
            hist_coplanarity_meas.Write();
            TH1D hist_yomega_kin                            = *rdf.Histo1D({("yomega_kin_"+ label).c_str(), ";y_{#omega};Counts", 200, 0.0, 2.0},"y_omega_kin","accidweight");
            hist_yomega_kin.Write();
            TH1D hist_yomega_meas                           = *rdf.Histo1D({("yomega_meas_"+ label).c_str(), ";y_{#omega};Counts", 200, 0.0, 2.0},"y_omega_meas","accidweight");
            hist_yomega_meas.Write();
            TH1D hist_phi_mass_kin                          = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c^{2});Counts", 100, 1.0, 2.0},"phi_mass_kin","accidweight");
            hist_phi_mass_kin.Write();
            TH1D hist_kinfit_fom_kin                        = *rdf.Histo1D({("kinfit_fom_kin_"+ label).c_str(), ";KinFit FOM;Counts", 100, 0.0, 1.0},"kinfit_fom_kin","accidweight");
            hist_kinfit_fom_kin.Write();
            TH2D hist_beam_energy_minust_kin                = *rdf.Histo2D({("beam_energy_minust_kin_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_kin","minust_kin","accidweight");
            hist_beam_energy_minust_kin.Write();

            TH2D hist_dalitz_s12_s23_kin                    = *rdf.Histo2D({("dalitz_s12_s23_kin_"+ label).c_str(), ";s(#pi^{+}#pi^{-}) (GeV^2/c);s(#pi^{-}#pi^{0}) (GeV^2/c)", 100, 0.0, 5.0, 100, 0.0, 5.0},"dalitz_s12_kin","dalitz_s23_kin","accidweight");
            hist_dalitz_s12_s23_kin.Write();

            TH1D hist_scatter_theta_com_kin                 = *rdf.Histo1D({("scatter_theta_com_kin_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"scatter_theta_com_kin","accidweight");
            hist_scatter_theta_com_kin.Write();
            TH2D hist_minust_scatter_theta_com_kin          = *rdf.Histo2D({("minust_scatter_theta_com_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 20, 0.0, 2.0, 40, 0.0, 40.0},"minust_kin","scatter_theta_com_kin","accidweight");
            hist_minust_scatter_theta_com_kin.Write();
            TH1D hist_polarization_phi_com_kin              = *rdf.Histo1D({("polarization_phi_com_kin_"+ label).c_str(), ";#phi_{com} (deg);Counts", 9, -180, 180.0},"polarization_phi_com_kin","accidweight");
            hist_polarization_phi_com_kin.Write();

            TH1D hist_decay_costheta_helicity_kin           = *rdf.Histo1D({("decay_costheta_helicity_kin_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"decay_costheta_helicity_kin","accidweight");
            hist_decay_costheta_helicity_kin.Write();
            TH1D hist_decay_phi_helicity_kin                = *rdf.Histo1D({("decay_phi_helicity_kin_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"decay_phi_helicity_kin","accidweight");
            hist_decay_phi_helicity_kin.Write();
            TH1D hist_psi_helicity_kin                      = *rdf.Histo1D({("psi_helicity_kin_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 18, -180.0, 180.0},"psi_helicity_kin","accidweight");
            hist_psi_helicity_kin.Write();


            if (reaction.find("sim") != string::npos)
            {
                TH1D hist_beam_energy_truth                 = *rdf.Histo1D({("beam_energy_truth_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_truth");
                hist_beam_energy_truth.Write();
                TH2D hist_pip_kinematics_truth              = *rdf.Histo2D({("pip_kinematics_truth_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_truth","pip_theta_truth");
                hist_pip_kinematics_truth.Write();
                TH2D hist_pim_kinematics_truth              = *rdf.Histo2D({("pim_kinematics_truth_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_truth","pim_theta_truth");
                hist_pim_kinematics_truth.Write();
                TH2D hist_d_kinematics_truth                = *rdf.Histo2D({("d_kinematics_truth_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_truth","d_theta_truth");
                hist_d_kinematics_truth.Write();
                TH1D hist_omega_mass_truth                  = *rdf.Histo1D({("omega_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c);Counts", 600, 0.0, 3.0},"omega_mass_truth");
                hist_omega_mass_truth.Write();
                TH2D hist_omega_kinematics_truth            = *rdf.Histo2D({("omega_kinematics_truth_"+ label).c_str(), ";P_{#omega} (GeV/c);#theta_{#omega} (deg)", 100, 0.0, 11.0, 180, 0.0, 180.0},"omega_momentum_truth","scatter_theta_truth");
                hist_omega_kinematics_truth.Write();
                TH2D hist_omega_d_theta_truth               = *rdf.Histo2D({("omega_d_theta_truth_"+ label).c_str(), ";#theta_{d} (deg);#theta_{#omega} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"d_theta_truth","scatter_theta_truth");
                hist_omega_d_theta_truth.Write();
                TH2D hist_omega_d_momentum_truth            = *rdf.Histo2D({("omega_d_momentum_truth_"+ label).c_str(), ";P_{d} (GeV/c);P_{#omega} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"d_momentum_truth","omega_momentum_truth");
                hist_omega_d_momentum_truth.Write();

                TH1D hist_miss_energy_truth                 = *rdf.Histo1D({("miss_energy_truth_"+ label).c_str(), ";E_{miss} - m_{A-2} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_balance_truth");
                hist_miss_energy_truth.Write();
                TH1D hist_miss_mass_truth                   = *rdf.Histo1D({("miss_mass_truth_"+ label).c_str(), ";m_{miss} - m_{A-2} (GeV/c^{2});Counts", 100, -0.5, 0.5},"miss_mass_balance_truth");
                hist_miss_mass_truth.Write();
                TH1D hist_miss_masssquared_truth            = *rdf.Histo1D({("miss_masssquared_truth_"+ label).c_str(), ";m_{miss}^{2} - m_{A-2}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 0.1},"miss_masssquared_balance_truth");
                hist_miss_masssquared_truth.Write();
                TH1D hist_miss_pminus_truth                 = *rdf.Histo1D({("miss_pminus_truth_"+ label).c_str(), ";P_{miss}^{-} - m_{A-2} (GeV/c);Counts", 100, -1.0, 1.0},"miss_pminus_truth");
                hist_miss_pminus_truth.Write();
                TH1D hist_miss_momentum_truth               = *rdf.Histo1D({("miss_momentum_truth_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_truth");
                hist_miss_momentum_truth.Write();
                TH2D hist_miss_momentum_pminus_truth        = *rdf.Histo2D({("miss_momentum_pminus_truth_"+ label).c_str(), ";P_{miss} (GeV/c);P_{miss}^{-} - m_{A-2} (GeV/c)", 200, 0.0, 2.0, 100, -1.0, 1.0},"miss_momentum_truth","miss_pminus_truth");
                hist_miss_momentum_pminus_truth.Write();
                TH2D hist_miss_momentum_energy_truth        = *rdf.Histo2D({("miss_momentum_energy_truth_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{A-2} (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_truth","miss_energy_truth");
                hist_miss_momentum_energy_truth.Write();

                TH1D hist_sqrts_truth                       = *rdf.Histo1D({("sqrts_truth_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 40, 4.0, 8.0},"sqrts_truth");
                hist_sqrts_truth.Write();
                TH1D hist_minust_truth                      = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 30, 0.0, 3.0},"minust_truth");
                hist_minust_truth.Write();
                TH1D hist_minusu_truth                      = *rdf.Histo1D({("minusu_truth_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 300, 15.0, 45.0},"minusu_truth");
                hist_minusu_truth.Write();
                TH1D hist_coplanarity_truth                 = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_truth");
                hist_coplanarity_truth.Write();
                TH1D hist_yomega_truth                      = *rdf.Histo1D({("yomega_truth_"+ label).c_str(), ";y_{#omega};Counts", 200, 0.0, 2.0},"y_omega_truth");
                hist_yomega_truth.Write();
                TH1D hist_phi_mass_truth                    = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c^{2});Counts", 100, 1.0, 2.0},"phi_mass_truth");
                hist_omega_mass_truth.Write();
                TH2D hist_beam_energy_minust_truth          = *rdf.Histo2D({("beam_energy_minust_truth_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_truth","minust_truth");
                hist_beam_energy_minust_truth.Write();

                TH1D hist_scatter_theta_com_truth           = *rdf.Histo1D({("scatter_theta_com_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"scatter_theta_com_truth");
                hist_scatter_theta_com_truth.Write();
                TH2D hist_minust_scatter_theta_com_truth    = *rdf.Histo2D({("minust_scatter_theta_com_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 30, 0.0, 3.0, 180, 0.0, 180.0},"minust_truth","scatter_theta_com_truth");
                hist_minust_scatter_theta_com_truth.Write();
                TH1D hist_polarization_phi_com_truth        = *rdf.Histo1D({("polarization_phi_com_truth_"+ label).c_str(), ";#omega_{com} (deg);Counts", 9, -180, 180.0},"polarization_phi_com_truth");
                hist_polarization_phi_com_truth.Write();

                TH1D hist_decay_costheta_helicity_truth     = *rdf.Histo1D({("decay_costheta_helicity_truth_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"decay_costheta_helicity_truth");
                hist_decay_costheta_helicity_truth.Write();
                TH1D hist_decay_phi_helicity_truth          = *rdf.Histo1D({("decay_phi_helicity_truth_"+ label).c_str(), ";#omega_{helicity} (deg);Counts", 9, -180.0, 180.0},"decay_phi_helicity_truth");
                hist_decay_phi_helicity_truth.Write();
                TH1D hist_psi_helicity_truth                = *rdf.Histo1D({("psi_helicity_truth_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 9, -270.0, 270.0},"psi_helicity_truth");
                hist_psi_helicity_truth.Write();

                TH2D hist_beam_energy_kin_truth             = *rdf.Histo2D({("beam_energy_kin_truth_"+ label).c_str(), ";E_{beam} (GeV);E_{beam} (GeV)", 60, 5.0, 11.0, 60, 5.0, 11.0},"beam_energy_kin","beam_energy_truth");
                hist_beam_energy_kin_truth.Write();
                TH2D hist_minust_kin_truth                  = *rdf.Histo2D({("minust_kin_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});-t (GeV^{2}/c^{2})", 30, 0.0, 3.0, 30, 0.0, 3.0},"minust_kin","minust_truth");
                hist_minust_kin_truth.Write();
            }
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}