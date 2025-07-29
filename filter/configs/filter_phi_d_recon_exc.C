#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

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

void filter_phi_d_recon_exc(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_recon_exc_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_phi_d_recon";
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());

    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    auto rdf_input = rdf_def
    .Define("target_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")
    .Define("sim_weight",                       "sim_weight_func(beam_p4_truth.E(), -(target_p4 - d_p4_truth).Mag2())")
    .Define("event_weight",                     "beam_accid_weight*combo_accid_weight*sim_weight")

    .Define("beam_energy_meas",                 "beam_p4_meas.E()")
    .Define("beam_energy_kin",                  "beam_p4_kin.E()")
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("beam_energy_diff",                 "beam_energy_kin - beam_energy_truth")
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
    .Define("phi_mass_diff",                    "phi_mass_kin - phi_mass_truth")
    .Define("phi_theta_meas",                   "phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin",                    "phi_p4_kin.Theta()*RadToDeg")
    .Define("phi_theta_truth",                  "phi_p4_truth.Theta()*RadToDeg")

    .Define("miss_p4_meas",                     "beam_p4_meas + target_p4 - phi_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                      "beam_p4_kin + target_p4 - phi_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                    "beam_p4_truth + target_p4 - phi_p4_truth - d_p4_truth")
    .Define("miss_energy_meas",                 "miss_p4_meas.E()")
    .Define("miss_energy_kin",                  "miss_p4_kin.E()")
    .Define("miss_energy_truth",                "miss_p4_truth.E()")
    .Define("miss_masssquared_meas",            "miss_p4_meas.M2()")
    .Define("miss_masssquared_kin",             "miss_p4_kin.M2()")
    .Define("miss_masssquared_truth",           "miss_p4_truth.M2()")
    .Define("miss_pminus_meas",                 "miss_p4_meas.Minus()")
    .Define("miss_pminus_kin",                  "miss_p4_kin.Minus()")
    .Define("miss_pminus_truth",                "miss_p4_truth.Minus()")
    .Define("miss_momentum_meas",               "miss_p4_meas.P()")
    .Define("miss_momentum_kin",                "miss_p4_kin.P()")
    .Define("miss_momentum_truth",              "miss_p4_truth.P()")

    .Define("total_p4_meas_initial",            "beam_p4_meas + target_p4")
    .Define("total_p4_meas_final",              "phi_p4_meas + d_p4_meas")
    .Define("total_p4_kin",                     "phi_p4_kin + d_p4_kin")
    .Define("total_p4_truth",                   "phi_p4_truth + d_p4_truth")
    .Define("minust_meas_target",               "-(d_p4_meas - target_p4).Mag2()")
    .Define("minust_meas_beam",                 "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minust_kin",                       "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minust_truth",                     "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minust_diff",                      "minust_kin - minust_truth")
    .Define("vertex_z_meas",                    "beam_x4_meas.Z()")
    .Define("vertex_z_kin",                     "beam_x4_kin.Z()")
    .Define("vertex_z_truth",                   "beam_x4_truth.Z()")
    .Define("vertex_x_meas",                    "beam_x4_meas.X()")
    .Define("vertex_x_kin",                     "beam_x4_kin.X()")
    .Define("vertex_x_truth",                   "beam_x4_truth.X()")
    .Define("vertex_y_meas",                    "beam_x4_meas.Y()")
    .Define("vertex_y_kin",                     "beam_x4_kin.Y()")
    .Define("vertex_y_truth",                   "beam_x4_truth.Y()")
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
    .Define("beam_p4_com_meas",                 "boost_lorentz_vector(beam_p4_meas, -total_p4_meas_final.BoostVector())")
    .Define("beam_p4_com_kin",                  "boost_lorentz_vector(beam_p4_kin, -total_p4_kin.BoostVector())")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_p4_truth.BoostVector())")
    .Define("phi_p4_com_meas",                  "boost_lorentz_vector(phi_p4_meas, -total_p4_meas_final.BoostVector())")
    .Define("phi_p4_com_kin",                   "boost_lorentz_vector(phi_p4_kin, -total_p4_kin.BoostVector())")
    .Define("phi_p4_com_truth",                 "boost_lorentz_vector(phi_p4_truth, -total_p4_truth.BoostVector())")
    .Define("kp_p4_com_meas",                   "boost_lorentz_vector(kp_p4_meas, -total_p4_meas_final.BoostVector())")
    .Define("kp_p4_com_kin",                    "boost_lorentz_vector(kp_p4_kin, -total_p4_kin.BoostVector())")
    .Define("kp_p4_com_truth",                  "boost_lorentz_vector(kp_p4_truth, -total_p4_truth.BoostVector())")
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
    .Define("polarization_phi_com_diff",        "polarization_phi_com_kin - polarization_phi_com_truth")
    .Define("scatter_theta_com_meas",           "phi_p4_com_meas.Vect().Angle(z_x3_com_meas)*RadToDeg")
    .Define("scatter_theta_com_kin",            "phi_p4_com_kin.Vect().Angle(z_x3_com_kin)*RadToDeg")
    .Define("scatter_theta_com_truth",          "phi_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")

    .Define("z_x3_helicity_meas",               "phi_p4_com_meas.Vect().Unit()")
    .Define("z_x3_helicity_kin",                "phi_p4_com_kin.Vect().Unit()")
    .Define("z_x3_helicity_truth",              "phi_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_meas",               "beam_p4_com_meas.Vect().Cross(phi_p4_com_meas.Vect()).Unit()")
    .Define("y_x3_helicity_kin",                "beam_p4_com_kin.Vect().Cross(phi_p4_com_kin.Vect()).Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_meas",               "y_x3_helicity_meas.Cross(z_x3_helicity_meas).Unit()")
    .Define("x_x3_helicity_kin",                "y_x3_helicity_kin.Cross(z_x3_helicity_kin).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_meas",              "boost_lorentz_vector(kp_p4_com_meas, -phi_p4_com_meas.BoostVector()).Vect().Unit()")
    .Define("pi_x3_helicity_kin",               "boost_lorentz_vector(kp_p4_com_kin, -phi_p4_com_kin.BoostVector()).Vect().Unit()")
    .Define("pi_x3_helicity_truth",             "boost_lorentz_vector(kp_p4_com_truth, -phi_p4_com_truth.BoostVector()).Vect().Unit()")
    .Define("decay_costheta_helicity_meas",     "pi_x3_helicity_meas.Dot(z_x3_helicity_meas)")
    .Define("decay_costheta_helicity_kin",      "pi_x3_helicity_kin.Dot(z_x3_helicity_kin)")
    .Define("decay_costheta_helicity_truth",    "pi_x3_helicity_truth.Dot(z_x3_helicity_truth)")
    .Define("decay_costheta_helicity_diff",     "decay_costheta_helicity_kin - decay_costheta_helicity_truth")
    .Define("decay_phi_helicity_meas",          "TMath::ATan2(-x_x3_helicity_meas.Dot(pi_x3_helicity_meas.Cross(z_x3_helicity_meas)), y_x3_helicity_meas.Dot(pi_x3_helicity_meas.Cross(z_x3_helicity_meas)))*RadToDeg")
    .Define("decay_phi_helicity_kin",           "TMath::ATan2(-x_x3_helicity_kin.Dot(pi_x3_helicity_kin.Cross(z_x3_helicity_kin)), y_x3_helicity_kin.Dot(pi_x3_helicity_kin.Cross(z_x3_helicity_kin)))*RadToDeg")
    .Define("decay_phi_helicity_truth",         "TMath::ATan2(-x_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)), y_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)))*RadToDeg")
    .Define("decay_phi_helicity_diff",          "decay_phi_helicity_kin - decay_phi_helicity_truth")
    .Define("psi_helicity_meas",                "fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0) >= 180 ? fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0) - 360 : fmod(polarization_phi_com_meas-decay_phi_helicity_meas+360, 360.0)")
    .Define("psi_helicity_kin",                 "fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0) >= 180 ? fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0) - 360 : fmod(polarization_phi_com_kin-decay_phi_helicity_kin+360, 360.0)")
    .Define("psi_helicity_truth",               "fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0)")
    .Define("psi_helicity_diff",                "psi_helicity_kin - psi_helicity_truth")
    .Define("Psi_helicity_meas",                "fmod(polarization_phi_com_meas+decay_phi_helicity_meas+360, 360.0) >= 180 ? fmod(polarization_phi_com_meas+decay_phi_helicity_meas+360, 360.0) - 360 : fmod(polarization_phi_com_meas+decay_phi_helicity_meas+360, 360.0)")
    .Define("Psi_helicity_kin",                 "fmod(polarization_phi_com_kin+decay_phi_helicity_kin+360, 360.0) >= 180 ? fmod(polarization_phi_com_kin+decay_phi_helicity_kin+360, 360.0) - 360 : fmod(polarization_phi_com_kin+decay_phi_helicity_kin+360, 360.0)")
    .Define("Psi_helicity_truth",               "fmod(polarization_phi_com_truth+decay_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth+decay_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth+decay_phi_helicity_truth+360, 360.0)")
    .Define("Psi_helicity_diff",                "Psi_helicity_kin - Psi_helicity_truth")
    ;

    cout << "Filtering events...\n";
    string KinematicsCut    = "kp_momentum_meas > 0.4 && km_momentum_meas > 0.4 && d_momentum_meas > 0.4 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
    string dEdxCut          = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.3*d_momentum_meas+4.1)+2.3) && d_dedx_st_keV_per_cm_meas > (TMath::Exp(-1.9*d_momentum_meas+2.8)+0.6)";
    string VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
    string KinFitChiSqCut   = "chisq_per_ndf_kin < 5.0";
    string PhiMassCut       = "phi_mass_kin > 1.005 && phi_mass_kin < 1.04";
    auto rdf_NoCut          = rdf_input;
    auto rdf_KinematicsCut  = rdf_input.Filter(dEdxCut.c_str()).Filter(VertexCut.c_str()).Filter(KinFitChiSqCut.c_str());
    auto rdf_dEdxCut        = rdf_input.Filter(KinematicsCut.c_str()).Filter(VertexCut.c_str()).Filter(KinFitChiSqCut.c_str());
    auto rdf_VertexCut      = rdf_input.Filter(KinematicsCut.c_str()).Filter(dEdxCut.c_str()).Filter(KinFitChiSqCut.c_str());
    auto rdf_KinFitChiSqCut = rdf_input.Filter(KinematicsCut.c_str()).Filter(dEdxCut.c_str()).Filter(VertexCut.c_str());
    auto rdf_NominalCut     = rdf_input.Filter(KinematicsCut.c_str()).Filter(dEdxCut.c_str()).Filter(VertexCut.c_str()).Filter(KinFitChiSqCut.c_str());
    auto rdf_AllCut         = rdf_input.Filter(KinematicsCut.c_str()).Filter(dEdxCut.c_str()).Filter(VertexCut.c_str()).Filter(KinFitChiSqCut.c_str()).Filter(PhiMassCut.c_str());
    RNode rdfs []           = {rdf_NoCut,   rdf_KinematicsCut,  rdf_dEdxCut,    rdf_VertexCut,  rdf_KinFitChiSqCut, rdf_NominalCut, rdf_AllCut};
    string labels []        = {"NoCut",     "KinematicsCut",    "dEdxCut",      "VertexCut",    "KinFitChiSqCut",   "NominalCut",   "AllCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_exc_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_phi_d_recon";
        rdf_AllCut.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_exc_%s.root",reaction.c_str());
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
            TH1D hist_combo_accid_weight                    = *rdf.Histo1D({("combo_accid_weight_"+ label).c_str(), ";Combo Accidental Weight;Counts", 2, -0.5, 1.5},"combo_accid_weight");
            hist_combo_accid_weight.Write();

            TH1D hist_beam_energy_meas                      = *rdf.Histo1D({("beam_energy_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_meas","event_weight");
            hist_beam_energy_meas.Write();
            TH1D hist_beam_DeltaT_meas                      = *rdf.Histo1D({("beam_DeltaT_"+ label).c_str(), ";#Delta t_{beam} (ns);Counts", 100, -20.0, 20.0},"beam_DeltaT_meas","beam_accid_weight");
            hist_beam_DeltaT_meas.Write();

            TH1D hist_kp_pidfom                             = *rdf.Histo1D({("kp_pidfom_"+ label).c_str(), ";kp_pidfom;Counts", 100, 0.0, 1.0},"kp_pidfom","event_weight");
            hist_kp_pidfom.Write();
            TH1D hist_kp_DeltaT_meas                        = *rdf.Histo1D({("kp_DeltaT_meas_"+ label).c_str(), ";#Delta t_{K^{+}} (ns);Counts", 100, -2.0, 2.0},"kp_DeltaT_meas","event_weight");
            hist_kp_DeltaT_meas.Write();
            TH2D hist_kp_DeltaT_momentum_meas               = *rdf.Histo2D({("kp_DeltaT_momentum_meas_"+ label).c_str(), ";p (GeV/c);#Delta t_{K^{+}} (ns)", 100, 0.0, 10.0, 100, -2.0, 2.0},"kp_momentum_meas","kp_DeltaT_meas","event_weight");
            hist_kp_DeltaT_momentum_meas.Write();
            TH2D hist_kp_dEdx_cdc_meas                      = *rdf.Histo2D({("kp_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_cdc_meas.Write();
            TH2D hist_kp_dEdx_fdc_meas                      = *rdf.Histo2D({("kp_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_fdc_meas.Write();
            TH2D hist_kp_dEdx_tof_meas                      = *rdf.Histo2D({("kp_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_tof_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_tof_meas.Write();
            TH2D hist_kp_dEdx_st_meas                       = *rdf.Histo2D({("kp_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"kp_momentum_meas","kp_dedx_st_keV_per_cm_meas","event_weight");
            hist_kp_dEdx_st_meas.Write();
            TH2D hist_kp_kinematics_meas                    = *rdf.Histo2D({("kp_kinematics_meas_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","event_weight");
            hist_kp_kinematics_meas.Write();
            TH2D hist_kp_kinematics_fdc_meas                = *rdf.Histo2D({("kp_kinematics_fdc_meas_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","kp_in_fdc_meas");
            hist_kp_kinematics_fdc_meas.Write();
            TH2D hist_kp_kinematics_fdc_cdc_meas            = *rdf.Histo2D({("kp_kinematics_fdc_cdc_meas_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","kp_in_fdc_cdc_meas");
            hist_kp_kinematics_fdc_cdc_meas.Write();
            TH2D hist_kp_kinematics_cdc_meas                = *rdf.Histo2D({("kp_kinematics_cdc_meas_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","kp_in_cdc_meas");
            hist_kp_kinematics_cdc_meas.Write();
            TH2D hist_kp_kinematics_neither_meas            = *rdf.Histo2D({("kp_kinematics_neither_meas_"+ label).c_str(), ";P_{K^{+}} (GeV/c);#theta_{K^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"kp_momentum_meas","kp_theta_meas","kp_in_neither_meas");
            hist_kp_kinematics_neither_meas.Write();

            TH1D hist_km_pidfom                             = *rdf.Histo1D({("km_pidfom_"+ label).c_str(), ";km_pidfom;Counts", 100, 0.0, 1.0},"km_pidfom","event_weight");
            hist_km_pidfom.Write();
            TH1D hist_km_DeltaT_meas                        = *rdf.Histo1D({("km_DeltaT_meas_"+ label).c_str(), ";#Delta t_{K^{-}} (ns);Counts", 100, -2.0, 2.0},"km_DeltaT_meas","event_weight");
            hist_km_DeltaT_meas.Write();
            TH2D hist_km_DeltaT_momentum_meas               = *rdf.Histo2D({("km_DeltaT_momentum_meas_"+ label).c_str(), ";p (GeV/c);#Delta t_{K^{-}} (ns)", 100, 0.0, 10.0, 100, -2.0, 2.0},"km_momentum_meas","km_DeltaT_meas","event_weight");
            hist_km_DeltaT_momentum_meas.Write();
            TH2D hist_km_dEdx_cdc_meas                      = *rdf.Histo2D({("km_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_km_dEdx_cdc_meas.Write();
            TH2D hist_km_dEdx_fdc_meas                      = *rdf.Histo2D({("km_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_km_dEdx_fdc_meas.Write();
            TH2D hist_km_dEdx_tof_meas                      = *rdf.Histo2D({("km_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_tof_keV_per_cm_meas","event_weight");
            hist_km_dEdx_tof_meas.Write();
            TH2D hist_km_dEdx_st_meas                       = *rdf.Histo2D({("km_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 100, 0.0, 40},"km_momentum_meas","km_dedx_st_keV_per_cm_meas","event_weight");
            hist_km_dEdx_st_meas.Write();
            TH2D hist_km_kinematics_meas                    = *rdf.Histo2D({("km_kinematics_meas_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","event_weight");
            hist_km_kinematics_meas.Write();
            TH2D hist_km_kinematics_fdc_meas                = *rdf.Histo2D({("km_kinematics_fdc_meas_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","km_in_fdc_meas");
            hist_km_kinematics_fdc_meas.Write();
            TH2D hist_km_kinematics_fdc_cdc_meas            = *rdf.Histo2D({("km_kinematics_fdc_cdc_meas_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","km_in_fdc_cdc_meas");
            hist_km_kinematics_fdc_cdc_meas.Write();
            TH2D hist_km_kinematics_cdc_meas                = *rdf.Histo2D({("km_kinematics_cdc_meas_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","km_in_cdc_meas");
            hist_km_kinematics_cdc_meas.Write();
            TH2D hist_km_kinematics_neither_meas            = *rdf.Histo2D({("km_kinematics_neither_meas_"+ label).c_str(), ";P_{K^{-}} (GeV/c);#theta_{K^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"km_momentum_meas","km_theta_meas","km_in_neither_meas");
            hist_km_kinematics_neither_meas.Write();

            TH1D hist_d_DeltaT_meas                         = *rdf.Histo1D({("d_DeltaT_meas_"+ label).c_str(), ";#Delta t_{d} (ns);Counts", 100, -5.0, 5.0},"d_DeltaT_meas","event_weight");
            hist_d_DeltaT_meas.Write();
            TH2D hist_d_DeltaT_momentum_meas                = *rdf.Histo2D({("d_DeltaT_momentum_meas_"+ label).c_str(), ";p (GeV/c);#Delta t_{d} (ns)", 200, 0.0, 2.0, 100, -5.0, 5.0},"d_momentum_meas","d_DeltaT_meas","event_weight");
            hist_d_DeltaT_momentum_meas.Write();
            TH2D hist_d_dEdx_cdc_meas                       = *rdf.Histo2D({("d_dEdx_cdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 200, 0.0, 2.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_cdc_keV_per_cm_meas","event_weight");
            hist_d_dEdx_cdc_meas.Write();
            TH2D hist_d_dEdx_fdc_meas                       = *rdf.Histo2D({("d_dEdx_fdc_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 200, 0.0, 2.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_fdc_keV_per_cm_meas","event_weight");
            hist_d_dEdx_fdc_meas.Write();
            TH2D hist_d_dEdx_tof_meas                       = *rdf.Histo2D({("d_dEdx_tof_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 200, 0.0, 2.0, 100, 0.0, 40},"d_momentum_meas","d_dedx_tof_keV_per_cm_meas","event_weight");
            hist_d_dEdx_tof_meas.Write();
            TH2D hist_d_dEdx_st_meas                        = *rdf.Histo2D({("d_dEdx_st_meas_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 200, 0.0, 2.0, 100, 0.0, 20},"d_momentum_meas","d_dedx_st_keV_per_cm_meas","event_weight");
            hist_d_dEdx_st_meas.Write();
            TH2D hist_d_kinematics_meas                     = *rdf.Histo2D({("d_kinematics_meas_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 200, 0.0, 2.0, 180, 0.0, 180.0},"d_momentum_meas","d_theta_meas","event_weight");
            hist_d_kinematics_meas.Write();

            TH1D hist_miss_energy_meas                      = *rdf.Histo1D({("miss_energy_meas_"+ label).c_str(), ";E_{miss} (GeV);Counts", 300, -3.0, 3.0},"miss_energy_meas","event_weight");
            hist_miss_energy_meas.Write();
            TH1D hist_miss_masssquared_meas                 = *rdf.Histo1D({("miss_masssquared_meas_"+ label).c_str(), ";m_{miss}^{2} (GeV^{2}/c^{4});Counts", 100, -0.1, 0.1},"miss_masssquared_meas","event_weight");
            hist_miss_masssquared_meas.Write();
            TH1D hist_miss_pminus_meas                      = *rdf.Histo1D({("miss_pminus_meas_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 400, -0.2, 0.2},"miss_pminus_meas","event_weight");
            hist_miss_pminus_meas.Write();
            TH1D hist_miss_momentum_meas                    = *rdf.Histo1D({("miss_momentum_meas_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 200, 0.0, 2.0},"miss_momentum_meas","event_weight");
            hist_miss_momentum_meas.Write();
            TH2D hist_miss_momentum_energy_meas             = *rdf.Histo2D({("miss_momentum_energy_meas_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} (GeV)", 200, 0.0, 2.0, 300, -3.0, 3.0},"miss_momentum_meas","miss_energy_meas","event_weight");
            hist_miss_momentum_energy_meas.Write();
            TH2D hist_miss_energy_phi_mass_meas             = *rdf.Histo2D({("miss_energy_phi_mass_meas_"+ label).c_str(), ";E_{miss} (GeV);m_{K^{+}K^{-}} (GeV/c)", 300, -3.0, 3.0, 500, 0.9, 1.9},"miss_energy_meas","phi_mass_meas","event_weight");
            hist_miss_energy_phi_mass_meas.Write();
            TH1D hist_coplanarity_meas                      = *rdf.Histo1D({("coplanarity_meas_"+ label).c_str(), ";Coplanarity (deg);Counts", 40, 160.0, 200.0},"coplanarity_meas","event_weight");
            hist_coplanarity_meas.Write();

            TH1D hist_phi_mass_kin                          = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 300, 0.9, 1.5},"phi_mass_kin","event_weight");
            hist_phi_mass_kin.Write();
            TH2D hist_phi_mass_chisq_kin                    = *rdf.Histo2D({("phi_mass_chisq_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);KinFit #Chi^2/NDF", 300, 0.9, 1.5, 100, 0.0, 10.0},"phi_mass_kin","chisq_per_ndf_kin","event_weight");
            hist_phi_mass_chisq_kin.Write();
            TH2D hist_phi_kinematics_kin                    = *rdf.Histo2D({("phi_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 110, 0.0, 11.0, 180, 0.0, 180.0},"phi_momentum_kin","phi_theta_kin","event_weight");
            hist_phi_kinematics_kin.Write();

            TH1D hist_minust_kin                            = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 100, 0.0, 2.0},"minust_kin","event_weight");
            hist_minust_kin.Write();
            TH1D hist_vertex_z_kin                          = *rdf.Histo1D({("vertex_z_kin_"+ label).c_str(), ";Z_{vertex} (cm);Counts", 100, 40.0, 90.0},"vertex_z_kin","event_weight");
            hist_vertex_z_kin.Write();
            TH2D hist_vertex_x_y_kin                        = *rdf.Histo2D({("vertex_x_y_kin_"+ label).c_str(), ";X_{vertex} (cm);Y_{vertex} (cm)", 100, -2.0, 2.0, 100, -2.0, 2.0},"vertex_x_kin","vertex_y_kin","event_weight");
            hist_vertex_x_y_kin.Write();
            TH1D hist_rho_mass_kin                          = *rdf.Histo1D({("rho_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 100, 0.0, 1.0},"rho_mass_kin","event_weight");
            hist_rho_mass_kin.Write();
            TH1D hist_kinfit_fom_kin                        = *rdf.Histo1D({("kinfit_fom_kin_"+ label).c_str(), ";log(KinFit FOM);Counts", 30, -15.0, 0},"log10_kinfit_fom_kin","event_weight");
            hist_kinfit_fom_kin.Write();
            TH1D hist_chisq_per_ndf_kin                     = *rdf.Histo1D({("chisq_per_ndf_kin_"+ label).c_str(), ";#chi^{2}/NDF;Counts", 20, 0.0, 10.0},"chisq_per_ndf_kin","event_weight");
            hist_chisq_per_ndf_kin.Write();
            TH2D hist_kinfit_fom_chisq_per_ndf_kin          = *rdf.Histo2D({("kinfit_fom_chisq_per_ndf_kin_"+ label).c_str(), ";log(KinFit FOM);#chi^{2}/NDF", 150, -15.0, 0.0, 100, 0.0, 10.0},"log10_kinfit_fom_kin","chisq_per_ndf_kin","event_weight");
            hist_kinfit_fom_chisq_per_ndf_kin.Write();
            TH2D hist_beam_energy_minust_kin                = *rdf.Histo2D({("beam_energy_minust_kin_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 20, 0.0, 2.0},"beam_energy_kin","minust_kin","event_weight");
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
            TH1D hist_Psi_helicity_kin                      = *rdf.Histo1D({("Psi_helicity_kin_"+ label).c_str(), ";#Psi_{helicity} (deg);Counts", 18, -180.0, 180.0},"Psi_helicity_kin","event_weight");
            hist_Psi_helicity_kin.Write();

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
                TH1D hist_phi_mass_truth                    = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 500, 0.9, 1.9},"phi_mass_truth","event_weight");
                hist_phi_mass_truth.Write();
                TH2D hist_phi_kinematics_truth              = *rdf.Histo2D({("phi_kinematics_truth_"+ label).c_str(), ";P_{#phi} (GeV/c);#theta_{#phi} (deg)", 100, 0.0, 11.0, 180, 0.0, 180.0},"phi_momentum_truth","phi_theta_truth","event_weight");
                hist_phi_kinematics_truth.Write();
                TH2D hist_phi_d_theta_truth                 = *rdf.Histo2D({("phi_d_theta_truth_"+ label).c_str(), ";#theta_{d} (deg);#theta_{#phi} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"d_theta_truth","phi_theta_truth","event_weight");
                hist_phi_d_theta_truth.Write();
                TH2D hist_phi_d_momentum_truth              = *rdf.Histo2D({("phi_d_momentum_truth_"+ label).c_str(), ";P_{d} (GeV/c);P_{#phi} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"d_momentum_truth","phi_momentum_truth","event_weight");
                hist_phi_d_momentum_truth.Write();
                TH1D hist_minust_truth                      = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 100, 0.0, 2.0},"minust_truth","event_weight");
                hist_minust_truth.Write();
                TH1D hist_rho_mass_truth                    = *rdf.Histo1D({("rho_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 100, 0.0, 1.0},"rho_mass_truth","event_weight");
                hist_rho_mass_truth.Write();
                TH2D hist_beam_energy_minust_truth          = *rdf.Histo2D({("beam_energy_minust_truth_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 20, 0.0, 2.0},"beam_energy_truth","minust_truth","event_weight");
                hist_beam_energy_minust_truth.Write();
                TH1D hist_scatter_theta_com_truth           = *rdf.Histo1D({("scatter_theta_com_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"scatter_theta_com_truth","event_weight");
                hist_scatter_theta_com_truth.Write();
                TH2D hist_minust_scatter_theta_com_truth    = *rdf.Histo2D({("minust_scatter_theta_com_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 20, 0.0, 2.0, 180, 0.0, 180.0},"minust_truth","scatter_theta_com_truth","event_weight");
                hist_minust_scatter_theta_com_truth.Write();
                TH1D hist_polarization_phi_com_truth        = *rdf.Histo1D({("polarization_phi_com_truth_"+ label).c_str(), ";#phi_{com} (deg);Counts", 9, -180, 180.0},"polarization_phi_com_truth","event_weight");
                hist_polarization_phi_com_truth.Write();
                TH1D hist_decay_costheta_helicity_truth     = *rdf.Histo1D({("decay_costheta_helicity_truth_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"decay_costheta_helicity_truth","event_weight");
                hist_decay_costheta_helicity_truth.Write();
                TH1D hist_decay_phi_helicity_truth          = *rdf.Histo1D({("decay_phi_helicity_truth_"+ label).c_str(), ";#phi_{helicity} (deg);Counts", 9, -180.0, 180.0},"decay_phi_helicity_truth","event_weight");
                hist_decay_phi_helicity_truth.Write();
                TH1D hist_psi_helicity_truth                = *rdf.Histo1D({("psi_helicity_truth_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 9, -270.0, 270.0},"psi_helicity_truth","event_weight");
                hist_psi_helicity_truth.Write();
                TH1D hist_Psi_helicity_truth                = *rdf.Histo1D({("Psi_helicity_truth_"+ label).c_str(), ";#Psi_{helicity} (deg);Counts", 9, -270.0, 270.0},"Psi_helicity_truth","event_weight");
                hist_Psi_helicity_truth.Write();

                TH2D hist_beam_energy_kin_truth             = *rdf.Histo2D({("beam_energy_kin_truth_"+ label).c_str(), ";E_{beam} (GeV);E_{beam} (GeV)", 60, 5.0, 11.0, 60, 5.0, 11.0},"beam_energy_kin","beam_energy_truth","event_weight");
                hist_beam_energy_kin_truth.Write();
                TH2D hist_minust_kin_truth                  = *rdf.Histo2D({("minust_kin_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});-t (GeV^{2}/c^{2})", 20, 0.0, 2.0, 20, 0.0, 2.0},"minust_kin","minust_truth","event_weight");
                hist_minust_kin_truth.Write();
                TH2D hist_beam_energy_diff                  = *rdf.Histo2D({("beam_energy_diff_"+ label).c_str(), ";E_{beam}^{truth} (GeV);E_{beam}^{kin} - E_{beam}^{truth} (GeV)", 60, 5.0, 11.0, 100, -0.5, 0.5},"beam_energy_truth","beam_energy_diff","event_weight");
                hist_beam_energy_diff.Write();
                TH2D hist_minust_diff                       = *rdf.Histo2D({("minust_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});-t_{kin} - -t_{truth} (GeV^{2}/c^{2})", 40, 0.0, 2.0, 80, -0.2, 0.2},"minust_truth","minust_diff","event_weight");
                hist_minust_diff.Write();
                TH2D hist_polarization_phi_com_diff         = *rdf.Histo2D({("polarization_phi_com_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});#phi_{com}^{kin} - #phi_{com}^{truth} (deg)", 40, 0.0, 2.0, 80, -2.0, 2.0},"minust_truth","polarization_phi_com_diff","event_weight");
                hist_polarization_phi_com_diff.Write();
                TH2D hist_decay_costheta_helicity_diff      = *rdf.Histo2D({("decay_costheta_helicity_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});cos(#theta_{helicity}^{kin}) - cos(#theta_{helicity}^{truth})", 40, 0.0, 2.0, 80, -1.0, 1.0},"minust_truth","decay_costheta_helicity_diff","event_weight");
                hist_decay_costheta_helicity_diff.Write();
                TH2D hist_decay_phi_helicity_diff           = *rdf.Histo2D({("decay_phi_helicity_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});#phi_{helicity}^{kin} - #phi_{helicity}^{truth} (deg)", 40, 0.0, 2.0, 80, -8.0, 8.0},"minust_truth","decay_phi_helicity_diff","event_weight");
                hist_decay_phi_helicity_diff.Write();
                TH2D hist_psi_helicity_diff                 = *rdf.Histo2D({("psi_helicity_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});#psi_{helicity}^{kin} - #psi_{helicity}^{truth} (deg)", 40, 0.0, 2.0, 80, -8.0, 8.0},"minust_truth","psi_helicity_diff","event_weight");
                hist_psi_helicity_diff.Write();
                TH2D hist_Psi_helicity_diff                 = *rdf.Histo2D({("Psi_helicity_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});#Psi_{helicity}^{kin} - #Psi_{helicity}^{truth} (deg)", 40, 0.0, 2.0, 80, -8.0, 8.0},"minust_truth","Psi_helicity_diff","event_weight");
                hist_Psi_helicity_diff.Write();
                TH2D hist_phi_mass_diff                     = *rdf.Histo2D({("phi_mass_diff_"+ label).c_str(), ";-t_{truth} (GeV^{2}/c^{2});m_{K^{+}K^{-}}^{kin} - m_{K^{+}K^{-}}^{truth} (GeV/c^{2})", 40, 0.0, 2.0, 80, -0.02, 0.02},"minust_truth","phi_mass_diff","event_weight");
                hist_phi_mass_diff.Write();
            }
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}