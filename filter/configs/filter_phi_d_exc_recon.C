#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double dxs_weight_func(double beam_energy_truth, double minust_truth, int sim_model_flag)
{
    double a1 = 10222.59309;
    double b1 = 19.33908;
    double a2 = 20.23543;
    double b2 = 3.33070;
    double normalization = 1.0;
    if (beam_energy_truth < 0.01 || sim_model_flag > 0)   // data with its truth variable set to zero as placeholder, or simulation with a cross section model already
        return 1.0;
    else                                                  // simulation with flat cross section, weighted by the double exponential function
        return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}

double psi_weight_func(double polarization_angle, double beam_energy_truth, double psi_helicity_truth)
{
    double polarization_degree[] = {                                0.0000,
                                    0.0614, 0.1087, 0.0971, 0.0763, 0.0629,
                                    0.0517, 0.0413, 0.1486, 0.1490, 0.2020,
                                    0.2395, 0.2812, 0.3042, 0.2247, 0.0528,
                                    0.0443, 0.0521, 0.0744, 0.0376, 0.0290,
                                    0.0132, 0.0359, 0.0282, 0.0576};  // 5.8-10.8 GeV, in steps of 0.2 GeV
    int energy_bin = (beam_energy_truth - 5.8) / 0.2;
    if (beam_energy_truth < 0.01)        // data with its truth variable set to zero as placeholder
        return 1.0;
    else if (polarization_angle < -0.01) // amorphous simulation runs
        return 1.0;
    else                                 // diamond simulation runs, to be weighted by the SCHC+NPE predictions
        return 1.0 + polarization_degree[energy_bin]*TMath::Cos(2*psi_helicity_truth/RadToDeg);
}

void filter_phi_d_exc_recon(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_exc_recon_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_phi_d_exc_recon";
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());

    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    if (reaction.find("model") != string::npos)
        rdf_def = rdf_def.Define("sim_model_flag", "1");
    else
        rdf_def = rdf_def.Define("sim_model_flag", "-1");
    auto rdf_input = rdf_def
    .Define("target_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_energy_meas",                 "beam_p4_meas .E()")
    .Define("beam_energy_kin",                  "beam_p4_kin  .E()")
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("beam_energy_diff",                 "beam_energy_kin - beam_energy_truth")
    .Define("beam_DeltaT_meas",                 "rftime + (beam_x4_meas.Z()-65.0)/29.9792458 - beam_x4_meas.T()")
    .Define("beam_DeltaT_kin",                  "rftime + (beam_x4_kin .Z()-65.0)/29.9792458 - beam_x4_kin .T()")

    .Define("kp_as_pion_p4_meas",               "TLorentzVector(kp_p4_meas .Vect(), TMath::Sqrt(kp_p4_meas .P()*kp_p4_meas .P() + mass_piplus*mass_piplus))")
    .Define("kp_as_pion_p4_kin",                "TLorentzVector(kp_p4_kin  .Vect(), TMath::Sqrt(kp_p4_kin  .P()*kp_p4_kin  .P() + mass_piplus*mass_piplus))")
    .Define("kp_as_pion_p4_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_energy_meas",                   "kp_p4_meas .E()")
    .Define("kp_energy_kin",                    "kp_p4_kin  .E()")
    .Define("kp_energy_truth",                  "kp_p4_truth.E()")
    .Define("kp_momentum_meas",                 "kp_p4_meas .P()")
    .Define("kp_momentum_kin",                  "kp_p4_kin  .P()")
    .Define("kp_momentum_truth",                "kp_p4_truth.P()")
    .Define("kp_theta_meas",                    "kp_p4_meas .Theta()*RadToDeg")
    .Define("kp_theta_kin",                     "kp_p4_kin  .Theta()*RadToDeg")
    .Define("kp_theta_truth",                   "kp_p4_truth.Theta()*RadToDeg")
    .Define("kp_DeltaT_meas",                   "rftime + (kp_x4_meas.Z()-65.0)/29.9792458 - kp_x4_meas.T()")
    .Define("kp_DeltaT_kin",                    "rftime + (kp_x4_kin .Z()-65.0)/29.9792458 - kp_x4_kin .T()")
    .Define("kp_dedx_fdc_keV_per_cm_meas",      "kp_dedx_fdc*1e6")
    .Define("kp_dedx_cdc_keV_per_cm_meas",      "kp_dedx_cdc*1e6")
    .Define("kp_dedx_st_keV_per_cm_meas",       "kp_dedx_st *1e3")
    .Define("kp_dedx_tof_keV_per_cm_meas",      "kp_dedx_tof*1e3")

    .Define("km_as_pion_p4_meas",               "TLorentzVector(km_p4_meas .Vect(), TMath::Sqrt(km_p4_meas .P()*km_p4_meas .P() + mass_piminus*mass_piminus))")
    .Define("km_as_pion_p4_kin",                "TLorentzVector(km_p4_kin  .Vect(), TMath::Sqrt(km_p4_kin  .P()*km_p4_kin  .P() + mass_piminus*mass_piminus))")
    .Define("km_as_pion_p4_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_energy_meas",                   "km_p4_meas .E()")
    .Define("km_energy_kin",                    "km_p4_kin  .E()")
    .Define("km_energy_truth",                  "km_p4_truth.E()")
    .Define("km_momentum_meas",                 "km_p4_meas .P()")
    .Define("km_momentum_kin",                  "km_p4_kin  .P()")
    .Define("km_momentum_truth",                "km_p4_truth.P()")
    .Define("km_theta_meas",                    "km_p4_meas .Theta()*RadToDeg")
    .Define("km_theta_kin",                     "km_p4_kin  .Theta()*RadToDeg")
    .Define("km_theta_truth",                   "km_p4_truth.Theta()*RadToDeg")
    .Define("km_DeltaT_meas",                   "rftime + (km_x4_meas.Z()-65.0)/29.9792458 - km_x4_meas.T()")
    .Define("km_DeltaT_kin",                    "rftime + (km_x4_kin .Z()-65.0)/29.9792458 - km_x4_kin .T()")
    .Define("km_dedx_fdc_keV_per_cm_meas",      "km_dedx_fdc*1e6")
    .Define("km_dedx_cdc_keV_per_cm_meas",      "km_dedx_cdc*1e6")
    .Define("km_dedx_st_keV_per_cm_meas",       "km_dedx_st *1e3")
    .Define("km_dedx_tof_keV_per_cm_meas",      "km_dedx_tof*1e3")

    .Define("d_energy_meas",                    "d_p4_meas .E()")
    .Define("d_energy_kin",                     "d_p4_kin  .E()")
    .Define("d_energy_truth",                   "d_p4_truth.E()")
    .Define("d_momentum_meas",                  "d_p4_meas .P()")
    .Define("d_momentum_kin",                   "d_p4_kin  .P()")
    .Define("d_momentum_truth",                 "d_p4_truth.P()")
    .Define("d_theta_meas",                     "d_p4_meas .Theta()*RadToDeg")
    .Define("d_theta_kin",                      "d_p4_kin  .Theta()*RadToDeg")
    .Define("d_theta_truth",                    "d_p4_truth.Theta()*RadToDeg")
    .Define("d_DeltaT_meas",                    "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_meas.T()")
    .Define("d_DeltaT_kin",                     "rftime + (d_x4_kin .Z()-65.0)/29.9792458 - d_x4_kin .T()")
    .Define("d_dedx_fdc_keV_per_cm_meas",       "d_dedx_fdc*1e6")
    .Define("d_dedx_cdc_keV_per_cm_meas",       "d_dedx_cdc*1e6")
    .Define("d_dedx_st_keV_per_cm_meas",        "d_dedx_st *1e3")
    .Define("d_dedx_tof_keV_per_cm_meas",       "d_dedx_tof*1e3")
    .Define("d_dedx_cdc_residue_meas",          "d_dedx_cdc*1e6-(TMath::Exp(-4.25*d_momentum_meas+5.20) + 3.89)")

    .Define("phi_p4_meas",                      "kp_p4_meas  + km_p4_meas")
    .Define("phi_p4_kin",                       "kp_p4_kin   + km_p4_kin")
    .Define("phi_p4_truth",                     "kp_p4_truth + km_p4_truth")
    .Define("phi_energy_meas",                  "phi_p4_meas .E()")
    .Define("phi_energy_kin",                   "phi_p4_kin  .E()")
    .Define("phi_energy_truth",                 "phi_p4_truth.E()")
    .Define("phi_momentum_meas",                "phi_p4_meas .P()")
    .Define("phi_momentum_kin",                 "phi_p4_kin  .P()")
    .Define("phi_momentum_truth",               "phi_p4_truth.P()")
    .Define("phi_mass_meas",                    "phi_p4_meas .M()")
    .Define("phi_mass_kin",                     "phi_p4_kin  .M()")
    .Define("phi_mass_truth",                   "phi_p4_truth.M()")
    .Define("phi_mass_diff",                    "phi_mass_kin - phi_mass_truth")
    .Define("phi_theta_meas",                   "phi_p4_meas .Theta()*RadToDeg")
    .Define("phi_theta_kin",                    "phi_p4_kin  .Theta()*RadToDeg")
    .Define("phi_theta_truth",                  "phi_p4_truth.Theta()*RadToDeg")

    .Define("miss_p4_meas",                     "beam_p4_meas  + target_p4 - phi_p4_meas  - d_p4_meas")
    .Define("miss_p4_kin",                      "beam_p4_kin   + target_p4 - phi_p4_kin   - d_p4_kin")
    .Define("miss_p4_truth",                    "beam_p4_truth + target_p4 - phi_p4_truth - d_p4_truth")
    .Define("miss_energy_meas",                 "miss_p4_meas. E()")
    .Define("miss_energy_kin",                  "miss_p4_kin  .E()")
    .Define("miss_energy_truth",                "miss_p4_truth.E()")
    .Define("miss_masssquared_meas",            "miss_p4_meas .M2()")
    .Define("miss_masssquared_kin",             "miss_p4_kin  .M2()")
    .Define("miss_masssquared_truth",           "miss_p4_truth.M2()")
    .Define("miss_pminus_meas",                 "miss_p4_meas .Minus()")
    .Define("miss_pminus_kin",                  "miss_p4_kin  .Minus()")
    .Define("miss_pminus_truth",                "miss_p4_truth.Minus()")
    .Define("miss_momentum_meas",               "miss_p4_meas .P()")
    .Define("miss_momentum_kin",                "miss_p4_kin  .P()")
    .Define("miss_momentum_truth",              "miss_p4_truth.P()")

    .Define("total_p4_meas_initial",            "beam_p4_meas + target_p4")
    .Define("total_p4_meas_final",              "phi_p4_meas  + d_p4_meas")
    .Define("total_p4_kin",                     "phi_p4_kin   + d_p4_kin")
    .Define("total_p4_truth",                   "phi_p4_truth + d_p4_truth")
    .Define("minust_meas_target",               "-(d_p4_meas     - target_p4   ).Mag2()")
    .Define("minust_meas_beam",                 "-(beam_p4_meas  - phi_p4_meas ).Mag2()")
    .Define("minust_kin",                       "-(beam_p4_kin   - phi_p4_kin  ).Mag2()")
    .Define("minust_truth",                     "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minust_diff",                      "minust_kin - minust_truth")
    .Define("minusu_meas_target",               "-(phi_p4_meas   - target_p4 ).Mag2()")
    .Define("minusu_meas_beam",                 "-(beam_p4_meas  - d_p4_meas ).Mag2()")
    .Define("minusu_kin",                       "-(beam_p4_kin   - d_p4_kin  ).Mag2()")
    .Define("minusu_truth",                     "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("minusu_diff",                      "minusu_kin - minusu_truth")
    .Define("vertex_z_meas",                    "beam_x4_meas .Z()")
    .Define("vertex_z_kin",                     "beam_x4_kin  .Z()")
    .Define("vertex_z_truth",                   "beam_x4_truth.Z()")
    .Define("vertex_x_meas",                    "beam_x4_meas .X()")
    .Define("vertex_x_kin",                     "beam_x4_kin  .X()")
    .Define("vertex_x_truth",                   "beam_x4_truth.X()")
    .Define("vertex_y_meas",                    "beam_x4_meas .Y()")
    .Define("vertex_y_kin",                     "beam_x4_kin  .Y()")
    .Define("vertex_y_truth",                   "beam_x4_truth.Y()")
    .Define("coplanarity_meas",                 "abs(phi_p4_meas.Phi()  - d_p4_meas.Phi()) *RadToDeg")
    .Define("coplanarity_kin",                  "abs(phi_p4_kin.Phi()   - d_p4_kin.Phi())  *RadToDeg")
    .Define("coplanarity_truth",                "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("rho_mass_meas",                    "(kp_as_pion_p4_meas  + km_as_pion_p4_meas) .M()")
    .Define("rho_mass_kin",                     "(kp_as_pion_p4_kin   + km_as_pion_p4_kin)  .M()")
    .Define("rho_mass_truth",                   "(kp_as_pion_p4_truth + km_as_pion_p4_truth).M()")
    .Define("rho_miss_pminus_meas",             "(beam_p4_meas  + target_p4 - kp_as_pion_p4_meas  - km_as_pion_p4_meas  - d_p4_meas) .Minus()")
    .Define("rho_miss_pminus_kin",              "(beam_p4_kin   + target_p4 - kp_as_pion_p4_kin   - km_as_pion_p4_kin   - d_p4_kin)  .Minus()")
    .Define("rho_miss_pminus_truth",            "(beam_p4_truth + target_p4 - kp_as_pion_p4_truth - km_as_pion_p4_truth - d_p4_truth).Minus()")
    .Define("chisq_per_ndf_kin",                "kin_chisq/kin_ndf")
    .Define("kinfit_fom_kin",                   "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("log10_kinfit_fom_kin",             "TMath::Log10(kinfit_fom_kin)")

    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_meas",                 "boost_lorentz_vector(beam_p4_meas,  -total_p4_meas_final.BoostVector())")
    .Define("beam_p4_com_kin",                  "boost_lorentz_vector(beam_p4_kin,   -total_p4_kin       .BoostVector())")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_p4_truth     .BoostVector())")
    .Define("phi_p4_com_meas",                  "boost_lorentz_vector(phi_p4_meas,   -total_p4_meas_final.BoostVector())")
    .Define("phi_p4_com_kin",                   "boost_lorentz_vector(phi_p4_kin,    -total_p4_kin       .BoostVector())")
    .Define("phi_p4_com_truth",                 "boost_lorentz_vector(phi_p4_truth,  -total_p4_truth     .BoostVector())")
    .Define("kp_p4_com_meas",                   "boost_lorentz_vector(kp_p4_meas,    -total_p4_meas_final.BoostVector())")
    .Define("kp_p4_com_kin",                    "boost_lorentz_vector(kp_p4_kin,     -total_p4_kin       .BoostVector())")
    .Define("kp_p4_com_truth",                  "boost_lorentz_vector(kp_p4_truth,   -total_p4_truth     .BoostVector())")
    .Define("z_x3_com_meas",                    "beam_p4_com_meas .Vect().Unit()")
    .Define("z_x3_com_kin",                     "beam_p4_com_kin  .Vect().Unit()")
    .Define("z_x3_com_truth",                   "beam_p4_com_truth.Vect().Unit()")
    .Define("y_x3_com_meas",                    "beam_p4_com_meas .Vect().Cross(phi_p4_com_meas .Vect()).Unit()")
    .Define("y_x3_com_kin",                     "beam_p4_com_kin  .Vect().Cross(phi_p4_com_kin  .Vect()).Unit()")
    .Define("y_x3_com_truth",                   "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_com_meas",                    "y_x3_com_meas .Cross(z_x3_com_meas ).Unit()")
    .Define("x_x3_com_kin",                     "y_x3_com_kin  .Cross(z_x3_com_kin  ).Unit()")
    .Define("x_x3_com_truth",                   "y_x3_com_truth.Cross(z_x3_com_truth).Unit()")
    .Define("polarization_phi_com_meas",        "TMath::ATan2(-x_x3_com_meas .Dot(epsilon_x3_com.Cross(z_x3_com_meas )), y_x3_com_meas .Dot(epsilon_x3_com.Cross(z_x3_com_meas )))*RadToDeg")
    .Define("polarization_phi_com_kin",         "TMath::ATan2(-x_x3_com_kin  .Dot(epsilon_x3_com.Cross(z_x3_com_kin  )), y_x3_com_kin  .Dot(epsilon_x3_com.Cross(z_x3_com_kin  )))*RadToDeg")
    .Define("polarization_phi_com_truth",       "TMath::ATan2(-x_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)), y_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)))*RadToDeg")
    .Define("polarization_phi_com_diff",        "polarization_phi_com_kin - polarization_phi_com_truth")
    .Define("scatter_theta_com_meas",           "phi_p4_com_meas .Vect().Angle(z_x3_com_meas )*RadToDeg")
    .Define("scatter_theta_com_kin",            "phi_p4_com_kin  .Vect().Angle(z_x3_com_kin  )*RadToDeg")
    .Define("scatter_theta_com_truth",          "phi_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")

    .Define("z_x3_helicity_meas",               "phi_p4_com_meas .Vect().Unit()")
    .Define("z_x3_helicity_kin",                "phi_p4_com_kin  .Vect().Unit()")
    .Define("z_x3_helicity_truth",              "phi_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_meas",               "beam_p4_com_meas .Vect().Cross(phi_p4_com_meas .Vect()).Unit()")
    .Define("y_x3_helicity_kin",                "beam_p4_com_kin  .Vect().Cross(phi_p4_com_kin  .Vect()).Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_meas",               "y_x3_helicity_meas .Cross(z_x3_helicity_meas ).Unit()")
    .Define("x_x3_helicity_kin",                "y_x3_helicity_kin  .Cross(z_x3_helicity_kin  ).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_meas",              "boost_lorentz_vector(kp_p4_com_meas,  -phi_p4_com_meas .BoostVector()).Vect().Unit()")
    .Define("pi_x3_helicity_kin",               "boost_lorentz_vector(kp_p4_com_kin,   -phi_p4_com_kin  .BoostVector()).Vect().Unit()")
    .Define("pi_x3_helicity_truth",             "boost_lorentz_vector(kp_p4_com_truth, -phi_p4_com_truth.BoostVector()).Vect().Unit()")
    .Define("decay_costheta_helicity_meas",     "pi_x3_helicity_meas .Dot(z_x3_helicity_meas )")
    .Define("decay_costheta_helicity_kin",      "pi_x3_helicity_kin  .Dot(z_x3_helicity_kin  )")
    .Define("decay_costheta_helicity_truth",    "pi_x3_helicity_truth.Dot(z_x3_helicity_truth)")
    .Define("decay_costheta_helicity_diff",     "decay_costheta_helicity_kin - decay_costheta_helicity_truth")
    .Define("decay_phi_helicity_meas",          "TMath::ATan2(-x_x3_helicity_meas .Dot(pi_x3_helicity_meas .Cross(z_x3_helicity_meas )), y_x3_helicity_meas .Dot(pi_x3_helicity_meas .Cross(z_x3_helicity_meas )))*RadToDeg")
    .Define("decay_phi_helicity_kin",           "TMath::ATan2(-x_x3_helicity_kin  .Dot(pi_x3_helicity_kin  .Cross(z_x3_helicity_kin  )), y_x3_helicity_kin  .Dot(pi_x3_helicity_kin  .Cross(z_x3_helicity_kin  )))*RadToDeg")
    .Define("decay_phi_helicity_truth",         "TMath::ATan2(-x_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)), y_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)))*RadToDeg")
    .Define("decay_phi_helicity_diff",          "decay_phi_helicity_kin - decay_phi_helicity_truth")
    .Define("psi_helicity_meas",                "fmod(polarization_phi_com_meas -decay_phi_helicity_meas +360, 360.0) >= 180 ? fmod(polarization_phi_com_meas -decay_phi_helicity_meas +360, 360.0) - 360 : fmod(polarization_phi_com_meas -decay_phi_helicity_meas +360, 360.0)")
    .Define("psi_helicity_kin",                 "fmod(polarization_phi_com_kin  -decay_phi_helicity_kin  +360, 360.0) >= 180 ? fmod(polarization_phi_com_kin  -decay_phi_helicity_kin  +360, 360.0) - 360 : fmod(polarization_phi_com_kin  -decay_phi_helicity_kin  +360, 360.0)")
    .Define("psi_helicity_truth",               "fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0)")
    .Define("psi_helicity_diff",                "psi_helicity_kin - psi_helicity_truth")

    .Define("dxs_weight",                       "dxs_weight_func(beam_energy_truth, minust_truth, sim_model_flag)")
    .Define("psi_weight",                       "psi_weight_func(polarization_angle, beam_energy_truth, psi_helicity_truth)")
    .Define("event_weight",                     "beam_accid_weight*combo_accid_weight*dxs_weight*psi_weight")

    .Define("if_kp_in_fdc",                     "event_weight*(kp_dedx_fdc >  0.0 && kp_dedx_cdc == 0.0)")
    .Define("if_kp_in_cdc",                     "event_weight*(kp_dedx_cdc >  0.0 && kp_dedx_fdc == 0.0)")
    .Define("if_kp_in_fdc_cdc",                 "event_weight*(kp_dedx_fdc >  0.0 && kp_dedx_cdc >  0.0)")
    .Define("if_kp_in_neither",                 "event_weight*(kp_dedx_fdc == 0.0 && kp_dedx_cdc == 0.0)")
    .Define("if_km_in_fdc",                     "event_weight*(km_dedx_fdc >  0.0 && km_dedx_cdc == 0.0)")
    .Define("if_km_in_cdc",                     "event_weight*(km_dedx_cdc >  0.0 && km_dedx_fdc == 0.0)")
    .Define("if_km_in_fdc_cdc",                 "event_weight*(km_dedx_fdc >  0.0 && km_dedx_cdc >  0.0)")
    .Define("if_km_in_neither",                 "event_weight*(km_dedx_fdc == 0.0 && km_dedx_cdc == 0.0)")
    .Define("if_best_combo",                    "(combo_accid_weight >  0.0)")
    .Define("if_not_best_combo",                "(combo_accid_weight == 0.0)")
    ;

    cout << "Filtering events...\n";
    string dEdxCut              = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.65*d_momentum_meas+4.47) + 2.57)";
    string dEdxCutSyst          = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.11*d_momentum_meas+3.90) + 1.83)";
    string MissPMinusCut        = "miss_pminus_meas > -0.02";
    string MissPMinusCutSyst    = "miss_pminus_meas > -0.04";
    string KinFitChiSqCut       = "chisq_per_ndf_kin < 5.0";
    string KinFitChiSqCutSyst   = "chisq_per_ndf_kin < 7.0";
    string KinematicsCut        = "kp_momentum_meas > 0.40 && km_momentum_meas > 0.40 && d_momentum_meas > 0.40 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
    string KinematicsCutSyst    = "kp_momentum_meas > 0.35 && km_momentum_meas > 0.35 && d_momentum_meas > 0.35 && kp_theta_meas > 1.0 && km_theta_meas > 1.0 && d_theta_meas > 1.0";
    string VertexCut            = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
    string VertexCutSyst        = "TMath::Abs(vertex_z_kin - 65.0) < 15.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.5";
    string PhiMassCut           = "phi_mass_kin > 1.00 && phi_mass_kin < 1.04";

    auto rdf_NoCut          = rdf_input;
    auto rdf_dEdxCut        = rdf_input                        .Filter(MissPMinusCut.c_str()).Filter(KinFitChiSqCut.c_str()).Filter(KinematicsCut.c_str()).Filter(VertexCut.c_str());
    auto rdf_MissPMinusCut  = rdf_input.Filter(dEdxCut.c_str())                              .Filter(KinFitChiSqCut.c_str()).Filter(KinematicsCut.c_str()).Filter(VertexCut.c_str());
    auto rdf_KinFitChiSqCut = rdf_input.Filter(dEdxCut.c_str()).Filter(MissPMinusCut.c_str())                               .Filter(KinematicsCut.c_str()).Filter(VertexCut.c_str());
    auto rdf_KinematicsCut  = rdf_input.Filter(dEdxCut.c_str()).Filter(MissPMinusCut.c_str()).Filter(KinFitChiSqCut.c_str())                              .Filter(VertexCut.c_str());
    auto rdf_VertexCut      = rdf_input.Filter(dEdxCut.c_str()).Filter(MissPMinusCut.c_str()).Filter(KinFitChiSqCut.c_str()).Filter(KinematicsCut.c_str());
    auto rdf_NominalCut     = rdf_input.Filter(dEdxCut.c_str()).Filter(MissPMinusCut.c_str()).Filter(KinFitChiSqCut.c_str()).Filter(KinematicsCut.c_str()).Filter(VertexCut.c_str());
    auto rdf_PlotCut        = rdf_NominalCut.Filter(PhiMassCut.c_str());
    auto rdf_SystCut        = rdf_input.Filter(dEdxCutSyst.c_str()).Filter(MissPMinusCutSyst.c_str()).Filter(KinFitChiSqCutSyst.c_str()).Filter(KinematicsCutSyst.c_str()).Filter(VertexCutSyst.c_str());
    RNode rdfs []           = {rdf_NoCut,   rdf_dEdxCut,    rdf_MissPMinusCut,  rdf_KinFitChiSqCut, rdf_KinematicsCut,  rdf_VertexCut,  rdf_NominalCut, rdf_PlotCut};
    string labels []        = {"NoCut",     "dEdxCut",      "MissPMinusCut",    "KinFitChiSqCut",   "KinematicsCut",    "VertexCut",    "NominalCut",   "PlotCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_exc_recon_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_phi_d_recon";
        rdf_SystCut.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_recon_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "--Processing " << label << " filter...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

            // cout << "----Processing observable plots..." << endl;
            // TH1D    hist_observable_phi_mass                            = *rdf.Histo1D({("observable_phi_mass_"+ label).c_str(),                ";m^{kin}_{K^{+}K^{-}} (GeV)            ;Counts",                                                           120, 0.98, 1.1},                        "phi_mass_kin",                                                     "event_weight");
            //         hist_observable_phi_mass.Write();
            // TH1D    hist_observable_Eg                                  = *rdf.Histo1D({("observable_Eg_"+ label).c_str(),                      ";E_{beam}^{kin} (GeV)                  ;Counts",                                                           60, 5.0, 11.0},                         "beam_energy_kin",                                                  "event_weight");
            //         hist_observable_Eg.Write();
            // TH1D    hist_observable_minust                              = *rdf.Histo1D({("observable_minust_"+ label).c_str(),                  ";-t^{kin} (GeV^{2})                    ;Counts",                                                           100, 0.0, 2.0},                         "minust_kin",                                                       "event_weight");
            //         hist_observable_minust.Write();
            // TH2D    hist_observable_Eg_minust                           = *rdf.Histo2D({("observable_Eg_minust_"+ label).c_str(),               ";E_{beam}^{kin} (GeV)                  ;-t^{kin} (GeV^{2})",                                               60, 5.0, 11.0, 20, 0.0, 2.0},           "beam_energy_kin",                  "minust_kin",                   "event_weight");
            //         hist_observable_Eg_minust.Write();
            // TH1D    hist_observable_scatter_theta                       = *rdf.Histo1D({("observable_scatter_theta_"+ label).c_str(),           ";#theta_{CM}^{kin} (deg)               ;Counts",                                                           180, 0.0, 180.0},                       "scatter_theta_com_kin",                                            "event_weight");
            //         hist_observable_scatter_theta.Write();
            // TH2D    hist_observable_minust_scatter_theta                = *rdf.Histo2D({("observable_minust_scatter_theta_"+ label).c_str(),    ";-t^{kin} (GeV^{2})                    ;#theta_{CM}^{kin} (deg)",                                          20, 0.0, 2.0, 40, 0.0, 40.0},           "minust_kin",                       "scatter_theta_com_kin",        "event_weight");
            //         hist_observable_minust_scatter_theta.Write();
            // TH2D    hist_observable_minust_d_momentum                   = *rdf.Histo2D({("observable_minust_d_momentum_"+ label).c_str(),       ";-t^{kin} (GeV^{2})                    ;p_{d}^{kin} (GeV)",                                                100, 0.0, 2.0, 200, 0.0, 2.0},          "minust_kin","d_momentum_kin",                                      "event_weight");
            //         hist_observable_minust_d_momentum.Write();
            // TH1D    hist_observable_decay_costheta                      = *rdf.Histo1D({("observable_decay_costheta_"+ label).c_str(),          ";cos(#theta_{helicity}^{kin})          ;Counts",                                                           20, -1.0, 1.0},                         "decay_costheta_helicity_kin",                                      "event_weight");
            //         hist_observable_decay_costheta.Write();
            // TH1D    hist_observable_decay_phi                           = *rdf.Histo1D({("observable_decay_phi_"+ label).c_str(),               ";#phi_{helicity}^{kin} (deg)           ;Counts",                                                           18, -180.0, 180.0},                     "decay_phi_helicity_kin",                                           "event_weight");
            //         hist_observable_decay_phi.Write();
            // TH1D    hist_observable_polarization_phi                    = *rdf.Histo1D({("observable_polarization_phi_"+ label).c_str(),        ";#phi_{com}^{kin} (deg)                ;Counts",                                                           18, -180, 180.0},                       "polarization_phi_com_kin",                                         "event_weight");
            //         hist_observable_polarization_phi.Write();
            // TH1D    hist_observable_psi                                 = *rdf.Histo1D({("observable_psi_"+ label).c_str(),                     ";#psi_{helicity}^{kin} (deg)           ;Counts",                                                           18, -180.0, 180.0},                     "psi_helicity_kin",                                                 "event_weight");
            //         hist_observable_psi.Write();

            // cout << "----Processing KinFit plots..." << endl;
            // TH1D    hist_kinfit_cut_chisq_per_ndf                       = *rdf.Histo1D({("kinfit_cut_chisq_per_ndf_"+ label).c_str(),           ";#chi^{2}/NDF                          ;Counts",                                                           100, 0.0, 10.0},                        "chisq_per_ndf_kin",                                                "event_weight");
            //         hist_kinfit_cut_chisq_per_ndf.Write();
            // TH1D    hist_kinfit_confidence_level                        = *rdf.Histo1D({("kinfit_confidence_level_"+ label).c_str(),            ";log(KinFit FOM)                       ;Counts",                                                           30, -15.0, 0},                          "log10_kinfit_fom_kin",                                             "event_weight");
            //         hist_kinfit_confidence_level.Write();
            // TH2D    hist_kinfit_cl_chisq_per_ndf                        = *rdf.Histo2D({("kinfit_cl_chisq_per_ndf_"+ label).c_str(),            ";log(KinFit FOM)                       ;#chi^{2}/NDF",                                                     150, -15.0, 0.0, 100, 0.0, 10.0},       "log10_kinfit_fom_kin",             "chisq_per_ndf_kin",            "event_weight");
            //         hist_kinfit_cl_chisq_per_ndf.Write();

            // cout << "----Processing PID plots..." << endl;
            // TH2D    hist_pid_cut_d_dEdx_cdc                             = *rdf.Histo2D({("pid_cut_d_dEdx_cdc_"+ label).c_str(),                 ";p_{d}^{meas} (GeV)                    ;(dE/dx)_{d}^{CDC} (keV/cm)",                                       200, 0.0, 2.0, 400, 0.0, 40},           "d_momentum_meas",                  "d_dedx_cdc_keV_per_cm_meas",   "event_weight");
            //         hist_pid_cut_d_dEdx_cdc.Write();
            // TH1D    hist_pid_kp_DeltaT                                  = *rdf.Histo1D({("pid_kp_DeltaT_"+ label).c_str(),                      ";#Delta t^{meas}_{K^{+}} (ns)          ;Counts",                                                           100, -2.0, 2.0},                        "kp_DeltaT_meas",                                                   "event_weight");
            //         hist_pid_kp_DeltaT.Write();
            // TH2D    hist_pid_kp_DeltaT_p                                = *rdf.Histo2D({("pid_kp_DeltaT_p_"+ label).c_str(),                    ";p_{K^{+}}^{meas} (GeV)                ;#Delta t_{K^{+}}^{meas} (ns)",                                     100, 0.0, 10.0, 100, -2.0, 2.0},        "kp_momentum_meas",                 "kp_DeltaT_meas",               "event_weight");
            //         hist_pid_kp_DeltaT_p.Write();
            // TH2D    hist_pid_kp_dEdx_cdc                                = *rdf.Histo2D({("pid_kp_dEdx_cdc_"+ label).c_str(),                    ";p_{K^{+}}^{meas} (GeV)                ;(dE/dx)_{K^{+}}^{CDC} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "kp_momentum_meas",                 "kp_dedx_cdc_keV_per_cm_meas",  "event_weight");
            //         hist_pid_kp_dEdx_cdc.Write();
            // TH2D    hist_pid_kp_dEdx_fdc                                = *rdf.Histo2D({("pid_kp_dEdx_fdc_"+ label).c_str(),                    ";p_{K^{+}}^{meas} (GeV)                ;(dE/dx)_{K^{+}}^{FDC} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "kp_momentum_meas",                 "kp_dedx_fdc_keV_per_cm_meas",  "event_weight");
            //         hist_pid_kp_dEdx_fdc.Write();
            // TH2D    hist_pid_kp_dEdx_tof                                = *rdf.Histo2D({("pid_kp_dEdx_tof_"+ label).c_str(),                    ";p_{K^{+}}^{meas} (GeV)                ;(dE/dx)_{K^{+}}^{TOF} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "kp_momentum_meas",                 "kp_dedx_tof_keV_per_cm_meas",  "event_weight");
            //         hist_pid_kp_dEdx_tof.Write();
            // TH2D    hist_pid_kp_dEdx_st                                 = *rdf.Histo2D({("pid_kp_dEdx_st_"+ label).c_str(),                     ";p_{K^{+}}^{meas} (GeV)                ;(dE/dx)_{K^{+}}^{ST} (keV/cm)",                                    100, 0.0, 10.0, 100, 0.0, 40},          "kp_momentum_meas",                 "kp_dedx_st_keV_per_cm_meas",   "event_weight");
            //         hist_pid_kp_dEdx_st.Write();
            // TH1D    hist_pid_km_DeltaT                                  = *rdf.Histo1D({("pid_km_DeltaT_"+ label).c_str(),                      ";#Delta t^{meas}_{K^{-}} (ns)          ;Counts",                                                           100, -2.0, 2.0},                        "km_DeltaT_meas",                                                   "event_weight");
            //         hist_pid_km_DeltaT.Write();
            // TH2D    hist_pid_km_DeltaT_p                                = *rdf.Histo2D({("pid_km_DeltaT_p_"+ label).c_str(),                    ";p_{K^{-}}^{meas} (GeV)                ;#Delta t_{K^{-}}^{meas} (ns)",                                     100, 0.0, 10.0, 100, -2.0, 2.0},        "km_momentum_meas",                 "km_DeltaT_meas",               "event_weight");
            //         hist_pid_km_DeltaT_p.Write();
            // TH2D    hist_pid_km_dEdx_cdc                                = *rdf.Histo2D({("pid_km_dEdx_cdc_"+ label).c_str(),                    ";p_{K^{-}}^{meas} (GeV)                ;(dE/dx)_{K^{-}}^{CDC} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "km_momentum_meas",                 "km_dedx_cdc_keV_per_cm_meas",  "event_weight");
            //         hist_pid_km_dEdx_cdc.Write();
            // TH2D    hist_pid_km_dEdx_fdc                                = *rdf.Histo2D({("pid_km_dEdx_fdc_"+ label).c_str(),                    ";p_{K^{-}}^{meas} (GeV)                ;(dE/dx)_{K^{-}}^{FDC} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "km_momentum_meas",                 "km_dedx_fdc_keV_per_cm_meas",  "event_weight");
            //         hist_pid_km_dEdx_fdc.Write();
            // TH2D    hist_pid_km_dEdx_tof                                = *rdf.Histo2D({("pid_km_dEdx_tof_"+ label).c_str(),                    ";p_{K^{-}}^{meas} (GeV)                ;(dE/dx)_{K^{-}}^{TOF} (keV/cm)",                                   100, 0.0, 10.0, 100, 0.0, 40},          "km_momentum_meas",                 "km_dedx_tof_keV_per_cm_meas",  "event_weight");
            //         hist_pid_km_dEdx_tof.Write();
            // TH2D    hist_pid_km_dEdx_st                                 = *rdf.Histo2D({("pid_km_dEdx_st_"+ label).c_str(),                     ";p_{K^{-}}^{meas} (GeV)                ;(dE/dx)_{K^{-}}^{ST} (keV/cm)",                                    100, 0.0, 10.0, 100, 0.0, 40},          "km_momentum_meas",                 "km_dedx_st_keV_per_cm_meas",   "event_weight");
            //         hist_pid_km_dEdx_st.Write();
            // TH1D    hist_pid_d_DeltaT                                   = *rdf.Histo1D({("pid_d_DeltaT_"+ label).c_str(),                       ";#Delta t_{d}^{meas} (ns)              ;Counts",                                                           100, -5.0, 5.0},                        "d_DeltaT_meas",                                                    "event_weight");
            //         hist_pid_d_DeltaT.Write();
            // TH2D    hist_pid_d_DeltaT_p                                 = *rdf.Histo2D({("pid_d_DeltaT_p_"+ label).c_str(),                     ";p_{d}^{meas} (GeV)                    ;#Delta t_{d}^{meas} (ns)",                                         200, 0.0, 2.0, 100, -5.0, 5.0},         "d_momentum_meas",                  "d_DeltaT_meas",                "event_weight");
            //         hist_pid_d_DeltaT_p.Write();
            // TH2D    hist_pid_d_dEdx_fdc                                 = *rdf.Histo2D({("pid_d_dEdx_fdc_"+ label).c_str(),                     ";p_{d}^{meas} (GeV)                    ;(dE/dx)_{d}^{FDC} (keV/cm)",                                       200, 0.0, 2.0, 100, 0.0, 40},           "d_momentum_meas",                  "d_dedx_fdc_keV_per_cm_meas",   "event_weight");
            //         hist_pid_d_dEdx_fdc.Write();
            // TH2D    hist_pid_d_dEdx_tof                                 = *rdf.Histo2D({("pid_d_dEdx_tof_"+ label).c_str(),                     ";p_{d}^{meas} (GeV)                    ;(dE/dx)_{d}^{TOF} (keV/cm)",                                       200, 0.0, 2.0, 100, 0.0, 40},           "d_momentum_meas",                  "d_dedx_tof_keV_per_cm_meas",   "event_weight");
            //         hist_pid_d_dEdx_tof.Write();
            // TH2D    hist_pid_d_dEdx_st                                  = *rdf.Histo2D({("pid_d_dEdx_st_"+ label).c_str(),                      ";p_{d}^{meas} (GeV)                    ;(dE/dx)_{d}^{ST} (keV/cm)",                                        200, 0.0, 2.0, 100, 0.0, 40},           "d_momentum_meas",                  "d_dedx_st_keV_per_cm_meas",    "event_weight");
            //         hist_pid_d_dEdx_st.Write();
            // TH2D    hist_pid_d_dEdx_cdc_st                              = *rdf.Histo2D({("pid_d_dEdx_cdc_st_"+ label).c_str(),                  ";(dE/dx)_{d}^{CDC} (keV/cm)            ;(dE/dx)_{d}^{ST} (keV/cm)",                                        100, 0.0, 40, 100, 0.0, 40},            "d_dedx_cdc_keV_per_cm_meas",       "d_dedx_st_keV_per_cm_meas",    "event_weight");
            //         hist_pid_d_dEdx_cdc_st.Write();

            cout << "----Processing exclusivity plots..." << endl;
            TH1D    hist_exclusivity_cut_miss_pminus                    = *rdf.Histo1D({("exclusivity_cut_miss_pminus_"+ label).c_str(),        ";(p_{miss}^{-})^{meas} (GeV)           ;Counts",                                                           400, -0.2, 0.2},                        "miss_pminus_meas",                                                 "event_weight");
                    hist_exclusivity_cut_miss_pminus.Write();
            TH1D    hist_exclusivity_miss_pminus_as_pion                = *rdf.Histo1D({("exclusivity_miss_pminus_as_pion_"+ label).c_str(),    ";(p_{miss}^{-})^{meas} (GeV)           ;Counts",                                                           400, -0.2, 0.2},                        "rho_miss_pminus_meas",                                             "event_weight");
                    hist_exclusivity_miss_pminus_as_pion.Write();
            TH2D    hist_exclusivity_miss_pminus_phi_mass               = *rdf.Histo2D({("exclusivity_miss_pminus_phi_mass_"+ label).c_str(),   ";(p_{miss}^{-})^{meas} (GeV)           ;m_{K^{+}K^{-}}^{meas} (GeV)",                                      400, -0.2, 0.2, 120, 0.98, 1.1},        "miss_pminus_meas",                 "phi_mass_meas",                "event_weight");
                    hist_exclusivity_miss_pminus_phi_mass.Write();
            TH1D    hist_exclusivity_miss_energy                        = *rdf.Histo1D({("exclusivity_miss_energy_"+ label).c_str(),            ";E_{miss}^{meas} (GeV)                 ;Counts",                                                           300, -3.0, 3.0},                        "miss_energy_meas",                                                 "event_weight");
                    hist_exclusivity_miss_energy.Write();
            TH1D    hist_exclusivity_miss_masssquared_meas              = *rdf.Histo1D({("exclusivity_miss_masssquared_"+ label).c_str(),       ";(m_{miss}^{2})^{meas} (GeV^{2})       ;Counts",                                                           100, -0.2, 0.2},                        "miss_masssquared_meas",                                            "event_weight");
                    hist_exclusivity_miss_masssquared_meas.Write();
            TH1D    hist_exclusivity_miss_momentum_meas                 = *rdf.Histo1D({("exclusivity_miss_momentum_"+ label).c_str(),          ";(p_{miss})^{meas} (GeV)               ;Counts",                                                           200, 0.0, 2.0},                         "miss_momentum_meas",                                               "event_weight");
                    hist_exclusivity_miss_momentum_meas.Write();
            TH1D    hist_exclusivity_unused_tracks                      = *rdf.Histo1D({("exclusivity_unused_tracks_"+ label).c_str(),          ";Number of unused tracks               ;Counts",                                                           10, 0.0, 10.0},                         "num_unused_tracks",                                                "event_weight");
                    hist_exclusivity_unused_tracks.Write();
            TH1D    hist_exclusivity_unused_showers                     = *rdf.Histo1D({("exclusivity_unused_showers_"+ label).c_str(),         ";Number of unused showers              ;Counts",                                                           10, 0.0, 10.0},                         "num_unused_showers",                                               "event_weight");
                    hist_exclusivity_unused_showers.Write();
            TH1D    hist_exclusivity_coplanarity_angle                  = *rdf.Histo1D({("exclusivity_coplanarity_angle_"+ label).c_str(),      ";(#phi_{#phi}-#phi_{d})^{meas} (deg)   ;Counts",                                                           41, 159.5, 200.5},                      "coplanarity_meas",                                                 "event_weight");
                    hist_exclusivity_coplanarity_angle.Write();

            cout << "----Processing correlation plots..." << endl;
            TH2D    hist_cut_correlation_chisq_piminus                  = *rdf.Histo2D({("cut_correlation_chisq_piminus_"+ label).c_str(),      ";#chi^{2}/NDF                          ;(p_{miss}^{-})^{meas} (GeV)",                                      100, 0.0, 10.0, 100, -0.2, 0.2},        "chisq_per_ndf_kin",                "miss_pminus_meas",             "event_weight");
                    hist_cut_correlation_chisq_piminus.Write();
            TH2D    hist_cut_correlation_chisq_dedx                     = *rdf.Histo2D({("cut_correlation_chisq_dedx_"+ label).c_str(),         ";#chi^{2}/NDF                          ;dE/dx Residue (keV/cm)",                                           100, 0.0, 10.0, 100, -20.0, 20.0},      "chisq_per_ndf_kin",                "d_dedx_cdc_residue_meas",      "event_weight");
                    hist_cut_correlation_chisq_dedx.Write();
            TH2D    hist_cut_correlation_pminus_dedx                    = *rdf.Histo2D({("cut_correlation_pminus_dedx_"+ label).c_str(),        ";(p_{miss}^{-})^{meas} (GeV)           ;dE/dx Residue (keV/cm)",                                           100, -0.2, 0.2, 100, -20.0, 20.0},      "missing_pminus_meas",              "d_dedx_cdc_residue_meas",      "event_weight");
                    hist_cut_correlation_pminus_dedx.Write();

            cout << "----Processing kinematics plots..." << endl;
            TH2D    hist_kinematics_cut_kp                              = *rdf.Histo2D({("kinematics_cut_kp_"+ label).c_str(),                  ";p_{K^{+}}^{meas} (GeV)                ;#theta_{K^{+}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "kp_momentum_meas",                 "kp_theta_meas",                "event_weight");
                    hist_kinematics_cut_kp.Write();
            TH2D    hist_kinematics_cut_km                              = *rdf.Histo2D({("kinematics_cut_km_"+ label).c_str(),                  ";p_{K^{-}}^{meas} (GeV)                ;#theta_{K^{-}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "km_momentum_meas",                 "km_theta_meas",                "event_weight");
                    hist_kinematics_cut_km.Write();
            TH2D    hist_kinematics_cut_d                               = *rdf.Histo2D({("kinematics_cut_d_"+ label).c_str(),                   ";p_{d}^{meas} (GeV)                    ;#theta_{d}^{meas} (deg)",                                          200, 0.0, 2.0, 180, 0.0, 180.0},        "d_momentum_meas",                  "d_theta_meas",                 "event_weight");
                    hist_kinematics_cut_d.Write();
            TH2D    hist_kinematics_kp_fdc                              = *rdf.Histo2D({("kinematics_kp_fdc_"+ label).c_str(),                  ";p_{K^{+}}^{meas} (GeV)                ;#theta_{K^{+}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "kp_momentum_meas",                 "kp_theta_meas",                "if_kp_in_fdc");
                    hist_kinematics_kp_fdc.Write();
            TH2D    hist_kinematics_kp_fdc_cdc                          = *rdf.Histo2D({("kinematics_kp_fdc_cdc_"+ label).c_str(),              ";p_{K^{+}}^{meas} (GeV)                ;#theta_{K^{+}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "kp_momentum_meas",                 "kp_theta_meas",                "if_kp_in_fdc_cdc");
                    hist_kinematics_kp_fdc_cdc.Write();
            TH2D    hist_kinematics_kp_cdc                              = *rdf.Histo2D({("kinematics_kp_cdc_"+ label).c_str(),                  ";p_{K^{+}}^{meas} (GeV)                ;#theta_{K^{+}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "kp_momentum_meas",                 "kp_theta_meas",                "if_kp_in_cdc");
                    hist_kinematics_kp_cdc.Write();
            TH2D    hist_kinematics_km_fdc                              = *rdf.Histo2D({("kinematics_km_fdc_"+ label).c_str(),                  ";p_{K^{-}}^{meas} (GeV)                ;#theta_{K^{-}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "km_momentum_meas",                 "km_theta_meas",                "if_km_in_fdc");
                    hist_kinematics_km_fdc.Write();
            TH2D    hist_kinematics_km_fdc_cdc                          = *rdf.Histo2D({("kinematics_km_fdc_cdc_"+ label).c_str(),              ";p_{K^{-}}^{meas} (GeV)                ;#theta_{K^{-}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "km_momentum_meas",                 "km_theta_meas",                "if_km_in_fdc_cdc");
                    hist_kinematics_km_fdc_cdc.Write();
            TH2D    hist_kinematics_km_cdc                              = *rdf.Histo2D({("kinematics_km_cdc_"+ label).c_str(),                  ";p_{K^{-}}^{meas} (GeV)                ;#theta_{K^{-}}^{meas} (deg)",                                      100, 0.0, 10.0, 180, 0.0, 180.0},       "km_momentum_meas",                 "km_theta_meas",                "if_km_in_cdc");
                    hist_kinematics_km_cdc.Write();
            TH2D    hist_kinematics_phi                                 = *rdf.Histo2D({("kinematics_phi_"+ label).c_str(),                     ";p_{#phi}^{kin} (GeV)                  ;#theta_{#phi}^{kin} (deg)",                                        110, 0.0, 11.0, 180, 0.0, 180.0},       "phi_momentum_kin",                 "phi_theta_kin",                "event_weight");
                    hist_kinematics_phi.Write();

            cout << "----Processing vertex plots..." << endl;
            TH1D    hist_vertex_cut_z                                   = *rdf.Histo1D({("vertex_cut_z_"+ label).c_str(),                       ";Z_{vertex}^{kin} (cm)                 ;Counts",                                                           100, 40.0, 90.0},                       "vertex_z_kin",                                                     "event_weight");
                    hist_vertex_cut_z.Write();
            TH2D    hist_vertex_cut_x_y                                 = *rdf.Histo2D({("vertex_cut_x_y_"+ label).c_str(),                     ";X_{vertex}^{kin} (cm)                 ;Y_{vertex}^{kin} (cm)",                                            100, -2.0, 2.0, 100, -2.0, 2.0},        "vertex_x_kin",                     "vertex_y_kin",                 "event_weight");
                    hist_vertex_cut_x_y.Write();

            cout << "----Processing accidental plots..." << endl;
            TH1D    hist_accidental_beam_weight                         = *rdf.Histo1D({("accidental_beam_weight_"+ label).c_str(),             ";Beam Accidental Weight                ;Counts",                                                           20, -0.5, 1.5},                         "beam_accid_weight");
                    hist_accidental_beam_weight.Write();
            TH1D    hist_accidental_beam_DeltaT                         = *rdf.Histo1D({("accidental_beam_DeltaT_"+ label).c_str(),             ";#Delta t_{beam}^ {meas} (ns)          ;Counts",                                                           1000, -25.0, 25.0},                     "beam_DeltaT_meas",                                                 "beam_accid_weight");
                    hist_accidental_beam_DeltaT.Write();
            TH1D    hist_accidental_combo_weight                        = *rdf.Histo1D({("accidental_combo_weight_"+ label).c_str(),            ";Combo Accidental Weight               ;Counts",                                                           2, -0.5, 1.5},                          "combo_accid_weight");
                    hist_accidental_combo_weight.Write();
            TH2D    hist_accidental_combo_kinematics_best               = *rdf.Histo2D({("accidental_combo_kinematics_best_"+ label).c_str(),   ";#theta_{K^{+}}^{meas} (deg)           ;#theta_{K^{-}}^{meas} (deg)",                                      200, 0.0, 20.0, 200, 0.0, 20.0},        "kp_theta_meas",                    "km_theta_meas",                "if_best_combo");
                    hist_accidental_combo_kinematics_best.Write();
            TH2D    hist_accidental_combo_kinematics_others             = *rdf.Histo2D({("accidental_combo_kinematics_others_"+ label).c_str(), ";#theta_{K^{+}}^{meas} (deg)           ;#theta_{K^{-}}^{meas} (deg)",                                      200, 0.0, 20.0, 200, 0.0, 20.0},        "kp_theta_meas",                    "km_theta_meas",                "if_not_best_combo");
                    hist_accidental_combo_kinematics_others.Write();

            if (reaction.find("sim") != string::npos)
            {
                cout << "----Processing truth observable plots..." << endl;
                TH1D    hist_truth_observable_phi_mass                  = *rdf.Histo1D({("truth_observable_phi_mass_"+ label).c_str(),          ";m_{K^{+}K^{-}}^{truth} (GeV)          ;Counts",                                                           500, 0.9, 1.9},                         "phi_mass_truth",                                                   "event_weight");
                        hist_truth_observable_phi_mass.Write();
                TH1D    hist_truth_observable_Eg                        = *rdf.Histo1D({("truth_observable_Eg_"+ label).c_str(),                ";E_{beam}^{truth} (GeV)                ;Counts",                                                           60, 5.0, 11.0},                         "beam_energy_truth",                                                "event_weight");
                        hist_truth_observable_Eg.Write();
                TH1D    hist_truth_observable_minust                    = *rdf.Histo1D({("truth_observable_minust_"+ label).c_str(),            ";-t^{truth} (GeV^{2})                  ;Counts",                                                           100, 0.0, 2.0},                         "minust_truth",                                                     "event_weight");
                        hist_truth_observable_minust.Write();
                TH2D    hist_truth_observable_Eg_minust                 = *rdf.Histo2D({("truth_observable_Eg_minust_"+ label).c_str(),         ";E_{beam}^{truth} (GeV)                ;-t^{truth} (GeV^{2})",                                             60, 5.0, 11.0, 20, 0.0, 2.0},           "beam_energy_truth",                "minust_truth",                 "event_weight");
                        hist_truth_observable_Eg_minust.Write();
                TH1D    hist_truth_observable_scatter_theta             = *rdf.Histo1D({("truth_observable_scatter_theta_"+ label).c_str(),     ";#theta_{CM}^{truth} (deg)             ;Counts",                                                           180, 0.0, 180.0},                       "scatter_theta_com_truth",                                          "event_weight");
                        hist_truth_observable_scatter_theta.Write();
                TH1D    hist_truth_observable_decay_costheta            = *rdf.Histo1D({("truth_observable_decay_costheta_"+ label).c_str(),    ";cos(#theta_{helicity}^{truth})        ;Counts",                                                           10, -1.0, 1.0},                         "decay_costheta_helicity_truth",                                    "event_weight");
                        hist_truth_observable_decay_costheta.Write();
                TH1D    hist_truth_observable_decay_phi                 = *rdf.Histo1D({("truth_observable_decay_phi_"+ label).c_str(),         ";#phi_{helicity}^{truth} (deg)         ;Counts",                                                           10, -180.0, 180.0},                     "decay_phi_helicity_truth",                                         "event_weight");
                        hist_truth_observable_decay_phi.Write();
                TH1D    hist_truth_observable_polarization_phi          = *rdf.Histo1D({("truth_observable_polarization_phi_"+ label).c_str(),  ";#phi_{com}^{truth} (deg)              ;Counts",                                                           10, -180, 180.0},                       "polarization_phi_com_truth",                                       "event_weight");
                        hist_truth_observable_polarization_phi.Write();
                TH1D    hist_truth_observable_psi                       = *rdf.Histo1D({("truth_observable_psi_"+ label).c_str(),               ";#psi_{helicity}^{truth} (deg)         ;Counts",                                                           10, -180.0, 180.0},                     "psi_helicity_truth",                                               "event_weight");
                        hist_truth_observable_psi.Write();

                cout << "----Processing truth kinematics plots..." << endl;
                TH2D    hist_truth_kinematics_kp                        = *rdf.Histo2D({("truth_kinematics_kp_"+ label).c_str(),                ";p_{K^{+}}^{truth} (GeV)               ;#theta_{K^{+}}^{truth} (deg)",                                     100, 0.0, 10.0, 180, 0.0, 180.0},       "kp_momentum_truth",                "kp_theta_truth",               "event_weight");
                        hist_truth_kinematics_kp.Write();
                TH2D    hist_truth_kinematics_km                        = *rdf.Histo2D({("truth_kinematics_km_"+ label).c_str(),                ";p_{K^{-}}^{truth} (GeV)               ;#theta_{K^{-}}^{truth} (deg)",                                     100, 0.0, 10.0, 180, 0.0, 180.0},       "km_momentum_truth",                "km_theta_truth",               "event_weight");
                        hist_truth_kinematics_km.Write();
                TH2D    hist_truth_kinematics_d                         = *rdf.Histo2D({("truth_kinematics_d_"+ label).c_str(),                 ";p_{d}^{truth} (GeV)                   ;#theta_{d}^{truth} (deg)",                                         100, 0.0, 10.0, 180, 0.0, 180.0},       "d_momentum_truth",                 "d_theta_truth",                "event_weight");
                        hist_truth_kinematics_d.Write();
                TH2D    hist_truth_kinematics_phi                       = *rdf.Histo2D({("truth_kinematics_phi_"+ label).c_str(),               ";p_{#phi}^{truth} (GeV)                ;#theta_{#phi}^{truth} (deg)",                                      100, 0.0, 11.0, 180, 0.0, 180.0},       "phi_momentum_truth",               "phi_theta_truth",              "event_weight");
                        hist_truth_kinematics_phi.Write();
                TH2D    hist_truth_kinematics_phi_d_theta               = *rdf.Histo2D({("truth_kinematics_phi_d_theta_"+ label).c_str(),       ";#theta_{d}^{truth} (deg)              ;#theta_{#phi}^{truth} (deg)",                                      180, 0.0, 180.0, 180, 0.0, 180.0},      "d_theta_truth",                    "phi_theta_truth",              "event_weight");
                        hist_truth_kinematics_phi_d_theta.Write();
                TH2D    hist_truth_kinematics_phi_d_p                   = *rdf.Histo2D({("truth_kinematics_phi_d_p_"+ label).c_str(),           ";p_{d}^{truth} (GeV)                   ;p_{#phi}^{truth} (GeV)",                                           100, 0.0, 10.0, 100, 0.0, 10.0},        "d_momentum_truth",                 "phi_momentum_truth",           "event_weight");
                        hist_truth_kinematics_phi_d_p.Write();

                cout << "----Processing bin migration plots..." << endl;
                TH2D    hist_bin_migration_Eg                           = *rdf.Histo2D({("bin_migration_Eg_"+ label).c_str(),                   ";E_{beam}^{truth} (GeV)                ;E_{beam}^{kin} (GeV)",                                             60, 5.0, 11.0, 60, 5.0, 11.0},          "beam_energy_truth",                "beam_energy_kin",              "event_weight");
                        hist_bin_migration_Eg.Write();
                TH2D    hist_bin_migration_minust                       = *rdf.Histo2D({("bin_migration_minust_"+ label).c_str(),               ";-t^{truth} (GeV^{2})                  ;-t^{kin} (GeV^{2})",                                               20, 0.0, 2.0, 20, 0.0, 2.0},            "minust_truth",                     "minust_kin",                   "event_weight");
                        hist_bin_migration_minust.Write();
                TH2D    hist_bin_migration_decay_costheta               = *rdf.Histo2D({("bin_migration_decay_costheta_"+ label).c_str(),       ";cos(#theta_{helicity}^{truth})        ;cos(#theta_{helicity}^{kin})",                                     10, -1.0, 1.0, 10, -1.0, 1.0},          "decay_costheta_helicity_truth",    "decay_costheta_helicity_kin",  "event_weight");
                        hist_bin_migration_decay_costheta.Write();
                TH2D    hist_bin_migration_decay_phi                    = *rdf.Histo2D({("bin_migration_decay_phi_"+ label).c_str(),            ";#phi_{helicity}^{truth} (deg)         ;#phi_{helicity}^{kin} (deg)",                                      10, -180.0, 180.0, 10, -180.0, 180.0},  "decay_phi_helicity_truth",         "decay_phi_helicity_kin",       "event_weight");
                        hist_bin_migration_decay_phi.Write();
                TH2D    hist_bin_migration_polarization_phi             = *rdf.Histo2D({("bin_migration_polarization_phi_"+ label).c_str(),     ";#phi_{com}^{truth} (deg)              ;#phi_{com}^{kin} (deg)",                                           10, -180.0, 180.0, 10, -180.0, 180.0},  "polarization_phi_com_truth",       "polarization_phi_com_kin",     "event_weight");
                        hist_bin_migration_polarization_phi.Write();
                TH2D    hist_bin_migration_psi                          = *rdf.Histo2D({("bin_migration_psi_"+ label).c_str(),                  ";#psi_{helicity}^{truth} (deg)         ;#psi_{helicity}^{kin} (deg)",                                      10, -180.0, 180.0, 10, -180.0, 180.0},  "psi_helicity_truth",               "psi_helicity_kin",             "event_weight");
                        hist_bin_migration_psi.Write();

                cout << "----Processing resolution plots..." << endl;
                TH2D    hist_resolution_Eg                              = *rdf.Histo2D({("resolution_Eg_"+ label).c_str(),                      ";E_{beam}^{truth} (GeV)                ;E_{beam}^{kin} - E_{beam}^{truth} (GeV)",                          60, 5.0, 11.0, 100, -0.5, 0.5},         "beam_energy_truth",                "beam_energy_diff",             "event_weight");
                        hist_resolution_Eg.Write();
                TH2D    hist_resolution_minust                          = *rdf.Histo2D({("resolution_minust_"+ label).c_str(),                  ";-t^{truth} (GeV^{2})                  ;-t^{kin} - -t^{truth} (GeV^{2})",                                  40, 0.0, 2.0, 80, -0.2, 0.2},           "minust_truth",                     "minust_diff",                  "event_weight");
                        hist_resolution_minust.Write();
                TH2D    hist_resolution_decay_costheta                  = *rdf.Histo2D({("resolution_decay_costheta_"+ label).c_str(),          ";cos(#theta_{helicity}^{truth})        ;cos(#theta_{helicity}^{kin}) - cos(#theta_{helicity}^{truth})",    40, 0.0, 2.0, 80, -1.0, 1.0},           "decay_costheta_helicity_truth",    "decay_costheta_helicity_diff", "event_weight");
                        hist_resolution_decay_costheta.Write();
                TH2D    hist_resolution_decay_phi                       = *rdf.Histo2D({("resolution_decay_phi_"+ label).c_str(),               ";#phi_{helicity}^{truth} (deg)         ;#phi_{helicity}^{kin} - #phi_{helicity}^{truth} (deg)",            40, 0.0, 2.0, 80, -8.0, 8.0},           "decay_phi_helicity_truth",         "decay_phi_helicity_diff",      "event_weight");
                        hist_resolution_decay_phi.Write();
                TH2D    hist_resolution_polarization_phi                = *rdf.Histo2D({("resolution_polarization_phi_"+ label).c_str(),        ";#phi_{com}^{truth} (deg)              ;#phi_{com}^{kin} - #phi_{com}^{truth} (deg)",                      40, 0.0, 2.0, 80, -2.0, 2.0},           "polarization_phi_com_truth",       "polarization_phi_com_diff",    "event_weight");
                        hist_resolution_polarization_phi.Write();
                TH2D    hist_resolution_psi                             = *rdf.Histo2D({("resolution_psi_"+ label).c_str(),                     ";#psi_{helicity}^{truth} (deg)         ;#psi_{helicity}^{kin} - #psi_{helicity}^{truth} (deg)",            40, 0.0, 2.0, 80, -8.0, 8.0},           "psi_helicity_truth",               "psi_helicity_diff",            "event_weight");
                        hist_resolution_psi.Write();

                cout << "----Processing resolution w.r.t. t plots..." << endl;
                TH2D    hist_resolution_t_phi_mass                      = *rdf.Histo2D({("resolution_t_phi_mass_"+ label).c_str(),              ";-t^{truth} (GeV^{2})                  ;m_{K^{+}K^{-}}^{kin} - m_{K^{+}K^{-}}^{truth} (GeV)",              40, 0.0, 2.0, 80, -0.02, 0.02},         "minust_truth",                     "phi_mass_diff",                "event_weight");
                        hist_resolution_t_phi_mass.Write();
                TH2D    hist_resolution_t_decay_costheta                = *rdf.Histo2D({("resolution_t_decay_costheta_"+ label).c_str(),        ";-t^{truth} (GeV^{2})                  ;cos(#theta_{helicity}^{kin}) - cos(#theta_{helicity}^{truth})",    40, 0.0, 2.0, 80, -1.0, 1.0},           "minust_truth",                     "decay_costheta_helicity_diff", "event_weight");
                        hist_resolution_t_decay_costheta.Write();
                TH2D    hist_resolution_t_decay_phi                     = *rdf.Histo2D({("resolution_t_decay_phi_"+ label).c_str(),             ";-t^{truth} (GeV^{2})                  ;#phi_{helicity}^{kin} - #phi_{helicity}^{truth} (deg)",            40, 0.0, 2.0, 80, -8.0, 8.0},           "minust_truth",                     "decay_phi_helicity_diff",      "event_weight");
                        hist_resolution_t_decay_phi.Write();
                TH2D    hist_resolution_t_polarization_phi              = *rdf.Histo2D({("resolution_t_polarization_phi_"+ label).c_str(),      ";-t^{truth} (GeV^{2})                  ;#phi_{com}^{kin} - #phi_{com}^{truth} (deg)",                      40, 0.0, 2.0, 80, -2.0, 2.0},           "minust_truth",                     "polarization_phi_com_diff",    "event_weight");
                        hist_resolution_t_polarization_phi.Write();
                TH2D    hist_resolution_t_psi                           = *rdf.Histo2D({("resolution_t_psi_"+ label).c_str(),                   ";-t^{truth} (GeV^{2})                  ;#psi_{helicity}^{kin} - #psi_{helicity}^{truth} (deg)",            40, 0.0, 2.0, 80, -8.0, 8.0},           "minust_truth",                     "psi_helicity_diff",            "event_weight");
                        hist_resolution_t_psi.Write();
            }
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}