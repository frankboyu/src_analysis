#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double get_polarization_angle(int run_number)
{
    int list_0[]    = { 90209, 90213, 90218, 90219, 90223, 90229, 90234, 90240, 90244, 90249, 90559, 90561, 90563, 90564, 90567,
                        90569, 90571, 90572, 90574, 90576, 90577, 90583, 90586, 90587, 90589, 90595, 90597, 90599, 90601};
    int list_45[]   = { 90211, 90215, 90221, 90226, 90231, 90242, 90247};
    int list_90[]   = { 90210, 90214, 90220, 90225, 90230, 90241, 90245, 90246, 90558, 90560, 90562, 90565, 90568, 90570, 90573,
                        90575, 90579, 90580, 90581, 90582, 90584, 90585, 90588, 90590, 90591, 90592, 90593, 90596, 90598, 90600};
    int list_135[]  = { 90208, 90212, 90217, 90222, 90227, 90232, 90233, 90243, 90248};
    int list_amo[]  = { 90228, 90578};
    double this_polarization_angle = 999;
    if      (std::find(std::begin(list_0),   std::end(list_0),   run_number) != std::end(list_0))
        this_polarization_angle = 0.0;
    else if (std::find(std::begin(list_45),  std::end(list_45),  run_number) != std::end(list_45))
        this_polarization_angle = 45.0;
    else if (std::find(std::begin(list_90),  std::end(list_90),  run_number) != std::end(list_90))
        this_polarization_angle = 90.0;
    else if (std::find(std::begin(list_135), std::end(list_135), run_number) != std::end(list_135))
        this_polarization_angle = 135.0;
    else if (std::find(std::begin(list_amo), std::end(list_amo), run_number) != std::end(list_amo))
        this_polarization_angle = -1.0;
    else
        cout << "Run number not found in the production list.";
    // cou << run_number << " is assigned with polarization angle = " << this_polarization_angle << " degree.\n";
    return this_polarization_angle;
}

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

void filter_phi_d_exc_thrown(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_exc_thrown_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_phi_d_exc_thrown";
    if (reaction.find("gen") != string::npos)
        input_tree_name = "genT";
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());

    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    if (reaction.find("model") != string::npos)
        rdf_def = rdf_def.Define("sim_model_flag", "1");
    else
        rdf_def = rdf_def.Define("sim_model_flag", "-1");
    if (reaction.find("gen") != string::npos)
    {
        rdf_def = rdf_def
        .Define("beam_p4_truth",            "pBeam")
        .Define("kp_p4_truth",              "pDecay1")
        .Define("km_p4_truth",              "pDecay2")
        .Define("d_p4_truth",               "pBaryon")
        .Define("polarization_angle",       "0.0")
        .Define("beam_x4_truth",            "TLorentzVector(0, 0, 0, 0)")
        ;
    }
    auto rdf_input = rdf_def
    .Define("target_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("kp_as_pion_p4_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_energy_truth",                  "kp_p4_truth.E()")
    .Define("kp_momentum_truth",                "kp_p4_truth.P()")
    .Define("kp_theta_truth",                   "kp_p4_truth.Theta()*RadToDeg")
    .Define("km_as_pion_p4_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_energy_truth",                  "km_p4_truth.E()")
    .Define("km_momentum_truth",                "km_p4_truth.P()")
    .Define("km_theta_truth",                   "km_p4_truth.Theta()*RadToDeg")
    .Define("d_energy_truth",                   "d_p4_truth.E()")
    .Define("d_momentum_truth",                 "d_p4_truth.P()")
    .Define("d_theta_truth",                    "d_p4_truth.Theta()*RadToDeg")
    .Define("phi_p4_truth",                     "kp_p4_truth + km_p4_truth")
    .Define("phi_energy_truth",                 "phi_p4_truth.E()")
    .Define("phi_momentum_truth",               "phi_p4_truth.P()")
    .Define("phi_mass_truth",                   "phi_p4_truth.M()")
    .Define("phi_theta_truth",                  "phi_p4_truth.Theta()*RadToDeg")
    .Define("miss_p4_truth",                    "beam_p4_truth + target_p4 - phi_p4_truth - d_p4_truth")
    .Define("miss_energy_truth",                "miss_p4_truth.E()")
    .Define("miss_masssquared_truth",           "miss_p4_truth.M2()")
    .Define("miss_pminus_truth",                "miss_p4_truth.Minus()")
    .Define("miss_momentum_truth",              "miss_p4_truth.P()")
    .Define("total_p4_truth",                   "phi_p4_truth + d_p4_truth")
    .Define("minust_truth",                     "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minusu_truth",                     "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("vertex_z_truth",                   "beam_x4_truth.Z()")
    .Define("vertex_x_truth",                   "beam_x4_truth.X()")
    .Define("vertex_y_truth",                   "beam_x4_truth.Y()")
    .Define("coplanarity_truth",                "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("rho_mass_truth",                   "(kp_as_pion_p4_truth + km_as_pion_p4_truth).M()")
    .Define("rho_miss_pminus_truth",            "(beam_p4_truth + target_p4 - kp_as_pion_p4_truth - km_as_pion_p4_truth - d_p4_truth).Minus()")
    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_p4_truth     .BoostVector())")
    .Define("phi_p4_com_truth",                 "boost_lorentz_vector(phi_p4_truth,  -total_p4_truth     .BoostVector())")
    .Define("kp_p4_com_truth",                  "boost_lorentz_vector(kp_p4_truth,   -total_p4_truth     .BoostVector())")
    .Define("z_x3_com_truth",                   "beam_p4_com_truth.Vect().Unit()")
    .Define("y_x3_com_truth",                   "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_com_truth",                   "y_x3_com_truth.Cross(z_x3_com_truth).Unit()")
    .Define("polarization_phi_com_truth",       "TMath::ATan2(-x_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)), y_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)))*RadToDeg")
    .Define("scatter_theta_com_truth",          "phi_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")
    .Define("z_x3_helicity_truth",              "phi_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(phi_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_truth",             "boost_lorentz_vector(kp_p4_com_truth, -phi_p4_com_truth.BoostVector()).Vect().Unit()")
    .Define("decay_costheta_helicity_truth",    "pi_x3_helicity_truth.Dot(z_x3_helicity_truth)")
    .Define("decay_phi_helicity_truth",         "TMath::ATan2(-x_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)), y_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)))*RadToDeg")
    .Define("psi_helicity_truth",               "fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth-decay_phi_helicity_truth+360, 360.0)")
    .Define("dxs_weight",                       "dxs_weight_func(beam_energy_truth, minust_truth, sim_model_flag)")
    .Define("psi_weight",                       "psi_weight_func(polarization_angle, beam_energy_truth, psi_helicity_truth)")
    .Define("event_weight",                     "dxs_weight*psi_weight")
    ;

    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_SystCut        = rdf_NoCut;
    RNode rdfs []           = {rdf_NoCut};
    string labels []        = {"NoCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_exc_thrown_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_phi_d_exc_thrown";
        rdf_SystCut.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_thrown_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "--Processing " << label << " filter...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

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
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}