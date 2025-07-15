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

void filter_phi_d_thrown_exc(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_thrown_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_phi_d_thrown";
    if (reaction.find("gen") != string::npos)
        input_tree_name = "genT";
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());

    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    if (reaction.find("gen") != string::npos)
    {
        rdf_def = rdf_def
        .Define("beam_p4_truth",            "pBeam")
        .Define("kp_p4_truth",              "pDecay1")
        // .Define("kp_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        .Define("km_p4_truth",              "pDecay2")
        // .Define("km_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        .Define("d_p4_truth",               "pBaryon")
        .Define("polarization_angle",       "-1")
        ;
    }
    auto rdf_input = rdf_def
    .Define("target_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")
    .Define("sim_weight",                       "sim_weight_func(beam_p4_truth.E(), -(target_p4 - d_p4_truth).Mag2())")
    .Define("event_weight",                     "sim_weight")
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
    .Define("coplanarity_truth",                "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("rho_mass_truth",                   "(kp_as_pion_p4_truth + km_as_pion_p4_truth).M()")
    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_p4_truth.BoostVector())")
    .Define("phi_p4_com_truth",                 "boost_lorentz_vector(phi_p4_truth, -total_p4_truth.BoostVector())")
    .Define("kp_p4_com_truth",                  "boost_lorentz_vector(kp_p4_truth, -total_p4_truth.BoostVector())")
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
    ;

    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_output         = rdf_NoCut;
    RNode rdfs []           = {rdf_NoCut};
    string labels []        = {"NoCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_thrown_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_phi_d_thrown";
        rdf_output.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_thrown_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

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
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}