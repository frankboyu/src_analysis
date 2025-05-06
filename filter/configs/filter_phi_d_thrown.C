#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;

void filter_phi_d_thrown(string reaction_name, string output_mode)
{
    // Define input and output names
    string input_tree;
    if (reaction_name.find("tagged") != string::npos)
        input_tree = "selectedtree_phi_d_thrown";
    else
        input_tree = "genT";
    string output_tree      = "filteredtree_phi_d_thrown";
    string input_tree_file  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_phi_d_thrown_%s.root",reaction_name.c_str());
    string output_hist_file = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_thrown_%s.root",reaction_name.c_str());
    string output_tree_file = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_thrown_%s.root",reaction_name.c_str());

    // Determine reaction specific parameters
    if      (reaction_name.find("2H")   != string::npos)
        mass_target = mass_2H;
    else if (reaction_name.find("4He")  != string::npos)
        mass_target = mass_4He;
    else if (reaction_name.find("12C")  != string::npos)
        mass_target = mass_12C;

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(input_tree.c_str());
    chain.Add(input_tree_file.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    if (reaction_name.find("gen") != string::npos)
    {
        rdf_def = rdf_def
        .Define("beam_p4_truth",            "pBeam")
        .Define("kp_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        .Define("km_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        .Define("d_p4_truth",               "pBaryon")
        ;
    }

    auto rdf_input = rdf_def
    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(kp_p4_truth + km_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")
    .Define("kp_p4pion_truth",              "TLorentzVector(kp_p4_truth.Vect(), TMath::Sqrt(kp_p4_truth.P()*kp_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("kp_p4helicity_truth",          "boost_lorentz_vector(kp_p4_truth, -(kp_p4_truth + km_p4_truth).BoostVector())")
    .Define("kp_energy_truth",              "kp_p4_truth.E()")
    .Define("kp_momentum_truth",            "kp_p4_truth.P()")
    .Define("kp_theta_truth",               "kp_p4_truth.Theta()*RadToDeg")
    .Define("km_p4pion_truth",              "TLorentzVector(km_p4_truth.Vect(), TMath::Sqrt(km_p4_truth.P()*km_p4_truth.P() + mass_piminus*mass_piminus))")
    .Define("km_p4helicity_truth",          "boost_lorentz_vector(km_p4_truth, -(kp_p4_truth + km_p4_truth).BoostVector())")
    .Define("km_energy_truth",              "km_p4_truth.E()")
    .Define("km_momentum_truth",            "km_p4_truth.P()")
    .Define("km_theta_truth",               "km_p4_truth.Theta()*RadToDeg")
    .Define("d_energy_truth",               "d_p4_truth.E()")
    .Define("d_momentum_truth",             "d_p4_truth.P()")
    .Define("d_theta_truth",                "d_p4_truth.Theta()*RadToDeg")
    .Define("phi_p4_truth",                 "kp_p4_truth + km_p4_truth")
    .Define("phi_p4com_truth",              "boost_lorentz_vector(phi_p4_truth, -(phi_p4_truth + d_p4_truth).BoostVector())")
    .Define("phi_energy_truth",             "phi_p4_truth.E()")
    .Define("phi_momentum_truth",           "phi_p4_truth.P()")
    .Define("phi_mass_truth",               "phi_p4_truth.M()")
    .Define("phi_proxymass_truth",          "TMath::Sqrt((kp_p4_truth.Minus() + km_p4_truth.Minus())*(2*beam_p4_truth.E() + mass_4He - d_p4_truth.Plus() - (mass_2H*mass_2H + (phi_p4_truth + d_p4_truth).Perp2())/(mass_4He - (phi_p4_truth + d_p4_truth).Minus())) - (kp_p4_truth.Px() + km_p4_truth.Px())*(kp_p4_truth.Px() + km_p4_truth.Px()) - (kp_p4_truth.Py() + km_p4_truth.Py())*(kp_p4_truth.Py() + km_p4_truth.Py()))")
    .Define("phi_theta_truth",              "phi_p4_truth.Theta()*RadToDeg")
    .Define("struck_p4_truth",              "phi_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("struck_energy_truth",          "struck_p4_truth.E()")
    .Define("struck_mass_truth",            "struck_p4_truth.M()")
    .Define("struck_masssquared_truth",     "struck_p4_truth.M2()")
    .Define("struck_momentum_truth",        "struck_p4_truth.P()")
    .Define("struck_pminus_truth",          "struck_p4_truth.Minus()")
    .Define("struck_theta_truth",           "struck_p4_truth.Theta()*RadToDeg")
    .Define("struck_energy_balance_truth",  "struck_energy_truth - mass_2H")
    .Define("miss_p4_truth",                "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target) - phi_p4_truth - d_p4_truth")
    .Define("miss_energy_truth",            "miss_p4_truth.E()")
    .Define("miss_mass_truth",              "miss_p4_truth.M()")
    .Define("miss_masssquared_truth",       "miss_p4_truth.M2()")
    .Define("miss_momentum_truth",          "miss_p4_truth.P()")
    .Define("miss_pminus_truth",            "miss_p4_truth.Minus()")
    .Define("miss_theta_truth",             "miss_p4_truth.Theta()*RadToDeg")
    .Define("miss_energy_balance_truth",    "miss_energy_truth - (mass_target - mass_2H)")
    .Define("sqrts_truth",                  "(phi_p4_truth + d_p4_truth).Mag()")
    .Define("minust_truth",                 "-(beam_p4_truth - phi_p4_truth).Mag2()")
    .Define("minusu_truth",                 "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("coplanarity_truth",            "abs(phi_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(phi_p4com_truth.Vect())*RadToDeg")
    .Define("y_phi_truth",                  "minust_truth/(2*mass_2H*(beam_p4_truth.E()-phi_p4_truth.E()))")
    .Define("rho_mass_truth",               "(kp_p4pion_truth + km_p4pion_truth).M()")
    .Define("pi_helicity_truth",            "kp_p4helicity_truth.Vect().Unit()")
    .Define("z_helicity_truth",             "phi_p4com_truth.Vect().Unit()")
    .Define("y_helicity_truth",             "beam_p4com_truth.Vect().Cross(phi_p4com_truth.Vect()).Unit()")
    .Define("x_helicity_truth",             "y_helicity_truth.Cross(z_helicity_truth).Unit()")
    .Define("costheta_helicity_truth",      "pi_helicity_truth.Dot(z_helicity_truth)")
    .Define("phi_helicity_truth",           "TMath::ATan2(-x_helicity_truth.Dot(pi_helicity_truth.Cross(z_helicity_truth)), y_helicity_truth.Dot(pi_helicity_truth.Cross(z_helicity_truth)))*RadToDeg")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_output         = rdf_NoCut;
    RNode rdfs []           = {rdf_NoCut};
    string labels []        = {"NoCut"};
    int N_filters           = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot(output_tree.c_str(),output_tree_file.c_str());
    }

    // Save histograms
    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        TFile * hist_file = new TFile(output_hist_file.c_str(), "RECREATE");
        hist_file->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = hist_file->mkdir(label.c_str());
            dir->cd();

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

                TH1D hist_struck_mass_truth                 = *rdf.Histo1D({("struck_mass_truth_"+ label).c_str(), ";m_{struck} (GeV/c^{2});Counts", 400, -4.0, 4.0},"struck_mass_truth");
                hist_struck_mass_truth.Write();
                TH1D hist_struck_masssquared_truth          = *rdf.Histo1D({("struck_masssquared_truth_"+ label).c_str(), ";m_{struck}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"struck_masssquared_truth");
                hist_struck_masssquared_truth.Write();
                TH1D hist_struck_momentum_truth             = *rdf.Histo1D({("struck_momentum_truth_"+ label).c_str(), ";P_{struck} (GeV/c);Counts", 200, 0.0, 2.0},"struck_momentum_truth");
                hist_struck_momentum_truth.Write();
                TH1D hist_struck_pminus_truth               = *rdf.Histo1D({("struck_pminus_truth_"+ label).c_str(), ";P_{struck}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"struck_pminus_truth");
                hist_struck_pminus_truth.Write();
                TH1D hist_struck_energy_balance_truth       = *rdf.Histo1D({("struck_energy_balance_truth_"+ label).c_str(), ";E_{struck} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"struck_energy_balance_truth");
                hist_struck_energy_balance_truth.Write();
                TH2D hist_struck_momentum_pminus_truth      = *rdf.Histo2D({("struck_momentum_pminus_truth_"+ label).c_str(), ";P_{struck} (GeV/c);P_{struck}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"struck_momentum_truth","struck_pminus_truth");
                hist_struck_momentum_pminus_truth.Write();
                TH2D hist_struck_momentum_energy_truth      = *rdf.Histo2D({("struck_momentum_energy_truth_"+ label).c_str(), ";P_{struck} (GeV/c);E_{struck} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"struck_momentum_truth","struck_energy_truth");
                hist_struck_momentum_energy_truth.Write();

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
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}