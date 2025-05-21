#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;
double mass_missing = 0.0;

void filter_rho_d_thrown(string reaction, string output_mode)
{
    cout << "Reading input files...\n";
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_rho_d_thrown_%s.root",reaction.c_str());
    string input_tree_name  = "selectedtree_rho_d_thrown";
    if (reaction.find("gen") != string::npos)
        input_tree_name = "genT";
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
    if (reaction.find("gen") != string::npos)
    {
        rdf_def = rdf_def
        .Define("beam_p4_truth",            "pBeam")
        // .Define("pip_p4_truth",              "pDecay1")
        .Define("pip_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        // .Define("pim_p4_truth",              "pDecay2")
        .Define("pim_p4_truth",              "TLorentzVector(0, 0, 0, 0)")
        .Define("d_p4_truth",               "pBaryon")
        .Define("polarization_angle",       "-1")
        ;
    }
    auto rdf_input = rdf_def
    .Define("beam_energy_truth",                "beam_p4_truth.E()")
    .Define("pip_as_kaon_p4_truth",             "TLorentzVector(pip_p4_truth.Vect(), TMath::Sqrt(pip_p4_truth.P()*pip_p4_truth.P() + mass_kplus*mass_kplus))")
    .Define("pip_energy_truth",                 "pip_p4_truth.E()")
    .Define("pip_momentum_truth",               "pip_p4_truth.P()")
    .Define("pip_theta_truth",                  "pip_p4_truth.Theta()*RadToDeg")
    .Define("pim_as_kaon_p4_truth",             "TLorentzVector(pim_p4_truth.Vect(), TMath::Sqrt(pim_p4_truth.P()*pim_p4_truth.P() + mass_kminus*mass_kminus))")
    .Define("pim_energy_truth",                 "pim_p4_truth.E()")
    .Define("pim_momentum_truth",               "pim_p4_truth.P()")
    .Define("pim_theta_truth",                  "pim_p4_truth.Theta()*RadToDeg")
    .Define("d_energy_truth",                   "d_p4_truth.E()")
    .Define("d_momentum_truth",                 "d_p4_truth.P()")
    .Define("d_theta_truth",                    "d_p4_truth.Theta()*RadToDeg")
    .Define("rho_p4_truth",                     "pip_p4_truth + pim_p4_truth")
    .Define("rho_energy_truth",                 "rho_p4_truth.E()")
    .Define("rho_momentum_truth",               "rho_p4_truth.P()")
    .Define("rho_mass_truth",                   "rho_p4_truth.M()")
    .Define("rho_proxymass_truth",              "TMath::Sqrt(rho_p4_truth.Minus()*(2*beam_energy_truth + mass_target - d_p4_truth.Plus() - (TMath::Sq(mass_missing) + (rho_p4_truth + d_p4_truth).Perp2())/(mass_target - (rho_p4_truth + d_p4_truth).Minus())) - rho_p4_truth.Perp2())")
    .Define("rho_theta_truth",                  "rho_p4_truth.Theta()*RadToDeg")
    .Define("struck_p4_truth",                  "rho_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("struck_energy_truth",              "struck_p4_truth.E()")
    .Define("struck_energy_balance_truth",      "struck_energy_truth - mass_2H")
    .Define("struck_mass_truth",                "struck_p4_truth.M()")
    .Define("struck_masssquared_truth",         "struck_p4_truth.M2()")
    .Define("struck_pminus_truth",              "struck_p4_truth.Minus()")
    .Define("struck_momentum_truth",            "struck_p4_truth.P()")
    .Define("struck_theta_truth",               "struck_p4_truth.Theta()*RadToDeg")
    .Define("miss_p4_truth",                    "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target) - rho_p4_truth - d_p4_truth")
    .Define("miss_energy_truth",                "miss_p4_truth.E()")
    .Define("miss_energy_balance_truth",        "miss_energy_truth - mass_missing")
    .Define("miss_mass_truth",                  "miss_p4_truth.M()")
    .Define("miss_mass_balance_truth",          "miss_mass_truth - mass_missing")
    .Define("miss_masssquared_truth",           "miss_p4_truth.M2()")
    .Define("miss_masssquared_balance_truth",   "miss_masssquared_truth - mass_missing*mass_missing")
    .Define("miss_pminus_truth",                "miss_p4_truth.Minus()")
    .Define("miss_pminus_balance_truth",        "miss_pminus_truth - mass_missing")
    .Define("miss_momentum_truth",              "miss_p4_truth.P()")
    .Define("miss_theta_truth",                 "miss_p4_truth.Theta()*RadToDeg")
    .Define("total_initial_p4_truth",           "beam_p4_truth + TLorentzVector(0, 0, 0, mass_target)")
    .Define("total_final_p4_truth",             "rho_p4_truth + d_p4_truth")
    .Define("sqrts_truth",                      "total_final_p4_truth.Mag()")
    .Define("minust_truth",                     "-(beam_p4_truth - rho_p4_truth).Mag2()")
    .Define("minusu_truth",                     "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("coplanarity_truth",                "abs(rho_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("y_rho_truth",                      "minust_truth/(2*mass_2H*(beam_p4_truth.E()-rho_p4_truth.E()))")
    .Define("phi_mass_truth",                   "(pip_as_kaon_p4_truth + pim_as_kaon_p4_truth).M()")
    .Define("epsilon_x3_com",                   "TVector3(TMath::Cos(polarization_angle/RadToDeg), TMath::Sin(polarization_angle/RadToDeg), 0)")
    .Define("beam_p4_com_truth",                "boost_lorentz_vector(beam_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("rho_p4_com_truth",                 "boost_lorentz_vector(rho_p4_truth, -total_final_p4_truth.BoostVector())")
    .Define("z_x3_com_truth",                   "beam_p4_com_truth.Vect().Unit()")
    .Define("y_x3_com_truth",                   "beam_p4_com_truth.Vect().Cross(rho_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_com_truth",                   "y_x3_com_truth.Cross(z_x3_com_truth).Unit()")
    .Define("polarization_phi_com_truth",       "TMath::ATan2(-x_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)), y_x3_com_truth.Dot(epsilon_x3_com.Cross(z_x3_com_truth)))*RadToDeg")
    .Define("rho_theta_com_truth",              "rho_p4_com_truth.Vect().Angle(z_x3_com_truth)*RadToDeg")
    .Define("z_x3_helicity_truth",              "rho_p4_com_truth.Vect().Unit()")
    .Define("y_x3_helicity_truth",              "beam_p4_com_truth.Vect().Cross(rho_p4_com_truth.Vect()).Unit()")
    .Define("x_x3_helicity_truth",              "y_x3_helicity_truth.Cross(z_x3_helicity_truth).Unit()")
    .Define("pi_x3_helicity_truth",             "boost_lorentz_vector(pip_p4_truth, -rho_p4_truth.BoostVector()).Vect().Unit()")
    .Define("pip_costheta_helicity_truth",      "pi_x3_helicity_truth.Dot(z_x3_helicity_truth)")
    .Define("pip_phi_helicity_truth",           "TMath::ATan2(-x_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)), y_x3_helicity_truth.Dot(pi_x3_helicity_truth.Cross(z_x3_helicity_truth)))*RadToDeg")
    .Define("psi_helicity_truth",               "fmod(polarization_phi_com_truth-pip_phi_helicity_truth+360, 360.0) >= 180 ? fmod(polarization_phi_com_truth-pip_phi_helicity_truth+360, 360.0) - 360 : fmod(polarization_phi_com_truth-pip_phi_helicity_truth+360, 360.0)")
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
        string output_treefile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_rho_d_thrown_%s.root",reaction.c_str());
        string output_tree_name = "filteredtree_rho_d_thrown";
        rdf_output.Snapshot(output_tree_name.c_str(), output_treefile_name.c_str());
    }

    if (output_mode == "hist" || output_mode == "both")
    {
        cout << "Plotting histograms...\n";
        string output_histfile_name = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_rho_d_thrown_%s.root",reaction.c_str());
        TFile * output_histfile = new TFile(output_histfile_name.c_str(), "RECREATE");
        output_histfile->cd();

        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";
            TDirectory * dir = output_histfile->mkdir(label.c_str());
            dir->cd();

                TH1D hist_beam_energy_truth                 = *rdf.Histo1D({("beam_energy_truth_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_truth");
                hist_beam_energy_truth.Write();
                TH2D hist_pip_kinematics_truth              = *rdf.Histo2D({("pip_kinematics_truth_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_truth","pip_theta_truth");
                hist_pip_kinematics_truth.Write();
                TH2D hist_pim_kinematics_truth              = *rdf.Histo2D({("pim_kinematics_truth_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_truth","pim_theta_truth");
                hist_pim_kinematics_truth.Write();
                TH2D hist_d_kinematics_truth                = *rdf.Histo2D({("d_kinematics_truth_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_truth","d_theta_truth");
                hist_d_kinematics_truth.Write();
                TH1D hist_rho_mass_truth                    = *rdf.Histo1D({("rho_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c);Counts", 500, 0.0, 2.0},"rho_mass_truth");
                hist_rho_mass_truth.Write();
                TH2D hist_rho_kinematics_truth              = *rdf.Histo2D({("rho_kinematics_truth_"+ label).c_str(), ";P_{#rho} (GeV/c);#theta_{#rho} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"rho_momentum_truth","rho_theta_truth");
                hist_rho_kinematics_truth.Write();
                TH2D hist_rho_d_theta_truth                 = *rdf.Histo2D({("rho_d_theta_truth_"+ label).c_str(), ";#theta_{d} (deg);#theta_{#rho} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"d_theta_truth","rho_theta_truth");
                hist_rho_d_theta_truth.Write();
                TH2D hist_rho_d_momentum_truth              = *rdf.Histo2D({("rho_d_momentum_truth_"+ label).c_str(), ";P_{d} (GeV/c);P_{#rho} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"d_momentum_truth","rho_momentum_truth");
                hist_rho_d_momentum_truth.Write();

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
                TH1D hist_yrho_truth                        = *rdf.Histo1D({("yrho_truth_"+ label).c_str(), ";y_{#rho};Counts", 200, 0.0, 2.0},"y_rho_truth");
                hist_yrho_truth.Write();
                TH1D hist_phi_mass_truth                    = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c^{2});Counts", 100, 0.9, 1.9},"phi_mass_truth");
                hist_rho_mass_truth.Write();
                TH2D hist_beam_energy_minust_truth          = *rdf.Histo2D({("beam_energy_minust_truth_"+ label).c_str(), ";E_{beam} (GeV);-t (GeV^{2}/c^{2})", 60, 5.0, 11.0, 30, 0.0, 3.0},"beam_energy_truth","minust_truth");
                hist_beam_energy_minust_truth.Write();

                TH1D hist_rho_theta_com_truth               = *rdf.Histo1D({("rho_theta_com_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"rho_theta_com_truth");
                hist_rho_theta_com_truth.Write();
                TH2D hist_minust_rho_theta_com_truth        = *rdf.Histo2D({("minust_rho_theta_com_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 30, 0.0, 3.0, 180, 0.0, 180.0},"minust_truth","rho_theta_com_truth");
                hist_minust_rho_theta_com_truth.Write();
                TH1D hist_polarization_rho_com_truth        = *rdf.Histo1D({("polarization_rho_com_truth_"+ label).c_str(), ";#rho_{com} (deg);Counts", 9, -180, 180.0},"polarization_rho_com_truth");
                hist_polarization_rho_com_truth.Write();

                TH1D hist_pip_costheta_helicity_truth       = *rdf.Histo1D({("pip_costheta_helicity_truth_"+ label).c_str(), ";cos(#theta_{helicity});Counts", 10, -1.0, 1.0},"pip_costheta_helicity_truth");
                hist_pip_costheta_helicity_truth.Write();
                TH1D hist_pip_rho_helicity_truth            = *rdf.Histo1D({("pip_rho_helicity_truth_"+ label).c_str(), ";#rho_{helicity} (deg);Counts", 9, -180.0, 180.0},"pip_rho_helicity_truth");
                hist_pip_rho_helicity_truth.Write();
                TH1D hist_psi_helicity_truth                = *rdf.Histo1D({("psi_helicity_truth_"+ label).c_str(), ";#psi_{helicity} (deg);Counts", 9, -270.0, 270.0},"psi_helicity_truth");
                hist_psi_helicity_truth.Write();
        }
        output_histfile->Close();
    }
    cout << "Done!\n";
}