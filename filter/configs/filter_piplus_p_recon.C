#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

void filter_piplus_p_recon(string reaction_name, string output_mode)
{
    string input_name   = Form("/work/halld2/home/boyu/src_analysis/selection/output/test/selectedtree_piplus_p_recon_%s.root",reaction_name.c_str());
    string hist_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredhist_piplus_p_recon_%s.root",reaction_name.c_str());
    string tree_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/test/filteredtree_piplus_p_recon_%s.root",reaction_name.c_str());

    // Determine reaction specific parameters
    double mass_target = 0.0;
    if (reaction_name.find("2H") != string::npos)
        mass_target = mass_2H;
    else if (reaction_name.find("4He") != string::npos)
        mass_target = mass_4He;
    else if (reaction_name.find("12C") != string::npos)
        mass_target = mass_12C;

    // Read input files
    cout << "Reading input files...\n";
    TChain chain("selectedtree_piplus_p_recon");
    chain.Add(input_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_input = rdf_raw
    .Define("kinfit_fom",                   "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4",                    "TLorentzVector(0, 0, 0, mass_target)")
    .Define("N2_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_p4com_meas",              "boost_lorentz_vector(beam_p4_meas, -(pip_p4_meas + p_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",               "boost_lorentz_vector(beam_p4_kin, -(pip_p4_kin + p_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(pip_p4_truth + p_p4_truth).BoostVector())")
    .Define("beam_energy_meas",             "beam_p4_meas.E()")
    .Define("beam_energy_kin",              "beam_p4_kin.E()")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",             "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",              "beam_x4_kin.T() - rftime")
    .Define("beam_DeltaT_truth",            "beam_x4_truth.T() - rftime")

    .Define("pip_p4com_meas",               "boost_lorentz_vector(pip_p4_meas, -(pip_p4_meas + p_p4_meas).BoostVector())")
    .Define("pip_p4com_kin",                "boost_lorentz_vector(pip_p4_kin, -(pip_p4_meas + p_p4_meas).BoostVector())")
    .Define("pip_p4com_truth",              "boost_lorentz_vector(pip_p4_truth, -(pip_p4_truth + p_p4_truth).BoostVector())")
    .Define("pip_energy_meas",              "pip_p4_meas.E()")
    .Define("pip_energy_kin",               "pip_p4_kin.E()")
    .Define("pip_energy_truth",             "pip_p4_truth.E()")
    .Define("pip_momentum_meas",            "pip_p4_meas.P()")
    .Define("pip_momentum_kin",             "pip_p4_kin.P()")
    .Define("pip_momentum_truth",           "pip_p4_truth.P()")
    .Define("pip_theta_meas",               "pip_p4_meas.Theta()*RadToDeg")
    .Define("pip_theta_kin",                "pip_p4_kin.Theta()*RadToDeg")
    .Define("pip_theta_truth",              "pip_p4_truth.Theta()*RadToDeg")
    .Define("pip_in_fdc",                   "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_in_cdc",                   "accidweight*(pip_dedx_cdc > 0.0 && pip_dedx_fdc == 0.0)")
    .Define("pip_in_fdc_cdc",               "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc > 0.0)")
    .Define("pip_in_neither",               "accidweight*(pip_dedx_fdc == 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_DeltaT_meas",              "rftime + (pip_x4_meas.Z()-65.0)/29.9792458 - pip_x4_meas.T()")
    .Define("pip_DeltaT_kin",               "rftime + (pip_x4_kin.Z()-65.0)/29.9792458 - pip_x4_kin.T()")
    .Define("pip_DeltaT_truth",             "rftime + (pip_x4_truth.Z()-65.0)/29.9792458 - pip_x4_truth.T()")
    .Define("pip_dedx_fdc_keV_per_cm",      "pip_dedx_fdc*1e6")
    .Define("pip_dedx_cdc_keV_per_cm",      "pip_dedx_cdc*1e6")
    .Define("pip_dedx_st_keV_per_cm",       "pip_dedx_st*1e6")
    .Define("pip_dedx_tof_keV_per_cm",      "pip_dedx_tof*1e6")

    .Define("p_p4pion_meas",                "TLorentzVector(p_p4_meas.Vect(), TMath::Sqrt(p_p4_meas.P()*p_p4_meas.P() + mass_piplus*mass_piplus))")
    .Define("p_p4pion_kin",                 "TLorentzVector(p_p4_kin.Vect(), TMath::Sqrt(p_p4_kin.P()*p_p4_kin.P() + mass_piplus*mass_piplus))")
    .Define("p_p4pion_truth",               "TLorentzVector(p_p4_truth.Vect(), TMath::Sqrt(p_p4_truth.P()*p_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("p_energy_meas",                "p_p4_meas.E()")
    .Define("p_energy_kin",                 "p_p4_kin.E()")
    .Define("p_energy_truth",               "p_p4_truth.E()")
    .Define("p_momentum_meas",              "p_p4_meas.P()")
    .Define("p_momentum_kin",               "p_p4_kin.P()")
    .Define("p_momentum_truth",             "p_p4_truth.P()")
    .Define("p_theta_meas",                 "p_p4_meas.Theta()*RadToDeg")
    .Define("p_theta_kin",                  "p_p4_kin.Theta()*RadToDeg")
    .Define("p_theta_truth",                "p_p4_truth.Theta()*RadToDeg")
    .Define("p_in_fdc",                     "accidweight*(p_dedx_fdc > 0.0 && p_dedx_cdc == 0.0)")
    .Define("p_in_cdc",                     "accidweight*(p_dedx_cdc > 0.0 && p_dedx_fdc == 0.0)")
    .Define("p_in_fdc_cdc",                 "accidweight*(p_dedx_fdc > 0.0 && p_dedx_cdc > 0.0)")
    .Define("p_in_neither",                 "accidweight*(p_dedx_fdc == 0.0 && p_dedx_cdc == 0.0)")
    .Define("p_DeltaT_meas",                "rftime + (p_x4_meas.Z()-65.0)/29.9792458 - p_x4_meas.T()")
    .Define("p_DeltaT_kin",                 "rftime + (p_x4_kin.Z()-65.0)/29.9792458 - p_x4_kin.T()")
    .Define("p_DeltaT_truth",               "rftime + (p_x4_truth.Z()-65.0)/29.9792458 - p_x4_truth.T()")
    .Define("p_dedx_fdc_keV_per_cm",        "p_dedx_fdc*1e6")
    .Define("p_dedx_cdc_keV_per_cm",        "p_dedx_cdc*1e6")
    .Define("p_dedx_st_keV_per_cm",         "p_dedx_st*1e6")
    .Define("p_dedx_tof_keV_per_cm",        "p_dedx_tof*1e6")

    .Define("n_p4_meas",                    "pip_p4_meas + p_p4_meas - beam_p4_meas")
    .Define("n_p4_kin",                     "pip_p4_kin + p_p4_kin - beam_p4_kin")
    .Define("n_p4_truth",                   "pip_p4_truth + p_p4_truth - beam_p4_truth")
    .Define("n_energy_meas",                "n_p4_meas.E()")
    .Define("n_energy_kin",                 "n_p4_kin.E()")
    .Define("n_energy_truth",               "n_p4_truth.E()")
    .Define("n_mass_meas",                  "n_p4_meas.M()")
    .Define("n_mass_kin",                   "n_p4_kin.M()")
    .Define("n_mass_truth",                 "n_p4_truth.M()")
    .Define("n_masssquared_meas",           "n_p4_meas.M2()")
    .Define("n_masssquared_kin",            "n_p4_kin.M2()")
    .Define("n_masssquared_truth",          "n_p4_truth.M2()")
    .Define("n_momentum_meas",              "n_p4_meas.P()")
    .Define("n_momentum_kin",               "n_p4_kin.P()")
    .Define("n_momentum_truth",             "n_p4_truth.P()")
    .Define("n_pminus_meas",                "n_p4_meas.Minus()")
    .Define("n_pminus_kin",                 "n_p4_kin.Minus()")
    .Define("n_pminus_truth",               "n_p4_truth.Minus()")
    .Define("n_energy_balance_meas",        "n_p4_meas.E() - mass_neutron")
    .Define("n_energy_balance_kin",         "n_p4_kin.E() - mass_neutron")
    .Define("n_energy_balance_truth",       "n_p4_truth.E() - mass_neutron")

    .Define("miss_p4_meas",                 "beam_p4_meas + target_p4 - pip_p4_meas - p_p4_meas")
    .Define("miss_p4_kin",                  "beam_p4_kin + target_p4 - pip_p4_kin - p_p4_kin")
    .Define("miss_p4_truth",                "beam_p4_truth + target_p4 - pip_p4_truth - p_p4_truth")
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

    .Define("N2miss_p4_meas",               "beam_p4_meas + N2_p4 - pip_p4_meas - p_p4_meas")
    .Define("N2miss_p4_kin",                "beam_p4_kin + N2_p4 - pip_p4_kin - p_p4_kin")
    .Define("N2miss_p4_truth",              "beam_p4_truth + N2_p4 - pip_p4_truth - p_p4_truth")
    .Define("N2miss_energy_meas",           "N2miss_p4_meas.E()")
    .Define("N2miss_energy_kin",            "N2miss_p4_kin.E()")
    .Define("N2miss_energy_truth",          "N2miss_p4_truth.E()")
    .Define("N2miss_mass_meas",             "N2miss_p4_meas.M()")
    .Define("N2miss_mass_kin",              "N2miss_p4_kin.M()")
    .Define("N2miss_mass_truth",            "N2miss_p4_truth.M()")
    .Define("N2miss_masssquared_meas",      "N2miss_p4_meas.M2()")
    .Define("N2miss_masssquared_kin",       "N2miss_p4_kin.M2()")
    .Define("N2miss_masssquared_truth",     "N2miss_p4_truth.M2()")
    .Define("N2miss_momentum_meas",         "N2miss_p4_meas.P()")
    .Define("N2miss_momentum_kin",          "N2miss_p4_kin.P()")
    .Define("N2miss_momentum_truth",        "N2miss_p4_truth.P()")
    .Define("N2miss_pminus_meas",           "N2miss_p4_meas.Minus()")
    .Define("N2miss_pminus_kin",            "N2miss_p4_kin.Minus()")
    .Define("N2miss_pminus_truth",          "N2miss_p4_truth.Minus()")
    .Define("N2miss_energy_balance_meas",   "N2miss_p4_meas.E() - mass_proton")
    .Define("N2miss_energy_balance_kin",    "N2miss_p4_kin.E() - mass_proton")
    .Define("N2miss_energy_balance_truth",  "N2miss_p4_truth.E() - mass_proton")

    .Define("sqrts_meas",                   "(pip_p4_meas + p_p4_meas).Mag()")
    .Define("sqrts_kin",                    "(pip_p4_kin + p_p4_kin).Mag()")
    .Define("sqrts_truth",                  "(pip_p4_truth + p_p4_truth).Mag()")
    .Define("minust_meas",                  "-(pip_p4_meas - beam_p4_meas).Mag2()")
    .Define("minust_kin",                   "-(pip_p4_kin - beam_p4_kin).Mag2()")
    .Define("minust_truth",                 "-(pip_p4_truth - beam_p4_truth).Mag2()")
    .Define("minusu_meas",                  "-(p_p4_meas - beam_p4_meas).Mag2()")
    .Define("minusu_kin",                   "-(p_p4_kin - beam_p4_kin).Mag2()")
    .Define("minusu_truth",                 "-(p_p4_truth - beam_p4_truth).Mag2()")
    .Define("coplanarity_meas",             "abs(pip_p4_meas.Phi() - p_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin",              "abs(pip_p4_kin.Phi() - p_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth",            "abs(pip_p4_truth.Phi() - p_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_meas",                 "beam_p4com_meas.Vect().Angle(pip_p4com_meas.Vect())*RadToDeg")
    .Define("thetaCM_kin",                  "beam_p4com_kin.Vect().Angle(pip_p4com_kin.Vect())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(pip_p4com_truth.Vect())*RadToDeg")
    .Define("rho_misid_mass_meas",          "(pip_p4_meas + p_p4pion_meas).M()")
    .Define("rho_misid_mass_kin",           "(pip_p4_kin + p_p4pion_kin).M()")
    .Define("rho_misid_mass_truth",         "(pip_p4_truth + p_p4pion_truth).M()")
    .Define("rho_misspion_mass_meas",       "(pip_p4_meas + TLorentzVector(0, 0, 0, mass_piplus)).M()")
    .Define("rho_misspion_mass_kin",        "(pip_p4_kin + TLorentzVector(0, 0, 0, mass_piplus)).M()")
    .Define("rho_misspion_mass_truth",      "(pip_p4_truth + TLorentzVector(0, 0, 0, mass_piplus)).M()")

    .Define("topology_pimprot",             "accidweight*(thrown_topology == 0)")
    .Define("topology_pippimprot",          "accidweight*(thrown_topology == 1)")
    .Define("topology_pippimn",             "accidweight*(thrown_topology == 2)")
    .Define("topology_pi0pimprot",          "accidweight*(thrown_topology == 3)")
    .Define("topology_others",              "accidweight*(thrown_topology == 9)")
    .Define("extrapion_momentum_truth",     "extrapion_p4_truth.P()")
    .Define("extrapion_theta_truth",        "extrapion_p4_truth.Theta()*RadToDeg")
    .Define("rho_thrown_p4_truth",          "extrapion_p4_truth + pip_p4_truth")
    .Define("rho_thrown_mass_truth",        "rho_thrown_p4_truth.M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    string miss_p_cut;
    if (reaction_name.find("inc") != string::npos)
        miss_p_cut = "0.50";
    else if (reaction_name.find("missprot") != string::npos)
        miss_p_cut = "0.20";
    else if (reaction_name.find("misshe3") != string::npos)
        miss_p_cut = "0.25";
    else if (reaction_name.find("missb11") != string::npos)
        miss_p_cut = "0.30";

    auto rdf_NoCut          = rdf_input;
    auto rdf_KinematicsCut  = rdf_NoCut.Filter("(minust_kin > 0.2) && (minusu_kin > 0.2)");
    auto rdf_KinFitFOMCut   = rdf_KinematicsCut.Filter("kinfit_fom > 0.01");
    auto rdf_PIDFOMCut      = rdf_KinFitFOMCut.Filter("(pip_pidfom > 0.01) && (p_pidfom > 0.01)");
    auto rdf_MissPCut       = rdf_PIDFOMCut.Filter("(n_momentum_kin < " + miss_p_cut + ")");
    auto rdf_MissPMinusCut  = rdf_MissPCut.Filter("(n_pminus_kin > 0.5) && (n_pminus_kin < 1.3)");
    auto rdf_output         = rdf_MissPMinusCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_KinematicsCut,  rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_MissPCut,   rdf_MissPMinusCut};
    string labels []    = {"NoCut",     "KinematicsCut",    "KinFitFOMCut",     "PIDFOMCut",    "MissPCut",     "MissPminusCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_piplus_p_recon",tree_name.c_str());
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

            TH1D hist_pip_pidfom                    = *rdf.Histo1D({("pip_pidfom_"+ label).c_str(), ";pip_pidfom;Counts", 100, 0.0, 1.0},"pip_pidfom","accidweight");
            hist_pip_pidfom.Write();
            TH2D hist_pip_dEdx_cdc_kin              = *rdf.Histo2D({("pip_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_cdc_keV_per_cm","accidweight");
            hist_pip_dEdx_cdc_kin.Write();
            TH2D hist_pip_dEdx_fdc_kin              = *rdf.Histo2D({("pip_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_fdc_keV_per_cm","accidweight");
            hist_pip_dEdx_fdc_kin.Write();
            TH2D hist_pip_dEdx_tof_kin              = *rdf.Histo2D({("pip_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_tof_keV_per_cm","accidweight");
            hist_pip_dEdx_tof_kin.Write();
            TH2D hist_pip_dEdx_st_kin               = *rdf.Histo2D({("pip_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_st_keV_per_cm","accidweight");
            hist_pip_dEdx_st_kin.Write();
            TH1D hist_pip_DeltaT_kin                = *rdf.Histo1D({("pip_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{+}} (ns);Counts", 400, -20.0, 20.0},"pip_DeltaT_kin","accidweight");
            hist_pip_DeltaT_kin.Write();
            TH2D hist_pip_kinematics_kin            = *rdf.Histo2D({("pip_kinematics_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","accidweight");
            hist_pip_kinematics_kin.Write();
            TH2D hist_pip_kinematics_fdc_kin        = *rdf.Histo2D({("pip_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc");
            hist_pip_kinematics_fdc_kin.Write();
            TH2D hist_pip_kinematics_fdc_cdc_kin    = *rdf.Histo2D({("pip_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc_cdc");
            hist_pip_kinematics_fdc_cdc_kin.Write();
            TH2D hist_pip_kinematics_cdc_kin        = *rdf.Histo2D({("pip_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_cdc");
            hist_pip_kinematics_cdc_kin.Write();
            TH2D hist_pip_kinematics_neither_kin    = *rdf.Histo2D({("pip_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_neither");
            hist_pip_kinematics_neither_kin.Write();

            TH1D hist_p_pidfom                      = *rdf.Histo1D({("p_pidfom_"+ label).c_str(), ";p_pidfom;Counts", 100, 0.0, 1.0},"p_pidfom","accidweight");
            hist_p_pidfom.Write();
            TH2D hist_p_dEdx_cdc_kin                = *rdf.Histo2D({("p_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"p_momentum_kin","p_dedx_cdc_keV_per_cm","accidweight");
            hist_p_dEdx_cdc_kin.Write();
            TH2D hist_p_dEdx_fdc_kin                = *rdf.Histo2D({("p_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"p_momentum_kin","p_dedx_fdc_keV_per_cm","accidweight");
            hist_p_dEdx_fdc_kin.Write();
            TH2D hist_p_dEdx_tof_kin                = *rdf.Histo2D({("p_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"p_momentum_kin","p_dedx_tof_keV_per_cm","accidweight");
            hist_p_dEdx_tof_kin.Write();
            TH2D hist_p_dEdx_st_kin                 = *rdf.Histo2D({("p_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 400, 0.0, 40},"p_momentum_kin","p_dedx_st_keV_per_cm","accidweight");
            hist_p_dEdx_st_kin.Write();
            TH1D hist_p_DeltaT_kin                  = *rdf.Histo1D({("p_DeltaT_kin_"+ label).c_str(), ";#Delta t_{p} (ns);Counts", 400, -20.0, 20.0},"p_DeltaT_kin","accidweight");
            hist_p_DeltaT_kin.Write();
            TH2D hist_p_kinematics_kin              = *rdf.Histo2D({("p_kinematics_kin_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_kin","p_theta_kin","accidweight");
            hist_p_kinematics_kin.Write();
            TH2D hist_p_kinematics_fdc_kin          = *rdf.Histo2D({("p_kinematics_fdc_kin_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_kin","p_theta_kin","p_in_fdc");
            hist_p_kinematics_fdc_kin.Write();
            TH2D hist_p_kinematics_fdc_cdc_kin      = *rdf.Histo2D({("p_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_kin","p_theta_kin","p_in_fdc_cdc");
            hist_p_kinematics_fdc_cdc_kin.Write();
            TH2D hist_p_kinematics_cdc_kin          = *rdf.Histo2D({("p_kinematics_cdc_kin_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_kin","p_theta_kin","p_in_cdc");
            hist_p_kinematics_cdc_kin.Write();
            TH2D hist_p_kinematics_neither_kin      = *rdf.Histo2D({("p_kinematics_neither_kin_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_kin","p_theta_kin","p_in_neither");
            hist_p_kinematics_neither_kin.Write();

            TH2D hist_pip_p_theta_kin               = *rdf.Histo2D({("pip_p_theta_kin_"+ label).c_str(), ";#theta_{#pi^{+}} (deg);#theta_{p} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"pip_theta_kin","p_theta_kin","accidweight");
            hist_pip_p_theta_kin.Write();
            TH2D hist_pip_p_momentum_kin            = *rdf.Histo2D({("pip_p_momentum_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);P_{p} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"pip_momentum_kin","p_momentum_kin","accidweight");
            hist_pip_p_momentum_kin.Write();

            TH1D hist_sqrts_kin                     = *rdf.Histo1D({("sqrts_kin_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_kin","accidweight");
            hist_sqrts_kin.Write();
            TH1D hist_minust_kin                    = *rdf.Histo1D({("minust_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minust_kin","accidweight");
            hist_minust_kin.Write();
            TH1D hist_coplanarity_kin               = *rdf.Histo1D({("coplanarity_kin_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_kin","accidweight");
            hist_coplanarity_kin.Write();
            TH1D hist_thetaCM_kin                   = *rdf.Histo1D({("thetaCM_kin_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_kin","accidweight");
            hist_thetaCM_kin.Write();
            TH2D hist_minust_thetaCM_kin            = *rdf.Histo2D({("minust_thetaCM_kin_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_kin","thetaCM_kin","accidweight");
            hist_minust_thetaCM_kin.Write();
            TH1D hist_rho_misid_mass_kin            = *rdf.Histo1D({("rho_misid_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misid_mass_kin","accidweight");
            hist_rho_misid_mass_kin.Write();

            TH1D hist_n_mass_kin                    = *rdf.Histo1D({("n_mass_kin_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 400, -4.0, 4.0},"n_mass_kin","accidweight");
            hist_n_mass_kin.Write();
            TH1D hist_n_masssquared_kin             = *rdf.Histo1D({("n_masssquared_kin_"+ label).c_str(), ";m_{n}^{2} (GeV^{2}/c^{4});Counts", 400, 0.0, 4.0},"n_masssquared_kin","accidweight");
            hist_n_masssquared_kin.Write();
            TH1D hist_n_momentum_kin                = *rdf.Histo1D({("n_momentum_kin_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 1.0},"n_momentum_kin","accidweight");
            hist_n_momentum_kin.Write();
            TH1D hist_n_pminus_kin                  = *rdf.Histo1D({("n_pminus_kin_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_kin","accidweight");
            hist_n_pminus_kin.Write();
            TH1D hist_n_energy_balance_kin          = *rdf.Histo1D({("n_energy_balance_kin_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"n_energy_balance_kin","accidweight");
            hist_n_energy_balance_kin.Write();
            TH2D hist_n_momentum_pminus_kin         = *rdf.Histo2D({("n_momentum_pminus_kin_"+ label).c_str(), ";P_{n} (GeV/c);P_{n}^{-} (GeV/c)", 100, 0.0, 1.0, 100, 0.4, 1.4},"n_momentum_kin","n_pminus_kin","accidweight");
            hist_n_momentum_pminus_kin.Write();
            TH1D hist_n_mass_meas                   = *rdf.Histo1D({("n_mass_meas_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 400, -4.0, 4.0},"n_mass_meas","accidweight");
            hist_n_mass_meas.Write();
            TH1D hist_n_masssquared_meas            = *rdf.Histo1D({("n_masssquared_meas_"+ label).c_str(), ";m_{n}^{2} (GeV^{2}/c^{4});Counts", 400, 0.0, 4.0},"n_masssquared_meas","accidweight");
            hist_n_masssquared_meas.Write();
            TH1D hist_n_momentum_meas               = *rdf.Histo1D({("n_momentum_meas_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 1.0},"n_momentum_meas","accidweight");
            hist_n_momentum_meas.Write();
            TH1D hist_n_pminus_meas                 = *rdf.Histo1D({("n_pminus_meas_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_meas","accidweight");
            hist_n_pminus_meas.Write();
            TH1D hist_n_energy_balance_meas         = *rdf.Histo1D({("n_energy_balance_meas_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"n_energy_balance_meas","accidweight");
            hist_n_energy_balance_meas.Write();
            TH2D hist_n_momentum_pminus_meas        = *rdf.Histo2D({("n_momentum_pminus_meas_"+ label).c_str(), ";P_{n} (GeV/c);P_{n}^{-} (GeV/c)", 100, 0.0, 1.0, 100, 0.4, 1.4},"n_momentum_meas","n_pminus_meas","accidweight");
            hist_n_momentum_pminus_meas.Write();

            TH1D hist_N2miss_mass_kin               = *rdf.Histo1D({("N2miss_mass_kin_"+ label).c_str(), ";m_{N2miss} (GeV/c^{2});Counts", 400, -4.0, 4.0},"N2miss_mass_kin","accidweight");
            hist_N2miss_mass_kin.Write();
            TH1D hist_N2miss_masssquared_kin        = *rdf.Histo1D({("N2miss_masssquared_kin_"+ label).c_str(), ";m_{N2miss}^{2} (GeV^{2}/c^{4});Counts", 400, 0.0, 4.0},"N2miss_masssquared_kin","accidweight");
            hist_N2miss_masssquared_kin.Write();
            TH1D hist_N2miss_momentum_kin           = *rdf.Histo1D({("N2miss_momentum_kin_"+ label).c_str(), ";P_{N2miss} (GeV/c);Counts", 100, 0.0, 1.0},"N2miss_momentum_kin","accidweight");
            hist_N2miss_momentum_kin.Write();
            TH1D hist_N2miss_pminus_kin             = *rdf.Histo1D({("N2miss_pminus_kin_"+ label).c_str(), ";P_{N2miss}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"N2miss_pminus_kin","accidweight");
            hist_N2miss_pminus_kin.Write();
            TH1D hist_N2miss_energy_balance_kin     = *rdf.Histo1D({("N2miss_energy_balance_kin_"+ label).c_str(), ";E_{N2miss} - m_{N2miss} (GeV);Counts", 400, -4.0, 4.0},"N2miss_energy_balance_kin","accidweight");
            hist_N2miss_energy_balance_kin.Write();
            TH2D hist_N2miss_momentum_pminus_kin    = *rdf.Histo2D({("N2miss_momentum_pminus_kin_"+ label).c_str(), ";P_{N2miss} (GeV/c);P_{N2miss}^{-} (GeV/c)", 100, 0.0, 1.0, 100, 0.4, 1.4},"N2miss_momentum_kin","N2miss_pminus_kin","accidweight");
            hist_N2miss_momentum_pminus_kin.Write();

            TH2D hist_n_N2miss_mass_kin             = *rdf.Histo2D({("n_N2miss_mass_kin_"+ label).c_str(), ";m_{n} (GeV/c^{2});m_{N2miss} (GeV/c^{2})", 400, -4.0, 4.0, 400, -4.0, 4.0},"n_mass_kin","N2miss_mass_kin","accidweight");
            hist_n_N2miss_mass_kin.Write();
            TH2D hist_n_N2miss_masssquared_kin      = *rdf.Histo2D({("n_N2miss_masssquared_kin_"+ label).c_str(), ";m_{n}^{2} (GeV^{2}/c^{4});m_{N2miss}^{2} (GeV^{2}/c^{4})", 400, 0.0, 4.0, 400, 0.0, 4.0},"n_masssquared_kin","N2miss_masssquared_kin","accidweight");
            hist_n_N2miss_masssquared_kin.Write();
            TH2D hist_n_N2miss_momentum_kin         = *rdf.Histo2D({("n_N2miss_momentum_kin_"+ label).c_str(), ";P_{n} (GeV/c);P_{N2miss} (GeV/c)", 100, 0.0, 1.0, 100, 0.0, 1.0},"n_momentum_kin","N2miss_momentum_kin","accidweight");
            hist_n_N2miss_momentum_kin.Write();
            TH2D hist_n_N2miss_pminus_kin           = *rdf.Histo2D({("n_N2miss_pminus_kin_"+ label).c_str(), ";P_{n}^{-} (GeV/c);P_{N2miss}^{-} (GeV/c)", 120, 0.3, 1.5, 120, 0.3, 1.5},"n_pminus_kin","N2miss_pminus_kin","accidweight");
            hist_n_N2miss_pminus_kin.Write();
            TH2D hist_n_N2miss_energy_balance_kin   = *rdf.Histo2D({("n_N2miss_energy_balance_kin_"+ label).c_str(), ";E_{n} - m_{n} (GeV);E_{N2miss} - m_{N2miss} (GeV)", 400, -4.0, 4.0, 400, -4.0, 4.0},"n_energy_balance_kin","N2miss_energy_balance_kin","accidweight");
            hist_n_N2miss_energy_balance_kin.Write();

            if (reaction_name.find("sim") != string::npos)
            {
                TH2D hist_pip_kinematics_truth          = *rdf.Histo2D({("pip_kinematics_truth_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_truth","pip_theta_truth");
                hist_pip_kinematics_truth.Write();

                TH2D hist_p_kinematics_truth            = *rdf.Histo2D({("p_kinematics_truth_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_truth","p_theta_truth");
                hist_p_kinematics_truth.Write();

                TH2D hist_pip_p_theta_truth             = *rdf.Histo2D({("pip_p_theta_truth_"+ label).c_str(), ";#theta_{#pi^{+}} (deg);#theta_{p} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"pip_theta_truth","p_theta_truth");
                hist_pip_p_theta_truth.Write();
                TH2D hist_pip_p_momentum_truth          = *rdf.Histo2D({("pip_p_momentum_truth_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);P_{p} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"pip_momentum_truth","p_momentum_truth");
                hist_pip_p_momentum_truth.Write();

                TH1D hist_sqrts_truth                   = *rdf.Histo1D({("sqrts_truth_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_truth");
                hist_sqrts_truth.Write();
                TH1D hist_minust_truth                  = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minust_truth");
                hist_minust_truth.Write();
                TH1D hist_coplanarity_truth             = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_truth");
                hist_coplanarity_truth.Write();
                TH1D hist_thetaCM_truth                 = *rdf.Histo1D({("thetaCM_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_truth");
                hist_thetaCM_truth.Write();
                TH2D hist_minust_thetaCM                = *rdf.Histo2D({("minust_thetaCM_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_truth","thetaCM_truth");
                hist_minust_thetaCM.Write();
                TH1D hist_rho_misid_mass_truth          = *rdf.Histo1D({("rho_misid_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misid_mass_truth");
                hist_rho_misid_mass_truth.Write();

                TH1D hist_n_mass_truth                  = *rdf.Histo1D({("n_mass_truth_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 100, 0.0, 4.0},"n_mass_truth");
                hist_n_mass_truth.Write();
                TH1D hist_n_momentum_truth              = *rdf.Histo1D({("n_momentum_truth_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 4.0},"n_momentum_truth");
                hist_n_momentum_truth.Write();
                TH1D hist_n_pminus_truth                = *rdf.Histo1D({("n_pminus_truth_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 100, 0.4, 1.4},"n_pminus_truth");
                hist_n_pminus_truth.Write();
                TH1D hist_n_energy_balance_truth        = *rdf.Histo1D({("n_energy_balance_truth_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"n_energy_balance_truth");
                hist_n_energy_balance_truth.Write();

                TH1D hist_N2miss_mass_truth             = *rdf.Histo1D({("N2miss_mass_truth_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 100, 0.0, 4.0},"N2miss_mass_truth");
                hist_N2miss_mass_truth.Write();
                TH1D hist_N2miss_momentum_truth         = *rdf.Histo1D({("N2miss_momentum_truth_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 4.0},"N2miss_momentum_truth");
                hist_N2miss_momentum_truth.Write();
                TH1D hist_N2miss_pminus_truth           = *rdf.Histo1D({("N2miss_pminus_truth_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 100, 0.4, 1.4},"N2miss_pminus_truth");
                hist_N2miss_pminus_truth.Write();
                TH1D hist_N2miss_energy_balance_truth   = *rdf.Histo1D({("N2miss_energy_balance_truth_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"N2miss_energy_balance_truth");
                hist_N2miss_energy_balance_truth.Write();

                TH2D hist_beam_energy_kin_truth         = *rdf.Histo2D({("beam_energy_kin_truth_"+ label).c_str(), ";E_{#gamma}^{kin} (deg);E_{#gamma}^{truth} (deg)", 55, 5.5, 11.0, 55, 5.5, 11.0},"beam_energy_kin","beam_energy_truth","accidweight");
                hist_beam_energy_kin_truth.Write();
                TH2D hist_thetaCM_kin_truth             = *rdf.Histo2D({("thetaCM_kin_truth_"+ label).c_str(), ";#theta_{CM}^{kin} (deg);#theta_{CM}^{truth} (deg)", 36, 0.0, 180.0, 36, 0.0, 180.0},"thetaCM_kin","thetaCM_truth","accidweight");
                hist_thetaCM_kin_truth.Write();
            }

            if (reaction_name.find("bggen") != string::npos)
            {
                TH1D hist_thrown_topology                           = *rdf.Histo1D({("thrown_topology_"+ label).c_str(), ";Thrown Topology;Counts", 10, -0.5, 9.5},"thrown_topology","accidweight");
                hist_thrown_topology.Write();
                TH1D hist_n_pminus_kin_pimprot                      = *rdf.Histo1D({("n_pminus_kin_pimprot_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_kin","topology_pimprot");
                hist_n_pminus_kin_pimprot.Write();
                TH1D hist_rho_misspion_mass_kin_pimprot             = *rdf.Histo1D({("rho_misspion_mass_kin_pimprot_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misspion_mass_kin","topology_pimprot");
                hist_rho_misspion_mass_kin_pimprot.Write();
                TH1D hist_n_pminus_kin_pippimprot                   = *rdf.Histo1D({("n_pminus_kin_pippimprot_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_kin","topology_pippimprot");
                hist_n_pminus_kin_pippimprot.Write();
                TH1D hist_rho_misspion_mass_kin_pippimprot          = *rdf.Histo1D({("rho_misspion_mass_kin_pippimprot_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misspion_mass_kin","topology_pippimprot");
                hist_rho_misspion_mass_kin_pippimprot.Write();
                TH2D hist_extrapion_kinematics_truth_pippimprot     = *rdf.Histo2D({("extrapion_kinematics_truth_pippimprot_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"extrapion_momentum_truth","extrapion_theta_truth","topology_pippimprot");
                hist_extrapion_kinematics_truth_pippimprot.Write();
                TH1D hist_rho_thrown_mass_truth_pippimprot          = *rdf.Histo1D({("rho_thrown_mass_truth_pippimprot_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_thrown_mass_truth","topology_pippimprot");
                hist_rho_thrown_mass_truth_pippimprot.Write();
                TH1D hist_n_pminus_kin_pi0pimprot                   = *rdf.Histo1D({("n_pminus_kin_pi0pimprot_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_kin","topology_pi0pimprot");
                hist_n_pminus_kin_pi0pimprot.Write();
                TH1D hist_rho_misspion_mass_kin_pi0pimprot          = *rdf.Histo1D({("rho_misspion_mass_kin_pi0pimprot_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misspion_mass_kin","topology_pi0pimprot");
                hist_rho_misspion_mass_kin_pi0pimprot.Write();
                TH2D hist_extrapion_kinematics_truth_pi0pimprot     = *rdf.Histo2D({("extrapion_kinematics_truth_pi0pimprot_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"extrapion_momentum_truth","extrapion_theta_truth","topology_pi0pimprot");
                hist_extrapion_kinematics_truth_pi0pimprot.Write();
                TH1D hist_rho_thrown_mass_truth_pi0pimprot          = *rdf.Histo1D({("rho_thrown_mass_truth_pi0pimprot_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_thrown_mass_truth","topology_pi0pimprot");
                hist_rho_thrown_mass_truth_pi0pimprot.Write();
                TH1D hist_n_pminus_kin_others                       = *rdf.Histo1D({("n_pminus_kin_others_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 120, 0.3, 1.5},"n_pminus_kin","topology_others");
                hist_n_pminus_kin_others.Write();
                TH1D hist_rho_misspion_mass_kin_others              = *rdf.Histo1D({("rho_misspion_mass_kin_others_"+ label).c_str(), ";m_{#pi^{+}#pi^{+}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_misspion_mass_kin","topology_others");
                hist_rho_misspion_mass_kin_others.Write();
            }
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}