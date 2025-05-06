#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

double mass_target = 0.0;

void filter_omega_d_recon(string reaction_name, string output_mode)
{
    string input_name   = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_omega_d_recon_%s.root",reaction_name.c_str());
    string hist_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_omega_d_recon_%s.root",reaction_name.c_str());
    string tree_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_omega_d_recon_%s.root",reaction_name.c_str());

    // Determine reaction specific parameters
    if (reaction_name.find("2H") != string::npos)
        mass_target = mass_2H;
    else if (reaction_name.find("4He") != string::npos)
        mass_target = mass_4He;
    else if (reaction_name.find("12C") != string::npos)
        mass_target = mass_12C;

    // Read input files
    cout << "Reading input files...\n";
    TChain chain("selectedtree_omega_d_recon");
    chain.Add(input_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_input = rdf_raw
    .Define("kinfit_fom",                   "TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4",                    "TLorentzVector(0, 0, 0, mass_target)")
    .Define("N2_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")
    .Define("total_p4_meas",                "pip_p4_meas + pim_p4_meas + d_p4_meas")
    .Define("total_p4_kin",                 "pip_p4_kin + pim_p4_kin + d_p4_kin")
    .Define("total_p4_truth",               "pip_p4_truth + pim_p4_truth + d_p4_truth")

    .Define("beam_p4com_meas",              "boost_lorentz_vector(beam_p4_meas, -(pip_p4_meas + pim_p4_meas + d_p4_meas).BoostVector())")
    .Define("beam_p4com_kin",               "boost_lorentz_vector(beam_p4_kin, -(pip_p4_kin + pim_p4_kin + d_p4_kin).BoostVector())")
    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(pip_p4_truth + pim_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_meas",             "beam_p4_meas.E()")
    .Define("beam_energy_kin",              "beam_p4_kin.E()")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")
    .Define("beam_DeltaT_meas",             "beam_x4_meas.T() - rftime")
    .Define("beam_DeltaT_kin",              "beam_x4_kin.T() - rftime")
    .Define("beam_DeltaT_truth",            "beam_x4_truth.T() - rftime")

    .Define("pip_p4kaon_meas",               "TLorentzVector(pip_p4_meas.Vect(), TMath::Sqrt(pip_p4_meas.P()*pip_p4_meas.P() + mass_kplus*mass_kplus))")
    .Define("pip_p4kaon_kin",                "TLorentzVector(pip_p4_kin.Vect(), TMath::Sqrt(pip_p4_kin.P()*pip_p4_kin.P() + mass_kplus*mass_kplus))")
    .Define("pip_p4kaon_truth",              "TLorentzVector(pip_p4_truth.Vect(), TMath::Sqrt(pip_p4_truth.P()*pip_p4_truth.P() + mass_kplus*mass_kplus))")
    .Define("pip_energy_meas",               "pip_p4_meas.E()")
    .Define("pip_energy_kin",                "pip_p4_kin.E()")
    .Define("pip_energy_truth",              "pip_p4_truth.E()")
    .Define("pip_momentum_meas",             "pip_p4_meas.P()")
    .Define("pip_momentum_kin",              "pip_p4_kin.P()")
    .Define("pip_momentum_truth",            "pip_p4_truth.P()")
    .Define("pip_theta_meas",                "pip_p4_meas.Theta()*RadToDeg")
    .Define("pip_theta_kin",                 "pip_p4_kin.Theta()*RadToDeg")
    .Define("pip_theta_truth",               "pip_p4_truth.Theta()*RadToDeg")
    .Define("pip_in_fdc",                    "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_in_cdc",                    "accidweight*(pip_dedx_cdc > 0.0 && pip_dedx_fdc == 0.0)")
    .Define("pip_in_fdc_cdc",                "accidweight*(pip_dedx_fdc > 0.0 && pip_dedx_cdc > 0.0)")
    .Define("pip_in_neither",                "accidweight*(pip_dedx_fdc == 0.0 && pip_dedx_cdc == 0.0)")
    .Define("pip_DeltaT_meas",               "rftime + (pip_x4_meas.Z()-65.0)/29.9792458 - pip_x4_meas.T()")
    .Define("pip_DeltaT_kin",                "rftime + (pip_x4_kin.Z()-65.0)/29.9792458 - pip_x4_kin.T()")
    .Define("pip_DeltaT_truth",              "rftime + (pip_x4_truth.Z()-65.0)/29.9792458 - pip_x4_truth.T()")
    .Define("pip_dedx_fdc_keV_per_cm",       "pip_dedx_fdc*1e6")
    .Define("pip_dedx_cdc_keV_per_cm",       "pip_dedx_cdc*1e6")
    .Define("pip_dedx_st_keV_per_cm",        "pip_dedx_st*1e6")
    .Define("pip_dedx_tof_keV_per_cm",       "pip_dedx_tof*1e6")

    .Define("pim_p4kaon_meas",               "TLorentzVector(pim_p4_meas.Vect(), TMath::Sqrt(pim_p4_meas.P()*pim_p4_meas.P() + mass_kminus*mass_kminus))")
    .Define("pim_p4kaon_kin",                "TLorentzVector(pim_p4_kin.Vect(), TMath::Sqrt(pim_p4_kin.P()*pim_p4_kin.P() + mass_kminus*mass_kminus))")
    .Define("pim_p4kaon_truth",              "TLorentzVector(pim_p4_truth.Vect(), TMath::Sqrt(pim_p4_truth.P()*pim_p4_truth.P() + mass_kminus*mass_kminus))")
    .Define("pim_energy_meas",               "pim_p4_meas.E()")
    .Define("pim_energy_kin",                "pim_p4_kin.E()")
    .Define("pim_energy_truth",              "pim_p4_truth.E()")
    .Define("pim_momentum_meas",             "pim_p4_meas.P()")
    .Define("pim_momentum_kin",              "pim_p4_kin.P()")
    .Define("pim_momentum_truth",            "pim_p4_truth.P()")
    .Define("pim_theta_meas",                "pim_p4_meas.Theta()*RadToDeg")
    .Define("pim_theta_kin",                 "pim_p4_kin.Theta()*RadToDeg")
    .Define("pim_theta_truth",               "pim_p4_truth.Theta()*RadToDeg")
    .Define("pim_in_fdc",                    "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc == 0.0)")
    .Define("pim_in_cdc",                    "accidweight*(pim_dedx_cdc > 0.0 && pim_dedx_fdc == 0.0)")
    .Define("pim_in_fdc_cdc",                "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc > 0.0)")
    .Define("pim_in_neither",                "accidweight*(pim_dedx_fdc == 0.0 && pim_dedx_cdc == 0.0)")
    .Define("pim_DeltaT_meas",               "rftime + (pim_x4_meas.Z()-65.0)/29.9792458 - pim_x4_meas.T()")
    .Define("pim_DeltaT_kin",                "rftime + (pim_x4_kin.Z()-65.0)/29.9792458 - pim_x4_kin.T()")
    .Define("pim_DeltaT_truth",              "rftime + (pim_x4_truth.Z()-65.0)/29.9792458 - pim_x4_truth.T()")
    .Define("pim_dedx_fdc_keV_per_cm",       "pim_dedx_fdc*1e6")
    .Define("pim_dedx_cdc_keV_per_cm",       "pim_dedx_cdc*1e6")
    .Define("pim_dedx_st_keV_per_cm",        "pim_dedx_st*1e6")
    .Define("pim_dedx_tof_keV_per_cm",       "pim_dedx_tof*1e6")

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
    .Define("d_DeltaT_truth",               "rftime + (d_x4_meas.Z()-65.0)/29.9792458 - d_x4_truth.T()")
    .Define("d_dedx_fdc_keV_per_cm",        "d_dedx_fdc*1e6")
    .Define("d_dedx_cdc_keV_per_cm",        "d_dedx_cdc*1e6")
    .Define("d_dedx_st_keV_per_cm",         "d_dedx_st*1e6")
    .Define("d_dedx_tof_keV_per_cm",        "d_dedx_tof*1e6")

    .Define("rho_p4_meas",                  "pip_p4_meas + pim_p4_meas")
    .Define("rho_p4_kin",                   "pip_p4_kin + pim_p4_kin")
    .Define("rho_p4_truth",                 "pip_p4_truth + pim_p4_truth")
    .Define("rho_p4com_meas",               "boost_lorentz_vector(rho_p4_meas, -(rho_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_p4com_kin",                "boost_lorentz_vector(rho_p4_kin, -(rho_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_p4com_truth",              "boost_lorentz_vector(rho_p4_truth, -(rho_p4_meas + d_p4_meas).BoostVector())")
    .Define("rho_energy_meas",              "rho_p4_meas.E()")
    .Define("rho_energy_kin",               "rho_p4_kin.E()")
    .Define("rho_energy_truth",             "rho_p4_truth.E()")
    .Define("rho_momentum_meas",            "rho_p4_meas.P()")
    .Define("rho_momentum_kin",             "rho_p4_kin.P()")
    .Define("rho_momentum_truth",           "rho_p4_truth.P()")
    .Define("rho_mass_meas",                "rho_p4_meas.M()")
    .Define("rho_mass_kin",                 "rho_p4_kin.M()")
    .Define("rho_mass_truth",               "rho_p4_truth.M()")
    .Define("rho_proxymass_meas",           "TMath::Sqrt((pip_p4_meas.Minus() + pim_p4_meas.Minus())*(2*beam_p4_meas.E() + mass_4He - d_p4_meas.Plus() - (mass_2H*mass_2H + total_p4_meas.Perp()*total_p4_meas.Perp())/(mass_4He - total_p4_meas.Minus())) - (pip_p4_meas.Px() + pim_p4_meas.Px())*(pip_p4_meas.Px() + pim_p4_meas.Px()) - (pip_p4_meas.Py() + pim_p4_meas.Py())*(pip_p4_meas.Py() + pim_p4_meas.Py()))")
    .Define("rho_proxymass_kin",            "TMath::Sqrt((pip_p4_kin.Minus() + pim_p4_kin.Minus())*(2*beam_p4_kin.E() + mass_4He - d_p4_kin.Plus() - (mass_2H*mass_2H + total_p4_kin.Perp()*total_p4_kin.Perp())/(mass_4He - total_p4_kin.Minus())) - (pip_p4_kin.Px() + pim_p4_kin.Px())*(pip_p4_kin.Px() + pim_p4_kin.Px()) - (pip_p4_kin.Py() + pim_p4_kin.Py())*(pip_p4_kin.Py() + pim_p4_kin.Py()))")
    .Define("rho_proxymass_truth",          "TMath::Sqrt((pip_p4_truth.Minus() + pim_p4_truth.Minus())*(2*beam_p4_truth.E() + mass_4He - d_p4_truth.Plus() - (mass_2H*mass_2H + total_p4_truth.Perp()*total_p4_truth.Perp())/(mass_4He - total_p4_truth.Minus())) - (pip_p4_truth.Px() + pim_p4_truth.Px())*(pip_p4_truth.Px() + pim_p4_truth.Px()) - (pip_p4_truth.Py() + pim_p4_truth.Py())*(pip_p4_truth.Py() + pim_p4_truth.Py()))")
    .Define("rho_theta_meas",               "rho_p4_meas.Theta()*RadToDeg")
    .Define("rho_theta_kin",                "rho_p4_kin.Theta()*RadToDeg")
    .Define("rho_theta_truth",              "rho_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_meas",               "rho_p4_meas + d_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin",                "rho_p4_kin + d_p4_kin - beam_p4_kin")
    .Define("struck_p4_truth",              "rho_p4_truth + d_p4_truth - beam_p4_truth")
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

    .Define("miss_p4_meas",                 "beam_p4_meas + target_p4 - rho_p4_meas - d_p4_meas")
    .Define("miss_p4_kin",                  "beam_p4_kin + target_p4 - rho_p4_kin - d_p4_kin")
    .Define("miss_p4_truth",                "beam_p4_truth + target_p4 - rho_p4_truth - d_p4_truth")
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

    .Define("sqrts_meas",                   "(rho_p4_meas + d_p4_meas).Mag()")
    .Define("sqrts_kin",                    "(rho_p4_kin + d_p4_kin).Mag()")
    .Define("sqrts_truth",                  "(rho_p4_truth + d_p4_truth).Mag()")
    .Define("minust_meas",                  "-(beam_p4_meas - rho_p4_meas).Mag2()")
    .Define("minust_kin",                   "-(beam_p4_kin - rho_p4_kin).Mag2()")
    .Define("minust_truth",                 "-(beam_p4_truth - rho_p4_truth).Mag2()")
    .Define("minusu_meas",                  "-(beam_p4_meas - d_p4_meas).Mag2()")
    .Define("minusu_kin",                   "-(beam_p4_kin - d_p4_kin).Mag2()")
    .Define("minusu_truth",                 "-(beam_p4_truth - d_p4_truth).Mag2()")
    .Define("coplanarity_meas",             "abs(rho_p4_meas.Phi() - d_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin",              "abs(rho_p4_kin.Phi() - d_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth",            "abs(rho_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_meas",                 "beam_p4com_meas.Vect().Angle(rho_p4com_meas.Vect())*RadToDeg")
    .Define("thetaCM_kin",                  "beam_p4com_kin.Vect().Angle(rho_p4com_kin.Vect())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(rho_p4com_truth.Vect())*RadToDeg")
    .Define("y_rho_meas",                   "minust_meas/(2*mass_2H*(beam_p4_meas.E()-rho_p4_meas.E()))")
    .Define("y_rho_kin",                    "minust_kin/(2*mass_2H*(beam_p4_kin.E()-rho_p4_kin.E()))")
    .Define("y_rho_truth",                  "minust_truth/(2*mass_2H*(beam_p4_truth.E()-rho_p4_truth.E()))")
    .Define("phi_mass_meas",                "(pip_p4kaon_meas + pim_p4kaon_meas).M()")
    .Define("phi_mass_kin",                 "(pip_p4kaon_kin + pim_p4kaon_kin).M()")
    .Define("phi_mass_truth",               "(pip_p4kaon_truth + pim_p4kaon_truth).M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_NoCut          = rdf_input;
    auto rdf_dEdxCut        = rdf_NoCut.Filter("(d_dedx_cdc_keV_per_cm > (TMath::Exp(-29.68353898*d_momentum_meas+13.50623694)+17.88279645*d_momentum_meas*d_momentum_meas-42.15473796*d_momentum_meas+28.83200736)) && (d_dedx_cdc_keV_per_cm < (TMath::Exp(-26.69276323*d_momentum_meas+15.92466317)+17.1164272*d_momentum_meas*d_momentum_meas-48.7542903*d_momentum_meas+40.25692313))");
    auto rdf_KinFitFOMCut   = rdf_dEdxCut.Filter("kinfit_fom > 0.01");
    auto rdf_PIDFOMCut      = rdf_KinFitFOMCut.Filter("(pip_pidfom > 0.01) && (pim_pidfom > 0.01)");
    auto rdf_MissPCut       = rdf_PIDFOMCut.Filter("(abs(struck_energy_balance_kin) < 1.0)");
    auto rdf_output         = rdf_MissPCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_dEdxCut,    rdf_KinFitFOMCut,   rdf_PIDFOMCut,  rdf_MissPCut};
    string labels []    = {"NoCut",     "dEdxCut",      "KinFitFOMCut",     "PIDFOMCut",    "MissPCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_omega_d_recon",tree_name);
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
            TH1D hist_beam_energy_kin               = *rdf.Histo1D({("beam_energy_"+ label).c_str(), ";E_{beam} (GeV);Counts", 60, 5.0, 11.0},"beam_energy_kin");
            hist_beam_energy_kin.Write();

            TH1D hist_pip_pidfom                    = *rdf.Histo1D({("pip_pidfom_"+ label).c_str(), ";pip_pidfom;Counts", 100, 0.0, 1.0},"pip_pidfom","accidweight");
            hist_pip_pidfom.Write();
            TH2D hist_pip_dEdx_cdc_kin              = *rdf.Histo2D({("pip_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_cdc_keV_per_cm","accidweight");
            hist_pip_dEdx_cdc_kin.Write();
            TH2D hist_pip_dEdx_fdc_kin              = *rdf.Histo2D({("pip_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_fdc_keV_per_cm","accidweight");
            hist_pip_dEdx_fdc_kin.Write();
            TH2D hist_pip_dEdx_tof_kin              = *rdf.Histo2D({("pip_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_tof_keV_per_cm","accidweight");
            hist_pip_dEdx_tof_kin.Write();
            TH2D hist_pip_dEdx_st_kin               = *rdf.Histo2D({("pip_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pip_momentum_kin","pip_dedx_st_keV_per_cm","accidweight");
            hist_pip_dEdx_st_kin.Write();
            TH1D hist_pip_DeltaT_kin                = *rdf.Histo1D({("pip_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{+}} (ns);Counts", 400, -20.0, 20.0},"pip_DeltaT_kin","accidweight");
            hist_pip_DeltaT_kin.Write();
            TH2D hist_pip_kinematics_kin            = *rdf.Histo2D({("pip_kinematics_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","accidweight");
            hist_pip_kinematics_kin.Write();
            TH2D hist_pip_kinematics_fdc_kin        = *rdf.Histo2D({("pip_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc");
            hist_pip_kinematics_fdc_kin.Write();
            TH2D hist_pip_kinematics_fdc_cdc_kin    = *rdf.Histo2D({("pip_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_fdc_cdc");
            hist_pip_kinematics_fdc_cdc_kin.Write();
            TH2D hist_pip_kinematics_cdc_kin        = *rdf.Histo2D({("pip_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_cdc");
            hist_pip_kinematics_cdc_kin.Write();
            TH2D hist_pip_kinematics_neither_kin    = *rdf.Histo2D({("pip_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{+}} (GeV/c);#theta_{#pi^{+}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pip_momentum_kin","pip_theta_kin","pip_in_neither");
            hist_pip_kinematics_neither_kin.Write();

            TH1D hist_pim_pidfom                    = *rdf.Histo1D({("pim_pidfom_"+ label).c_str(), ";pim_pidfom;Counts", 100, 0.0, 1.0},"pim_pidfom","accidweight");
            hist_pim_pidfom.Write();
            TH2D hist_pim_dEdx_cdc_kin              = *rdf.Histo2D({("pim_dEdx_cdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pim_momentum_kin","pim_dedx_cdc_keV_per_cm","accidweight");
            hist_pim_dEdx_cdc_kin.Write();
            TH2D hist_pim_dEdx_fdc_kin              = *rdf.Histo2D({("pim_dEdx_fdc_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pim_momentum_kin","pim_dedx_fdc_keV_per_cm","accidweight");
            hist_pim_dEdx_fdc_kin.Write();
            TH2D hist_pim_dEdx_tof_kin              = *rdf.Histo2D({("pim_dEdx_tof_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pim_momentum_kin","pim_dedx_tof_keV_per_cm","accidweight");
            hist_pim_dEdx_tof_kin.Write();
            TH2D hist_pim_dEdx_st_kin               = *rdf.Histo2D({("pim_dEdx_st_kin_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 500, 0.0, 10.0, 400, 0.0, 40},"pim_momentum_kin","pim_dedx_st_keV_per_cm","accidweight");
            hist_pim_dEdx_st_kin.Write();
            TH1D hist_pim_DeltaT_kin                = *rdf.Histo1D({("pim_DeltaT_kin_"+ label).c_str(), ";#Delta t_{#pi^{-}} (ns);Counts", 400, -20.0, 20.0},"pim_DeltaT_kin","accidweight");
            hist_pim_DeltaT_kin.Write();
            TH2D hist_pim_kinematics_kin            = *rdf.Histo2D({("pim_kinematics_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","accidweight");
            hist_pim_kinematics_kin.Write();
            TH2D hist_pim_kinematics_fdc_kin        = *rdf.Histo2D({("pim_kinematics_fdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_fdc");
            hist_pim_kinematics_fdc_kin.Write();
            TH2D hist_pim_kinematics_fdc_cdc_kin    = *rdf.Histo2D({("pim_kinematics_fdc_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_fdc_cdc");
            hist_pim_kinematics_fdc_cdc_kin.Write();
            TH2D hist_pim_kinematics_cdc_kin        = *rdf.Histo2D({("pim_kinematics_cdc_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_cdc");
            hist_pim_kinematics_cdc_kin.Write();
            TH2D hist_pim_kinematics_neither_kin    = *rdf.Histo2D({("pim_kinematics_neither_kin_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_kin","pim_theta_kin","pim_in_neither");
            hist_pim_kinematics_neither_kin.Write();

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
            TH2D hist_d_kinematics_kin              = *rdf.Histo2D({("d_kinematics_kin_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"d_momentum_kin","d_theta_kin","accidweight");
            hist_d_kinematics_kin.Write();

            TH1D hist_rho_mass_kin                  = *rdf.Histo1D({("rho_mass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 600, 0.0, 1.5},"rho_mass_kin","accidweight");
            hist_rho_mass_kin.Write();
            TH1D hist_rho_mass_meas                 = *rdf.Histo1D({("rho_mass_meas_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 600, 0.0, 1.5},"rho_mass_meas","accidweight");
            hist_rho_mass_meas.Write();
            TH1D hist_rho_proxymass_kin             = *rdf.Histo1D({("rho_proxymass_kin_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 600, 0.0, 1.5},"rho_proxymass_kin","accidweight");
            hist_rho_proxymass_kin.Write();
            TH1D hist_rho_proxymass_meas            = *rdf.Histo1D({("rho_proxymass_meas_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 600, 0.0, 1.5},"rho_proxymass_meas","accidweight");
            hist_rho_proxymass_meas.Write();
            TH2D hist_rho_kinematics_kin            = *rdf.Histo2D({("rho_kinematics_kin_"+ label).c_str(), ";p (GeV/c);#theta (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"rho_momentum_kin","rho_theta_kin","accidweight");
            hist_rho_kinematics_kin.Write();

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
            TH1D hist_phi_mass_kin                  = *rdf.Histo1D({("phi_mass_kin_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"phi_mass_kin","accidweight");
            hist_phi_mass_kin.Write();
            TH1D hist_yrho_kin                      = *rdf.Histo1D({("yrho_kin_"+ label).c_str(), ";y_{#rho};Counts", 200, 0.0, 2.0},"y_rho_kin","accidweight");
            hist_yrho_kin.Write();

            TH1D hist_struck_mass_kin               = *rdf.Histo1D({("struck_mass_kin_"+ label).c_str(), ";m_{struck} (GeV/c^{2});Counts", 400, -4.0, 4.0},"struck_mass_kin","accidweight");
            hist_struck_mass_kin.Write();
            TH1D hist_struck_mass_meas              = *rdf.Histo1D({("struck_mass_meas_"+ label).c_str(), ";m_{struck} (GeV/c^{2});Counts", 400, -4.0, 4.0},"struck_mass_meas","accidweight");
            hist_struck_mass_meas.Write();
            TH1D hist_struck_masssquared_kin        = *rdf.Histo1D({("struck_masssquared_kin_"+ label).c_str(), ";m_{struck}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"struck_masssquared_kin","accidweight");
            hist_struck_masssquared_kin.Write();
            TH1D hist_struck_masssquared_meas       = *rdf.Histo1D({("struck_masssquared_meas_"+ label).c_str(), ";m_{struck}^{2} (GeV^{2}/c^{4});Counts", 100, 0.0, 10.0},"struck_masssquared_meas","accidweight");
            hist_struck_masssquared_meas.Write();
            TH1D hist_struck_momentum_kin           = *rdf.Histo1D({("struck_momentum_kin_"+ label).c_str(), ";P_{struck} (GeV/c);Counts", 200, 0.0, 2.0},"struck_momentum_kin","accidweight");
            hist_struck_momentum_kin.Write();
            TH1D hist_struck_momentum_meas          = *rdf.Histo1D({("struck_momentum_meas_"+ label).c_str(), ";P_{struck} (GeV/c);Counts", 200, 0.0, 2.0},"struck_momentum_meas","accidweight");
            hist_struck_momentum_meas.Write();
            TH1D hist_struck_pminus_kin             = *rdf.Histo1D({("struck_pminus_kin_"+ label).c_str(), ";P_{struck}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"struck_pminus_kin","accidweight");
            hist_struck_pminus_kin.Write();
            TH1D hist_struck_pminus_meas            = *rdf.Histo1D({("struck_pminus_meas_"+ label).c_str(), ";P_{struck}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"struck_pminus_meas","accidweight");
            hist_struck_pminus_meas.Write();
            TH1D hist_struck_energy_balance_kin     = *rdf.Histo1D({("struck_energy_balance_kin_"+ label).c_str(), ";E_{struck} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"struck_energy_balance_kin","accidweight");
            hist_struck_energy_balance_kin.Write();
            TH1D hist_struck_energy_balance_meas    = *rdf.Histo1D({("struck_energy_balance_meas_"+ label).c_str(), ";E_{struck} - m_{2H} (GeV);Counts", 400, -4.0, 4.0},"struck_energy_balance_meas","accidweight");
            hist_struck_energy_balance_meas.Write();
            TH2D hist_struck_momentum_pminus_kin    = *rdf.Histo2D({("struck_momentum_pminus_kin_"+ label).c_str(), ";P_{struck} (GeV/c);P_{struck}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"struck_momentum_kin","struck_pminus_kin","accidweight");
            hist_struck_momentum_pminus_kin.Write();
            TH2D hist_struck_momentum_pminus_meas   = *rdf.Histo2D({("struck_momentum_pminus_meas_"+ label).c_str(), ";P_{struck} (GeV/c);P_{struck}^{-} (GeV/c)", 200, 0.0, 2.0, 200, 1.0, 3.0},"struck_momentum_meas","struck_pminus_meas","accidweight");
            hist_struck_momentum_pminus_meas.Write();
            TH2D hist_struck_momentum_energy_kin    = *rdf.Histo2D({("struck_momentum_energy_kin_"+ label).c_str(), ";P_{struck} (GeV/c);E_{struck} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"struck_momentum_kin","struck_energy_kin","accidweight");
            hist_struck_momentum_energy_kin.Write();
            TH2D hist_struck_momentum_energy_meas   = *rdf.Histo2D({("struck_momentum_energy_meas_"+ label).c_str(), ";P_{struck} (GeV/c);E_{struck} (GeV)", 200, 0.0, 2.0, 250, 0.5, 3.0},"struck_momentum_meas","struck_energy_meas","accidweight");
            hist_struck_momentum_energy_meas.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}