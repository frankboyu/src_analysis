#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <cstring>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_deuteron    = 1.875612;
double mass_proton      = 0.938272;
double mass_neutron     = 0.939565;
double mass_pion        = 0.139570;
double RadToDeg         = 180.0 / 3.14159265;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_piminus_p_recon(string Tag, string InputMode, string OutputMode)
{
    string InputFile  = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_piminus_p_recon_%s.root",Tag.c_str());
    string HistFile   = Form("output/filteredhist_piminus_p_recon_%s.root",Tag.c_str());
    string TreeFile   = Form("output/filteredtree_piminus_p_recon_%s.root",Tag.c_str());

    // Read input files
    cout << "Reading input files...\n";
    TChain chain("selectedtree_piminus_p_recon");
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_input = RNode(rdf_raw);
    if (InputMode == "one" && Reaction.find("2H") != string::npos)
        rdf_input = rdf_input.Filter("run == 90213");
    else if (InputMode == "one" && Reaction.find("4He") != string::npos)
        rdf_input = rdf_input.Filter("run == 90061");
    else if (InputMode == "one" && Reaction.find("12C") != string::npos)
        rdf_input = rdf_input.Filter("run == 90291");

    auto rdf_def = rdf_input
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")

    .Define("beam_p4com_meas","boost_lorentz_vector(beam_p4_meas, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("beam_p4com_kin","boost_lorentz_vector(beam_p4_kin, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("beam_p4com_truth","boost_lorentz_vector(beam_p4_truth, -(pim_p4_truth + p_p4_truth).BoostVector())")
    .Define("beam_E_meas","beam_p4_meas.E()")
    .Define("beam_E_kin","beam_p4_kin.E()")
    .Define("beam_E_truth","beam_p4_truth.E()")

    .Define("miss_p4_meas","pim_p4_meas + p_p4_meas - beam_p4_meas")
    .Define("miss_p4_kin","pim_p4_kin + p_p4_kin - beam_p4_kin")
    .Define("miss_p4_truth","pim_p4_truth + p_p4_truth - beam_p4_truth")
    .Define("miss_mass_meas","miss_p4_meas.M()")
    .Define("miss_mass_kin","miss_p4_kin.M()")
    .Define("miss_mass_truth","miss_p4_truth.M()")
    .Define("miss_p_meas","miss_p4_meas.P()")
    .Define("miss_p_kin","miss_p4_kin.P()")
    .Define("miss_p_truth","miss_p4_truth.P()")
    .Define("miss_pminus_meas","miss_p4_meas.Minus()")
    .Define("miss_pminus_kin","miss_p4_kin.Minus()")
    .Define("miss_pminus_truth","miss_p4_truth.Minus()")
    .Define("miss_energy_balance_meas","miss_p4_meas.E() - mass_neutron")
    .Define("miss_energy_balance_kin","miss_p4_kin.E() - mass_neutron")
    .Define("miss_energy_balance_truth","miss_p4_truth.E() - mass_neutron")

    .Define("pim_p4com_meas","boost_lorentz_vector(pim_p4_meas, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("pim_p4com_kin","boost_lorentz_vector(pim_p4_kin, -(pim_p4_meas + p_p4_meas).BoostVector())")
    .Define("pim_p4com_truth","boost_lorentz_vector(pim_p4_truth, -(pim_p4_truth + p_p4_truth).BoostVector())")
    .Define("pim_p_meas","pim_p4_meas.P()")
    .Define("pim_p_kin","pim_p4_kin.P()")
    .Define("pim_p_truth","pim_p4_truth.P()")
    .Define("pim_theta_meas","pim_p4_meas.Theta()*RadToDeg")
    .Define("pim_theta_kin","pim_p4_kin.Theta()*RadToDeg")
    .Define("pim_theta_truth","pim_p4_truth.Theta()*RadToDeg")
    .Define("pim_in_fdc", "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc == 0.0)")
    .Define("pim_in_cdc", "accidweight*(pim_dedx_cdc > 0.0 && pim_dedx_fdc == 0.0)")
    .Define("pim_in_fdc_cdc", "accidweight*(pim_dedx_fdc > 0.0 && pim_dedx_cdc > 0.0)")
    .Define("pim_in_neither", "accidweight*(pim_dedx_fdc == 0.0 && pim_dedx_cdc == 0.0)")

    .Define("p_as_pip_p4_meas","TLorentzVector(p_p4_meas.Vect(), TMath::Sqrt(p_p4_meas.P()*p_p4_meas.P() + mass_pion*mass_pion))")
    .Define("p_as_pip_p4_kin","TLorentzVector(p_p4_kin.Vect(), TMath::Sqrt(p_p4_kin.P()*p_p4_kin.P() + mass_pion*mass_pion))")
    .Define("p_as_pip_p4_truth","TLorentzVector(p_p4_truth.Vect(), TMath::Sqrt(p_p4_truth.P()*p_p4_truth.P() + mass_pion*mass_pion))")
    .Define("p_p_meas","p_p4_meas.P()")
    .Define("p_p_kin","p_p4_kin.P()")
    .Define("p_p_truth","p_p4_truth.P()")
    .Define("p_theta_meas","p_p4_meas.Theta()*RadToDeg")
    .Define("p_theta_kin","p_p4_kin.Theta()*RadToDeg")
    .Define("p_theta_truth","p_p4_truth.Theta()*RadToDeg")
    .Define("p_in_fdc", "accidweight*(p_dedx_fdc > 0.0 && p_dedx_cdc == 0.0)")
    .Define("p_in_cdc", "accidweight*(p_dedx_cdc > 0.0 && p_dedx_fdc == 0.0)")
    .Define("p_in_fdc_cdc", "accidweight*(p_dedx_fdc > 0.0 && p_dedx_cdc > 0.0)")
    .Define("p_in_neither", "accidweight*(p_dedx_fdc == 0.0 && p_dedx_cdc == 0.0)")

    .Define("sqrt_s_meas", "(pim_p4_meas + p_p4_meas).Mag()")
    .Define("sqrt_s_kin", "(pim_p4_meas + p_p4_meas).Mag()")
    .Define("sqrt_s_truth", "(pim_p4_truth + p_p4_truth).Mag()")
    .Define("minus_t_meas", "-(pim_p4_meas - beam_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(pim_p4_kin - beam_p4_kin).Mag2()")
    .Define("minus_t_truth", "-(pim_p4_truth - beam_p4_truth).Mag2()")
    .Define("minus_u_meas", "-(p_p4_meas - beam_p4_meas).Mag2()")
    .Define("minus_u_kin", "-(p_p4_kin - beam_p4_kin).Mag2()")
    .Define("minus_u_truth", "-(p_p4_truth - beam_p4_truth).Mag2()")
    .Define("coplanarity_meas","abs(pim_p4_meas.Phi() - p_p4_meas.Phi())*RadToDeg")
    .Define("coplanarity_kin","abs(pim_p4_kin.Phi() - p_p4_kin.Phi())*RadToDeg")
    .Define("coplanarity_truth","abs(pim_p4_truth.Phi() - p_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_meas", "beam_p4com_meas.Vect().Angle(pim_p4com_meas.Vect())*RadToDeg")
    .Define("thetaCM_kin", "beam_p4com_kin.Vect().Angle(pim_p4com_kin.Vect())*RadToDeg")
    .Define("thetaCM_truth", "beam_p4com_truth.Vect().Angle(pim_p4com_truth.Vect())*RadToDeg")

    .Define("rho_mass_meas","(pim_p4_meas + p_as_pip_p4_meas).M()")
    .Define("rho_mass_kin","(pim_p4_kin + p_as_pip_p4_kin).M()")
    .Define("rho_mass_truth","(pim_p4_truth + p_as_pip_p4_truth).M()")
    .Define("coherent_2pi_missing_mass_meas","(beam_p4_meas + TLorentzVector(0.0, 0.0, 0.0, mass_deuteron) - pim_p4_meas - p_as_pip_p4_meas).M()")
    .Define("coherent_2pi_missing_mass_kin","(beam_p4_kin + TLorentzVector(0.0, 0.0, 0.0, mass_deuteron) - pim_p4_kin - p_as_pip_p4_kin).M()")
    .Define("coherent_2pi_missing_mass_truth","(beam_p4_truth + TLorentzVector(0.0, 0.0, 0.0, mass_deuteron) - pim_p4_truth - p_as_pip_p4_truth).M()")

    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    string miss_p_cut; // Default value
    if (Tag.find("inc") != string::npos)
        miss_p_cut = "0.5";
    else if (Tag.find("missprot") != string::npos)
        miss_p_cut = "0.2";
    else if (Tag.find("misstri") != string::npos)
        miss_p_cut = "0.25";
    else if (Tag.find("missb11") != string::npos)
        miss_p_cut = "0.3";

    auto rdf_no_filter                  = rdf_def.Filter("(minus_t_kin > 0.5) && (minus_u_kin > 0.5)");
    auto rdf_cl_filtered                = rdf_no_filter.Filter("kin_cl > 0.01");
    auto rdf_pidfom_filtered            = rdf_cl_filtered.Filter("(pim_pidfom > 0.01) && (p_pidfom > 0.01)");
    auto rdf_miss_p_filtered            = rdf_pidfom_filtered.Filter("(miss_p_kin < " + miss_p_cut + ")");
    auto rdf_miss_pminus_filtered       = rdf_miss_p_filtered.Filter("(miss_pminus_kin > 0.5) && (miss_pminus_kin < 1.3)");
    auto rdf_kinematics_filtered        = rdf_miss_pminus_filtered.Filter("(minus_t_kin > 0.5) && (minus_u_kin > 0.5)");
    auto rdf_output                     = rdf_kinematics_filtered;

    // Save tree
    if (OutputMode == "tree" || OutputMode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_piminus_p_recon",TreeFile);
    }

    // Save histograms
    if (OutputMode == "hist" || OutputMode == "both")
    {
        cout << "Plotting histograms...\n";
        TFile * histFile = new TFile(HistFile.c_str(),"RECREATE");
        histFile->cd();
        vector<TH1*> hist_list;

        int N_filters = 6;
        RNode rdfs [] = {rdf_no_filter, rdf_cl_filtered, rdf_pidfom_filtered, rdf_miss_p_filtered, rdf_miss_pminus_filtered, rdf_kinematics_filtered};
        string labels [] = {"NoCut", "CLCut", "PIDFOMCut", "MissPCut", "MissPminusCut", "KinematicsCut"};


        for (int i = 0; i < N_filters; i++)
        {
            auto rdf = rdfs[i];
            string label = labels[i];
            cout << "Processing " << label << "...\n";

            TDirectory * dir = histFile->mkdir(label.c_str());
            dir->cd();

            // Data columns
            TH1D hist_accidweight = *rdf.Histo1D({("accidweight_"+ label).c_str(), ";accidweight;Counts", 20, -0.5, 1.5},"accidweight");
            hist_accidweight.Write();
            TH1D hist_kin_cl = *rdf.Histo1D({("kin_cl_"+ label).c_str(), ";kin_cl;Counts", 100, 0.0, 1.0},"kin_cl","accidweight");
            hist_kin_cl.Write();
            TH1D hist_run = *rdf.Histo1D({("run_"+ label).c_str(), ";run;Counts", 700, 90000, 90700},"run","accidweight");
            hist_run.Write();

            // PID
            TH1D hist_pim_pidfom = *rdf.Histo1D({("pim_pidfom_"+ label).c_str(), ";pim_pidfom;Counts", 100, 0.0, 1.0},"pim_pidfom","accidweight");
            hist_pim_pidfom.Write();
            TH2D hist_pim_dEdx_cdc = *rdf.Histo2D({("pim_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_cdc","accidweight");
            hist_pim_dEdx_cdc.Write();
            TH2D hist_pim_dEdx_fdc = *rdf.Histo2D({("pim_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_fdc","accidweight");
            hist_pim_dEdx_fdc.Write();
            TH2D hist_pim_dEdx_fdc_cdc = *rdf.Histo2D({("pim_dEdx_fdc_cdc_"+ label).c_str(), ";FDC dE/dx (keV/cm);CDC dE/dx (keV/cm)", 200, 0.0, 2e-5, 200, 0.0, 2e-5},"pim_dedx_fdc","pim_dedx_cdc","accidweight");
            hist_pim_dEdx_fdc_cdc.Write();
            TH2D hist_pim_dEdx_tof = *rdf.Histo2D({("pim_dEdx_tof_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_tof","accidweight");
            hist_pim_dEdx_tof.Write();
            TH2D hist_pim_dEdx_st = *rdf.Histo2D({("pim_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"pim_p_kin","pim_dedx_st","accidweight");
            hist_pim_dEdx_st.Write();

            TH1D hist_p_pidfom = *rdf.Histo1D({("p_pidfom_"+ label).c_str(), ";p_pidfom;Counts", 100, 0.0, 1.0},"p_pidfom","accidweight");
            hist_p_pidfom.Write();
            TH2D hist_p_dEdx_cdc = *rdf.Histo2D({("p_dEdx_cdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_cdc","accidweight");
            hist_p_dEdx_cdc.Write();
            TH2D hist_p_dEdx_fdc = *rdf.Histo2D({("p_dEdx_fdc_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_fdc","accidweight");
            hist_p_dEdx_fdc.Write();
            TH2D hist_p_dEdx_fdc_cdc = *rdf.Histo2D({("p_dEdx_fdc_cdc_"+ label).c_str(), ";FDC dE/dx (keV/cm);CDC dE/dx (keV/cm)", 200, 0.0, 2e-5, 200, 0.0, 2e-5},"p_dedx_fdc","p_dedx_cdc","accidweight");
            hist_p_dEdx_fdc_cdc.Write();
            TH2D hist_p_dEdx_tof = *rdf.Histo2D({("p_dEdx_tof_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_tof","accidweight");
            hist_p_dEdx_tof.Write();
            TH2D hist_p_dEdx_st = *rdf.Histo2D({("p_dEdx_st_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"p_p_kin","p_dedx_st","accidweight");
            hist_p_dEdx_st.Write();

            // Background
            TH1D hist_rho_mass = *rdf.Histo1D({("rho_mass_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_mass_kin","accidweight");
            hist_rho_mass.Write();
            TH1D hist_coherent_2pi = *rdf.Histo1D({("coherent_2pi_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}X} (GeV/c^{2});Counts", 400, 0.0, 4.0},"coherent_2pi_missing_mass_kin","accidweight");
            hist_coherent_2pi.Write();

            // Initial neutron
            TH1D hist_miss_pminus = *rdf.Histo1D({("miss_pminus_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 100, 0.4, 1.4},"miss_pminus_kin","accidweight");
            hist_miss_pminus.Write();
            TH1D hist_miss_mass = *rdf.Histo1D({("miss_mass_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 100, 0.0, 4.0},"miss_mass_kin","accidweight");
            hist_miss_mass.Write();
            TH1D hist_miss_p = *rdf.Histo1D({("miss_p_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 100, 0.0, 4.0},"miss_p_kin","accidweight");
            hist_miss_p.Write();
            TH2D hist_miss_pminus_miss_p = *rdf.Histo2D({("miss_pminus_miss_p_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);P_{miss} (GeV/c)", 100, 0.4, 1.4, 50, 0.0, 1.0},"miss_pminus_kin","miss_p_kin","accidweight");
            hist_miss_pminus_miss_p.Write();
            TH2D hist_miss_p_energy_balance = *rdf.Histo2D({("miss_p_energy_balance_"+ label).c_str(), ";P_{miss} (GeV/c);E_{miss} - m_{n} (GeV)", 100, 0.0, 4.0, 400, -4.0, 4.0},"miss_p_kin","miss_energy_balance_kin","accidweight");
            hist_miss_p_energy_balance.Write();

            // Kinematics
            TH1D hist_sqrt_s = *rdf.Histo1D({("sqrt_s_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrt_s_kin","accidweight");
            hist_sqrt_s.Write();
            TH1D hist_minus_t = *rdf.Histo1D({("minus_t_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minus_t_kin","accidweight");
            hist_minus_t.Write();
            TH1D hist_thetaCM = *rdf.Histo1D({("thetaCM_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_kin","accidweight");
            hist_thetaCM.Write();
            TH2D hist_minus_t_thetaCM = *rdf.Histo2D({("minus_t_thetaCM_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minus_t_kin","thetaCM_kin","accidweight");
            hist_minus_t_thetaCM.Write();
            TH2D hist_theta_pim_p = *rdf.Histo2D({("theta_pim_p_"+ label).c_str(), ";#theta_{#pi^{-}} (deg);#theta_{p} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"pim_theta_kin","p_theta_kin","accidweight");
            hist_theta_pim_p.Write();
            TH1D hist_coplanarity = *rdf.Histo1D({("coplanarity_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_kin","accidweight");
            hist_coplanarity.Write();

            TH2D hist_kinematics_pim = *rdf.Histo2D({("kinematics_pim_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_p_kin","pim_theta_kin","accidweight");
            hist_kinematics_pim.Write();
            TH2D hist_kinematics_pim_fdc = *rdf.Histo2D({("kinematics_pim_fdc_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_p_kin","pim_theta_kin","pim_in_fdc");
            hist_kinematics_pim_fdc.Write();
            TH2D hist_kinematics_pim_fdc_cdc = *rdf.Histo2D({("kinematics_pim_fdc_cdc_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_p_kin","pim_theta_kin","pim_in_fdc_cdc");
            hist_kinematics_pim_fdc_cdc.Write();
            TH2D hist_kinematics_pim_cdc = *rdf.Histo2D({("kinematics_pim_cdc_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_p_kin","pim_theta_kin","pim_in_cdc");
            hist_kinematics_pim_cdc.Write();
            TH2D hist_kinematics_pim_neither = *rdf.Histo2D({("kinematics_pim_neither_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_p_kin","pim_theta_kin","pim_in_neither");
            hist_kinematics_pim_neither.Write();

            TH2D hist_kinematics_p = *rdf.Histo2D({("kinematics_p_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_p_kin","p_theta_kin","accidweight");
            hist_kinematics_p.Write();
            TH2D hist_kinematics_p_fdc = *rdf.Histo2D({("kinematics_p_fdc_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_p_kin","p_theta_kin","p_in_fdc");
            hist_kinematics_p_fdc.Write();
            TH2D hist_kinematics_p_fdc_cdc = *rdf.Histo2D({("kinematics_p_fdc_cdc_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_p_kin","p_theta_kin","p_in_fdc_cdc");
            hist_kinematics_p_fdc_cdc.Write();
            TH2D hist_kinematics_p_cdc = *rdf.Histo2D({("kinematics_p_cdc_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_p_kin","p_theta_kin","p_in_cdc");
            hist_kinematics_p_cdc.Write();
            TH2D hist_kinematics_p_neither = *rdf.Histo2D({("kinematics_p_neither_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_p_kin","p_theta_kin","p_in_neither");
            hist_kinematics_p_neither.Write();

            // Thrown information
            TH1D hist_miss_p_truth = *rdf.Histo1D({("miss_p_truth_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 100, 0.0, 4.0},"miss_p_truth","accidweight");
            hist_miss_p_truth.Write();
            TH2D hist_thetaCM_wrt_truth = *rdf.Histo2D({("thetaCM_wrt_truth_"+ label).c_str(), ";#theta_{CM} (deg);#theta_{CM}^{truth} (deg)", 36, 0.0, 180.0, 36, 0.0, 180.0},"thetaCM_kin","thetaCM_truth","accidweight");
            hist_thetaCM_wrt_truth.Write();

        }
        histFile->Close();
    }

    cout << "Done!\n";
}