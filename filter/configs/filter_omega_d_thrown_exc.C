#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include </work/halld2/home/boyu/src_analysis/filter/configs/const.h>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_target = 0.0;

TLorentzVector boost_lorentz_vector(TLorentzVector p4, TVector3 boost_vector)
{
    TLorentzVector p4_boosted(p4);
    p4_boosted.Boost(boost_vector);
    return p4_boosted;
}

void filter_rho_d_thrown_exc(string reaction_name, string output_mode)
{
    // Read input files
    cout << "Reading input files...\n";
    string input_name, input_tree_name;
    if (reaction_name.find("tagged") != string::npos)
    {
        input_name      = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_rho_d_thrown_%s.root",reaction_name.c_str());
        input_tree_name = "selectedtree_rho_d_thrown";
    }
    else
    {
        input_tree_name = "genT";
        if (reaction_name == "gen_2H")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/rho_d_2H_ver01/root/generator/*.root";
        else if (reaction_name == "gen_4He")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/rho_d_4He_ver01/root/generator/*.root";
        else if (reaction_name == "gen_12C")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/rho_d_12C_ver01/root/generator/*.root";
    }
    string hist_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_rho_d_thrown_%s.root",reaction_name.c_str());
    string tree_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_rho_d_thrown_%s.root",reaction_name.c_str());

    TChain chain(input_tree_name.c_str());
    chain.Add(input_name.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    if      (reaction_name.find("2H") != string::npos)
        mass_target = mass_2H;
    else if (reaction_name.find("4He") != string::npos)
        mass_target = mass_4He;
    else if (reaction_name.find("12C") != string::npos)
        mass_target = mass_12C;

    RDataFrame rdf_raw(chain);
    auto rdf_def = RNode(rdf_raw);
    if (reaction_name.find("gen") != string::npos)
    {
        rdf_def = rdf_def
        .Define("beam_p4_truth", "pBeam")
        .Define("pip_p4_truth",  "TLorentzVector(0, 0, 0, 0)")
        .Define("pip_p4kaon_truth", "TLorentzVector(0, 0, 0, 0)")
        .Define("pim_p4_truth",  "TLorentzVector(0, 0, 0, 0)")
        .Define("pim_p4kaon_truth", "TLorentzVector(0, 0, 0, 0)")
        .Define("rho_p4_truth", "pMeson")
        .Define("d_p4_truth", "pBaryon")
        ;
    }
    else if (reaction_name.find("tagged") != string::npos)
    {
        rdf_def = rdf_def
        .Define("pip_p4kaon_truth",         "TLorentzVector(pip_p4_truth.Vect(), TMath::Sqrt(pip_p4_truth.P()*pip_p4_truth.P() + mass_kplus*mass_kplus))")
        .Define("pim_p4kaon_truth",         "TLorentzVector(pim_p4_truth.Vect(), TMath::Sqrt(pim_p4_truth.P()*pim_p4_truth.P() + mass_kminus*mass_kminus))")
        .Define("rho_p4_truth",             "pip_p4_truth + pim_p4_truth")
        ;
    }

    auto rdf_input = rdf_def
    .Define("target_p4",                    "TLorentzVector(0, 0, 0, mass_target)")
    .Define("N2_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(rho_p4_truth + d_p4_truth).BoostVector())")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")

    .Define("pip_energy_truth",             "pip_p4_truth.E()")
    .Define("pip_momentum_truth",           "pip_p4_truth.P()")
    .Define("pip_theta_truth",              "pip_p4_truth.Theta()*RadToDeg")

    .Define("pim_energy_truth",             "pim_p4_truth.E()")
    .Define("pim_momentum_truth",           "pim_p4_truth.P()")
    .Define("pim_theta_truth",              "pim_p4_truth.Theta()*RadToDeg")

    .Define("d_energy_truth",               "d_p4_truth.E()")
    .Define("d_momentum_truth",             "d_p4_truth.P()")
    .Define("d_theta_truth",                "d_p4_truth.Theta()*RadToDeg")

    .Define("rho_p4com_truth",              "boost_lorentz_vector(rho_p4_truth, -(rho_p4_truth + d_p4_truth).BoostVector())")
    .Define("rho_energy_truth",             "rho_p4_truth.E()")
    .Define("rho_momentum_truth",           "rho_p4_truth.P()")
    .Define("rho_mass_truth",               "rho_p4_truth.M()")
    .Define("rho_theta_truth",              "rho_p4_truth.Theta()*RadToDeg")

    .Define("struck_p4_truth",              "rho_p4_truth + d_p4_truth - beam_p4_truth")
    .Define("struck_energy_truth",          "struck_p4_truth.E()")
    .Define("struck_mass_truth",            "struck_p4_truth.M()")
    .Define("struck_masssquared_truth",     "struck_p4_truth.M2()")
    .Define("struck_momentum_truth",        "struck_p4_truth.P()")
    .Define("struck_pminus_truth",          "struck_p4_truth.Minus()")
    .Define("struck_theta_truth",           "struck_p4_truth.Theta()*RadToDeg")
    .Define("struck_energy_balance_truth",  "struck_p4_truth.E() - mass_2H")

    .Define("miss_p4_truth",                "beam_p4_truth + target_p4 - rho_p4_truth - d_p4_truth")
    .Define("miss_energy_truth",            "miss_p4_truth.E()")
    .Define("miss_mass_truth",              "miss_p4_truth.M()")
    .Define("miss_masssquared_truth",       "miss_p4_truth.M2()")
    .Define("miss_momentum_truth",          "miss_p4_truth.P()")
    .Define("miss_pminus_truth",            "miss_p4_truth.Minus()")
    .Define("miss_theta_truth",             "miss_p4_truth.Theta()*RadToDeg")
    .Define("miss_energy_balance_truth",    "miss_p4_truth.E() - mass_2H")

    .Define("sqrts_truth",                  "(rho_p4_truth + d_p4_truth).Mag()")
    .Define("minust_truth",                 "-(rho_p4_truth - beam_p4_truth).Mag2()")
    .Define("minusu_truth",                 "-(d_p4_truth - beam_p4_truth).Mag2()")
    .Define("coplanarity_truth",            "abs(rho_p4_truth.Phi() - d_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(rho_p4com_truth.Vect())*RadToDeg")
    .Define("y_rho_truth",                  "minust_truth/(2*mass_2H*(beam_p4_truth.E()-rho_p4_truth.E()))")
    .Define("phi_mass_truth",               "(pip_p4kaon_truth + pim_p4kaon_truth).M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";

    auto rdf_NoCut          = rdf_input;
    auto rdf_output         = rdf_NoCut;

    RNode rdfs []       = {rdf_NoCut};
    string labels []    = {"NoCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_rho_d_thrown",tree_name.c_str());
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

            TH2D hist_rho_kinematics_truth = *rdf.Histo2D({("rho_kinematics_truth_"+ label).c_str(), ";P_{#rho} (GeV/c);#theta_{#rho} (deg)", 60, 5.0, 11.0, 20, 0.0, 20.0},"rho_momentum_truth","rho_theta_truth");
            hist_rho_kinematics_truth.Write();

            TH2D hist_d_kinematics_truth = *rdf.Histo2D({("d_kinematics_truth_"+ label).c_str(), ";P_{d} (GeV/c);#theta_{d} (deg)", 200, 0.0, 2.0, 90, 0.0, 90.0},"d_momentum_truth","d_theta_truth");
            hist_d_kinematics_truth.Write();

            TH2D hist_rho_d_theta_truth = *rdf.Histo2D({("rho_d_theta_truth_"+ label).c_str(), ";#theta_{#rho} (deg);#theta_{d} (deg)", 20, 0.0, 20.0, 90, 0.0, 90.0},"rho_theta_truth","d_theta_truth");
            hist_rho_d_theta_truth.Write();
            TH2D hist_rho_d_momentum_truth = *rdf.Histo2D({("rho_d_momentum_truth_"+ label).c_str(), ";P_{#rho} (GeV/c);P_{d} (GeV/c)", 60, 5.0, 11.0, 200, 0.0, 2.0},"rho_momentum_truth","d_momentum_truth");
            hist_rho_d_momentum_truth.Write();

            TH1D hist_sqrts_truth = *rdf.Histo1D({("sqrts_truth_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_truth");
            hist_sqrts_truth.Write();
            TH1D hist_minust_truth = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 300, 0.0, 3.0},"minust_truth");
            hist_minust_truth.Write();
            TH1D hist_minusu_truth = *rdf.Histo1D({("minusu_truth_"+ label).c_str(), ";-u (GeV^{2}/c^{2});Counts", 250, 15.0, 40.0},"minusu_truth");
            hist_minusu_truth.Write();
            TH1D hist_coplanarity_truth = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_truth");
            hist_coplanarity_truth.Write();
            TH1D hist_thetaCM_truth = *rdf.Histo1D({("thetaCM_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_truth");
            hist_thetaCM_truth.Write();
            TH2D hist_minust_thetaCM = *rdf.Histo2D({("minust_thetaCM_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_truth","thetaCM_truth");
            hist_minust_thetaCM.Write();
            TH1D hist_phi_mass_truth = *rdf.Histo1D({("phi_mass_truth_"+ label).c_str(), ";m_{#K^{+}#K^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"phi_mass_truth");
            hist_phi_mass_truth.Write();

            TH1D hist_struck_mass_truth = *rdf.Histo1D({("struck_mass_truth_"+ label).c_str(), ";m_{struck} (GeV/c^{2});Counts", 100, 0.0, 4.0},"struck_mass_truth");
            hist_struck_mass_truth.Write();
            TH1D hist_struck_momentum_truth = *rdf.Histo1D({("struck_momentum_truth_"+ label).c_str(), ";P_{struck} (GeV/c);Counts", 100, 0.0, 4.0},"struck_momentum_truth");
            hist_struck_momentum_truth.Write();
            TH1D hist_struck_pminus_truth = *rdf.Histo1D({("struck_pminus_truth_"+ label).c_str(), ";P_{struck}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"struck_pminus_truth");
            hist_struck_pminus_truth.Write();
            TH1D hist_struck_energy_balance_truth = *rdf.Histo1D({("struck_energy_balance_truth_"+ label).c_str(), ";E_{struck} - m_{struck} (GeV);Counts", 400, -4.0, 4.0},"struck_energy_balance_truth");
            hist_struck_energy_balance_truth.Write();

            TH1D hist_miss_mass_truth = *rdf.Histo1D({("miss_mass_truth_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts", 100, 0.0, 4.0},"miss_mass_truth");
            hist_miss_mass_truth.Write();
            TH1D hist_miss_momentum_truth = *rdf.Histo1D({("miss_momentum_truth_"+ label).c_str(), ";P_{miss} (GeV/c);Counts", 100, 0.0, 4.0},"miss_momentum_truth");
            hist_miss_momentum_truth.Write();
            TH1D hist_miss_pminus_truth = *rdf.Histo1D({("miss_pminus_truth_"+ label).c_str(), ";P_{miss}^{-} (GeV/c);Counts", 200, 1.0, 3.0},"miss_pminus_truth");
            hist_miss_pminus_truth.Write();
            TH1D hist_miss_energy_balance_truth = *rdf.Histo1D({("miss_energy_balance_truth_"+ label).c_str(), ";E_{miss} - m_{miss} (GeV);Counts", 400, -4.0, 4.0},"miss_energy_balance_truth");
            hist_miss_energy_balance_truth.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}