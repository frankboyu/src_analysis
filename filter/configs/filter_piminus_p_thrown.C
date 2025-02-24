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

void filter_piminus_p_thrown(string reaction_name, string input_mode, string output_mode)
{
    // Read input files
    cout << "Reading input files...\n";
    string input_name, input_tree_name;
    if (reaction_name.find("tagged") != string::npos)
    {
        input_name      = Form("/work/halld2/home/boyu/src_analysis/selection/output/selectedtree_piminus_p_thrown_%s.root",reaction_name.c_str());
        input_tree_name = "selectedtree_piminus_p_thrown";
    }
    else
    {
        input_tree_name = "genT";
        if (reaction_name == "gen_2H_model")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/generator/*.root";
        else if (reaction_name == "gen_2H_flat")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/generator/*.root";
        else if (reaction_name == "gen_4He_model")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/generator/*.root";
        else if (reaction_name == "gen_4He_flat")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/generator/*.root";
        else if (reaction_name == "gen_12C_model")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/generator/*.root";
        else if (reaction_name == "gen_12C_flat")
            input_name = "/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/generator/*.root";
    }
    string hist_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_piminus_p_thrown_%s.root",reaction_name.c_str());
    string tree_name    = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_piminus_p_thrown_%s.root",reaction_name.c_str());

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
        .Define("pim_p4_truth", "pMeson")
        .Define("p_p4_truth", "pBaryon")
        ;
    }

    auto rdf_input = rdf_def
    .Define("target_p4",                    "TLorentzVector(0, 0, 0, mass_target)")
    .Define("N2_p4",                        "TLorentzVector(0, 0, 0, mass_2H)")

    .Define("beam_p4com_truth",             "boost_lorentz_vector(beam_p4_truth, -(pim_p4_truth + p_p4_truth).BoostVector())")
    .Define("beam_energy_truth",            "beam_p4_truth.E()")

    .Define("pim_p4com_truth",              "boost_lorentz_vector(pim_p4_truth, -(pim_p4_truth + p_p4_truth).BoostVector())")
    .Define("pim_energy_truth",             "pim_p4_truth.E()")
    .Define("pim_momentum_truth",           "pim_p4_truth.P()")
    .Define("pim_theta_truth",              "pim_p4_truth.Theta()*RadToDeg")

    .Define("p_p4pion_truth",               "TLorentzVector(p_p4_truth.Vect(), TMath::Sqrt(p_p4_truth.P()*p_p4_truth.P() + mass_piplus*mass_piplus))")
    .Define("p_energy_truth",               "p_p4_truth.E()")
    .Define("p_momentum_truth",             "p_p4_truth.P()")
    .Define("p_theta_truth",                "p_p4_truth.Theta()*RadToDeg")

    .Define("n_p4_truth",                   "pim_p4_truth + p_p4_truth - beam_p4_truth")
    .Define("n_energy_truth",               "n_p4_truth.E()")
    .Define("n_mass_truth",                 "n_p4_truth.M()")
    .Define("n_momentum_truth",             "n_p4_truth.P()")
    .Define("n_pminus_truth",               "n_p4_truth.Minus()")
    .Define("n_energy_balance_truth",       "n_p4_truth.E() - mass_neutron")

    .Define("miss_p4_truth",                "beam_p4_truth + target_p4 - pim_p4_truth - p_p4_truth")
    .Define("miss_energy_truth",            "miss_p4_truth.E()")
    .Define("miss_mass_truth",              "miss_p4_truth.M()")
    .Define("miss_momentum_truth",          "miss_p4_truth.P()")
    .Define("miss_pminus_truth",            "miss_p4_truth.Minus()")

    .Define("N2miss_p4_truth",              "beam_p4_truth + N2_p4 - pim_p4_truth - p_p4_truth")
    .Define("N2miss_energy_truth",          "N2miss_p4_truth.E()")
    .Define("N2miss_mass_truth",            "N2miss_p4_truth.M()")
    .Define("N2miss_momentum_truth",        "N2miss_p4_truth.P()")
    .Define("N2miss_pminus_truth",          "N2miss_p4_truth.Minus()")
    .Define("N2miss_energy_balance_truth",  "N2miss_p4_truth.E() - mass_proton")

    .Define("sqrts_truth",                  "(pim_p4_truth + p_p4_truth).Mag()")
    .Define("minust_truth",                 "-(pim_p4_truth - beam_p4_truth).Mag2()")
    .Define("minusu_truth",                 "-(p_p4_truth - beam_p4_truth).Mag2()")
    .Define("coplanarity_truth",            "abs(pim_p4_truth.Phi() - p_p4_truth.Phi())*RadToDeg")
    .Define("thetaCM_truth",                "beam_p4com_truth.Vect().Angle(pim_p4com_truth.Vect())*RadToDeg")
    .Define("rho_mass_truth",               "(pim_p4_truth + p_p4pion_truth).M()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    string miss_p_cut;
    if (reaction_name.find("2H") != string::npos)
        miss_p_cut = "0.20";
    else if (reaction_name.find("4He") != string::npos)
        miss_p_cut = "0.25";
    else if (reaction_name.find("12C") != string::npos)
        miss_p_cut = "0.30";

    auto rdf_NoCut          = rdf_input;
    auto rdf_MissPCut       = rdf_NoCut.Filter("(n_momentum_truth < " + miss_p_cut + ")");
    auto rdf_output         = rdf_MissPCut;

    RNode rdfs []       = {rdf_NoCut,   rdf_MissPCut};
    string labels []    = {"NoCut",     "MissPCut"};
    int N_filters = sizeof(labels) / sizeof(labels[0]);

    // Save tree
    if (output_mode == "tree" || output_mode == "both")
    {
        cout << "Saving to new tree...\n";
        rdf_output.Snapshot("filteredtree_piminus_p_thrown",tree_name.c_str());
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

            TH2D hist_pim_kinematics_truth = *rdf.Histo2D({("pim_kinematics_truth_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);#theta_{#pi^{-}} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"pim_momentum_truth","pim_theta_truth");
            hist_pim_kinematics_truth.Write();

            TH2D hist_p_kinematics_truth = *rdf.Histo2D({("p_kinematics_truth_"+ label).c_str(), ";P_{p} (GeV/c);#theta_{p} (deg)", 100, 0.0, 10.0, 180, 0.0, 180.0},"p_momentum_truth","p_theta_truth");
            hist_p_kinematics_truth.Write();

            TH2D hist_pim_p_theta_truth = *rdf.Histo2D({("pim_p_theta_truth_"+ label).c_str(), ";#theta_{#pi^{-}} (deg);#theta_{p} (deg)", 180, 0.0, 180.0, 180, 0.0, 180.0},"pim_theta_truth","p_theta_truth");
            hist_pim_p_theta_truth.Write();
            TH2D hist_pim_p_momentum_truth = *rdf.Histo2D({("pim_p_momentum_truth_"+ label).c_str(), ";P_{#pi^{-}} (GeV/c);P_{p} (GeV/c)", 100, 0.0, 10.0, 100, 0.0, 10.0},"pim_momentum_truth","p_momentum_truth");
            hist_pim_p_momentum_truth.Write();

            TH1D hist_sqrts_truth = *rdf.Histo1D({("sqrts_truth_"+ label).c_str(), ";#sqrt{s} (GeV);Counts", 100, 0.0, 10.0},"sqrts_truth");
            hist_sqrts_truth.Write();
            TH1D hist_minust_truth = *rdf.Histo1D({("minust_truth_"+ label).c_str(), ";-t (GeV^{2}/c^{2});Counts", 200, 0.0, 20.0},"minust_truth");
            hist_minust_truth.Write();
            TH1D hist_coplanarity_truth = *rdf.Histo1D({("coplanarity_truth_"+ label).c_str(), ";Coplanarity (deg);Counts", 360, 0.0, 360.0},"coplanarity_truth");
            hist_coplanarity_truth.Write();
            TH1D hist_thetaCM_truth = *rdf.Histo1D({("thetaCM_truth_"+ label).c_str(), ";#theta_{CM} (deg);Counts", 180, 0.0, 180.0},"thetaCM_truth");
            hist_thetaCM_truth.Write();
            TH2D hist_minust_thetaCM = *rdf.Histo2D({("minust_thetaCM_"+ label).c_str(), ";-t (GeV^{2}/c^{2});#theta_{CM} (deg)", 200, 0.0, 20.0, 180, 0.0, 180.0},"minust_truth","thetaCM_truth");
            hist_minust_thetaCM.Write();
            TH1D hist_rho_mass_truth = *rdf.Histo1D({("rho_mass_truth_"+ label).c_str(), ";m_{#pi^{+}#pi^{-}} (GeV/c^{2});Counts", 400, 0.0, 4.0},"rho_mass_truth");
            hist_rho_mass_truth.Write();

            TH1D hist_n_mass_truth = *rdf.Histo1D({("n_mass_truth_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 100, 0.0, 4.0},"n_mass_truth");
            hist_n_mass_truth.Write();
            TH1D hist_n_momentum_truth = *rdf.Histo1D({("n_momentum_truth_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 4.0},"n_momentum_truth");
            hist_n_momentum_truth.Write();
            TH1D hist_n_pminus_truth = *rdf.Histo1D({("n_pminus_truth_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 100, 0.4, 1.4},"n_pminus_truth");
            hist_n_pminus_truth.Write();
            TH1D hist_n_energy_balance_truth = *rdf.Histo1D({("n_energy_balance_truth_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"n_energy_balance_truth");
            hist_n_energy_balance_truth.Write();

            TH1D hist_N2miss_mass_truth = *rdf.Histo1D({("N2miss_mass_truth_"+ label).c_str(), ";m_{n} (GeV/c^{2});Counts", 100, 0.0, 4.0},"N2miss_mass_truth");
            hist_N2miss_mass_truth.Write();
            TH1D hist_N2miss_momentum_truth = *rdf.Histo1D({("N2miss_momentum_truth_"+ label).c_str(), ";P_{n} (GeV/c);Counts", 100, 0.0, 4.0},"N2miss_momentum_truth");
            hist_N2miss_momentum_truth.Write();
            TH1D hist_N2miss_pminus_truth = *rdf.Histo1D({("N2miss_pminus_truth_"+ label).c_str(), ";P_{n}^{-} (GeV/c);Counts", 100, 0.4, 1.4},"N2miss_pminus_truth");
            hist_N2miss_pminus_truth.Write();
            TH1D hist_N2miss_energy_balance_truth = *rdf.Histo1D({("N2miss_energy_balance_truth_"+ label).c_str(), ";E_{n} - m_{n} (GeV);Counts", 400, -4.0, 4.0},"N2miss_energy_balance_truth");
            hist_N2miss_energy_balance_truth.Write();
        }
        hist_file->Close();
    }
    cout << "Done!\n";
}