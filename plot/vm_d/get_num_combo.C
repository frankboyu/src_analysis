#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <algorithm>
using namespace std;

int get_test()
{
    // TFile *input_treefile = new TFile("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_exc_data_2H.root", "UPDATE");
    TFile *input_treefile = new TFile("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_exc_sim_2H.root", "UPDATE");
    TTree *input_tree = (TTree*) input_treefile->Get("filteredtree_phi_d_recon");

    int nevents = 0;
    int nmultiple = 0;
    vector<Int_t> run_list;

    UInt_t run, run_before;
    ULong64_t event, event_before;
    Int_t beam_id, beam_id_before;
    Int_t kp_id, kp_id_before;
    Int_t km_id, km_id_before;
    Int_t d_id, d_id_before;
    double phi_mass_meas, phi_mass_meas_before;
    double kp_momentum_meas, kp_momentum_meas_before;
    double km_momentum_meas, km_momentum_meas_before;
    double d_momentum_meas, d_momentum_meas_before;
    double beam_energy_meas, beam_energy_meas_before;
    double kp_theta_meas, kp_theta_meas_before;
    double km_theta_meas, km_theta_meas_before;
    double d_theta_meas, d_theta_meas_before;
    double combo_accid_weight, combo_accid_weight_before;;

    input_tree->SetBranchAddress("run", &run);
    input_tree->SetBranchAddress("event", &event);
    input_tree->SetBranchAddress("phi_mass_meas", &phi_mass_meas);
    input_tree->SetBranchAddress("beam_id", &beam_id);
    input_tree->SetBranchAddress("kp_id", &kp_id);
    input_tree->SetBranchAddress("km_id", &km_id);
    input_tree->SetBranchAddress("d_id", &d_id);
    input_tree->SetBranchAddress("kp_momentum_meas", &kp_momentum_meas);
    input_tree->SetBranchAddress("km_momentum_meas", &km_momentum_meas);
    input_tree->SetBranchAddress("d_momentum_meas", &d_momentum_meas);
    input_tree->SetBranchAddress("beam_energy_meas", &beam_energy_meas);
    input_tree->SetBranchAddress("kp_theta_meas", &kp_theta_meas);
    input_tree->SetBranchAddress("km_theta_meas", &km_theta_meas);
    input_tree->SetBranchAddress("d_theta_meas", &d_theta_meas);
    input_tree->SetBranchAddress("combo_accid_weight", &combo_accid_weight);

    map<vector<int>, int> NumComboMap;
    map<vector<int>, int> HasBestComboMap;
    for (int i = 0; i < input_tree->GetEntries(); i++)
    {
        input_tree->GetEntry(i);

        vector<int> combo = {static_cast<int>(run), static_cast<int>(event), static_cast<int>(beam_id)};
        NumComboMap[combo]      += 1;
        HasBestComboMap[combo]  += combo_accid_weight;
    }

    int NumCombos, HasBestCombo;
    TBranch *NumCombos_branch       = input_tree->Branch("NumCombos", &NumCombos, "NumCombos/I");
    TBranch *HasBestCombo_branch    = input_tree->Branch("HasBestCombo", &HasBestCombo, "HasBestCombo/I");
    for (int i = 0; i < input_tree->GetEntries(); i++)
    {
        input_tree->GetEntry(i);
        vector<int> combo = {static_cast<int>(run), static_cast<int>(event), static_cast<int>(beam_id)};
        NumCombos = NumComboMap[combo];
        HasBestCombo = HasBestComboMap[combo];
        NumCombos_branch->Fill();
        HasBestCombo_branch->Fill();
    }
    input_tree->Write("", TObject::kOverwrite);

    return 0;
}