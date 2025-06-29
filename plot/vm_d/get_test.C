#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
using namespace std;

typedef enum
{
    KPlus = 0,
    KMinus = 1,
    Deuteron = 2,
    Unknown = 3
} Particle_t;

int get_test()
{
    TFile *input_treefile = new TFile("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_phi_d_recon_data_2H_exc.root", "read");
    TTree *input_tree = (TTree*) input_treefile->Get("filteredtree_phi_d_recon");

    int nevents = 0;
    int nmultiple = 0;

    UInt_t run, run_before;
    ULong64_t event, event_before;
    Int_t beam_id, beam_id_before;
    Int_t kp_id, kp_id_before;
    Int_t km_id, km_id_before;
    Int_t d_id, d_id_before;
    double phi_mass_meas, phi_mass_meas_before;

    input_tree->SetBranchAddress("run", &run);
    input_tree->SetBranchAddress("event", &event);
    input_tree->SetBranchAddress("phi_mass_meas", &phi_mass_meas);
    input_tree->SetBranchAddress("beam_id", &beam_id);
    input_tree->SetBranchAddress("kp_id", &kp_id);
    input_tree->SetBranchAddress("km_id", &km_id);
    input_tree->SetBranchAddress("d_id", &d_id);

    set<map<Particle_t, set<Int_t> > > locUsedSoFar_Mass;

    // for (int i = 0; i < 100; i++)
    for (int i = 0; i < input_tree->GetEntries(); i++)
    {
        run_before = run;
        event_before = event;
        phi_mass_meas_before = phi_mass_meas;
        beam_id_before = beam_id;
        kp_id_before = kp_id;
        km_id_before = km_id;
        d_id_before = d_id;

        input_tree->GetEntry(i);

        map<Particle_t, set<Int_t> > locUsedThisCombo_Mass;
		locUsedThisCombo_Mass[KPlus].insert(kp_id);
		locUsedThisCombo_Mass[KMinus].insert(km_id);
		locUsedThisCombo_Mass[Deuteron].insert(d_id);
        locUsedThisCombo_Mass[Unknown].insert(beam_id);

        // cout << "Run: " << run << ", Event: " << event << ", Phi Mass Meas: " << phi_mass_meas << endl;
        if (run_before != run || event_before != event)
        {
            nevents++;
            // if (locUsedSoFar_Mass.size() > 1)
            // {
                // nmultiple++;
                // cout << "Multiple entries found for Run: " << run_before << ", Event: " << event_before << endl;
                // cout << locUsedSoFar_Mass.size() << " entries found." << endl;
                // for (const auto& entry : locUsedSoFar_Mass)
                // {
                //     cout << "Combination: ";
                //     for (const auto& particle : entry)
                //     {
                //         cout << particle.first << ": ";
                //         for (const auto& id : particle.second)
                //         {
                //             cout << id << " ";
                //         }
                //         cout << "; ";
                //     }
                //     cout << endl;
                // }
            // }
            locUsedSoFar_Mass.clear();
            locUsedSoFar_Mass.insert(locUsedThisCombo_Mass);
        }
        else
        {
            if(locUsedSoFar_Mass.find(locUsedThisCombo_Mass) == locUsedSoFar_Mass.end())
                locUsedSoFar_Mass.insert(locUsedThisCombo_Mass);
            else
                nmultiple++;
        }
    }

    cout << "Total number of events: " << nevents << endl;
    cout << "Number of multiple entries: " << nmultiple << endl;

    return 0;

}