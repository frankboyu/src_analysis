using namespace std;
using namespace ROOT;
using namespace RooFit;

int plot_tree()
{
    // Open the file and get the tree
    TChain *chain = new TChain("flattree_phi_c_2H_data");
    chain->Add("/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_c_2H_data_ver01/*.root");
    TFile *output_file = new TFile("output/plots_phi_c_2H_data.root", "RECREATE");

    // Set branch addresses
    TLorentzVector *KPlusP4 = new TLorentzVector();
    TLorentzVector *KMinusP4 = new TLorentzVector();
    double WeightFactor;

    chain->SetBranchAddress("KPlusP4", &KPlusP4);
    chain->SetBranchAddress("KMinusP4", &KMinusP4);
    chain->SetBranchAddress("WeightFactor", &WeightFactor);

    // Create histogram
    TH2F *hist_MPhi_ThetaPhi = new TH2F("hist_MPhi_ThetaPhi", "hist_MPhi_ThetaPhi", 100, 0.9, 1.4, 180, 0.0, 30.0);

    // Loop over tree entries
    for (Long64_t i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);

        // Calculate variables
        double phi_mass = (*KPlusP4 + *KMinusP4).M();
        double phi_theta = (*KPlusP4 + *KMinusP4).Theta() * 180.0 / TMath::Pi();

        // Fill histogram
        hist_MPhi_ThetaPhi->Fill(phi_mass, phi_theta, WeightFactor);
    }

    // Save histograms to file
    hist_MPhi_ThetaPhi->Write();
    output_file->Close();

    return 0;
}
