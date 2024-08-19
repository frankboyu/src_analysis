using namespace std;
using namespace ROOT;
using namespace RooFit;

double mass_proton = 0.938272;
double mass_deuteron = 1.875612;
double mass_pion = 0.13957018;

int filter_phi_p_2H_data()
{
    // Open the file and get the tree
    TChain *chain = new TChain("flattree_phi_p_2H_data");
    chain->Add("/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_p_2H_data_test/flattree_phi_p_2H_data_1_initial.root");
    TFile *output_file = new TFile("output/plots_phi_p_2H_data.root", "RECREATE");

    // Set branch addresses
    TLorentzVector *BeamP4      = new TLorentzVector();
    TLorentzVector *KPlusP4     = new TLorentzVector();
    TLorentzVector *KMinusP4    = new TLorentzVector();
    TLorentzVector *ProtonP4    = new TLorentzVector();
    TLorentzVector *MissingP4   = new TLorentzVector();
    double WeightFactor;

    chain->SetBranchAddress("BeamP4", &BeamP4);
    chain->SetBranchAddress("MissingP4", &MissingP4);
    chain->SetBranchAddress("KPlusP4", &KPlusP4);
    chain->SetBranchAddress("KMinusP4", &KMinusP4);
    chain->SetBranchAddress("ProtonP4", &ProtonP4);
    chain->SetBranchAddress("WeightFactor", &WeightFactor);

    // Create histogram
    TH1F *hist_Mkk                  = new TH1F("hist_Mkk", "hist_Mkk", 400, 0.9, 1.3);
    TH1F *hist_Mpkk                 = new TH1F("hist_Mpkk", "hist_Mpkk", 400, 1.9, 2.3);
    TH2F *hist_Mkk_Mpkk             = new TH2F("hist_Mkk_Mpkk", "hist_Mkk_Mpkk", 400, 0.9, 1.3, 400, 1.9, 2.3);

    // Loop over tree entries
    for (Long64_t i = 0; i < chain->GetEntries(); i++)
    {
        chain->GetEntry(i);

        // Calculate variables
        // TLorentzVector *TargetP4 = new TLorentzVector(0.0, 0.0, 0.0, mass_deuteron);
        TLorentzVector *PhiP4 = new TLorentzVector(*KPlusP4 + *KMinusP4);
        TLorentzVector *KPlusP4AsPion = new TLorentzVector(KPlusP4->X(), KPlusP4->Y(), KPlusP4->Z(), sqrt(KPlusP4->P()*KPlusP4->P() + mass_pion*mass_pion));
        TLorentzVector *KMinusP4AsPion = new TLorentzVector(KMinusP4->X(), KMinusP4->Y(), KMinusP4->Z(), sqrt(KMinusP4->P()*KMinusP4->P() + mass_pion*mass_pion));

        double phi_mass     = (*KPlusP4 + *KMinusP4).M();
        double rho_mass     = (*KPlusP4AsPion + *KMinusP4AsPion).M();
        double phi_theta    = (*KPlusP4 + *KMinusP4).Theta() * 180.0 / TMath::Pi();
        double sqrt_s       = (*PhiP4 + *ProtonP4).Mag();
        double minus_t      = -(*BeamP4 - *PhiP4).Mag2();
        double minus_u      = -(*BeamP4 - *ProtonP4).Mag2();

        // Fill histograms: before filter

        // Filter events
        if ((*ProtonP4 + *KPlusP4).M() > 1.5) continue;
        if ((*ProtonP4 + *KMinusP4).M() > 1.5) continue;
        if ((*KPlusP4 + *KMinusP4).M() > 1.05) continue;
        cout << (*KPlusP4 + *KMinusP4).M() << " " << (*ProtonP4 + *KPlusP4 + *KMinusP4).M() << endl;

        // Fill histograms: after filter
        hist_Mkk->Fill((*KPlusP4 + *KMinusP4).M(), WeightFactor);
        hist_Mpkk->Fill((*ProtonP4 + *KPlusP4 + *KMinusP4).M(), WeightFactor);
        hist_Mkk_Mpkk->Fill((*KPlusP4 + *KMinusP4).M(), (*ProtonP4 + *KPlusP4 + *KMinusP4).M(), WeightFactor);
    }

    // Save histograms to file
    hist_Mkk->Write();
    hist_Mpkk->Write();
    hist_Mkk_Mpkk->Write();

    output_file->Close();

    return 0;
}
