using namespace std;
using namespace ROOT;
using namespace RooFit;

double mass_proton      = 0.938272;
double mass_deuteron    = 1.875612;
double mass_pion        = 0.13957018;
double mass_phi         = 1.019461;

int filter_phi_c_2H_data()
{
    // Open the file and get the tree
    TChain *chain = new TChain("flattree_phi_c_2H_data");
    chain->Add("/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_c_2H_data_ver01/*.root");
    TFile *output_file = new TFile("output/plots_phi_c_2H_data.root", "RECREATE");

    // Set branch addresses
    TLorentzVector *BeamP4      = new TLorentzVector();
    TLorentzVector *KPlusP4     = new TLorentzVector();
    TLorentzVector *KMinusP4    = new TLorentzVector();
    TLorentzVector *MissingP4   = new TLorentzVector();
    double WeightFactor;
    int RunNumber, Entry, Combo;

    chain->SetBranchAddress("BeamP4", &BeamP4);
    chain->SetBranchAddress("KPlusP4", &KPlusP4);
    chain->SetBranchAddress("KMinusP4", &KMinusP4);
    chain->SetBranchAddress("MissingP4", &MissingP4);
    chain->SetBranchAddress("WeightFactor", &WeightFactor);
    chain->SetBranchAddress("RunNumber", &RunNumber);
    chain->SetBranchAddress("Entry", &Entry);
    chain->SetBranchAddress("Combo", &Combo);

    // Create histogram
    TH1F *hist_MPhi                 = new TH1F("hist_MPhi", "hist_MPhi", 400, 0.9, 1.3);
    TH1F *hist_MissingMass          = new TH1F("hist_MissingMass", "hist_MissingMass", 200, 1.0, 3.0);
    TH2F *hist_MPhi_ThetaPhi        = new TH2F("hist_MPhi_ThetaPhi", "hist_MPhi_ThetaPhi", 400, 0.9, 1.3, 200, 0.0, 20.0);
    TH2F *hist_MPhi_MinusT          = new TH2F("hist_MPhi_MinusT", "hist_MPhi_MinusT", 400, 0.9, 1.3, 200, 0.0, 2.0);
    TH2F *hist_MPhi_yPhi            = new TH2F("hist_MPhi_yPhi", "hist_MPhi_yPhi", 400, 0.9, 1.3, 200, 0.0, 2.0);
    TH2F *hist_MPhi_MissingMass     = new TH2F("hist_MPhi_MissingMass", "hist_MPhi_MissingMass", 400, 0.9, 1.3, 200, 1.0, 3.0);
    TH2F *hist_MPhi_DeltaE          = new TH2F("hist_MPhi_DeltaE", "hist_MPhi_DeltaE", 400, 0.9, 1.3, 200, -2.0, 2.0);
    TH2F *hist_MPhi_MRho            = new TH2F("hist_MPhi_MRho", "hist_MPhi_MRho", 400, 0.9, 1.3, 800, 0.3, 1.1);
    TH2F *hist_ThetaPhi_yPhi        = new TH2F("hist_ThetaPhi_yPhi", "hist_ThetaPhi_yPhi", 200, 0.0, 20.0, 200, 0.0, 2.0);
    TH2F *hist_MissingMass_yPhi     = new TH2F("hist_MissingMass_yPhi", "hist_MissingMass_yPhi", 200, 1.0, 3.0, 200, 0.0, 2.0);
    TH2F *hist_DeltaE_yPhi          = new TH2F("hist_DeltaE_yPhi", "hist_DeltaE_yPhi", 200, -2.0, 2.0, 200, 0.0, 2.0);
    TH2F *hist_DeltaE_MissingMass   = new TH2F("hist_DeltaE_MissingMass", "hist_DeltaE_MissingMass", 200, -2.0, 2.0, 200, 1.0, 3.0);

    // Loop over tree entries
    for (Long64_t i = 0; i < chain->GetEntries(); i++)
    // for (Long64_t i = 0; i < 10000; i++)
    {
        chain->GetEntry(i);

        // Calculate variables
        TLorentzVector *TargetP4 = new TLorentzVector(0.0, 0.0, 0.0, mass_deuteron);
        TLorentzVector *PhiP4 = new TLorentzVector(*KPlusP4 + *KMinusP4);
        TLorentzVector *PhiP4COM = new TLorentzVector(*KPlusP4 + *KMinusP4);
        PhiP4COM->Boost(-(*BeamP4+*TargetP4).BoostVector());
        TLorentzVector *KPlusP4AsPion = new TLorentzVector(KPlusP4->X(), KPlusP4->Y(), KPlusP4->Z(), sqrt(KPlusP4->P()*KPlusP4->P() + mass_pion*mass_pion));
        TLorentzVector *KMinusP4AsPion = new TLorentzVector(KMinusP4->X(), KMinusP4->Y(), KMinusP4->Z(), sqrt(KMinusP4->P()*KMinusP4->P() + mass_pion*mass_pion));

        double sqrt_s               = (*BeamP4 + *MissingP4).Mag();
        double minus_t              = -(*BeamP4 - *PhiP4).Mag2();
        double minus_u              = -(*BeamP4 - *MissingP4).Mag2();
        double phi_mass             = (*KPlusP4 + *KMinusP4).M();
        double phi_energy_expected  = (pow(sqrt_s, 2) - pow(TargetP4->M(), 2) + pow(mass_phi, 2)) / (2. * sqrt_s);
        double rho_mass             = (*KPlusP4AsPion + *KMinusP4AsPion).M();
        double phi_theta            = (*KPlusP4 + *KMinusP4).Theta() * 180.0 / TMath::Pi();
        double y_phi                = minus_t/(2*mass_deuteron*(BeamP4->E()-PhiP4->E()));

        // Fill histograms: before filter

        // Filter events
        // if (y_phi < 0.6) continue;
        hist_MPhi->Fill(phi_mass, WeightFactor);
        // if (phi_mass < 1.01 || phi_mass > 1.03) continue;

        // Fill histograms: after filter
        hist_MissingMass->Fill(MissingP4->M(), WeightFactor);
        hist_MPhi_ThetaPhi->Fill(phi_mass, phi_theta, WeightFactor);
        hist_MPhi_MinusT->Fill(phi_mass, minus_t, WeightFactor);
        hist_MPhi_yPhi->Fill(phi_mass, y_phi, WeightFactor);
        hist_MPhi_MissingMass->Fill(phi_mass, MissingP4->M(), WeightFactor);
        hist_MPhi_DeltaE->Fill(phi_mass, phi_energy_expected - PhiP4COM->E(), WeightFactor);
        hist_MPhi_MRho->Fill(phi_mass, rho_mass, WeightFactor);
        hist_ThetaPhi_yPhi->Fill(phi_theta, y_phi, WeightFactor);
        hist_MissingMass_yPhi->Fill(MissingP4->M(), y_phi, WeightFactor);
        hist_DeltaE_yPhi->Fill(phi_energy_expected - PhiP4COM->E(), y_phi, WeightFactor);
        hist_DeltaE_MissingMass->Fill(phi_energy_expected - PhiP4COM->E(), MissingP4->M(), WeightFactor);
    }

    // Save histograms to file
    hist_MPhi->Write();
    hist_MissingMass->Write();
    hist_MPhi_ThetaPhi->Write();
    hist_MPhi_MinusT->Write();
    hist_MPhi_yPhi->Write();
    hist_MPhi_MissingMass->Write();
    hist_MPhi_DeltaE->Write();
    hist_MPhi_MRho->Write();
    hist_ThetaPhi_yPhi->Write();
    hist_MissingMass_yPhi->Write();
    hist_DeltaE_yPhi->Write();
    hist_DeltaE_MissingMass->Write();
    output_file->Close();

    return 0;
}
