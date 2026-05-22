#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

double rad_to_deg = 180.0 / TMath::Pi();

class DSelector_phi_d_exc_test_chisqcut : public DSelector
{
public:

    DSelector_phi_d_exc_test_chisqcut(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_phi_d_exc_test_chisqcut(){}

    void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);

private:

    void Get_ComboWrappers(void);
    void Finalize(void);

    // BEAM POLARIZATION INFORMATION
    UInt_t  dPreviousRunNumber;
    bool    dIsPolarizedFlag;
    bool    dIsPARAFlag;
    int     dPolarizationAngle;

    // FLAGS
    bool dIsMC;

    // PARTICLE WRAPPERS
    DParticleComboStep*         dStep0Wrapper;
    DBeamParticle*              dComboBeamWrapper;
    DChargedTrackHypothesis*    dKPlusWrapper;
    DChargedTrackHypothesis*    dKMinusWrapper;
    DChargedTrackHypothesis*    dDeuteronWrapper;

    // CUSTOM HISTOGRAMS
    TH1D* dHist_Chisq_ndf_Before;
    TH1D* dHist_Chisq_ndf_After;
    TH1D* dHist_PhiMass_Before;
    TH1D* dHist_PhiMass_After;

    ClassDef(DSelector_phi_d_exc_test_chisqcut, 0);
};

void DSelector_phi_d_exc_test_chisqcut::Get_ComboWrappers(void)
{
	dStep0Wrapper       = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper   = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper       = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
    dKMinusWrapper      = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dDeuteronWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
}

void DSelector_phi_d_exc_test_chisqcut::Init(TTree *locTree)
{
    // SET OUTPUT FILE NAME
    dOutputFileName          = "selectedhist_phi_d_exc_test_chisqcut.root";
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = "";
    dFlatTreeName            = "";
    dSaveDefaultFlatBranches = true;
    //dSkipNoTriggerEvents     = false;

	// INITIALIZE THE TREE INTERFACE
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;
	Get_ComboWrappers();
    Initialize_Actions();

    // CUSTOM HISTOGRAMS
    dHist_Chisq_ndf_Before = new TH1D("chisq_ndf_before", ";#chi^{2}/NDF;Counts", 100, 0.0, 50.0);
    dHist_Chisq_ndf_After  = new TH1D("chisq_ndf_after", ";#chi^{2}/NDF;Counts", 100, 0.0, 50.0);
    dHist_PhiMass_Before = new TH1D("phi_mass_before", ";m_{K^{+}K^{-}} (GeV/c);Counts", 120, 0.98, 1.1);
    dHist_PhiMass_After  = new TH1D("phi_mass_after", ";m_{K^{+}K^{-}} (GeV/c);Counts", 120, 0.98, 1.1);
}
// END OF INITIALIZATION

Bool_t DSelector_phi_d_exc_test_chisqcut::Process(Long64_t locEntry)
{
	// CALL THIS FIRST
	DSelector::Process(locEntry); // gets the data from the tree for the entry
    Reset_Actions_NewEvent();

    // GET BEAM POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag   = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
        dAnalysisUtilities.Get_PolarizationAngle(locRunNumber, dPolarizationAngle);
		dPreviousRunNumber = locRunNumber;
	}

    // MC INFORMATION
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

    // PRELOOP THE COMBOS TO SORT THE COMBOS BY CHISQ
    set<Int_t> locUsedSoFar_BeamID;
    std::vector<std::pair<UInt_t, Double_t>> loc_combos;
    for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
    {
        dComboWrapper->Set_ComboIndex(loc_i);
        Double_t locChiSq = dComboWrapper->Get_ChiSq_KinFit("");
        loc_combos.push_back(std::make_pair(loc_i, locChiSq));
    }
    std::sort(loc_combos.begin(), loc_combos.end(), [](const std::pair<UInt_t, Double_t>& a, const std::pair<UInt_t, Double_t>& b) { return a.second < b.second;});

    // LOOP OVER COMBOS
    for(const auto& loc_combo : loc_combos)
	{
		// INITIALIZE THE COMBO
        UInt_t loc_i = loc_combo.first;
		dComboWrapper->Set_ComboIndex(loc_i);  // set branch array indices
		if(dComboWrapper->Get_IsComboCut())    // check whether the combo has been cut
			continue;                          // combo has been cut previously

		// GET PARTICLE INDICES
		Int_t locBeamID             = dComboBeamWrapper->Get_BeamID();
        Int_t locKPlusTrackID       = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID      = dKMinusWrapper->Get_TrackID();
		Int_t locDeuteronTrackID    = dDeuteronWrapper->Get_TrackID();

		// GET RECONSTRUCTED P4
        TLorentzVector locBeamP4        = dComboBeamWrapper->Get_P4();
        TLorentzVector locKPlusP4       = dKPlusWrapper->Get_P4();
        TLorentzVector locKMinusP4      = dKMinusWrapper->Get_P4();
        TLorentzVector locDeuteronP4    = dDeuteronWrapper->Get_P4();

        // GET RECONSTRUCTED MEASURED P4
        TLorentzVector locBeamP4_Measured       = dComboBeamWrapper->Get_P4_Measured();
        TLorentzVector locKPlusP4_Measured      = dKPlusWrapper->Get_P4_Measured();
        TLorentzVector locKMinusP4_Measured     = dKMinusWrapper->Get_P4_Measured();
        TLorentzVector locDeuteronP4_Measured   = dDeuteronWrapper->Get_P4_Measured();

        // GET RECONSTRUCTED MEASURED P4 WITH PION MASS HYPOTHESIS
        TLorentzVector locKPlusP4_Measured_PionHypothesis      = TLorentzVector(locKPlusP4_Measured.X(), locKPlusP4_Measured.Y(), locKPlusP4_Measured.Z(), sqrt(locKPlusP4_Measured.P()*locKPlusP4_Measured.P() + 0.13957*0.13957));
        TLorentzVector locKMinusP4_Measured_PionHypothesis     = TLorentzVector(locKMinusP4_Measured.X(), locKMinusP4_Measured.Y(), locKMinusP4_Measured.Z(), sqrt(locKMinusP4_Measured.P()*locKMinusP4_Measured.P() + 0.13957*0.13957));

        double locMissPMinus = (locBeamP4_Measured + TLorentzVector(0,0,0,1.875612859) - locKPlusP4_Measured - locKMinusP4_Measured - locDeuteronP4_Measured).Minus();
        double locVertexR = TMath::Sqrt(dComboBeamWrapper->Get_X4_Measured().X()*dComboBeamWrapper->Get_X4_Measured().X() + dComboBeamWrapper->Get_X4_Measured().Y()*dComboBeamWrapper->Get_X4_Measured().Y());

        // PERFORM CUTS
        if(locBeamP4_Measured.E() < 5.84)                                                                       dComboWrapper->Set_IsComboCut(true);
        if(dDeuteronWrapper->Get_dEdx_CDC()     == 0.0  || dDeuteronWrapper->Get_dEdx_ST()  == 0.0)             dComboWrapper->Set_IsComboCut(true);
        if(dDeuteronWrapper->Get_dEdx_CDC()*1e6 < (TMath::Exp(-3.65*locDeuteronP4_Measured.P())+4.47) + 2.57)   dComboWrapper->Set_IsComboCut(true);
        if((locKPlusP4+locKMinusP4).M()         > 1.1)                                                          dComboWrapper->Set_IsComboCut(true);
        if(locMissPMinus < -0.02)                                                                               dComboWrapper->Set_IsComboCut(true);
        if(locKPlusP4_Measured.P() < 0.40)                                                                      dComboWrapper->Set_IsComboCut(true);
        if(locKMinusP4_Measured.P() < 0.40)                                                                     dComboWrapper->Set_IsComboCut(true);
        if(locDeuteronP4_Measured.P() < 0.40)                                                                   dComboWrapper->Set_IsComboCut(true);
        if(locKPlusP4_Measured.Theta()*rad_to_deg < 2.0)                                                        dComboWrapper->Set_IsComboCut(true);
        if(locKMinusP4_Measured.Theta()*rad_to_deg < 2.0)                                                       dComboWrapper->Set_IsComboCut(true);
        if(locDeuteronP4_Measured.Theta()*rad_to_deg < 2.0)                                                     dComboWrapper->Set_IsComboCut(true);
        if(locVertexR > 1.0)                                                                                    dComboWrapper->Set_IsComboCut(true);
        if(TMath::Abs(dComboBeamWrapper->Get_X4_Measured().Z() - 65.0) > 14.0)                                  dComboWrapper->Set_IsComboCut(true);

        if(dComboWrapper->Get_IsComboCut())  continue;

		// GET THE BEAM ACCIDENTAL WEIGHT FACTOR
		TLorentzVector locBeamX4                       = dComboBeamWrapper->Get_X4_Measured();
		Double_t       locBunchPeriod                  = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t       locDeltaT_RF                    = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4, dComboWrapper);
		Int_t          locRelBeamBucket                = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4, dComboWrapper);          // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t          locNumOutOfTimeBunchesInTree    = 4;                                                                                             // Number of out-of-time beam bunches in tree on a single side
		Bool_t         locSkipNearestOutOfTimeBunch    = false;                                                                                         // true: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t          locNumOutOfTimeBunchesToUse     = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree;
		Double_t       locAccidentalScalingFactor      = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC);         // ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t       locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E());           // ideal value would be 1, but deviations observed, need added factor.
		Double_t       locHistAccidWeightFactor        = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ;        // weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
            locHistAccidWeightFactor = 0.0;
        if(abs(locRelBeamBucket)==5)    // skip the 5th off-time bunch, still saved for systematics studies
            locHistAccidWeightFactor = 0.0;

        // GET THE COMBO ACCIDENTAL WEIGHT FACTOR
        Double_t locComboAccidWeightFactor = 0.0; // default value
        if(locUsedSoFar_BeamID.find(locBeamID) == locUsedSoFar_BeamID.end())
        {
            locComboAccidWeightFactor = 1.0; // combo with best chisq for this beam photon
            locUsedSoFar_BeamID.insert(locBeamID);
        }

        dHist_Chisq_ndf_Before->Fill(dComboWrapper->Get_ChiSq_KinFit("")/dComboWrapper->Get_NDF_KinFit(""), locHistAccidWeightFactor*locComboAccidWeightFactor);
        dHist_PhiMass_Before->Fill((locKPlusP4+locKMinusP4).M(), locHistAccidWeightFactor*locComboAccidWeightFactor);

        if(dComboWrapper->Get_ChiSq_KinFit("")/dComboWrapper->Get_NDF_KinFit("") > 5.0)   // CHISQ CUT
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }

        dHist_Chisq_ndf_After->Fill(dComboWrapper->Get_ChiSq_KinFit("")/dComboWrapper->Get_NDF_KinFit(""), locHistAccidWeightFactor*locComboAccidWeightFactor);
        dHist_PhiMass_After->Fill((locKPlusP4+locKMinusP4).M(), locHistAccidWeightFactor*locComboAccidWeightFactor);

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        // FILL FLAT TREE
        // Fill_FlatTree();
	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_phi_d_exc_test_chisqcut::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION