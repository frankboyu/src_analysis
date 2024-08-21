#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

class DSelector_phi_c_2H_data : public DSelector
{
	public:

		DSelector_phi_c_2H_data(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_phi_c_2H_data(){}

		void   Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool   dIsPolarizedFlag;
		bool   dIsPARAFlag;

		// MONTE CARLO INFORMATION
        bool dIsMC;

		// PARTICLE WRAPPERS
		DParticleComboStep*      dStep0Wrapper;
		DBeamParticle*           dComboBeamWrapper;
        DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dKMinusWrapper;

		// DECLARE HISTOGRAMS
        TH1D* dHist_NumUnusedTracks_Before;
        TH1D* dHist_NumUnusedShowers_Before;
        TH1F* dHist_PhotonEnergy_Before;
        TH1F* dHist_VertexZ_Before;
        TH2F* dHist_VertexXY_Before;
        TH1F* dHist_InvariantMassPhi_Before;

        TH1F* dHist_ConfidenceLevel_After;
        TH1F* dHist_KPlusPIDFOM_After;
        TH1F* dHist_KMinusPIDFOM_After;

        TH1F* dHist_PhotonTiming_Raw;
        TH1F* dHist_PhotonTiming_Weighted;

	ClassDef(DSelector_phi_c_2H_data, 0);
};

void DSelector_phi_c_2H_data::Get_ComboWrappers(void)
{
	dStep0Wrapper     = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper     = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
    dKMinusWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
}

void DSelector_phi_c_2H_data::Init(TTree *locTree)
{
	// SET OUTPUT FILE NAME
	dOutputFileName          = "";
	dOutputTreeFileName      = "";
	dFlatTreeFileName        = "flattree_phi_c_2H_data.root";
	dFlatTreeName            = "flattree_phi_c_2H_data";
    dSaveDefaultFlatBranches = true;

	// INITIALIZE THE TREE INTERFACE AND WRAPPERS
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    Initialize_Actions();

    // DEFINE HISTOGRAMS
    dHist_NumUnusedTracks_Before       = new TH1D("NumUnusedTracks_Before",     ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_Before      = new TH1D("NumUnusedShowers_Before",    ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PhotonEnergy_Before          = new TH1F("PhotonEnergy_Before",        ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_Before               = new TH1F("VertexZ_Before",             ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_Before              = new TH2F("VertexXY_Before",            ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_InvariantMassPhi_Before      = new TH1F("InvariantMassPhi_Before",    ";M_{K^{+}K^{-}} (GeV)        ;Events/0.01 GeV",        500,  0.0,   5.0);

    dHist_ConfidenceLevel_After        = new TH1F("ConfidenceLevel_After",      ";Confidence Level            ;Events/0.001",          1000,  0.0,  0.0001);
    dHist_KPlusPIDFOM_After            = new TH1F("KPlusPIDFOM_After",          ";PIDFOM_{K^{+}}              ;Events/0.001",          1000,  0.0,  0.0001);
    dHist_KMinusPIDFOM_After           = new TH1F("KMinusPIDFOM_After",         ";PIDFOM_{K^{-}}              ;Events/0.001",          1000,  0.0,  0.0001);

    dHist_PhotonTiming_Raw             = new TH1F("PhotonTiming_Raw",           ";#Delta t_{Beam-RF} (ns)     ;Events/0.1 ns",          360,-18.0,  18.0);
    dHist_PhotonTiming_Weighted        = new TH1F("PhotonTiming_Weighted",      ";#Delta t_{Beam-RF} (ns)     ;Events/0.1 ns",          360,-18.0,  18.0);

    // FLAT TREE BRANCHES
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");
}
// END OF INITIALIZATION

Bool_t DSelector_phi_c_2H_data::Process(Long64_t locEntry)
{
	// CALL THIS FIRST
	DSelector::Process(locEntry); // gets the data from the tree for the entry
    Reset_Actions_NewEvent();

    // GET BEAM POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag   = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

    // MC INFORMATION
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

	// LOOP OVER COMBOS
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		// INITIALIZE THE COMBO
		dComboWrapper->Set_ComboIndex(loc_i);  // set branch array indices
		if(dComboWrapper->Get_IsComboCut())    // check whether the combo has been cut
			continue;                          // combo has been cut previously

		// GET PARTICLE INDICES
		Int_t locBeamID         = dComboBeamWrapper->Get_BeamID();
        Int_t locKPlusTrackID   = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID  = dKMinusWrapper->Get_TrackID();

		// GET RECONSTRUCTED P4
        TLorentzVector locBeamP4     = dComboBeamWrapper->Get_P4_Measured();
        TLorentzVector locKPlusP4    = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4   = dKMinusWrapper->Get_P4_Measured();

        // FILL HISTOGRAMS
        dHist_NumUnusedTracks_Before    ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_Before   ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_Before       ->Fill(locBeamP4.E());
        dHist_VertexZ_Before            ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_Before           ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_InvariantMassPhi_Before   ->Fill((locKPlusP4+locKMinusP4).M());

        // PERFORM CUTS
        if(dComboWrapper->Get_NumUnusedTracks()        > 0)                                                 dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_NumUnusedShowers()       > 0)                                                 dComboWrapper->Set_IsComboCut(true);
        if(locBeamP4.E()                               < 5.8  || locBeamP4.E()                   > 10.7)    dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Z()             < 51.0 || dComboBeamWrapper->Get_X4().Z() > 79.0)    dComboWrapper->Set_IsComboCut(true);
        if(sqrt(pow(dComboBeamWrapper->Get_X4().X(),2) + pow(dComboBeamWrapper->Get_X4().Y(),2)) > 1.0)     dComboWrapper->Set_IsComboCut(true);
        if((locKPlusP4+locKMinusP4).M()                > 1.20)                                              dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < 1e-5)                                              dComboWrapper->Set_IsComboCut(true);
        if(dKPlusWrapper->Get_PIDFOM()                 < 1e-5)                                              dComboWrapper->Set_IsComboCut(true);
        if(dKMinusWrapper->Get_PIDFOM()                < 1e-5)                                              dComboWrapper->Set_IsComboCut(true);

        if(dComboWrapper->Get_IsComboCut())  continue;

        // FILL HISTOGRAMS
        dHist_ConfidenceLevel_After ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        dHist_KPlusPIDFOM_After     ->Fill(dKPlusWrapper->Get_PIDFOM());
        dHist_KMinusPIDFOM_After    ->Fill(dKMinusWrapper->Get_PIDFOM());

		// GET THE ACCIDENTAL WEIGHT FACTOR
		TLorentzVector locBeamX4                       = dComboBeamWrapper->Get_X4_Measured();
		Double_t       locBunchPeriod                  = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t       locDeltaT_RF                    = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4, dComboWrapper);
		Int_t          locRelBeamBucket                = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t          locNumOutOfTimeBunchesInTree    = 4;                                                                                             // Number of out-of-time beam bunches in tree on a single side
		Bool_t         locSkipNearestOutOfTimeBunch    = true;                                                                                          // true: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t          locNumOutOfTimeBunchesToUse     = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree;
		Double_t       locAccidentalScalingFactor      = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC);         // ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t       locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E());           // ideal value would be 1, but deviations observed, need added factor.
		Double_t       locHistAccidWeightFactor        = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ;        // weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time

        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)                                                                                    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
        {
            dComboWrapper->Set_IsComboCut(true);
			continue;
        }

        dHist_PhotonTiming_Raw->Fill(locDeltaT_RF);
        dHist_PhotonTiming_Weighted->Fill(locDeltaT_RF, locHistAccidWeightFactor);

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight", locHistAccidWeightFactor);
        Fill_FlatTree(); //for the active combo
	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_phi_c_2H_data::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION