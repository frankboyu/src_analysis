#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

double RadToDeg     = 180.0/3.1415926;

class DSelector_piminus_p_recon : public DSelector
{
	public:

		DSelector_piminus_p_recon(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_piminus_p_recon(){}

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
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// CUSTOM HISTOGRAMS
        TH1F* dHist_NumUnusedTracks_Before;
        TH1F* dHist_NumUnusedShowers_Before;
        TH1F* dHist_PhotonEnergy_Before;
        TH1F* dHist_VertexZ_Before;
        TH2F* dHist_VertexXY_Before;
        TH1F* dHist_ConfidenceLevel_Before;
        TH1F* dHist_PiMinusPIDFOM_Before;
        TH1F* dHist_ProtonPIDFOM_Before;
        TH1F* dHist_MissingPMinus_Before;

        TH1F* dHist_NumUnusedTracks_After;
        TH1F* dHist_NumUnusedShowers_After;
        TH1F* dHist_PhotonEnergy_After;
        TH1F* dHist_VertexZ_After;
        TH2F* dHist_VertexXY_After;
        TH1F* dHist_ConfidenceLevel_After;
        TH1F* dHist_PiMinusPIDFOM_After;
        TH1F* dHist_ProtonPIDFOM_After;
        TH1F* dHist_MissingPMinus_After;

        TH1F* dHist_PhotonTiming_Raw;
        TH1F* dHist_PhotonTiming_Weighted;

        TH1F* dHist_NumUnusedTracks_Weighted;
        TH1F* dHist_NumUnusedShowers_Weighted;
        TH1F* dHist_PhotonEnergy_Weighted;
        TH1F* dHist_VertexZ_Weighted;
        TH2F* dHist_VertexXY_Weighted;
        TH1F* dHist_ConfidenceLevel_Weighted;
        TH1F* dHist_PiMinusPIDFOM_Weighted;
        TH1F* dHist_ProtonPIDFOM_Weighted;
        TH1F* dHist_MissingPMinus_Weighted;

	ClassDef(DSelector_piminus_p_recon, 0);
};

void DSelector_piminus_p_recon::Get_ComboWrappers(void)
{
	dStep0Wrapper     = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiMinusWrapper   = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dProtonWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
}

void DSelector_piminus_p_recon::Init(TTree *locTree)
{
	// SET OUTPUT FILE NAME
	dOutputFileName          = "";
	dOutputTreeFileName      = "";
	dFlatTreeName            = "flattree_piminus_p_recon";
    dSaveDefaultFlatBranches = true;

    if      (locTree->GetBranch("MCWeight") == NULL && !strncmp(locTree->GetName(), "gd_pimprotmissprot__B4_F4_T1_S4_Tree",   strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_data_2H_missprot.root";
    else if (locTree->GetBranch("MCWeight") == NULL && !strncmp(locTree->GetName(), "gd_pimprotinc__B4_F4_T1_S4_Tree",        strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_data_2H_inc.root";
    else if (locTree->GetBranch("MCWeight") == NULL && !strncmp(locTree->GetName(), "ghe_pimprotinc__B4_F4_T1_S4_Tree",       strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_data_4He_inc.root";
    else if (locTree->GetBranch("MCWeight") == NULL && !strncmp(locTree->GetName(), "gc12_pimprotinc__B4_F4_T2_S5_Tree",      strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_data_12C_inc.root";
    else if (locTree->GetBranch("MCWeight") != NULL && !strncmp(locTree->GetName(), "gd_pimprotmissprot__B4_F4_T1_S4_Tree",   strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_sim_2H_missprot.root";
    else if (locTree->GetBranch("MCWeight") != NULL && !strncmp(locTree->GetName(), "gd_pimprotinc__B4_F4_T1_S4_Tree",        strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_sim_2H_inc.root";
    else if (locTree->GetBranch("MCWeight") != NULL && !strncmp(locTree->GetName(), "ghe_pimprotinc__B4_F4_T1_S4_Tree",       strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_sim_4He_inc.root";
    else if (locTree->GetBranch("MCWeight") != NULL && !strncmp(locTree->GetName(), "gc12_pimprotinc__B4_F4_T2_S5_Tree",      strlen(locTree->GetName())))
        dFlatTreeFileName = "flattree_piminus_p_recon_sim_12C_inc.root";

	// INITIALIZE THE TREE INTERFACE AND WRAPPERS
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    Initialize_Actions();

    // CUSTOM HISTOGRAMS
    dHist_NumUnusedTracks_Before    = new TH1F("NumUnusedTracks_Before",    ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_Before   = new TH1F("NumUnusedShowers_Before",   ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PhotonEnergy_Before       = new TH1F("PhotonEnergy_Before",       ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_Before            = new TH1F("VertexZ_Before",            ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_Before           = new TH2F("VertexXY_Before",           ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_ConfidenceLevel_Before    = new TH1F("ConfidenceLevel_Before",    ";Confidence Level            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_PiMinusPIDFOM_Before      = new TH1F("PiMinusPIDFOM_Before",      ";PIDFOM_{#pi^{-}}            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_Before       = new TH1F("ProtonPIDFOM_Before",       ";PIDFOM_{p}                  ;Events/0.001",          1000,  0.0,   1.0);
    dHist_MissingPMinus_Before      = new TH1F("MissingPMinus_Before",      ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);

    dHist_NumUnusedTracks_After     = new TH1F("NumUnusedTracks_After",     ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_After    = new TH1F("NumUnusedShowers_After",    ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PhotonEnergy_After        = new TH1F("PhotonEnergy_After",        ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_After             = new TH1F("VertexZ_After",             ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_After            = new TH2F("VertexXY_After",            ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_ConfidenceLevel_After     = new TH1F("ConfidenceLevel_After",     ";Confidence Level            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_PiMinusPIDFOM_After       = new TH1F("PiMinusPIDFOM_After",       ";PIDFOM_{#pi^{-}}            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_After        = new TH1F("ProtonPIDFOM_After",        ";PIDFOM_{p}                  ;Events/0.001",          1000,  0.0,   1.0);
    dHist_MissingPMinus_After       = new TH1F("MissingPMinus_After",       ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);

    dHist_PhotonTiming_Raw          = new TH1F("PhotonTiming_Raw",          ";#Delta t_{Beam-RF} (ns)     ;Events/0.1 ns",          360,-18.0,  18.0);
    dHist_PhotonTiming_Weighted     = new TH1F("PhotonTiming_Weighted",     ";#Delta t_{Beam-RF} (ns)     ;Events/0.1 ns",          360,-18.0,  18.0);

    dHist_NumUnusedTracks_Weighted  = new TH1F("NumUnusedTracks_Weighted",  ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_Weighted = new TH1F("NumUnusedShowers_Weighted", ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PhotonEnergy_Weighted     = new TH1F("PhotonEnergy_Weighted",     ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_Weighted          = new TH1F("VertexZ_Weighted",          ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_Weighted         = new TH2F("VertexXY_Weighted",         ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_ConfidenceLevel_Weighted  = new TH1F("ConfidenceLevel_Weighted",  ";Confidence Level            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_PiMinusPIDFOM_Weighted    = new TH1F("PiMinusPIDFOM_Weighted",    ";PIDFOM_{#pi^{-}}            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_Weighted     = new TH1F("ProtonPIDFOM_Weighted",     ";PIDFOM_{p}                  ;Events/0.001",          1000,  0.0,   1.0);
    dHist_MissingPMinus_Weighted    = new TH1F("MissingPMinus_Weighted",    ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pim_pidfom");  // the PIDFOM in the default flat branches pim_pid_fom is corrupted and always 0
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("p_pidfom");    // the PIDFOM in the default flat branches p_pid_fom is corrupted and always 0
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("p_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("p_p4_truth");

}
// END OF INITIALIZATION

Bool_t DSelector_piminus_p_recon::Process(Long64_t locEntry)
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

        // DISCARD EVENTS WITH NO L1 TRIGGER BITS
        if (dComboWrapper->Get_L1TriggerBits() == 0)
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }

		// GET PARTICLE INDICES
		Int_t locBeamID         = dComboBeamWrapper->Get_BeamID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID  = dProtonWrapper->Get_TrackID();

		// GET RECONSTRUCTED P4
		TLorentzVector locBeamP4    = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4  = dProtonWrapper->Get_P4();

        //GET THROWN P4
        TLorentzVector locBeamX4_Thrown, locPiMinusX4_Thrown, locProtonX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locPiMinusP4_Thrown, locProtonP4_Thrown;
        if(dThrownBeam != NULL)
        {
            locBeamX4_Thrown = dThrownBeam->Get_X4();
            locBeamP4_Thrown = dThrownBeam->Get_P4();
        }
        for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
        {
            dThrownWrapper->Set_ArrayIndex(loc_i);
            if (dThrownWrapper->Get_PID() == PiMinus)
            {
                locPiMinusX4_Thrown = dThrownWrapper->Get_X4();
                locPiMinusP4_Thrown = dThrownWrapper->Get_P4();
            }
            else if (dThrownWrapper->Get_PID() == Proton)
            {
                locProtonX4_Thrown = dThrownWrapper->Get_X4();
                locProtonP4_Thrown = dThrownWrapper->Get_P4();
            }
        }

        // FILL HISTOGRAMS BEFORE CUTS
        dHist_NumUnusedTracks_Before  ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_Before ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_Before     ->Fill(locBeamP4.E());
        dHist_VertexZ_Before          ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_Before         ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_ConfidenceLevel_Before  ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        dHist_PiMinusPIDFOM_Before    ->Fill(dPiMinusWrapper->Get_PIDFOM());
        dHist_ProtonPIDFOM_Before     ->Fill(dProtonWrapper->Get_PIDFOM());
        dHist_MissingPMinus_Before    ->Fill((locPiMinusP4 + locProtonP4 - locBeamP4).Minus());

        // PERFORM CUTS
        if(dComboWrapper->Get_NumUnusedTracks()             > 0)                                                                dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_NumUnusedShowers()            > 0)                                                                dComboWrapper->Set_IsComboCut(true);
        if(locBeamP4.E()                                    < 5.8  || locBeamP4.E()                   > 10.7)                   dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Z()                  < 51.0 || dComboBeamWrapper->Get_X4().Z() > 79.0)                   dComboWrapper->Set_IsComboCut(true);
        if(sqrt(pow(dComboBeamWrapper->Get_X4().X(),2)      + pow(dComboBeamWrapper->Get_X4().Y(),2)) > 1.0)                    dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_ConfidenceLevel_KinFit()      < 1e-4)                                                             dComboWrapper->Set_IsComboCut(true);
        if(dPiMinusWrapper->Get_PIDFOM()                    < 1e-4)                                                             dComboWrapper->Set_IsComboCut(true);
        if(dProtonWrapper->Get_PIDFOM()                     < 1e-4)                                                             dComboWrapper->Set_IsComboCut(true);
        if((locPiMinusP4 + locProtonP4 - locBeamP4).Minus() < 0.4  || (locPiMinusP4 + locProtonP4 - locBeamP4).Minus() > 1.4)   dComboWrapper->Set_IsComboCut(true);

        if(dComboWrapper->Get_IsComboCut())  continue;

        // FILL HISTOGRAMS AFTER CUTS
        dHist_NumUnusedTracks_After  ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_After ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_After     ->Fill(locBeamP4.E());
        dHist_VertexZ_After          ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_After         ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_ConfidenceLevel_After  ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        dHist_PiMinusPIDFOM_After    ->Fill(dPiMinusWrapper->Get_PIDFOM());
        dHist_ProtonPIDFOM_After     ->Fill(dProtonWrapper->Get_PIDFOM());
        dHist_MissingPMinus_After    ->Fill((locPiMinusP4 + locProtonP4 - locBeamP4).Minus());

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

        // FILL WEIGHTED HISTOGRAMS
        dHist_NumUnusedTracks_Weighted  ->Fill(dComboWrapper->Get_NumUnusedTracks(), locHistAccidWeightFactor);
        dHist_NumUnusedShowers_Weighted ->Fill(dComboWrapper->Get_NumUnusedShowers(), locHistAccidWeightFactor);
        dHist_PhotonEnergy_Weighted     ->Fill(locBeamP4.E(), locHistAccidWeightFactor);
        dHist_VertexZ_Weighted          ->Fill(dComboBeamWrapper->Get_X4().Z(), locHistAccidWeightFactor);
        dHist_VertexXY_Weighted         ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y(), locHistAccidWeightFactor);
        dHist_ConfidenceLevel_Weighted  ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(), locHistAccidWeightFactor);
        dHist_PiMinusPIDFOM_Weighted    ->Fill(dPiMinusWrapper->Get_PIDFOM(), locHistAccidWeightFactor);
        dHist_ProtonPIDFOM_Weighted     ->Fill(dProtonWrapper->Get_PIDFOM(), locHistAccidWeightFactor);
        dHist_MissingPMinus_Weighted    ->Fill((locPiMinusP4 + locProtonP4 - locBeamP4).Minus(), locHistAccidWeightFactor);


		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight", locHistAccidWeightFactor);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("pim_pidfom", dPiMinusWrapper->Get_PIDFOM());
		dFlatTreeInterface->Fill_Fundamental<Double_t>("p_pidfom", dProtonWrapper->Get_PIDFOM());
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_x4_truth", locBeamX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_x4_truth", locPiMinusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("p_x4_truth", locProtonX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_p4_truth", locBeamP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_p4_truth", locPiMinusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("p_p4_truth", locProtonP4_Thrown);
        Fill_FlatTree(); //for the active combo
	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_piminus_p_recon::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION


