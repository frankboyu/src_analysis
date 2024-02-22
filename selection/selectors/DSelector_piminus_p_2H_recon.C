#include "DSelector_piminus_p_recon.h"

double RadToDeg = 180.0/3.1415926;

void DSelector_piminus_p_recon::Init(TTree *locTree)
{
	// SET OUTPUT FILE NAME
	dOutputFileName          = "";                                 // "" for none
	dOutputTreeFileName      = "";                                 // "" for none
	dFlatTreeFileName        = "flattree.root";                    // output flat tree (one combo per tree entry), "" for none
	dFlatTreeName            = "piminus_p_recon";                  // if blank, default name will be chosen
    dSaveDefaultFlatBranches = false;                              // False: don't save default branches, reduce disk footprint.

	// INITIALIZE THE TREE INTERFACE AND WRAPPERS
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    Initialize_Actions();
    
    // MC INFORMATION
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

    // CUSTOM HISTOGRAMS: BEFORE THE CUTS
    dHist_NumUnusedTracks_Measured_Before         = new TH1F("NumUnusedTracks_Measured_Before",      ";Unused Tracks                  ;Events/1",             10,  0.0,  10.0);
    dHist_NumUnusedShowers_Measured_Before        = new TH1F("NumUnusedShowers_Measured_Before",     ";Unused Showers                 ;Events/1",             10,  0.0,  10.0);
    dHist_PiMinusPIDFOM_Measured_Before           = new TH1F("PiMinusPIDFOM_Measured_Before",        ";PID FOM                        ;Events/0.001",       1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_Measured_Before            = new TH1F("ProtonPIDFOM_Measured_Before",         ";PID FOM                        ;Events/0.001",       1000,  0.0,   1.0);
    dHist_PhotonEnergy_Measured_Before            = new TH1F("PhotonEnergy_Before",                  ";Photon Energy (GeV)            ;Events/0.01 GeV",     900,  3.0,  12.0);
    dHist_ConfidenceLevel_KinFit_Before           = new TH1F("ConfidenceLevel_KinFit_Before",        ";Confidence Level               ;Events/0.001",       1000,  0.0,   1.0);
    dHist_VertexZ_KinFit_Before                   = new TH1F("VertexZ_KinFit_Before",                ";Vertex Z (cm)                  ;Events/1 cm",         200,  0.0, 200.0);
    dHist_VertexXY_KinFit_Before                  = new TH2F("VertexXY_KinFit_Before",               ";Vertex X (cm)                  ;Vertex Y (cm)",       100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_MissingPMinus_KinFit_Before             = new TH1F("MissingPMinus_KinFit_Before",          ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);        
    dHist_MissingMomentum_KinFit_Before           = new TH1F("MissingMomentum_KinFit_Before",        ";Missing Momentum (GeV)         ;Events/0.01 GeV",     200,  0.0,   2.0);
    
    // CUSTOM HISTOGRAMS: CUT EFFECTS
    dHist_CutEffect_NumUnusedTracks               = new TH1F("CutEffect_NumUnusedTracks",            ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_NumUnusedShowers              = new TH1F("CutEffect_NumUnusedShowers",           ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_TrackPIDFOM                   = new TH1F("CutEffect_TrackPIDFOM",                ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_PhotonEnergy                  = new TH1F("CutEffect_PhotonEnergy",               ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_ConfidenceLevel               = new TH1F("CutEffect_ConfidenceLevel",            ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_CommonVertex                  = new TH1F("CutEffect_CommonVertex",               ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_MissingMomentum               = new TH1F("CutEffect_MissingMomentum",            ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     200,  0.0,   2.0);
    
    // CUSTOM HISTOGRAMS: AFTER THE CUTS
    dHist_NumUnusedTracks_Measured_After          = new TH1F("NumUnusedTracks_Measured_After",       ";Unused Tracks                  ;Events/1",              10,  0.0,  10.0);
    dHist_NumUnusedShowers_Measured_After         = new TH1F("NumUnusedShowers_Measured_After",      ";Unused Showers                 ;Events/1",              10,  0.0,  10.0);
    dHist_PiMinusPIDFOM_Measured_After            = new TH1F("PiMinusPIDFOM_Measured_After",         ";PID FOM                        ;Events/0.001",        1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_Measured_After             = new TH1F("ProtonPIDFOM_Measured_After",          ";PID FOM                        ;Events/0.001",        1000,  0.0,   1.0);
    dHist_PiMinusPVsdEdx_Measured_After           = new TH2F("PiMinusPVsdEdx_Measured_After",        ";P (GeV)                        ;dE/dx (keV/cm)",       100,  0.0,  10.0,  250,  0.0,  25.0);
    dHist_ProtonPVsdEdx_Measured_After            = new TH2F("ProtonPVsdEdx_Measured_After",         ";P (GeV)                        ;dE/dx (keV/cm)",       100,  0.0,  10.0,  250,  0.0,  25.0);    
    dHist_PhotonEnergy_Measured_After             = new TH1F("PhotonEnergy_After",                   ";Photon Energy (GeV)            ;Events/0.01 GeV",      900,  3.0,  12.0);
    dHist_PhotonTiming_Measured_After             = new TH1F("PhotonTiming_After",                   ";#Delta t_{Beam-RF} (ns)        ;Events/0.1 ns",        360,-18.0,  18.0);
    dHist_ConfidenceLevel_KinFit_After            = new TH1F("ConfidenceLevel_KinFit_After",         ";Confidence Level               ;Events/0.001",        1000,  0.0,   1.0);
    dHist_ChiSquarePerNDF_KinFit_After            = new TH1F("ChiSquarePerNDF_KinFit_After",         ";#chi^{2}/NDF;                  ;Events/0.01",         1000,  0.0,  10.0);
    dHist_VertexZ_KinFit_After                    = new TH1F("VertexZ_KinFit_After",                 ";Vertex Z (cm)                  ;Events/1 cm",          200,  0.0, 200.0);
    dHist_VertexXY_KinFit_After                   = new TH2F("VertexXY_KinFit_After",                ";Vertex X (cm)                  ;Vertex Y (cm)",        100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_MissingPMinus_Measured_After            = new TH1F("MissingPMinus_Measured_After",         ";Missing P^{-} (GeV)            ;Events/0.01 GeV",      200,  0.0,   2.0);
    dHist_MissingPMinus_KinFit_After              = new TH1F("MissingPMinus_KinFit_After",           ";Missing P^{-} (GeV)            ;Events/0.01 GeV",      200,  0.0,   2.0);        
    dHist_MissingMassSquared_Measured_After       = new TH1F("MissingMassSquared_Measured_After",    ";Missing Mass Squared (GeV^{2}) ;Events/0.01 GeV^{2}",  600,  0.0,   6.0);
    dHist_MissingMassSquared_KinFit_After         = new TH1F("MissingMassSquared_KinFit_After",      ";Missing Mass Squared (GeV^{2}) ;Events/0.01 GeV^{2}",  600,  0.0,   6.0);
    dHist_MissingMomentum_Measured_After          = new TH1F("MissingMomentum_Measured_After",       ";Missing Momentum (GeV)         ;Events/0.01 GeV",      200,  0.0,   2.0);
    dHist_MissingMomentum_KinFit_After            = new TH1F("MissingMomentum_KinFit_After",         ";Missing Momentum (GeV)         ;Events/0.01 GeV",      200,  0.0,   2.0);
    dHist_SqrtS_Measured_After                    = new TH1F("SqrtS_Measured_After",                 ";#sqrt(s) (GeV)                 ;Events/0.01 GeV",     1000,  0.0,  10.0);
    dHist_SqrtS_KinFit_After                      = new TH1F("SqrtS_KinFit_After",                   ";#sqrt(s) (GeV)                 ;Events/0.01 GeV",     1000,  0.0,  10.0);
    dHist_MinusT_Measured_After                   = new TH1F("MinusT_Measured_After",                ";-t (GeV^{2})                   ;Events/0.01 GeV^{2}", 1500,  0.0,  15.0);
    dHist_MinusT_KinFit_After                     = new TH1F("MinusT_KinFit_After",                  ";-t (GeV^{2})                   ;Events/0.01 GeV^{2}", 1500,  0.0,  15.0);
    dHist_MinusU_Measured_After                   = new TH1F("MinusU_Measured_After",                ";-u (GeV^{2})                   ;Events/0.01 GeV^{2}", 1500,  0.0,  15.0);
    dHist_MinusU_KinFit_After                     = new TH1F("MinusU_KinFit_After",                  ";-u (GeV^{2})                   ;Events/0.01 GeV^{2}", 1500,  0.0,  15.0);
    dHist_ThetaCM_Measured_After                  = new TH1F("ThetaCM_Measured_After",               ";#theta_{c.m.} (deg)            ;Events/1 deg",         180,  0.0, 180.0);
    dHist_ThetaCM_KinFit_After                    = new TH1F("ThetaCM_KinFit_After",                 ";#theta_{c.m.} (deg)            ;Events/1 deg",         180,  0.0, 180.0);
    dHist_Coplanarity_Measured_After              = new TH1F("Coplanarity_Measured_After",           ";Coplanarity Angle (deg)        ;Events/1 deg",         360,  0.0, 360.0);
    dHist_Coplanarity_KinFit_After                = new TH1F("Coplanarity_KinFit_After",             ";Coplanarity Angle (deg)        ;Events/1 deg",         360,  0.0, 360.0);
    dHist_PiMinusPVsTheta_Measured_After          = new TH2F("PiMinusPVsTheta_Measured_After",       ";P (GeV)                        ;#theta (deg)",        1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_PiMinusPVsTheta_KinFit_After            = new TH2F("PiMinusPVsTheta_KinFit_After",         ";P (GeV)                        ;#theta (deg)",        1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ProtonPVsTheta_Measured_After           = new TH2F("ProtonPVsTheta_Measured_After",        ";P (GeV)                        ;#theta (deg)",        1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ProtonPVsTheta_KinFit_After             = new TH2F("ProtonPVsTheta_KinFit_After",          ";P (GeV)                        ;#theta (deg)",        1000,  0.0,  10.0,  180,  0.0, 180.0);
    
    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("RunNumber");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("Entry");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("Combo");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("L1TriggerBits");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("HistAccidWeightFactor");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_Measured");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_Measured");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_Measured");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_Measured");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_KinFit");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_KinFit");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_KinFit");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_KinFit");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_Thrown");
} 
// END OF INITIALIZATION

Bool_t DSelector_piminus_p_recon::Process(Long64_t locEntry)
{
	// CALL THIS FIRST
	DSelector::Process(locEntry); // gets the data from the tree for the entry
    Reset_Actions_NewEvent();

    // GET PHOTON POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag   = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	// LOOP OVER COMBOS
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		// INITIALIZE THE COMBO
		dComboWrapper->Set_ComboIndex(loc_i);  // set branch array indices
		if(dComboWrapper->Get_IsComboCut())    // check whether the combo has been cut
			continue;                          // combo has been cut previously

		// GET PARTICLE INDICES
		Int_t locBeamID         = dComboBeamWrapper->Get_BeamID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID  = dProtonWrapper->Get_TrackID();
        
		// GET MEASURED P4
		TLorentzVector locBeamP4_Measured    = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured  = dProtonWrapper->Get_P4_Measured();
        TLorentzVector locMissingP4_Measured = locPiMinusP4_Measured + locProtonP4_Measured - locBeamP4_Measured;
		
        // GET KINFIT P4
        TLorentzVector locBeamP4_KinFit      = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiMinusP4_KinFit   = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4_KinFit    = dProtonWrapper->Get_P4();
        TLorentzVector locMissingP4_KinFit   = locPiMinusP4_KinFit + locProtonP4_KinFit - locBeamP4_KinFit;
        
        // GET THROWN P4
        TLorentzVector locBeamP4_Thrown, locPiMinusP4_Thrown, locProtonP4_Thrown, locMissingP4_Thrown;
        if(dThrownBeam != NULL)
            locBeamP4_Thrown = dThrownBeam->Get_P4();
        for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
        {                                                                                                                                                                                                                                                     
            dThrownWrapper->Set_ArrayIndex(loc_i);
            if (dThrownWrapper->Get_PID() == PiMinus)
                locPiMinusP4_Thrown = dThrownWrapper->Get_P4();
            else if (dThrownWrapper->Get_PID() == Proton)
                locProtonP4_Thrown  = dThrownWrapper->Get_P4();
        }
        locMissingP4_Thrown = locPiMinusP4_Thrown + locProtonP4_Thrown - locBeamP4_Thrown;

		// GET RF TIMING INFO
		TLorentzVector locBeamX4_Measured              = dComboBeamWrapper->Get_X4_Measured();
		Double_t       locBunchPeriod                  = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t       locDeltaT_RF                    = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t          locRelBeamBucket                = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t          locNumOutOfTimeBunchesInTree    = 4;                                                                                             // Number of out-of-time beam bunches in tree on a single side
		Bool_t         locSkipNearestOutOfTimeBunch    = true;                                                                                          // true: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t          locNumOutOfTimeBunchesToUse     = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		Double_t       locAccidentalScalingFactor      = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4_KinFit.E(), dIsMC);  // ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t       locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4_KinFit.E());    // ideal value would be 1, but deviations observed, need added factor.
		Double_t       locHistAccidWeightFactor        = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ;        // weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)                                                                                    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
        { 
		 	dComboWrapper->Set_IsComboCut(true); 
			continue;
        }

        // VECTORS in C.M. FRAME
        TVector3       boostCM_Measured        = (locPiMinusP4_Measured + locProtonP4_Measured).BoostVector();
        TVector3       boostCM_KinFit          = (locPiMinusP4_KinFit   + locProtonP4_KinFit).BoostVector();
        TLorentzVector locBeamP4CM_Measured    = locBeamP4_Measured;
        TLorentzVector locBeamP4CM_KinFit      = locBeamP4_KinFit;
        TLorentzVector locPiMinusP4CM_Measured = locPiMinusP4_Measured;
        TLorentzVector locPiMinusP4CM_KinFit   = locPiMinusP4_KinFit;
        locBeamP4CM_Measured.Boost(-boostCM_Measured);
        locBeamP4CM_KinFit.Boost(-boostCM_KinFit);
        locPiMinusP4CM_Measured.Boost(-boostCM_Measured); 
        locPiMinusP4CM_KinFit.Boost(-boostCM_KinFit);

        // CALCULATE CUSTOM VARIABLES
        double locVertexR_KinFit               = sqrt(pow(dComboBeamWrapper->Get_X4().X(),2) + pow(dComboBeamWrapper->Get_X4().Y(),2));
        double locChiSquarePerNDF_KinFit       = dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit();
        double locSqrtS_Measured               = (locPiMinusP4_Measured + locProtonP4_Measured).Mag();
        double locSqrtS_KinFit                 = (locPiMinusP4_KinFit   + locProtonP4_KinFit).Mag();
        double locMinusT_Measured              = -(locBeamP4_Measured - locPiMinusP4_Measured).Mag2();
        double locMinusT_KinFit                = -(locBeamP4_KinFit   - locPiMinusP4_KinFit).Mag2();
        double locMinusU_Measured              = -(locBeamP4_Measured - locProtonP4_Measured).Mag2();
        double locMinusU_KinFit                = -(locBeamP4_KinFit   - locProtonP4_KinFit).Mag2();
        double locThetaCM_Measured             = locBeamP4CM_Measured.Vect().Angle(locPiMinusP4CM_Measured.Vect())*RadToDeg;
        double locThetaCM_KinFit               = locBeamP4CM_KinFit.Vect().Angle(locPiMinusP4CM_KinFit.Vect())*RadToDeg;
        double locCoplanarity_Measured         = abs(locPiMinusP4_Measured.Phi() - locProtonP4_Measured.Phi())*RadToDeg;        
        double locCoplanarity_KinFit           = abs(locPiMinusP4_KinFit.Phi()   - locProtonP4_KinFit.Phi())*RadToDeg;
        
        // FILL CUSTOM HISTOGRAMS: BEFORE CUTS
        dHist_NumUnusedTracks_Measured_Before  ->Fill(dComboWrapper->Get_NumUnusedTracks(),                                         locHistAccidWeightFactor);
        dHist_NumUnusedShowers_Measured_Before ->Fill(dComboWrapper->Get_NumUnusedShowers(),                                        locHistAccidWeightFactor);
        dHist_PiMinusPIDFOM_Measured_Before    ->Fill(dPiMinusWrapper->Get_PIDFOM(),                                                locHistAccidWeightFactor);
        dHist_ProtonPIDFOM_Measured_Before     ->Fill(dProtonWrapper->Get_PIDFOM(),                                                 locHistAccidWeightFactor);
        dHist_PhotonEnergy_Measured_Before     ->Fill(locBeamP4_Measured.E(),                                                       locHistAccidWeightFactor);       
        dHist_ConfidenceLevel_KinFit_Before    ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(),                                  locHistAccidWeightFactor);
        dHist_VertexZ_KinFit_Before            ->Fill(dComboBeamWrapper->Get_X4().Z(),                                              locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_Before           ->Fill(dComboBeamWrapper->Get_X4().X(),             dComboBeamWrapper->Get_X4().Y(), locHistAccidWeightFactor);
        dHist_MissingPMinus_KinFit_Before      ->Fill(locMissingP4_KinFit.Minus(),                                                  locHistAccidWeightFactor);  
        dHist_MissingMomentum_KinFit_Before    ->Fill(locMissingP4_KinFit.P(),                                                      locHistAccidWeightFactor); 

        // INITIALIZE THE CUT FLAGS
        bool locCutFlag_NumUnusedTracks  = false;
        bool locCutFlag_NumUnusedShowers = false;
        bool locCutFlag_TrackPIDFOM      = false;
        bool locCutFlag_PhotonEnergy     = false;
        bool locCutFlag_ConfidenceLevel  = false;
        bool locCutFlag_CommonVertex     = false;
        bool locCutFlag_MissingMomentum  = false;
        
        // PERFORM THE CUTS
        if(dComboWrapper->Get_NumUnusedTracks()        > 5)                                                                         locCutFlag_NumUnusedTracks  = true;
        if(dComboWrapper->Get_NumUnusedShowers()       > 5)                                                                         locCutFlag_NumUnusedShowers = true;
        if(dPiMinusWrapper->Get_PIDFOM()               < 0.1  || dProtonWrapper->Get_PIDFOM()    < 0.1)                             locCutFlag_TrackPIDFOM      = true;
        if(locBeamP4_Measured.E()                      < 6.0  || locBeamP4_Measured.E()          > 10.8)                            locCutFlag_PhotonEnergy     = true;
        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < 0.01)                                                                      locCutFlag_ConfidenceLevel  = true;
        if(dComboBeamWrapper->Get_X4().Z()             < 51.0 || dComboBeamWrapper->Get_X4().Z() > 79.0 || locVertexR_KinFit > 1.0) locCutFlag_CommonVertex     = true;                                                 
        if(locMissingP4_KinFit.P()                     > 0.6)                                                                       locCutFlag_MissingMomentum  = true;
        
        // FILL CUSTOM HISTOGRAMS: CUT EFFECTS
        if(                               !locCutFlag_NumUnusedShowers && !locCutFlag_TrackPIDFOM && !locCutFlag_PhotonEnergy && !locCutFlag_ConfidenceLevel && !locCutFlag_CommonVertex && !locCutFlag_MissingMomentum)  dHist_CutEffect_NumUnusedTracks  ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);
        if(!locCutFlag_NumUnusedTracks &&                                 !locCutFlag_TrackPIDFOM && !locCutFlag_PhotonEnergy && !locCutFlag_ConfidenceLevel && !locCutFlag_CommonVertex && !locCutFlag_MissingMomentum)  dHist_CutEffect_NumUnusedShowers ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);
        if(!locCutFlag_NumUnusedTracks && !locCutFlag_NumUnusedShowers &&                            !locCutFlag_PhotonEnergy && !locCutFlag_ConfidenceLevel && !locCutFlag_CommonVertex && !locCutFlag_MissingMomentum)  dHist_CutEffect_TrackPIDFOM      ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);        
        if(!locCutFlag_NumUnusedTracks && !locCutFlag_NumUnusedShowers && !locCutFlag_TrackPIDFOM &&                             !locCutFlag_ConfidenceLevel && !locCutFlag_CommonVertex && !locCutFlag_MissingMomentum)  dHist_CutEffect_PhotonEnergy     ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);
        if(!locCutFlag_NumUnusedTracks && !locCutFlag_NumUnusedShowers && !locCutFlag_TrackPIDFOM && !locCutFlag_PhotonEnergy &&                                !locCutFlag_CommonVertex && !locCutFlag_MissingMomentum)  dHist_CutEffect_ConfidenceLevel  ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);
        if(!locCutFlag_NumUnusedTracks && !locCutFlag_NumUnusedShowers && !locCutFlag_TrackPIDFOM && !locCutFlag_PhotonEnergy && !locCutFlag_ConfidenceLevel &&                             !locCutFlag_MissingMomentum)  dHist_CutEffect_CommonVertex     ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);    
        if(!locCutFlag_NumUnusedTracks && !locCutFlag_NumUnusedShowers && !locCutFlag_TrackPIDFOM && !locCutFlag_PhotonEnergy && !locCutFlag_ConfidenceLevel && !locCutFlag_CommonVertex                               )  dHist_CutEffect_MissingMomentum  ->Fill(locMissingP4_KinFit.Minus(), locHistAccidWeightFactor);
      
        // SKIP THE REST OF THE LOOP IF CUT
        if( locCutFlag_NumUnusedTracks ||  locCutFlag_NumUnusedShowers ||  locCutFlag_TrackPIDFOM ||  locCutFlag_PhotonEnergy ||  locCutFlag_ConfidenceLevel ||  locCutFlag_CommonVertex ||  locCutFlag_MissingMomentum)
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }

        if (dComboWrapper->Get_L1TriggerBits() == 0)   
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }     
        
        // FILL CUSTOM HISTOGRAMS: AFTER CUTS
        dHist_NumUnusedTracks_Measured_After    ->Fill(dComboWrapper->Get_NumUnusedTracks(),                                                 locHistAccidWeightFactor);
        dHist_NumUnusedShowers_Measured_After   ->Fill(dComboWrapper->Get_NumUnusedShowers(),                                                locHistAccidWeightFactor);
        dHist_PiMinusPIDFOM_Measured_After      ->Fill(dPiMinusWrapper->Get_PIDFOM(),                                                        locHistAccidWeightFactor);
        dHist_ProtonPIDFOM_Measured_After       ->Fill(dProtonWrapper->Get_PIDFOM(),                                                         locHistAccidWeightFactor);
        dHist_PiMinusPVsdEdx_Measured_After     ->Fill(locPiMinusP4_Measured.P(),                   dPiMinusWrapper->Get_dEdx_CDC()*1000000, locHistAccidWeightFactor);
        dHist_ProtonPVsdEdx_Measured_After      ->Fill(locProtonP4_Measured.P(),                    dProtonWrapper->Get_dEdx_CDC()*1000000,  locHistAccidWeightFactor);
        dHist_PhotonEnergy_Measured_After       ->Fill(locBeamP4_Measured.E(),                                                               locHistAccidWeightFactor);       
        dHist_PhotonTiming_Measured_After       ->Fill(locDeltaT_RF,                                                                         locHistAccidWeightFactor);       
        dHist_ConfidenceLevel_KinFit_After      ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(),                                          locHistAccidWeightFactor);
        dHist_ChiSquarePerNDF_KinFit_After      ->Fill(locChiSquarePerNDF_KinFit,                                                            locHistAccidWeightFactor);
        dHist_VertexZ_KinFit_After              ->Fill(dComboBeamWrapper->Get_X4().Z(),                                                      locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_After             ->Fill(dComboBeamWrapper->Get_X4().X(),             dComboBeamWrapper->Get_X4().Y(),         locHistAccidWeightFactor);
        dHist_MissingPMinus_Measured_After      ->Fill(locMissingP4_Measured.Minus(),                                                        locHistAccidWeightFactor); 
        dHist_MissingPMinus_KinFit_After        ->Fill(locMissingP4_KinFit.Minus(),                                                          locHistAccidWeightFactor);  
        dHist_MissingMassSquared_Measured_After ->Fill(locMissingP4_Measured.M2(),                                                           locHistAccidWeightFactor); 
        dHist_MissingMassSquared_KinFit_After   ->Fill(locMissingP4_KinFit.M2(),                                                             locHistAccidWeightFactor);
        dHist_MissingMomentum_Measured_After    ->Fill(locMissingP4_Measured.P(),                                                            locHistAccidWeightFactor); 
        dHist_MissingMomentum_KinFit_After      ->Fill(locMissingP4_KinFit.P(),                                                              locHistAccidWeightFactor);
        dHist_SqrtS_Measured_After              ->Fill(locSqrtS_Measured,                                                                    locHistAccidWeightFactor);    
        dHist_SqrtS_KinFit_After                ->Fill(locSqrtS_KinFit,                                                                      locHistAccidWeightFactor);
        dHist_MinusT_Measured_After             ->Fill(locMinusT_Measured,                                                                   locHistAccidWeightFactor);    
        dHist_MinusT_KinFit_After               ->Fill(locMinusT_KinFit,                                                                     locHistAccidWeightFactor);    
        dHist_MinusU_Measured_After             ->Fill(locMinusU_Measured,                                                                   locHistAccidWeightFactor);      
        dHist_MinusU_KinFit_After               ->Fill(locMinusU_KinFit,                                                                     locHistAccidWeightFactor);    
        dHist_ThetaCM_Measured_After            ->Fill(locThetaCM_Measured,                                                                  locHistAccidWeightFactor);
        dHist_ThetaCM_KinFit_After              ->Fill(locThetaCM_KinFit,                                                                    locHistAccidWeightFactor);
        dHist_Coplanarity_Measured_After        ->Fill(locCoplanarity_Measured,                                                              locHistAccidWeightFactor);      
        dHist_Coplanarity_KinFit_After          ->Fill(locCoplanarity_KinFit,                                                                locHistAccidWeightFactor);
        dHist_PiMinusPVsTheta_Measured_After    ->Fill(locPiMinusP4_Measured.P(),                   locPiMinusP4_Measured.Theta()*RadToDeg,  locHistAccidWeightFactor);     
        dHist_PiMinusPVsTheta_KinFit_After      ->Fill(locPiMinusP4_KinFit.P(),                     locPiMinusP4_KinFit.Theta()*RadToDeg,    locHistAccidWeightFactor);       
        dHist_ProtonPVsTheta_Measured_After     ->Fill(locProtonP4_Measured.P(),                    locProtonP4_Measured.Theta()*RadToDeg,   locHistAccidWeightFactor);    
        dHist_ProtonPVsTheta_KinFit_After       ->Fill(locProtonP4_KinFit.P(),                      locProtonP4_KinFit.Theta()*RadToDeg,     locHistAccidWeightFactor);
        
		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        //FILL CUSTOM BRANCHES: FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Int_t>("RunNumber",             locRunNumber);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("Entry",                 locEntry);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("Combo",                 loc_i);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("L1TriggerBits",         dComboWrapper->Get_L1TriggerBits());
        dFlatTreeInterface->Fill_Fundamental<Double_t>("HistAccidWeightFactor", locHistAccidWeightFactor);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("BeamP4_Measured",     locBeamP4_Measured);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("PiMinusP4_Measured",  locPiMinusP4_Measured);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("ProtonP4_Measured",   locProtonP4_Measured);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("MissingP4_Measured",  locMissingP4_Measured);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("BeamP4_KinFit",       locBeamP4_KinFit);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("PiMinusP4_KinFit",    locPiMinusP4_KinFit);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("ProtonP4_KinFit",     locProtonP4_KinFit);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("MissingP4_KinFit",    locMissingP4_KinFit);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("BeamP4_Thrown",       locBeamP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("PiMinusP4_Thrown",    locPiMinusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("ProtonP4_Thrown",     locProtonP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("MissingP4_Thrown",    locMissingP4_Thrown);
        
        // FILL FLAT TREE
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


