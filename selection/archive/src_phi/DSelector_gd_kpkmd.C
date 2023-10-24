#include "DSelector_gd_kpkmd.h"

double RadToDeg = 180.0/3.1415926;

void DSelector_gd_kpkmd::Init(TTree *locTree)
{    
	//SET OUTPUT FILE NAME
	dOutputFileName     = "hist_gd_kpkmd.root"; //"" for none
	dOutputTreeFileName = "";                       //"" for none
	dFlatTreeFileName   = "";                       //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName       = "";                       //if blank, default name will be chosen

	//INITIALIZE THE TREE INTERFACE AND WRAPPERS
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree);                        //This must be called to initialize wrappers for each new TTree, gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return;                                      //have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    Initialize_Actions();
    
    //MC INFORMATION
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);
    
    //CUSTOM HISTOGRAMS: BEFORE THE CUTS
    dHist_NumUnusedTracks_Measured_Before         = new TH1F("NumUnusedTracks_Measured_Before",        ";Unused Tracks                  ;Events/1",             10,  0.0,  10.0);
    dHist_NumUnusedShowers_Measured_Before        = new TH1F("NumUnusedShowers_Measured_Before",       ";Unused Showers                 ;Events/1",             10,  0.0,  10.0);
    dHist_ConfidenceLevel_KinFit_Before           = new TH1F("ConfidenceLevel_KinFit_Before",          ";Confidence Level               ;Events/0.001",       1000,  0.0,   1.0);
    dHist_InvariantMass_Measured_Before           = new TH1F("InvariantMass_Measured_Before",          ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.8,   2.8);
    dHist_MissingPMinus_Measured_Before           = new TH1F("MissingPMinus_Measured_Before",          ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     400,  0.0,   4.0);
    dHist_MissingMassSquared_Measured_Before      = new TH1F("MissingMassSquared_Measured_Before",     ";Missing Mass Squared (GeV^{2}) ;Events/0.01 GeV^{2}", 600,  0.0,   6.0);
    dHist_MissingMomentum_Measured_Before         = new TH1F("MissingMomentum_Measured_Before",        ";Missing Momentum (GeV)         ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_Coplanarity_Measured_Before             = new TH1F("Coplanarity_Measured_Before",            ";Coplanarity Angle (deg)        ;Events/1 deg",        360,  0.0, 360.0);
    dHist_PhotonEnergy_Measured_Before            = new TH1F("PhotonEnergy_Before",                    ";Photon Energy (GeV)            ;Events/0.01 GeV",     900,  3.0,  12.0);
    dHist_VertexZ_KinFit_Before                   = new TH1F("VertexZ_KinFit_Before",                  ";Vertex Z (cm)                  ;Events/1 cm",         200,  0.0, 200.0);
    dHist_VertexXY_KinFit_Before                  = new TH2F("VertexXY_KinFit_Before",                 ";Vertex X (cm)                  ;Vertex Y (cm)",       100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_MinusT_Measured_Before                  = new TH1F("MinusT_Measured_Before",                 ";-t (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_MinusU_Measured_Before                  = new TH1F("MinusU_Measured_Before",                 ";-u (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_KPlusPVsdEdx_Measured_Before            = new TH2F("KPlusPVsdEdx_Measured_Before",           ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);
    dHist_KMinusPVsdEdx_Measured_Before           = new TH2F("KMinusPVsdEdx_Measured_Before",          ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);    
    dHist_DeuteronPVsdEdx_Measured_Before         = new TH2F("DeuteronPVsdEdx_Measured_Before",        ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);
    
    //CUSTOM HISTOGRAMS: CUT EFFECTS
    dHist_CutEffect_NumUnusedTracks               = new TH1F("CutEffect_NumUnusedTracks",              ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_NumUnusedShowers              = new TH1F("CutEffect_NumUnusedShowers",             ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_ConfidenceLevel               = new TH1F("CutEffect_ConfidenceLevel",              ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_MissingPMinus                 = new TH1F("CutEffect_MissingPMinus",                ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_Coplanarity                   = new TH1F("CutEffect_Coplanarity",                  ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_CommonVertex                  = new TH1F("CutEffect_CommonVertex",                 ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutEffect_TAndU                         = new TH1F("CutEffect_TAndU",                        ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    
    //FILL CUSTOM HISTOGRAMS: EXPERIMENTING WITH CUTS
    dHist_CutTest_NoShower                        = new TH1F("CutTest_NoShower",                       ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutTest_OneShower                       = new TH1F("CutTest_OneShower",                      ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutTest_TwoShower                       = new TH1F("CutTest_TwoShower",                      ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutTest_ThreeShower                     = new TH1F("CutTest_ThreeShower",                    ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_CutTest_AnyShower                       = new TH1F("CutTest_AnyShower",                      ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.0,   2.0);
    
    //CUSTOM HISTOGRAMS: AFTER THE CUTS
    dHist_ConfidenceLevel_KinFit_After            = new TH1F("ConfidenceLevel_KinFit_After",           ";Confidence Level               ;Events/0.001",       1000,  0.0,   1.0);
    dHist_ChiSquarePerNDF_KinFit_After            = new TH1F("ChiSquarePerNDF_KinFit_After",           ";#chi^{2}/NDF;                  ;Events/0.01",        1000,  0.0,  10.0);
    dHist_InvariantMass_Measured_After            = new TH1F("InvariantMass_Measured_After",           ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.8,   2.8);
    dHist_InvariantMass_KinFit_After              = new TH1F("InvariantMass_KinFit_After",             ";Invariant Mass (K^+K^-)        ;Events/0.01 GeV",     200,  0.8,   2.8);
    dHist_MissingPMinus_Measured_After            = new TH1F("MissingPMinus_Measured_After",           ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     400,  0.0,   4.0);
    dHist_MissingPMinus_KinFit_After              = new TH1F("MissingPMinus_KinFit_After",             ";Missing P^{-} (GeV)            ;Events/0.01 GeV",     400,  0.0,   4.0);        
    dHist_MissingMassSquared_Measured_After       = new TH1F("MissingMassSquared_Measured_After",      ";Missing Mass Squared (GeV^{2}) ;Events/0.01 GeV^{2}", 600,  0.0,   6.0);
    dHist_MissingMassSquared_KinFit_After         = new TH1F("MissingMassSquared_KinFit_After",        ";Missing Mass Squared (GeV^{2}) ;Events/0.01 GeV^{2}", 600,  0.0,   6.0);
    dHist_MissingMomentum_Measured_After          = new TH1F("MissingMomentum_Measured_After",         ";Missing Momentum (GeV)         ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_MissingMomentum_KinFit_After            = new TH1F("MissingMomentum_KinFit_After",           ";Missing Momentum (GeV)         ;Events/0.01 GeV",     200,  0.0,   2.0);
    dHist_SqrtS_Measured_After                    = new TH1F("SqrtS_Measured_After",                   ";#sqrt(s) (GeV)                 ;Events/0.01 GeV",    1000,  0.0,  10.0);
    dHist_SqrtS_KinFit_After                      = new TH1F("SqrtS_KinFit_After",                     ";#sqrt(s) (GeV)                 ;Events/0.01 GeV",    1000,  0.0,  10.0);
    dHist_MinusT_Measured_After                   = new TH1F("MinusT_Measured_After",                  ";-t (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_MinusT_KinFit_After                     = new TH1F("MinusT_KinFit_After",                    ";-t (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_MinusU_Measured_After                   = new TH1F("MinusU_Measured_After",                  ";-u (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_MinusU_KinFit_After                     = new TH1F("MinusU_KinFit_After",                    ";-u (GeV^{2})                   ;Events/0.01 GeV^{2}", 500,  0.0,   5.0);
    dHist_Coplanarity_Measured_After              = new TH1F("Coplanarity_Measured_After",             ";Coplanarity Angle (deg)        ;Events/1 deg",        360,  0.0, 360.0);
    dHist_Coplanarity_KinFit_After                = new TH1F("Coplanarity_KinFit_After",               ";Coplanarity Angle (deg)        ;Events/1 deg",        360,  0.0, 360.0);
    dHist_ThetaCM_Measured_After                  = new TH1F("ThetaCM_Measured_After",                 ";#theta_{c.m.} (deg)            ;Events/1 deg",        180,  0.0, 180.0);
    dHist_ThetaCM_KinFit_After                    = new TH1F("ThetaCM_KinFit_After",                   ";#theta_{c.m.} (deg)            ;Events/1 deg",        180,  0.0, 180.0);
    dHist_PhotonEnergy_Measured_After             = new TH1F("PhotonEnergy_After",                     ";Photon Energy (GeV)            ;Events/0.01 GeV",     900,  3.0,  12.0);
    dHist_PhotonTiming_Measured_After             = new TH1F("PhotonTiming_After",                     ";#Delta t_{Beam-RF} (ns)        ;Events/0.1 ns",       360,-18.0,  18.0);
    dHist_VertexZ_KinFit_After                    = new TH1F("VertexZ_KinFit_After",                   ";Vertex Z (cm)                  ;Events/1 cm",         200,  0.0, 200.0);
    dHist_VertexXY_KinFit_After                   = new TH2F("VertexXY_KinFit_After",                  ";Vertex X (cm)                  ;Vertex Y (cm)",       100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_KPlusPVsdEdx_Measured_After             = new TH2F("KPlusPVsdEdx_Measured_After",            ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);
    dHist_KMinusPVsdEdx_Measured_After            = new TH2F("KMinusPVsdEdx_Measured_After",           ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);    
    dHist_DeuteronPVsdEdx_Measured_After          = new TH2F("DeuteronPVsdEdx_Measured_After",         ";P (GeV)                        ;dE/dx (keV/cm)",      100,  0.0,  10.0,  250,  0.0,  25.0);
    dHist_KPlusPVsTheta_Measured_After            = new TH2F("KPlusPVsTheta_Measured_After",           ";P (GeV)                        ;#theta (deg)",       1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_KPlusPVsTheta_KinFit_After              = new TH2F("KPlusPVsTheta_KinFit_After",             ";P (GeV)                        ;#theta (deg)",       1000,  0.0, 10.0,   180,  0.0, 180.0);
    dHist_KMinusPVsTheta_Measured_After           = new TH2F("KMinusPVsTheta_Measured_After",          ";P (GeV)                        ;#theta (deg)",       1000,  0.0, 10.0,   180,  0.0, 180.0);
    dHist_KMinusPVsTheta_KinFit_After             = new TH2F("KMinusPVsTheta_KinFit_After",            ";P (GeV)                        ;#theta (deg)",       1000,  0.0, 10.0,   180,  0.0, 180.0);
    dHist_DeuteronPVsTheta_Measured_After         = new TH2F("DeuteronPVsTheta_Measured_After",        ";P (GeV)                        ;#theta (deg)",       1000,  0.0, 10.0,   180,  0.0, 180.0);
    dHist_DeuteronPVsTheta_KinFit_After           = new TH2F("DeuteronPVsTheta_KinFit_After",          ";P (GeV)                        ;#theta (deg)",       1000,  0.0, 10.0,   180,  0.0, 180.0);

} 
//END OF INITIALIZATION

Bool_t DSelector_gd_kpkmd::Process(Long64_t locEntry)
{    
	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
    Reset_Actions_NewEvent();
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

    //GET PHOTON POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	//LOOP OVER COMBOS
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
    {
		//INITIALIZE THE COMBO
		dComboWrapper->Set_ComboIndex(loc_i);  //Set branch array indices
		if(dComboWrapper->Get_IsComboCut()) //check whether the combo has been cut
			continue; // Combo has been cut previously

        //GET PARTICLE INDICES
		Int_t locBeamID          = dComboBeamWrapper->Get_BeamID();
		Int_t locKPlusTrackID    = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID   = dKMinusWrapper->Get_TrackID();
        Int_t locDeuteronTrackID = dDeuteronWrapper->Get_TrackID();
        
        //GET KINFIT P4
		TLorentzVector locBeamP4                     = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4                    = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4                   = dKMinusWrapper->Get_P4();
        TLorentzVector locDeuteronP4                 = dDeuteronWrapper->Get_P4();
        
        // Get MEASURED P4
		TLorentzVector locBeamP4_Measured            = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured           = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured          = dKMinusWrapper->Get_P4_Measured();
        TLorentzVector locDeuteronP4_Measured        = dDeuteronWrapper->Get_P4_Measured();

		//GET RF TIMING INFO
		TLorentzVector locBeamX4_Measured              = dComboBeamWrapper->Get_X4_Measured();
		Double_t       locBunchPeriod                  = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t       locDeltaT_RF                    = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t          locRelBeamBucket                = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t          locNumOutOfTimeBunchesInTree    = 4; //YOU need to specify this number. Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 
		Bool_t         locSkipNearestOutOfTimeBunch    = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t          locNumOutOfTimeBunchesToUse     = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		Double_t       locAccidentalScalingFactor      = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t       locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		Double_t       locHistAccidWeightFactor        = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)  // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
        { 
		 	dComboWrapper->Set_IsComboCut(true); 
			continue;
        }
        
        //VECTORS in C.M. FRAME
        TVector3       boostCM_Measured        = (locBeamP4_Measured + dTargetP4).BoostVector();
        TVector3       boostCM_KinFit          = (locBeamP4 + dTargetP4).BoostVector();
        TLorentzVector locBeamP4CM_Measured    = locBeamP4_Measured;
        TLorentzVector locBeamP4CM_KinFit      = locBeamP4;
        TLorentzVector locKPlusP4CM_Measured   = locKPlusP4_Measured;
        TLorentzVector locKPlusP4CM_KinFit     = locKPlusP4;
        TLorentzVector locKMinusP4CM_Measured  = locKMinusP4_Measured;
        TLorentzVector locKMinusP4CM_KinFit    = locKMinusP4;        
        locBeamP4CM_Measured.Boost(-boostCM_Measured);
        locBeamP4CM_KinFit.Boost(-boostCM_KinFit);
        locKPlusP4CM_Measured.Boost(-boostCM_Measured); 
        locKPlusP4CM_KinFit.Boost(-boostCM_KinFit);
        locKMinusP4CM_Measured.Boost(-boostCM_Measured); 
        locKMinusP4CM_KinFit.Boost(-boostCM_KinFit);
        
        //CALCULATE CUSTOM VARIABLES
        double locConfidenceLevel_KinFit       = dComboWrapper->Get_ConfidenceLevel_KinFit();
        double locChiSquarePerNDF_KinFit       = dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit();
        double locVertexR_Measured             = sqrt(pow(locBeamX4_Measured.X(),2) + pow(locBeamX4_Measured.Y(),2));
        double locInvariantMass_Measured       = (locKPlusP4_Measured + locKMinusP4_Measured).M();
        double locInvariantMass_KinFit         = (locKPlusP4 + locKMinusP4).M();
        double locMissingPMinus_Measured       = (locBeamP4_Measured + dTargetP4 - locKPlusP4_Measured - locKMinusP4_Measured).Minus();
        double locMissingPMinus_KinFit         = (locBeamP4 + dTargetP4 - locKPlusP4 - locKMinusP4).Minus(); 
        double locMissingMassSquared_Measured  = (locBeamP4_Measured + dTargetP4 - locKPlusP4_Measured - locKMinusP4_Measured).M2();
        double locMissingMassSquared_KinFit    = (locBeamP4 + dTargetP4 - locKPlusP4 - locKMinusP4).M2();
        double locMissingMomentum_Measured     = (locBeamP4_Measured + dTargetP4 - locKPlusP4_Measured - locKMinusP4_Measured).P();
        double locMissingMomentum_KinFit       = (locBeamP4 + dTargetP4 - locKPlusP4 - locKMinusP4).P();
        double locSqrtS_Measured               = (locBeamP4_Measured + dTargetP4).Mag();
        double locSqrtS_KinFit                 = (locBeamP4 + dTargetP4).Mag();
        double locMinusT_Measured              = -(locBeamP4_Measured - locKPlusP4_Measured - locKMinusP4_Measured).Mag2();
        double locMinusT_KinFit                = -(locBeamP4 - locKPlusP4 - locKMinusP4).Mag2();
        double locMinusU_Measured              = -(locBeamP4_Measured - locDeuteronP4_Measured).Mag2();
        double locMinusU_KinFit                = -(locBeamP4 - locDeuteronP4).Mag2();
        double locCoplanarity_Measured         = abs((locKPlusP4_Measured + locKMinusP4_Measured).Phi() - locDeuteronP4_Measured.Phi())*RadToDeg;        
        double locCoplanarity_KinFit           = abs((locKPlusP4 + locKMinusP4).Phi() - locDeuteronP4.Phi())*RadToDeg;
        double locThetaCM_Measured             = locBeamP4CM_Measured.Vect().Angle((locKPlusP4CM_Measured + locKMinusP4CM_Measured).Vect())*RadToDeg;
        double locThetaCM_KinFit               = locBeamP4CM_KinFit.Vect().Angle((locKPlusP4CM_KinFit + locKMinusP4CM_KinFit).Vect())*RadToDeg;
        
        //FILL CUSTOM HISTOGRAMS: BEFORE CUTS
        dHist_NumUnusedTracks_Measured_Before     ->Fill(dComboWrapper->Get_NumUnusedTracks(),                                     locHistAccidWeightFactor);
        dHist_NumUnusedShowers_Measured_Before    ->Fill(dComboWrapper->Get_NumUnusedShowers(),                                    locHistAccidWeightFactor);
        dHist_ConfidenceLevel_KinFit_Before       ->Fill(locConfidenceLevel_KinFit,                                                locHistAccidWeightFactor);
        dHist_InvariantMass_Measured_Before       ->Fill(locInvariantMass_Measured,                                                locHistAccidWeightFactor);
        dHist_MissingPMinus_Measured_Before       ->Fill(locMissingPMinus_Measured,                                                locHistAccidWeightFactor);
        dHist_MissingMassSquared_Measured_Before  ->Fill(locMissingMassSquared_Measured,                                           locHistAccidWeightFactor);
        dHist_MissingMomentum_Measured_Before     ->Fill(locMissingMomentum_Measured,                                              locHistAccidWeightFactor); 
        dHist_Coplanarity_Measured_Before         ->Fill(locCoplanarity_Measured,                                                  locHistAccidWeightFactor);      
        dHist_PhotonEnergy_Measured_Before        ->Fill(locBeamP4_Measured.E(),                                                   locHistAccidWeightFactor);       
        dHist_VertexZ_KinFit_Before               ->Fill(locBeamX4_Measured.Z(),                                                   locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_Before              ->Fill(locBeamX4_Measured.X(),         locBeamX4_Measured.Y(),                   locHistAccidWeightFactor);
        dHist_MinusT_Measured_Before              ->Fill(locMinusT_Measured,                                                       locHistAccidWeightFactor);    
        dHist_MinusU_Measured_Before              ->Fill(locMinusU_Measured,                                                       locHistAccidWeightFactor);      
        dHist_KPlusPVsdEdx_Measured_Before        ->Fill(locKPlusP4_Measured.P(),        dKPlusWrapper->Get_dEdx_CDC()*1000000,    locHistAccidWeightFactor);
        dHist_KMinusPVsdEdx_Measured_Before       ->Fill(locKMinusP4_Measured.P(),       dKMinusWrapper->Get_dEdx_CDC()*1000000,   locHistAccidWeightFactor);
        dHist_DeuteronPVsdEdx_Measured_Before     ->Fill(locDeuteronP4_Measured.P(),     dDeuteronWrapper->Get_dEdx_CDC()*1000000, locHistAccidWeightFactor);
        
        //EXECUTE MANUAL CUTS
        //Confidence level cut
        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < 0.01){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        //if((locMissingMassSquared_Measured < 2.5) || (locMissingMassSquared_Measured > 4.5)){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}
        //if((locCoplanarityAngle_Measured < 170.0) || (locCoplanarityAngle_Measured > 190.0)){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}
        if(dComboWrapper->Get_NumUnusedTracks() > 0){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        if(dComboWrapper->Get_NumUnusedShowers() > 0){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        //if((locBeamX4_Measured.Z() < 51.0) || (locBeamX4_Measured.Z() > 76.0) || sqrt(pow(locBeamX4_Measured.X(),2) + pow(locBeamX4_Measured.Y(),2)) > 1.0){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}    
        //if((locBeamP4_Measured.E() < 6.0) || (locBeamP4_Measured.E() > 11.0)){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}   
        //if((locInvariantMass_Measured < 0.9) || (locInvariantMass_Measured > 1.1)){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}
        
        //FILL CUSTOM HISTOGRAMS: AFTER CUTS
        dHist_ConfidenceLevel_KinFit_After           ->Fill(locConfidenceLevel_KinFit,                                                         locHistAccidWeightFactor);
        dHist_ChiSquarePerNDF_KinFit_After           ->Fill(locChiSquarePerNDF_KinFit,                                                         locHistAccidWeightFactor);
        dHist_InvariantMass_Measured_After           ->Fill(locInvariantMass_Measured,                                                         locHistAccidWeightFactor);
        dHist_InvariantMass_KinFit_After             ->Fill(locInvariantMass_KinFit,                                                           locHistAccidWeightFactor);
        dHist_MissingPMinus_Measured_After           ->Fill(locMissingPMinus_Measured,                                                         locHistAccidWeightFactor); 
        dHist_MissingPMinus_KinFit_After             ->Fill(locMissingPMinus_KinFit,                                                           locHistAccidWeightFactor);  
        dHist_MissingMassSquared_Measured_After      ->Fill(locMissingMassSquared_Measured,                                                    locHistAccidWeightFactor); 
        dHist_MissingMassSquared_KinFit_After        ->Fill(locMissingMassSquared_KinFit,                                                      locHistAccidWeightFactor);
        dHist_MissingMomentum_Measured_After         ->Fill(locMissingMomentum_Measured,                                                       locHistAccidWeightFactor); 
        dHist_MissingMomentum_KinFit_After           ->Fill(locMissingMomentum_KinFit,                                                         locHistAccidWeightFactor);
        dHist_SqrtS_Measured_After                   ->Fill(locSqrtS_Measured,                                                                 locHistAccidWeightFactor);    
        dHist_SqrtS_KinFit_After                     ->Fill(locSqrtS_KinFit,                                                                   locHistAccidWeightFactor);
        dHist_MinusT_Measured_After                  ->Fill(locMinusT_Measured,                                                                locHistAccidWeightFactor);    
        dHist_MinusT_KinFit_After                    ->Fill(locMinusT_KinFit,                                                                  locHistAccidWeightFactor);    
        dHist_MinusU_Measured_After                  ->Fill(locMinusU_Measured,                                                                locHistAccidWeightFactor);      
        dHist_MinusU_KinFit_After                    ->Fill(locMinusU_KinFit,                                                                  locHistAccidWeightFactor);    
        dHist_Coplanarity_Measured_After             ->Fill(locCoplanarity_Measured,                                                           locHistAccidWeightFactor);      
        dHist_Coplanarity_KinFit_After               ->Fill(locCoplanarity_KinFit,                                                             locHistAccidWeightFactor);
        dHist_ThetaCM_Measured_After                 ->Fill(locThetaCM_Measured,                                                               locHistAccidWeightFactor);
        dHist_ThetaCM_KinFit_After                   ->Fill(locThetaCM_KinFit,                                                                 locHistAccidWeightFactor);
        dHist_PhotonEnergy_Measured_After            ->Fill(locBeamP4_Measured.E(),                                                            locHistAccidWeightFactor);       
        dHist_PhotonTiming_Measured_After            ->Fill(locDeltaT_RF,                                                                      locHistAccidWeightFactor);       
        dHist_VertexZ_KinFit_After                   ->Fill(locBeamX4_Measured.Z(),                                                            locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_After                  ->Fill(locBeamX4_Measured.X(),            locBeamX4_Measured.Y(),                         locHistAccidWeightFactor);
        dHist_KPlusPVsdEdx_Measured_After            ->Fill(locKPlusP4_Measured.P(),           dKPlusWrapper->Get_dEdx_CDC()*1000000,          locHistAccidWeightFactor);
        dHist_KMinusPVsdEdx_Measured_After           ->Fill(locKMinusP4_Measured.P(),          dKMinusWrapper->Get_dEdx_CDC()*1000000,         locHistAccidWeightFactor);
        dHist_DeuteronPVsdEdx_Measured_After         ->Fill(locDeuteronP4_Measured.P(),        dDeuteronWrapper->Get_dEdx_CDC()*1000000,       locHistAccidWeightFactor);
        dHist_KPlusPVsTheta_Measured_After           ->Fill(locKPlusP4_Measured.P(),           locKPlusP4_Measured.Theta()*RadToDeg,           locHistAccidWeightFactor);     
        dHist_KPlusPVsTheta_KinFit_After             ->Fill(locKPlusP4.P(),                    locKPlusP4.Theta()*RadToDeg,                    locHistAccidWeightFactor);       
        dHist_KMinusPVsTheta_Measured_After          ->Fill(locKMinusP4_Measured.P(),          locKMinusP4_Measured.Theta()*RadToDeg,          locHistAccidWeightFactor);    
        dHist_KMinusPVsTheta_KinFit_After            ->Fill(locKMinusP4.P(),                   locKMinusP4.Theta()*RadToDeg,                   locHistAccidWeightFactor);
        dHist_DeuteronPVsTheta_Measured_After        ->Fill(locDeuteronP4_Measured.P(),        locDeuteronP4_Measured.Theta()*RadToDeg,        locHistAccidWeightFactor);    
        dHist_DeuteronPVsTheta_KinFit_After          ->Fill(locDeuteronP4.P(),                 locDeuteronP4.Theta()*RadToDeg,                 locHistAccidWeightFactor);

		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;
    }  
    //END OF COMBO LOOP

	return kTRUE;
}  
//end of processing

void DSelector_gd_kpkmd::Finalize(void)
{
	//CALL THIS LAST
	DSelector::Finalize();
} //end of finalization