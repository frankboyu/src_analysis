#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"

double RadToDeg = 180.0/3.1415926;
double mass_neutron = 0.939565421;

class DSelector_piminus_p_12C_sim : public DSelector
{
	public:
    
		DSelector_piminus_p_12C_sim(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_piminus_p_12C_sim(){}

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

		// CUSTOM HISTOGRAMS: BEFORE THE CUTS
        TH1F* dHist_NumUnusedTracks_Before;
        TH1F* dHist_NumUnusedShowers_Before;
        TH1F* dHist_PiMinusPIDFOM_Before;
        TH1F* dHist_ProtonPIDFOM_Before;
        TH2F* dHist_PiMinusPVsTheta_Before;
        TH2F* dHist_ProtonPVsTheta_Before;
        TH1F* dHist_ConfidenceLevel_Before;
        TH1F* dHist_MissingMomentum_Before;
        TH1F* dHist_EnergyBalance_Before;
        TH1F* dHist_PhotonEnergy_Before;
        TH1F* dHist_VertexZ_Before;
        TH2F* dHist_VertexXY_Before;
        TH1F* dHist_MissingPMinus_Before;
        TH1F* dHist_InvariantMassRho_Before;
        
        // CUSTOM HISTOGRAMS: CUT EFFECTS
        TH1F* dHist_CutEffect_NumUnusedTracks;
        TH1F* dHist_CutEffect_NumUnusedShowers;
        TH1F* dHist_CutEffect_TrackPIDFOM;
        TH1F* dHist_CutEffect_ParticleKinematics;
        TH1F* dHist_CutEffect_ConfidenceLevel;
        TH1F* dHist_CutEffect_MissingMomentum;
        TH1F* dHist_CutEffect_PhotonEnergy;
        TH1F* dHist_CutEffect_CommonVertex;
        TH1F* dHist_CutEffect_MissingPMinus;
        
        // CUSTOM HISTOGRAMS: AFTER THE CUTS
        TH1F* dHist_NumUnusedTracks_After;
        TH1F* dHist_NumUnusedShowers_After;
        TH1F* dHist_PiMinusPIDFOM_After;
        TH1F* dHist_ProtonPIDFOM_After;
        TH2F* dHist_PiMinusPVsTheta_After;
        TH2F* dHist_ProtonPVsTheta_After;
        TH1F* dHist_ConfidenceLevel_After;
        TH1F* dHist_MissingMomentum_After;
        TH1F* dHist_PhotonEnergy_After;
        TH1F* dHist_VertexZ_After;
        TH2F* dHist_VertexXY_After;
        TH1F* dHist_MissingPMinus_After;
        TH1F* dHist_PhotonTiming_After;      
        TH1F* dHist_ChiSquarePerNDF_After;
        TH1F* dHist_MissingMassSquared_After;
        TH1F* dHist_SqrtS_After;
        TH1F* dHist_MinusT_After;
        TH1F* dHist_MinusU_After;
        TH1F* dHist_ThetaCM_After;
        TH1F* dHist_Coplanarity_After;
        TH1F* dHist_InvariantMassRho_After;


	ClassDef(DSelector_piminus_p_12C_sim, 0);
};

void DSelector_piminus_p_12C_sim::Get_ComboWrappers(void)
{
	dStep0Wrapper     = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiMinusWrapper   = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dProtonWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
}

void DSelector_piminus_p_12C_sim::Init(TTree *locTree)
{
	// SET OUTPUT FILE NAME
	dOutputFileName          = "";                                 
	dOutputTreeFileName      = "";                                 
	dFlatTreeFileName        = "flattree_piminus_p_12C_sim.root";
	dFlatTreeName            = "flattree_piminus_p_12C_sim";              
    dSaveDefaultFlatBranches = false;                              

	// INITIALIZE THE TREE INTERFACE AND WRAPPERS
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    Initialize_Actions();

    // CUSTOM HISTOGRAMS: BEFORE THE CUTS
    dHist_NumUnusedTracks_Before       = new TH1F("NumUnusedTracks_Before",         ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_Before      = new TH1F("NumUnusedShowers_Before",        ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PiMinusPIDFOM_Before         = new TH1F("PiMinusPIDFOM_Before",           ";PIDFOM_{#pi^{-}}            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_Before          = new TH1F("ProtonPIDFOM_Before",            ";PIDFOM_{p}                  ;Events/0.001",          1000,  0.0,   1.0);
    dHist_PiMinusPVsTheta_Before       = new TH2F("PiMinusPVsTheta_Before",         ";P_{#pi^{-}} (GeV)           ;#theta (deg)",          1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ProtonPVsTheta_Before        = new TH2F("ProtonPVsTheta_Before",          ";P_{p} (GeV)                 ;#theta (deg)",          1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ConfidenceLevel_Before       = new TH1F("ConfidenceLevel_Before",         ";Confidence Level            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_MissingMomentum_Before       = new TH1F("MissingMomentum_Before",         ";P_{miss} (GeV)              ;Events/0.01 GeV",       1000,  0.0,  10.0);    
    dHist_EnergyBalance_Before         = new TH1F("EnergyBalance_Before",           ";Missing Energy (GeV)        ;Events/0.01 GeV",        800, -4.0,   4.0);
    dHist_PhotonEnergy_Before          = new TH1F("PhotonEnergy_Before",            ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_Before               = new TH1F("VertexZ_Before",                 ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_Before              = new TH2F("VertexXY_Before",                ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_MissingPMinus_Before         = new TH1F("MissingPMinus_Before",           ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);        
    dHist_InvariantMassRho_Before      = new TH1F("InvariantMassRho_Before",        ";M_{'#pi^{+}'#pi^{-}} (GeV)  ;Events/0.01 GeV",        500,  0.0,   5.0);

    // CUSTOM HISTOGRAMS: CUT EFFECTS
    dHist_CutEffect_NumUnusedTracks    = new TH1F("CutEffect_NumUnusedTracks",      ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_NumUnusedShowers   = new TH1F("CutEffect_NumUnusedShowers",     ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_TrackPIDFOM        = new TH1F("CutEffect_TrackPIDFOM",          ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_ParticleKinematics = new TH1F("CutEffect_ParticleKinematics",   ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_ConfidenceLevel    = new TH1F("CutEffect_ConfidenceLevel",      ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_MissingMomentum    = new TH1F("CutEffect_MissingMomentum",      ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_PhotonEnergy       = new TH1F("CutEffect_PhotonEnergy",         ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_CommonVertex       = new TH1F("CutEffect_CommonVertex",         ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    dHist_CutEffect_MissingPMinus      = new TH1F("CutEffect_MissingPMinus",        ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);
    
    // CUSTOM HISTOGRAMS: AFTER THE CUTS
    dHist_NumUnusedTracks_After        = new TH1F("NumUnusedTracks_After",     ";Unused Tracks               ;Events/1",                10,  0.0,  10.0);
    dHist_NumUnusedShowers_After       = new TH1F("NumUnusedShowers_After",    ";Unused Showers              ;Events/1",                10,  0.0,  10.0);
    dHist_PiMinusPIDFOM_After          = new TH1F("PiMinusPIDFOM_After",       ";PIDFOM_{#pi^{-}}            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_ProtonPIDFOM_After           = new TH1F("ProtonPIDFOM_After",        ";PIDFOM_{p}                  ;Events/0.001",          1000,  0.0,   1.0);
    dHist_PiMinusPVsTheta_After        = new TH2F("PiMinusPVsTheta_After",     ";P_{#pi^{-}} (GeV)           ;#theta (deg)",          1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ProtonPVsTheta_After         = new TH2F("ProtonPVsTheta_After",      ";P_{p} (GeV)                 ;#theta (deg)",          1000,  0.0,  10.0,  180,  0.0, 180.0);
    dHist_ConfidenceLevel_After        = new TH1F("ConfidenceLevel_After",     ";Confidence Level            ;Events/0.001",          1000,  0.0,   1.0);
    dHist_MissingMomentum_After        = new TH1F("MissingMomentum_After",     ";P_{miss} (GeV)              ;Events/0.01 GeV",       1000,  0.0,  10.0);    
    dHist_PhotonEnergy_After           = new TH1F("PhotonEnergy_After",        ";Photon Energy (GeV)         ;Events/0.01 GeV",        900,  3.0,  12.0);
    dHist_VertexZ_After                = new TH1F("VertexZ_After",             ";Vertex Z (cm)               ;Events/1 cm",            200,  0.0, 200.0);
    dHist_VertexXY_After               = new TH2F("VertexXY_After",            ";Vertex X (cm)               ;Vertex Y (cm)",          100, -5.0,   5.0,  100, -5.0,   5.0);
    dHist_MissingPMinus_After          = new TH1F("MissingPMinus_After",       ";P^{-}_{miss} (GeV)          ;Events/0.01 GeV",        200,  0.0,   2.0);        
    dHist_PhotonTiming_After           = new TH1F("PhotonTiming_After",        ";#Delta t_{Beam-RF} (ns)     ;Events/0.1 ns",          360,-18.0,  18.0);
    dHist_ChiSquarePerNDF_After        = new TH1F("ChiSquarePerNDF_After",     ";#chi^{2}/NDF;               ;Events/0.01",           1000,  0.0,  10.0);
    dHist_MissingMassSquared_After     = new TH1F("MissingMassSquared_After",  ";MM^{2}_{miss} (GeV^{2})     ;Events/0.01 GeV^{2}",    600,  0.0,   6.0);
    dHist_SqrtS_After                  = new TH1F("SqrtS_After",               ";#sqrt(s) (GeV)              ;Events/0.01 GeV",       1000,  0.0,  10.0);
    dHist_MinusT_After                 = new TH1F("MinusT_After",              ";-t (GeV^{2})                ;Events/0.01 GeV^{2}",   1500,  0.0,  15.0);
    dHist_MinusU_After                 = new TH1F("MinusU_After",              ";-u (GeV^{2})                ;Events/0.01 GeV^{2}",   1500,  0.0,  15.0);
    dHist_ThetaCM_After                = new TH1F("ThetaCM_After",             ";#theta_{c.m.} (deg)         ;Events/1 deg",           180,  0.0, 180.0);
    dHist_Coplanarity_After            = new TH1F("Coplanarity_After",         ";Coplanarity Angle (deg)     ;Events/1 deg",           360,  0.0, 360.0);
    dHist_InvariantMassRho_After       = new TH1F("InvariantMassRho_After",    ";M_{'#pi^{+}'#pi^{-}} (GeV)  ;Events/0.01 GeV",        500,  0.0,   5.0);

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("RunNumber");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("Entry");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("Combo");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("WeightFactor");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("MissingPMinus");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("BeamEnergy");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("MinusT");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaCM");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_Thrown");
} 
// END OF INITIALIZATION

Bool_t DSelector_piminus_p_12C_sim::Process(Long64_t locEntry)
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
        TLorentzVector locMissingP4 = locPiMinusP4 + locProtonP4 - locBeamP4;
        
        // GET THROWN P4
        TLorentzVector locBeamP4_Thrown, locPiMinusP4_Thrown, locProtonP4_Thrown, locMissingP4_Thrown;
        if(dThrownBeam != NULL)
            locBeamP4_Thrown = dThrownBeam->Get_P4();
        for(UInt_t loc_j = 0; loc_j < Get_NumThrown(); ++loc_j)
        {                                                                                                                                                                                                                                                     
            dThrownWrapper->Set_ArrayIndex(loc_j);
            if (dThrownWrapper->Get_PID() == PiMinus)
                locPiMinusP4_Thrown = dThrownWrapper->Get_P4();
            else if (dThrownWrapper->Get_PID() == Proton)
                locProtonP4_Thrown  = dThrownWrapper->Get_P4();
        }
        locMissingP4_Thrown = locPiMinusP4_Thrown + locProtonP4_Thrown - locBeamP4_Thrown;

		// GET RF TIMING INFO
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
        
        // DISCARD EVENTS FROM NEAREST OUT-OF-TIME BUNCH
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)                                                                                    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
        { 
		 	dComboWrapper->Set_IsComboCut(true); 
			continue;
        }

        // VECTORS in C.M. FRAME
        TVector3       boostCM          = (locPiMinusP4   + locProtonP4).BoostVector();
        TLorentzVector locBeamP4CM      = locBeamP4;
        TLorentzVector locPiMinusP4CM   = locPiMinusP4;
        locBeamP4CM.Boost(-boostCM);
        locPiMinusP4CM.Boost(-boostCM);

        TLorentzVector locProtonP4AsPion;
        locProtonP4AsPion.SetXYZM(locProtonP4.X(), locProtonP4.Y(), locProtonP4.Z(), 0.13957);

        // CALCULATE CUSTOM VARIABLES
        double locVertexR               = sqrt(pow(dComboBeamWrapper->Get_X4().X(),2) + pow(dComboBeamWrapper->Get_X4().Y(),2));
        double locChiSquarePerNDF       = dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit();
        double locSqrtS                 = (locPiMinusP4   + locProtonP4).Mag();
        double locMinusT                = -(locBeamP4   - locPiMinusP4).Mag2();
        double locMinusU                = -(locBeamP4   - locProtonP4).Mag2();
        double locThetaCM               = locBeamP4CM.Vect().Angle(locPiMinusP4CM.Vect())*RadToDeg;
        double locCoplanarity           = abs(locPiMinusP4.Phi()   - locProtonP4.Phi())*RadToDeg;

        // FILL CUSTOM HISTOGRAMS: BEFORE CUTS
        dHist_NumUnusedTracks_Before  ->Fill(dComboWrapper->Get_NumUnusedTracks(),                                         locHistAccidWeightFactor);
        dHist_NumUnusedShowers_Before ->Fill(dComboWrapper->Get_NumUnusedShowers(),                                        locHistAccidWeightFactor);
        dHist_PiMinusPIDFOM_Before    ->Fill(dPiMinusWrapper->Get_PIDFOM(),                                                locHistAccidWeightFactor);
        dHist_ProtonPIDFOM_Before     ->Fill(dProtonWrapper->Get_PIDFOM(),                                                 locHistAccidWeightFactor);
        dHist_PiMinusPVsTheta_Before  ->Fill(locPiMinusP4.P(),                            locPiMinusP4.Theta()*RadToDeg,   locHistAccidWeightFactor);     
        dHist_ProtonPVsTheta_Before   ->Fill(locProtonP4.P(),                             locProtonP4.Theta()*RadToDeg,    locHistAccidWeightFactor);    
        dHist_ConfidenceLevel_Before  ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(),                                  locHistAccidWeightFactor);
        dHist_MissingMomentum_Before  ->Fill(locMissingP4.P(),                                                             locHistAccidWeightFactor);
        dHist_EnergyBalance_Before    ->Fill(mass_neutron-locMissingP4.E(),                                                locHistAccidWeightFactor);
        dHist_PhotonEnergy_Before     ->Fill(locBeamP4.E(),                                                                locHistAccidWeightFactor);       
        dHist_VertexZ_Before          ->Fill(dComboBeamWrapper->Get_X4().Z(),                                              locHistAccidWeightFactor);
        dHist_VertexXY_Before         ->Fill(dComboBeamWrapper->Get_X4().X(),             dComboBeamWrapper->Get_X4().Y(), locHistAccidWeightFactor);
        dHist_MissingPMinus_Before    ->Fill(locMissingP4.Minus(),                                                         locHistAccidWeightFactor);  
        dHist_InvariantMassRho_Before ->Fill((locPiMinusP4 + locProtonP4AsPion).M(),                                       locHistAccidWeightFactor);

        // SET CUT FLAGS
        bool locCutFlags[9] = {false, false, false, false, false, false, false, false, false};
        int locIsComboCut = 0;
        if(dComboWrapper->Get_NumUnusedTracks()        > 0)                                                                  locCutFlags[0] = true;
        if(dComboWrapper->Get_NumUnusedShowers()       > 0)                                                                  locCutFlags[1] = true;
        if(dPiMinusWrapper->Get_PIDFOM()               < 0.01 || dProtonWrapper->Get_PIDFOM()    < 0.01)                     locCutFlags[2] = true;
        if((locPiMinusP4 + locProtonP4AsPion).M()      < 1.0)                                                                locCutFlags[3] = true;        
        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < 0.001)                                                              locCutFlags[4] = true;
        if(locMissingP4.P()                            > 2.0)                                                                locCutFlags[5] = true;
        if(locBeamP4.E()                               < 5.5  || locBeamP4.E()                   > 11.0)                     locCutFlags[6] = true;
        if(dComboBeamWrapper->Get_X4().Z()             < 51.0 || dComboBeamWrapper->Get_X4().Z() > 79.0 || locVertexR > 1.0) locCutFlags[7] = true;
        if(locMissingP4.Minus()                        < 0.4  || locMissingP4.Minus()            > 1.4)                      locCutFlags[8] = true;
        for(int loc_j = 0; loc_j < 9; ++loc_j)
        {
            locIsComboCut += locCutFlags[loc_j];
        }
            
        // FILL CUSTOM HISTOGRAMS: CUT EFFECTS
        if((locIsComboCut-locCutFlags[0]) == 0) dHist_CutEffect_NumUnusedTracks   ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[1]) == 0) dHist_CutEffect_NumUnusedShowers  ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[2]) == 0) dHist_CutEffect_TrackPIDFOM       ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[3]) == 0) dHist_CutEffect_ParticleKinematics->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[4]) == 0) dHist_CutEffect_ConfidenceLevel   ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[5]) == 0) dHist_CutEffect_MissingMomentum   ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[6]) == 0) dHist_CutEffect_PhotonEnergy      ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[7]) == 0) dHist_CutEffect_CommonVertex      ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);
        if((locIsComboCut-locCutFlags[8]) == 0) dHist_CutEffect_MissingPMinus     ->Fill(locMissingP4.Minus(), locHistAccidWeightFactor);

        // DISCARD EVENTS THAT FAIL CUTS
        if(locIsComboCut > 0)
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }
        
        // FILL CUSTOM HISTOGRAMS: AFTER CUTS
        dHist_NumUnusedTracks_After    ->Fill(dComboWrapper->Get_NumUnusedTracks(),                                                 locHistAccidWeightFactor);
        dHist_NumUnusedShowers_After   ->Fill(dComboWrapper->Get_NumUnusedShowers(),                                                locHistAccidWeightFactor);
        dHist_PiMinusPIDFOM_After      ->Fill(dPiMinusWrapper->Get_PIDFOM(),                                                        locHistAccidWeightFactor);
        dHist_ProtonPIDFOM_After       ->Fill(dProtonWrapper->Get_PIDFOM(),                                                         locHistAccidWeightFactor);
        dHist_PiMinusPVsTheta_After    ->Fill(locPiMinusP4.P(),                            locPiMinusP4.Theta()*RadToDeg,           locHistAccidWeightFactor);     
        dHist_ProtonPVsTheta_After     ->Fill(locProtonP4.P(),                             locProtonP4.Theta()*RadToDeg,            locHistAccidWeightFactor);    
        dHist_ConfidenceLevel_After    ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(),                                          locHistAccidWeightFactor);        
        dHist_MissingMomentum_After    ->Fill(locMissingP4.P(),                                                                     locHistAccidWeightFactor); 
        dHist_PhotonEnergy_After       ->Fill(locBeamP4.E(),                                                                        locHistAccidWeightFactor);       
        dHist_VertexZ_After            ->Fill(dComboBeamWrapper->Get_X4().Z(),                                                      locHistAccidWeightFactor);
        dHist_VertexXY_After           ->Fill(dComboBeamWrapper->Get_X4().X(),             dComboBeamWrapper->Get_X4().Y(),         locHistAccidWeightFactor);
        dHist_MissingPMinus_After      ->Fill(locMissingP4.Minus(),                                                                 locHistAccidWeightFactor); 
        dHist_PhotonTiming_After       ->Fill(locDeltaT_RF,                                                                         locHistAccidWeightFactor);       
        dHist_ChiSquarePerNDF_After    ->Fill(locChiSquarePerNDF,                                                                   locHistAccidWeightFactor);
        dHist_MissingMassSquared_After ->Fill(locMissingP4.M2(),                                                                    locHistAccidWeightFactor); 
        dHist_SqrtS_After              ->Fill(locSqrtS,                                                                             locHistAccidWeightFactor);    
        dHist_MinusT_After             ->Fill(locMinusT,                                                                            locHistAccidWeightFactor);    
        dHist_MinusU_After             ->Fill(locMinusU,                                                                            locHistAccidWeightFactor);      
        dHist_ThetaCM_After            ->Fill(locThetaCM,                                                                           locHistAccidWeightFactor);
        dHist_Coplanarity_After        ->Fill(locCoplanarity,                                                                       locHistAccidWeightFactor);      
        dHist_InvariantMassRho_After   ->Fill((locPiMinusP4 + locProtonP4AsPion).M(),                                               locHistAccidWeightFactor);

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        //FILL CUSTOM BRANCHES: FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Int_t>("RunNumber",                locRunNumber);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("Entry",                    locEntry);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("Combo",                    loc_i);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("WeightFactor",          locHistAccidWeightFactor);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("MissingPMinus",         locMissingP4.Minus());
        dFlatTreeInterface->Fill_Fundamental<Double_t>("BeamEnergy",            locBeamP4.E());
        dFlatTreeInterface->Fill_Fundamental<Double_t>("MinusT",                locMinusT);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaCM",               locThetaCM);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("BeamP4",              locBeamP4);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("PiMinusP4",           locPiMinusP4);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("ProtonP4",            locProtonP4);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("MissingP4",           locMissingP4);
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

void DSelector_piminus_p_12C_sim::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION


