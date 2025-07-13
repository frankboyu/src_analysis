#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

double rad_to_deg = 180.0 / TMath::Pi();

class DSelector_phi_d_recon : public DSelector
{
public:

    DSelector_phi_d_recon(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_phi_d_recon(){}

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
    string dTag;

    // PARTICLE WRAPPERS
    DParticleComboStep*         dStep0Wrapper;
    DBeamParticle*              dComboBeamWrapper;
    DChargedTrackHypothesis*    dKPlusWrapper;
    DChargedTrackHypothesis*    dKMinusWrapper;
    DChargedTrackHypothesis*    dDeuteronWrapper;

    // CUSTOM HISTOGRAMS
    TH1D* dHist_NumUnusedTracks_Before;
    TH1D* dHist_NumUnusedShowers_Before;
    TH1F* dHist_PhotonEnergy_Before;
    TH1F* dHist_VertexZ_Before;
    TH2F* dHist_VertexXY_Before;
    TH2F* dHist_KPlusKinematics_Before;
    TH2F* dHist_KMinusKinematics_Before;
    TH2F* dHist_DeuteronKinematics_Before;
    TH1F* dHist_KPlusPIDFOM_Before;
    TH1F* dHist_KMinusPIDFOM_Before;
    TH2F* dHist_DeuterondEdxCDC_Before;
    TH1F* dHist_InvariantMassPhi_Before;
    TH2F* dHist_ChiSqPerNDF_Before;
    TH1F* dHist_ThrownTopology_Before;

    TH1D* dHist_NumUnusedTracks_After;
    TH1D* dHist_NumUnusedShowers_After;
    TH1F* dHist_PhotonEnergy_After;
    TH1F* dHist_VertexZ_After;
    TH2F* dHist_VertexXY_After;
    TH2F* dHist_KPlusKinematics_After;
    TH2F* dHist_KMinusKinematics_After;
    TH2F* dHist_DeuteronKinematics_After;
    TH1F* dHist_KPlusPIDFOM_After;
    TH1F* dHist_KMinusPIDFOM_After;
    TH2F* dHist_DeuterondEdxCDC_After;
    TH1F* dHist_InvariantMassPhi_After;
    TH2F* dHist_ChiSqPerNDF_After;
    TH1F* dHist_ThrownTopology_After;

    ClassDef(DSelector_phi_d_recon, 0);
};

void DSelector_phi_d_recon::Get_ComboWrappers(void)
{
	dStep0Wrapper       = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper   = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper       = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
    dKMinusWrapper      = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dDeuteronWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
}

void DSelector_phi_d_recon::Init(TTree *locTree)
{
    // DETERMINE THE TAG NAME
    if      (string(locTree->GetName()).find("kpkmd_") != string::npos)
        dTag += "exc";
    else if (string(locTree->GetName()).find("kpkmdinc_") != string::npos)
        dTag += "inc";

    if (locTree->GetBranch("MCWeight") == NULL)
        dTag += "_data";
    else
        dTag += "_sim";

    if      (string(locTree->GetName()).find("gd") != string::npos)
        dTag += "_2H";
    else if (string(locTree->GetName()).find("ghe") != string::npos)
        dTag += "_4He";
    else if (string(locTree->GetName()).find("gc12") != string::npos)
        dTag += "_12C";

    // SET OUTPUT FILE NAME
    dOutputFileName          = Form("selectedhist_phi_d_recon_%s.root", dTag.c_str());
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = Form("selectedtree_phi_d_recon_%s.root", dTag.c_str());
    dFlatTreeName            = "selectedtree_phi_d_recon";
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
    dHist_NumUnusedTracks_Before    = new TH1D("NumUnusedTracks_Before",    ";Unused Tracks             ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers_Before   = new TH1D("NumUnusedShowers_Before",   ";Unused Showers            ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy_Before       = new TH1F("PhotonEnergy_Before",       ";Photon Energy (GeV)       ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_VertexZ_Before            = new TH1F("VertexZ_Before",            ";Vertex Z (cm)             ;Events/1 cm",          200,    0.0,    200.0);
    dHist_VertexXY_Before           = new TH2F("VertexXY_Before",           ";Vertex X (cm)             ;Vertex Y (cm)",        100,    -5.0,   5.0,    100,    -5.0,   5.0);
    dHist_KPlusKinematics_Before    = new TH2F("KPlusKinematics_Before",    ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_KMinusKinematics_Before   = new TH2F("KMinusKinematics_Before",   ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_DeuteronKinematics_Before = new TH2F("DeuteronKinematics_Before", ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_KPlusPIDFOM_Before        = new TH1F("KPlusPIDFOM_Before",        ";log(PIDFOM)               ;Events/1",             50,     -50.0,  0);
    dHist_KMinusPIDFOM_Before       = new TH1F("KMinusPIDFOM_Before",       ";log(PIDFOM)               ;Events/1",             50,     -50.0,  0);
    dHist_DeuterondEdxCDC_Before    = new TH2F("DeuterondEdxCDC_Before",    ";P (GeV/c)                 ;dE/dx_{CDC} (keV/cm)", 100,    0.0,    10.0,   400,    0.0,    40.0);
    dHist_InvariantMassPhi_Before   = new TH1F("InvariantMassPhi_Before",	";M_{K^{+}K^{-}} (GeV)      ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ChiSqPerNDF_Before        = new TH2F("ChiSqPerNDF_Before",        ";M_{K^{+}K^{-}} (GeV)      ;#Chi^2/NDF",           500,    0.0,    5.0,    20,     0.0,    20);
    dHist_ThrownTopology_Before     = new TH1F("ThrownTopology_Before",     ";Thrown Topology           ;Events/1",             20,     0.5,    20.5);

    dHist_NumUnusedTracks_After     = new TH1D("NumUnusedTracks_After",     ";Unused Tracks             ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers_After    = new TH1D("NumUnusedShowers_After",    ";Unused Showers            ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy_After        = new TH1F("PhotonEnergy_After",        ";Photon Energy (GeV)       ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_VertexZ_After             = new TH1F("VertexZ_After",             ";Vertex Z (cm)             ;Events/1 cm",          200,    0.0,    200.0);
    dHist_VertexXY_After            = new TH2F("VertexXY_After",            ";Vertex X (cm)             ;Vertex Y (cm)",        100,    -5.0,   5.0,    100,    -5.0,   5.0);
    dHist_KPlusKinematics_After     = new TH2F("KPlusKinematics_After",     ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_KMinusKinematics_After    = new TH2F("KMinusKinematics_After",    ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_DeuteronKinematics_After  = new TH2F("DeuteronKinematics_After",  ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   90,     0.0,    180.0);
    dHist_KPlusPIDFOM_After         = new TH1F("KPlusPIDFOM_After",         ";log(PIDFOM)               ;Events/1",             50,     -50.0,  0);
    dHist_KMinusPIDFOM_After        = new TH1F("KMinusPIDFOM_After",        ";log(PIDFOM)               ;Events/1",             50,     -50.0,  0);
    dHist_DeuterondEdxCDC_After     = new TH2F("DeuterondEdxCDC_After",     ";P (GeV/c)                 ;dE/dx_{CDC} (keV/cm)", 100,    0.0,    10.0,   400,    0.0,    40.0);
    dHist_InvariantMassPhi_After    = new TH1F("InvariantMassPhi_After",	";M_{K^{+}K^{-}} (GeV)      ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ChiSqPerNDF_After         = new TH2F("ChiSqPerNDF_After",         ";M_{K^{+}K^{-}} (GeV)      ;#Chi^2/NDF",           500,    0.0,    5.0,    20,     0.0,    20);
    dHist_ThrownTopology_After      = new TH1F("ThrownTopology_After",      ";Thrown Topology           ;Events/1",             20,     0.5,    20.5);

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("num_unused_tracks");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("num_unused_showers");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("beam_id");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("kp_id");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("km_id");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("d_id");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidental_weight");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("kp_pidfom");  // the PIDFOM in the default flat branches kp_pid_fom is corrupted and always 0
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("km_pidfom");  // the PIDFOM in the default flat branches km_pid_fom is corrupted and always 0
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("thrown_topology");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("polarization_angle");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("kp_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("kp_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("km_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("km_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_p4_truth");
}
// END OF INITIALIZATION

Bool_t DSelector_phi_d_recon::Process(Long64_t locEntry)
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

        //GET THROWN P4 AND TOPOLOGY
        TLorentzVector locBeamX4_Thrown, locKPlusX4_Thrown, locKMinusX4_Thrown, locDeuteronX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locKPlusP4_Thrown, locKMinusP4_Thrown, locDeuteronP4_Thrown;
        TString locThrownTopology = Get_ThrownTopologyString();
        Int_t locThrownTopologyFlag = -1;
        if (dIsMC)
        {
            locBeamX4_Thrown = dThrownBeam->Get_X4();
            locBeamP4_Thrown = dThrownBeam->Get_P4();
            for(UInt_t loc_j = 0; loc_j < Get_NumThrown(); ++loc_j)
            {
                dThrownWrapper->Set_ArrayIndex(loc_j);
                if (dThrownWrapper->Get_PID() == KPlus)
                {
                    locKPlusX4_Thrown = dThrownWrapper->Get_X4();
                    locKPlusP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == KMinus)
                {
                    locKMinusX4_Thrown = dThrownWrapper->Get_X4();
                    locKMinusP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == Deuteron)
                {
                    locDeuteronX4_Thrown = dThrownWrapper->Get_X4();
                    locDeuteronP4_Thrown = dThrownWrapper->Get_P4();
                }
            }
            locDeuteronX4_Thrown = locKPlusX4_Thrown;  // workaround for the missing deuteron info in the tree
            locDeuteronP4_Thrown = locBeamP4_Thrown + TLorentzVector(0, 0, 0, 1.875612859) - locKPlusP4_Thrown - locKMinusP4_Thrown; // workaround for the missing deuteron info in the tree
        }

        // FILL HISTOGRAMS BEFORE CUTS
        dHist_NumUnusedTracks_Before    ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_Before   ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_Before       ->Fill(locBeamP4.E());
        dHist_VertexZ_Before            ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_Before           ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_KPlusKinematics_Before    ->Fill(locKPlusP4_Measured.P(), locKPlusP4_Measured.Theta()*rad_to_deg);
        dHist_KMinusKinematics_Before   ->Fill(locKMinusP4_Measured.P(), locKMinusP4_Measured.Theta()*rad_to_deg);
        dHist_DeuteronKinematics_Before ->Fill(locDeuteronP4_Measured.P(), locDeuteronP4_Measured.Theta()*rad_to_deg);
        dHist_KPlusPIDFOM_Before        ->Fill(TMath::Log10(dKPlusWrapper->Get_PIDFOM()));
        dHist_KMinusPIDFOM_Before       ->Fill(TMath::Log10(dKMinusWrapper->Get_PIDFOM()));
        dHist_DeuterondEdxCDC_Before    ->Fill(locDeuteronP4_Measured.P(), dDeuteronWrapper->Get_dEdx_CDC()*1e6);
        dHist_InvariantMassPhi_Before   ->Fill((locKPlusP4+locKMinusP4).M());
        dHist_ChiSqPerNDF_Before        ->Fill((locKPlusP4+locKMinusP4).M(), dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit());
        dHist_ThrownTopology_Before     ->Fill(locThrownTopology.Data(), 1);

        // PERFORM CUTS
        if(locBeamP4.E()                        < 5.8   || locBeamP4.E()                                > 10.7) dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Z()      < 51.0  || dComboBeamWrapper->Get_X4().Z()              > 79.0) dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Perp()   > 1.0)                                                          dComboWrapper->Set_IsComboCut(true);
        if(locKPlusP4_Measured.P()              < 0.4   || locKPlusP4_Measured.Theta()*rad_to_deg       < 2)    dComboWrapper->Set_IsComboCut(true);
        if(locKMinusP4_Measured.P()             < 0.4   || locKMinusP4_Measured.Theta()*rad_to_deg      < 2)    dComboWrapper->Set_IsComboCut(true);
        if(locDeuteronP4_Measured.P()           < 0.4   || locDeuteronP4_Measured.Theta()*rad_to_deg    < 2)    dComboWrapper->Set_IsComboCut(true);
        if(dDeuteronWrapper->Get_dEdx_CDC()     == 0.0)                                                         dComboWrapper->Set_IsComboCut(true);
        if((locKPlusP4+locKMinusP4).M()         > 1.5)                                                          dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit()                            > 10)   dComboWrapper->Set_IsComboCut(true);

        if(dComboWrapper->Get_IsComboCut())  continue;

        // FILL HISTOGRAMS AFTER CUTS
        dHist_NumUnusedTracks_After     ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_After    ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_After        ->Fill(locBeamP4.E());
        dHist_VertexZ_After             ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_After            ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_KPlusKinematics_After     ->Fill(locKPlusP4_Measured.P(), locKPlusP4_Measured.Theta()*rad_to_deg);
        dHist_KMinusKinematics_After    ->Fill(locKMinusP4_Measured.P(), locKMinusP4_Measured.Theta()*rad_to_deg);
        dHist_DeuteronKinematics_After  ->Fill(locDeuteronP4_Measured.P(), locDeuteronP4_Measured.Theta()*rad_to_deg);
        dHist_KPlusPIDFOM_After         ->Fill(TMath::Log10(dKPlusWrapper->Get_PIDFOM()));
        dHist_KMinusPIDFOM_After        ->Fill(TMath::Log10(dKMinusWrapper->Get_PIDFOM()));
        dHist_DeuterondEdxCDC_After     ->Fill(locDeuteronP4_Measured.P(), dDeuteronWrapper->Get_dEdx_CDC()*1e6);
        dHist_InvariantMassPhi_After    ->Fill((locKPlusP4+locKMinusP4).M());
        dHist_ChiSqPerNDF_After         ->Fill((locKPlusP4+locKMinusP4).M(), dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit());
        dHist_ThrownTopology_After      ->Fill(locThrownTopology.Data(), 1);

		// GET THE ACCIDENTAL WEIGHT FACTOR
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

        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)                                                                                    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
        {
            dComboWrapper->Set_IsComboCut(true);
			continue;
        }

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Int_t>("num_unused_tracks", dComboWrapper->Get_NumUnusedTracks());
        dFlatTreeInterface->Fill_Fundamental<Int_t>("num_unused_showers", dComboWrapper->Get_NumUnusedShowers());
        dFlatTreeInterface->Fill_Fundamental<Int_t>("beam_id", locBeamID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("kp_id", locKPlusTrackID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("km_id", locKMinusTrackID);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("d_id", locDeuteronTrackID);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("accidental_weight", locHistAccidWeightFactor);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("kp_pidfom", dKPlusWrapper->Get_PIDFOM());
        dFlatTreeInterface->Fill_Fundamental<Double_t>("km_pidfom", dKMinusWrapper->Get_PIDFOM());
        dFlatTreeInterface->Fill_Fundamental<Int_t>("thrown_topology", locThrownTopologyFlag);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("polarization_angle", dPolarizationAngle);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_x4_truth", locBeamX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_p4_truth", locBeamP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("kp_x4_truth", locKPlusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("kp_p4_truth", locKPlusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("km_x4_truth", locKMinusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("km_p4_truth", locKMinusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_x4_truth", locDeuteronX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_p4_truth", locDeuteronP4_Thrown);
        Fill_FlatTree();
	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_phi_d_recon::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION