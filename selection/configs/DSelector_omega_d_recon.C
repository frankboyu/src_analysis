#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

class DSelector_omega_d_recon : public DSelector
{
public:

    DSelector_omega_d_recon(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_omega_d_recon(){}

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
    DChargedTrackHypothesis*    dPiPlusWrapper;
    DChargedTrackHypothesis*    dPiMinusWrapper;
    DChargedTrackHypothesis*    dDeuteronWrapper;
    DParticleComboStep*         dStep1Wrapper;
    DKinematicData*             dDecayingPi0Wrapper;
    DNeutralParticleHypothesis* dPhoton1Wrapper;
    DNeutralParticleHypothesis* dPhoton2Wrapper;

    // CUSTOM HISTOGRAMS
    TH1F* dHist_NumUnusedTracks_Before;
    TH1F* dHist_NumUnusedShowers_Before;
    TH1F* dHist_PhotonEnergy_Before;
    TH1F* dHist_VertexZ_Before;
    TH2F* dHist_VertexXY_Before;
    TH1F* dHist_ConfidenceLevel_Before;
    TH1F* dHist_ShowerQuality_Before;
    TH1F* dHist_PiPlusPIDFOM_Before;
    TH1F* dHist_PiMinusPIDFOM_Before;
    TH1F* dHist_DeuterondEdxCDC_Before;
    TH1F* dHist_InvariantMassOmega_Before;
    TH1F* dHist_ThrownTopology_Before;

    TH1F* dHist_NumUnusedTracks_After;
    TH1F* dHist_NumUnusedShowers_After;
    TH1F* dHist_PhotonEnergy_After;
    TH1F* dHist_VertexZ_After;
    TH2F* dHist_VertexXY_After;
    TH1F* dHist_ConfidenceLevel_After;
    TH1F* dHist_ShowerQuality_After;
    TH1F* dHist_PiPlusPIDFOM_After;
    TH1F* dHist_PiMinusPIDFOM_After;
    TH1F* dHist_DeuterondEdxCDC_After;
    TH1F* dHist_InvariantMassOmega_After;
    TH1F* dHist_ThrownTopology_After;

    TH1F* dHist_PhotonTiming_Raw;
    TH1F* dHist_PhotonTiming_Weighted;

    TH1F* dHist_NumUnusedTracks_Weighted;
    TH1F* dHist_NumUnusedShowers_Weighted;
    TH1F* dHist_PhotonEnergy_Weighted;
    TH1F* dHist_VertexZ_Weighted;
    TH2F* dHist_VertexXY_Weighted;
    TH1F* dHist_ConfidenceLevel_Weighted;
    TH1F* dHist_ShowerQuality_Weighted;
    TH1F* dHist_PiPlusPIDFOM_Weighted;
    TH1F* dHist_PiMinusPIDFOM_Weighted;
    TH1F* dHist_DeuterondEdxCDC_Weighted;
    TH1F* dHist_InvariantMassOmega_Weighted;
    TH1F* dHist_ThrownTopology_Weighted;

	ClassDef(DSelector_omega_d_recon, 0);
};

void DSelector_omega_d_recon::Get_ComboWrappers(void)
{
	dStep0Wrapper       = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper   = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiPlusWrapper      = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
    dPiMinusWrapper     = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dDeuteronWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(3));
    dStep1Wrapper       = dComboWrapper->Get_ParticleComboStep(1);
	dDecayingPi0Wrapper = dStep1Wrapper->Get_InitialParticle();
	dPhoton1Wrapper     = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(0));
	dPhoton2Wrapper     = static_cast<DNeutralParticleHypothesis*>(dStep1Wrapper->Get_FinalParticle(1));
}

void DSelector_omega_d_recon::Init(TTree *locTree)
{
    // DETERMINE THE TAG NAME
    if (locTree->GetBranch("MCWeight") == NULL)
        dTag += "data";
    else
        dTag += "sim";

    if      (string(locTree->GetName()).find("gd") != string::npos)
        dTag += "_2H";
    else if (string(locTree->GetName()).find("ghe") != string::npos)
        dTag += "_4He";
    else if (string(locTree->GetName()).find("gc12") != string::npos)
        dTag += "_12C";

    if      (string(locTree->GetName()).find("pi0pippimd_") != string::npos)
        dTag += "_exc";
    else if (string(locTree->GetName()).find("pi0pippimdmiss") != string::npos)
        dTag += "_IA";
    else if (string(locTree->GetName()).find("pi0pippimdinc_") != string::npos)
        dTag += "_inc";

    // SET OUTPUT FILE NAME
    dOutputFileName          = Form("selectedhist_omega_d_recon_%s.root", dTag.c_str());
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = Form("selectedtree_omega_d_recon_%s.root", dTag.c_str());
    dFlatTreeName            = "selectedtree_omega_d_recon";
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
    dHist_NumUnusedTracks_Before        = new TH1F("NumUnusedTracks_Before",        ";Unused Tracks                     ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers_Before       = new TH1F("NumUnusedShowers_Before",       ";Unused Showers                    ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy_Before           = new TH1F("PhotonEnergy_Before",           ";Photon Energy (GeV)               ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_VertexZ_Before                = new TH1F("VertexZ_Before",                ";Vertex Z (cm)                     ;Events/1 cm",          200,    0.0,    200.0);
    dHist_VertexXY_Before               = new TH2F("VertexXY_Before",               ";Vertex X (cm)                     ;Vertex Y (cm)",        100,    -5.0,   5.0,    100,    -5.0,   5.0);
    dHist_ConfidenceLevel_Before        = new TH1F("ConfidenceLevel_Before",        ";Confidence Level                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_ShowerQuality_Before          = new TH1F("ShowerQuality_Before",          ";Shower Quality                    ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiPlusPIDFOM_Before           = new TH1F("PiPlusPIDFOM_Before",           ";PIDFOM_{#pi^{+}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiMinusPIDFOM_Before          = new TH1F("PiMinusPIDFOM_Before",          ";PIDFOM_{#pi^{-}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_DeuterondEdxCDC_Before        = new TH1F("DeuterondEdxCDC_Before",        ";dE/dx_{CDC} (keV/cm)              ;Events/0.1 keV/cm",    400,    0.0,    40.0);
    dHist_InvariantMassOmega_Before     = new TH1F("InvariantMassOmega_Before",	    ";M_{#pi^{0}#pi^{+}#pi^{-}} (GeV)   ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ThrownTopology_Before         = new TH1F("ThrownTopology_Before",         ";Thrown Topology                   ;Events/1",             20,     0.5,    20.5);

    dHist_NumUnusedTracks_After         = new TH1F("NumUnusedTracks_After",         ";Unused Tracks                     ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers_After        = new TH1F("NumUnusedShowers_After",        ";Unused Showers                    ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy_After            = new TH1F("PhotonEnergy_After",            ";Photon Energy (GeV)               ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_VertexZ_After                 = new TH1F("VertexZ_After",                 ";Vertex Z (cm)                     ;Events/1 cm",          200,    0.0,    200.0);
    dHist_VertexXY_After                = new TH2F("VertexXY_After",                ";Vertex X (cm)                     ;Vertex Y (cm)",        100,    -5.0,   5.0,    100,    -5.0,   5.0);
    dHist_ConfidenceLevel_After         = new TH1F("ConfidenceLevel_After",         ";Confidence Level                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_ShowerQuality_After           = new TH1F("ShowerQuality_After",           ";Shower Quality                    ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiPlusPIDFOM_After            = new TH1F("PiPlusPIDFOM_After",            ";PIDFOM_{#pi^{+}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiMinusPIDFOM_After           = new TH1F("PiMinusPIDFOM_After",           ";PIDFOM_{#pi^{-}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_DeuterondEdxCDC_After         = new TH1F("DeuterondEdxCDC_After",         ";dE/dx_{CDC} (keV/cm)              ;Events/0.1 keV/cm",    400,    0.0,    40.0);
    dHist_InvariantMassOmega_After      = new TH1F("InvariantMassOmega_After",	    ";M_{#pi^{0}#pi^{+}#pi^{-}} (GeV)   ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ThrownTopology_After          = new TH1F("ThrownTopology_After",          ";Thrown Topology                   ;Events/1",             20,     0.5,    20.5);

    dHist_PhotonTiming_Raw              = new TH1F("PhotonTiming_Raw",              ";#Delta t_{Beam-RF} (ns)           ;Events/0.1 ns",        360,    -18.0,  18.0);
    dHist_PhotonTiming_Weighted         = new TH1F("PhotonTiming_Weighted",         ";#Delta t_{Beam-RF} (ns)           ;Events/0.1 ns",        360,    -18.0,  18.0);

    dHist_NumUnusedTracks_Weighted      = new TH1F("NumUnusedTracks_Weighted",      ";Unused Tracks                     ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers_Weighted     = new TH1F("NumUnusedShowers_Weighted",     ";Unused Showers                    ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy_Weighted         = new TH1F("PhotonEnergy_Weighted",         ";Photon Energy (GeV)               ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_VertexZ_Weighted              = new TH1F("VertexZ_Weighted",              ";Vertex Z (cm)                     ;Events/1 cm",          200,    0.0,    200.0);
    dHist_VertexXY_Weighted             = new TH2F("VertexXY_Weighted",             ";Vertex X (cm)                     ;Vertex Y (cm)",        100,    -5.0,   5.0,    100,    -5.0,   5.0);
    dHist_ConfidenceLevel_Weighted      = new TH1F("ConfidenceLevel_Weighted",      ";Confidence Level                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_ShowerQuality_Weighted        = new TH1F("ShowerQuality_Weighted",        ";Shower Quality                    ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiPlusPIDFOM_Weighted         = new TH1F("PiPlusPIDFOM_Weighted",         ";PIDFOM_{#pi^{+}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_PiMinusPIDFOM_Weighted        = new TH1F("PiMinusPIDFOM_Weighted",        ";PIDFOM_{#pi^{-}}                  ;Events/0.001",         1000,   0.0,    1.0);
    dHist_DeuterondEdxCDC_Weighted      = new TH1F("DeuterondEdxCDC_Weighted",      ";dE/dx_{CDC} (keV/cm)              ;Events/0.1 keV/cm",    400,    0.0,    40.0);
    dHist_InvariantMassOmega_Weighted   = new TH1F("InvariantMassOmega_Weighted",   ";M_{#pi^{0}#pi^{+}#pi^{-}} (GeV)   ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ThrownTopology_Weighted       = new TH1F("ThrownTopology_Weighted",       ";Thrown Topology                   ;Events/1",             20,     0.5,    20.5);

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("thrown_topology");
    dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("polarization_angle");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pip_pidfom");  // the PIDFOM in the default flat branches pip_pid_fom is corrupted and always 0
	dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("pim_pidfom");  // the PIDFOM in the default flat branches pim_pid_fom is corrupted and always 0
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("g1_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("g1_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("g2_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("g2_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("decaypi0_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("decaypi0_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pip_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pip_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_p4_truth");
}
// END OF INITIALIZATION

Bool_t DSelector_omega_d_recon::Process(Long64_t locEntry)
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
        Int_t locPiPlusTrackID      = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID     = dPiMinusWrapper->Get_TrackID();
		Int_t locDeuteronTrackID    = dDeuteronWrapper->Get_TrackID();
        Int_t locPhoton1NeutralID   = dPhoton1Wrapper->Get_NeutralID();
        Int_t locPhoton2NeutralID   = dPhoton2Wrapper->Get_NeutralID();

		// GET RECONSTRUCTED P4
        TLorentzVector locBeamP4        = dComboBeamWrapper->Get_P4();
        TLorentzVector locPiPlusP4      = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4     = dPiMinusWrapper->Get_P4();
		TLorentzVector locDeuteronP4    = dDeuteronWrapper->Get_P4();
        TLorentzVector locDecayingPi0P4 = dDecayingPi0Wrapper->Get_P4();
		TLorentzVector locPhoton1P4     = dPhoton1Wrapper->Get_P4();
		TLorentzVector locPhoton2P4     = dPhoton2Wrapper->Get_P4();

        //GET THROWN P4 AND TOPOLOGY
        TLorentzVector locBeamX4_Thrown, locDecayingPi0X4_Thrown, locPhoton1X4_Thrown, locPhoton2X4_Thrown, locPiPlusX4_Thrown, locPiMinusX4_Thrown, locDeuteronX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locDecayingPi0P4_Thrown, locPhoton1P4_Thrown, locPhoton2P4_Thrown, locPiPlusP4_Thrown, locPiMinusP4_Thrown, locDeuteronP4_Thrown;
        TString locThrownTopology = Get_ThrownTopologyString();
        Int_t locThrownTopologyFlag = -1;
        if (dIsMC)
        {
            locBeamX4_Thrown = dThrownBeam->Get_X4();
            locBeamP4_Thrown = dThrownBeam->Get_P4();
            for(UInt_t loc_j = 0; loc_j < Get_NumThrown(); ++loc_j)
            {
                dThrownWrapper->Set_ArrayIndex(loc_j);
                if (dThrownWrapper->Get_PID() == Pi0)
                {
                    locDecayingPi0X4_Thrown = dThrownWrapper->Get_X4();
                    locDecayingPi0P4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == PiPlus)
                {
                    locPiPlusX4_Thrown = dThrownWrapper->Get_X4();
                    locPiPlusP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == PiMinus)
                {
                    locPiMinusX4_Thrown = dThrownWrapper->Get_X4();
                    locPiMinusP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == Deuteron)
                {
                    locDeuteronX4_Thrown = dThrownWrapper->Get_X4();
                    locDeuteronP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == Gamma)
                {
                    if (dPhoton1Wrapper->Get_ThrownIndex() == Int_t(loc_j))
                    {
                        locPhoton1X4_Thrown = dThrownWrapper->Get_X4();
                        locPhoton1P4_Thrown = dThrownWrapper->Get_P4();
                        cout << "Found thrown photon 1 with ID: " << locPhoton1NeutralID << endl;
                    }
                    else if (dPhoton1Wrapper->Get_ThrownIndex() == Int_t(loc_j))
                    {
                        locPhoton2X4_Thrown = dThrownWrapper->Get_X4();
                        locPhoton2P4_Thrown = dThrownWrapper->Get_P4();
                        cout << "Found thrown photon 2 with ID: " << locPhoton2NeutralID << endl;
                    }
                    else
                    {
                        cout << "Unexpected thrown photon with ID: " << dThrownWrapper->Get_PID() << endl;
                    }
                }
                else
                    cout << "Unexpected PID: " << dThrownWrapper->Get_PID() << endl;
            }
        }

        // FILL HISTOGRAMS BEFORE CUTS
        dHist_NumUnusedTracks_Before    ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_Before   ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_Before       ->Fill(locBeamP4.E());
        dHist_VertexZ_Before            ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_Before           ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_ConfidenceLevel_Before    ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        dHist_ShowerQuality_Before      ->Fill(dPhoton1Wrapper->Get_Shower_Quality());
        dHist_ShowerQuality_Before      ->Fill(dPhoton2Wrapper->Get_Shower_Quality());
        dHist_PiPlusPIDFOM_Before       ->Fill(dPiPlusWrapper->Get_PIDFOM());
        dHist_PiMinusPIDFOM_Before      ->Fill(dPiMinusWrapper->Get_PIDFOM());
        dHist_DeuterondEdxCDC_Before    ->Fill(dDeuteronWrapper->Get_dEdx_CDC()*1e6);
        dHist_InvariantMassOmega_Before ->Fill((locDecayingPi0P4+locPiPlusP4+locPiMinusP4).M());
        dHist_ThrownTopology_Before     ->Fill(locThrownTopology.Data(), 1);

        // PERFORM CUTS
        if(dComboWrapper->Get_NumUnusedTracks()         > 0)                                                dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_NumUnusedShowers()        > 0)                                                dComboWrapper->Set_IsComboCut(true);
        if(locBeamP4.E()                                < 5.8   || locBeamP4.E()                    > 10.7) dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Z()              < 51.0  || dComboBeamWrapper->Get_X4().Z()  > 79.0) dComboWrapper->Set_IsComboCut(true);
        if(dComboBeamWrapper->Get_X4().Perp()           > 1.0)                                              dComboWrapper->Set_IsComboCut(true);
        if(dComboWrapper->Get_ConfidenceLevel_KinFit()  < 1e-4)                                             dComboWrapper->Set_IsComboCut(true);
        if(dPhoton1Wrapper->Get_Shower_Quality()        < 0.5)                                              dComboWrapper->Set_IsComboCut(true);
        if(dPhoton2Wrapper->Get_Shower_Quality()        < 0.5)                                              dComboWrapper->Set_IsComboCut(true);
        if(dPiPlusWrapper->Get_PIDFOM()                 < 1e-4)                                             dComboWrapper->Set_IsComboCut(true);
        if(dPiMinusWrapper->Get_PIDFOM()                < 1e-4)                                             dComboWrapper->Set_IsComboCut(true);
        if(dDeuteronWrapper->Get_dEdx_CDC()             == 0.0)                                             dComboWrapper->Set_IsComboCut(true);
        if((locPiPlusP4+locPiMinusP4).M()               <= 0.0)                                             dComboWrapper->Set_IsComboCut(true);

        if(dComboWrapper->Get_IsComboCut())  continue;

        // FILL HISTOGRAMS AFTER CUTS
        dHist_NumUnusedTracks_After     ->Fill(dComboWrapper->Get_NumUnusedTracks());
        dHist_NumUnusedShowers_After    ->Fill(dComboWrapper->Get_NumUnusedShowers());
        dHist_PhotonEnergy_After        ->Fill(locBeamP4.E());
        dHist_VertexZ_After             ->Fill(dComboBeamWrapper->Get_X4().Z());
        dHist_VertexXY_After       	    ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y());
        dHist_ConfidenceLevel_After     ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit());
        dHist_ShowerQuality_After       ->Fill(dPhoton1Wrapper->Get_Shower_Quality());
        dHist_ShowerQuality_After       ->Fill(dPhoton2Wrapper->Get_Shower_Quality());
        dHist_PiPlusPIDFOM_After        ->Fill(dPiPlusWrapper->Get_PIDFOM());
        dHist_PiMinusPIDFOM_After       ->Fill(dPiMinusWrapper->Get_PIDFOM());
        dHist_DeuterondEdxCDC_After     ->Fill(dDeuteronWrapper->Get_dEdx_CDC()*1e6);
        dHist_InvariantMassOmega_After  ->Fill((locDecayingPi0P4+locPiPlusP4+locPiMinusP4).M());
        dHist_ThrownTopology_After      ->Fill(locThrownTopology.Data(), 1);

		// GET THE ACCIDENTAL WEIGHT FACTOR
		TLorentzVector locBeamX4                       = dComboBeamWrapper->Get_X4_Measured();
		Double_t       locBunchPeriod                  = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t       locDeltaT_RF                    = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4, dComboWrapper);
		Int_t          locRelBeamBucket                = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4, dComboWrapper);          // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t          locNumOutOfTimeBunchesInTree    = 4;                                                                                             // Number of out-of-time beam bunches in tree on a single side
		Bool_t         locSkipNearestOutOfTimeBunch    = true;                                                                                          // true: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t          locNumOutOfTimeBunchesToUse     = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree;
		Double_t       locAccidentalScalingFactor      = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC);         // ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t       locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E());           // ideal value would be 1, but deviations observed, need added factor.
		Double_t       locHistAccidWeightFactor        = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ;        // weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time

        dHist_PhotonTiming_Raw->Fill(locDeltaT_RF);
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1)                                                                                    // skip nearest out-of-time bunch: tails of in-time distribution also leak in
        {
            dComboWrapper->Set_IsComboCut(true);
			continue;
        }
        dHist_PhotonTiming_Weighted->Fill(locDeltaT_RF, locHistAccidWeightFactor);

		// FILL WEIGHTED HISTOGRAMS
		dHist_NumUnusedTracks_Weighted      ->Fill(dComboWrapper->Get_NumUnusedTracks(), locHistAccidWeightFactor);
		dHist_NumUnusedShowers_Weighted     ->Fill(dComboWrapper->Get_NumUnusedShowers(), locHistAccidWeightFactor);
		dHist_PhotonEnergy_Weighted         ->Fill(locBeamP4.E(), locHistAccidWeightFactor);
		dHist_VertexZ_Weighted              ->Fill(dComboBeamWrapper->Get_X4().Z(), locHistAccidWeightFactor);
		dHist_VertexXY_Weighted             ->Fill(dComboBeamWrapper->Get_X4().X(), dComboBeamWrapper->Get_X4().Y(), locHistAccidWeightFactor);
		dHist_ConfidenceLevel_Weighted      ->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(), locHistAccidWeightFactor);
        dHist_ShowerQuality_Weighted        ->Fill(dPhoton1Wrapper->Get_Shower_Quality(), locHistAccidWeightFactor);
        dHist_ShowerQuality_Weighted        ->Fill(dPhoton2Wrapper->Get_Shower_Quality(), locHistAccidWeightFactor);
		dHist_PiPlusPIDFOM_Weighted         ->Fill(dPiPlusWrapper->Get_PIDFOM(), locHistAccidWeightFactor);
		dHist_PiMinusPIDFOM_Weighted        ->Fill(dPiMinusWrapper->Get_PIDFOM(), locHistAccidWeightFactor);
        dHist_DeuterondEdxCDC_Weighted      ->Fill(dDeuteronWrapper->Get_dEdx_CDC()*1e6, locHistAccidWeightFactor);
        dHist_InvariantMassOmega_Weighted   ->Fill((locDecayingPi0P4+locPiPlusP4+locPiMinusP4).M(), locHistAccidWeightFactor);
        dHist_ThrownTopology_Weighted       ->Fill(locThrownTopology.Data(), locHistAccidWeightFactor);

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Int_t>("thrown_topology", locThrownTopologyFlag);
        dFlatTreeInterface->Fill_Fundamental<Int_t>("polarization_angle", dPolarizationAngle);
        dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight", locHistAccidWeightFactor);
		dFlatTreeInterface->Fill_Fundamental<Double_t>("pip_pidfom", dPiPlusWrapper->Get_PIDFOM());
		dFlatTreeInterface->Fill_Fundamental<Double_t>("pim_pidfom", dPiMinusWrapper->Get_PIDFOM());
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_x4_truth", locBeamX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_p4_truth", locBeamP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("g1_x4_truth", locPhoton1X4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("g1_p4_truth", locPhoton1P4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("g2_x4_truth", locPhoton2X4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("g2_p4_truth", locPhoton2P4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("decaypi0_x4_truth", locDecayingPi0X4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("decaypi0_p4_truth", locDecayingPi0P4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pip_x4_truth", locPiPlusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pip_p4_truth", locPiPlusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_x4_truth", locPiMinusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_p4_truth", locPiMinusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_x4_truth", locDeuteronX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_p4_truth", locDeuteronP4_Thrown);
        Fill_FlatTree();
	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_omega_d_recon::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION