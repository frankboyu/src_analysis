#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>

using namespace std;

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

double rad_to_deg = 180.0 / TMath::Pi();

double sigma_jpsi_p(double Ephoton, double t) {
    std::ifstream ifs("/work/halld2/home/boyu/src_analysis/selection/configs/dvmp_xsec_proton2025.csv");
    if (Ephoton < 8.2) return 0.0; // table only goes down to 8.2 GeV, so return 0 below that

    if (!ifs.is_open()) return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    double bestE=0, bestT=0, bestXS=std::numeric_limits<double>::quiet_NaN();
    double bestDist = std::numeric_limits<double>::infinity();
    bool any = false;

    while (std::getline(ifs, line)) {
        // trim leading spaces
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos==std::string::npos) continue;
        if (line[pos]=='#' || line[pos]=='%' || line[pos]=='q') continue; // skip comments

        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')) {
            // trim token
            size_t a = token.find_first_not_of(" \t\r\n");
            size_t b = token.find_last_not_of(" \t\r\n");
            if (a==std::string::npos) cols.push_back(""); else cols.push_back(token.substr(a, b-a+1));
        }

        try {
            double e = std::stod(cols[0]);
            double tt = std::stod(cols[2]);
            double xs = std::stod(cols[3]);

            any = true;

            // keep nearest (Euclidean) if no exact match
            double dx = e - Ephoton;
            double dt = tt + t;  // note: table has t as positive, but function argument is negative, so add to get distance
            double dist = std::sqrt(dx*dx + dt*dt);
            if (dist < bestDist) {
                bestDist = dist;
                bestE = e; bestT = tt; bestXS = xs;
            }
        } catch (...) {
            continue;
        }
    }

    if (!any) return std::numeric_limits<double>::quiet_NaN();
    return bestXS;
}

class DSelector_jpsi_p_recon : public DSelector
{
public:

    DSelector_jpsi_p_recon(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_jpsi_p_recon(){}

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
    DChargedTrackHypothesis*    dPositronWrapper;
    DChargedTrackHypothesis*    dElectronWrapper;
    DChargedTrackHypothesis*    dProtonWrapper;

    // CUSTOM HISTOGRAMS
    TH1D* dHist_NumUnusedTracks;
    TH1D* dHist_NumUnusedShowers;
    TH1D* dHist_PhotonEnergy;
    TH2D* dHist_PositronKinematics;
    TH2D* dHist_PositronKinematics_Thrown;
    TH2D* dHist_ElectronKinematics;
    TH2D* dHist_ElectronKinematics_Thrown;
    TH2D* dHist_ProtonKinematics;
    TH2D* dHist_ProtonKinematics_Thrown;
    TH2D* dHist_ProtondEdxCDC;
    TH1D* dHist_InvariantMassJpsi;
    TH2D* dHist_ChiSqPerNDF;
    TH1D* dHist_Minust;

    ClassDef(DSelector_jpsi_p_recon, 0);
};

void DSelector_jpsi_p_recon::Get_ComboWrappers(void)
{
	dStep0Wrapper       = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper   = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPositronWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
    dElectronWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper      = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
}

void DSelector_jpsi_p_recon::Init(TTree *locTree)
{

    // SET OUTPUT FILE NAME
    dOutputFileName          = "selectedhist_jpsi_p_recon.root";
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
    dHist_NumUnusedTracks    = new TH1D("NumUnusedTracks",    ";Unused Tracks             ;Events/1",             10,     0.0,    10.0);
    dHist_NumUnusedShowers   = new TH1D("NumUnusedShowers",   ";Unused Showers            ;Events/1",             10,     0.0,    10.0);
    dHist_PhotonEnergy       = new TH1D("PhotonEnergy",       ";Photon Energy (GeV)       ;Events/0.01 GeV",      900,    3.0,    12.0);
    dHist_PositronKinematics    = new TH2D("PositronKinematics",    ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_PositronKinematics_Thrown    = new TH2D("PositronKinematics_Thrown",    ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_ElectronKinematics   = new TH2D("ElectronKinematics",   ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_ElectronKinematics_Thrown   = new TH2D("ElectronKinematics_Thrown",   ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_ProtonKinematics = new TH2D("ProtonKinematics", ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_ProtonKinematics_Thrown = new TH2D("ProtonKinematics_Thrown", ";P (GeV/c)                 ;#theta (deg)",         100,    0.0,    10.0,   180,    0.0,    180.0);
    dHist_ProtondEdxCDC    = new TH2D("ProtondEdxCDC", ";P (GeV/c)                 ;dE/dx CDC (MeV/cm)", 100,    0.0,    10.0,   100,    0.0,    40.0);
    dHist_InvariantMassJpsi   = new TH1D("InvariantMassJpsi",	";M_{K^{+}K^{-}} (GeV)      ;Events/0.01 GeV",      500,    0.0,    5.0);
    dHist_ChiSqPerNDF        = new TH2D("ChiSqPerNDF",        ";M_{K^{+}K^{-}} (GeV)      ;log(#chi^{2}/NDF)",    400,    0.9,    4.9,    100,    0.0,    5);
    dHist_Minust              = new TH1D("Minust",              ";-t (GeV^{2})             ;Events/0.01 GeV^{2}", 100,    0.0,    5.0);
}
// END OF INITIALIZATION

Bool_t DSelector_jpsi_p_recon::Process(Long64_t locEntry)
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

        // DISCARD SIMULATION EVENTS WITH NO L1 TRIGGER BITS
        if (dComboWrapper->Get_L1TriggerBits() == 0)
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }

		// GET PARTICLE INDICES
		Int_t locBeamID             = dComboBeamWrapper->Get_BeamID();
        Int_t locPositronTrackID       = dPositronWrapper->Get_TrackID();
		Int_t locElectronTrackID      = dElectronWrapper->Get_TrackID();
		Int_t locProtonTrackID    = dProtonWrapper->Get_TrackID();

		// GET RECONSTRUCTED P4
        TLorentzVector locBeamP4        = dComboBeamWrapper->Get_P4();
        TLorentzVector locPositronP4       = dPositronWrapper->Get_P4();
        TLorentzVector locElectronP4      = dElectronWrapper->Get_P4();
        TLorentzVector locProtonP4    = dProtonWrapper->Get_P4();

        // GET RECONSTRUCTED MEASURED P4
        TLorentzVector locBeamP4_Measured       = dComboBeamWrapper->Get_P4_Measured();
        TLorentzVector locPositronP4_Measured      = dPositronWrapper->Get_P4_Measured();
        TLorentzVector locElectronP4_Measured     = dElectronWrapper->Get_P4_Measured();
        TLorentzVector locProtonP4_Measured   = dProtonWrapper->Get_P4_Measured();

        //GET THROWN P4 AND TOPOLOGY
        TLorentzVector locBeamX4_Thrown, locPositronX4_Thrown, locElectronX4_Thrown, locProtonX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locPositronP4_Thrown, locElectronP4_Thrown, locProtonP4_Thrown;
        TString locThrownTopology = Get_ThrownTopologyString();
        Int_t locThrownTopologyFlag = -1;
        if (dIsMC)
        {
            locBeamX4_Thrown = dThrownBeam->Get_X4();
            locBeamP4_Thrown = dThrownBeam->Get_P4();
            for(UInt_t loc_j = 0; loc_j < Get_NumThrown(); ++loc_j)
            {
                dThrownWrapper->Set_ArrayIndex(loc_j);
                if (dThrownWrapper->Get_PID() == Positron)
                {
                    locPositronX4_Thrown = dThrownWrapper->Get_X4();
                    locPositronP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == Electron)
                {
                    locElectronX4_Thrown = dThrownWrapper->Get_X4();
                    locElectronP4_Thrown = dThrownWrapper->Get_P4();
                }
                else if (dThrownWrapper->Get_PID() == Proton)
                {
                    locProtonX4_Thrown = dThrownWrapper->Get_X4();
                    locProtonP4_Thrown = dThrownWrapper->Get_P4();
                }
            }
        }

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

        double event_weight = locHistAccidWeightFactor*locComboAccidWeightFactor*sigma_jpsi_p(locBeamP4_Thrown.E(), (locBeamP4_Thrown - locPositronP4_Thrown - locElectronP4_Thrown - locProtonP4_Thrown).M2());

        if ((locPositronP4 + locElectronP4).M() < 3.0 || (locPositronP4 + locElectronP4).M() > 3.2)
        {
            dComboWrapper->Set_IsComboCut(true);
            continue;
        }
        // FILL HISTOGRAMS BEFORE CUTS
        dHist_NumUnusedTracks    ->Fill(dComboWrapper->Get_NumUnusedTracks(), event_weight);
        dHist_NumUnusedShowers   ->Fill(dComboWrapper->Get_NumUnusedShowers(), event_weight);
        dHist_PhotonEnergy       ->Fill(locBeamP4_Measured.E(), event_weight);
        dHist_PositronKinematics    ->Fill(locPositronP4_Measured.P(), locPositronP4_Measured.Theta()*rad_to_deg, event_weight);
        dHist_PositronKinematics_Thrown    ->Fill(locPositronP4_Thrown.P(), locPositronP4_Thrown.Theta()*rad_to_deg, event_weight);
        dHist_ElectronKinematics   ->Fill(locElectronP4_Measured.P(), locElectronP4_Measured.Theta()*rad_to_deg, event_weight);
        dHist_ElectronKinematics_Thrown   ->Fill(locElectronP4_Thrown.P(), locElectronP4_Thrown.Theta()*rad_to_deg, event_weight);
        dHist_ProtonKinematics ->Fill(locProtonP4_Measured.P(), locProtonP4_Measured.Theta()*rad_to_deg, event_weight);
        dHist_ProtonKinematics_Thrown ->Fill(locProtonP4_Thrown.P(), locProtonP4_Thrown.Theta()*rad_to_deg, event_weight);
        dHist_ProtondEdxCDC    ->Fill(locProtonP4_Measured.P(), dProtonWrapper->Get_dEdx_CDC()*1e6, event_weight);
        dHist_InvariantMassJpsi   ->Fill((locPositronP4+locElectronP4).M(), event_weight);
        dHist_ChiSqPerNDF        ->Fill((locPositronP4+locElectronP4).M(), TMath::Log10(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit()), event_weight);
        dHist_Minust              ->Fill(-1*(locBeamP4 - locPositronP4 - locElectronP4).M2(), event_weight);

		// EXECUTE ANALYSIS ACTIONS
        if(!Execute_Actions()) // if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

	}
    // END OF COMBO LOOP

	return kTRUE;
}
// END OF PROCESSING

void DSelector_jpsi_p_recon::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION