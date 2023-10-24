#include "DSelector_gd_kpkmp.h"

double RadToDeg = 180.0/3.1415926;
double mass_proton = 0.938272;

void DSelector_gd_kpkmp::Init(TTree *locTree){
    
	//SET OUTPUT FILE NAME AND FORMAT
	dOutputFileName = "test.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

	//INITIALIZE THE TREE INTERFACE AND WRAPPERS
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree, gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit
	Get_ComboWrappers();
	dPreviousRunNumber = 0;
    
    //ANALYSIS ACTIONS
    std::deque<Particle_t> MyPhi;
    MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);
    dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false, ""));
    dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false, ""));
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions(dAnalysisActions, dComboWrapper, false, 0, MyPhi, 200, 0.5, 2.5, "CutActionEffect");

	//INITIALIZE ANALYSIS ACTIONS
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); //If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well

    dHist_MultiCombo       = new TH1I("MultiCombo", ";Number of combos", 10, 0, 10);
    
	//CUSTOM HISTOGRAMS: BEFORE THE CUTS
    //Combo
    dHist_InvariantMass_Measured_Before            = new TH1F("InvariantMass_Measured_Before", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_KinFit_Before              = new TH1F("InvariantMass_KinFit_Before", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_ConfidenceLevel_KinFit_Before            = new TH1F("ConfidenceLevel_KinFit_Before", ";Confidence Level", 1000, 0.0, 1.0);
    dHist_ChiSq_KinFit_Before                      = new TH1F("ChiSq_KinFit_Before", ";ChiSq/NDF", 1000, 0.0, 10.0);
    dHist_ChiSqComparison_KinFit_Before            = new TH2F("ChiSqComparison_KinFit_Before", ";ChiSqKaon;ChiSqPion", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	dHist_MissingMassSquared_Measured_Before       = new TH1F("MissingMassSquared_Measured_Before", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_MissingMassSquared_KinFit_Before         = new TH1F("MissingMassSquared_KinFit_Before", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_EnergyBalance_Measured_Before            = new TH1F("EnergyBalance_Measured_Before", ";Energy Balance (GeV)", 800, -4.0, 4.0);
    dHist_EnergyBalance_KinFit_Before              = new TH1F("EnergyBalance_KinFit_Before", ";Energy Balance (GeV)", 800, -4.0, 4.0);
    dHist_VertexZ_KinFit_Before                    = new TH1F("VertexZ_KinFit_Before", ";Vertex Z (cm)", 200, 0.0, 200.0);
    dHist_VertexXY_KinFit_Before                   = new TH2F("VertexXY_KinFit_Before", ";Vertex X (cm);Vertex Y (cm)", 100, -5.0, 5.0, 100, -5.0, 5.0);
    dHist_MinusT_Measured_Before                   = new TH1F("MinusT_Measured_Before", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusT_KinFit_Before                     = new TH1F("MinusT_KinFit_Before", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_Measured_Before                   = new TH1F("MinusU_Measured_Before", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_KinFit_Before                     = new TH1F("MinusU_KinFit_Before", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_Coplanarity_Measured_Before              = new TH1F("Coplanarity_Measured_Before", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    dHist_Coplanarity_KinFit_Before                = new TH1F("Coplanarity_KinFit_Before", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    //Beam photon
    dHist_PhotonEnergy_Measured_Before             = new TH1F("PhotonEnergy_Before", ";Beam Energy (GeV)", 600, 0.0, 12.0);
    dHist_PhotonTiming_Measured_Before             = new TH1F("PhotonTiming_Before", ";#Delta t_{Beam-RF} (ns)", 360, -18.0, 18.0);
    //Detected particles
    dHist_KPlusPVsTheta_Measured_Before            = new TH2F("KPlusPVsTheta_Measured_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KPlusPVsTheta_KinFit_Before              = new TH2F("KPlusPVsTheta_KinFit_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_Measured_Before           = new TH2F("KMinusPVsTheta_Measured_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_KinFit_Before             = new TH2F("KMinusPVsTheta_KinFit_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_Measured_Before           = new TH2F("ProtonPVsTheta_Measured_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_KinFit_Before             = new TH2F("ProtonPVsTheta_KinFit_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    //Undetected particles
    dHist_MissingNeutronPVsTheta_Measured_Before   = new TH2F("MissingNeutronPVsTheta_Measured_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_MissingNeutronPVsTheta_KinFit_Before     = new TH2F("MissingNeutronPVsTheta_KinFit_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    dHist_InitialProtonPVsTheta_Measured_Before    = new TH2F("InitialProtonPVsTheta_Measured_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_InitialProtonPVsTheta_KinFit_Before      = new TH2F("InitialProtonPVsTheta_KinFit_Before", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    
    //CUSTOM HISTOGRAMS: DURING THE CUTS
    dHist_InvariantMass_ConfidenceLevelCut         = new TH1F("InvariantMass_ConfidenceLevelCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_MissingMassSquaredCut      = new TH1F("InvariantMass_MissingMassSquaredCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_EnergyBalanceCut           = new TH1F("InvariantMass_EnergyBalanceCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_CoplanarityCut             = new TH1F("InvariantMass_CoplanarityCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_UnusedTracksCut            = new TH1F("InvariantMass_UnusedTracksCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_UnusedShowersCut           = new TH1F("InvariantMass_UnusedShowersCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_VertexCut                  = new TH1F("InvariantMass_VertexCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_PhotonEnergyCut            = new TH1F("InvariantMass_PhotonEnergyCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_InvariantMassCut           = new TH1F("InvariantMass_InvariantMassCut", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    
    //CUSTOM HISTOGRAMS: AFTER THE CUTS
    //Combo
    dHist_InvariantMass_Measured_After             = new TH1F("InvariantMass_Measured_After", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_KinFit_After               = new TH1F("InvariantMass_KinFit_After", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_ConfidenceLevel_KinFit_After             = new TH1F("ConfidenceLevel_KinFit_After", ";Confidence Level", 1000, 0.0, 1.0);
    dHist_ChiSq_KinFit_After                       = new TH1F("ChiSq_KinFit_After", ";ChiSq/NDF", 1000, 0.0, 10.0);
    dHist_ChiSqComparison_KinFit_After             = new TH2F("ChiSqComparison_KinFit_After", ";ChiSqKaon;ChiSqPion", 1000, 0.0, 10.0, 1000, 0.0, 10.0);
	dHist_MissingMassSquared_Measured_After        = new TH1F("MissingMassSquared_Measured_After", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_MissingMassSquared_KinFit_After          = new TH1F("MissingMassSquared_KinFit_After", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_EnergyBalance_Measured_After             = new TH1F("EnergyBalance_Measured_After", ";Energy Balance (GeV)", 800, -4.0, 4.0);
    dHist_EnergyBalance_KinFit_After               = new TH1F("EnergyBalance_KinFit_After", ";Energy Balance (GeV)", 800, -4.0, 4.0);
    dHist_VertexZ_KinFit_After                     = new TH1F("VertexZ_KinFit_After", ";Vertex Z (cm)", 200, 0.0, 200.0);
    dHist_VertexXY_KinFit_After                    = new TH2F("VertexXY_KinFit_After", ";Vertex X (cm);Vertex Y (cm)", 100, -5.0, 5.0, 100, -5.0, 5.0);
    dHist_MinusT_Measured_After                    = new TH1F("MinusT_Measured_After", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusT_KinFit_After                      = new TH1F("MinusT_KinFit_After", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_Measured_After                    = new TH1F("MinusU_Measured_After", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_KinFit_After                      = new TH1F("MinusU_KinFit_After", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_Coplanarity_Measured_After               = new TH1F("Coplanarity_Measured_After", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    dHist_Coplanarity_KinFit_After                 = new TH1F("Coplanarity_KinFit_After", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    //Beam photon
    dHist_PhotonEnergy_Measured_After              = new TH1F("PhotonEnergy_After", ";Beam Energy (GeV)", 600, 0.0, 12.0);
    dHist_PhotonTiming_Measured_After              = new TH1F("PhotonTiming_After", ";#Delta t_{Beam-RF} (ns)", 360, -18.0, 18.0);
    //Detected particles
    dHist_KPlusPVsTheta_Measured_After             = new TH2F("KPlusPVsTheta_Measured_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KPlusPVsTheta_KinFit_After               = new TH2F("KPlusPVsTheta_KinFit_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_Measured_After            = new TH2F("KMinusPVsTheta_Measured_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_KinFit_After              = new TH2F("KMinusPVsTheta_KinFit_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_Measured_After            = new TH2F("ProtonPVsTheta_Measured_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_KinFit_After              = new TH2F("ProtonPVsTheta_KinFit_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    //Undetected particles
    dHist_MissingNeutronPVsTheta_Measured_After    = new TH2F("MissingNeutronPVsTheta_Measured_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_MissingNeutronPVsTheta_KinFit_After      = new TH2F("MissingNeutronPVsTheta_KinFit_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    dHist_InitialProtonPVsTheta_Measured_After     = new TH2F("InitialProtonPVsTheta_Measured_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_InitialProtonPVsTheta_KinFit_After       = new TH2F("InitialProtonPVsTheta_KinFit_After", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    
    //CUSTOM HISTOGRAMS: WITH MULTI-COMBO WEIGHT
    //Combo
    dHist_InvariantMass_Measured_Weighted          = new TH1F("InvariantMass_Measured_Weighted", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_InvariantMass_KinFit_Weighted            = new TH1F("InvariantMass_KinFit_Weighted", ";Invariant Mass (GeV)", 200, 0.5, 2.5);
    dHist_ConfidenceLevel_KinFit_Weighted          = new TH1F("ConfidenceLevel_KinFit_Weighted", ";Confidence Level", 1000, 0.0, 1.0);
	dHist_MissingMassSquared_Measured_Weighted     = new TH1F("MissingMassSquared_Measured_Weighted", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_MissingMassSquared_KinFit_Weighted       = new TH1F("MissingMassSquared_KinFit_Weighted", ";Missing Mass Squared (GeV^{2})", 400, 0.0, 4.0);
    dHist_EnergyBalance_Measured_Weighted          = new TH1F("EnergyBalance_Measured_Weighted", ";Energy Balance (GeV)", 800, -4.0, 4.0);
    dHist_EnergyBalance_KinFit_Weighted            = new TH1F("EnergyBalance_KinFit_Weighted", ";Energy Balance (GeV)", 800, -4.0, 4.0);    
    dHist_VertexZ_KinFit_Weighted                  = new TH1F("VertexZ_KinFit_Weighted", ";Vertex Z (cm)", 200, 0.0, 200.0);
    dHist_VertexXY_KinFit_Weighted                 = new TH2F("VertexXY_KinFit_Weighted", ";Vertex X (cm);Vertex Y (cm)", 100, -5.0, 5.0, 100, -5.0, 5.0);
    dHist_MinusT_Measured_Weighted                 = new TH1F("MinusT_Measured_Weighted", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusT_KinFit_Weighted                   = new TH1F("MinusT_KinFit_Weighted", ";-t (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_Measured_Weighted                 = new TH1F("MinusU_Measured_Weighted", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_MinusU_KinFit_Weighted                   = new TH1F("MinusU_KinFit_Weighted", ";-u (GeV)", 500, 0.0, 5.0);
    dHist_Coplanarity_Measured_Weighted            = new TH1F("Coplanarity_Measured_Weighted", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    dHist_Coplanarity_KinFit_Weighted              = new TH1F("Coplanarity_KinFit_Weighted", ";Coplanarity Angle (deg)", 360, 0.0, 360.0);
    //Beam photon
    dHist_PhotonEnergy_Measured_Weighted           = new TH1F("PhotonEnergy_Weighted", ";Beam Energy (GeV)", 600, 0.0, 12.0);
    dHist_PhotonTiming_Measured_Weighted           = new TH1F("PhotonTiming_Weighted", ";#Delta t_{Beam-RF} (ns)", 360, -18.0, 18.0);
    //Detected particles
    dHist_KPlusPVsTheta_Measured_Weighted          = new TH2F("KPlusPVsTheta_Measured_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KPlusPVsTheta_KinFit_Weighted            = new TH2F("KPlusPVsTheta_KinFit_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_Measured_Weighted         = new TH2F("KMinusPVsTheta_Measured_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_KMinusPVsTheta_KinFit_Weighted           = new TH2F("KMinusPVsTheta_KinFit_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_Measured_Weighted         = new TH2F("ProtonPVsTheta_Measured_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_ProtonPVsTheta_KinFit_Weighted           = new TH2F("ProtonPVsTheta_KinFit_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    //Undetected particles
    dHist_MissingNeutronPVsTheta_Measured_Weighted = new TH2F("MissingNeutronPVsTheta_Measured_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_MissingNeutronPVsTheta_KinFit_Weighted   = new TH2F("MissingNeutronPVsTheta_KinFit_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    dHist_InitialProtonPVsTheta_Measured_Weighted  = new TH2F("InitialProtonPVsTheta_Measured_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);
    dHist_InitialProtonPVsTheta_KinFit_Weighted    = new TH2F("InitialProtonPVsTheta_KinFit_Weighted", ";Momentum (GeV);Theta (deg)", 1000, 0.0, 10.0, 180, 0.0, 180.0);    
    
	//CUSTOM OUTPUT BRANCHES IN THE MAIN TREE
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	//dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	//dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	//dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	//dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	//dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");

	//CUSTOM OUTPUT BRANCHES IN THE FLAT TREE
	//dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");
    //dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("multicombo");
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	//dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	//dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	//dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	//dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");

	//CUSTOM BRANCHES TO READ
	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	//DETERMINE IF ANALYZING SIMULATED DATA
	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);} //end of initialization

Bool_t DSelector_gd_kpkmp::Process(Long64_t locEntry){
    
	//GET THE EVENT
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;

	//GET POLARIZATION ORIENTATION
	UInt_t locRunNumber = Get_RunNumber();  //Only if the run number changes. RCDB environment must be setup in order for this to work! (Will return false otherwise)
	if(locRunNumber != dPreviousRunNumber){
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;}

	//RESET ANALYSIS ACTIONS
	Reset_Actions_NewEvent();
	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
    int locMultiCombo = 0;
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	//EXAMPLE 1: Particle-specific info:
	set<Int_t> locUsedSoFar_BeamEnergy; //Int_t: Unique ID for beam particles. set: easy to use, fast to search

	//EXAMPLE 2: Combo-specific info:
		//In general: Could have multiple particles with the same PID: Use a set of Int_t's
		//In general: Multiple PIDs, so multiple sets: Contain within a map
		//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;

	//FILL CUSTOM OUTPUT BRANCHES FOR THIS EVENT
	//Int_t locMyInt = 7;
	//dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);
	//TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	//dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);
	//for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
	//	dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index

	//LOOP OVER COMBOS
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i){
        
		//GET THE COMBO
		dComboWrapper->Set_ComboIndex(loc_i);
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		//GET PARTICLE INFO    
        //Tracking ID
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID = dKMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();      
        //P4's: is KinFit if KinFit performed, else is measured
		TLorentzVector locBeamP4                    = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4                   = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4                  = dKMinusWrapper->Get_P4();
		TLorentzVector locProtonP4                  = dProtonWrapper->Get_P4();
        TLorentzVector locMissingNeutronP4          = locBeamP4 + dTargetP4 - locKPlusP4 - locKMinusP4 - locProtonP4;
        TLorentzVector locInitialProtonP4           = locKPlusP4 + locKMinusP4 + locProtonP4 - locBeamP4;
        //Measured P4's:
		TLorentzVector locBeamP4_Measured           = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured          = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured         = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured         = dProtonWrapper->Get_P4_Measured();
        TLorentzVector locMissingNeutronP4_Measured = locBeamP4_Measured + dTargetP4 - locKPlusP4_Measured - locKMinusP4_Measured - locProtonP4_Measured;
        TLorentzVector locInitialProtonP4_Measured = locKPlusP4 + locKMinusP4 + locProtonP4 - locBeamP4;

		//GET RF TIMING INFO
		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
		Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		Int_t locNumOutOfTimeBunchesInTree = 4; //YOU need to specify this number. Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 
		Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		Double_t locHistAccidWeightFactor = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
        if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		 	dComboWrapper->Set_IsComboCut(true); 
			continue;}
        
        //CALCULATE CUSTOM VARIABLES
        double locInvariantMass_Measured      = (locKPlusP4_Measured + locKMinusP4_Measured).M();
        double locInvariantMass_KinFit        = (locKPlusP4 + locKMinusP4).M();
		double locMissingMassSquared_Measured = (locBeamP4_Measured + dTargetP4 - locKPlusP4_Measured - locKMinusP4_Measured - locProtonP4_Measured).M2();
        double locMissingMassSquared_KinFit   = (locBeamP4 + dTargetP4 - locKPlusP4 - locKMinusP4 - locProtonP4).M2();
        double locEnergyBalance_Measured      = (locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured).E() - locBeamP4_Measured.E() - mass_proton;
        double locEnergyBalance_KinFit        = (locKPlusP4 + locKMinusP4 + locProtonP4).E() - locBeamP4.E() - mass_proton;   
        double locMinusT_Measured             = -(locBeamP4_Measured - locKPlusP4_Measured - locKMinusP4_Measured).Mag2();
        double locMinusT                      = -(locBeamP4 - locKPlusP4 - locKMinusP4).Mag2();
        double locMinusU_Measured             = -(locBeamP4_Measured - locProtonP4_Measured).Mag2();
        double locMinusU                      = -(locBeamP4 - locProtonP4).Mag2();
        double locCoplanarityAngle_Measured   = abs((locKPlusP4_Measured + locKMinusP4_Measured).Phi() - locProtonP4_Measured.Phi())*RadToDeg;        
        double locCoplanarityAngle_KinFit     = abs((locKPlusP4 + locKMinusP4).Phi() - locProtonP4.Phi())*RadToDeg;
        
        //FILL CUSTOM HISTOGRAMS: BEFORE CUTS
        //Combo
        dHist_InvariantMass_Measured_Before->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        dHist_InvariantMass_KinFit_Before->Fill(locInvariantMass_KinFit, locHistAccidWeightFactor);
        dHist_ConfidenceLevel_KinFit_Before->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(), locHistAccidWeightFactor);
        dHist_ChiSq_KinFit_Before->Fill(dComboWrapper->Get_ChiSq_KinFit("gd")/dComboWrapper->Get_NDF_KinFit("gd"), locHistAccidWeightFactor);
        if (dComboWrapper->Get_ChiSq_KinFit("gd") > 0)
            dHist_ChiSqComparison_KinFit_Before->Fill(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit(), dComboWrapper->Get_ChiSq_KinFit("gd")/dComboWrapper->Get_NDF_KinFit("gd"), locHistAccidWeightFactor);
        dHist_MissingMassSquared_Measured_Before->Fill(locMissingMassSquared_Measured, locHistAccidWeightFactor); 
        dHist_MissingMassSquared_KinFit_Before->Fill(locMissingMassSquared_KinFit, locHistAccidWeightFactor);
        dHist_EnergyBalance_Measured_Before->Fill(locEnergyBalance_Measured, locHistAccidWeightFactor);
        dHist_EnergyBalance_KinFit_Before->Fill(locEnergyBalance_KinFit, locHistAccidWeightFactor);
        dHist_VertexZ_KinFit_Before->Fill(locBeamX4_Measured.Z(), locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_Before->Fill(locBeamX4_Measured.X(), locBeamX4_Measured.Y(), locHistAccidWeightFactor);     
        dHist_MinusT_Measured_Before->Fill(locMinusT_Measured, locHistAccidWeightFactor);    
        dHist_MinusT_KinFit_Before->Fill(locMinusT, locHistAccidWeightFactor);    
        dHist_MinusU_Measured_Before->Fill(locMinusU_Measured, locHistAccidWeightFactor);      
        dHist_MinusU_KinFit_Before->Fill(locMinusU, locHistAccidWeightFactor);    
        dHist_Coplanarity_Measured_Before->Fill(locCoplanarityAngle_Measured, locHistAccidWeightFactor);      
        dHist_Coplanarity_KinFit_Before->Fill(locCoplanarityAngle_KinFit, locHistAccidWeightFactor);
        //Beam photon
        dHist_PhotonEnergy_Measured_Before->Fill(locBeamP4_Measured.E(), locHistAccidWeightFactor);       
        dHist_PhotonTiming_Measured_Before->Fill(locDeltaT_RF, locHistAccidWeightFactor);       
        //Detected particles
        dHist_KPlusPVsTheta_Measured_Before->Fill(locKPlusP4_Measured.P(), locKPlusP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);     
        dHist_KPlusPVsTheta_KinFit_Before->Fill(locKPlusP4.P(), locKPlusP4.Theta()*RadToDeg, locHistAccidWeightFactor);       
        dHist_KMinusPVsTheta_Measured_Before->Fill(locKMinusP4_Measured.P(), locKMinusP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);    
        dHist_KMinusPVsTheta_KinFit_Before->Fill(locKMinusP4.P(), locKMinusP4.Theta()*RadToDeg, locHistAccidWeightFactor);      
        dHist_ProtonPVsTheta_Measured_Before->Fill(locProtonP4_Measured.P(), locProtonP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);    
        dHist_ProtonPVsTheta_KinFit_Before->Fill(locProtonP4.P(), locProtonP4.Theta()*RadToDeg, locHistAccidWeightFactor);
        
        //EXECUTE MANUAL CUTS
        //Confidence level cut
        if(dComboWrapper->Get_ConfidenceLevel_KinFit() < 0.01){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        dHist_InvariantMass_ConfidenceLevelCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Missing mass squared cut
        if((locMissingMassSquared_Measured < 0.80) || (locMissingMassSquared_Measured > 1.00)){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        dHist_InvariantMass_MissingMassSquaredCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Energy balance cut
        if((locEnergyBalance_Measured < -1.0) || (locEnergyBalance_Measured > 1.0)){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        dHist_InvariantMass_EnergyBalanceCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Coplanarity cut
        if((locCoplanarityAngle_Measured < 160.0) || (locCoplanarityAngle_Measured > 200.0)){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        dHist_InvariantMass_CoplanarityCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        if(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit() < dComboWrapper->Get_ChiSq_KinFit("gd")/dComboWrapper->Get_NDF_KinFit("gd")){
			dComboWrapper->Set_IsComboCut(true);
			continue;}
        dHist_InvariantMass_VertexCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Unused tracks cut
        //if(dComboWrapper->Get_NumUnusedTracks() > 0){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}
        //dHist_InvariantMass_UnusedTracksCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Unused showers cut
        //if(dComboWrapper->Get_NumUnusedShowers() > 0){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}
        //dHist_InvariantMass_UnusedShowersCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);        
        //Vertex cut
        //if((locBeamX4_Measured.Z() < 51.0) || (locBeamX4_Measured.Z() > 76.0) || sqrt(pow(locBeamX4_Measured.X(),2) + pow(locBeamX4_Measured.Y(),2)) > 1.0){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}    
        //dHist_InvariantMass_VertexCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Photon energy cut    
        //if((locBeamP4_Measured.E() < 6.0) || (locBeamP4_Measured.E() > 11.0)){
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;}   
        //dHist_InvariantMass_PhotonEnergyCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        //Invariant mass cut
        if((locInvariantMass_Measured < 1.00) || (locInvariantMass_Measured > 1.04)){
			dComboWrapper->Set_IsComboCut(true);
			continue;}  
        dHist_InvariantMass_InvariantMassCut->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);       

        //FILL CUSTOM HISTOGRAMS: AFTER CUTS
        //Combo
        dHist_InvariantMass_Measured_After->Fill(locInvariantMass_Measured, locHistAccidWeightFactor);
        dHist_InvariantMass_KinFit_After->Fill(locInvariantMass_KinFit, locHistAccidWeightFactor);
        dHist_ConfidenceLevel_KinFit_After->Fill(dComboWrapper->Get_ConfidenceLevel_KinFit(), locHistAccidWeightFactor);
        dHist_ChiSq_KinFit_After->Fill(dComboWrapper->Get_ChiSq_KinFit("gd")/dComboWrapper->Get_NDF_KinFit("gd"), locHistAccidWeightFactor);
        if (dComboWrapper->Get_ChiSq_KinFit("gd") < 0)
            dHist_ChiSqComparison_KinFit_After->Fill(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit(), 0.0, locHistAccidWeightFactor);
        else
            dHist_ChiSqComparison_KinFit_After->Fill(dComboWrapper->Get_ChiSq_KinFit()/dComboWrapper->Get_NDF_KinFit(), dComboWrapper->Get_ChiSq_KinFit("gd")/dComboWrapper->Get_NDF_KinFit("gd"), locHistAccidWeightFactor);
        dHist_MissingMassSquared_Measured_After->Fill(locMissingMassSquared_Measured, locHistAccidWeightFactor); 
        dHist_MissingMassSquared_KinFit_After->Fill(locMissingMassSquared_KinFit, locHistAccidWeightFactor); 
        dHist_EnergyBalance_Measured_After->Fill(locEnergyBalance_Measured, locHistAccidWeightFactor);
        dHist_EnergyBalance_KinFit_After->Fill(locEnergyBalance_KinFit, locHistAccidWeightFactor);
        dHist_VertexZ_KinFit_After->Fill(locBeamX4_Measured.Z(), locHistAccidWeightFactor);
        dHist_VertexXY_KinFit_After->Fill(locBeamX4_Measured.X(), locBeamX4_Measured.Y(), locHistAccidWeightFactor);     
        dHist_MinusT_Measured_After->Fill(locMinusT_Measured, locHistAccidWeightFactor);    
        dHist_MinusT_KinFit_After->Fill(locMinusT, locHistAccidWeightFactor);    
        dHist_MinusU_Measured_After->Fill(locMinusU_Measured, locHistAccidWeightFactor);      
        dHist_MinusU_KinFit_After->Fill(locMinusU, locHistAccidWeightFactor);    
        dHist_Coplanarity_Measured_After->Fill(locCoplanarityAngle_Measured, locHistAccidWeightFactor);      
        dHist_Coplanarity_KinFit_After->Fill(locCoplanarityAngle_KinFit, locHistAccidWeightFactor);
        //Beam photon
        dHist_PhotonEnergy_Measured_After->Fill(locBeamP4_Measured.E(), locHistAccidWeightFactor);       
        dHist_PhotonTiming_Measured_After->Fill(locDeltaT_RF, locHistAccidWeightFactor);       
        //Detected particles
        dHist_KPlusPVsTheta_Measured_After->Fill(locKPlusP4_Measured.P(), locKPlusP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);     
        dHist_KPlusPVsTheta_KinFit_After->Fill(locKPlusP4.P(), locKPlusP4.Theta()*RadToDeg, locHistAccidWeightFactor);       
        dHist_KMinusPVsTheta_Measured_After->Fill(locKMinusP4_Measured.P(), locKMinusP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);    
        dHist_KMinusPVsTheta_KinFit_After->Fill(locKMinusP4.P(), locKMinusP4.Theta()*RadToDeg, locHistAccidWeightFactor);      
        dHist_ProtonPVsTheta_Measured_After->Fill(locProtonP4_Measured.P(), locProtonP4_Measured.Theta()*RadToDeg, locHistAccidWeightFactor);    
        dHist_ProtonPVsTheta_KinFit_After->Fill(locProtonP4.P(), locProtonP4.Theta()*RadToDeg, locHistAccidWeightFactor); 
        
		//EXECUTE ANALYSIS ACTIONS
		dAnalyzeCutActions->Perform_Action(); // Loop through the analysis actions, executing them in order for the active particle combo. Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;
        
        locMultiCombo += 1;
        
        //Beam energy
		//if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end()){
		//	dHist_BeamEnergy->Fill(locBeamP4.E(),locHistAccidWeightFactor);
		//	locUsedSoFar_BeamEnergy.insert(locBeamID);}
		//Missing Mass Squared
        //TLorentzVector locMissingP4_Measured = (locBeamP4_Measured + dTargetP4) - (locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured);
		//double locMissingMassSquared = locMissingP4_Measured.M2();
		//map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		//locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		//locUsedThisCombo_MissingMass[KPlus].insert(locKPlusTrackID);
		//locUsedThisCombo_MissingMass[KMinus].insert(locKMinusTrackID);
		//locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		//if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())  //unique missing mass combo: histogram it, and register this combo of particles
		//{
		//	dHist_MissingMassSquared->Fill(locMissingMassSquared,locHistAccidWeightFactor);
		//	locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		//}
        
        //FILL CUSTOM OUTPUT BRANCHES FOR THIS COMBO
		//TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
		//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		//dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		//dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
        
		//FILL FLAT TREE
		//dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);
		//Int_t locMyInt_Flat = 7;
		//dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);
		//TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		//dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);
		//for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		//{
		//	dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
		//	TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
		//	dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		//}

		//WRITE TO THE FLAT TREE
		//Fill_FlatTree();
    }  //end of combo loop
    
    if (locMultiCombo > 0)
        dHist_MultiCombo->Fill(locMultiCombo);

	//FILL HISTOGRAMS OF SURVIVING COMBOS/EVENTS
	//Fill_NumCombosSurvivedHists();

	//LOOP OVER THROWN DATA (OPTIONAL)
    //Thrown beam: just use directly
	//if(dThrownBeam != NULL)
	//	double locEnergy = dThrownBeam->Get_P4().E();
	//Loop over throwns
	//for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i){
		//Set branch array indices corresponding to this particle
	//	dThrownWrapper->Set_ArrayIndex(loc_i);
		//Do stuff with the wrapper here ...}
    
	//LOOP OVER OTHER ARRAYS (OPTIONAL)
    /*
	//Loop over beam particles (note, only those appearing in combos are present)
	for(UInt_t loc_i = 0; loc_i < Get_NumBeam(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dBeamWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
	//Loop over charged track hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumChargedHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dChargedHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
	//Loop over neutral particle hypotheses
	for(UInt_t loc_i = 0; loc_i < Get_NumNeutralHypos(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dNeutralHypoWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
    */

	//FILL CLONE OF TTREE HERE WITH CUTS APPLIED
	//Bool_t locIsEventCut = true;
	//for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
	//	//Set branch array indices for combo and all combo particles
	//	dComboWrapper->Set_ComboIndex(loc_i);
	//	// Is used to indicate when combos have been cut
	//	if(dComboWrapper->Get_IsComboCut())
	//		continue;
	//	locIsEventCut = false; // At least one combo succeeded
	//	break;
	//}
	//if(!locIsEventCut && dOutputTreeFileName != "")
	//	Fill_OutputTree();

	return kTRUE;}  //end of processing

void DSelector_gd_kpkmp::Finalize(void){
	//CALL THIS LAST
	DSelector::Finalize();} //end of finalization