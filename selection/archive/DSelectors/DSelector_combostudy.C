#include "DSelector_combostudy.h"

void DSelector_combostudy::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "combostudy.root"; //"" for none
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen

	//Because this function gets called for each TTree in the TChain, we must be careful:
		//We need to re-initialize the tree interface & branch wrappers, but don't want to recreate histograms
	bool locInitializedPriorFlag = dInitializedFlag; //save whether have been initialized previously
	DSelector::Init(locTree); //This must be called to initialize wrappers for each new TTree
	//gDirectory now points to the output file with name dOutputFileName (if any)
	if(locInitializedPriorFlag)
		return; //have already created histograms, etc. below: exit

	Get_ComboWrappers();
	dPreviousRunNumber = 0;

	/*********************************** EXAMPLE USER INITIALIZATION: ANALYSIS ACTIONS **********************************/

	// EXAMPLE: Create deque for histogramming particle masses:
	// // For histogramming the phi mass in phi -> K+ K-
	// // Be sure to change this and dAnalyzeCutActions to match reaction
	std::deque<Particle_t> Lambda1520;
	Lambda1520.push_back(Proton); Lambda1520.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data
	
// 	dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper,1E-7));
// 	dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.2, 8.8));
// 	dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//PID
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));

	//MASSES
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, Lambda1520, 1000, 1.4, 2.4, "Lambda_Measured"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, Lambda1520, 1000, 1.4, 2.4, "Lambda_KinFit"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, 0, {KPlus,KMinus}, 1000, 0.8, 2.8, "Phi_Measured"));
	dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, true, 0, {KPlus,KMinus}, 1000, 0.8, 2.8, "Phi_KinFit"));
	dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change Lambda1520 to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, Lambda1520, 1000, 1.4, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/
	dHist_MultComboWeight = new TH1F("MultComboWeight", ";# combos NOT from beam", 10, 0, 10);
	
	//MEASURED HISTOGRAMS:
	dHist_BeamDeltaTRF_Measured = new TH2F("BeamDeltaTRF_Measured", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_Measured = new TH1F("BeamEnergy_Measured", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_Measured = new TH1F("MissingMassSquared_Measured", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_Measured = new TH1F("PKMinus_InvMass_Measured", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_Measured = new TH2F("BeamDeltaTRF_PKMinus_InvMass_Measured", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);

	//KINFIT HISTOGRAMS:
	dHist_BeamDeltaTRF_Kinfit = new TH2F("BeamDeltaTRF_Kinfit", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_Kinfit = new TH1F("BeamEnergy_Kinfit", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_Kinfit = new TH1F("MissingMassSquared_Kinfit", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_Kinfit = new TH1F("PKMinus_InvMass_Kinfit", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_Kinfit = new TH2F("BeamDeltaTRF_PKMinus_InvMass_Kinfit", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);

	//CORRECT COMBO HISTOGRAMS:
	dHist_BeamDeltaTRF_CorrectCombo = new TH2F("BeamDeltaTRF_CorrectCombo", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_CorrectCombo = new TH1F("BeamEnergy_CorrectCombo", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_CorrectCombo = new TH1F("MissingMassSquared_CorrectCombo", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_CorrectCombo = new TH1F("PKMinus_InvMass_CorrectCombo", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_CorrectCombo = new TH2F("BeamDeltaTRF_PKMinus_InvMass_CorrectCombo", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);

	//WRONG COMBO HISTOGRAMS:
	dHist_BeamDeltaTRF_WrongCombo = new TH2F("BeamDeltaTRF_WrongCombo", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_WrongCombo = new TH1F("BeamEnergy_WrongCombo", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_WrongCombo = new TH1F("MissingMassSquared_WrongCombo", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_WrongCombo = new TH1F("PKMinus_InvMass_WrongCombo", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_WrongCombo = new TH2F("BeamDeltaTRF_PKMinus_InvMass_WrongCombo", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);

	//THROWN HISTOGRAMS:
	dHist_BeamEnergy_Thrown = new TH1F("BeamEnergy_Thrown", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_Thrown = new TH1F("MissingMassSquared_Thrown", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_Thrown = new TH1F("PKMinus_InvMass_Thrown", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);

	//WEIGHTED HISTOGRAMS:
	dHist_BeamDeltaTRF_Measured_Weighted = new TH2F("BeamDeltaTRF_Measured_Weighted", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_Measured_Weighted = new TH1F("BeamEnergy_Measured_Weighted", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_Measured_Weighted = new TH1F("MissingMassSquared_Measured_Weighted", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_Measured_Weighted = new TH1F("PKMinus_InvMass_Measured_Weighted", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_Measured_Weighted = new TH2F("BeamDeltaTRF_PKMinus_InvMass_Measured_Weighted", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);
	
	dHist_BeamDeltaTRF_Kinfit_Weighted = new TH2F("BeamDeltaTRF_Kinfit_Weighted", ";Beam DeltaTRF (ns)", 200, -20, 20, 2, 0, 2);
	dHist_BeamEnergy_Kinfit_Weighted = new TH1F("BeamEnergy_Kinfit_Weighted", ";Beam Energy (GeV)", 900, 3, 12);
	dHist_MissingMassSquared_Kinfit_Weighted = new TH1F("MissingMassSquared_Kinfit_Weighted", ";Missing Mass Squared (GeV/c^{2})^{2}", 200, -0.1, 0.1);
	dHist_PKMinus_InvMass_Kinfit_Weighted = new TH1F("PKMinus_InvMass_Kinfit_Weighted", ";pK^{-} Inv. Mass (GeV/c^{2})", 500, 1.4, 6.4);
	dHist_BeamDeltaTRF_PKMinus_InvMass_Kinfit_Weighted = new TH2F("BeamDeltaTRF_PKMinus_InvMass_Kinfit_Weighted", ";Beam DeltaTRF (ns);pK^{-} Inv. Mass (GeV/c^{2})", 200, -20, 20, 500, 1.4, 6.4);
}

Bool_t DSelector_combostudy::Process(Long64_t locEntry)
{
	// The Process() function is called for each entry in the tree. The entry argument
	// specifies which entry in the currently loaded tree is to be processed.
	//
	// This function should contain the "body" of the analysis. It can contain
	// simple or elaborate selection criteria, run algorithms on the data
	// of the event and typically fill histograms.
	//
	// The processing can be stopped by calling Abort().
	// Use fStatus to set the return value of TTree::Process().
	// The return value is currently not used.

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry
	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
// 		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
		cout << "Run " << locRunNumber << endl;
	}
	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
		//Then for each combo, just compare to what you used before, and make sure it's unique

	set<Int_t>  locUsedSoFar_BeamEnergy_Measured;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass_Measured;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PKMinus_InvMass_Measured;
	
	set<Int_t>  locUsedSoFar_BeamEnergy_Kinfit;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass_Kinfit;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_PKMinus_InvMass_Kinfit;

	/************************** FOR WEIGHTS CALCULATIONS ********************************************/
	vector<Int_t> locCombosPassedCuts;
	set<map<Particle_t, set<Int_t> > > AdditionalCombos;  //similar to uniqueness tracking to make sure additional combo doesn't come from beam photons, we want to keep those and subtract with side peaks
	Double_t MultComboWeight = 0; //start at 0 because we will add at least 1 to it -> weight 1

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i)
	{
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locKPlusTrackID = dKPlusWrapper->Get_TrackID();
		Int_t locKMinusTrackID = dKMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locBeamX4 = dComboBeamWrapper->Get_X4();
		TLorentzVector locKPlusX4 = dKPlusWrapper->Get_X4();
		TLorentzVector locKMinusX4 = dKMinusWrapper->Get_X4();
		TLorentzVector locProtonX4 = dProtonWrapper->Get_X4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		TLorentzVector locKPlusX4_Measured = dKPlusWrapper->Get_X4_Measured();
		TLorentzVector locKMinusX4_Measured = dKMinusWrapper->Get_X4_Measured();
		TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		Execute_Actions();
// 		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
// 			continue;
		
		/************************************ HISTOGRAM MEASURED ************************************/

		//Evaluate iIsPrompt_Measured
		Double_t dRFTime_Measured = dComboWrapper->Get_RFTime_Measured();
		Double_t dTargetCenterZ = dComboWrapper->Get_TargetCenter().Z();
		Double_t locDeltaTRF_Measured = locBeamX4_Measured.T() - (dRFTime_Measured + (locBeamX4_Measured.Z() - dTargetCenterZ)/29.9792458);
		dHist_BeamDeltaTRF_Measured->Fill(locDeltaTRF_Measured,0);
		Double_t iIsPrompt_Measured=0;
		if((locDeltaTRF_Measured > -18.036) && (locDeltaTRF_Measured < -2.004)){
			iIsPrompt_Measured = 2;
			dHist_BeamDeltaTRF_Measured->Fill(locDeltaTRF_Measured,1,-0.125);
		}
		else if((locDeltaTRF_Measured >= -2.004) && (locDeltaTRF_Measured <= 2.004)){
			iIsPrompt_Measured = 1;
			dHist_BeamDeltaTRF_Measured->Fill(locDeltaTRF_Measured,1);
		}
		else if((locDeltaTRF_Measured > 2.004) && (locDeltaTRF_Measured < 18.036)){
			iIsPrompt_Measured = 2;
			dHist_BeamDeltaTRF_Measured->Fill(locDeltaTRF_Measured,1,-0.125);
		}
		else
			continue;
// 		iIsPrompt_Measured=1; //TODO deactivate accidental subtraction
		
		//Beam Energy
		if(locUsedSoFar_BeamEnergy_Measured.find(locBeamID) == locUsedSoFar_BeamEnergy_Measured.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			Fill1DHisto(dHist_BeamEnergy_Measured,locBeamP4_Measured.E(),iIsPrompt_Measured);
			locUsedSoFar_BeamEnergy_Measured.insert(locBeamID);
		}
		
		//Missing Mass Squared
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured;
		double locMissingMassSquared = locMissingP4_Measured.M2();
		
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass_Measured;
		locUsedThisCombo_MissingMass_Measured[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass_Measured[KPlus].insert(locKPlusTrackID);
		locUsedThisCombo_MissingMass_Measured[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_MissingMass_Measured[Proton].insert(locProtonTrackID);
		
		if(locUsedSoFar_MissingMass_Measured.find(locUsedThisCombo_MissingMass_Measured) == locUsedSoFar_MissingMass_Measured.end())
		{
			Fill1DHisto(dHist_MissingMassSquared_Measured,locMissingMassSquared,iIsPrompt_Measured);
			locUsedSoFar_MissingMass_Measured.insert(locUsedThisCombo_MissingMass_Measured);
		}
		
		//Invariant Mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_PKMinus_InvMass_Measured;
		locUsedThisCombo_PKMinus_InvMass_Measured[Unknown].insert(locBeamID); //beam for accidental subtraction
		locUsedThisCombo_PKMinus_InvMass_Measured[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_PKMinus_InvMass_Measured[Proton].insert(locProtonTrackID);
		
		if(locUsedSoFar_PKMinus_InvMass_Measured.find(locUsedThisCombo_PKMinus_InvMass_Measured) == locUsedSoFar_PKMinus_InvMass_Measured.end())
		{
			Fill1DHisto(dHist_PKMinus_InvMass_Measured,(locKMinusP4_Measured+locProtonP4_Measured).M(),iIsPrompt_Measured);
			Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_Measured,locDeltaTRF_Measured,(locKMinusP4_Measured+locProtonP4_Measured).M(),1);
			locUsedSoFar_PKMinus_InvMass_Measured.insert(locUsedThisCombo_PKMinus_InvMass_Measured);
		}
		
		
		if(dComboWrapper->Get_IsTrueCombo()){
			dHist_BeamDeltaTRF_CorrectCombo->Fill(locDeltaTRF_Measured,0);
			Fill1DHisto(dHist_BeamEnergy_CorrectCombo,locBeamP4_Measured.E(),iIsPrompt_Measured);
			Fill1DHisto(dHist_MissingMassSquared_CorrectCombo,locMissingMassSquared,iIsPrompt_Measured);
			Fill1DHisto(dHist_PKMinus_InvMass_CorrectCombo,(locKMinusP4_Measured+locProtonP4_Measured).M(),iIsPrompt_Measured);
			Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_CorrectCombo,locDeltaTRF_Measured,(locKMinusP4_Measured+locProtonP4_Measured).M(),1);
		}
		else{
			dHist_BeamDeltaTRF_WrongCombo->Fill(locDeltaTRF_Measured,0);
			Fill1DHisto(dHist_BeamEnergy_WrongCombo,locBeamP4_Measured.E(),iIsPrompt_Measured);
			Fill1DHisto(dHist_MissingMassSquared_WrongCombo,locMissingMassSquared,iIsPrompt_Measured);
			Fill1DHisto(dHist_PKMinus_InvMass_WrongCombo,(locKMinusP4_Measured+locProtonP4_Measured).M(),iIsPrompt_Measured);
			Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_WrongCombo,locDeltaTRF_Measured,(locKMinusP4_Measured+locProtonP4_Measured).M(),1);
		}
		
		/************************************ HISTOGRAM KINFIT ************************************/
		//Evaluate iIsPrompt_Kinfit
		Double_t dRFTime_Kinfit = dComboWrapper->Get_RFTime();
		Double_t locDeltaTRF_Kinfit = locBeamX4.T() - (dRFTime_Kinfit + (locBeamX4.Z() - dTargetCenterZ)/29.9792458);
		dHist_BeamDeltaTRF_Kinfit->Fill(locDeltaTRF_Kinfit,0);
		Double_t iIsPrompt_Kinfit=0;
		if((locDeltaTRF_Kinfit > -18.036) && (locDeltaTRF_Kinfit < -2.004)){
			iIsPrompt_Kinfit = 2;
			dHist_BeamDeltaTRF_Kinfit->Fill(locDeltaTRF_Kinfit,1,-0.125);
		}
		else if((locDeltaTRF_Kinfit >= -2.004) && (locDeltaTRF_Kinfit <= 2.004)){
			iIsPrompt_Kinfit = 1;
			dHist_BeamDeltaTRF_Kinfit->Fill(locDeltaTRF_Kinfit,1);
		}
		else if((locDeltaTRF_Kinfit > 2.004) && (locDeltaTRF_Kinfit < 18.036)){
			iIsPrompt_Kinfit = 2;
			dHist_BeamDeltaTRF_Kinfit->Fill(locDeltaTRF_Kinfit,1,-0.125);
		}
		else
			continue;
// 		iIsPrompt_Kinfit=1; //TODO deactivate accidental subtraction
		
		//Beam Energy
		if(locUsedSoFar_BeamEnergy_Kinfit.find(locBeamID) == locUsedSoFar_BeamEnergy_Kinfit.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			Fill1DHisto(dHist_BeamEnergy_Kinfit,locBeamP4.E(),iIsPrompt_Kinfit);
			locUsedSoFar_BeamEnergy_Kinfit.insert(locBeamID);
		}

		//Missing Mass Squared
		TLorentzVector locMissingP4_Kinfit = locBeamP4 + dTargetP4;
		locMissingP4_Kinfit -= locKPlusP4+ locKMinusP4 + locProtonP4;
		double locMissingMassSquared_Kinfit = locMissingP4_Kinfit.M2();
		
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass_Kinfit;
		locUsedThisCombo_MissingMass_Kinfit[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass_Kinfit[KPlus].insert(locKPlusTrackID);
		locUsedThisCombo_MissingMass_Kinfit[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_MissingMass_Kinfit[Proton].insert(locProtonTrackID);
		
		if(locUsedSoFar_MissingMass_Kinfit.find(locUsedThisCombo_MissingMass_Kinfit) == locUsedSoFar_MissingMass_Kinfit.end())
		{
			Fill1DHisto(dHist_MissingMassSquared_Kinfit,locMissingMassSquared_Kinfit,iIsPrompt_Kinfit);
			locUsedSoFar_MissingMass_Kinfit.insert(locUsedThisCombo_MissingMass_Kinfit);
		}
		
		//Invariant Mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_PKMinus_InvMass_Kinfit;
		locUsedThisCombo_PKMinus_InvMass_Kinfit[Unknown].insert(locBeamID); //beam for accidental subtraction
		locUsedThisCombo_PKMinus_InvMass_Kinfit[KMinus].insert(locKMinusTrackID);
		locUsedThisCombo_PKMinus_InvMass_Kinfit[Proton].insert(locProtonTrackID);
		
		if(locUsedSoFar_PKMinus_InvMass_Kinfit.find(locUsedThisCombo_PKMinus_InvMass_Kinfit) == locUsedSoFar_PKMinus_InvMass_Kinfit.end())
		{
			Fill1DHisto(dHist_PKMinus_InvMass_Kinfit,(locKMinusP4+locProtonP4).M(),iIsPrompt_Kinfit);
			Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_Kinfit,locDeltaTRF_Kinfit,(locKMinusP4+locProtonP4).M(),1);
			locUsedSoFar_PKMinus_InvMass_Kinfit.insert(locUsedThisCombo_PKMinus_InvMass_Kinfit);
		}
		
		/************************** CALCULATE WEIGHTS ********************************************/
		if( false ){ // don't place any cuts for now
			dComboWrapper->Set_IsComboCut(true);
		}
		else{
			locCombosPassedCuts.push_back(loc_i);

			map<Particle_t, set<Int_t> > AdditionalComboMap; // only for non beam particles
			AdditionalComboMap[KPlus].insert(locKPlusTrackID);
			AdditionalComboMap[KMinus].insert(locKMinusTrackID);
			AdditionalComboMap[Proton].insert(locProtonTrackID);
			if(AdditionalCombos.find(AdditionalComboMap) == AdditionalCombos.end()){
				MultComboWeight++;
				AdditionalCombos.insert(AdditionalComboMap);
			}
		}
		
	
	} // end of combo loop
	
	dHist_MultComboWeight->Fill(MultComboWeight);
	
	/***************************** SECOND LOOP *********************************************/
	//Loop over combos
	for(UInt_t loc_i:locCombosPassedCuts)
	{
		
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		
		/*********************************************** GET FOUR-MOMENTUM **********************************************/
		
		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locKPlusP4 = dKPlusWrapper->Get_P4();
		TLorentzVector locKMinusP4 = dKMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locBeamX4 = dComboBeamWrapper->Get_X4();
		TLorentzVector locKPlusX4 = dKPlusWrapper->Get_X4();
		TLorentzVector locKMinusX4 = dKMinusWrapper->Get_X4();
		TLorentzVector locProtonX4 = dProtonWrapper->Get_X4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locKPlusP4_Measured = dKPlusWrapper->Get_P4_Measured();
		TLorentzVector locKMinusP4_Measured = dKMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		TLorentzVector locKPlusX4_Measured = dKPlusWrapper->Get_X4_Measured();
		TLorentzVector locKMinusX4_Measured = dKMinusWrapper->Get_X4_Measured();
		TLorentzVector locProtonX4_Measured = dProtonWrapper->Get_X4_Measured();
		

		/************************************ HISTOGRAM WEIGHTED MEASURED ************************************/

		//Evaluate iIsPrompt_Measured
		Double_t dRFTime_Measured = dComboWrapper->Get_RFTime_Measured();
		Double_t dTargetCenterZ = dComboWrapper->Get_TargetCenter().Z();
		Double_t locDeltaTRF_Measured = locBeamX4_Measured.T() - (dRFTime_Measured + (locBeamX4_Measured.Z() - dTargetCenterZ)/29.9792458);
		dHist_BeamDeltaTRF_Measured_Weighted->Fill(locDeltaTRF_Measured,0);
		Double_t iIsPrompt_Measured=0;
		if((locDeltaTRF_Measured > -18.036) && (locDeltaTRF_Measured < -2.004)){
			iIsPrompt_Measured = 2;
			dHist_BeamDeltaTRF_Measured_Weighted->Fill(locDeltaTRF_Measured,1,-0.125);
		}
		else if((locDeltaTRF_Measured >= -2.004) && (locDeltaTRF_Measured <= 2.004)){
			iIsPrompt_Measured = 1;
			dHist_BeamDeltaTRF_Measured_Weighted->Fill(locDeltaTRF_Measured,1);
		}
		else if((locDeltaTRF_Measured > 2.004) && (locDeltaTRF_Measured < 18.036)){
			iIsPrompt_Measured = 2;
			dHist_BeamDeltaTRF_Measured_Weighted->Fill(locDeltaTRF_Measured,1,-0.125);
		}
		else
			continue;
// 		iIsPrompt_Measured=1; //TODO deactivate accidental subtraction
		
		
		//Beam Energy
		Fill1DHisto(dHist_BeamEnergy_Measured_Weighted,locBeamP4_Measured.E(),1./MultComboWeight,iIsPrompt_Measured);
		
		//Missing Mass Squared
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locKPlusP4_Measured + locKMinusP4_Measured + locProtonP4_Measured;
		double locMissingMassSquared = locMissingP4_Measured.M2();
		Fill1DHisto(dHist_MissingMassSquared_Measured_Weighted,locMissingMassSquared,1./MultComboWeight,iIsPrompt_Measured);
		
		//Invariant Mass
		Fill1DHisto(dHist_PKMinus_InvMass_Measured_Weighted,(locKMinusP4_Measured+locProtonP4_Measured).M(),1./MultComboWeight,iIsPrompt_Measured);
		Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_Measured_Weighted,locDeltaTRF_Measured,(locKMinusP4_Measured+locProtonP4_Measured).M(),1./MultComboWeight,1);
		
		
		/************************************ HISTOGRAM WEIGHTED KINFIT ************************************/
		//Evaluate iIsPrompt_Measured
		Double_t dRFTime_Kinfit = dComboWrapper->Get_RFTime();
		Double_t locDeltaTRF_Kinfit = locBeamX4.T() - (dRFTime_Kinfit + (locBeamX4.Z() - dTargetCenterZ)/29.9792458);
		dHist_BeamDeltaTRF_Kinfit_Weighted->Fill(locDeltaTRF_Kinfit,0);
		Double_t iIsPrompt_Kinfit=0;
		if((locDeltaTRF_Kinfit > -18.036) && (locDeltaTRF_Kinfit < -2.004)){
			iIsPrompt_Kinfit = 2;
			dHist_BeamDeltaTRF_Kinfit_Weighted->Fill(locDeltaTRF_Kinfit,1,-0.125);
		}
		else if((locDeltaTRF_Kinfit >= -2.004) && (locDeltaTRF_Kinfit <= 2.004)){
			iIsPrompt_Kinfit = 1;
			dHist_BeamDeltaTRF_Kinfit_Weighted->Fill(locDeltaTRF_Kinfit,1);
		}
		else if((locDeltaTRF_Kinfit > 2.004) && (locDeltaTRF_Kinfit < 18.036)){
			iIsPrompt_Kinfit = 2;
			dHist_BeamDeltaTRF_Kinfit_Weighted->Fill(locDeltaTRF_Kinfit,1,-0.125);
		}
		else
			continue;
// 		iIsPrompt_Kinfit=1; //TODO deactivate accidental subtraction
		
		//Beam Energy
		Fill1DHisto(dHist_BeamEnergy_Kinfit_Weighted,locBeamP4.E(),1./MultComboWeight,iIsPrompt_Kinfit);
		
		//Missing Mass Squared
		TLorentzVector locMissingP4_Kinfit = locBeamP4 + dTargetP4;
		locMissingP4_Kinfit -= locKPlusP4+ locKMinusP4 + locProtonP4;
		double locMissingMassSquared_Kinfit = locMissingP4_Kinfit.M2();
		Fill1DHisto(dHist_MissingMassSquared_Kinfit_Weighted,locMissingMassSquared_Kinfit,1./MultComboWeight,iIsPrompt_Kinfit);
		
		//Invariant Mass
		Fill1DHisto(dHist_PKMinus_InvMass_Kinfit_Weighted,(locKMinusP4+locProtonP4).M(),1./MultComboWeight,iIsPrompt_Kinfit);
		Fill2DHisto(dHist_BeamDeltaTRF_PKMinus_InvMass_Kinfit_Weighted,locDeltaTRF_Kinfit,(locKMinusP4+locProtonP4).M(),1./MultComboWeight,1);
		
	}
	
	/******************************************* HISTOGRAM THROWN ***************************************/
	
	if(dThrownBeam != NULL){
		
		TLorentzVector genBeamP4 = dThrownBeam->Get_P4();
			TLorentzVector genBeamX4 = dThrownBeam->Get_X4();
			dThrownWrapper->Set_ArrayIndex(0);
			if(!(dThrownWrapper->Get_PID()==11)){//KPlus
				dThrownWrapper->Set_ArrayIndex(1);
				if(!(dThrownWrapper->Get_PID()==11)){
					dThrownWrapper->Set_ArrayIndex(2);
				}
			}
			TLorentzVector genKPlusP4 = dThrownWrapper->Get_P4();
			TLorentzVector genKPlusX4 = dThrownWrapper->Get_X4();
			
			dThrownWrapper->Set_ArrayIndex(0);
			if(!(dThrownWrapper->Get_PID()==12)){//KMinus
				dThrownWrapper->Set_ArrayIndex(1);
				if(!(dThrownWrapper->Get_PID()==12)){
					dThrownWrapper->Set_ArrayIndex(2);
				}
			}
			TLorentzVector genKMinusP4 = dThrownWrapper->Get_P4();
			TLorentzVector genKMinusX4 = dThrownWrapper->Get_X4();
			
			dThrownWrapper->Set_ArrayIndex(0);
			if(!(dThrownWrapper->Get_PID()==14)){//Proton
				dThrownWrapper->Set_ArrayIndex(1);
				if(!(dThrownWrapper->Get_PID()==14)){
					dThrownWrapper->Set_ArrayIndex(2);
				}
			}
			TLorentzVector genProtonP4 = dThrownWrapper->Get_P4();
			TLorentzVector genProtonX4 = dThrownWrapper->Get_X4();
			
			TLorentzVector locMissingP4_Thrown = genBeamP4 + dTargetP4;
			locMissingP4_Thrown -= genKPlusP4 + genKMinusP4 + genProtonP4;
			
			
			dHist_BeamEnergy_Thrown->Fill(genBeamP4.E());
			dHist_MissingMassSquared_Thrown->Fill(locMissingP4_Thrown.M2());
			dHist_PKMinus_InvMass_Thrown->Fill((genKMinusP4+genProtonP4).M());
	}
	
	
	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	return kTRUE;
}

void DSelector_combostudy::Finalize(void)
{
	//Save anything to output here that you do not want to be in the default DSelector output ROOT file.

	//Otherwise, don't do anything else (especially if you are using PROOF).
		//If you are using PROOF, this function is called on each thread,
		//so anything you do will not have the combined information from the various threads.
		//Besides, it is best-practice to do post-processing (e.g. fitting) separately, in case there is a problem.

	//DO YOUR STUFF HERE

	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
