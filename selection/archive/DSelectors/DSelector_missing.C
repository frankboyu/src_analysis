#include "DSelector_missing.h"

void DSelector_missing::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	//dOutputFileName = "root/D_1_14__8_9_14_m13.root"; //"" not working
	dOutputFileName = "Rhocandidate.root"; //"" not working
	dOutputTreeFileName = ""; //"" for none
	dFlatTreeFileName = ""; //output flat tree (one combo per tree entry), "" for none
	dFlatTreeName = ""; //if blank, default name will be chosen
	//dSaveDefaultFlatBranches = true; // False: don't save default branches, reduce disk footprint.
	//dSaveTLorentzVectorsAsFundamentaFlatTree = false; // Default (or false): save particles as TLorentzVector objects. True: save as four doubles instead.

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
	std::deque<Particle_t> MyPhi;
	MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data

	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//PIDFOM (for charged tracks)
	dAnalysisActions.push_back(new DHistogramAction_PIDFOM(dComboWrapper));
	//dAnalysisActions.push_back(new DCutAction_PIDFOM(dComboWrapper, KPlus, 0.1));
	//dAnalysisActions.push_back(new DCutAction_EachPIDFOM(dComboWrapper, 0.1));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//CUT ON SHOWER QUALITY
	//dAnalysisActions.push_back(new DCutAction_ShowerQuality(dComboWrapper, SYS_FCAL, 0.5));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 6.5, 13));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));
	// dAnalysisActions.push_back(new DCutAction_KinFitFOM(dComboWrapper, 0.02, ""));


	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	//dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()
	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	//Beam and Vertex 
       dHist_BeamDeltaT = new TH1F("BeamDeltaT", "; t_{Tagger} - t_{RF} (ns)", 600, -20., 20.);
       dHist_BeamDeltaT_After = new TH1F("BeamDeltaT_After", "; t_{Tagger} - t_{RF} (ns)", 600, -20., 20.);
       dHist_BeamDeltaT_After_Weight = new TH1F("dHist_BeamDeltaT_After_Weight", "; t_{Tagger} - t_{RF} (ns)", 600, -20., 20.);
       dHist_vertex = new TH1F("VertexZ", "; Vertex Z (cm); # Combos / 0.1 cm", 1000, 30, 130);
       dHist_vertex_cut = new TH1F("VertexZ_cut", "; Vertex Z (cm); # Combos / 0.1 cm", 1000, 30, 130);
       dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	
	 
    // dHist_Bunch = new TH1F("BeamBunch", "; t_{Tagger} - t_{RF} (ns)", 600, -20., 20.);
     


	
//1.1  Missing Mass for different combination for measured value without any cut
	 dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -1, 4);
	 dHist_MissingMass = new TH1I("MissingMass", ";Missing Mass  (GeV/)^", 600, 0.5, 2);
	 dHist_MissingMass_piplus_piminus = new TH1I("MissingMass_piplus_piminus", ";Missing Mass  (GeV/)^", 600, 1.7, 4);
	 dHist_MissingMass_proton = new TH1I("MissingMass_proton", ";Missing Mass  (GeV/)^", 600, 1, 4);
 //1.2 Missing Mass for different Combination using Kinfit Value without  any cut
	dHist_MissingMassSquaredkinfit = new TH1I("MissingMassSquaredkinfit", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, 0.7, 0.9);
	dHist_MissingMasskinfit = new TH1I("MissingMasskinfit", ";Missing Mass  (GeV/)^", 600, 0.86, 9);
	dHist_MissingMass_piplus_piminuskinfit = new TH1I("MissingMass_piplus_piminuskinfit", ";Missing Mass  (GeV/)^", 600, 1.7, 4);
	dHist_MissingMass_protonkinfit = new TH1I("MissingMass_protonkinfit", ";Missing Mass  (GeV/)^", 600, 1.3, 4);
//1.3  Invariant Mass for measured values without cut.
	 dHist_PiPlus_Proton_Mass = new TH1F("PiPlus_Proton", ";#PiPlusProton Mass(GeV)" , 600 , 1 ,2);
      dHist_PiPlus_PiMinus_Mass =  new TH1F("PiPlus_PiMinus", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
	 dHist_PiMinus_Proton_Mass = new TH1F("PiMinus_Proton", ";#PiMinusProton Mass(GeV)" , 600 , 1 ,3);
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
	 dHist_PiPlus_Proton_vs_PiMinus_Proton_mass = new  TH2F("PiPlus_protonvs_Piminus_Proton",";#PiPlus_Proton(GeV); PiMinus_Proton Mass(GeV)" , 100 , 1 ,3 , 100 , 1 ,3);	           
	 dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass = new TH2F("PiPlus_PiMinus_vs_PiMinus_Proton",";#PiPlus_PiMinus(GeV); PiMinusProton Mass(GeV)" , 100 , 0.2 ,2 , 100 , 1 ,3);
	        //Kinhistogram
	 dHist_PiPlus_Proton_Mass_kin = new TH1F("PiPlus_Proton_kin", ";#PiPlusProton Mass(GeV)" , 600 , 1 ,2);
      dHist_PiPlus_PiMinus_Mass_kin =  new TH1F("PiPlus_PiMinus_kin", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
	 dHist_PiMinus_Proton_Mass_kin = new TH1F("PiMinus_Proton_kin", ";#PiMinusProton Mass(GeV)" , 600 , 1 ,3);
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_kin = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus_kin",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
	 dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_kin = new  TH2F("PiPlus_protonvs_Piminus_Proton_kin",";#PiPlus_Proton(GeV); PiMinus_Proton Mass(GeV)" , 100 , 1 ,3 , 100 , 1 ,3);	           
	 dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_kin = new TH2F("PiPlus_PiMinus_vs_PiMinus_Proton_kin",";#PiPlus_PiMinus(GeV); PiMinusProton Mass(GeV)" , 100 , 0.2 ,2 , 100 , 1 ,3);
	


   ///////////////Cut on Kinematic Confidence Level:
	 dHist_KinFitCL = new TH1F("KinFit_CL", ";Kinematic Fit Confidence Level", 1200, 0., 1.);
	 dHist_KinFitCL_cut = new TH1F("KinFit_CL_cut", ";Kinematic Fit Confidence Level", 1200, 0., 1.);

//2.1 Missing mass and invariant mass after Confidence Level for both Measured and Kinfit Vaue;
	 dHist_MissingMassSquared_CL_cut = new TH1I("MissingMassSquared_CL", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -1, 4);
	 dHist_MissingMass_CL_cut = new TH1I("MissingMass_CL", ";Missing Mass  (GeV/)^", 600, 0.5, 2);					
	 dHist_MissingMassSquaredkinfit_CL_cut = new TH1I("MissingMassSquaredkinfit_CL", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, 0.7, 0.9);
	 dHist_MissingMasskinfit_CL_cut = new TH1I("MissingMasskinfit_CL", ";Missing Mass  (GeV/)^", 600, 0.86, 9);		          
      dHist_PiPlus_PiMinus_Mass_CL_cut = new TH1F("PiPlus_PiMinus_CL_cut", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_CL_cut = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus_CL_cut",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
	    //kin
	 dHist_PiPlus_PiMinus_Mass_CL_cut_kin = new TH1F("PiPlus_PiMinus_CL_cut_kin", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_CL_cut_kin = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus_CL_cut_kin",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
		   	   
//3.1 Missing Mass for different values after Missing Mass Cut
	
	dHist_MissingMassSquared_cut = new TH1I("MissingMassSquared_cut", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -1, 4);
	dHist_MissingMass_cut = new TH1I("MissingMass_cut", ";Missing Mass  (GeV/)^", 600, 0.5, 2);
	dHist_MissingMass_piplus_piminus_cut = new TH1I("MissingMass_piplus_piminus_cut", ";Missing Mass  (GeV/)^", 600, 1.7, 4);
	dHist_MissingMass_proton_cut = new TH1I("MissingMass_proton_cut", ";Missing Mass  (GeV/)^", 600, 1, 4);
			
	dHist_MissingMassSquaredkinfit_cut = new TH1I("MissingMassSquaredkinfit_cut", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, 0.7, 0.9);
	dHist_MissingMasskinfit_cut = new TH1I("MissingMasskinfit_cut", ";Missing Mass  (GeV/)^", 600, 0.86, 0.9);
	dHist_MissingMass_piplus_piminuskinfit_cut = new TH1I("MissingMass_piplus_piminuskinfit_cut", ";Missing Mass  (GeV/)^", 600, 1.7, 4);
	dHist_MissingMass_protonkinfit_cut = new TH1I("MissingMass_protonkinfit_cut", ";Missing Mass  (GeV/)^", 600, 1.3, 4);
// Invariant Mass for Measured value after Cut
   	dHist_PiPlus_Proton_Mass_cut = new TH1F("PiPlus_Proton_cut", ";#PiPlusProton Mass(GeV)" , 600 , 1 ,2);
     dHist_PiPlus_PiMinus_Mass_cut =  new TH1F("PiPlus_PiMinus_cut", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
	 dHist_PiMinus_Proton_Mass_cut = new TH1F("PiMinus_Proton_cut", ";#PiMinusProton Mass(GeV)" , 600 , 1 ,3);
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_cut = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus_cut",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
	 dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_cut = new  TH2F("PiPlus_protonvs_Piminus_Proton_cut",";#PiPlus_Proton(GeV); PiMinus_Proton Mass(GeV)" , 100 , 1 ,3 , 100 , 1 ,3);	           
	 dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_cut = new TH2F("PiPlus_PiMinus_vs_PiMinus_Proton_cut",";#PiPlus_PiMinus(GeV); PiMinusProton Mass(GeV)" , 100 , 0.2 ,2 , 100 , 1 ,3);
		
		//kin
	 dHist_PiPlus_Proton_Mass_cut_kin = new TH1F("PiPlus_Proton_cut_kin", ";#PiPlusProton Mass(GeV)" , 600 , 1 ,2);
     dHist_PiPlus_PiMinus_Mass_cut_kin =  new TH1F("PiPlus_PiMinus_cut_kin", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5); 
	 dHist_PiMinus_Proton_Mass_cut_kin = new TH1F("PiMinus_Proton_cut_kin", ";#PiMinusProton Mass(GeV)" , 600 , 1 ,3);
      dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_cut_kin = new TH2F("PiPlus_Proton_vs_PiPlus_PiMinus_cut_kin",";#PiPlus_Proton Mass(GeV);PiPlus_PiMinus Mass(GeV)" , 100 , 1 ,3 , 100 , 0.25 ,2.5);
	 dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_cut_kin = new  TH2F("PiPlus_protonvs_Piminus_Proton_cut_kin",";#PiPlus_Proton(GeV); PiMinus_Proton Mass(GeV)" , 100 , 1 ,3 , 100 , 1 ,3);	           
	 dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_cut_kin = new TH2F("PiPlus_PiMinus_vs_PiMinus_Proton_cut_kin",";#PiPlus_PiMinus(GeV); PiMinusProton Mass(GeV)" , 100 , 0.2 ,2 , 100 , 1 ,3);
		

			   ////t and u distribution
     dHist_minus_t = new TH1F("t_dist",";-t (GeV/c)^2",40,0,20.0);
	dHist_minus_t1 = new TH1F("t1_dist",";-t (GeV/c)^2",40,0,20.0);
	dHist_minus_u = new TH1F("u_dist",";-u (GeV/c)^2",40,-2.0,18.0);
	dHist_minus_u1 = new TH1F("u1_dist",";-u (GeV/c)^2",40,-2.0,18.0);   

	       //////t and u cut 
     dHist_minus_tcut = new TH1F("t_distcut",";-t (GeV/c)^2",40,0,20.0);	     
     dHist_minus_ucut = new TH1F("u_distcut",";-u (GeV/c)^2",40,-2.0,18.0);	       
     dHist_pippim =  new TH1F("pippim", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5);   
	       ///rho mass cut
     dHist_minus_tcut2= new TH1F("t_distcut2",";-t (GeV/c)^2",40,0,20.0);
  	dHist_minus_ucut2 = new TH1F("u_distcut2",";-u (GeV/c)^2",40,-2.0,18.0);
	dHist_pippim2 =  new TH1F("pippim2", ";#PiPlus_PiMinus(GeV)" , 600 , 0.1 ,2.5);   








	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - MAIN TREE *************************/

	//EXAMPLE MAIN TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dTreeInterface->Create_Branch_Fundamental<Int_t>("my_int"); //fundamental = char, int, float, double, etc.
	dTreeInterface->Create_Branch_FundamentalArray<Int_t>("my_int_array", "my_int");
	dTreeInterface->Create_Branch_FundamentalArray<Float_t>("my_combo_array", "NumCombos");
	dTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("my_p4");
	dTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("my_p4_array");
	*/

	/************************** EXAMPLE USER INITIALIZATION: CUSTOM OUTPUT BRANCHES - FLAT TREE *************************/

	// RECOMMENDED: CREATE ACCIDENTAL WEIGHT BRANCH
	// dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("accidweight");

	//EXAMPLE FLAT TREE CUSTOM BRANCHES (OUTPUT ROOT FILE NAME MUST FIRST BE GIVEN!!!! (ABOVE: TOP)):
	//The type for the branch must be included in the brackets
	//1st function argument is the name of the branch
	//2nd function argument is the name of the branch that contains the size of the array (for fundamentals only)
	/*
	dFlatTreeInterface->Create_Branch_Fundamental<Int_t>("flat_my_int"); //fundamental = char, int, float, double, etc.
	dFlatTreeInterface->Create_Branch_FundamentalArray<Int_t>("flat_my_int_array", "flat_my_int");
	dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("flat_my_p4");
	dFlatTreeInterface->Create_Branch_ClonesArray<TLorentzVector>("flat_my_p4_array");
	*/

	/************************************* ADVANCED EXAMPLE: CHOOSE BRANCHES TO READ ************************************/

	//TO SAVE PROCESSING TIME
		//If you know you don't need all of the branches/data, but just a subset of it, you can speed things up
		//By default, for each event, the data is retrieved for all branches
		//If you know you only need data for some branches, you can skip grabbing data from the branches you don't need
		//Do this by doing something similar to the commented code below

	//dTreeInterface->Clear_GetEntryBranches(); //now get none
	//dTreeInterface->Register_GetEntryBranch("Proton__P4"); //manually set the branches you want

	/************************************** DETERMINE IF ANALYZING SIMULATED DATA *************************************/

	//dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_missing::Process(Long64_t locEntry)
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
/*	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}
*/
	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	Reset_Actions_NewEvent();
	//dAnalyzeCutActions->Reset_NewEvent(); // manual action, must call Reset_NewEvent()

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
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
	
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass_cut;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_CLcut;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_tdistribution;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_tcutdistribution;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_rhocutdistribution;

	//INSERT USER ANALYSIS UNIQUENESS TRACKING HERE

	/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

	/*
	Int_t locMyInt = 7;
	dTreeInterface->Fill_Fundamental<Int_t>("my_int", locMyInt);

	TLorentzVector locMyP4(4.0, 3.0, 2.0, 1.0);
	dTreeInterface->Fill_TObject<TLorentzVector>("my_p4", locMyP4);

	for(int loc_i = 0; loc_i < locMyInt; ++loc_i)
		dTreeInterface->Fill_Fundamental<Int_t>("my_int_array", 3*loc_i, loc_i); //2nd argument = value, 3rd = array index
	*/

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
		Int_t locPiPlusTrackID = dPiPlusWrapper->Get_TrackID();
		Int_t locPiMinusTrackID = dPiMinusWrapper->Get_TrackID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locPiPlusP4 = dPiPlusWrapper->Get_P4();
		TLorentzVector locPiMinusP4 = dPiMinusWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		TLorentzVector locMissingNeutronP4 = dMissingNeutronWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locPiPlusP4_Measured = dPiPlusWrapper->Get_P4_Measured();
		TLorentzVector locPiMinusP4_Measured = dPiMinusWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();


		/********************************************* GET COMBO RF TIMING INFO *****************************************/

		TLorentzVector locBeamX4_Measured = dComboBeamWrapper->Get_X4_Measured();
		Double_t locBunchPeriod = dAnalysisUtilities.Get_BeamBunchPeriod(Get_RunNumber());
	//	dHist_Bunch->Fill(locBunchPeriod);
		 Double_t locDeltaT_RF = dAnalysisUtilities.Get_DeltaT_RF(Get_RunNumber(), locBeamX4_Measured, dComboWrapper);
		  dHist_BeamDeltaT->Fill(locDeltaT_RF);
		 Int_t locRelBeamBucket = dAnalysisUtilities.Get_RelativeBeamBucket(Get_RunNumber(), locBeamX4_Measured, dComboWrapper); // 0 for in-time events, non-zero integer for out-of-time photons
		 Int_t locNumOutOfTimeBunchesInTree = 4; //YOU need to specify this number
			//Number of out-of-time beam bunches in tree (on a single side, so that total number out-of-time bunches accepted is 2 times this number for left + right bunches) 

		 Bool_t locSkipNearestOutOfTimeBunch = true; // True: skip events from nearest out-of-time bunch on either side (recommended).
		 Int_t locNumOutOfTimeBunchesToUse = locSkipNearestOutOfTimeBunch ? locNumOutOfTimeBunchesInTree-1:locNumOutOfTimeBunchesInTree; 
		 Double_t locAccidentalScalingFactor = dAnalysisUtilities.Get_AccidentalScalingFactor(Get_RunNumber(), locBeamP4.E(), dIsMC); // Ideal value would be 1, but deviations require added factor, which is different for data and MC.
		 Double_t locAccidentalScalingFactorError = dAnalysisUtilities.Get_AccidentalScalingFactorError(Get_RunNumber(), locBeamP4.E()); // Ideal value would be 1, but deviations observed, need added factor.
		 Double_t locAccWeight = locRelBeamBucket==0 ? 1 : -locAccidentalScalingFactor/(2*locNumOutOfTimeBunchesToUse) ; // Weight by 1 for in-time events, ScalingFactor*(1/NBunches) for out-of-time
		 if(locSkipNearestOutOfTimeBunch && abs(locRelBeamBucket)==1) { // Skip nearest out-of-time bunch: tails of in-time distribution also leak in
		 	dComboWrapper->Set_IsComboCut(true); 
		 	continue; 
		 } 

         dHist_BeamDeltaT_After->Fill(locDeltaT_RF);
         dHist_BeamDeltaT_After_Weight->Fill(locDeltaT_RF,locAccWeight);

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locPiPlusP4_Measured + locPiMinusP4_Measured + locProtonP4_Measured;
		
		
		TLorentzVector locPiPlus_Proton_P4_Meas = locProtonP4_Measured+locPiPlusP4_Measured;
		TLorentzVector locPiPlus_PiMinus_P4_Meas = locPiPlusP4_Measured +locPiMinusP4_Measured;
		TLorentzVector locPiMinus_Proton_P4_Meas = locPiMinusP4_Measured +locProtonP4_Measured;


///////Kinfit
		TLorentzVector locPiPlus_Proton_P4 = locProtonP4+ locPiPlusP4;
		TLorentzVector locPiPlus_PiMinus_P4 = locPiPlusP4 +locPiMinusP4;
		TLorentzVector locPiMinus_Proton_P4 = locPiMinusP4 +locProtonP4;
		
		
		
		
		 double locPiPlus_Proton_Mass_Meas = locPiPlus_Proton_P4_Meas.M();
		 double locPiPlus_PiMinus_Mass_Meas = locPiPlus_PiMinus_P4_Meas.M();
		 double locPiMinus_Proton_Mass_Meas = locPiMinus_Proton_P4_Meas.M();
               


		 double locPiPlus_Proton_Mass = locPiPlus_Proton_P4.M();
		 double locPiPlus_PiMinus_Mass = locPiPlus_PiMinus_P4.M();
		 double locPiMinus_Proton_Mass = locPiMinus_Proton_P4.M();
               



		//********************************************************////////////////////////
		TLorentzVector locMissingP4_Measured_piplus_piminus = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured_piplus_piminus -=  locPiPlusP4_Measured + locPiMinusP4_Measured;
		
		TLorentzVector locMissingP4_Measured_proton = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured_proton -= locProtonP4_Measured; 
		
		
		//Kinfit
		TLorentzVector locMissingP4 = locBeamP4 + dTargetP4;
		locMissingP4 -= locPiPlusP4 + locPiMinusP4 + locProtonP4;   
		
	        TLorentzVector locMissingP4_piplus_piminus = locBeamP4 + dTargetP4;
		locMissingP4_piplus_piminus -=  locPiPlusP4 + locPiMinusP4;
		
		TLorentzVector locMissingP4_proton = locBeamP4 + dTargetP4;
		locMissingP4_proton -= locProtonP4; 


		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/

		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;

		//if you manually execute any actions, and it fails a cut, be sure to call:
			//dComboWrapper->Set_IsComboCut(true);

		/**************************************** EXAMPLE: FILL CUSTOM OUTPUT BRANCHES **************************************/

		/*
		TLorentzVector locMyComboP4(8.0, 7.0, 6.0, 5.0);
		//for arrays below: 2nd argument is value, 3rd is array index
		//NOTE: By filling here, AFTER the cuts above, some indices won't be updated (and will be whatever they were from the last event)
			//So, when you draw the branch, be sure to cut on "IsComboCut" to avoid these.
		dTreeInterface->Fill_Fundamental<Float_t>("my_combo_array", -2*loc_i, loc_i);
		dTreeInterface->Fill_TObject<TLorentzVector>("my_p4_array", locMyComboP4, loc_i);
		*/
          
           /********************************************Cut 1: Cut on Beam Energy. ******************************************************/
		    
                if (locBeamP4.E() < 7.5)
                {
                	dComboWrapper->Set_IsComboCut(true);
                	continue;
                }
                  
             
		/**************************************** EXAMPLE: HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			//dHist_BeamEnergy->Fill(locBeamP4.E()); // Fills in-time and out-of-time beam photon combos
			dHist_BeamEnergy->Fill(locBeamP4.E(),locAccWeight); // Alternate version with accidental subtraction
              
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}
              
              dHist_vertex->Fill(locBeamX4_Measured.Z());
		    if (fabs(locBeamX4_Measured.Z() - 65) > 13){
		     dComboWrapper->Set_IsComboCut(true);
		     continue;
		}
		dHist_vertex_cut->Fill(locBeamX4_Measured.Z());
		
		/************************************ EXAMPLE: HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();
		double  locMissingMass    = locMissingP4_Measured.M();
		double  locMissingEnergy = locMissingP4_Measured.E();                                                                                                                                         
       	double locMissingMass_proton      = locMissingP4_Measured_proton.M();

		double locMissingMass_piplus_piminus    =locMissingP4_Measured_piplus_piminus.M();
		double locMissingMass_piplus_piminuskinfit      =locMissingP4_piplus_piminus.M();
		double locMissingMass_protonkinfit      =locMissingP4_proton.M();
		double locMissingMassSquaredkinfit = locMissingP4.M2();
		double locMissingMasskinfit    = locMissingP4.M();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{

			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared,locAccWeight);
			dHist_MissingMass->Fill(locMissingMass,locAccWeight);
			dHist_MissingMass_piplus_piminus->Fill(locMissingMass_piplus_piminus ,locAccWeight);
			dHist_MissingMass_proton->Fill(locMissingMass_proton,locAccWeight);
			
			
			dHist_MissingMassSquaredkinfit->Fill(locMissingMassSquaredkinfit,locAccWeight);
			dHist_MissingMasskinfit->Fill(locMissingMasskinfit,locAccWeight); // Fills in-time and out-of-time beam photon combos
			dHist_MissingMass_piplus_piminuskinfit->Fill(locMissingMass_piplus_piminuskinfit ,locAccWeight);
			dHist_MissingMass_protonkinfit->Fill(locMissingMass_protonkinfit,locAccWeight);
			//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
			
                        dHist_PiPlus_Proton_Mass->Fill(locPiPlus_Proton_Mass_Meas,locAccWeight);
                        dHist_PiPlus_PiMinus_Mass->Fill(locPiPlus_PiMinus_Mass_Meas,locAccWeight);	        
		                dHist_PiMinus_Proton_Mass->Fill(locPiMinus_Proton_Mass_Meas,locAccWeight);
		          
		                dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass->Fill(locPiPlus_Proton_Mass_Meas, locPiPlus_PiMinus_Mass_Meas,locAccWeight);
		                dHist_PiPlus_Proton_vs_PiMinus_Proton_mass->Fill(locPiPlus_Proton_Mass_Meas, locPiMinus_Proton_Mass_Meas,locAccWeight);
		           
		                dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass->Fill(locPiPlus_PiMinus_Mass_Meas, locPiMinus_Proton_Mass_Meas,locAccWeight); 

		                 dHist_PiPlus_Proton_Mass_kin->Fill(locPiPlus_Proton_Mass,locAccWeight);
                        dHist_PiPlus_PiMinus_Mass_kin->Fill(locPiPlus_PiMinus_Mass,locAccWeight);	        
		                dHist_PiMinus_Proton_Mass_kin->Fill(locPiMinus_Proton_Mass,locAccWeight);
		          
		                dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_kin->Fill(locPiPlus_Proton_Mass, locPiPlus_PiMinus_Mass,locAccWeight);
		                dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_kin->Fill(locPiPlus_Proton_Mass, locPiMinus_Proton_Mass,locAccWeight);
		           
		                dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_kin->Fill(locPiPlus_PiMinus_Mass, locPiMinus_Proton_Mass,locAccWeight);


			//unique missing mass combo: histogram it, and register this combo of particles
			//dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
			
                        
                        
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//E.g. Cut
		//if((locMissingMassSquared < -0.04) || (locMissingMassSquared > 0.04))
		//{
		//	dComboWrapper->Set_IsComboCut(true);
		//	continue;
		//}
               double locKinFit_CL = dComboWrapper->Get_ConfidenceLevel_KinFit("");
		      		        dHist_KinFitCL->Fill(locKinFit_CL);

			             //	Cut locKinFit_CLon KinFit CL
			      	if(locKinFit_CL < 0.001)  //0.002
			      	{
			      		dComboWrapper->Set_IsComboCut(true);
			      		continue; 
			      	}

       dHist_KinFitCL_cut->Fill(locKinFit_CL);

		


        //   if((locMissingMass_piplus_piminus < 1.8) || (locMissingMass_piplus_piminus> 1.9))
	//	{
	  //     	dComboWrapper->Set_IsComboCut(true);
	//		continue;
	//	}

		map<Particle_t, set<Int_t> > locUsedThisCombo_CLcut;
		locUsedThisCombo_CLcut[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_CLcut[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_CLcut[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_CLcut[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_CLcut.find(locUsedThisCombo_CLcut) == locUsedSoFar_CLcut.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared_CL_cut->Fill(locMissingMassSquared,locAccWeight);
			dHist_MissingMass_CL_cut->Fill(locMissingMass,locAccWeight);
						
			
			dHist_MissingMassSquaredkinfit_CL_cut->Fill(locMissingMassSquaredkinfit,locAccWeight);
			dHist_MissingMasskinfit_CL_cut->Fill(locMissingMasskinfit,locAccWeight);// Fills in-time and out-of-time beam photon combos			          
               dHist_PiPlus_PiMinus_Mass_CL_cut->Fill(locPiPlus_PiMinus_Mass_Meas,locAccWeight);      
		    	          
		    dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_CL_cut->Fill(locPiPlus_Proton_Mass_Meas, locPiPlus_PiMinus_Mass_Meas,locAccWeight);
		   

		    dHist_PiPlus_PiMinus_Mass_CL_cut_kin->Fill(locPiPlus_PiMinus_Mass,locAccWeight);     		    	          
		    dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_CL_cut_kin->Fill(locPiPlus_Proton_Mass, locPiPlus_PiMinus_Mass,locAccWeight);
		   
		       
		     //dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_CL_->Fill(locPiPlus_PiMinus_Mass_Meas, locPiMinus_Proton_Mass_Meas,locAccWeight);
			locUsedSoFar_CLcut.insert(locUsedThisCombo_CLcut);
		}

                   // if((locMissingMass < 0.8) || (locMissingMass > 1.2))
		            if((locMissingMassSquared < 0.5) || (locMissingMassSquared > 1.2))
		   {
			      dComboWrapper->Set_IsComboCut(true);
			      continue;
		    }

               map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass_cut;
		locUsedThisCombo_MissingMass_cut[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass_cut[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_MissingMass_cut[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_MissingMass_cut[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass_cut.find(locUsedThisCombo_MissingMass_cut) == locUsedSoFar_MissingMass_cut.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared_cut->Fill(locMissingMassSquared,locAccWeight);
			dHist_MissingMass_cut->Fill(locMissingMass,locAccWeight);
			dHist_MissingMass_piplus_piminus_cut->Fill(locMissingMass_piplus_piminus,locAccWeight);
			dHist_MissingMass_proton_cut->Fill(locMissingMass_proton,locAccWeight);
			
			
			dHist_MissingMassSquaredkinfit_cut->Fill(locMissingMassSquaredkinfit,locAccWeight);
			dHist_MissingMasskinfit_cut->Fill(locMissingMasskinfit,locAccWeight); // Fills in-time and out-of-time beam photon combos
			dHist_MissingMass_piplus_piminuskinfit_cut->Fill(locMissingMass_piplus_piminuskinfit ,locAccWeight);
			dHist_MissingMass_protonkinfit_cut->Fill(locMissingMass_protonkinfit,locAccWeight); // Fills in-time and out-of-time beam photon combos
			
            dHist_PiPlus_Proton_Mass_cut->Fill(locPiPlus_Proton_Mass_Meas,locAccWeight);
            dHist_PiPlus_PiMinus_Mass_cut->Fill(locPiPlus_PiMinus_Mass_Meas,locAccWeight);	        
		  dHist_PiMinus_Proton_Mass_cut->Fill(locPiMinus_Proton_Mass_Meas,locAccWeight);		          
		  dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_cut->Fill(locPiPlus_Proton_Mass_Meas, locPiPlus_PiMinus_Mass_Meas,locAccWeight);
		  dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_cut->Fill(locPiPlus_Proton_Mass_Meas, locPiMinus_Proton_Mass_Meas,locAccWeight);		       
		  dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass->Fill(locPiPlus_PiMinus_Mass_Meas, locPiMinus_Proton_Mass_Meas,locAccWeight);
                    ///kinfit mass
		  dHist_PiPlus_Proton_Mass_cut_kin->Fill(locPiPlus_Proton_Mass,locAccWeight);
            dHist_PiPlus_PiMinus_Mass_cut_kin->Fill(locPiPlus_PiMinus_Mass,locAccWeight);	        
		  dHist_PiMinus_Proton_Mass_cut_kin->Fill(locPiMinus_Proton_Mass,locAccWeight);		          
		  dHist_PiPlus_Proton_vs_PiPlus_PiMinus_mass_cut_kin->Fill(locPiPlus_Proton_Mass, locPiPlus_PiMinus_Mass,locAccWeight);
		  dHist_PiPlus_Proton_vs_PiMinus_Proton_mass_cut_kin->Fill(locPiPlus_Proton_Mass, locPiMinus_Proton_Mass,locAccWeight);		       
		  dHist_PiPlus_PiMinus_vs_PiMinus_Proton_mass_kin->Fill(locPiPlus_PiMinus_Mass, locPiMinus_Proton_Mass,locAccWeight);

			locUsedSoFar_MissingMass_cut.insert(locUsedThisCombo_MissingMass_cut);
		}

		double minus_t = -((locBeamP4_Measured - locPiPlus_PiMinus_P4_Meas).M2());                
           double minus_t1 = -((locProtonP4_Measured -dTargetP4).M2());

           double minus_u = -((locBeamP4_Measured - locProtonP4_Measured).M2());
           double minus_u1 = -((locPiPlus_PiMinus_P4_Meas- dTargetP4).M2());


		
			map<Particle_t, set<Int_t> > locUsedThisCombo_tdistribution;
		locUsedThisCombo_tdistribution[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_tdistribution[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_tdistribution[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_tdistribution[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_tdistribution.find(locUsedThisCombo_tdistribution) == locUsedSoFar_tdistribution.end())
			//unique missing mass combo: histogram it, and register this combo of particles
			//dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
					        
                          
          {    dHist_minus_t->Fill(minus_t,locAccWeight);
               dHist_minus_t1->Fill(minus_t1,locAccWeight);
               dHist_minus_u->Fill(minus_u,locAccWeight);
               dHist_minus_u1->Fill(minus_u1,locAccWeight);                           
			  locUsedSoFar_tdistribution.insert(locUsedThisCombo_tdistribution);
		}

		// if((fabs(minus_t) < 1))
		         if((fabs(minus_t) < 1) || (fabs(minus_u) < 1))
	    {
	       	dComboWrapper->Set_IsComboCut(true);
			continue;
		}
   
			map<Particle_t, set<Int_t> > locUsedThisCombo_tcutdistribution;
		locUsedThisCombo_tcutdistribution[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_tcutdistribution[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_tcutdistribution[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_tcutdistribution[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_tcutdistribution.find(locUsedThisCombo_tcutdistribution) == locUsedSoFar_tcutdistribution.end())
			//unique missing mass combo: histogram it, and register this combo of particles
			//dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
					        
                          
          {    dHist_minus_tcut->Fill(minus_t,locAccWeight);                         
               dHist_minus_ucut->Fill(minus_u,locAccWeight);
               dHist_pippim->Fill(locPiPlus_PiMinus_Mass,locAccWeight); 
                         
                          
			  locUsedSoFar_tcutdistribution.insert(locUsedThisCombo_tcutdistribution);
		}

       
                    if((locPiPlus_Proton_Mass_Meas < 1.4) || (locPiMinus_Proton_Mass_Meas < 1.4))
		   {
			      dComboWrapper->Set_IsComboCut(true);
			      continue;
		    }



			map<Particle_t, set<Int_t> > locUsedThisCombo_rhocutdistribution;
		locUsedThisCombo_rhocutdistribution[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_rhocutdistribution[PiPlus].insert(locPiPlusTrackID);
		locUsedThisCombo_rhocutdistribution[PiMinus].insert(locPiMinusTrackID);
		locUsedThisCombo_rhocutdistribution[Proton].insert(locProtonTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_rhocutdistribution.find(locUsedThisCombo_rhocutdistribution) == locUsedSoFar_rhocutdistribution.end())
			//unique missing mass combo: histogram it, and register this combo of particles
			//dHist_MissingMassSquared->Fill(locMissingMassSquared); // Fills in-time and out-of-time beam photon combos
					        
                          
          {    dHist_minus_tcut2->Fill(minus_t,locAccWeight);                         
               dHist_minus_ucut2->Fill(minus_u,locAccWeight);
               dHist_pippim2->Fill(locPiPlus_PiMinus_Mass,locAccWeight);                     
                           
			  locUsedSoFar_rhocutdistribution.insert(locUsedThisCombo_rhocutdistribution);
		}
		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

		// RECOMMENDED: FILL ACCIDENTAL WEIGHT
		// dFlatTreeInterface->Fill_Fundamental<Double_t>("accidweight",locHistAccidWeightFactor);

		/*
		//FILL ANY CUSTOM BRANCHES FIRST!!
		Int_t locMyInt_Flat = 7;
		dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int", locMyInt_Flat);

		TLorentzVector locMyP4_Flat(4.0, 3.0, 2.0, 1.0);
		dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4", locMyP4_Flat);

		for(int loc_j = 0; loc_j < locMyInt_Flat; ++loc_j)
		{
			dFlatTreeInterface->Fill_Fundamental<Int_t>("flat_my_int_array", 3*loc_j, loc_j); //2nd argument = value, 3rd = array index
			TLorentzVector locMyComboP4_Flat(8.0, 7.0, 6.0, 5.0);
			dFlatTreeInterface->Fill_TObject<TLorentzVector>("flat_my_p4_array", locMyComboP4_Flat, loc_j);
		}
		*/

		//FILL FLAT TREE
		//Fill_FlatTree(); //for the active combo
	} // end of combo loop

	//FILL HISTOGRAMS: Num combos / events surviving actions
	Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/
/*
	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

		//Do stuff with the wrapper here ...
	}
*/
	/****************************************** LOOP OVER OTHER ARRAYS (OPTIONAL) ***************************************/
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

	/************************************ EXAMPLE: FILL CLONE OF TTREE HERE WITH CUTS APPLIED ************************************/
/*
	Bool_t locIsEventCut = true;
	for(UInt_t loc_i = 0; loc_i < Get_NumCombos(); ++loc_i) {
		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);
		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut())
			continue;
		locIsEventCut = false; // At least one combo succeeded
		break;
	}
	if(!locIsEventCut && dOutputTreeFileName != "")
		Fill_OutputTree();
*/

	return kTRUE;
}

void DSelector_missing::Finalize(void)
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
