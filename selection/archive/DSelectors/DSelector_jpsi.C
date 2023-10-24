#include "DSelector_jpsi.h"

const double m_Jpsi = 3.0969;
const double m_p = 0.9382720;
const double sthr = pow(m_Jpsi+m_p,2.);
const double Ethr = m_Jpsi + pow(m_Jpsi,2.)/2/m_p;



bool CheckEPCut(double BCAL_EP, double FCAL_EP, double BCAL_Epre)
{
	BCAL_EP = 1./BCAL_EP;
	FCAL_EP = 1./FCAL_EP;
	//cout << "  FCAL: " << FCAL_EP << endl;
	//cout << "  BCAL: " << BCAL_EP << " " << BCAL_Epre << endl;
	if(BCAL_EP > 0.8413 && BCAL_EP < 1.1827) {
		//if( BCAL_Epre < 0.03 && BCAL_Epre > 0.001 )
		if( BCAL_Epre < 0.03 )
			return false;
		else
			return true;
	} else if(FCAL_EP > 0.9562 && FCAL_EP < 1.1521) {
		return true;
	}
	
	return false;

		/*
        if( (FCAL_EP>0.8) || (BCAL_EP>0.8) ) {
                return true;   // NOTE!

                if(BCAL_Epre>0.05)
                        return true;
                else
                        return false;
        } else {
                return false;
        }
		*/
		
		/*
		if( BCAL_EP < 0.001 && FCAL_EP < 0.001)
			return false;
        
        if(BCAL_EP > 0.) {
            if( BCAL_Epre < 0.03 )
                   return false;
           if( BCAL_EP < 0.8413 )
                  return false;
           if( BCAL_EP > 1.1827 )
                  return false;
        }
        if(FCAL_EP > 0.) {
            if( FCAL_EP < 0.9562 )
                 return false;
            if( FCAL_EP > 1.1521 )
                 return false;
                     
        }
        return true;

        */
        
        
}

static bool PRINT_CUTS = false;

void DSelector_jpsi::Init(TTree *locTree)
{
	// USERS: IN THIS FUNCTION, ONLY MODIFY SECTIONS WITH A "USER" OR "EXAMPLE" LABEL. LEAVE THE REST ALONE.

	rand = new TRandom3(0);
	
	// The Init() function is called when the selector needs to initialize a new tree or chain.
	// Typically here the branch addresses and branch pointers of the tree will be set.
	// Init() will be called many times when running on PROOF (once per file to be processed).

	//USERS: SET OUTPUT FILE NAME //can be overriden by user in PROOF
	dOutputFileName = "jpsi.root"; //"" for none
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
	//std::deque<Particle_t> MyPhi;
	//MyPhi.push_back(KPlus); MyPhi.push_back(KMinus);

	//ANALYSIS ACTIONS: //Executed in order if added to dAnalysisActions
	//false/true below: use measured/kinfit data
	/*
	//PID
	dAnalysisActions.push_back(new DHistogramAction_ParticleID(dComboWrapper, false));
	//below: value: +/- N ns, Unknown: All PIDs, SYS_NULL: all timing systems
	//dAnalysisActions.push_back(new DCutAction_PIDDeltaT(dComboWrapper, false, 0.5, KPlus, SYS_BCAL));

	//MASSES
	//dAnalysisActions.push_back(new DHistogramAction_InvariantMass(dComboWrapper, false, Lambda, 1000, 1.0, 1.2, "Lambda"));
	//dAnalysisActions.push_back(new DHistogramAction_MissingMassSquared(dComboWrapper, false, 1000, -0.1, 0.1));

	//KINFIT RESULTS
	dAnalysisActions.push_back(new DHistogramAction_KinFitResults(dComboWrapper));

	//CUT MISSING MASS
	//dAnalysisActions.push_back(new DCutAction_MissingMassSquared(dComboWrapper, false, -0.03, 0.02));

	//BEAM ENERGY
	dAnalysisActions.push_back(new DHistogramAction_BeamEnergy(dComboWrapper, false));
	//dAnalysisActions.push_back(new DCutAction_BeamEnergy(dComboWrapper, false, 8.4, 9.05));

	//KINEMATICS
	dAnalysisActions.push_back(new DHistogramAction_ParticleComboKinematics(dComboWrapper, false));

	// ANALYZE CUT ACTIONS
	// // Change MyPhi to match reaction
	dAnalyzeCutActions = new DHistogramAction_AnalyzeCutActions( dAnalysisActions, dComboWrapper, false, 0, MyPhi, 1000, 0.9, 2.4, "CutActionEffect" );

	//INITIALIZE ACTIONS
	//If you create any actions that you want to run manually (i.e. don't add to dAnalysisActions), be sure to initialize them here as well
	Initialize_Actions();
	dAnalyzeCutActions->Initialize(); // manual action, must call Initialize()
	*/
	
	dBinRange.push_back( pair<double,double>(8.2,8.38) );
	dBinRange.push_back( pair<double,double>(8.38,8.56) );
	dBinRange.push_back( pair<double,double>(8.56,8.74) );
	dBinRange.push_back( pair<double,double>(8.74,8.92) );
	dBinRange.push_back( pair<double,double>(8.92,9.10) );
	dBinRange.push_back( pair<double,double>(9.10,9.28) );
	dBinRange.push_back( pair<double,double>(9.28,9.46) );
	dBinRange.push_back( pair<double,double>(9.46,9.64) );
	dBinRange.push_back( pair<double,double>(9.64,9.82) );
	dBinRange.push_back( pair<double,double>(9.82,10.) );
	dBinRange.push_back( pair<double,double>(10.00,10.18) );
	dBinRange.push_back( pair<double,double>(10.18,10.36) );
	dBinRange.push_back( pair<double,double>(10.36,10.54) );
	dBinRange.push_back( pair<double,double>(10.54,10.72) );
	dBinRange.push_back( pair<double,double>(10.72,10.90) );
	dBinRange.push_back( pair<double,double>(10.90,11.08) );
	dBinRange.push_back( pair<double,double>(11.08,11.26) );
	dBinRange.push_back( pair<double,double>(11.26,11.44) );
	dBinRange.push_back( pair<double,double>(11.44,11.62) );
	dBinRange.push_back( pair<double,double>(11.62,11.8) );

	/******************************** EXAMPLE USER INITIALIZATION: STAND-ALONE HISTOGRAMS *******************************/

	//EXAMPLE MANUAL HISTOGRAMS:
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", ";Missing Mass Squared (GeV/c^{2})^{2}", 600, -0.06, 0.06);
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 12.0);
	dHist_jpsi_BeamEnergy = new TH1I("BeamEnergy_jpsi", ";Beam Energy (GeV)", 600, 0.0, 12.0);

	dHist_RFTime = new TH1F("RFTime", "; RF Time (ns)", 400, -20.0, 20.0);

	dHist_MJpsi = new TH1F("MJpsi", ";M(e^{+}e^{-}) (GeV)", 500, 0.9, 3.4);
	dHist_MJpsi_withcombos = new TH1F("MJpsi_withcombos", ";M(e^{+}e^{-}) (GeV)", 500, 0.9, 3.4);
	dHist_MJpsi_acc = new TH1F("MJpsi_acc", ";M(e^{+}e^{-}) (GeV)", 500, 0.9, 3.4);

	dHist_Ebeam_t = new TH2F("Ebeam_t", ";E(#gamma) (GeV);-t (GeV^{2})", 50, 8.2, 12, 50, 0, 8);


	dHist_MJpsi_final = new TH1F("MJpsi_final", ";M(e^{+}e^{-}) (GeV)", 500, 0.9, 3.4);
	dHist_MJpsi_acc_final = new TH1F("MJpsi_acc_final", ";M(e^{+}e^{-}) (GeV)", 500, 0.9, 3.4);

	dHist_NumGoodEvents = new TH1F("num_good_evts", "", 5, 0, 5);

    hthrown_egang = new TH1F("egang_thrown", ";#theta(e,#gamma)", 2000, 0, 10.);
    hthrown_egang_theta = new TH2F("egang_theta", ";#theta(e,#gamma);#theta(e)", 300, 0, 6., 300, 0, 6.);
    hrecon_egang = new TH1F("egang_recon", ";#theta(e,#gamma)", 2000, 0, 10.);
    hrecon_egang_theta = new TH2F("egang_theta_recon", ";#theta(e,#gamma);#theta(e)", 300, 0, 6., 300, 0, 6.);
    hrecon_egang_all = new TH1F("egang_recon_all", ";#theta(e,#gamma)", 2000, 0, 10.);
    hrecon_egang_theta_all = new TH2F("egang_theta_recon_all", ";#theta(e,#gamma);#theta(e)", 300, 0, 6., 300, 0, 6.);

	int bin=0;
	for(auto range : dBinRange) {
		bin++;
		char name[200], title[100];

		sprintf(name, "MJpsi_bin%d",bin);
		sprintf(title, "e^{+}e^{-} Mass E(#gamma)=%3.1f-%3.1f;M(e^{+}e^{-}) (GeV)",range.first,range.second);
		dHist_MJpsi_bins.push_back( new TH1F(name, title, 70, 2.9, 3.25) );

		sprintf(name, "MJpsi_withcombos_bin%d",bin);
		sprintf(title, "e^{+}e^{-} Mass E(#gamma)=%3.1f-%3.1f;M(e^{+}e^{-}) (GeV)",range.first,range.second);
		dHist_MJpsi_bins_withcombos.push_back( new TH1F(name, title, 70, 2.9, 3.25) );

		sprintf(name, "MJpsi_final_bin%d",bin);
		sprintf(title, "e^{+}e^{-} Mass E(#gamma)=%3.1f-%3.1f;M(e^{+}e^{-}) (GeV)",range.first,range.second);
		dHist_MJpsi_bins_final.push_back( new TH1F(name, title, 70, 2.9, 3.25) );

	}


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

	dIsMC = (dTreeInterface->Get_Branch("MCWeight") != NULL);

}

Bool_t DSelector_jpsi::Process(Long64_t locEntry)
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

	int num_good_evts = 0;
	
	double best_chisq = 1e100;
	double best_mass = -10.;
	double best_t = -10.;
	double best_beam_en = -10.;

	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry

	//if( Get_EventNumber() != 6635852 &&  Get_EventNumber() != 8407933 )  // RUN 11145
	//	return kFALSE;;
	//if( Get_EventNumber() != 16626266 && Get_EventNumber() != 54319305 && Get_EventNumber() != 82009578 && Get_EventNumber() != 82135454 
	//  && Get_EventNumber() != 95854679 && Get_EventNumber() != 122967382 ) // RUN 11448
	//	return kFALSE;;

	//cout << "RUN " << Get_RunNumber() << ", EVENT " << Get_EventNumber() << endl;
	//TLorentzVector locProductionX4 = Get_X4_Production();

	/******************************************** GET POLARIZATION ORIENTATION ******************************************/

	//Only if the run number changes
	//RCDB environment must be setup in order for this to work! (Will return false otherwise)
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//ANALYSIS ACTIONS: Reset uniqueness tracking for each action
	//For any actions that you are executing manually, be sure to call Reset_NewEvent() on them here
	//Reset_Actions_NewEvent();
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
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Mass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Mass_t;

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

		if(PRINT_CUTS) cout << "COMBO #" << loc_i << endl;

		//Set branch array indices for combo and all combo particles
		dComboWrapper->Set_ComboIndex(loc_i);

		// Is used to indicate when combos have been cut
		if(dComboWrapper->Get_IsComboCut()) // Is false when tree originally created
			continue; // Combo has been cut previously

		/********************************************** GET PARTICLE INDICES *********************************************/

		//Used for tracking uniqueness when filling histograms, and for determining unused particles

		//Step 0
		Int_t locBeamID = dComboBeamWrapper->Get_BeamID();
		Int_t locProtonTrackID = dProtonWrapper->Get_TrackID();

		//Step 1
		Int_t locElectronTrackID = dElectronWrapper->Get_TrackID();
		Int_t locPositronTrackID = dPositronWrapper->Get_TrackID();

		/*********************************************** GET FOUR-MOMENTUM **********************************************/

		// Get P4's: //is kinfit if kinfit performed, else is measured
		//dTargetP4 is target p4
		//Step 0
		TLorentzVector locBeamX4 = dComboBeamWrapper->Get_X4();
		TLorentzVector locBeamP4 = dComboBeamWrapper->Get_P4();
		TLorentzVector locProtonP4 = dProtonWrapper->Get_P4();
		//Step 1
		TLorentzVector locElectronP4 = dElectronWrapper->Get_P4();
		TLorentzVector locPositronP4 = dPositronWrapper->Get_P4();

		// Get Measured P4's:
		//Step 0
		TLorentzVector locBeamP4_Measured = dComboBeamWrapper->Get_P4_Measured();
		TLorentzVector locProtonP4_Measured = dProtonWrapper->Get_P4_Measured();
		//Step 1
		TLorentzVector locElectronP4_Measured = dElectronWrapper->Get_P4_Measured();
		TLorentzVector locPositronP4_Measured = dPositronWrapper->Get_P4_Measured();

		/********************************************* COMBINE FOUR-MOMENTUM ********************************************/

		// DO YOUR STUFF HERE

		// Combine 4-vectors
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4;
		locMissingP4_Measured -= locProtonP4_Measured + locElectronP4_Measured + locPositronP4_Measured;

        TLorentzVector locJpsiP4_Measured = locPositronP4_Measured + locElectronP4_Measured;
        TLorentzVector locJpsiP4 = locPositronP4 + locElectronP4;

		/******************************************** EXECUTE ANALYSIS ACTIONS *******************************************/
		/*
		// Loop through the analysis actions, executing them in order for the active particle combo
		dAnalyzeCutActions->Perform_Action(); // Must be executed before Execute_Actions()
		if(!Execute_Actions()) //if the active combo fails a cut, IsComboCutFlag automatically set
			continue;
		*/
	
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

		/**************************************** EVENT SELECTIONS *****************************************/

		double locRFCut = 2.;

		// momentum cuts
		if( locElectronP4.P() < 0.4 ) {
			if(PRINT_CUTS) cout << "ELECTRON P4 CUT: " << locElectronP4.P() << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		if( locPositronP4.P() < 0.4 ) {
			if(PRINT_CUTS) cout << "POSITRON P4 CUT: " << locPositronP4.P() << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		if( locProtonP4.P() < 0.4 ) {
			if(PRINT_CUTS) cout << "PROTON P4 CUT: " << locProtonP4.P() << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		
		// angle cuts
		if( locElectronP4.Theta()*180./TMath::Pi() < 2. ) {
			if(PRINT_CUTS) cout << "ELECTRON THETA CUT: " << locElectronP4.Theta()*180./TMath::Pi() << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		if( locPositronP4.Theta()*180./TMath::Pi() < 2. ) {
			if(PRINT_CUTS) cout << "ELECTRON THETA CUT: " << locPositronP4.Theta()*180./TMath::Pi()  << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}

		// RF time selections
		double locDeltaTRF = locBeamX4.T() - (dComboWrapper->Get_RFTime() 
						  + (locBeamX4.Z() - dComboWrapper->Get_TargetCenter().Z() )/29.9792458);
		dHist_RFTime->Fill(locDeltaTRF);

		double weight_acc = 0.;
		if(fabs(locDeltaTRF) > locRFCut) {
			if(PRINT_CUTS) cout << "RF cut: " << fabs(locDeltaTRF)  << endl;
			weight_acc = -1./8.;   // depends on the number of beam bunches used 
			//dComboWrapper->Set_IsComboCut(true);
			//continue;
		} else {    // in-time
			weight_acc = 1.;
		}
		
		
		// Kin fit chi^2
		if(dComboWrapper->Get_ChiSq_KinFit("") < 0. || dComboWrapper->Get_ChiSq_KinFit("") > 5000.) { 
			if(PRINT_CUTS) cout << "KIN FIT CHISQ CUT: " << dComboWrapper->Get_ChiSq_KinFit("")  << endl;
			dComboWrapper->Set_IsComboCut(true);
			continue;
		}
		

		// ELECTRON PID 
		double Positron_BCAL_EP = (dPositronWrapper->Get_Energy_BCAL()/locPositronP4.P());
		double Positron_FCAL_EP = (dPositronWrapper->Get_Energy_FCAL()/locPositronP4.P());
		double Positron_BCAL_EpreP =(dPositronWrapper->Get_Energy_BCALPreshower()) * sin(locPositronP4.Theta());

		double Electron_BCAL_EP = (dElectronWrapper->Get_Energy_BCAL()/locElectronP4.P());
		double Electron_FCAL_EP = (dElectronWrapper->Get_Energy_FCAL()/locElectronP4.P());
		double Electron_BCAL_EpreP =(dElectronWrapper->Get_Energy_BCALPreshower()) * sin(locElectronP4.Theta());

		//cout << "ELECTRON:  "; CheckEPCut(Electron_BCAL_EP, Electron_FCAL_EP, Electron_BCAL_EpreP);  cout << endl;
		//cout << "POSITRON:  "; CheckEPCut(Positron_BCAL_EP, Positron_FCAL_EP, Positron_BCAL_EpreP);  cout << endl;E
		if(!CheckEPCut(Positron_BCAL_EP, Positron_FCAL_EP, Positron_BCAL_EpreP)
			|| !CheckEPCut(Electron_BCAL_EP, Electron_FCAL_EP, Electron_BCAL_EpreP)) {
			dComboWrapper->Set_IsComboCut(true);
			
			if(PRINT_CUTS) {
				if(!CheckEPCut(Positron_BCAL_EP, Positron_FCAL_EP, Positron_BCAL_EpreP))
					cout << "POSITRON EP CUT (EP FCAL/EP BCAL/EpreP BCAL): " << Positron_FCAL_EP
					     << " / " << Positron_BCAL_EP << " / " << Positron_BCAL_EpreP << endl;
				if(!CheckEPCut(Electron_BCAL_EP, Electron_FCAL_EP, Electron_BCAL_EpreP))
					cout << "ELECTRON EP CUT (EP FCAL/EP BCAL/EpreP BCAL): " << Electron_FCAL_EP
					     << " / " << Electron_BCAL_EP << " / " << Electron_BCAL_EpreP << endl;
			}

			continue;
		}

		/**************************************** HISTOGRAM BEAM ENERGY *****************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamID) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4.E());
			locUsedSoFar_BeamEnergy.insert(locBeamID);
		}

		/************************************ HISTOGRAM MISSING MASS SQUARED ************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Uniqueness tracking: Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy).
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamID); //beam
		locUsedThisCombo_MissingMass[Proton].insert(locProtonTrackID);
		locUsedThisCombo_MissingMass[Electron].insert(locElectronTrackID);
		locUsedThisCombo_MissingMass[Positron].insert(locPositronTrackID);

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}
		
		// calculate kinematic variables
		double t = (locProtonP4 - dTargetP4).M2();
		
		double s = m_p*(m_p + 2.*locBeamP4.E());
		double q_g = (s - m_p*m_p)/2./sqrt(s);
		double q_Jpsi = 1./2./sqrt(s) * sqrt((s - pow(m_Jpsi + m_p,2.))*(s - pow(m_Jpsi - m_p,2.)));
		double t0 = pow(m_Jpsi*m_Jpsi,2.)/4./s - pow(q_g - q_Jpsi,2.);
		double t1 = t0 - 4*q_g*q_Jpsi;
		double tp = t - t0;

		// TEST
		if(locBeamP4.E() > 9. && locBeamP4.E() < 9.5) {
			dHist_ptheta_t->Fill(locProtonP4.Theta()*180./TMath::Pi(), fabs(t));
		}
		
		//locUsedSoFar_Mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_Mass;
		locUsedThisCombo_Mass[Electron].insert(locElectronTrackID);
		locUsedThisCombo_Mass[Positron].insert(locPositronTrackID);

		map<Particle_t, set<Int_t> > locUsedThisCombo_Mass_t;
		locUsedThisCombo_Mass_t[Proton].insert(locProtonTrackID);
		locUsedThisCombo_Mass_t[Electron].insert(locElectronTrackID);
		locUsedThisCombo_Mass_t[Positron].insert(locPositronTrackID);
			
		// weighted histograms
		if(locUsedSoFar_Mass.find(locUsedThisCombo_Mass) == locUsedSoFar_Mass.end()) {
			dHist_MJpsi_acc->Fill(locJpsiP4.M(), weight_acc);
		}
		
		// all in time events
		if(weight_acc > 0.) {
			dHist_jpsi_BeamEnergy->Fill(locBeamP4.E());

			// include all combos
			dHist_MJpsi_withcombos->Fill(locJpsiP4.M());
			int bin=0;
			for(auto range : dBinRange) {
				if( locBeamP4.E()>range.first && locBeamP4.E()<range.second ) {
					dHist_MJpsi_bins_withcombos[bin]->Fill(locJpsiP4.M());
					break;
				}
				bin++;
			}
		
			// mass
			if(locUsedSoFar_Mass.find(locUsedThisCombo_Mass) == locUsedSoFar_Mass.end()) {
				locUsedSoFar_Mass.insert(locUsedThisCombo_Mass);
				
				if(best_chisq > dComboWrapper->Get_ChiSq_KinFit("")) {
					best_chisq = dComboWrapper->Get_ChiSq_KinFit("");
					best_mass = locJpsiP4.M();
					best_t = t;
					best_beam_en = locBeamP4.E();
				}
				
				// Fill masses
				dHist_MJpsi->Fill(locJpsiP4.M());

				int bin=0;
				for(auto range : dBinRange) {
					if( locBeamP4.E()>range.first && locBeamP4.E()<range.second ) {
						dHist_MJpsi_bins[bin]->Fill(locJpsiP4.M());
					}
					bin++;
				}
			}
			
			// mass vs. t
			if(locUsedSoFar_Mass_t.find(locUsedThisCombo_Mass_t) == locUsedSoFar_Mass_t.end()) {
				locUsedSoFar_Mass_t.insert(locUsedThisCombo_Mass_t);
		
				if(locJpsiP4.M() > 3.04 && locJpsiP4.M() < 3.13) {
					dHist_Ebeam_t->Fill(locBeamP4.E(), fabs(t));
				}
			}

		}

		/****************************************** FILL FLAT TREE (IF DESIRED) ******************************************/

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
	//Fill_NumCombosSurvivedHists();

	/******************************************* LOOP OVER THROWN DATA (OPTIONAL) ***************************************/

	//Thrown beam: just use directly
	if(dThrownBeam != NULL)
		double locEnergy = dThrownBeam->Get_P4().E();

    TLorentzVector ep4v(0.,0.,0.,0.);
    TLorentzVector em4v(0.,0.,0.,0.);
    vector<TLorentzVector> gam4vs;

	//Loop over throwns
	for(UInt_t loc_i = 0; loc_i < Get_NumThrown(); ++loc_i)
	{
		//Set branch array indices corresponding to this particle
		dThrownWrapper->Set_ArrayIndex(loc_i);

        //cout << dThrownWrapper->Get_PID() << endl;

		//Do stuff with the wrapper here ...
        if((int)dThrownWrapper->Get_PID() == 2)
            em4v = dThrownWrapper->Get_P4();
        if((int)dThrownWrapper->Get_PID() == 3)
            ep4v = dThrownWrapper->Get_P4();

        if((int)dThrownWrapper->Get_PID() == 1)
            gam4vs.push_back(dThrownWrapper->Get_P4());
	}

    /*
    cout << "electron:  "; em4v.Print(); cout << endl;
    cout << "positron:  "; ep4v.Print(); cout << endl;
    cout << "photons = " << gam4vs.size() << endl;
    */

    for(int i=0; i<gam4vs.size(); i++) {
        double angle = (180./TMath::Pi())*em4v.Vect().Angle(gam4vs[i].Vect());
        double angle2 = (180./TMath::Pi())*ep4v.Vect().Angle(gam4vs[i].Vect());
        //cout << "angle = " << angle << endl;
        double theta = (180./TMath::Pi())*em4v.Vect().Theta();
        double theta2 = (180./TMath::Pi())*ep4v.Vect().Theta();


        if(angle < angle2) {
            hthrown_egang->Fill(angle);
            hthrown_egang_theta->Fill(theta, angle);
        } else {
            hthrown_egang->Fill(angle2);
            hthrown_egang_theta->Fill(theta2, angle2);

        }
    }

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

	// fill overall histograms
	dHist_NumGoodEvents->Fill(num_good_evts);
	
	if(best_chisq < 5000.) {

		dHist_MJpsi_final->Fill(best_mass);

		int bin=0;
		for(auto range : dBinRange) {
			if( best_mass>range.first && best_mass<range.second ) {
				dHist_MJpsi_bins_final[bin]->Fill(best_mass);
			}
			bin++;
		}

	}

	return kTRUE;
}

void DSelector_jpsi::Finalize(void)
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
