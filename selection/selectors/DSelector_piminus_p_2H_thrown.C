#include "DSelector_piminus_p_2H_thrown.h"

void DSelector_piminus_p_2H_thrown::Init(TTree *locTree)
{
    //SET OUTPUT FILE NAME 
    dOutputFileName          = "";                                  //"" for none
    dOutputTreeFileName      = "";                                  //"" for none
    dFlatTreeFileName        = "flattree_piminus_p_2H_thrown.root"; //Output flat tree (one combo per tree entry), "" for none
    dFlatTreeName            = "piminus_p_2H_thrown";               //If blank, default name will be chosen
    dSaveDefaultFlatBranches = false;                               //False: don't save default branches, reduce disk footprint.

	//INITIALIZE THE TREE INTERFACE
	bool locInitializedPriorFlag = dInitializedFlag; //Save whether have been initialized previously
	DSelector::Init(locTree);                        //This must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                      //Have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;

    //CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_Thrown");
}
//END OF INITIALIZATION

Bool_t DSelector_piminus_p_2H_thrown::Process(Long64_t locEntry)
{
	//CALL THIS FIRST
	DSelector::Process(locEntry); //Gets the data from the tree for the entry

    //GET PHOTON POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

	//GET THROWN P4
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

    //FILL CUSTOM BRANCHES: FLAT TREE
    dFlatTreeInterface->Fill_TObject<TLorentzVector>("BeamP4_Thrown",    locBeamP4_Thrown);
    dFlatTreeInterface->Fill_TObject<TLorentzVector>("PiMinusP4_Thrown", locPiMinusP4_Thrown);
    dFlatTreeInterface->Fill_TObject<TLorentzVector>("ProtonP4_Thrown",  locProtonP4_Thrown);
    dFlatTreeInterface->Fill_TObject<TLorentzVector>("MissingP4_Thrown", locMissingP4_Thrown);

    // FILL FLAT TREE
    Fill_FlatTree();

	return kTRUE;
}
//END OF PROCESSING

void DSelector_piminus_p_2H_thrown::Finalize(void)
{
	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
//END OF FINALIZATION
