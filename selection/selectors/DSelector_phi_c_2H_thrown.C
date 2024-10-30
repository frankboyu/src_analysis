#include <iostream>
#include <string>

#include "DSelector/DSelector.h"

double RadToDeg = 180.0/3.1415926;

class DSelector_phi_c_2H_thrown : public DSelector
{
	public:

		DSelector_phi_c_2H_thrown(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_phi_c_2H_thrown(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Finalize(void);

		//BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag;      //Else is AMO
		bool dIsPARAFlag;           //Else is PERP or AMO

	ClassDef(DSelector_phi_c_2H_thrown, 0);
};

void DSelector_phi_c_2H_thrown::Init(TTree *locTree)
{
    // SET OUTPUT FILE NAME
    dOutputFileName          = "";
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = "flattree_phi_c_2H_thrown.root";
    dFlatTreeName            = "flattree_phi_c_2H_thrown";
    dSaveDefaultFlatBranches = false;
    dSkipNoTriggerEvents     = false;

	// INITIALIZE THE TREE INTERFACE
	bool locInitializedPriorFlag = dInitializedFlag; //Save whether have been initialized previously
	DSelector::Init(locTree);                        //This must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                      //Have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;

}
// END OF INITIALIZATION

Bool_t DSelector_phi_c_2H_thrown::Process(Long64_t locEntry)
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

    // FILL FLAT TREE
    Fill_FlatTree();

	return kTRUE;
}
//END OF PROCESSING

void DSelector_phi_c_2H_thrown::Finalize(void)
{
	//CALL THIS LAST
	DSelector::Finalize(); //Saves results to the output file
}
//END OF FINALIZATION
