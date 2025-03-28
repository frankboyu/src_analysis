#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

class DSelector_rho_d_thrown : public DSelector
{
public:

    DSelector_rho_d_thrown(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_rho_d_thrown(){}

    void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);

private:

    void Finalize(void);

    // BEAM POLARIZATION INFORMATION
    UInt_t  dPreviousRunNumber;
    bool    dIsPolarizedFlag;
    bool    dIsPARAFlag;

    ClassDef(DSelector_rho_d_thrown, 0);
};

void DSelector_rho_d_thrown::Init(TTree *locTree)
{
    // SET OUTPUT FILE NAME
    dOutputFileName          = "";
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = "selectedtree_rho_d_thrown.root";
    dFlatTreeName            = "selectedtree_rho_d_thrown";
    dSaveDefaultFlatBranches = true;
    dSkipNoTriggerEvents     = false;

	// INITIALIZE THE TREE INTERFACE
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pip_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_x4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("beam_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pip_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("pim_p4_truth");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("d_p4_truth");
}
// END OF INITIALIZATION

Bool_t DSelector_rho_d_thrown::Process(Long64_t locEntry)
{
	// CALL THIS FIRST
	DSelector::Process(locEntry); // gets the data from the tree for the entry

    // GET BEAM POLARIZATION INFO. RCDB ENVIRONMENT REQUIRED
	UInt_t locRunNumber = Get_RunNumber();
	if(locRunNumber != dPreviousRunNumber)
	{
		dIsPolarizedFlag   = dAnalysisUtilities.Get_IsPolarizedBeam(locRunNumber, dIsPARAFlag);
		dPreviousRunNumber = locRunNumber;
	}

        //GET THROWN P4
        TLorentzVector locBeamX4_Thrown, locPiPlusX4_Thrown, locPiMinusX4_Thrown, locDeuteronX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locPiPlusP4_Thrown, locPiMinusP4_Thrown, locDeuteronP4_Thrown;
        locBeamX4_Thrown = dThrownBeam->Get_X4();
        locBeamP4_Thrown = dThrownBeam->Get_P4();
        dThrownWrapper->Set_ArrayIndex(0);
        locPiPlusX4_Thrown = dThrownWrapper->Get_X4();
        locPiPlusP4_Thrown = dThrownWrapper->Get_P4();
        dThrownWrapper->Set_ArrayIndex(1);
        locPiMinusX4_Thrown = dThrownWrapper->Get_X4();
        locPiMinusP4_Thrown = dThrownWrapper->Get_P4();
        dThrownWrapper->Set_ArrayIndex(2);
        locDeuteronX4_Thrown = dThrownWrapper->Get_X4();
        locDeuteronP4_Thrown = dThrownWrapper->Get_P4();

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_x4_truth", locBeamX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pip_x4_truth", locPiPlusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_x4_truth", locPiMinusX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_x4_truth", locDeuteronX4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("beam_p4_truth", locBeamP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pip_p4_truth", locPiPlusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("pim_p4_truth", locPiMinusP4_Thrown);
        dFlatTreeInterface->Fill_TObject<TLorentzVector>("d_p4_truth", locDeuteronP4_Thrown);
        Fill_FlatTree();

	return kTRUE;
}
// END OF PROCESSING

void DSelector_rho_d_thrown::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION