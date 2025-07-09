#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

class DSelector_phi_d_thrown : public DSelector
{
public:

    DSelector_phi_d_thrown(TTree* locTree = NULL) : DSelector(locTree){}
    virtual ~DSelector_phi_d_thrown(){}

    void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);

private:

    void Finalize(void);

    // BEAM POLARIZATION INFORMATION
    UInt_t  dPreviousRunNumber;
    bool    dIsPolarizedFlag;
    bool    dIsPARAFlag;
    int     dPolarizationAngle;

    // FLAGS
    bool dIsMC;

    ClassDef(DSelector_phi_d_thrown, 0);
};

void DSelector_phi_d_thrown::Init(TTree *locTree)
{
    // SET OUTPUT FILE NAME
    dOutputFileName          = "";
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = "selectedtree_phi_d_thrown.root";
    dFlatTreeName            = "selectedtree_phi_d_thrown";
    dSaveDefaultFlatBranches = true;
    dSkipNoTriggerEvents     = false;

	// INITIALIZE THE TREE INTERFACE
    bool locInitializedPriorFlag = dInitializedFlag;               // save whether have been initialized previously
	DSelector::Init(locTree);                                      // this must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                                    // have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;

    // CUSTOM OUTPUT BRACHES: FLAT TREE
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

Bool_t DSelector_phi_d_thrown::Process(Long64_t locEntry)
{
	// CALL THIS FIRST
	DSelector::Process(locEntry); // gets the data from the tree for the entry

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

        //GET THROWN P4 AND TOPOLOGY
        TLorentzVector locBeamX4_Thrown, locKPlusX4_Thrown, locKMinusX4_Thrown, locDeuteronX4_Thrown;
        TLorentzVector locBeamP4_Thrown, locKPlusP4_Thrown, locKMinusP4_Thrown, locDeuteronP4_Thrown;
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
                else
                    cout << "Unexpected PID: " << dThrownWrapper->Get_PID() << endl;
            }
            locDeuteronX4_Thrown = locKPlusX4_Thrown;  // workaround for the missing deuteron info in the tree
            locDeuteronP4_Thrown = locBeamP4_Thrown + TLorentzVector(0, 0, 0, 1.875612859) - locKPlusP4_Thrown - locKMinusP4_Thrown; // workaround for the missing deuteron info in the tree
        }

        // FILL FLAT TREE
        dFlatTreeInterface->Fill_Fundamental<Int_t>("thrown_topology", Get_ThrownTopologyString().Data());
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

	return kTRUE;
}
// END OF PROCESSING

void DSelector_phi_d_thrown::Finalize(void)
{
	// CALL THIS LAST
	DSelector::Finalize(); // saves results to the output file
}
// END OF FINALIZATION