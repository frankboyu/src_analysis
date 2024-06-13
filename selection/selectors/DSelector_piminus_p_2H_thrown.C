#include <iostream>
#include <string>

#include "DSelector/DSelector.h"

double RadToDeg = 180.0/3.1415926;

class DSelector_piminus_p_2H_thrown : public DSelector
{
	public:

		DSelector_piminus_p_2H_thrown(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_piminus_p_2H_thrown(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Finalize(void);

		//BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag;      //Else is AMO
		bool dIsPARAFlag;           //Else is PERP or AMO

	ClassDef(DSelector_piminus_p_2H_thrown, 0);
};

void DSelector_piminus_p_2H_thrown::Init(TTree *locTree)
{
    // SET OUTPUT FILE NAME 
    dOutputFileName          = "";
    dOutputTreeFileName      = "";
    dFlatTreeFileName        = "flattree_piminus_p_2H_thrown.root";
    dFlatTreeName            = "flattree_piminus_p_2H_thrown";
    dSaveDefaultFlatBranches = false;
    dSkipNoTriggerEvents     = false;

	// INITIALIZE THE TREE INTERFACE
	bool locInitializedPriorFlag = dInitializedFlag; //Save whether have been initialized previously
	DSelector::Init(locTree);                        //This must be called to initialize wrappers for each new TTree
	if(locInitializedPriorFlag)
		return;                                      //Have already created histograms, etc. below: exit
	dPreviousRunNumber = 0;

    // CUSTOM OUTPUT BRACHES: FLAT TREE
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("WeightFactor");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("MissingPMinus");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("BeamEnergy");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("MinusT");
    dFlatTreeInterface->Create_Branch_Fundamental<Double_t>("thetaCM");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("BeamP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("PiMinusP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("ProtonP4_Thrown");
    dFlatTreeInterface->Create_Branch_NoSplitTObject<TLorentzVector>("MissingP4_Thrown");
}
// END OF INITIALIZATION

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

    TVector3       boostCM          = (locPiMinusP4_Thrown   + locProtonP4_Thrown).BoostVector();
    TLorentzVector locBeamP4CM_Thrown      = locBeamP4_Thrown;
    TLorentzVector locPiMinusP4CM_Thrown   = locPiMinusP4_Thrown;
    locBeamP4CM_Thrown.Boost(-boostCM);
    locPiMinusP4CM_Thrown.Boost(-boostCM);

    double locMinusT_Thrown          = -(locBeamP4_Thrown   - locPiMinusP4_Thrown).Mag2();
    double locThetaCM_Thrown         = locBeamP4CM_Thrown.Vect().Angle(locPiMinusP4CM_Thrown.Vect())*RadToDeg;

    //FILL CUSTOM BRANCHES: FLAT TREE
    dFlatTreeInterface->Fill_Fundamental<Double_t>("WeightFactor",          1.0);
    dFlatTreeInterface->Fill_Fundamental<Double_t>("MissingPMinus",         locMissingP4_Thrown.Minus());
    dFlatTreeInterface->Fill_Fundamental<Double_t>("BeamEnergy",            locBeamP4_Thrown.E());
    dFlatTreeInterface->Fill_Fundamental<Double_t>("MinusT",                locMinusT_Thrown);
    dFlatTreeInterface->Fill_Fundamental<Double_t>("thetaCM",               locThetaCM_Thrown);
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
