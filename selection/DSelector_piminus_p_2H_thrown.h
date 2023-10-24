#ifndef DSelector_piminus_p_2H_h
#define DSelector_piminus_p_2H_h

#include <iostream>
#include <string>

#include "DSelector/DSelector.h"

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

#endif
