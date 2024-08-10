#include <iostream>
#include <fstream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TProof.h"
#include "TProofDebug.h"

using namespace std;
R__LOAD_LIBRARY(libDSelector)

int process_chain(string locInputFileName, string locTreeName, string locSelectorName, unsigned int locNThreads)
{
	//tell it to compile selector (if user did not)
	if(locSelectorName[locSelectorName.size() - 1] != '+')
		locSelectorName += '+';

	cout << "file name, tree name, selector name, #threads = " << locInputFileName << ", " << locTreeName << ", " << locSelectorName << ", " << locNThreads << endl;

	//process chain directly
	gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");
	TChain *locChain = new TChain(locTreeName.c_str());
	locChain->Add(locInputFileName.c_str());

	Long64_t locStatus = locChain->Process(locSelectorName.c_str());
	cout << "chain status = " << locStatus << endl;
	return ((locStatus >= Long64_t(0)) ? 0 : 999); //0 = success
}


