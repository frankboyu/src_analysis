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

void test_selector()
{
    gROOT->ProcessLine(".x $(ROOT_ANALYSIS_HOME)/scripts/Load_DSelector.C");

    TChain *chain = new TChain("gd_pimprotinc__B4_F4_T1_S4_Tree");

	chain->Add("/cache/halld/home/boyu/src_analysis/sim/piminus_p_12C/ver01/tree_thrown_gen_MF/merged/*.root");

    chain->Process("selectors/DSelector_piminus_p_12C_thrown.C+","");
}
