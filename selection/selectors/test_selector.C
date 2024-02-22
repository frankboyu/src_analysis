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
    
    //TChain *chain = new TChain("gd_pimprotmissprot__B4_F4_T1_S4_Tree");
    TChain *chain = new TChain("ghe_pimprotinc__B4_F4_T1_S4_Tree");
	// TChain *chain = new TChain("Thrown_Tree");
    //TChain *chain = new TChain("gd_kpkmprotmissn__B4_F4_T2_S5_Tree");

    chain->Add("/cache/halld/RunPeriod-2021-11/analysis/ver06/tree_ghe_pimprotinc__B4_F4_T1_S4/merged/tree_ghe_pimprotinc__B4_F4_T1_S4_090034.root");
    //chain->Add("/cache/halld/RunPeriod-2021-11/analysis/ver03/tree_gd_pimprotmissprot__B4_F4_T1_S4/merged/*.root");
    // chain->Add("/volatile/halld/home/boyu/src_analysis/sim/piminus_p_4He/ver01/root/trees/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF_090167_000.root");
    //chain->Add("/volatile/halld/home/boyu/src_sim/piminus_p_2H_MF_temp/flat_hist_0.0cut_nobkg_2M/root/trees/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF_*.root");
	//chain->Add("/volatile/halld/home/boyu/src_sim/piminus_p_2H_MF_temp/flat_hist_0.0cut_nobkg_2M/root/thrown/*.root");    
    //chain->Add("/volatile/halld/home/boyu/src_analysis/sim/piminus_p_2H_MF/root/trees/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF_*.root");
    // chain->Add("/volatile/halld/home/boyu/src_analysis/sim/piminus_p_2H_MF/root/thrown/*.root");
    //chain->Add("/volatile/halld/home/boyu/src_analysis/skim/tree_gd_pimprotmissprot__B4_F4_T1_S4/090220/*.root");
    
    chain->Process("DSelector_piminus_p_recon.C+","");
    // chain->Process("selectors/DSelector_piminus_p_2H_thrown.C+","");
    //chain->Process("DSelector_phi_p_2H_recon.C+","");
}
