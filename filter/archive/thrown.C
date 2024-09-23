#include <TString.h>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include <iostream>

#include "TLorentzVector.h"

void thrown(){

    TChain *chain = new TChain("genT");
    chain->Add("/work/halld2/home/boyu/src_sim/SRC-CT_Simulation/output/pim_p_2H_MF/model_hist_0.1cut_nobkg_2M/root/generator/*.root");
    
    TLorentzVector *pBeam = 0;
    TLorentzVector *pMeson = 0;
    TLorentzVector *pBaryon = 0;
    
    double mandelstam_s = 0;
    double mandelstam_t = 0;
    
    int thrown_events[5][10];
    
    for (int i=0; i<5; i++)
    {
        for (int j=0; j<10; j++)
        {
            thrown_events[i][j] = 0;
        }
    }

    chain->SetBranchAddress("pBeam", &pBeam);
    chain->SetBranchAddress("pMeson", &pMeson);
    chain->SetBranchAddress("pBaryon", &pBaryon);
  
    for (int index=0; index<chain->GetEntries(); index++)
    {
        chain->GetEntry(index);
        mandelstam_s = (*pMeson + *pBaryon).Mag2();
        mandelstam_t = -(*pBeam - *pMeson).Mag2();
        
        double bin_S[6] = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
        double bin_T[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0};
        for(int i = 0; i < 5; i++)
        {
            for(int j = 0; j < 10; j++)
            {
                if(mandelstam_s >=  bin_S[i] && mandelstam_s < bin_S[i+1] && mandelstam_t >= bin_T[j] && mandelstam_t < bin_T[j+1])
                {
                    thrown_events[i][j]++;
                    continue;
                }
            }
        }

    }
    
    for (int i=0; i<5; i++)
    {
        for (int j=0; j<10; j++)
        {
            cout << thrown_events[i][j] << endl;
        }
    }

  }


 
  
  
