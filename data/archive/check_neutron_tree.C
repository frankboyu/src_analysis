

{

    TFile f("tree_kpkm__F2_B4_T0.root");
    //TFile f("/work/halld2/home/boyu/src_phi_neutron/tree_kpkm__F2_B4_T0.root");
    TTree *T = (TTree*)f.Get("kpkm__F2_B4_T0_Tree");
    
    cout << T->GetEntries() << endl;

    ULong64_t numevent;
    UInt_t numNeutrals;
    TClonesArray *NeutralP4s = new TClonesArray();
    
    T->SetBranchAddress("EventNumber", &numevent);
    T->SetBranchAddress("NumNeutralHypos", &numNeutrals);
    T->SetBranchAddress("NeutralHypo__P4_Measured", &NeutralP4s);

    for (Long64_t i=0;i<T->GetEntries();i++) {
        T->GetEntry(i);

        if(numNeutrals != NeutralP4s->GetEntries())
            cout << numevent << " " << numNeutrals << " " << NeutralP4s->GetEntries() << endl;
    }

}
