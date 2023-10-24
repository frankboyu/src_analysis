
#include "protonrho_candidates.h"
#include "rootalias.h"
#include <dirent.h>

void protonrho_candidates(
	TString infile 		  	  = "tree_1p2pi",
	//string infile_folder 	  = "/work/halld2/home/nathaly/test/HallD_SRC-CT_Analysis/analysis_scripts/proton_rho0/sim/02/",
	string infile_folder     = "/work/halld2/home/nathaly/simresutls/MFrho",
	//string infile_folder     = "/work/halld2/home/nathaly/simresutls/SRCrho",
	//string infile_folder     = "./",
	//string infile_folder     = "./",
	//string infile_folder     = "/w/halld-scshelf2001/halld2/home/sns/test/MF",
	//string infile_folder     = "/w/halld-scshelf2001/halld2/home/sns/test/SRC",
	
	TString outfile 	  	  = "rhocandidates1",
	TString outfile_folder 	  = "./",
	int maxNEvent        = -0,
    	Int_t verbosity           = 0,
	//double weight	= 6.59575e+00  //MF
	//double weight  = 2.03753e+01  //SRC
	double weight  = 1  //SRC
) 
{  
  
  //===== Infile ==========================================
  const char * dataPath = infile_folder.c_str();
  vector<string> inputFiles;
  
  DIR *dirp;
  dirp = opendir(dataPath);
  regex extension ("(.*)(.root)");
 int numberruns  = 0; 
  if (dirp)
    {
      string dataPathString = dataPath;
      if (dataPathString.back()!='/')
	dataPathString = dataPathString + "/";
      
      struct dirent *directory;
      while ((directory = readdir(dirp)))
	{
	  string fileName = directory->d_name;
	  string fullName = dataPathString + fileName;
	  
          if (!regex_match(fileName,regex("(tree_1p2pi)(.*)")))
            continue;
	  
	  inputFiles.push_back(fullName);
	  
	}
      
      closedir(dirp);
      
    }
  else if (ifstream(dataPath))
    {
      string fileName = dataPath;
      
      if (!regex_match(fileName,extension))
        {
          cout << "Invalid input. Exiting...\n";
          return;
        }
      
      inputFiles.push_back(fileName);
    }
  else
    {
      cout << "Input data not found. Exiting...\n";
      return;
    }
  
  TChain *T = new TChain("tree_1p2pi");
  //T->Add(Form(""))
  
  sort(inputFiles.begin(),inputFiles.end());
  
  int ilt=0;
  for (const auto& fullName: inputFiles)
    {
      
      cout << fullName << "\n";
      
      TFile *f = new TFile(fullName.c_str());      
      
      if (f->GetListOfKeys()->Contains("tree_1p2pi"))
	{
	  //inTree = (TTree *)f->Get("tree_0p1pi1pi0");
	  T->Add(Form("%s",fullName.c_str()));
	  numberruns++;
	}
      else
	{
	  cout << "Data tree not found, abandoning file\n";
	  continue;
	}
      ilt++;
    }
  
  /*
    if(!(infile.EndsWith(".root"))) infile.Append(".root");
    if(infile_folder.EndsWith("/")) infile_folder.Remove(infile_folder.Length()-0);
    TString inputfile = Form("%s/%s",infile_folder.Data(),infile.Data());
    
    if(!gSystem->AccessPathName(inputfile)){
    cout << "Reading file: " << inputfile << endl;
    } else {
    cout << "File: " << inputfile << " does not exists" << endl;
    return 0;
    }*/
  
  //===== Outfile ==========================================
  
  if(!(outfile.EndsWith(".root"))) outfile.Append(".root");
  TString outputfile = Form("%s/%s",outfile_folder.Data(),outfile.Data());
  
  if(!gSystem->AccessPathName(outfile_folder)){
    cout << "Output file: " << outputfile << endl;
  } else {
    cout << "File: " << outfile_folder << " path does not exists" << endl;
    return 0;
  }
  
  
  TFile *MyFile = new TFile(Form("%s",outputfile.Data()),"RECREATE");
  
  TTree *EvTree = new TTree( "EvTree", "EvTree" );
  
  SetBranchAddressesTree(T);	
  SetBranchAddressesOutTree(EvTree);
  
  
  Long64_t tentries = T->GetEntries();
  cout << "Number of events in trees: " << tentries << endl;
  if(maxNEvent>0 && maxNEvent<tentries){
    tentries = maxNEvent;
    cout << "Will only analyze first " << maxNEvent << " events in the trees..." << endl;
  }
  
  
  Double_t CL_cut = 0.00002; 
  Int_t norecon = 0;
  Double_t commontime = 0;
  Double_t t0_protcand = 0;
  
  TLorentzVector target(0.0, 0.0, 0.0, 0.938272088);
  Double_t CL;
  
  Double_t thetamin = 0;
  Double_t massmin1 = 0.5; Double_t massmax1 = 1.4; 
  Double_t phimin   = 170;  Double_t phimax = 190; 
  double preltheta_p, preltheta_pip;
  double prelPmag_p, prelPmag_pip; 

  for(Long64_t i = 0; i < tentries; i++  ){
    
    T->GetEntry(i);

    if (_nShower>10) continue; 

    TLorentzVector hyp1pim(_pX_piminus[0],_pY_piminus[0],_pZ_piminus[0], _E_piminus[0]);
    TLorentzVector hyp1pip(_pX_piplus[0],_pY_piplus[0],_pZ_piplus[0], _E_piplus[0]);
    TLorentzVector hyp1p(_pX_proton[0],_pY_proton[0],_pZ_proton[0], _E_proton[0]);
    TLorentzVector hyp1rho = hyp1pim + hyp1pip;
    
    
    pX_pim = _pX_piminus[0]; pY_pim = _pY_piminus[0]; pZ_pim = _pZ_piminus[0]; E_pim = _E_piminus[0];
    pX_pip = _pX_piplus[0]; pY_pip = _pY_piplus[0]; pZ_pip = _pZ_piplus[0]; E_pip = _E_piplus[0];
    pX_p = _pX_proton[0]; pY_p = _pY_proton[0]; pZ_p = _pZ_proton[0]; E_p = _E_proton[0];
    preltheta_p = hyp1p.Theta(); preltheta_pip = hyp1pip.Theta(); 
    prelPmag_p = hyp1p.Vect().Mag();  prelPmag_pip = hyp1pip.Vect().Mag(); 
    
    X_vtx = _X_vertex[0]; Y_vtx = _Y_vertex[0]; Z_vtx = _Z_vertex[0]; 
    t0_protcand = _t_proton[0]; CL = _CLKinFit[0]; t_vtx = _T_vertex[0];
    commontime = _Common_Time[0]; 
  
   FOM_piminuscand = _FOM_piminuscand[0]; NDF_piminuscand = _NDF_piminuscand[0]; ChiSq_piminuscand = _ChiSq_piminuscand[0];
piminus_dedx_dc_NDF = _piminus_dedx_dc_NDF[0]; piminus_dedx_dc_ChiSq = _piminus_dedx_dc_ChiSq[0];  Beta_Timing_piminuscand = _Beta_Timing_piminuscand[0]; ChiSq_Timing_piminuscand = _ChiSq_Timing_piminuscand[0];
   dEdx_TOF_piminuscand = _dEdx_TOF_piminuscand[0]; Energy_BCAL_piminuscand = _Energy_BCAL_piminuscand[0]; Energy_BCALPreshower_piminuscand=_Energy_BCALPreshower_piminuscand[0]; SigLong_BCAL_piminuscand =_SigLong_BCAL_piminuscand[0]; SigTheta_BCAL_piminuscand = _SigTheta_BCAL_piminuscand[0]; RMSTime_BCAL_piminuscand = _RMSTime_BCAL_piminuscand[0];  E_BCAL_piminuscand = _E_BCAL_piminuscand[0];  Energy_FCAL_piminuscand=_Energy_FCAL_piminuscand[0]; E1E9_FCAL_piminuscand=_E1E9_FCAL_piminuscand[0]; E9E25_FCAL_piminuscand=_E9E25_FCAL_piminuscand[0]; SumU_FCAL_piminuscand=_SumU_FCAL_piminuscand[0]; SumV_FCAL_piminuscand=_SumV_FCAL_piminuscand[0]; TrackBCAL_DeltaPhi_piminuscand = _TrackBCAL_DeltaPhi_piminuscand[0]; TrackBCAL_DeltaZ_piminuscand=_TrackBCAL_DeltaZ_piminuscand[0]; dEdx_ST_piminuscand=_dEdx_ST_piminuscand[0]; dEdx_CDC_piminuscand=_dEdx_CDC_piminuscand[0];

   FOM_pipluscand = _FOM_pipluscand[0]; NDF_pipluscand = _NDF_pipluscand[0]; ChiSq_pipluscand = _ChiSq_pipluscand[0];
piplus_dedx_dc_NDF = _piplus_dedx_dc_NDF[0]; piplus_dedx_dc_ChiSq = _piplus_dedx_dc_ChiSq[0];  Beta_Timing_pipluscand = _Beta_Timing_pipluscand[0]; ChiSq_Timing_pipluscand = _ChiSq_Timing_pipluscand[0];
dEdx_TOF_pipluscand = _dEdx_TOF_pipluscand[0]; Energy_BCAL_pipluscand = _Energy_BCAL_pipluscand[0]; Energy_BCALPreshower_pipluscand=_Energy_BCALPreshower_pipluscand[0]; SigLong_BCAL_pipluscand =_SigLong_BCAL_pipluscand[0]; SigTheta_BCAL_pipluscand = _SigTheta_BCAL_pipluscand[0]; RMSTime_BCAL_pipluscand = _RMSTime_BCAL_pipluscand[0];  E_BCAL_pipluscand = _E_BCAL_pipluscand[0];  Energy_FCAL_pipluscand=_Energy_FCAL_pipluscand[0]; E1E9_FCAL_pipluscand=_E1E9_FCAL_pipluscand[0]; E9E25_FCAL_pipluscand=_E9E25_FCAL_pipluscand[0]; SumU_FCAL_pipluscand=_SumU_FCAL_pipluscand[0]; SumV_FCAL_pipluscand=_SumV_FCAL_pipluscand[0]; TrackBCAL_DeltaPhi_pipluscand = _TrackBCAL_DeltaPhi_pipluscand[0]; TrackBCAL_DeltaZ_pipluscand=_TrackBCAL_DeltaZ_pipluscand[0]; dEdx_ST_pipluscand=_dEdx_ST_pipluscand[0]; dEdx_CDC_pipluscand=_dEdx_CDC_pipluscand[0];  


   FOM_protcand = _FOM_protcand[0]; NDF_protcand = _NDF_protcand[0]; ChiSq_protcand = _ChiSq_protcand[0];
prot_dedx_dc_NDF = _prot_dedx_dc_NDF[0]; prot_dedx_dc_ChiSq = _prot_dedx_dc_ChiSq[0];  Beta_Timing_protcand = _Beta_Timing_protcand[0]; ChiSq_Timing_protcand = _ChiSq_Timing_protcand[0];
   dEdx_TOF_protcand = _dEdx_TOF_protcand[0]; Energy_BCAL_protcand = _Energy_BCAL_protcand[0]; Energy_BCALPreshower_protcand=_Energy_BCALPreshower_protcand[0]; SigLong_BCAL_protcand =_SigLong_BCAL_protcand[0]; SigTheta_BCAL_protcand = _SigTheta_BCAL_protcand[0]; RMSTime_BCAL_protcand = _RMSTime_BCAL_protcand[0];  E_BCAL_protcand = _E_BCAL_protcand[0];  Energy_FCAL_protcand=_Energy_FCAL_protcand[0]; E1E9_FCAL_protcand=_E1E9_FCAL_protcand[0]; E9E25_FCAL_protcand=_E9E25_FCAL_protcand[0]; SumU_FCAL_protcand=_SumU_FCAL_protcand[0]; SumV_FCAL_protcand=_SumV_FCAL_protcand[0]; TrackBCAL_DeltaPhi_protcand = _TrackBCAL_DeltaPhi_protcand[0]; TrackBCAL_DeltaZ_protcand=_TrackBCAL_DeltaZ_protcand[0]; dEdx_ST_protcand=_dEdx_ST_protcand[0]; dEdx_CDC_protcand=_dEdx_CDC_protcand[0];
   locDeltaBCAL_piplus = _locDeltaBCAL_piplus[0]; locDeltaTOF_piplus = _locDeltaTOF_piplus[0];
locDeltaBCAL_piminus = _locDeltaBCAL_piminus[0]; locDeltaTOF_piminus = _locDeltaTOF_piminus[0];
locDeltaBCAL_prot = _locDeltaBCAL_prot[0]; locDeltaTOF_prot = _locDeltaTOF_prot[0];

    hypnumber = 0;
 
    _numberruns = numberruns;   
//    if (_nHyp>0){
      
      if  (   (hyp1rho.M()<massmin1 ||  hyp1rho.M()>massmax1
	    || fabs(hyp1rho.Phi()-hyp1p.Phi())*RadToDeg < phimin ||  fabs(hyp1rho.Phi()-hyp1p.Phi())*RadToDeg > phimax)
            || ( dEdx_CDC_protcand*1e6 < exp(-4*prelPmag_p + 2.25) +1 && dEdx_CDC_protcand!=0  && preltheta_p*RadToDeg>11) 
	    || ( dEdx_CDC_pipluscand*1e6 > exp(-7.0*prelPmag_pip + 3.0) + 6.2 && dEdx_CDC_pipluscand!=0 && preltheta_pip*RadToDeg>11)  
	    || ( fabs(locDeltaBCAL_prot)>1 &&  preltheta_p*RadToDeg>11 ) 
	    || ( fabs(locDeltaTOF_prot)>1 &&  preltheta_p*RadToDeg<11 )
	    || ( fabs(locDeltaBCAL_piplus)>1 &&  preltheta_pip*RadToDeg>11 )
            || ( fabs(locDeltaTOF_piplus)>1 &&  preltheta_pip*RadToDeg<11 )
){
	
        TLorentzVector hyp2pim(_pX_piminus[0],_pY_piminus[0],_pZ_piminus[0], _E_piminus[0]);
	double Ep = sqrt((_pX_piplus[0]*_pX_piplus[0] + _pY_piplus[0]*_pY_piplus[0] + _pZ_piplus[0]*_pZ_piplus[0]) + mp*mp); 
	TLorentzVector hyp2p(_pX_piplus[0],_pY_piplus[0],_pZ_piplus[0], Ep);
	double Epip = sqrt((_pX_proton[0]*_pX_proton[0] + _pY_proton[0]*_pY_proton[0] + _pZ_proton[0]*_pZ_proton[0]) + mpi*mpi); 
	TLorentzVector hyp2pip(_pX_proton[0],_pY_proton[0],_pZ_proton[0], Epip);

	TLorentzVector hyp2rho = hyp2pim + hyp2pip;
	preltheta_p = hyp2p.Theta(); preltheta_pip = hyp2pip.Theta();
	prelPmag_p  = hyp2p.Vect().Mag();  prelPmag_pip = hyp2pip.Vect().Mag();

	pX_pim = _pX_piminus[0]; pY_pim = _pY_piminus[0]; pZ_pim = _pZ_piminus[0]; E_pim = _E_piminus[0];
        pX_pip = _pX_proton[0]; pY_pip = _pY_proton[0]; pZ_pip = _pZ_proton[0]; E_pip = Epip;
        pX_p = _pX_piplus[0]; pY_p = _pY_piplus[0]; pZ_p = _pZ_piplus[0]; E_p = Ep;
	
	X_vtx = _X_vertex[0]; Y_vtx = _Y_vertex[0]; Z_vtx = _Z_vertex[0];	
	t0_protcand = _t_piplus[0]; CL = _CLKinFit[0]; t_vtx = _T_vertex[0];
	commontime = _Common_Time[0]; 

	FOM_piminuscand = _FOM_piminuscand[0]; NDF_piminuscand = _NDF_piminuscand[0]; ChiSq_piminuscand = _ChiSq_piminuscand[0];
	piminus_dedx_dc_NDF = _piminus_dedx_dc_NDF[0]; piminus_dedx_dc_ChiSq = _piminus_dedx_dc_ChiSq[0];  Beta_Timing_piminuscand = _Beta_Timing_piminuscand[0]; ChiSq_Timing_piminuscand = _ChiSq_Timing_piminuscand[0];
	dEdx_TOF_piminuscand = _dEdx_TOF_piminuscand[0]; Energy_BCAL_piminuscand = _Energy_BCAL_piminuscand[0]; Energy_BCALPreshower_piminuscand=_Energy_BCALPreshower_piminuscand[0]; SigLong_BCAL_piminuscand =_SigLong_BCAL_piminuscand[0]; SigTheta_BCAL_piminuscand = _SigTheta_BCAL_piminuscand[0]; RMSTime_BCAL_piminuscand = _RMSTime_BCAL_piminuscand[0];  E_BCAL_piminuscand = _E_BCAL_piminuscand[0];  Energy_FCAL_piminuscand=_Energy_FCAL_piminuscand[0]; E1E9_FCAL_piminuscand=_E1E9_FCAL_piminuscand[0]; E9E25_FCAL_piminuscand=_E9E25_FCAL_piminuscand[0]; SumU_FCAL_piminuscand=_SumU_FCAL_piminuscand[0]; SumV_FCAL_piminuscand=_SumV_FCAL_piminuscand[0]; TrackBCAL_DeltaPhi_piminuscand = _TrackBCAL_DeltaPhi_piminuscand[0]; TrackBCAL_DeltaZ_piminuscand=_TrackBCAL_DeltaZ_piminuscand[0]; dEdx_ST_piminuscand=_dEdx_ST_piminuscand[0]; dEdx_CDC_piminuscand=_dEdx_CDC_piminuscand[0];

	FOM_pipluscand = _FOM_protcand[0]; NDF_pipluscand = _NDF_protcand[0]; ChiSq_pipluscand = _ChiSq_protcand[0];
  piplus_dedx_dc_NDF = _prot_dedx_dc_NDF[0]; piplus_dedx_dc_ChiSq = _prot_dedx_dc_ChiSq[0];  Beta_Timing_pipluscand = _Beta_Timing_protcand[0]; ChiSq_Timing_pipluscand = _ChiSq_Timing_protcand[0];
  dEdx_TOF_pipluscand = _dEdx_TOF_protcand[0]; Energy_BCAL_pipluscand = _Energy_BCAL_protcand[0]; Energy_BCALPreshower_pipluscand=_Energy_BCALPreshower_protcand[0]; SigLong_BCAL_pipluscand =_SigLong_BCAL_protcand[0]; 
  SigTheta_BCAL_pipluscand = _SigTheta_BCAL_protcand[0]; RMSTime_BCAL_pipluscand = _RMSTime_BCAL_protcand[0];  E_BCAL_pipluscand = _E_BCAL_protcand[0];  Energy_FCAL_pipluscand=_Energy_FCAL_protcand[0]; 
  E1E9_FCAL_pipluscand=_E1E9_FCAL_protcand[0]; E9E25_FCAL_pipluscand=_E9E25_FCAL_protcand[0]; SumU_FCAL_pipluscand=_SumU_FCAL_protcand[0]; SumV_FCAL_pipluscand=_SumV_FCAL_protcand[0]; 
  TrackBCAL_DeltaPhi_pipluscand = _TrackBCAL_DeltaPhi_protcand[0]; TrackBCAL_DeltaZ_pipluscand=_TrackBCAL_DeltaZ_protcand[0]; dEdx_ST_pipluscand=_dEdx_ST_protcand[0]; dEdx_CDC_pipluscand=_dEdx_CDC_protcand[0];

   FOM_protcand = _FOM_pipluscand[0]; NDF_protcand = _NDF_pipluscand[0]; ChiSq_protcand = _ChiSq_pipluscand[0];
  prot_dedx_dc_NDF = _piplus_dedx_dc_NDF[0]; prot_dedx_dc_ChiSq = _piplus_dedx_dc_ChiSq[0];  Beta_Timing_protcand = _Beta_Timing_pipluscand[0]; 
  ChiSq_Timing_protcand = _ChiSq_Timing_pipluscand[0];
  dEdx_TOF_protcand = _dEdx_TOF_pipluscand[0]; Energy_BCAL_pipluscand = _Energy_BCAL_pipluscand[0]; 
  Energy_BCALPreshower_protcand=_Energy_BCALPreshower_pipluscand[0]; SigLong_BCAL_protcand =_SigLong_BCAL_pipluscand[0]; SigTheta_BCAL_protcand = _SigTheta_BCAL_pipluscand[0]; 
  RMSTime_BCAL_protcand = _RMSTime_BCAL_pipluscand[0];  E_BCAL_protcand = _E_BCAL_pipluscand[0];  Energy_FCAL_protcand=_Energy_FCAL_pipluscand[0]; E1E9_FCAL_protcand=_E1E9_FCAL_pipluscand[0]; 
  E9E25_FCAL_protcand=_E9E25_FCAL_pipluscand[0]; SumU_FCAL_protcand=_SumU_FCAL_pipluscand[0]; SumV_FCAL_protcand=_SumV_FCAL_pipluscand[0]; TrackBCAL_DeltaPhi_protcand = _TrackBCAL_DeltaPhi_pipluscand[0]; 
  TrackBCAL_DeltaZ_protcand=_TrackBCAL_DeltaZ_pipluscand[0]; dEdx_ST_protcand=_dEdx_ST_pipluscand[0]; dEdx_CDC_protcand=_dEdx_CDC_pipluscand[0];

	locDeltaBCAL_piplus = _locDeltaBCAL_prot[0]; locDeltaTOF_piplus = _locDeltaTOF_prot[0];
	locDeltaBCAL_piminus = _locDeltaBCAL_piminus[0]; locDeltaTOF_piminus = _locDeltaTOF_piminus[0];
	locDeltaBCAL_prot = _locDeltaBCAL_piplus[0]; locDeltaTOF_prot = _locDeltaTOF_piplus[0];
	
	hypnumber = 2;
	
}
   // }
   
    if( dEdx_CDC_protcand*1e6 < exp(-4*prelPmag_p + 2.25) +1 && dEdx_CDC_protcand!=0) continue;
    if( dEdx_CDC_pipluscand*1e6 > exp(-7*prelPmag_pip + 3) +6.2 && dEdx_CDC_pipluscand!=0) continue; 
    //if( fabs(locDeltaBCAL_prot)>1 &&  preltheta_p*RadToDeg>11 ) continue;
    //if( fabs(locDeltaTOF_prot)>0.6 &&  preltheta_p*RadToDeg<11 ) continue;
    //if( fabs(locDeltaBCAL_piplus)>1 &&  preltheta_pip*RadToDeg>11 ) continue;
    //if( fabs(locDeltaTOF_piplus)>0.4 &&  preltheta_pip*RadToDeg<11 ) continue;

   //if( fabs(locDeltaBCAL_piplus)>0|| fabs(locDeltaTOF_piplus)>0.5) continue;
    TLorentzVector pPim(pX_pim,pY_pim,pZ_pim,E_pim);
    TLorentzVector pPip(pX_pip,pY_pip,pZ_pip,E_pip);
    TLorentzVector pP(pX_p,pY_p,pZ_p,E_p);
     
    TLorentzVector rho = pPim + pPip;
    TLorentzVector Pimp = pPim + pP; 
    TLorentzVector Pipp = pPip + pP;
    
    if( fabs(pP.M() - 0.938272)> 0.0001) continue; 
   
    if  (rho.M()<massmin1 || rho.M()>massmax1
	 ||  fabs(rho.Phi()   - pP.Phi())*RadToDeg < phimin ||  fabs(rho.Phi() - pP.Phi())*RadToDeg > phimax  ) continue;
    
    if ( (rho.E() + pP.E())<6 ) continue;
    
    if (CL< CL_cut) continue;
 
    mass_pimp = Pimp.M(); mass_pipp = Pipp.M();
    theta_pip = pPip.Theta(); theta_pim = pPim.Theta(); theta_p = pP.Theta();
    phi_pip = pPip.Phi(); phi_pim = pPim.Phi(); phi_p = pP.Phi();
    theta_rho = rho.Theta(); theta_pimp = Pimp.Theta(); theta_pipp = Pipp.Theta();
    phi_rho = rho.Phi(); phi_pimp = Pimp.Phi(); phi_pipp = Pipp.Phi();
    Pmag_pip = pPip.Vect().Mag(); Pmag_pim = pPim.Vect().Mag(); Pmag_p = pP.Vect().Mag();
    
    mass_rho  = rho.M();  
    s =  (pPim+pPip+pP)*(pPim+pPip+pP);
    
    //---------------- Information with Photon Energy
    
    InTimePhotons = OutTimePhotons = allphotons = 0;
  
    TLorentzVector Ptot	 = pPim + pPip + pP;
    InTimePhotons = 1;  
    TVector3 boost_rhop =  -Ptot.BoostVector();
    TLorentzVector pPim_CM  =  pPim;
    TLorentzVector pPip_CM  =  pPip;
    TLorentzVector pP_CM    =  pP;
    TLorentzVector tot_CM   =  Ptot;
    
    pPim_CM.Boost(boost_rhop);
    pPip_CM.Boost(boost_rhop);
    pP_CM.Boost(boost_rhop);
    tot_CM.Boost(boost_rhop);
    
    Double_t dEt = 0;  
    
    TLorentzVector Ep(0,0,_bmE[0],_bmE[0]);
      
    dEt = (Ep + target - rho - pP).E();
    if( fabs(dEt) >1. ) continue;
    if( _bmE[0] < 6. && _bmE[0] > 10.5) continue;
    TLorentzVector v4Rho_CM = rho;   
    TLorentzVector v4Beam_CM = Ep;
    v4Rho_CM.Boost(boost_rhop);
    v4Beam_CM.Boost(boost_rhop);
    double cosTheta_CM = cos(v4Rho_CM.Vect().Angle(v4Beam_CM.Vect())); 
  
    _w = weight; 
    TLorentzVector Pmissvec = pP + rho - Ep;
      
    Double_t deltat = t_vtx - (_bmtime[0]  + (Z_vtx - 65. )/29.9792458);
    E_photon[0] = _bmE[0];
    massinv2[0] = (rho + pP - Ep).M2(); 
    kmiss[0] = mp*sqrt(((Pmissvec.Perp2() + mp*mp)/((Pmissvec.E() - Pmissvec.Pz())*(2*mp-(Pmissvec.E() - Pmissvec.Pz()))))-1);
    Mmiss[0] = sqrt(pow((2*mp -(Pmissvec.E() - Pmissvec.Pz())),2) - pow((Pmissvec.Vect()).Mag(),2));
    dE[0] =dEt;                 
    t[0] = -(Ep-rho)*(Ep-rho);
    u[0] = -(Ep-pP)*(Ep-pP);
    Emiss[0] = Pmissvec.E();
    Pmiss[0] =(Pmissvec.Vect()).Mag();
    Pmissmin[0] = Pmissvec.E() - Pmissvec.Pz();
    Pmissplus[0] = Pmissvec.E() + Pmissvec.Pz();
    Pmissperp[0] = Pmissvec.Perp();     
    q[0] = sqrt(pPim_CM.Z()*pPim_CM.Z() + pPip_CM.Z()*pPip_CM.Z() +  pP_CM.Z()*pP_CM.Z());
    omega[0]  = TMath::ASin((pPim_CM.Z())/(sqrt(2./3.)*q[0]));
    xcm[0] = q[0]*TMath::Cos(omega[0]);
    ycm[0] = q[0]*TMath::Sin(omega[0]);         
    dt[0] = deltat;
    cosThetaCM[0] = cosTheta_CM;        

    //mass_rho[0] = rho.M(); 

    /* 
    for (Int_t j=0; j<_nPhotonCandidates; j++ ){
      TLorentzVector Ep(0,0,_bmE[j],_bmE[j]);
      
      dEt = (Ep + target - rho - pP).E();
      if( fabs(dEt) >1. ) continue;
      if( _bmE[j] < 6. && _bmE[j] > 10.5) continue;
      // CM Scattering Angle
      TLorentzVector v4Rho_CM = rho;   
      TLorentzVector v4Beam_CM = Ep;
      v4Rho_CM.Boost(boost_rhop);
      v4Beam_CM.Boost(boost_rhop);
      double cosTheta_CM = cos(v4Rho_CM.Vect().Angle(v4Beam_CM.Vect()));      
 
      _w = weight; 
      TLorentzVector Pmissvec = pP + rho - Ep;
      
      Double_t deltat = t_vtx - (_bmtime[j]  + (Z_vtx - 65. )/29.9792458); 
      
      if( fabs(deltat) < 2.){
	E_photon[InTimePhotons] = _bmE[j];
	massinv2[InTimePhotons] = (rho + pP - Ep).M2();	
	kmiss[InTimePhotons] = mp*sqrt(((Pmissvec.Perp2() + mp*mp)/((Pmissvec.E() - Pmissvec.Pz())*(2*mp-(Pmissvec.E() - Pmissvec.Pz()))))-1);
	Mmiss[InTimePhotons] = sqrt(pow((2*mp -(Pmissvec.E() - Pmissvec.Pz())),2) - pow((Pmissvec.Vect()).Mag(),2));
	dE[InTimePhotons] =dEt;  		
	t[InTimePhotons] = -(Ep-rho)*(Ep-rho);
	u[InTimePhotons] = -(Ep-pP)*(Ep-pP);
	Emiss[InTimePhotons] = Pmissvec.E();
	Pmiss[InTimePhotons] =(Pmissvec.Vect()).Mag();
	Pmissmin[InTimePhotons] = Pmissvec.E() - Pmissvec.Pz();
	Pmissplus[InTimePhotons] = Pmissvec.E() + Pmissvec.Pz();
	Pmissperp[InTimePhotons] = Pmissvec.Perp();	
	q[InTimePhotons] = sqrt(pPim_CM.Z()*pPim_CM.Z() + pPip_CM.Z()*pPip_CM.Z() +  pP_CM.Z()*pP_CM.Z());
	omega[InTimePhotons]  = TMath::ASin((pPim_CM.Z())/(sqrt(2./3.)*q[InTimePhotons]));
	xcm[InTimePhotons] = q[InTimePhotons]*TMath::Cos(omega[InTimePhotons]);
	ycm[InTimePhotons] = q[InTimePhotons]*TMath::Sin(omega[InTimePhotons]);		
	dt[InTimePhotons] = deltat;
	cosThetaCM[InTimePhotons] = cosTheta_CM; 	
	//mass_rho[InTimePhotons] = rho.M(); 
	InTimePhotons++;
      }
      
      
      if( fabs(deltat) > 6. && fabs(deltat)< 18. ){
	massinv2bgd[OutTimePhotons] = (rho + pP - Ep).M2();	
	kmissbgd[OutTimePhotons] = mp*sqrt(((Pmissvec.Perp2() + mp*mp)/((Pmissvec.E() - Pmissvec.Pz())*(2*mp-(Pmissvec.E() - Pmissvec.Pz()))))-0);
	Mmissbgd[OutTimePhotons] = sqrt(pow((2*mp -(Pmissvec.E() - Pmissvec.Pz())),2) - pow((Pmissvec.Vect()).Mag(),2));
	dEbgd[OutTimePhotons] = dEt; 
	tbgd[OutTimePhotons] = -(Ep-rho)*(Ep-rho);
	ubgd[OutTimePhotons] = -(Ep-pP)*(Ep-pP);
	Emissbgd[OutTimePhotons] = Pmissvec.E();
	Pmissbgd[OutTimePhotons] =(Pmissvec.Vect()).Mag();
	Pmissminbgd[OutTimePhotons] = Pmissvec.E() - Pmissvec.Pz();
	Pmissplusbgd[OutTimePhotons] = Pmissvec.E() + Pmissvec.Pz();
	Pmissperpbgd[OutTimePhotons] = Pmissvec.Perp();	
	qbgd[OutTimePhotons] = sqrt(pPim_CM.Z()*pPim_CM.Z() + pPip_CM.Z()*pPip_CM.Z() +  pP_CM.Z()*pP_CM.Z());
	omegabgd[OutTimePhotons]  = TMath::ASin((pPim_CM.Z())/(sqrt(2./3.)*qbgd[OutTimePhotons]));
	xcmbgd[OutTimePhotons] = qbgd[OutTimePhotons]*TMath::Cos(omegabgd[OutTimePhotons]);
	ycmbgd[OutTimePhotons] = qbgd[OutTimePhotons]*TMath::Sin(omegabgd[OutTimePhotons]);		
	dtbgd[OutTimePhotons] = deltat;
	cosThetaCMbgd[OutTimePhotons] = cosTheta_CM; 
	//mass_rhobgd[OutTimePhotons] = rho.M();	
	OutTimePhotons++;
      }
      
      
      if( fabs(deltat)< 18. ){
	dttot[allphotons] = deltat;
	allphotons++;
      }
      
    }
   
    if (Z_vtx>zmin && Z_vtx<zmax ) EvTree->Fill(); // Change this cut for production data*/
 
     EvTree->Fill();

  }
  
  cout << tentries << endl;
  
  EvTree->Write();
  MyFile->Close();
  
}
