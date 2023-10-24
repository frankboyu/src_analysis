#include "phi_selector.h"

const double shower_number_max_cut     = 10;
const double hypothesis_number_min_cut = 1;
const double photon_number_min_cut     = 1;

const double confidence_level_min_cut  = 0.0001;
const double phi_mass_min_cut          = 1.00;    
const double phi_mass_max_cut          = 1.04; 
const double coplanarity_min_cut       = 170;   
const double coplanarity_max_cut       = 190; 
const double vertex_z_min_cut          = 51;   
const double vertex_z_max_cut          = 79;
const double CDC_kaon_max_cut[3]       = {-7.0, 3.0, 6.2};
const double CDC_proton_min_cut[3]     = {-4.0, 2.25, 1.0};
const double BCAL_kaon_max_cut         = 0.75;
const double BCAL_proton_max_cut       = 1.0;
const double TOF_kaon_max_cut          = 0.3;
const double TOF_proton_max_cut        = 0.6;

const double intime_photon_max_cut     = 2.0;
const double outtime_photon_min_cut    = 6.0;
const double outtime_photon_max_cut    = 18.0;

const double energy_balance_max_cut    = 1.0;
const double photon_energy_min_cut     = 6.0; 
const double photon_energy_max_cut     = 10.5;

void phi_selector(int runnumb) {

  TString infile_folder = "/volatile/halld/home/sns/1p2K_hddm/0" + std::to_string(runnumb) + "/";	
  TString outputfile = "results/phi_candidates_" + std::to_string(runnumb) + ".root";
  //TString infile_folder = "results/kinfit/";	
  //TString outputfile = "results/phi_candidates/phi_candidates_4He_MF_100k.root";
  
  const char * dataPath = infile_folder;
  vector<string> inputFiles;
  DIR *dirp;
  dirp = opendir(dataPath);
  regex extension ("(.*)(.root)");
  if (dirp){
    string dataPathString = dataPath;
    struct dirent *directory;
    while ((directory = readdir(dirp))){
      string fileName = directory->d_name;
      string fullName = dataPathString + fileName;
      if (!regex_match(fileName,regex("(tree_1p2k)(.*)")))
	continue;	  
      inputFiles.push_back(fullName);
    }
    closedir(dirp);
  }
  else if (ifstream(dataPath)){
    string fileName = dataPath;  
    if (!regex_match(fileName,extension)){
      cout << "Invalid input. Exiting...\n";
      return;
    }
    inputFiles.push_back(fileName);
  }
  else{
    cout << "Input data not found. Exiting...\n";
    return;
  }

  TChain *T = new TChain("tree_1p2k");
  sort(inputFiles.begin(),inputFiles.end());
  for (const auto& fullName: inputFiles){
    cout << fullName << "\n"; 
    TFile *f = new TFile(fullName.c_str());      
    if (f->GetListOfKeys()->Contains("tree_1p2k")){
      T->Add(Form("%s",fullName.c_str()));
    }
    else{
      cout << "Data tree not found, abandoning file\n";
      continue;
    }
  }
  
  TFile *MyFile = new TFile(Form("%s",outputfile.Data()),"RECREATE");
  TTree *EvTree = new TTree("EvTree","EvTree");
  SetBranchAddressesTree(T);	
  SetBranchAddressesOutTree(EvTree);  
  
  Long64_t tentries = T->GetEntries();
  cout << "Number of events in trees: " << tentries << endl;  
  int pass = 0;
  
  for(Long64_t i=0; i<tentries; i++){

    if (i%1000000==0) cout << "events processed: " << int(i/1000000) << "M" << endl; 

    T->GetEntry(i);
    
    shower_number     = _nShower;
    hypothesis_number = _nHyp;
    photon_number     = _nPhotonCandidates;

    if (shower_number     > shower_number_max_cut     || 
	hypothesis_number < hypothesis_number_min_cut ||
	photon_number     < photon_number_min_cut    
	)
      continue; 

    pass = 0;
    for(Int_t j=0; j<hypothesis_number; j++){
      pX_km = _pX_kminus[j]; pY_km = _pY_kminus[j]; pZ_km = _pZ_kminus[j]; E_km = _E_kminus[j];
      pX_kp = _pX_kplus[j];  pY_kp = _pY_kplus[j];  pZ_kp = _pZ_kplus[j];  E_kp = _E_kplus[j];
      pX_p  = _pX_proton[j]; pY_p  = _pY_proton[j]; pZ_p  = _pZ_proton[j]; E_p  = _E_proton[j];
    
      pkm.SetPxPyPzE(pX_km,pY_km,pZ_km,E_km);
      pkp.SetPxPyPzE(pX_kp,pY_kp,pZ_kp,E_kp);
      pP.SetPxPyPzE(pX_p,pY_p,pZ_p,E_p);
      phi = pkm + pkp;
      kmp = pkm + pP; 
      kpp = pkp + pP;    
      ppim.SetXYZM(pX_km,pY_km,pZ_km,mpi);
      ppip.SetXYZM(pX_kp,pY_kp,pZ_kp,mpi);
      rho = ppim + ppip;

      theta_km  = pkm.Theta(); phi_km  = pkm.Phi(); Pmag_km  = pkm.Vect().Mag();
      theta_kp  = pkp.Theta(); phi_kp  = pkp.Phi(); Pmag_kp  = pkp.Vect().Mag();
      theta_p   = pP.Theta();  phi_p   = pP.Phi();  Pmag_p   = pP.Vect().Mag();    
      theta_phi = phi.Theta(); phi_phi = phi.Phi(); mass_phi = phi.M();
      theta_kmp = kmp.Theta(); phi_kmp = kmp.Phi(); mass_kmp = kmp.M();
      theta_kpp = kpp.Theta(); phi_kpp = kpp.Phi(); mass_kpp = kpp.M();
      s = (pkm+pkp+pP)*(pkm+pkp+pP);
      
      mass_rho = rho.M();
    
      X_vtx       = _X_vertex[j]; 
      Y_vtx       = _Y_vertex[j]; 
      Z_vtx       = _Z_vertex[j]; 
      T_vtx       = _T_vertex[j];    
      CL          = _CLKinFit[j]; 
      t0_protcand = _t_proton[j]; 
      commontime  = _Common_Time[j]; 
  
      FOM_kminuscand                  = _FOM_kminuscand[j]; 
      NDF_kminuscand                  = _NDF_kminuscand[j]; 
      ChiSq_kminuscand                = _ChiSq_kminuscand[j];
      kminus_dedx_dc_NDF              = _kminus_dedx_dc_NDF[j]; 
      kminus_dedx_dc_ChiSq            = _kminus_dedx_dc_ChiSq[j];  
      Beta_Timing_kminuscand          = _Beta_Timing_kminuscand[j]; 
      ChiSq_Timing_kminuscand         = _ChiSq_Timing_kminuscand[j];
      dEdx_TOF_kminuscand             = _dEdx_TOF_kminuscand[j]; 
      Energy_BCAL_kminuscand          = _Energy_BCAL_kminuscand[j]; 
      Energy_BCALPreshower_kminuscand = _Energy_BCALPreshower_kminuscand[j]; 
      SigLong_BCAL_kminuscand         = _SigLong_BCAL_kminuscand[j]; 
      SigTheta_BCAL_kminuscand        = _SigTheta_BCAL_kminuscand[j]; 
      RMSTime_BCAL_kminuscand         = _RMSTime_BCAL_kminuscand[j];  
      E_BCAL_kminuscand               = _E_BCAL_kminuscand[j];  
      Energy_FCAL_kminuscand          = _Energy_FCAL_kminuscand[j]; 
      E1E9_FCAL_kminuscand            = _E1E9_FCAL_kminuscand[j]; 
      E9E25_FCAL_kminuscand           = _E9E25_FCAL_kminuscand[j]; 
      SumU_FCAL_kminuscand            = _SumU_FCAL_kminuscand[j]; 
      SumV_FCAL_kminuscand            = _SumV_FCAL_kminuscand[j]; 
      TrackBCAL_DeltaPhi_kminuscand   = _TrackBCAL_DeltaPhi_kminuscand[j]; 
      TrackBCAL_DeltaZ_kminuscand     = _TrackBCAL_DeltaZ_kminuscand[j]; 
      dEdx_ST_kminuscand              = _dEdx_ST_kminuscand[j]; 
      dEdx_CDC_kminuscand             = _dEdx_CDC_kminuscand[j];
      locDeltaBCAL_kminus             = _locDeltaBCAL_kminus[j]; 
      locDeltaTOF_kminus              = _locDeltaTOF_kminus[j];

      FOM_kpluscand                   = _FOM_kpluscand[j]; 
      NDF_kpluscand                   = _NDF_kpluscand[j]; 
      ChiSq_kpluscand                 = _ChiSq_kpluscand[j];
      kplus_dedx_dc_NDF               = _kplus_dedx_dc_NDF[j]; 
      kplus_dedx_dc_ChiSq             = _kplus_dedx_dc_ChiSq[j];  
      Beta_Timing_kpluscand           = _Beta_Timing_kpluscand[j]; 
      ChiSq_Timing_kpluscand          = _ChiSq_Timing_kpluscand[j];
      dEdx_TOF_kpluscand              = _dEdx_TOF_kpluscand[j]; 
      Energy_BCAL_kpluscand           = _Energy_BCAL_kpluscand[j]; 
      Energy_BCALPreshower_kpluscand  = _Energy_BCALPreshower_kpluscand[j]; 
      SigLong_BCAL_kpluscand          = _SigLong_BCAL_kpluscand[j]; 
      SigTheta_BCAL_kpluscand         = _SigTheta_BCAL_kpluscand[j]; 
      RMSTime_BCAL_kpluscand          = _RMSTime_BCAL_kpluscand[j];  
      E_BCAL_kpluscand                = _E_BCAL_kpluscand[j];  
      Energy_FCAL_kpluscand           = _Energy_FCAL_kpluscand[j]; 
      E1E9_FCAL_kpluscand             = _E1E9_FCAL_kpluscand[j]; 
      E9E25_FCAL_kpluscand            = _E9E25_FCAL_kpluscand[j]; 
      SumU_FCAL_kpluscand             = _SumU_FCAL_kpluscand[j]; 
      SumV_FCAL_kpluscand             = _SumV_FCAL_kpluscand[j]; 
      TrackBCAL_DeltaPhi_kpluscand    = _TrackBCAL_DeltaPhi_kpluscand[j]; 
      TrackBCAL_DeltaZ_kpluscand      = _TrackBCAL_DeltaZ_kpluscand[j]; 
      dEdx_ST_kpluscand               = _dEdx_ST_kpluscand[j]; 
      dEdx_CDC_kpluscand              = _dEdx_CDC_kpluscand[j];  
      locDeltaBCAL_kplus              = _locDeltaBCAL_kplus[j]; 
      locDeltaTOF_kplus               = _locDeltaTOF_kplus[j];

      FOM_protcand                    = _FOM_protcand[j]; 
      NDF_protcand                    = _NDF_protcand[j]; 
      ChiSq_protcand                  = _ChiSq_protcand[j];
      prot_dedx_dc_NDF                = _prot_dedx_dc_NDF[j]; 
      prot_dedx_dc_ChiSq              = _prot_dedx_dc_ChiSq[j];  
      Beta_Timing_protcand            = _Beta_Timing_protcand[j]; 
      ChiSq_Timing_protcand           = _ChiSq_Timing_protcand[j];
      dEdx_TOF_protcand               = _dEdx_TOF_protcand[j]; 
      Energy_BCAL_protcand            = _Energy_BCAL_protcand[j]; 
      Energy_BCALPreshower_protcand   = _Energy_BCALPreshower_protcand[j]; 
      SigLong_BCAL_protcand           = _SigLong_BCAL_protcand[j]; 
      SigTheta_BCAL_protcand          = _SigTheta_BCAL_protcand[j]; 
      RMSTime_BCAL_protcand           = _RMSTime_BCAL_protcand[j];  
      E_BCAL_protcand                 = _E_BCAL_protcand[j];  
      Energy_FCAL_protcand            = _Energy_FCAL_protcand[j]; 
      E1E9_FCAL_protcand              = _E1E9_FCAL_protcand[j]; 
      E9E25_FCAL_protcand             = _E9E25_FCAL_protcand[j]; 
      SumU_FCAL_protcand              = _SumU_FCAL_protcand[j]; 
      SumV_FCAL_protcand              = _SumV_FCAL_protcand[j]; 
      TrackBCAL_DeltaPhi_protcand     = _TrackBCAL_DeltaPhi_protcand[j]; 
      TrackBCAL_DeltaZ_protcand       = _TrackBCAL_DeltaZ_protcand[j]; 
      dEdx_ST_protcand                = _dEdx_ST_protcand[j]; 
      dEdx_CDC_protcand               = _dEdx_CDC_protcand[j];
      locDeltaBCAL_prot               = _locDeltaBCAL_prot[j]; 
      locDeltaTOF_prot                = _locDeltaTOF_prot[j];

      if(CL                                   > confidence_level_min_cut
	 && phi.M()                           > phi_mass_min_cut  
	 && phi.M()                           < phi_mass_max_cut                                                                                       
	 && fabs(phi.Phi()-pP.Phi())*RadToDeg > coplanarity_min_cut
	 && fabs(phi.Phi()-pP.Phi())*RadToDeg < coplanarity_max_cut
	 && Z_vtx                             > vertex_z_min_cut
	 && Z_vtx                             < vertex_z_max_cut
	 && (dEdx_CDC_kminuscand*1e6          < exp(CDC_kaon_max_cut[0]*Pmag_kp+CDC_kaon_max_cut[1])+CDC_kaon_max_cut[2]      || !dEdx_CDC_kminuscand)
	 && (dEdx_CDC_kpluscand*1e6           < exp(CDC_kaon_max_cut[0]*Pmag_kp+CDC_kaon_max_cut[1])+CDC_kaon_max_cut[2]      || !dEdx_CDC_kpluscand)
	 && (dEdx_CDC_protcand*1e6            > exp(CDC_proton_min_cut[0]*Pmag_p+CDC_proton_min_cut[1])+CDC_proton_min_cut[2] || !dEdx_CDC_protcand)   
	 && (fabs(locDeltaBCAL_kminus)        < BCAL_kaon_max_cut   ||  locDeltaBCAL_kminus == 999.0)
	 && (fabs(locDeltaBCAL_kplus)         < BCAL_kaon_max_cut   ||  locDeltaBCAL_kplus  == 999.0)
	 && (fabs(locDeltaBCAL_prot)          < BCAL_proton_max_cut ||  locDeltaBCAL_prot   == 999.0)
	 && (fabs(locDeltaTOF_kminus)         < TOF_kaon_max_cut    ||  locDeltaTOF_kminus  == 999.0)
	 && (fabs(locDeltaTOF_kplus)          < TOF_kaon_max_cut    ||  locDeltaTOF_kplus   == 999.0)
	 && (fabs(locDeltaTOF_prot)           < TOF_proton_max_cut  ||  locDeltaTOF_prot    == 999.0)
	 )
	pass++;
    }
    if(pass == 0) continue;

    //photon info    
    allphotons = 0;
    for (Int_t j=0; j<photon_number; j++){
      TVector3 boost_phip =  -(pkm+pkp+pP).BoostVector();      
      TLorentzVector Ep(0,0,_bmE[j],_bmE[j]);
      TLorentzVector target(0.0, 0.0, 0.0, mp);
      TLorentzVector Pmissvec = pP + phi - Ep;
      TLorentzVector pkm_CM     = pkm;
      TLorentzVector pkp_CM     = pkp;
      TLorentzVector pP_CM      = pP;
      TLorentzVector pphi_CM    = phi;
      TLorentzVector pphoton_CM = Ep;
      pkm_CM.Boost(boost_phip);
      pkp_CM.Boost(boost_phip);
      pP_CM.Boost(boost_phip);
      pphi_CM.Boost(boost_phip);   
      pphoton_CM.Boost(boost_phip);

      E_photon[allphotons]       = _bmE[j];
      massinv2[allphotons]       = (phi + pP - Ep).M2();	
      kmiss[allphotons]          = mp*sqrt(((Pmissvec.Perp2() + mp*mp)/((Pmissvec.E() - Pmissvec.Pz())*(2*mp-(Pmissvec.E() - Pmissvec.Pz()))))-1);
      Mmiss[allphotons]          = sqrt(pow((2*mp -(Pmissvec.E() - Pmissvec.Pz())),2) - pow((Pmissvec.Vect()).Mag(),2));
      dE[allphotons]             = (Ep + target - phi - pP).E();  		
      t[allphotons]              = -(Ep-phi)*(Ep-phi);
      u[allphotons]              = -(Ep-pP)*(Ep-pP);
      Emiss[allphotons]          = Pmissvec.E();
      Pmiss[allphotons]          = (Pmissvec.Vect()).Mag();
      Pmissmin[allphotons]       = Pmissvec.E() - Pmissvec.Pz();
      Pmissplus[allphotons]      = Pmissvec.E() + Pmissvec.Pz();
      Pmissperp[allphotons]      = Pmissvec.Perp();	
      q[allphotons]              = sqrt(pkm_CM.Z()*pkm_CM.Z() + pkp_CM.Z()*pkp_CM.Z() +  pP_CM.Z()*pP_CM.Z());
      omega[allphotons]          = TMath::ASin((pkm_CM.Z())/(sqrt(2./3.)*q[allphotons]));
      xcm[allphotons]            = q[allphotons]*TMath::Cos(omega[allphotons]);
      ycm[allphotons]            = q[allphotons]*TMath::Sin(omega[allphotons]);		
      dt[allphotons]             = T_vtx - (_bmtime[j]  + (Z_vtx - 65. )/29.9792458);
      cosThetaCM[allphotons]     = cos(pphi_CM.Vect().Angle(pphoton_CM.Vect())); 	
     
      if (fabs(dt[allphotons]) < intime_photon_max_cut) 
        weight[allphotons] = 1;
      else if (fabs(dt[allphotons]) > outtime_photon_min_cut && fabs(dt[allphotons]) < outtime_photon_max_cut)   
	weight[allphotons] = -4.0/(2.0*(outtime_photon_max_cut - outtime_photon_min_cut));
      else 
	continue;

      if(fabs(dE[allphotons]) < energy_balance_max_cut && 
	 E_photon[allphotons] > photon_energy_min_cut  &&
	 E_photon[allphotons] < photon_energy_max_cut
	 ) 
	allphotons++;
    }
    if (allphotons == 0) continue;
    
    EvTree->Fill();
  }
  EvTree->Write();
  MyFile->Close();
}
