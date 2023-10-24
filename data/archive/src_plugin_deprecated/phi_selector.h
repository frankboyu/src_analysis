#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector> 
#include <array>
#include <stdio.h> 
#include <dirent.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TMath.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TFile.h>
using namespace std;

const Double_t RadToDeg = 180.0 / 3.14159265358979323846;
const Double_t mp  = 0.938272;
const Double_t mk  = 0.493677;
const Double_t mpi = 0.139570;

//input tree
////hyperparameters
const int maxnHyp = 3;
const int maxnPhotonCandidates	= 100;

////event info
int _eventNumber;
int _L1TriggerBits;
int _nShower;
int _nHyp;
int _nPhotonCandidates;

////kinematic fit
double* _pX_kminus			  = new double[maxnHyp];
double* _pY_kminus			  = new double[maxnHyp];
double* _pZ_kminus			  = new double[maxnHyp];
double* _E_kminus			  = new double[maxnHyp];
double* _t_kminus			  = new double[maxnHyp];
double* _pX_kplus			  = new double[maxnHyp];
double* _pY_kplus			  = new double[maxnHyp];
double* _pZ_kplus			  = new double[maxnHyp];
double* _E_kplus			  = new double[maxnHyp];
double* _t_kplus			  = new double[maxnHyp];
double* _pX_proton			  = new double[maxnHyp];
double* _pY_proton			  = new double[maxnHyp];
double* _pZ_proton			  = new double[maxnHyp];
double* _E_proton			  = new double[maxnHyp];
double* _t_proton			  = new double[maxnHyp];
double* _X_vertex			  = new double[maxnHyp];
double* _Y_vertex			  = new double[maxnHyp];
double* _Z_vertex			  = new double[maxnHyp];
double* _T_vertex			  = new double[maxnHyp];
double* _CLKinFit			  = new double[maxnHyp];
double* _NDF				  = new double[maxnHyp];
double* _ChiSqFit			  = new double[maxnHyp];
double* _Common_Time	 	          = new double[maxnHyp];

////KMinus PID
double* _FOM_kminuscand                   = new double[maxnHyp];
double* _NDF_kminuscand                   = new double[maxnHyp];
double* _ChiSq_kminuscand                 = new double[maxnHyp];
double* _kminus_dedx_dc_NDF               = new double[maxnHyp];
double* _kminus_dedx_dc_ChiSq             = new double[maxnHyp];
double* _Beta_Timing_kminuscand           = new double[maxnHyp];
double* _ChiSq_Timing_kminuscand          = new double[maxnHyp];
double* _dEdx_TOF_kminuscand              = new double[maxnHyp];
double* _Energy_BCAL_kminuscand           = new double[maxnHyp];
double* _Energy_BCALPreshower_kminuscand  = new double[maxnHyp];
double* _SigLong_BCAL_kminuscand          = new double[maxnHyp];
double* _SigTheta_BCAL_kminuscand         = new double[maxnHyp];
double* _RMSTime_BCAL_kminuscand          = new double[maxnHyp];
double* _E_BCAL_kminuscand                = new double[maxnHyp];
double* _Energy_FCAL_kminuscand           = new double[maxnHyp];
double* _E1E9_FCAL_kminuscand             = new double[maxnHyp];
double* _E9E25_FCAL_kminuscand            = new double[maxnHyp];
double* _SumU_FCAL_kminuscand             = new double[maxnHyp];
double* _SumV_FCAL_kminuscand             = new double[maxnHyp];
double* _TrackBCAL_DeltaPhi_kminuscand    = new double[maxnHyp];
double* _TrackBCAL_DeltaZ_kminuscand      = new double[maxnHyp];
double* _dEdx_ST_kminuscand               = new double[maxnHyp];
double* _dEdx_CDC_kminuscand              = new double[maxnHyp];
double* _locDeltaBCAL_kminus              = new double[maxnHyp];
double* _locDeltaTOF_kminus               = new double[maxnHyp];

////KPlus PID
double* _FOM_kpluscand                    = new double[maxnHyp];
double* _NDF_kpluscand                    = new double[maxnHyp];
double* _ChiSq_kpluscand                  = new double[maxnHyp];
double* _kplus_dedx_dc_NDF                = new double[maxnHyp];
double* _kplus_dedx_dc_ChiSq              = new double[maxnHyp];
double* _Beta_Timing_kpluscand            = new double[maxnHyp];
double* _ChiSq_Timing_kpluscand           = new double[maxnHyp];
double* _dEdx_TOF_kpluscand               = new double[maxnHyp];
double* _Energy_BCAL_kpluscand            = new double[maxnHyp];
double* _Energy_BCALPreshower_kpluscand   = new double[maxnHyp];
double* _SigLong_BCAL_kpluscand           = new double[maxnHyp];
double* _SigTheta_BCAL_kpluscand          = new double[maxnHyp];
double* _RMSTime_BCAL_kpluscand           = new double[maxnHyp];
double* _E_BCAL_kpluscand                 = new double[maxnHyp];
double* _Energy_FCAL_kpluscand            = new double[maxnHyp];
double* _E1E9_FCAL_kpluscand              = new double[maxnHyp];
double* _E9E25_FCAL_kpluscand             = new double[maxnHyp];
double* _SumU_FCAL_kpluscand              = new double[maxnHyp];
double* _SumV_FCAL_kpluscand              = new double[maxnHyp];
double* _TrackBCAL_DeltaPhi_kpluscand     = new double[maxnHyp];
double* _TrackBCAL_DeltaZ_kpluscand       = new double[maxnHyp];
double* _dEdx_ST_kpluscand                = new double[maxnHyp];
double* _dEdx_CDC_kpluscand               = new double[maxnHyp];
double* _locDeltaBCAL_kplus               = new double[maxnHyp];
double* _locDeltaTOF_kplus                = new double[maxnHyp];

////proton PID
double* _FOM_protcand                     = new double[maxnHyp];
double* _NDF_protcand                     = new double[maxnHyp];
double* _ChiSq_protcand                   = new double[maxnHyp];
double* _prot_dedx_dc_NDF                 = new double[maxnHyp];
double* _prot_dedx_dc_ChiSq               = new double[maxnHyp];
double* _Beta_Timing_protcand             = new double[maxnHyp];
double* _ChiSq_Timing_protcand            = new double[maxnHyp];
double* _dEdx_TOF_protcand                = new double[maxnHyp];
double* _Energy_BCAL_protcand             = new double[maxnHyp];
double* _Energy_BCALPreshower_protcand    = new double[maxnHyp];
double* _SigLong_BCAL_protcand            = new double[maxnHyp];
double* _SigTheta_BCAL_protcand           = new double[maxnHyp];
double* _RMSTime_BCAL_protcand            = new double[maxnHyp];
double* _E_BCAL_protcand                  = new double[maxnHyp];
double* _Energy_FCAL_protcand             = new double[maxnHyp];
double* _E1E9_FCAL_protcand               = new double[maxnHyp];
double* _E9E25_FCAL_protcand              = new double[maxnHyp];
double* _SumU_FCAL_protcand               = new double[maxnHyp];
double* _SumV_FCAL_protcand               = new double[maxnHyp];
double* _TrackBCAL_DeltaPhi_protcand      = new double[maxnHyp];
double* _TrackBCAL_DeltaZ_protcand        = new double[maxnHyp];
double* _dEdx_ST_protcand                 = new double[maxnHyp];
double* _dEdx_CDC_protcand                = new double[maxnHyp];
double* _locDeltaBCAL_prot		  = new double[maxnHyp];
double* _locDeltaTOF_prot		  = new double[maxnHyp];

////photon info
double* _bmE				  = new double[maxnPhotonCandidates];
double* _bmtime				  = new double[maxnPhotonCandidates];

////set branches
void SetBranchAddressesTree(TTree* T){
  T->SetBranchStatus("*",0); 
  
  T->SetBranchAddress("eventNumber",&_eventNumber);
  T->SetBranchAddress("L1TriggerBits",&_L1TriggerBits);
  T->SetBranchAddress("nShower",&_nShower);
  T->SetBranchAddress("nHyp",&_nHyp);
  T->SetBranchAddress("nPhotonCandidates",&_nPhotonCandidates);
  
  T->SetBranchAddress("pX_kminus",_pX_kminus);
  T->SetBranchAddress("pY_kminus",_pY_kminus);
  T->SetBranchAddress("pZ_kminus",_pZ_kminus);
  T->SetBranchAddress("E_kminus",_E_kminus);
  T->SetBranchAddress("t_kminus",_t_kminus);
  T->SetBranchAddress("pX_kplus",_pX_kplus);
  T->SetBranchAddress("pY_kplus",_pY_kplus);
  T->SetBranchAddress("pZ_kplus",_pZ_kplus);
  T->SetBranchAddress("E_kplus",_E_kplus);
  T->SetBranchAddress("t_kplus",_t_kplus);
  T->SetBranchAddress("pX_proton",_pX_proton);
  T->SetBranchAddress("pY_proton",_pY_proton);
  T->SetBranchAddress("pZ_proton",_pZ_proton);
  T->SetBranchAddress("E_proton",_E_proton);
  T->SetBranchAddress("t_proton",_t_proton);
  T->SetBranchAddress("X_vertex",_X_vertex);
  T->SetBranchAddress("Y_vertex",_Y_vertex);
  T->SetBranchAddress("Z_vertex",_Z_vertex);
  T->SetBranchAddress("T_vertex",_T_vertex);
  T->SetBranchAddress("CLKinFit",_CLKinFit);
  T->SetBranchAddress("NDF",_NDF);
  T->SetBranchAddress("ChiSqFit",_ChiSqFit);
  T->SetBranchAddress("Common_Time",_Common_Time);
  
  T->SetBranchAddress("FOM_kminuscand",_FOM_kminuscand );
  T->SetBranchAddress("NDF_kminuscand",_NDF_kminuscand);
  T->SetBranchAddress("ChiSq_kminuscand",_ChiSq_kminuscand);
  T->SetBranchAddress("kminus_dedx_dc_NDF",_kminus_dedx_dc_NDF);
  T->SetBranchAddress("kminus_dedx_dc_ChiSq",_kminus_dedx_dc_ChiSq);
  T->SetBranchAddress("Beta_Timing_kminuscand",_Beta_Timing_kminuscand);
  T->SetBranchAddress("ChiSq_Timing_kminuscand",_ChiSq_Timing_kminuscand);
  T->SetBranchAddress("dEdx_TOF_kminuscand",_dEdx_TOF_kminuscand );
  T->SetBranchAddress("Energy_BCAL_kminuscand",_Energy_BCAL_kminuscand);
  T->SetBranchAddress("Energy_BCALPreshower_kminuscand",_Energy_BCALPreshower_kminuscand);
  T->SetBranchAddress("SigLong_BCAL_kminuscand",_SigLong_BCAL_kminuscand);
  T->SetBranchAddress("SigTheta_BCAL_kminuscand",_SigTheta_BCAL_kminuscand);
  T->SetBranchAddress("RMSTime_BCAL_kminuscand",_RMSTime_BCAL_kminuscand);
  T->SetBranchAddress("E_BCAL_kminuscand",_E_BCAL_kminuscand);
  T->SetBranchAddress("Energy_FCAL_kminuscand",_Energy_FCAL_kminuscand);
  T->SetBranchAddress("E1E9_FCAL_kminuscand",_E1E9_FCAL_kminuscand);
  T->SetBranchAddress("E9E25_FCAL_kminuscand",_E9E25_FCAL_kminuscand);
  T->SetBranchAddress("SumU_FCAL_kminuscand",_SumU_FCAL_kminuscand);
  T->SetBranchAddress("SumV_FCAL_kminuscand",_SumV_FCAL_kminuscand);
  T->SetBranchAddress("TrackBCAL_DeltaPhi_kminuscand",_TrackBCAL_DeltaPhi_kminuscand);
  T->SetBranchAddress("TrackBCAL_DeltaZ_kminuscand",_TrackBCAL_DeltaZ_kminuscand);
  T->SetBranchAddress("dEdx_ST_kminuscand",_dEdx_ST_kminuscand);
  T->SetBranchAddress("dEdx_CDC_kminuscand",_dEdx_CDC_kminuscand);
  T->SetBranchAddress("locDeltaBCAL_kminus",_locDeltaBCAL_kminus);
  T->SetBranchAddress("locDeltaTOF_kminus",_locDeltaTOF_kminus);

  T->SetBranchAddress("FOM_kpluscand",_FOM_kpluscand );
  T->SetBranchAddress("NDF_kpluscand",_NDF_kpluscand);
  T->SetBranchAddress("ChiSq_kpluscand",_ChiSq_kpluscand);
  T->SetBranchAddress("kplus_dedx_dc_NDF",_kplus_dedx_dc_NDF);
  T->SetBranchAddress("kplus_dedx_dc_ChiSq",_kplus_dedx_dc_ChiSq);
  T->SetBranchAddress("Beta_Timing_kpluscand",_Beta_Timing_kpluscand);
  T->SetBranchAddress("ChiSq_Timing_kpluscand",_ChiSq_Timing_kpluscand);
  T->SetBranchAddress("dEdx_TOF_kpluscand",_dEdx_TOF_kpluscand );
  T->SetBranchAddress("Energy_BCAL_kpluscand",_Energy_BCAL_kpluscand);
  T->SetBranchAddress("Energy_BCALPreshower_kpluscand",_Energy_BCALPreshower_kpluscand);
  T->SetBranchAddress("SigLong_BCAL_kpluscand",_SigLong_BCAL_kpluscand);
  T->SetBranchAddress("SigTheta_BCAL_kpluscand",_SigTheta_BCAL_kpluscand);
  T->SetBranchAddress("RMSTime_BCAL_kpluscand",_RMSTime_BCAL_kpluscand);
  T->SetBranchAddress("E_BCAL_kpluscand",_E_BCAL_kpluscand);
  T->SetBranchAddress("Energy_FCAL_kpluscand",_Energy_FCAL_kpluscand);
  T->SetBranchAddress("E1E9_FCAL_kpluscand",_E1E9_FCAL_kpluscand);
  T->SetBranchAddress("E9E25_FCAL_kpluscand",_E9E25_FCAL_kpluscand);
  T->SetBranchAddress("SumU_FCAL_kpluscand",_SumU_FCAL_kpluscand);
  T->SetBranchAddress("SumV_FCAL_kpluscand",_SumV_FCAL_kpluscand);
  T->SetBranchAddress("TrackBCAL_DeltaPhi_kpluscand",_TrackBCAL_DeltaPhi_kpluscand);
  T->SetBranchAddress("TrackBCAL_DeltaZ_kpluscand",_TrackBCAL_DeltaZ_kpluscand);
  T->SetBranchAddress("dEdx_ST_kpluscand",_dEdx_ST_kpluscand);
  T->SetBranchAddress("dEdx_CDC_kpluscand",_dEdx_CDC_kpluscand);
  T->SetBranchAddress("locDeltaBCAL_kplus",_locDeltaBCAL_kplus);
  T->SetBranchAddress("locDeltaTOF_kplus", _locDeltaTOF_kplus);

  T->SetBranchAddress("FOM_protcand",_FOM_protcand );
  T->SetBranchAddress("NDF_protcand",_NDF_protcand);
  T->SetBranchAddress("ChiSq_protcand",_ChiSq_protcand);
  T->SetBranchAddress("prot_dedx_dc_NDF",_prot_dedx_dc_NDF);
  T->SetBranchAddress("prot_dedx_dc_ChiSq",_prot_dedx_dc_ChiSq);
  T->SetBranchAddress("Beta_Timing_protcand",_Beta_Timing_protcand);
  T->SetBranchAddress("ChiSq_Timing_protcan",_ChiSq_Timing_protcand);
  T->SetBranchAddress("dEdx_TOF_protcand",_dEdx_TOF_protcand );
  T->SetBranchAddress("Energy_BCAL_protcand",_Energy_BCAL_protcand);
  T->SetBranchAddress("Energy_BCALPreshower_protcand",_Energy_BCALPreshower_protcand);
  T->SetBranchAddress("SigLong_BCAL_protcand",_SigLong_BCAL_protcand);
  T->SetBranchAddress("SigTheta_BCAL_protcand",_SigTheta_BCAL_protcand);
  T->SetBranchAddress("RMSTime_BCAL_protcand",_RMSTime_BCAL_protcand);
  T->SetBranchAddress("E_BCAL_protcand",_E_BCAL_protcand);
  T->SetBranchAddress("Energy_FCAL_protcand",_Energy_FCAL_protcand);
  T->SetBranchAddress("E1E9_FCAL_protcand",_E1E9_FCAL_protcand);
  T->SetBranchAddress("E9E25_FCAL_protcand",_E9E25_FCAL_protcand);
  T->SetBranchAddress("SumU_FCAL_protcand",_SumU_FCAL_protcand);
  T->SetBranchAddress("SumV_FCAL_protcand",_SumV_FCAL_protcand);
  T->SetBranchAddress("dEdx_ST_protcand",_dEdx_ST_protcand);
  T->SetBranchAddress("dEdx_CDC_protcand",_dEdx_CDC_protcand);
  T->SetBranchAddress("locDeltaBCAL_prot",_locDeltaBCAL_prot);
  T->SetBranchAddress("locDeltaTOF_prot", _locDeltaTOF_prot);

  T->SetBranchAddress("bmE",_bmE);
  T->SetBranchAddress("bmtime",_bmtime);
}

//output tree
////hyperparameters
int Nmax = 100;

////event info
int allphotons;
int shower_number;
int hypothesis_number;
int photon_number;

////kinematic fit
TLorentzVector pkm, pkp, pP, phi, kmp, kpp;
TLorentzVector ppim, ppip, rho;
double pX_km, pY_km, pZ_km, E_km; 
double pX_kp, pY_kp, pZ_kp, E_kp;
double pX_p, pY_p, pZ_p, E_p; 
double theta_km, phi_km, Pmag_km;
double theta_kp, phi_kp, Pmag_kp;
double theta_p, phi_p, Pmag_p;
double theta_phi, phi_phi, mass_phi;
double theta_kmp, phi_kmp, mass_kmp;
double theta_kpp, phi_kpp, mass_kpp;
double s;
double mass_rho;
double X_vtx, Y_vtx, Z_vtx, T_vtx, CL, t0_protcand, commontime;

////KMinus PID
double FOM_kminuscand; 
double NDF_kminuscand; 
double ChiSq_kminuscand; 
double kminus_dedx_dc_NDF; 
double kminus_dedx_dc_ChiSq; 
double Beta_Timing_kminuscand; 
double ChiSq_Timing_kminuscand; 
double dEdx_TOF_kminuscand; 
double Energy_BCAL_kminuscand; 
double Energy_BCALPreshower_kminuscand; 
double SigLong_BCAL_kminuscand; 
double SigTheta_BCAL_kminuscand; 
double RMSTime_BCAL_kminuscand; 
double E_BCAL_kminuscand; 
double Energy_FCAL_kminuscand; 
double E1E9_FCAL_kminuscand; 
double E9E25_FCAL_kminuscand; 
double SumU_FCAL_kminuscand; 
double SumV_FCAL_kminuscand; 
double TrackBCAL_DeltaPhi_kminuscand; 
double TrackBCAL_DeltaZ_kminuscand; 
double dEdx_ST_kminuscand; 
double dEdx_CDC_kminuscand; 
double locDeltaBCAL_kminus; 
double locDeltaTOF_kminus;

////KPlus PID
double FOM_kpluscand; 
double NDF_kpluscand; 
double ChiSq_kpluscand;
double kplus_dedx_dc_NDF;
double kplus_dedx_dc_ChiSq;
double Beta_Timing_kpluscand;
double ChiSq_Timing_kpluscand; 
double dEdx_TOF_kpluscand; 
double Energy_BCAL_kpluscand;
double Energy_BCALPreshower_kpluscand;
double SigLong_BCAL_kpluscand;
double SigTheta_BCAL_kpluscand;
double RMSTime_BCAL_kpluscand;
double E_BCAL_kpluscand;
double Energy_FCAL_kpluscand;
double E1E9_FCAL_kpluscand;
double E9E25_FCAL_kpluscand;
double SumU_FCAL_kpluscand;
double SumV_FCAL_kpluscand;
double TrackBCAL_DeltaPhi_kpluscand;
double TrackBCAL_DeltaZ_kpluscand;
double dEdx_ST_kpluscand;
double dEdx_CDC_kpluscand;
double locDeltaBCAL_kplus;
double locDeltaTOF_kplus;

////proton PID
double FOM_protcand;
double NDF_protcand;
double ChiSq_protcand;
double prot_dedx_dc_NDF;
double prot_dedx_dc_ChiSq;
double Beta_Timing_protcand;
double ChiSq_Timing_protcand;
double dEdx_TOF_protcand;
double Energy_BCAL_protcand;
double Energy_BCALPreshower_protcand;
double SigLong_BCAL_protcand;
double SigTheta_BCAL_protcand;
double RMSTime_BCAL_protcand;
double E_BCAL_protcand;
double Energy_FCAL_protcand;
double E1E9_FCAL_protcand;
double E9E25_FCAL_protcand;
double SumU_FCAL_protcand;
double SumV_FCAL_protcand;
double TrackBCAL_DeltaPhi_protcand;
double TrackBCAL_DeltaZ_protcand;
double dEdx_ST_protcand;
double dEdx_CDC_protcand;
double locDeltaBCAL_prot;
double locDeltaTOF_prot;

////photon info
double* E_photon				= new double[Nmax];
double* massinv2				= new double[Nmax];
double* kmiss					= new double[Nmax];
double* Mmiss					= new double[Nmax];
double* t					= new double[Nmax];
double* u					= new double[Nmax];
double* Emiss					= new double[Nmax];
double* Pmiss					= new double[Nmax];
double* Pmissmin				= new double[Nmax];
double* Pmissplus				= new double[Nmax];
double* Pmissperp				= new double[Nmax];
double* dE					= new double[Nmax];
double* q					= new double[Nmax];
double* omega					= new double[Nmax];
double* dt					= new double[Nmax];
double* xcm					= new double[Nmax];
double* ycm					= new double[Nmax];
double* cosThetaCM 				= new double[Nmax];
double* weight                                  = new double[Nmax];

////set branches
void SetBranchAddressesOutTree(TTree* EvTree){

  EvTree->Branch("shower_number",&shower_number,"shower_number/I");
  EvTree->Branch("hypothesis_number",&hypothesis_number,"hypothesis_number/I");
  EvTree->Branch("photon_number",&photon_number,"photon_number/I");

  EvTree->Branch("pX_km",&pX_km,"pX_km/D");
  EvTree->Branch("pY_km",&pY_km,"pY_km/D");
  EvTree->Branch("pZ_km",&pZ_km,"pZ_km/D");
  EvTree->Branch("E_km",&E_km,"E_km/D");
  EvTree->Branch("pX_kp",&pX_kp,"pX_kp/D");
  EvTree->Branch("pY_kp",&pY_kp,"pY_kp/D");
  EvTree->Branch("pZ_kp",&pZ_kp,"pZ_kp/D");
  EvTree->Branch("E_kp",&E_kp,"E_kp/D");
  EvTree->Branch("pX_p",&pX_p,"pX_p/D");
  EvTree->Branch("pY_p",&pY_p,"pY_p/D");
  EvTree->Branch("pZ_p",&pZ_p,"pZ_p/D");
  EvTree->Branch("E_p",&E_p,"E_p/D");
  EvTree->Branch("theta_km",&theta_km,"theta_km/D");
  EvTree->Branch("phi_km",&phi_km,"phi_km/D");
  EvTree->Branch("Pmag_km",&Pmag_km,"Pmag_km/D");
  EvTree->Branch("theta_kp",&theta_kp,"theta_kp/D");
  EvTree->Branch("phi_kp",&phi_kp,"phi_kp/D");
  EvTree->Branch("Pmag_kp",&Pmag_kp,"Pmag_kp/D");
  EvTree->Branch("theta_p",&theta_p,"theta_p/D");
  EvTree->Branch("phi_p",&phi_p,"phi_p/D");
  EvTree->Branch("Pmag_p",&Pmag_p,"Pmag_p/D");
  EvTree->Branch("theta_phi",&theta_phi,"theta_phi/D");
  EvTree->Branch("phi_phi",&phi_phi,"phi_phi/D");
  EvTree->Branch("mass_phi",&mass_phi,"mass_phi/D");
  EvTree->Branch("theta_kmp",&theta_kmp,"theta_kmp/D");
  EvTree->Branch("phi_kmp",&phi_kmp,"phi_kmp/D");
  EvTree->Branch("mass_kmp",&mass_kmp,"mass_kmp/D");
  EvTree->Branch("theta_kpp",&theta_kpp,"theta_kpp/D");
  EvTree->Branch("phi_kpp",&phi_kpp,"phi_kpp/D");
  EvTree->Branch("mass_kpp",&mass_kpp,"mass_kpp/D");
  EvTree->Branch("s",&s,"s/D");
  EvTree->Branch("mass_rho",&mass_rho,"mass_rho/D");
  EvTree->Branch("X_vtx",&X_vtx,"X_vtx/D");
  EvTree->Branch("Y_vtx",&Y_vtx,"Y_vtx/D");
  EvTree->Branch("Z_vtx",&Z_vtx,"Z_vtx/D");
  EvTree->Branch("T_vtx",&T_vtx,"T_vtx/D");
  EvTree->Branch("CL",&CL,"CL/D");
  EvTree->Branch("t0_protcand",&t0_protcand,"t0_protcand/D");
  EvTree->Branch("commontime",&commontime,"commontime/D");

  EvTree->Branch("FOM_kminuscand",&FOM_kminuscand,"FOM_kminuscand/D");
  EvTree->Branch("NDF_kminuscand",&NDF_kminuscand,"NDF_kminuscand/D");
  EvTree->Branch("ChiSq_kminuscand",&ChiSq_kminuscand,"ChiSq_kminuscand/D");
  EvTree->Branch("kminus_dedx_dc_NDF",&kminus_dedx_dc_NDF,"kminus_dedx_dc_NDF/D");
  EvTree->Branch("kminus_dedx_dc_ChiSq",&kminus_dedx_dc_ChiSq,"kminus_dedx_dc_ChiSq/D");
  EvTree->Branch("Beta_Timing_kminuscand",&Beta_Timing_kminuscand,"Beta_Timing_kminuscand/D");
  EvTree->Branch("ChiSq_Timing_kminuscand",&ChiSq_Timing_kminuscand,"ChiSq_Timing_kminuscand/D");
  EvTree->Branch("dEdx_TOF_kminuscand",&dEdx_TOF_kminuscand,"dEdx_TOF_kminuscand/D");
  EvTree->Branch("Energy_BCAL_kminuscand",&Energy_BCAL_kminuscand,"Energy_BCAL_kminuscand/D");
  EvTree->Branch("Energy_BCALPreshower_kminuscand",&Energy_BCALPreshower_kminuscand,"Energy_BCALPreshower_kminuscand/D");
  EvTree->Branch("SigLong_BCAL_kminuscand",&SigLong_BCAL_kminuscand,"SigLong_BCAL_kminuscand/D");
  EvTree->Branch("SigTheta_BCAL_kminuscand",&SigTheta_BCAL_kminuscand,"SigTheta_BCAL_kminuscand/D");
  EvTree->Branch("RMSTime_BCAL_kminuscand",&RMSTime_BCAL_kminuscand,"RMSTime_BCAL_kminuscand/D");
  EvTree->Branch("E_BCAL_kminuscand",&E_BCAL_kminuscand,"E_BCAL_kminuscand/D");
  EvTree->Branch("Energy_FCAL_kminuscand",&Energy_FCAL_kminuscand,"Energy_FCAL_kminuscand/D");
  EvTree->Branch("E1E9_FCAL_kminuscand",&E1E9_FCAL_kminuscand,"E1E9_FCAL_kminuscand/D");
  EvTree->Branch("E9E25_FCAL_kminuscand",&E9E25_FCAL_kminuscand,"E9E25_FCAL_kminuscand/D");
  EvTree->Branch("SumU_FCAL_kminuscand",&SumU_FCAL_kminuscand,"SumU_FCAL_kminuscand/D");
  EvTree->Branch("SumV_FCAL_kminuscand",&SumV_FCAL_kminuscand,"SumV_FCAL_kminuscand/D");
  EvTree->Branch("TrackBCAL_DeltaPhi_kminuscand",&TrackBCAL_DeltaPhi_kminuscand,"TrackBCAL_DeltaPhi_kminuscand/D");
  EvTree->Branch("TrackBCAL_DeltaZ_kminuscand",&TrackBCAL_DeltaZ_kminuscand,"TrackBCAL_DeltaZ_kminuscand/D");
  EvTree->Branch("dEdx_ST_kminuscand",&dEdx_ST_kminuscand,"dEdx_ST_kminuscand/D");
  EvTree->Branch("dEdx_CDC_kminuscand",&dEdx_CDC_kminuscand,"dEdx_CDC_kminuscand/D");
  EvTree->Branch("locDeltaBCAL_kminus",&locDeltaBCAL_kminus,"locDeltaBCAL_kminus/D");
  EvTree->Branch("locDeltaTOF_kminus",&locDeltaTOF_kminus,"locDeltaTOF_kminus/D");       

  EvTree->Branch("FOM_kpluscand",&FOM_kpluscand,"FOM_kpluscand/D");
  EvTree->Branch("NDF_kpluscand",&NDF_kpluscand,"NDF_kpluscand/D");
  EvTree->Branch("ChiSq_kpluscand",&ChiSq_kpluscand,"ChiSq_kpluscand/D");
  EvTree->Branch("kplus_dedx_dc_NDF",&kplus_dedx_dc_NDF,"kplus_dedx_dc_NDF/D");
  EvTree->Branch("kplus_dedx_dc_ChiSq",&kplus_dedx_dc_ChiSq,"kplus_dedx_dc_ChiSq/D");
  EvTree->Branch("Beta_Timing_kpluscand",&Beta_Timing_kpluscand,"Beta_Timing_kpluscand/D");
  EvTree->Branch("ChiSq_Timing_kpluscand",&ChiSq_Timing_kpluscand,"ChiSq_Timing_kpluscand/D");
  EvTree->Branch("dEdx_TOF_kpluscand",&dEdx_TOF_kpluscand ,"dEdx_TOF_kpluscand/D");
  EvTree->Branch("Energy_BCAL_kpluscand",&Energy_BCAL_kpluscand,"Energy_BCAL_kpluscand/D");
  EvTree->Branch("Energy_BCALPreshower_kpluscand",&Energy_BCALPreshower_kpluscand,"Energy_BCALPreshower_kpluscand/D");
  EvTree->Branch("SigLong_BCAL_kpluscand",&SigLong_BCAL_kpluscand,"SigLong_BCAL_kpluscand/D");
  EvTree->Branch("SigTheta_BCAL_kpluscand",&SigTheta_BCAL_kpluscand,"SigTheta_BCAL_kpluscand/D");
  EvTree->Branch("RMSTime_BCAL_kpluscand",&RMSTime_BCAL_kpluscand,"RMSTime_BCAL_kpluscand/D");
  EvTree->Branch("E_BCAL_kpluscand",&E_BCAL_kpluscand,"E_BCAL_kpluscand/D");
  EvTree->Branch("Energy_FCAL_kpluscand",&Energy_FCAL_kpluscand,"Energy_FCAL_kpluscand/D");
  EvTree->Branch("E1E9_FCAL_kpluscand",&E1E9_FCAL_kpluscand,"E1E9_FCAL_kpluscand/D");
  EvTree->Branch("E9E25_FCAL_kpluscand",&E9E25_FCAL_kpluscand,"E9E25_FCAL_kpluscand/D");
  EvTree->Branch("SumU_FCAL_kpluscand",&SumU_FCAL_kpluscand,"SumU_FCAL_kpluscand/D");
  EvTree->Branch("SumV_FCAL_kpluscand",&SumV_FCAL_kpluscand,"SumU_FCAL_kpluscand/D"); 
  EvTree->Branch("TrackBCAL_DeltaPhi_kpluscand",&TrackBCAL_DeltaPhi_kpluscand,"TrackBCAL_DeltaPhi_kpluscand/D");
  EvTree->Branch("TrackBCAL_DeltaZ_kpluscand",&TrackBCAL_DeltaZ_kpluscand,"TrackBCAL_DeltaZ_kpluscand/D");
  EvTree->Branch("dEdx_ST_kpluscand",&dEdx_ST_kpluscand,"dEdx_ST_kpluscand/D");
  EvTree->Branch("dEdx_CDC_kpluscand",&dEdx_CDC_kpluscand,"dEdx_CDC_kpluscand/D");
  EvTree->Branch("locDeltaBCAL_kplus",&locDeltaBCAL_kplus,"locDeltaBCAL_kplus/D");
  EvTree->Branch("locDeltaTOF_kplus",&locDeltaTOF_kplus,"locDeltaTOF_kplus/D");

  EvTree->Branch("FOM_protcand",&FOM_protcand,"FOM_protcand/D");
  EvTree->Branch("NDF_protcand",&NDF_protcand,"NDF_protcand/D");
  EvTree->Branch("ChiSq_protcand",&ChiSq_protcand,"ChiSq_protcand/D");
  EvTree->Branch("prot_dedx_dc_NDF",&prot_dedx_dc_NDF,"prot_dedx_dc_NDF/D");
  EvTree->Branch("prot_dedx_dc_ChiSq",&prot_dedx_dc_ChiSq,"prot_dedx_dc_ChiSq/D");
  EvTree->Branch("Beta_Timing_protcand",&Beta_Timing_protcand,"Beta_Timing_protcand/D");
  EvTree->Branch("ChiSq_Timing_protcand",&ChiSq_Timing_protcand,"ChiSq_Timing_protcand/D");
  EvTree->Branch("dEdx_TOF_protcand",&dEdx_TOF_protcand,"dEdx_TOF_protcand/D");
  EvTree->Branch("Energy_BCAL_protcand",&Energy_BCAL_protcand,"Energy_BCAL_protcand/D");
  EvTree->Branch("Energy_BCALPreshower_protcand",&Energy_BCALPreshower_protcand,"Energy_BCALPreshower_protcand/D");
  EvTree->Branch("SigLong_BCAL_protcand",&SigLong_BCAL_protcand,"SigLong_BCAL_protcand/D");
  EvTree->Branch("SigTheta_BCAL_protcand",&SigTheta_BCAL_protcand,"SigTheta_BCAL_protcand/D");
  EvTree->Branch("RMSTime_BCAL_protcand",&RMSTime_BCAL_protcand,"RMSTime_BCAL_protcand/D");
  EvTree->Branch("E_BCAL_protcand",&E_BCAL_protcand,"E_BCAL_protcand/D"); 
  EvTree->Branch("Energy_FCAL_protcand",&Energy_FCAL_protcand,"Energy_FCAL_protcand/D");
  EvTree->Branch("E1E9_FCAL_protcand",&E1E9_FCAL_protcand,"E1E9_FCAL_protcand/D");
  EvTree->Branch("E9E25_FCAL_protcand",&E9E25_FCAL_protcand,"E9E25_FCAL_protcand/D");
  EvTree->Branch("SumU_FCAL_protcand",&SumU_FCAL_protcand,"SumU_FCAL_protcand/D");
  EvTree->Branch("SumV_FCAL_protcand",&SumV_FCAL_protcand,"SumV_FCAL_protcand/D");
  EvTree->Branch("TrackBCAL_DeltaPhi_protcand",&TrackBCAL_DeltaPhi_protcand,"TrackBCAL_DeltaPhi_protcand/D");
  EvTree->Branch("TrackBCAL_DeltaZ_protcand",&TrackBCAL_DeltaZ_protcand,"TrackBCAL_DeltaZ_protcand/D");
  EvTree->Branch("dEdx_ST_protcand",&dEdx_ST_protcand,"dEdx_ST_protcand/D");
  EvTree->Branch("dEdx_CDC_protcand",&dEdx_CDC_protcand,"dEdx_CDC_protcand/D");
  EvTree->Branch("locDeltaBCAL_prot",&locDeltaBCAL_prot,"locDeltaBCAL_prot/D");
  EvTree->Branch("locDeltaTOF_prot",&locDeltaTOF_prot,"locDeltaTOF_prot/D");

  EvTree->Branch("allphotons",&allphotons,"allphotons/I");
  EvTree->Branch("E_photon",E_photon,"E_photon[allphotons]/D");
  EvTree->Branch("massinv2",massinv2,"massinv2[allphotons]/D");
  EvTree->Branch("kmiss",kmiss,"kmiss[allphotons]/D");
  EvTree->Branch("Mmiss",Mmiss,"Mmiss[allphotons]/D");
  EvTree->Branch("t",t,"t[allphotons]/D");
  EvTree->Branch("u",u,"u[allphotons]/D");
  EvTree->Branch("Emiss",Emiss,"Emiss[allphotons]/D");
  EvTree->Branch("Pmiss",Pmiss,"Pmiss[allphotons]/D");
  EvTree->Branch("Pmissmin",Pmissmin,"Pmissmin[allphotons]/D");
  EvTree->Branch("Pmissplus",Pmissplus,"Pmissplus[allphotons]/D");
  EvTree->Branch("Pmissperp",Pmissperp,"Pmissperp[allphotons]/D");
  EvTree->Branch("dE",dE,"dE[allphotons]/D");
  EvTree->Branch("q",q,"q[allphotons]/D");
  EvTree->Branch("omega",omega,"omega[allphotons]/D");
  EvTree->Branch("xcm",xcm,"xcm[allphotons]/D");
  EvTree->Branch("ycm",ycm,"ycm[allphotons]/D");
  EvTree->Branch("dt",dt,"dt[allphotons]/D");
  EvTree->Branch("cosThetaCM",cosThetaCM,"cosThetaCM[allphotons]/D");
  EvTree->Branch("weight",weight,"weight[allphotons]/D");  
}
