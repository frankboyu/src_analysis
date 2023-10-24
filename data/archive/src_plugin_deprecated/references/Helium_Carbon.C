#include <TString.h>
void Helium_Carbon(){
       gStyle->SetPalette(1,0);
       gStyle->SetOptStat(1000011);
       gStyle->SetOptFit(1);
       gStyle->SetOptTitle(0);
       gStyle->SetTextFont(132);
       gStyle->SetTextSize(0.05);
       gStyle->SetLabelFont(132,"xyz");
       gStyle->SetTitleSize(0.07,"xyz");
       gStyle->SetLabelSize(0.05,"xyz");
       gStyle->SetTitleOffset(0.7,"Y");
       gStyle->SetTitleOffset(0.7,"X");
      //  gStyle->SetPadLeftMargin(0.12);
       gStyle->SetOptStat(0);
       gStyle->SetMarkerStyle(20);
       gStyle->SetMarkerSize(0.75);
       gStyle->SetHistLineWidth(2);
       gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
       // get rid of X error bars and y error bar caps
        gStyle->SetErrorX(0.001);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gROOT->ForceStyle();

	   
         TH1F *hmassrho_he = new TH1F("hmassrho_he","Inv_mass_He", 100, 0, 2);          ///  Mass of Rho of helium   
         TH1F *hmassrhoacc_he = new TH1F("hmassrhoacc_he","Inv_mass_acc_He", 100, 0, 2);/// Accidental Subtraction for Rho Mass        
         TH1F *htrho_he = new TH1F("htrho_he","|t| He", 20, 0.1, 15);
         htrho_he->Sumw2();        
         TH1F *htrhoacc_he = new TH1F("htrhoacc_he","|t| acc", 20, 0.1, 15);
         htrhoacc_he->Sumw2();
         TH1F *htrhom_he = new TH1F("htrhom_he","|t| He", 20, 0.1, 15);
         htrhom_he->Sumw2();        
         TH1F *htrhomacc_he = new TH1F("htrhomacc_he","|t| acc", 20, 0.1, 15);
         htrhomacc_he->Sumw2();
       
///////////////////////Helium///////////////////////////////////////////////////
  TChain *tt_he = new TChain("EvTree");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090076/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090077/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090078/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090080/rhocandidates.root");
    //   tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090127/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090152/rhocandidates.root");
   tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090153/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090154/rhocandidates.root");
    tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090155/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090158/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090161/rhocandidates.root");
  //  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090162/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090163/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090164/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090165/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090167/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090168/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090169/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090170/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090171/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090174/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090175/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090176/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090177/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090178/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090179/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090181/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090182/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090183/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090184/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090185/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090190/rhocandidates.root");
  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090191/rhocandidates.root");
  //  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090192/rhocandidates.root");
  //  tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090193/rhocandidates.root");
  // tt_he->Add("/volatile/halld/home/src-ct/analysis/1p2pi/ver02/090191/rhocandidates.root");

  Double_t dEdx_CDC_protcand;
  Double_t Pmag_p;
  Double_t locDeltaBCAL_prot;
  Double_t locDeltaTOF_prot;
  Double_t dEdx_CDC_pipluscand;
  Double_t Pmag_pip;
  Double_t locDeltaBCAL_piplus;
  Double_t locDeltaTOF_piplus;
  Double_t Z_vtx;
  Double_t mass_rho;

  //===== From here down are arrays =========//
  // The size of the Intime Variables is given by InTimePhotons //
  // The size of the OUtime (or bgs) is given by OutTimePhotons//
  
 
  Int_t InTimePhotons, OutTimePhotons;
  Int_t maxarr =20; 
  Double_t Pmissmin[maxarr];
  Double_t trho[maxarr];
  Double_t u[maxarr];
  Double_t tbgd[maxarr]; 
  Double_t ubgd[maxarr];
  Double_t Pmissminbgd[maxarr];
  Double_t massinv2bgd[maxarr];         
  
  tt_he->SetBranchAddress("dEdx_CDC_protcand" , &dEdx_CDC_protcand);
  tt_he->SetBranchAddress("Pmag_p", &Pmag_p);
  tt_he->SetBranchAddress("locDeltaBCAL_prot", &locDeltaBCAL_prot);
  tt_he->SetBranchAddress("locDeltaTOF_prot", &locDeltaTOF_prot);
  tt_he->SetBranchAddress("dEdx_CDC_pipluscand" , &dEdx_CDC_pipluscand);
  tt_he->SetBranchAddress("Pmag_pip" , &Pmag_pip);
  tt_he->SetBranchAddress("locDeltaBCAL_piplus" , &locDeltaBCAL_piplus);
  tt_he->SetBranchAddress("Z_vtx", &Z_vtx );
  tt_he->SetBranchAddress("locDeltaTOF_piplus", &locDeltaTOF_piplus );
  tt_he->SetBranchAddress("InTimePhotons", &InTimePhotons );
  tt_he->SetBranchAddress("OutTimePhotons", &OutTimePhotons );
  tt_he->SetBranchAddress("mass_rho" , &mass_rho);  
  tt_he->SetBranchAddress("Pmissmin",  Pmissmin);
  tt_he->SetBranchAddress("t",trho);
  tt_he->SetBranchAddress("u" ,u);
  tt_he->SetBranchAddress("tbgd",  tbgd );
  tt_he->SetBranchAddress("ubgd", ubgd);
  tt_he->SetBranchAddress("Pmissminbgd", Pmissminbgd);
  tt_he->SetBranchAddress("massinv2bgd",massinv2bgd);

   
  Bool_t targetcut;
  Bool_t datacut;
  Bool_t datacut2;
  Bool_t tcut1;
  Bool_t datacut2acc;
  Bool_t tacccut1; 
  Bool_t Mrhocut; 
  Bool_t Mrhocutacc; 

  Long64_t nentriesHe = tt_he->GetEntries();
  for (int kk=0; kk<nentriesHe;  kk++){
       tt_he->GetEntry(kk);

	tcut1     = false;
        tacccut1  = false;        
	targetcut = false;
	datacut   = false;
	datacut2  = false;
	datacut2acc = false; 
        Mrhocut = mass_rho > 0.62 && mass_rho < 0.92;
	
        // Target Cut
       targetcut =  Z_vtx > 51 && Z_vtx < 79;
         // Data Cut 
   datacut =((dEdx_CDC_protcand!=0 && dEdx_CDC_protcand*1e6>exp(-4*Pmag_p + 2.25 ) + 1 && fabs(locDeltaBCAL_prot)<1 ) || (dEdx_CDC_protcand ==0    && fabs(locDeltaTOF_prot)<0.6 ) )  && Pmag_p < 10 && ((dEdx_CDC_pipluscand*1e6 < exp(-7*Pmag_pip+3)+6.2 && dEdx_CDC_pipluscand!=0 && fabs(locDeltaBCAL_piplus)<1 ) || (dEdx_CDC_pipluscand==0 && fabs(locDeltaTOF_piplus)<0.5 ));


          if (!datacut || !targetcut) continue;

//	cout << "============== InTimePhotons " << InTimePhotons << endl; 

	for(int j=0; j<InTimePhotons; j++){
	   if(trho[j]>1 && u[j] >1) {
	    tcut1=true ;
	    //	   if ((u[j]>1.5*Pmissmin[j] +1) && (trho[j] > -0.2*u[j] +2)) {
	    //  datacut2=true;
	for (int ii=0; ii<13; ii++){
	      htrho_he->Fill(trho[j]);
              hmassrho_he->Fill(mass_rho);
	   
           if(Mrhocut)htrhom_he->Fill(trho[j]);
	}
	      //}
	   }
	}
 
       	for(int k=0; k<OutTimePhotons; k++){
          if(tbgd[k]>1 && ubgd[k] > 1) {
	   tacccut1=true;
	   //	  if ((ubgd[k]>1.5*Pmissminbgd[k] +1) && (tbgd[k] > -0.2*ubgd[k] +2)) {
	   //  datacut2acc=true;
	for (int ii=0; ii<13; ii++){
             htrhoacc_he->Fill(tbgd[k]);
             hmassrhoacc_he->Fill(mass_rho);
	   if(Mrhocut)htrhomacc_he->Fill(tbgd[k]);
	  	}
	   // }
	  }       
      }
  }//entries
     htrho_he->Add(htrhoacc_he,-1/6.);
     htrhom_he->Add(htrhomacc_he,-1/6.);
     hmassrho_he->Add(hmassrhoacc_he,-1/6.);
   

////////////////////////////////////////////////Deuterium////////////////////////////////////////////

        TH1F *hmassrho_D = new TH1F("hmassrho_C","Inv_mass_C", 100, 0, 2);          ///  Mass of Rho of helium   
        TH1F *hmassrhoacc_D = new TH1F("hmassrhoacc_C","Inv_mass_acc_C", 100, 0, 2);/// Accidental Subtraction for Rho Mass       
        TH1F *htrho_D = new TH1F("htrho_C","|t| C", 20, 0.1, 15);
        htrho_D->Sumw2();      
        TH1F *htrhoacc_D = new TH1F("htrhoacc_C","|t| acc", 20, 0.1, 15);
        htrhoacc_D->Sumw2();      
         TH1F *htrhom_D = new TH1F("htrhom_C","|t| C", 20, 0.1, 15);
         htrhom_D->Sumw2();        
         TH1F *htrhomacc_D = new TH1F("htrhomacc_C","|t| acc", 20, 0.1, 15);
         htrhomacc_D->Sumw2();


  TChain *tt_D = new TChain("EvTree");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090262/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090263/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090315/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090267/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090270/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090281/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090282/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090283/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090284/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090285/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090286/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090287/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090314/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090291/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090292/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090293/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090294/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090295/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090296/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090297/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090298/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090299/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090301/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090313/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090305/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090306/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090307/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090312/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090309/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090310/rhocandidates.root");
    tt_D->Add("/volatile/halld/home/src-ct/analysis/C/090311/rhocandidates.root");


  Double_t dEdx_CDC_protcand_D;
  Double_t Pmag_p_D;
  Double_t locDeltaBCAL_prot_D;
  Double_t locDeltaTOF_prot_D;
  Double_t dEdx_CDC_pipluscand_D;
  Double_t Pmag_pip_D;
  Double_t locDeltaBCAL_piplus_D;
  Double_t locDeltaTOF_piplus_D;
  Double_t Z_vtx_D;
  Double_t mass_rho_D;

  //===== From here down are arrays =========//
  // The size of the Intime Variables is given by InTimePhotons //
  // The size of the OUtime (or bgs) is given by OutTimePhotons//
  
 
  Int_t InTimePhotons_D, OutTimePhotons_D;
  Double_t Pmissmin_D[maxarr];
  Double_t trho_D[maxarr];
  Double_t u_D[maxarr];
  Double_t tbgd_D[maxarr]; 
  Double_t ubgd_D[maxarr];
  Double_t Pmissminbgd_D[maxarr];
  Double_t massinv2bgd_D[maxarr];  
  Double_t massinv2_D;     
  
  tt_D->SetBranchAddress("dEdx_CDC_protcand" , &dEdx_CDC_protcand_D);
  tt_D->SetBranchAddress("Pmag_p", &Pmag_p_D);
  tt_D->SetBranchAddress("locDeltaBCAL_prot", &locDeltaBCAL_prot_D);
  tt_D->SetBranchAddress("locDeltaTOF_prot", &locDeltaTOF_prot_D);
  tt_D->SetBranchAddress("dEdx_CDC_pipluscand" , &dEdx_CDC_pipluscand_D);
  tt_D->SetBranchAddress("Pmag_pip" , &Pmag_pip_D);
  tt_D->SetBranchAddress("locDeltaBCAL_piplus" , &locDeltaBCAL_piplus_D);
  tt_D->SetBranchAddress("Z_vtx", &Z_vtx_D );
  tt_D->SetBranchAddress("locDeltaTOF_piplus", &locDeltaTOF_piplus_D );  
  tt_D->SetBranchAddress("InTimePhotons", &InTimePhotons_D );
  tt_D->SetBranchAddress("OutTimePhotons", &OutTimePhotons_D );
  tt_D->SetBranchAddress("mass_rho" , &mass_rho_D);
  tt_D->SetBranchAddress("Pmissmin",  Pmissmin_D);
  tt_D->SetBranchAddress("t",trho_D);
  tt_D->SetBranchAddress("u" ,u_D);
  tt_D->SetBranchAddress("tbgd",  tbgd_D );
  tt_D->SetBranchAddress("ubgd", ubgd_D);
  tt_D->SetBranchAddress("Pmissminbgd", Pmissminbgd_D);
  tt_D->SetBranchAddress("massinv2bgd",massinv2bgd_D);

   
  Bool_t targetcut_D;
  Bool_t datacut_D;
  Bool_t datacut2_D;
  Bool_t tcut1_D;
  Bool_t datacut2acc_D;
  Bool_t tacccut1_D;
  Bool_t Mrhocut_D;
  Bool_t Mrhocutacc_D;
  
  Long64_t nentriesD = tt_D->GetEntries();
  for (int ll=0; ll<nentriesD;  ll++){
       tt_D->GetEntry(ll);

	   targetcut_D = false;
           tcut1_D     = false;       
           tacccut1_D  = false;
	   datacut_D   = false;
	   datacut2_D  = false;
	   datacut2acc_D = false; 
	   Mrhocut_D = mass_rho_D > 0.62 && mass_rho_D < 0.92;
	
        // Target Cut
       targetcut_D =  Z_vtx_D > 51 && Z_vtx_D < 79;
         // Data Cut 
   datacut_D =((dEdx_CDC_protcand_D!=0 && dEdx_CDC_protcand_D*1e6>exp(-4*Pmag_p_D + 2.25 ) + 1 && fabs(locDeltaBCAL_prot_D)<1 ) || (dEdx_CDC_protcand_D ==0    && fabs(locDeltaTOF_prot_D)<0.6 ) )  && Pmag_p_D < 10 && ((dEdx_CDC_pipluscand_D*1e6 < exp(-7*Pmag_pip_D+3)+6.2 && dEdx_CDC_pipluscand_D!=0 && fabs(locDeltaBCAL_piplus_D)<1 ) || (dEdx_CDC_pipluscand_D==0 && fabs(locDeltaTOF_piplus_D)<0.5 ));


          if (!datacut_D || !targetcut_D) continue;	 
//	cout << "============== InTimePhotons " << InTimePhotons << endl; 
	for(int j=0; j<InTimePhotons_D; j++){
	   if(trho_D[j]>1 && u_D[j] >1) {
	    tcut1_D=true ;
	    //	    if ((u_D[j]>1.5*Pmissmin_D[j] +1) && (trho_D[j] > -0.2*u_D[j] +2)) {
	    //	      datacut2_D=true;
	for (int ii=0; ii<25; ii++){
	      htrho_D->Fill(trho_D[j]);
              hmassrho_D->Fill(mass_rho_D);
	      if(Mrhocut_D)htrhom_D->Fill(trho_D[j]);
	}
	      // }
	   } 
	}
        
	for(int k=0; k<OutTimePhotons_D; k++){
          if(tbgd_D[k]>1 && ubgd_D[k] > 1) {
	   tacccut1_D=true;
	   //	   if ((ubgd_D[k]>1.5*Pmissminbgd_D[k] +1) && (tbgd_D[k] > -0.2*ubgd_D[k] +2)) {
	   // datacut2acc_D=true;
	for (int ii=0; ii<25; ii++){
             htrhoacc_D->Fill(tbgd_D[k]);
             hmassrhoacc_D->Fill(mass_rho_D);
	     if(Mrhocut_D)htrhomacc_D->Fill(tbgd_D[k]);
	}
	     // }
          }
        }      
   }
  
   htrho_D->Add(htrhoacc_D,-1/6.);
   htrhom_D->Add(htrhomacc_D,-1/6.);
   hmassrho_D->Add(hmassrhoacc_D,-1/6.);
  
 

/////////////////////////////////////////////////////Ratio ////////////////////////////////////
//  TH1F *hRatio_HeD = new TH1F (*htrho_he);
  TH1F *hRatio_HeD = new TH1F("hRatio_HeC","|t|", 20, 0.1, 15);
  hRatio_HeD -> Divide(htrho_he, htrho_D,1.0,1.35);	  	
  //  TH1F *hRatiom_HeD = new TH1F (*htrhom_he);
  TH1F *hRatiom_HeD = new TH1F("hRatiom_HeC","|t|", 20, 0.1, 15);
  hRatiom_HeD -> Divide(htrhom_D, htrhom_he, 1.0/27963, 1.0/12902);	  	


	TCanvas *c1 = new TCanvas("c1","c1", 1000,800);
        c1->Divide(2,3);
        c1->cd(1);
	
	 hmassrho_he->Draw("HIST");
	 hmassrho_he->GetXaxis()->SetTitle("Invariant Mass (GeV)");
	 hmassrho_he->GetYaxis()->SetTitle("He Counts");
        c1->Update();
        c1->cd(2);
        hmassrho_D->GetXaxis()->SetTitle("Invariant Mass (GeV)");
        hmassrho_D->GetYaxis()->SetTitle("C Counts");
        hmassrho_D->Draw("HIST");
        c1->Update();
        c1->cd(3);
	gPad->SetLogy();
        htrho_he-> SetMarkerStyle(kFullCircle);
        htrho_he->SetMarkerSize(0.8);
        htrho_he->GetXaxis()->SetTitle("|t| (GeV^{1})");
        htrho_he->GetYaxis()->SetTitle("He Counts");
        htrho_he->Draw();
        c1->Update();
        c1->cd(4);
	gPad->SetLogy();
        htrho_D-> SetMarkerStyle(kFullCircle);
        htrho_D->SetMarkerSize(0.8);
        htrho_D->GetXaxis()->SetTitle("|t| (GeV^{1})");
        htrho_D->GetYaxis()->SetTitle("C Counts");
        htrho_D->Draw();	
        c1->Update();
        c1->cd(5);
	gPad->SetLogy();
        htrhom_he-> SetMarkerStyle(kFullCircle);
        htrhom_he->SetMarkerSize(0.8);
        htrhom_he->GetXaxis()->SetTitle("|t| (GeV^{1})");
        htrhom_he->GetYaxis()->SetTitle("H Counts (with mass cut)");
        htrhom_he->Draw();
        c1->Update();
        c1->cd(6);
	gPad->SetLogy();
        htrhom_D-> SetMarkerStyle(kFullCircle);
        htrhom_D->SetMarkerSize(0.8);
        htrhom_D->GetXaxis()->SetTitle("|t| (GeV^{1})");
        htrhom_D->GetYaxis()->SetTitle("C Counts (with inv. mass cut)");
        htrhom_D->Draw();	
        c1->Update();

        c1->SaveAs("HeD_mass and t plots.pdf");

      //  gStyle->SetOptStat(1000011);
	

        gStyle->SetOptStat(0);
	TCanvas *cHeD = new TCanvas("cHeC","cHeC", 1000,800);
	//         cHeD->Divide(1,2);
	// cHeD->cd(1);
	// hRatio_HeD-> SetMarkerStyle(kFullCircle);
	//hRatio_HeD ->SetMarkerSize(1.0);
	//hRatio_HeD->GetXaxis()->SetTitle("|t| (GeV^{1})");
	//hRatio_HeD->GetYaxis()->SetTitle("Ratio He/D (Full)");
	//hRatio_HeD->GetYaxis()->SetRangeUser(0.,1.0); 
	//hRatio_HeD->Draw("P E1");
      
        //cHeD->Update();

	// cHeD->cd(2);
       hRatiom_HeD-> SetMarkerStyle(kFullCircle);
       hRatiom_HeD ->SetMarkerSize(0.8);
       hRatiom_HeD->GetXaxis()->SetTitle("|t| (GeV^{1})");
       hRatiom_HeD->GetYaxis()->SetTitle("Ratio C/He");
       hRatiom_HeD->GetYaxis()->SetRangeUser(0.,1.5);
       hRatiom_HeD->Draw("P E1");
      
  cHeD->Update();

  cHeD->SaveAs("HeD_ratioplots.pdf");
  }//void


 
  
  
