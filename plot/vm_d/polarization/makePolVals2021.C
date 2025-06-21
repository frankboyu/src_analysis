#include <fstream>
#include <string>
#include <vector>
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TH1.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"




//SimFit
int nParPhi0 = 5;
int iparPhi0[5] = { 0,
                    1,
                    2,//Shared (First phi offset)
                    3,//Shared (Amplitude of cosPhi)
                    4 //Shared (Second phi offset)
};
int nParPhi45 = 5;
int iparPhi45[5] = { 5,
                     6,
                     2,//Shared
                     3,//Shared
                     4 //Shared
};
int nParPhi90 = 5;
int iparPhi90[5] = { 7,
                     8,                                   
                     2,//Shared
                     3,//Shared
                     4 //Shared
};

int nParPhi135 = 5;
int iparPhi135[5] = { 9,
                      10,                                
		      2,//Shared
                      3,//Shared
                      4 //Shared
};
const int Npar = 11;

double gAmp1 = 160000;
double gAsymm      = 0.07;
double gPhiOffset1 = 6;
double gF1         = 0.1;
double gPhiOffset2 = -89;
double gAmp2 = 160000;
double gAsymm2 = 0.06;

double gAmp3 = 160000;
double gAsymm3 = 0.06;

double gAmp4 = 160000;
double gAsymm4 = 0.06;

double simPars[Npar] = {gAmp1,gAsymm,gPhiOffset1,gF1,gPhiOffset2,gAmp2,gAsymm2,gAmp3,gAsymm3,gAmp4,gAsymm4};

// Create the GlobalCHi2 structure                                                                                   

struct GlobalChi2 {
  GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
               ROOT::Math::IMultiGenFunction & f2,
               ROOT::Math::IMultiGenFunction & f3,
               ROOT::Math::IMultiGenFunction & f4) :
    fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3), fChi2_4(&f4) {}

  double operator() (const double *par) const {
    double p1[nParPhi0];
    for (int i = 0; i < nParPhi0; ++i) p1[i] = par[iparPhi0[i] ];

    double p2[nParPhi45];
    for (int i = 0; i < nParPhi45; ++i) p2[i] = par[iparPhi45[i] ];

    double p3[nParPhi90];
    for (int i = 0; i < nParPhi90; ++i) p3[i] = par[iparPhi90[i] ];

    double p4[nParPhi135];
    for (int i = 0; i < nParPhi135; ++i) p4[i] = par[iparPhi135[i] ];

    return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3) + (*fChi2_4)(p4);
  }

  const  ROOT::Math::IMultiGenFunction * fChi2_1;
  const  ROOT::Math::IMultiGenFunction * fChi2_2;
  const  ROOT::Math::IMultiGenFunction * fChi2_3;
  const  ROOT::Math::IMultiGenFunction * fChi2_4;
};



//Structures
struct tpolToGo_t {
  Double_t polVal0;
  Double_t polErr0;
  Double_t polVal45;
  Double_t polErr45;
  Double_t polVal90;
  Double_t polErr90;
  Double_t polVal135;
  Double_t polErr135;
};
struct fluxToGo_t {
  TH1D *hNorm0;
  TH1D *hNorm45;
  TH1D *hNorm90;
  TH1D *hNorm135;
};

struct phiHistos_t {
  TH1D *hPhi0;
  TH1D *hPhi45;
  TH1D *hPhi90;
  TH1D *hPhi135;
};
struct phiPars_t {
  double pol0;
  double err0;
  double pol45;
  double err45;
  double pol90;
  double err90;
  double pol135;
  double err135;
  double offset;
  double aOffset;
  double aCenter;

};

//Function prototypes
phiPars_t phiSimFit(phiHistos_t phiHistos);
int readVec(char *fToGo,vector<int> *vec);
int writeRuns(char *fTogo,TH1D *hMask);
int setMask2d(TH2D *hMask, vector<int> *vec);
int setMask1d(TH1D *hMask, vector<int> *vec);
void setHistoParsE(TH1D *histo);
void setHistoParsEAmo(TH1D *histo);
void phi1D(TH2D *hPhiVsRun,TH2D *hMask,TH1D *hPhi,char *hString, TF1 *fitFun,bool qVal);
tpolToGo_t makePolValsTmp2021(int year,int part,int userSet,int radTh,bool mixMode);

void makePolVals2021(int userSet){

  bool mixMode = false;
  tpolToGo_t tpolToGo75;
  tpolToGo_t tpolToGo750;

  int part = 1;
  int radTh = 75;
  int year = 2021;

  
  makePolValsTmp2021(year,part,userSet,radTh,mixMode);


  return;
}


tpolToGo_t makePolValsTmp2021(int year,int part,int userSet,int radTh,bool mixMode){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(5);
  char *s = new char[1];
  int fitShow = -9000;
  //fitShow = 0;
  tpolToGo_t tpolToGo;
  tpolToGo.polVal0   = 0;
  tpolToGo.polErr0   = 0;
  tpolToGo.polVal45  = 0;
  tpolToGo.polErr45  = 0;
  tpolToGo.polVal90  = 0;
  tpolToGo.polErr90  = 0;
  tpolToGo.polVal135 = 0;
  tpolToGo.polErr135 = 0;
  //bool userSet = false;
  bool sp16Set = false;
  bool amoPlot = false;
  bool sp17Low = false;
  bool sp17High = false;
  bool sp18Set  = false;
  bool fa18Set  = false;
  bool sp20Set  = false;
  bool y21Set  = true;
  bool rad750   = false;

  double anaPower = 0.1904;

  if (userSet != 0 && userSet != 1){
      cout<<"The user defined set definition must be 0 or 1"<<endl;
      cout<<"User defined set definition: 0 -> Standard run set"<<endl;
      cout<<"                             1 -> User defined run set"<<endl;
      cout<<"You have the user defined set definition set to "<<userSet<<endl;
      return tpolToGo;
  }

  char  fString[120];
  vector<int> veco0;
  vector<int> veco45;
  vector<int> veco90;
  vector<int> veco135;
  vector<int> vecAmo;
  vector<int> vecAmoAll;
  vector<int> vecAmoGood;
  vector<int> vecAmoBad1;
  vector<int> vecAmoBad2;
  vector<int> vecAmoBad3;
  vector<int> vecUser;

  //READ AMO
  //sprintf(fString,"./inFiles/amoSet.txt");
  //readVec(fString,&vecAmo);

  //READ 0
  sprintf(fString,"./inFiles/runList21-0-75Full.txt");
  readVec(fString,&veco0);

  //READ 45

  sprintf(fString,"./inFiles/runList21-45-75Full.txt");
  readVec(fString,&veco45);

  //READ 90
  sprintf(fString,"./inFiles/runList21-90-75Full.txt");
  readVec(fString,&veco90);

  //READ 135
  sprintf(fString,"./inFiles/runList21-135-75Full.txt");
  readVec(fString,&veco135);

  sprintf(fString,"./inFiles/userSet.txt");


  readVec(fString,&vecUser);
  
  sprintf(fString,"./inFiles/all2021.root");

  TFile *hFile1 = new TFile(fString);
  
  TH2D *hA0RunTag; hFile1->GetObject("hA0RunTag",hA0RunTag);
  int nBinsX  = hA0RunTag->GetNbinsX();
  double xMin = hA0RunTag->GetXaxis()->GetXmin();
  double xMax = hA0RunTag->GetXaxis()->GetXmax();
  TH1D *hMaskAmo  = new TH1D("hMaskAmo","",nBinsX,xMin,xMax);
  TH1D *hMask0    = new TH1D("hMask0","",nBinsX,xMin,xMax);
  TH1D *hMask45   = new TH1D("hMask45","",nBinsX,xMin,xMax);
  TH1D *hMask90   = new TH1D("hMask90","",nBinsX,xMin,xMax);
  TH1D *hMask135  = new TH1D("hMask135","",nBinsX,xMin,xMax);
  TH1D *hMaskUser = new TH1D("hMaskUser","",nBinsX,xMin,xMax);
  TH1D *hMaskTmp  = new TH1D("hMaskTmp","",nBinsX,xMin,xMax);
  TH1D *hMaskAll  = new TH1D("hMaskAll","",nBinsX,xMin,xMax);
  TH1D *hMaskSub  = new TH1D("hMaskSub","",nBinsX,xMin,xMax);
  
  //setMask1d(hMaskAmo,&vecAmo);
  setMask1d(hMask0,&veco0);
  setMask1d(hMask45,&veco45);
  setMask1d(hMask90,&veco90);
  setMask1d(hMask135,&veco135);
  setMask1d(hMaskUser,&vecUser);

  
  //hMaskAll->Add(hMaskAmo,hMask0);
  hMaskAll->Add(hMaskAll,hMask45);
  hMaskAll->Add(hMaskAll,hMask90);
  hMaskAll->Add(hMaskAll,hMask135);
  
  hMaskSub->Add(hMaskAll,hMaskUser,1.0,-1.0);


  bool badRun = false;
  if (userSet == true && mixMode == false){ 
    for (int run = 0; run < vecUser.size(); run++) {
      if (vecUser.at(run) > xMax || vecUser.at(run) < xMin){
	cout<<"!!!!! Run number "<<vecUser.at(run)<<" is outside range of TPOL master list for this run period!!!!!"<<endl;
	badRun = true;
      }
    }


    for (int iTest = 1; iTest <= nBinsX; iTest++) {
      double testVal = hMaskSub->GetBinContent(iTest);
      if (testVal < 0){
	int runNumber = (int)hMaskSub->GetBinCenter(iTest);
	cout<<"!!!!! Run number "<<runNumber<<" is NOT in TPOL master list !!!!!"<<endl;
	badRun = true;
      }
    }
  }

  if (badRun == true){
    cout<<""<<endl;
    cout<<"Remove the above run number(s) from ./inFiles/userSet.txt and try again"<<endl;
    cout<<""<<endl;
    return tpolToGo;
  }

  if (userSet == true) {
    //hMaskAmo->Multiply(hMaskAmo,hMaskUser);
    hMask0->Multiply(hMask0,hMaskUser);
    hMask45->Multiply(hMask45,hMaskUser);
    hMask90->Multiply(hMask90,hMaskUser);
    hMask135->Multiply(hMask135,hMaskUser);
  }

  if (sp18Set == true || fa18Set == true || sp20Set == true || y21Set == true) {

    TH2D *hPhiCutVsRun; hFile1->GetObject("hPhiCutVsRun",hPhiCutVsRun);

    TH2D *hPhiCutVsRunAtEg7100; hFile1->GetObject("hPhiCutVsRunAtEg7",hPhiCutVsRunAtEg7100);
    TH2D *hPhiCutVsRunAtEg7300; hFile1->GetObject("hPhiCutVsRunAtEg8",hPhiCutVsRunAtEg7300);
    TH2D *hPhiCutVsRunAtEg7500; hFile1->GetObject("hPhiCutVsRunAtEg9",hPhiCutVsRunAtEg7500);
    TH2D *hPhiCutVsRunAtEg7700; hFile1->GetObject("hPhiCutVsRunAtEg10",hPhiCutVsRunAtEg7700);
    TH2D *hPhiCutVsRunAtEg7900; hFile1->GetObject("hPhiCutVsRunAtEg11",hPhiCutVsRunAtEg7900);

    TH2D *hPhiCutVsRunAtEg8100; hFile1->GetObject("hPhiCutVsRunAtEg12",hPhiCutVsRunAtEg8100);
    TH2D *hPhiCutVsRunAtEg8300; hFile1->GetObject("hPhiCutVsRunAtEg13",hPhiCutVsRunAtEg8300);
    TH2D *hPhiCutVsRunAtEg8500; hFile1->GetObject("hPhiCutVsRunAtEg14",hPhiCutVsRunAtEg8500);
    TH2D *hPhiCutVsRunAtEg8700; hFile1->GetObject("hPhiCutVsRunAtEg15",hPhiCutVsRunAtEg8700);
    TH2D *hPhiCutVsRunAtEg8900; hFile1->GetObject("hPhiCutVsRunAtEg16",hPhiCutVsRunAtEg8900);

    TH2D *hPhiCutVsRunAtEg9100; hFile1->GetObject("hPhiCutVsRunAtEg17",hPhiCutVsRunAtEg9100);
    TH2D *hPhiCutVsRunAtEg9300; hFile1->GetObject("hPhiCutVsRunAtEg18",hPhiCutVsRunAtEg9300);
    TH2D *hPhiCutVsRunAtEg9500; hFile1->GetObject("hPhiCutVsRunAtEg19",hPhiCutVsRunAtEg9500);
    TH2D *hPhiCutVsRunAtEg9700; hFile1->GetObject("hPhiCutVsRunAtEg20",hPhiCutVsRunAtEg9700);
    TH2D *hPhiCutVsRunAtEg9900; hFile1->GetObject("hPhiCutVsRunAtEg21",hPhiCutVsRunAtEg9900);

    TH2D *hPhiCutVsRunAtEg10100; hFile1->GetObject("hPhiCutVsRunAtEg22",hPhiCutVsRunAtEg10100);
    TH2D *hPhiCutVsRunAtEg10300; hFile1->GetObject("hPhiCutVsRunAtEg23",hPhiCutVsRunAtEg10300);
    TH2D *hPhiCutVsRunAtEg10500; hFile1->GetObject("hPhiCutVsRunAtEg24",hPhiCutVsRunAtEg10500);
    TH2D *hPhiCutVsRunAtEg10700; hFile1->GetObject("hPhiCutVsRunAtEg25",hPhiCutVsRunAtEg10700);
    TH2D *hPhiCutVsRunAtEg10900; hFile1->GetObject("hPhiCutVsRunAtEg26",hPhiCutVsRunAtEg10900);

    TH2D *hPhiCutVsRunAtEg11100; hFile1->GetObject("hPhiCutVsRunAtEg27",hPhiCutVsRunAtEg11100);
    TH2D *hPhiCutVsRunAtEg11300; hFile1->GetObject("hPhiCutVsRunAtEg28",hPhiCutVsRunAtEg11300);
    TH2D *hPhiCutVsRunAtEg11500; hFile1->GetObject("hPhiCutVsRunAtEg29",hPhiCutVsRunAtEg11500);
    TH2D *hPhiCutVsRunAtEg11700; hFile1->GetObject("hPhiCutVsRunAtEg30",hPhiCutVsRunAtEg11700);



    hPhiCutVsRun->Add(hPhiCutVsRunAtEg8300,hPhiCutVsRunAtEg8500);
    hPhiCutVsRun->Add(hPhiCutVsRun,hPhiCutVsRunAtEg8700);
    if (sp20Set == true || y21Set == true) {
      hPhiCutVsRun->Add(hPhiCutVsRunAtEg8100,hPhiCutVsRunAtEg8300);
      hPhiCutVsRun->Add(hPhiCutVsRun,hPhiCutVsRunAtEg8500);
      
    }
    //TH2D *hMask2dEAmo  = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dEAmo");
    TH2D *hMask2dE0    = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dE0");
    TH2D *hMask2dE45   = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dE45");
    TH2D *hMask2dE90   = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dE90");
    TH2D *hMask2dE135  = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dE135");
    TH2D *hMask2dEUser = (TH2D*)hPhiCutVsRunAtEg8500->Clone("hMask2dEUser");
    
    //setMask2d(hMask2dEAmo,&vecAmo);
    setMask2d(hMask2dE0,&veco0);
    setMask2d(hMask2dE45,&veco45);
    setMask2d(hMask2dE90,&veco90);
    setMask2d(hMask2dE135,&veco135);
    setMask2d(hMask2dEUser,&vecUser);
    
    if (userSet == true) {
      //hMask2dEAmo->Multiply(hMask2dEAmo,hMask2dEUser);
      hMask2dE0->Multiply(hMask2dE0,hMask2dEUser);
      hMask2dE45->Multiply(hMask2dE45,hMask2dEUser);
      hMask2dE90->Multiply(hMask2dE90,hMask2dEUser);
      hMask2dE135->Multiply(hMask2dE135,hMask2dEUser);
    }

    char hString[120];
    double par0=0,par1=0,par2=0,par3=0,par4=0,par5=0;
    double err0=0,err1=0,err2=0,err3=0,err4=0,err5=0;
    int nBy = hPhiCutVsRun->GetNbinsY();
    double yLow = hPhiCutVsRun->GetYaxis()->GetXmin();
    double yHigh = hPhiCutVsRun->GetYaxis()->GetXmax();

    double centerShift = 0.1;
    TH1D *hPol0   = new TH1D("hPol0","",25,7.0-centerShift,12.0-centerShift);   hPol0->Sumw2();
    TH1D *hPol45  = new TH1D("hPol45","",25,7.0-centerShift,12.0-centerShift);  hPol45->Sumw2();
    TH1D *hPol90  = new TH1D("hPol90","",25,7.0-centerShift,12.0-centerShift);  hPol90->Sumw2();
    TH1D *hPol135 = new TH1D("hPol135","",25,7.0-centerShift,12.0-centerShift); hPol135->Sumw2();

    hPol0->SetLineColor(1);
    hPol45->SetLineColor(2);
    hPol90->SetLineColor(3);
    hPol135->SetLineColor(4);
 

    TF1 *fitFun = new TF1("fitFun", "[0]*(1-[1]*cos(2*(3.14/180)*(x-[2]))+[3]*cos(3.14*(x-[4])/180))",0.0,360.0);

    double phiOffset = 5;
    //START 0
    sprintf(hString,"hPhiCut0");
    TH1D *hPhiCut0 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0->Sumw2();
    //fitFun->FixParameter(2,10.0);
    fitFun->SetParameter(2,phiOffset);

    if (fitShow == 0) {
      phi1D(hPhiCutVsRun,hMask2dE0,hPhiCut0,hString,fitFun,false);
      //hPhiCut0->SetMinimum(45000);hPhiCut0->SetMaximum(60000);
      hPhiCut0->Draw();
      return tpolToGo;
    }
    phi1D(hPhiCutVsRun,hMask2dE0,hPhiCut0,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    par2 = fitFun->GetParameter(2); err2 = fitFun->GetParError(2);
    double polVal0 = par1/anaPower;
    double polErr0 = err1/anaPower;

    sprintf(hString,"hPhiCut45");
    TH1D *hPhiCut45 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45->Sumw2();
    //fitFun->FixParameter(2,45.0+10.0);
    fitFun->SetParameter(2,45.0+phiOffset);
    if (fitShow == 45) {
      phi1D(hPhiCutVsRun,hMask2dE45,hPhiCut45,hString,fitFun,false);
      hPhiCut45->Draw();
      return tpolToGo;
    }
    phi1D(hPhiCutVsRun,hMask2dE45,hPhiCut45,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    par2 = fitFun->GetParameter(2); err2 = fitFun->GetParError(2);
    double polVal45 = par1/anaPower;
    double polErr45 = err1/anaPower;

    sprintf(hString,"hPhiCut90");
    //fitFun->FixParameter(2,90+10);
    TH1D *hPhiCut90 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90->Sumw2();
    fitFun->SetParameter(2,90+phiOffset);
    if (fitShow == 90) {
      phi1D(hPhiCutVsRun,hMask2dE90,hPhiCut90,hString,fitFun,false);
      //hPhiCut90->SetMinimum(70000);hPhiCut90->SetMaximum(95000);
      hPhiCut90->Draw();
      return tpolToGo;
    }
    phi1D(hPhiCutVsRun,hMask2dE90,hPhiCut90,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    par2 = fitFun->GetParameter(2); err2 = fitFun->GetParError(2);
    double polVal90 = par1/anaPower;
    double polErr90 = err1/anaPower;

    sprintf(hString,"hPhiCut135");
    TH1D *hPhiCut135 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135->Sumw2();
    //fitFun->FixParameter(2,135+10.0);
    fitFun->SetParameter(2,135+phiOffset-180);
    if (fitShow == 135) {
      phi1D(hPhiCutVsRun,hMask2dE135,hPhiCut135,hString,fitFun,false);
      //hPhiCut135->SetMinimum(70000);hPhiCut135->SetMaximum(95000);
      hPhiCut135->Draw();
      return tpolToGo;
    }
    phi1D(hPhiCutVsRun,hMask2dE135,hPhiCut135,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    par2 = fitFun->GetParameter(2); err2 = fitFun->GetParError(2);
    //cout<<"par2(135) = "<<par2<<" +/- "<<err2<<endl;
    double polVal135 = par1/anaPower;
    double polErr135 = err1/anaPower;

    phiPars_t phiPars;
    bool simFit = true;
    if (simFit == true && (sp20Set == true || y21Set == true)) {
      cout<<""<<endl;
      cout<<"Starting SimFit"<<endl;

      //E_gamma between 8.0 and 8.6 GeV
      phiHistos_t phiHistos;
      phiHistos.hPhi0   = hPhiCut0;
      phiHistos.hPhi45  = hPhiCut45;
      phiHistos.hPhi90  = hPhiCut90;
      phiHistos.hPhi135 = hPhiCut135;

      phiPars = phiSimFit(phiHistos);

      cout<<"Ended SimFit"<<endl;
      cout<<""<<endl;

      fitFun->FixParameter(3,phiPars.aCenter);
      fitFun->FixParameter(4,phiPars.aOffset);

    }




    //START 0
    if (simFit == true)fitFun->FixParameter(2,phiPars.offset+0);
    sprintf(hString,"hPhiCut0AtEg7100");
    TH1D *hPhiCut0AtEg7100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg7100->Sumw2();
    phi1D(hPhiCutVsRunAtEg7100,hMask2dE0,hPhiCut0AtEg7100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(1,fabs(par1)); 
    hPol0->SetBinError(1,err1);
    sprintf(hString,"hPhiCut0AtEg7300");
    TH1D *hPhiCut0AtEg7300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg7300->Sumw2();
    phi1D(hPhiCutVsRunAtEg7300,hMask2dE0,hPhiCut0AtEg7300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(2,fabs(par1)); 
    hPol0->SetBinError(2,err1);
    sprintf(hString,"hPhiCut0AtEg7500");
    TH1D *hPhiCut0AtEg7500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg7500->Sumw2();
    phi1D(hPhiCutVsRunAtEg7500,hMask2dE0,hPhiCut0AtEg7500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(3,fabs(par1)); 
    hPol0->SetBinError(3,err1);
    sprintf(hString,"hPhiCut0AtEg7700");
    TH1D *hPhiCut0AtEg7700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg7700->Sumw2();
    phi1D(hPhiCutVsRunAtEg7700,hMask2dE0,hPhiCut0AtEg7700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(4,fabs(par1)); 
    hPol0->SetBinError(4,err1);
    sprintf(hString,"hPhiCut0AtEg7900");
    TH1D *hPhiCut0AtEg7900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg7900->Sumw2();
    phi1D(hPhiCutVsRunAtEg7900,hMask2dE0,hPhiCut0AtEg7900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(5,fabs(par1)); 
    hPol0->SetBinError(5,err1);
    
    sprintf(hString,"hPhiCut0AtEg8100");
    TH1D *hPhiCut0AtEg8100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg8100->Sumw2();
    phi1D(hPhiCutVsRunAtEg8100,hMask2dE0,hPhiCut0AtEg8100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(6,fabs(par1)); 
    hPol0->SetBinError(6,err1);
    sprintf(hString,"hPhiCut0AtEg8300");
    TH1D *hPhiCut0AtEg8300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg8300->Sumw2();
    phi1D(hPhiCutVsRunAtEg8300,hMask2dE0,hPhiCut0AtEg8300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(7,fabs(par1)); 
    hPol0->SetBinError(7,err1);
    sprintf(hString,"hPhiCut0AtEg8500");
    TH1D *hPhiCut0AtEg8500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg8500->Sumw2();
    phi1D(hPhiCutVsRunAtEg8500,hMask2dE0,hPhiCut0AtEg8500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(8,fabs(par1)); 
    hPol0->SetBinError(8,err1);
    sprintf(hString,"hPhiCut0AtEg8700");
    TH1D *hPhiCut0AtEg8700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg8700->Sumw2();
    phi1D(hPhiCutVsRunAtEg8700,hMask2dE0,hPhiCut0AtEg8700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(9,fabs(par1)); 
    hPol0->SetBinError(9,err1);
    sprintf(hString,"hPhiCut0AtEg8900");
    TH1D *hPhiCut0AtEg8900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg8900->Sumw2();
    phi1D(hPhiCutVsRunAtEg8900,hMask2dE0,hPhiCut0AtEg8900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(10,fabs(par1)); 
    hPol0->SetBinError(10,err1);
    
    sprintf(hString,"hPhiCut0AtEg9100");
    TH1D *hPhiCut0AtEg9100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg9100->Sumw2();
    phi1D(hPhiCutVsRunAtEg9100,hMask2dE0,hPhiCut0AtEg9100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(11,fabs(par1)); 
    hPol0->SetBinError(11,err1);
    sprintf(hString,"hPhiCut0AtEg9300");
    TH1D *hPhiCut0AtEg9300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg9300->Sumw2();
    phi1D(hPhiCutVsRunAtEg9300,hMask2dE0,hPhiCut0AtEg9300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(12,fabs(par1)); 
    hPol0->SetBinError(12,err1);
    sprintf(hString,"hPhiCut0AtEg9500");
    TH1D *hPhiCut0AtEg9500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg9500->Sumw2();
    phi1D(hPhiCutVsRunAtEg9500,hMask2dE0,hPhiCut0AtEg9500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(13,fabs(par1)); 
    hPol0->SetBinError(13,err1);
    sprintf(hString,"hPhiCut0AtEg9700");
    TH1D *hPhiCut0AtEg9700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg9700->Sumw2();
    phi1D(hPhiCutVsRunAtEg9700,hMask2dE0,hPhiCut0AtEg9700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(14,fabs(par1)); 
    hPol0->SetBinError(14,err1);
    sprintf(hString,"hPhiCut0AtEg9900");
    TH1D *hPhiCut0AtEg9900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg9900->Sumw2();
    phi1D(hPhiCutVsRunAtEg9900,hMask2dE0,hPhiCut0AtEg9900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(15,fabs(par1)); 
    hPol0->SetBinError(15,err1);
    
    sprintf(hString,"hPhiCut0AtEg10100");
    TH1D *hPhiCut0AtEg10100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg10100->Sumw2();
    phi1D(hPhiCutVsRunAtEg10100,hMask2dE0,hPhiCut0AtEg10100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(16,fabs(par1)); 
    hPol0->SetBinError(16,err1);
    sprintf(hString,"hPhiCut0AtEg10300");
    TH1D *hPhiCut0AtEg10300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg10300->Sumw2();
    phi1D(hPhiCutVsRunAtEg10300,hMask2dE0,hPhiCut0AtEg10300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(17,fabs(par1)); 
    hPol0->SetBinError(17,err1);
    sprintf(hString,"hPhiCut0AtEg10500");
    TH1D *hPhiCut0AtEg10500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg10500->Sumw2();
    phi1D(hPhiCutVsRunAtEg10500,hMask2dE0,hPhiCut0AtEg10500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(18,fabs(par1)); 
    hPol0->SetBinError(18,err1);
    sprintf(hString,"hPhiCut0AtEg10700");
    TH1D *hPhiCut0AtEg10700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg10700->Sumw2();
    phi1D(hPhiCutVsRunAtEg10700,hMask2dE0,hPhiCut0AtEg10700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(19,fabs(par1)); 
    hPol0->SetBinError(19,err1);
    sprintf(hString,"hPhiCut0AtEg10900");
    TH1D *hPhiCut0AtEg10900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg10900->Sumw2();
    phi1D(hPhiCutVsRunAtEg10900,hMask2dE0,hPhiCut0AtEg10900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(20,fabs(par1)); 
    hPol0->SetBinError(20,err1);
    
    sprintf(hString,"hPhiCut0AtEg11100");
    TH1D *hPhiCut0AtEg11100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg11100->Sumw2();
    phi1D(hPhiCutVsRunAtEg11100,hMask2dE0,hPhiCut0AtEg11100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(21,fabs(par1)); 
    hPol0->SetBinError(21,err1);
    sprintf(hString,"hPhiCut0AtEg11300");
    TH1D *hPhiCut0AtEg11300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut0AtEg11300->Sumw2();
    phi1D(hPhiCutVsRunAtEg11300,hMask2dE0,hPhiCut0AtEg11300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol0->SetBinContent(22,fabs(par1)); 
    hPol0->SetBinError(22,err1);
    hPol0->GetXaxis()->SetRangeUser(7.8,11.4);
    
      
    //START 45
    if (simFit == true)fitFun->FixParameter(2,phiPars.offset+45);
    sprintf(hString,"hPhiCut45AtEg7100");
    TH1D *hPhiCut45AtEg7100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg7100->Sumw2();
    phi1D(hPhiCutVsRunAtEg7100,hMask2dE45,hPhiCut45AtEg7100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(1,fabs(par1)); 
    hPol45->SetBinError(1,err1);
    sprintf(hString,"hPhiCut45AtEg7300");
    TH1D *hPhiCut45AtEg7300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg7300->Sumw2();
    phi1D(hPhiCutVsRunAtEg7300,hMask2dE45,hPhiCut45AtEg7300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(2,fabs(par1)); 
    hPol45->SetBinError(2,err1);
    sprintf(hString,"hPhiCut45AtEg7500");
    TH1D *hPhiCut45AtEg7500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg7500->Sumw2();
    phi1D(hPhiCutVsRunAtEg7500,hMask2dE45,hPhiCut45AtEg7500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(3,fabs(par1)); 
    hPol45->SetBinError(3,err1);
    sprintf(hString,"hPhiCut45AtEg7700");
    TH1D *hPhiCut45AtEg7700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg7700->Sumw2();
    phi1D(hPhiCutVsRunAtEg7700,hMask2dE45,hPhiCut45AtEg7700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(4,fabs(par1)); 
    hPol45->SetBinError(4,err1);
    sprintf(hString,"hPhiCut45AtEg7900");
    TH1D *hPhiCut45AtEg7900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg7900->Sumw2();
    phi1D(hPhiCutVsRunAtEg7900,hMask2dE45,hPhiCut45AtEg7900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(5,fabs(par1)); 
    hPol45->SetBinError(5,err1);
    
    sprintf(hString,"hPhiCut45AtEg8100");
    TH1D *hPhiCut45AtEg8100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg8100->Sumw2();
    phi1D(hPhiCutVsRunAtEg8100,hMask2dE45,hPhiCut45AtEg8100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(6,fabs(par1)); 
    hPol45->SetBinError(6,err1);
    sprintf(hString,"hPhiCut45AtEg8300");
    TH1D *hPhiCut45AtEg8300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg8300->Sumw2();
    phi1D(hPhiCutVsRunAtEg8300,hMask2dE45,hPhiCut45AtEg8300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(7,fabs(par1)); 
    hPol45->SetBinError(7,err1);
    sprintf(hString,"hPhiCut45AtEg8500");
    TH1D *hPhiCut45AtEg8500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg8500->Sumw2();
    phi1D(hPhiCutVsRunAtEg8500,hMask2dE45,hPhiCut45AtEg8500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(8,fabs(par1)); 
    hPol45->SetBinError(8,err1);
    sprintf(hString,"hPhiCut45AtEg8700");
    TH1D *hPhiCut45AtEg8700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg8700->Sumw2();
    phi1D(hPhiCutVsRunAtEg8700,hMask2dE45,hPhiCut45AtEg8700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(9,fabs(par1)); 
    hPol45->SetBinError(9,err1);
    sprintf(hString,"hPhiCut45AtEg8900");
    TH1D *hPhiCut45AtEg8900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg8900->Sumw2();
    phi1D(hPhiCutVsRunAtEg8900,hMask2dE45,hPhiCut45AtEg8900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(10,fabs(par1)); 
    hPol45->SetBinError(10,err1);
    
    sprintf(hString,"hPhiCut45AtEg9100");
    TH1D *hPhiCut45AtEg9100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg9100->Sumw2();
    phi1D(hPhiCutVsRunAtEg9100,hMask2dE45,hPhiCut45AtEg9100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(11,fabs(par1)); 
    hPol45->SetBinError(11,err1);
    sprintf(hString,"hPhiCut45AtEg9300");
    TH1D *hPhiCut45AtEg9300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg9300->Sumw2();
    phi1D(hPhiCutVsRunAtEg9300,hMask2dE45,hPhiCut45AtEg9300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(12,fabs(par1)); 
    hPol45->SetBinError(12,err1);
    sprintf(hString,"hPhiCut45AtEg9500");
    TH1D *hPhiCut45AtEg9500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg9500->Sumw2();
    phi1D(hPhiCutVsRunAtEg9500,hMask2dE45,hPhiCut45AtEg9500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(13,fabs(par1)); 
    hPol45->SetBinError(13,err1);
    sprintf(hString,"hPhiCut45AtEg9700");
    TH1D *hPhiCut45AtEg9700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg9700->Sumw2();
    phi1D(hPhiCutVsRunAtEg9700,hMask2dE45,hPhiCut45AtEg9700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(14,fabs(par1)); 
    hPol45->SetBinError(14,err1);
    sprintf(hString,"hPhiCut45AtEg9900");
    TH1D *hPhiCut45AtEg9900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg9900->Sumw2();
    phi1D(hPhiCutVsRunAtEg9900,hMask2dE45,hPhiCut45AtEg9900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(15,fabs(par1)); 
    hPol45->SetBinError(15,err1);
    
    sprintf(hString,"hPhiCut45AtEg10100");
    TH1D *hPhiCut45AtEg10100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg10100->Sumw2();
    phi1D(hPhiCutVsRunAtEg10100,hMask2dE45,hPhiCut45AtEg10100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(16,fabs(par1)); 
    hPol45->SetBinError(16,err1);
    sprintf(hString,"hPhiCut45AtEg10300");
    TH1D *hPhiCut45AtEg10300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg10300->Sumw2();
    phi1D(hPhiCutVsRunAtEg10300,hMask2dE45,hPhiCut45AtEg10300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(17,fabs(par1)); 
    hPol45->SetBinError(17,err1);
    sprintf(hString,"hPhiCut45AtEg10500");
    TH1D *hPhiCut45AtEg10500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg10500->Sumw2();
    phi1D(hPhiCutVsRunAtEg10500,hMask2dE45,hPhiCut45AtEg10500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(18,fabs(par1)); 
    hPol45->SetBinError(18,err1);
    sprintf(hString,"hPhiCut45AtEg10700");
    TH1D *hPhiCut45AtEg10700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg10700->Sumw2();
    phi1D(hPhiCutVsRunAtEg10700,hMask2dE45,hPhiCut45AtEg10700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(19,fabs(par1)); 
    hPol45->SetBinError(19,err1);
    sprintf(hString,"hPhiCut45AtEg10900");
    TH1D *hPhiCut45AtEg10900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg10900->Sumw2();
    phi1D(hPhiCutVsRunAtEg10900,hMask2dE45,hPhiCut45AtEg10900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(20,fabs(par1)); 
    hPol45->SetBinError(20,err1);
    
    sprintf(hString,"hPhiCut45AtEg11100");
    TH1D *hPhiCut45AtEg11100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg11100->Sumw2();
    phi1D(hPhiCutVsRunAtEg11100,hMask2dE45,hPhiCut45AtEg11100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(21,fabs(par1)); 
    hPol45->SetBinError(21,err1);
    sprintf(hString,"hPhiCut45AtEg11300");
    TH1D *hPhiCut45AtEg11300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut45AtEg11300->Sumw2();
    phi1D(hPhiCutVsRunAtEg11300,hMask2dE45,hPhiCut45AtEg11300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol45->SetBinContent(22,fabs(par1)); 
    hPol45->SetBinError(22,err1);
    
    hPol45->GetXaxis()->SetRangeUser(7.5,11.4);
    
    //START 90
    if (simFit == true)fitFun->FixParameter(2,phiPars.offset+90);  
    sprintf(hString,"hPhiCut90AtEg7100");
    TH1D *hPhiCut90AtEg7100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg7100->Sumw2();
    phi1D(hPhiCutVsRunAtEg7100,hMask2dE90,hPhiCut90AtEg7100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(1,fabs(par1)); 
    hPol90->SetBinError(1,err1);
    sprintf(hString,"hPhiCut90AtEg7300");
    TH1D *hPhiCut90AtEg7300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg7300->Sumw2();
    phi1D(hPhiCutVsRunAtEg7300,hMask2dE90,hPhiCut90AtEg7300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(2,fabs(par1)); 
    hPol90->SetBinError(2,err1);
    sprintf(hString,"hPhiCut90AtEg7500");
    TH1D *hPhiCut90AtEg7500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg7500->Sumw2();
    phi1D(hPhiCutVsRunAtEg7500,hMask2dE90,hPhiCut90AtEg7500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(3,fabs(par1)); 
    hPol90->SetBinError(3,err1);
    sprintf(hString,"hPhiCut90AtEg7700");
    TH1D *hPhiCut90AtEg7700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg7700->Sumw2();
    phi1D(hPhiCutVsRunAtEg7700,hMask2dE90,hPhiCut90AtEg7700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(4,fabs(par1)); 
    hPol90->SetBinError(4,err1);
    sprintf(hString,"hPhiCut90AtEg7900");
    TH1D *hPhiCut90AtEg7900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg7900->Sumw2();
    phi1D(hPhiCutVsRunAtEg7900,hMask2dE90,hPhiCut90AtEg7900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(5,fabs(par1)); 
    hPol90->SetBinError(5,err1);
    hPol90->GetXaxis()->SetRangeUser(7.6,11.4);
    
    sprintf(hString,"hPhiCut90AtEg8100");
    TH1D *hPhiCut90AtEg8100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg8100->Sumw2();
    phi1D(hPhiCutVsRunAtEg8100,hMask2dE90,hPhiCut90AtEg8100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(6,fabs(par1)); 
    hPol90->SetBinError(6,err1);
    sprintf(hString,"hPhiCut90AtEg8300");
    TH1D *hPhiCut90AtEg8300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg8300->Sumw2();
    phi1D(hPhiCutVsRunAtEg8300,hMask2dE90,hPhiCut90AtEg8300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(7,fabs(par1)); 
    hPol90->SetBinError(7,err1);
    sprintf(hString,"hPhiCut90AtEg8500");
    TH1D *hPhiCut90AtEg8500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg8500->Sumw2();
    phi1D(hPhiCutVsRunAtEg8500,hMask2dE90,hPhiCut90AtEg8500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(8,fabs(par1)); 
    hPol90->SetBinError(8,err1);
    sprintf(hString,"hPhiCut90AtEg8700");
    TH1D *hPhiCut90AtEg8700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg8700->Sumw2();
    phi1D(hPhiCutVsRunAtEg8700,hMask2dE90,hPhiCut90AtEg8700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(9,fabs(par1)); 
    hPol90->SetBinError(9,err1);
    sprintf(hString,"hPhiCut90AtEg8900");
    TH1D *hPhiCut90AtEg8900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg8900->Sumw2();
    phi1D(hPhiCutVsRunAtEg8900,hMask2dE90,hPhiCut90AtEg8900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(10,fabs(par1)); 
    hPol90->SetBinError(10,err1);
    
    sprintf(hString,"hPhiCut90AtEg9100");
    TH1D *hPhiCut90AtEg9100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg9100->Sumw2();
    phi1D(hPhiCutVsRunAtEg9100,hMask2dE90,hPhiCut90AtEg9100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(11,fabs(par1)); 
    hPol90->SetBinError(11,err1);
    sprintf(hString,"hPhiCut90AtEg9300");
    TH1D *hPhiCut90AtEg9300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg9300->Sumw2();
    phi1D(hPhiCutVsRunAtEg9300,hMask2dE90,hPhiCut90AtEg9300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(12,fabs(par1)); 
    hPol90->SetBinError(12,err1);
    sprintf(hString,"hPhiCut90AtEg9500");
    TH1D *hPhiCut90AtEg9500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg9500->Sumw2();
    phi1D(hPhiCutVsRunAtEg9500,hMask2dE90,hPhiCut90AtEg9500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(13,fabs(par1)); 
    hPol90->SetBinError(13,err1);
    sprintf(hString,"hPhiCut90AtEg9700");
    TH1D *hPhiCut90AtEg9700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg9700->Sumw2();
    phi1D(hPhiCutVsRunAtEg9700,hMask2dE90,hPhiCut90AtEg9700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(14,fabs(par1)); 
    hPol90->SetBinError(14,err1);
    sprintf(hString,"hPhiCut90AtEg9900");
    TH1D *hPhiCut90AtEg9900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg9900->Sumw2();
    phi1D(hPhiCutVsRunAtEg9900,hMask2dE90,hPhiCut90AtEg9900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(15,fabs(par1)); 
    hPol90->SetBinError(15,err1);
    
    sprintf(hString,"hPhiCut90AtEg10100");
    TH1D *hPhiCut90AtEg10100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg10100->Sumw2();
    phi1D(hPhiCutVsRunAtEg10100,hMask2dE90,hPhiCut90AtEg10100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(16,fabs(par1)); 
    hPol90->SetBinError(16,err1);
    sprintf(hString,"hPhiCut90AtEg10300");
    TH1D *hPhiCut90AtEg10300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg10300->Sumw2();
    phi1D(hPhiCutVsRunAtEg10300,hMask2dE90,hPhiCut90AtEg10300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(17,fabs(par1)); 
    hPol90->SetBinError(17,err1);
    sprintf(hString,"hPhiCut90AtEg10500");
    TH1D *hPhiCut90AtEg10500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg10500->Sumw2();
    phi1D(hPhiCutVsRunAtEg10500,hMask2dE90,hPhiCut90AtEg10500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(18,fabs(par1)); 
    hPol90->SetBinError(18,err1);
    sprintf(hString,"hPhiCut90AtEg10700");
    TH1D *hPhiCut90AtEg10700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg10700->Sumw2();
    phi1D(hPhiCutVsRunAtEg10700,hMask2dE90,hPhiCut90AtEg10700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(19,fabs(par1)); 
    hPol90->SetBinError(19,err1);
    sprintf(hString,"hPhiCut90AtEg10900");
    TH1D *hPhiCut90AtEg10900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg10900->Sumw2();
    phi1D(hPhiCutVsRunAtEg10900,hMask2dE90,hPhiCut90AtEg10900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(20,fabs(par1)); 
    hPol90->SetBinError(20,err1);
    
    sprintf(hString,"hPhiCut90AtEg11100");
    TH1D *hPhiCut90AtEg11100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg11100->Sumw2();
    phi1D(hPhiCutVsRunAtEg11100,hMask2dE90,hPhiCut90AtEg11100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(21,fabs(par1)); 
    hPol90->SetBinError(21,err1);
    sprintf(hString,"hPhiCut90AtEg11300");
    TH1D *hPhiCut90AtEg11300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut90AtEg11300->Sumw2();
    phi1D(hPhiCutVsRunAtEg11300,hMask2dE90,hPhiCut90AtEg11300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol90->SetBinContent(22,fabs(par1)); 
    hPol90->SetBinError(22,err1);
    hPol90->GetXaxis()->SetRangeUser(7.6,11.4);
    
    //START 135
    if (simFit == true)fitFun->FixParameter(2,phiPars.offset+135);   
    sprintf(hString,"hPhiCut135AtEg7100");
    TH1D *hPhiCut135AtEg7100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg7100->Sumw2();
    phi1D(hPhiCutVsRunAtEg7100,hMask2dE135,hPhiCut135AtEg7100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(1,fabs(par1)); 
    hPol135->SetBinError(1,err1);
    sprintf(hString,"hPhiCut135AtEg7300");
    TH1D *hPhiCut135AtEg7300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg7300->Sumw2();
    phi1D(hPhiCutVsRunAtEg7300,hMask2dE135,hPhiCut135AtEg7300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(2,fabs(par1)); 
    hPol135->SetBinError(2,err1);
    sprintf(hString,"hPhiCut135AtEg7500");
    TH1D *hPhiCut135AtEg7500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg7500->Sumw2();
    phi1D(hPhiCutVsRunAtEg7500,hMask2dE135,hPhiCut135AtEg7500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(3,fabs(par1)); 
    hPol135->SetBinError(3,err1);
    sprintf(hString,"hPhiCut135AtEg7700");
    TH1D *hPhiCut135AtEg7700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg7700->Sumw2();
    phi1D(hPhiCutVsRunAtEg7700,hMask2dE135,hPhiCut135AtEg7700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(4,fabs(par1)); 
    hPol135->SetBinError(4,err1);
    sprintf(hString,"hPhiCut135AtEg7900");
    TH1D *hPhiCut135AtEg7900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg7900->Sumw2();
    phi1D(hPhiCutVsRunAtEg7900,hMask2dE135,hPhiCut135AtEg7900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(5,fabs(par1)); 
    hPol135->SetBinError(5,err1);
    
    sprintf(hString,"hPhiCut135AtEg8100");
    TH1D *hPhiCut135AtEg8100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg8100->Sumw2();
    phi1D(hPhiCutVsRunAtEg8100,hMask2dE135,hPhiCut135AtEg8100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(6,fabs(par1)); 
    hPol135->SetBinError(6,err1);
    sprintf(hString,"hPhiCut135AtEg8300");
    TH1D *hPhiCut135AtEg8300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg8300->Sumw2();
    phi1D(hPhiCutVsRunAtEg8300,hMask2dE135,hPhiCut135AtEg8300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(7,fabs(par1)); 
    hPol135->SetBinError(7,err1);
    sprintf(hString,"hPhiCut135AtEg8500");
    TH1D *hPhiCut135AtEg8500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg8500->Sumw2();
    phi1D(hPhiCutVsRunAtEg8500,hMask2dE135,hPhiCut135AtEg8500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(8,fabs(par1)); 
    hPol135->SetBinError(8,err1);
    sprintf(hString,"hPhiCut135AtEg8700");
    TH1D *hPhiCut135AtEg8700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg8700->Sumw2();
    phi1D(hPhiCutVsRunAtEg8700,hMask2dE135,hPhiCut135AtEg8700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(9,fabs(par1)); 
    hPol135->SetBinError(9,err1);
    sprintf(hString,"hPhiCut135AtEg8900");
    TH1D *hPhiCut135AtEg8900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg8900->Sumw2();
    phi1D(hPhiCutVsRunAtEg8900,hMask2dE135,hPhiCut135AtEg8900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(10,fabs(par1)); 
    hPol135->SetBinError(10,err1);
    
    sprintf(hString,"hPhiCut135AtEg9100");
    TH1D *hPhiCut135AtEg9100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg9100->Sumw2();
    phi1D(hPhiCutVsRunAtEg9100,hMask2dE135,hPhiCut135AtEg9100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(11,fabs(par1)); 
    hPol135->SetBinError(11,err1);
    sprintf(hString,"hPhiCut135AtEg9300");
    TH1D *hPhiCut135AtEg9300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg9300->Sumw2();
    phi1D(hPhiCutVsRunAtEg9300,hMask2dE135,hPhiCut135AtEg9300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(12,fabs(par1)); 
    hPol135->SetBinError(12,err1);
    sprintf(hString,"hPhiCut135AtEg9500");
    TH1D *hPhiCut135AtEg9500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg9500->Sumw2();
    phi1D(hPhiCutVsRunAtEg9500,hMask2dE135,hPhiCut135AtEg9500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(13,fabs(par1)); 
    hPol135->SetBinError(13,err1);
    sprintf(hString,"hPhiCut135AtEg9700");
    TH1D *hPhiCut135AtEg9700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg9700->Sumw2();
    phi1D(hPhiCutVsRunAtEg9700,hMask2dE135,hPhiCut135AtEg9700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(14,fabs(par1)); 
    hPol135->SetBinError(14,err1);
    sprintf(hString,"hPhiCut135AtEg9900");
    TH1D *hPhiCut135AtEg9900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg9900->Sumw2();
    phi1D(hPhiCutVsRunAtEg9900,hMask2dE135,hPhiCut135AtEg9900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(15,fabs(par1)); 
    hPol135->SetBinError(15,err1);
      
    sprintf(hString,"hPhiCut135AtEg10100");
    TH1D *hPhiCut135AtEg10100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg10100->Sumw2();
    phi1D(hPhiCutVsRunAtEg10100,hMask2dE135,hPhiCut135AtEg10100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(16,fabs(par1)); 
    hPol135->SetBinError(16,err1);
    sprintf(hString,"hPhiCut135AtEg10300");
    TH1D *hPhiCut135AtEg10300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg10300->Sumw2();
    phi1D(hPhiCutVsRunAtEg10300,hMask2dE135,hPhiCut135AtEg10300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(17,fabs(par1)); 
    hPol135->SetBinError(17,err1);
    sprintf(hString,"hPhiCut135AtEg10500");
    TH1D *hPhiCut135AtEg10500 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg10500->Sumw2();
    phi1D(hPhiCutVsRunAtEg10500,hMask2dE135,hPhiCut135AtEg10500,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(18,fabs(par1)); 
    hPol135->SetBinError(18,err1);
    sprintf(hString,"hPhiCut135AtEg10700");
    TH1D *hPhiCut135AtEg10700 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg10700->Sumw2();
    phi1D(hPhiCutVsRunAtEg10700,hMask2dE135,hPhiCut135AtEg10700,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(19,fabs(par1)); 
    hPol135->SetBinError(19,err1);
    sprintf(hString,"hPhiCut135AtEg10900");
    TH1D *hPhiCut135AtEg10900 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg10900->Sumw2();
    phi1D(hPhiCutVsRunAtEg10900,hMask2dE135,hPhiCut135AtEg10900,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(20,fabs(par1)); 
    hPol135->SetBinError(20,err1);
    
    sprintf(hString,"hPhiCut135AtEg11100");
    TH1D *hPhiCut135AtEg11100 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg11100->Sumw2();
    phi1D(hPhiCutVsRunAtEg11100,hMask2dE135,hPhiCut135AtEg11100,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(21,fabs(par1)); 
    hPol135->SetBinError(21,err1);
    sprintf(hString,"hPhiCut135AtEg11300");
    TH1D *hPhiCut135AtEg11300 = new TH1D(hString,"",nBy,yLow,yHigh); hPhiCut135AtEg11300->Sumw2();
    phi1D(hPhiCutVsRunAtEg11300,hMask2dE135,hPhiCut135AtEg11300,hString,fitFun,true);
    par1 = fitFun->GetParameter(1); err1 = fitFun->GetParError(1);
    hPol135->SetBinContent(22,fabs(par1)); 
    hPol135->SetBinError(22,err1);
    hPol135->GetXaxis()->SetRangeUser(7.6,11.4);
    //TH2D *hZo135 = (TH2D*)hZVal->Clone("hZo135");
    //hZAmo->Multiply(hZAmo,hMask2dEAmo);

    hPol0->Scale(1.0/anaPower);
    hPol45->Scale(1.0/anaPower);
    hPol90->Scale(1.0/anaPower);
    hPol135->Scale(1.0/anaPower);

    if (sp20Set == true){
      sprintf(fString,"./outFiles/runList0Sp20.txt");
      writeRuns(fString,hMask0);
      sprintf(fString,"./outFiles/runList45Sp20.txt");
      writeRuns(fString,hMask45);
      sprintf(fString,"./outFiles/runList90Sp20.txt");
      writeRuns(fString,hMask90);
      sprintf(fString,"./outFiles/runList135Sp20.txt");
      writeRuns(fString,hMask135);

      char  pString[120];
      cout<<""<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<"                  Run lists 2020                         "<<endl;
      cout<<"---------------------------------------------------------"<<endl;
      cout<<"\tBeam orientation\tFile"<<endl;
      cout<<"\t  0 degrees:\t./outFiles/runList0Sp20.txt"<<endl;
      cout<<"\t 45 degrees:\t./outFiles/runList45Sp20.txt"<<endl;
      cout<<"\t 90 degrees:\t./outFiles/runList90Sp20.txt"<<endl;
      cout<<"\t135 degrees:\t./outFiles/runList135Sp20.txt"<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<" Polarization values for E_gamma between 8.0 and 8.6 GeV"<<endl;
      cout<<"---------------------------------------------------------"<<endl;
      cout<<"\tBeam orientation\tPolarization"<<endl;
      sprintf(pString,"\t  0 degrees:\t     %5.4f +/- %5.4f",polVal0,polErr0);
      cout<<pString<<endl;
      sprintf(pString,"\t 45 degrees:\t     %5.4f +/- %5.4f",polVal45,polErr45);
      cout<<pString<<endl;
      sprintf(pString,"\t 90 degrees:\t     %5.4f +/- %5.4f",polVal90,polErr90);
      cout<<pString<<endl;
      sprintf(pString,"\t135 degrees:\t     %5.4f +/- %5.4f",polVal135,polErr135);
      cout<<pString<<endl;
      cout<<"*********************************************************"<<endl;
      if (userSet == true) {
	cout<<"NOTE: Run set defined by user in file ./inFiles/userSet.txt"<<endl;
      }else{
	cout<<"NOTE: Run set = all runs"<<endl;
      }
      cout<<"NOTE: Histograms written to ./outFiles/sp20TPol.root"<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<""<<endl;
    }

    if (y21Set == true){
      sprintf(fString,"./outFiles/runList0_2021-11.txt");
      writeRuns(fString,hMask0);
      sprintf(fString,"./outFiles/runList45_2021-11.txt");
      writeRuns(fString,hMask45);
      sprintf(fString,"./outFiles/runList90_2021-11.txt");
      writeRuns(fString,hMask90);
      sprintf(fString,"./outFiles/runList135_2021-11.txt");
      writeRuns(fString,hMask135);

      char  pString[120];
      cout<<""<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<"                  Run lists 2021-11                      "<<endl;
      cout<<"---------------------------------------------------------"<<endl;
      cout<<"\tBeam orientation\tFile"<<endl;
      cout<<"\t  0 degrees:\t./outFiles/runList0_2021-11.txt"<<endl;
      cout<<"\t 45 degrees:\t./outFiles/runList45_2021-11.txt"<<endl;
      cout<<"\t 90 degrees:\t./outFiles/runList90_2021-11.txt"<<endl;
      cout<<"\t135 degrees:\t./outFiles/runList135_2021-11.txt"<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<" Polarization values for E_gamma between 7.9 and 8.5 GeV"<<endl;
      cout<<"---------------------------------------------------------"<<endl;
      cout<<"\tBeam orientation\tPolarization"<<endl;
      sprintf(pString,"\t  0 degrees:\t     %5.4f +/- %5.4f",fabs(polVal0),polErr0);
      cout<<pString<<endl;
      sprintf(pString,"\t 45 degrees:\t     %5.4f +/- %5.4f",fabs(polVal45),polErr45);
      cout<<pString<<endl;
      sprintf(pString,"\t 90 degrees:\t     %5.4f +/- %5.4f",fabs(polVal90),polErr90);
      cout<<pString<<endl;
      sprintf(pString,"\t135 degrees:\t     %5.4f +/- %5.4f",fabs(polVal135),polErr135);
      cout<<pString<<endl;
      cout<<"*********************************************************"<<endl;
      if (userSet == true) {
	cout<<"NOTE: Run set defined by user in file ./inFiles/userSet.txt"<<endl;
      }else{
	cout<<"NOTE: Run set = all runs"<<endl;
      }
      cout<<"NOTE: Histograms written to ./outFiles/2021-11_TPol.root"<<endl;
      cout<<"*********************************************************"<<endl;
      cout<<""<<endl;
    }

    if (mixMode == false){
      TCanvas *c1 = new TCanvas();
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1)->SetPad(0.0,0.5,0.5,1.0);
      c1->cd(2)->SetPad(0.5,0.5,1.0,1.0);
      c1->cd(3)->SetPad(0.0,0.0,0.5,0.5);
      c1->cd(4)->SetPad(0.5,0.0,1.0,0.5);
      c1->cd(1)->SetBottomMargin(0.13);
      c1->cd(2)->SetBottomMargin(0.13);
      c1->cd(3)->SetBottomMargin(0.13);
      c1->cd(4)->SetBottomMargin(0.13);
      c1->cd(1)->SetGridy(1);
      c1->cd(2)->SetGridy(1);
      c1->cd(3)->SetGridy(1);
      c1->cd(4)->SetGridy(1);
      
      TLatex *t1 = new TLatex();
      t1->SetNDC();
      t1->SetTextSize(0.1);
      
      setHistoParsE(hPol0);
      setHistoParsE(hPol45);
      setHistoParsE(hPol90);
      setHistoParsE(hPol135);
      hPol0->GetXaxis()->SetRangeUser(7.6,10.6);
      hPol45->GetXaxis()->SetRangeUser(7.6,10.6);
      hPol90->GetXaxis()->SetRangeUser(7.6,10.6);
      hPol135->GetXaxis()->SetRangeUser(7.6,10.6);

      hPol0->GetXaxis()->SetRangeUser(7.5,10.6);
      hPol45->GetXaxis()->SetRangeUser(7.5,10.6);
      hPol90->GetXaxis()->SetRangeUser(7.5,10.6);
      hPol135->GetXaxis()->SetRangeUser(7.5,10.6);

      c1->cd(1); hPol0->Draw("e1"); t1->DrawLatex(0.75,0.8,"0^{o}");
      c1->cd(2); hPol45->Draw("e1"); t1->DrawLatex(0.75,0.8,"45^{o}");
      c1->cd(3); hPol90->Draw("e1"); t1->DrawLatex(0.75,0.8,"90^{o}");
      c1->cd(4); hPol135->Draw("e1"); t1->DrawLatex(0.75,0.8,"-45^{o}");

      sprintf(fString,"./outFiles/2021-11_TPol.root");
      TFile *fOut = new TFile(fString,"RECREATE");
      hPol0->Write();
      hPol45->Write();
      hPol90->Write();
      hPol135->Write();
      fOut->Close();
    }
  } 
  return tpolToGo;
}

void setHistoParsE(TH1D *histo){
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetYaxis()->SetTitleFont(42);

  histo->SetLineColor(4);
  histo->SetLineWidth(2);  
  histo->SetMarkerStyle(21);  
  histo->SetMarkerColor(4);  
  histo->SetMinimum(-0.05);
  histo->SetMaximum(0.4);

  histo->SetLineWidth(1);  
  histo->SetLineColor(1);  
  histo->SetMarkerColor(1);  
  histo->SetMarkerSize(0.6);  
  
  histo->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  histo->GetYaxis()->SetTitle("P");
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetTitleOffset(0.65);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetRangeUser(7.5,11.4);

}

void setHistoParsEAmo(TH1D *histo){
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetYaxis()->SetTitleFont(42);

  histo->SetLineColor(4);
  histo->SetLineWidth(2);  
  histo->SetMarkerStyle(21);  
  histo->SetMarkerColor(4);  
  //histo->SetMinimum(-0.1);
  //histo->SetMaximum(0.5);

  histo->SetLineWidth(1);  
  histo->SetLineColor(1);  
  histo->SetMarkerColor(1);  
  histo->SetMarkerSize(0.6);  
  
  histo->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  histo->GetYaxis()->SetTitle("P");
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetYaxis()->SetTitleSize(0.07);
  histo->GetYaxis()->SetTitleOffset(0.65);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  //histo->GetXaxis()->SetRangeUser(7.5,11.5);

}

int readVec(char *fToGo,vector<int> *vec)
{
  string line;
  stringstream os(line);
  string temp;
  ifstream myfile (fToGo);
  if (myfile.is_open())
    {
      int j = 0;
      while ( myfile.good() )
        {
          int i = 0;
          getline(myfile,line);
          stringstream os(line);
          while (os >> temp) {
	    vec->push_back(atoi(temp.c_str()));
            i++;
          }
          j++;
        }
      myfile.close();
    }
  else cout << "Unable to open file";
  return 0;
}

int writeRuns(char *fToGo,TH1D *hMask)
{
  ofstream myfile;
  myfile.open (fToGo);
  int vLength = hMask->GetNbinsX();
  for (int j=1; j<=vLength; j++) {
    double maskVal = hMask->GetBinContent(j);
    if (maskVal > 0) {
      int runNumber = (int)hMask->GetBinCenter(j);
      myfile<<runNumber;
      endl(myfile);
    }
  }
  myfile.close();
  return 0;
}

int setMask1d(TH1D *hMask, vector<int> *vec){
  hMask->Reset();
  for (int run = 0; run < vec->size(); run++) {
    int rBin = hMask->FindBin(vec->at(run)*1.0);
    hMask->SetBinContent(rBin,1.0);
    hMask->SetBinError(rBin,0.0);
  }
  return 0;
}
int setMask2d(TH2D *hMask, vector<int> *vec){
  hMask->Reset();
  int nBinsY = hMask->GetNbinsY();
  for (int run = 0; run < vec->size(); run++) {
    int rBin = hMask->GetXaxis()->FindBin(vec->at(run)*1.0);
    for (int yBin = 1; yBin <= nBinsY; yBin++) {
      hMask->SetBinContent(rBin,yBin,1.0);
      hMask->SetBinError(rBin,yBin,0.0);
    }
  }
  return 0;
}

void phi1D(TH2D *hPhiVsRun,TH2D *hMask,TH1D *hPhi,char *hString, TF1 *fitFun, bool qVal){
  TH2D *hPhiVsRunM   = (TH2D*)hPhiVsRun->Clone("hPhiVsRun0M");
  hPhiVsRunM->Multiply(hPhiVsRunM,hMask);
  int nBinsX = hPhiVsRunM->GetNbinsX();
  hPhi   = (TH1D*)hPhiVsRunM->ProjectionY(hString,1,nBinsX);

  int iKillLow = 14;//Actual value
  int iKillHigh = 16;//Actual value


  int iKillLow2 = 8;
  int iKillHigh2 = 8;

  double nPhi = hPhi->Integral();


  if (nPhi == 0){
    return;
  } 
  
  hPhi->Fit("fitFun","QN");

  for (int iKill = iKillLow; iKill<= iKillHigh; iKill++){
    hPhi->SetBinError(iKill,0);
  }
  for (int iKill = iKillLow2; iKill<= iKillHigh2; iKill++){
    hPhi->SetBinError(iKill,0); 
  }

  double maxVal = hPhi->GetMaximum();
  double minVal = hPhi->GetMinimum();
  hPhi->SetMaximum(maxVal*1.1);
  hPhi->SetMinimum(maxVal*0.75);

  hPhi->GetXaxis()->SetTitleFont(62);
  hPhi->GetYaxis()->SetTitleFont(62);
  hPhi->GetXaxis()->SetLabelFont(62);
  hPhi->GetYaxis()->SetLabelFont(62);
  hPhi->GetXaxis()->SetTitle("#phi (degree)");
  hPhi->GetYaxis()->SetTitle("Counts");

  if (qVal == true){
    hPhi->Fit("fitFun","QN");
  }else{
    hPhi->Fit("fitFun","");
  }

  double par0 = fitFun->GetParameter(0);
  double par1 = fitFun->GetParameter(1);
  double par3 = fitFun->GetParameter(3);

  double f0Val = par3/par0;
  //cout<<"f0Val = "<<f0Val<<endl;

  double pTmp = par1/0.1904;

  //cout<<"pTmp = "<<pTmp<<endl;

}
phiPars_t phiSimFit(phiHistos_t phiHistos){

  double rangeLow = 0;
  double rangeHigh = 360;
  
  TF1 * fPhi0   = new TF1("fPhi0","[0]*(1-[1]*cos(2*(3.14/180)*(x-[2]))+[3]*cos(3.14*(x-[4])/180))",rangeLow,rangeHigh);
  TF1 * fPhi45  = new TF1("fPhi45","[0]*(1-[1]*cos(2*(3.14/180)*(x-[2]-45))+[3]*cos(3.14*(x-[4])/180))",rangeLow,rangeHigh);
  TF1 * fPhi90  = new TF1("fPhi90","[0]*(1-[1]*cos(2*(3.14/180)*(x-[2]-90))+[3]*cos(3.14*(x-[4])/180))",rangeLow,rangeHigh);
  TF1 * fPhi135 = new TF1("fPhi135","[0]*(1-[1]*cos(2*(3.14/180)*(x-[2]-135))+[3]*cos(3.14*(x-[4])/180))",rangeLow,rangeHigh);
  
  ROOT::Math::WrappedMultiTF1 wfPhi0(*fPhi0,1);
  ROOT::Math::WrappedMultiTF1 wfPhi45(*fPhi45,1);
  ROOT::Math::WrappedMultiTF1 wfPhi90(*fPhi90,1);
  ROOT::Math::WrappedMultiTF1 wfPhi135(*fPhi135,1);
  
  ROOT::Fit::DataOptions opt;
  
  ROOT::Fit::DataRange rangePhi0;
  rangePhi0.SetRange(rangeLow,rangeHigh);
  ROOT::Fit::BinData dataPhi0(opt,rangePhi0);
  ROOT::Fit::FillData(dataPhi0, phiHistos.hPhi0);

  ROOT::Fit::DataRange rangePhi45;
  rangePhi45.SetRange(rangeLow,rangeHigh);
  ROOT::Fit::BinData dataPhi45(opt,rangePhi45);
  ROOT::Fit::FillData(dataPhi45, phiHistos.hPhi45);
  
  ROOT::Fit::DataRange rangePhi90;
  rangePhi90.SetRange(rangeLow,rangeHigh);
  ROOT::Fit::BinData dataPhi90(opt,rangePhi90);
  ROOT::Fit::FillData(dataPhi90, phiHistos.hPhi90);
  
  ROOT::Fit::DataRange rangePhi135;
  rangePhi135.SetRange(rangeLow,rangeHigh);
  ROOT::Fit::BinData dataPhi135(opt,rangePhi135);
  ROOT::Fit::FillData(dataPhi135, phiHistos.hPhi135);
  
  ROOT::Fit::Chi2Function chi2_Phi0(dataPhi0, wfPhi0);
  ROOT::Fit::Chi2Function chi2_Phi45(dataPhi45, wfPhi45);
  ROOT::Fit::Chi2Function chi2_Phi90(dataPhi90, wfPhi90);
  ROOT::Fit::Chi2Function chi2_Phi135(dataPhi135, wfPhi135);
  
  GlobalChi2 globalChi2(chi2_Phi0, chi2_Phi45, chi2_Phi90, chi2_Phi135);
  
  ROOT::Fit::Fitter fitter;
  
  // create before the parameter settings in order to fix or set range on them          
  fitter.Config().SetParamsSettings(Npar,simPars);
  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","Migrad");
  fitter.FitFCN(Npar,globalChi2,0,dataPhi0.Size()+dataPhi45.Size()+dataPhi90.Size()+dataPhi135.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);
  
  fPhi0->SetFitResult( result, iparPhi0);
  fPhi45->SetFitResult( result, iparPhi45);
  fPhi90->SetFitResult( result, iparPhi90);
  fPhi135->SetFitResult( result, iparPhi135);


  double pol0   = fPhi0->GetParameter(1)/0.1904;
  double err0   = fPhi0->GetParError(1)/0.1904;
  double pol45  = fPhi45->GetParameter(1)/0.1904;
  double err45  = fPhi45->GetParError(1)/0.1904;
  double pol90  = fPhi90->GetParameter(1)/0.1904;
  double err90  = fPhi90->GetParError(1)/0.1904;
  double pol135 = fPhi135->GetParameter(1)/0.1904;
  double err135 = fPhi135->GetParError(1)/0.1904;


  cout<<"simPars = "<<*simPars<<endl;
  cout<<"npar = "<<fPhi0->GetNumberFreeParameters()<<endl;
  cout<<"tmp0   = "<<fPhi0->GetParError(0)<<endl;
  cout<<"tmp1   = "<<fPhi0->GetParError(1)<<endl;
  cout<<"tmp2   = "<<fPhi0->GetParError(2)<<endl;
  cout<<"tmp3   = "<<fPhi0->GetParError(3)<<endl;
  cout<<"tmp4   = "<<fPhi0->GetParError(4)<<endl;


  cout<<"pol0   = "<<pol0<<" +/- "<<err0<<endl;
  cout<<"pol45  = "<<pol45<<" +/- "<<err45<<endl;
  cout<<"pol90  = "<<pol90<<" +/- "<<err90<<endl;
  cout<<"pol135 = "<<pol135<<" +/- "<<err135<<endl;
  
  double sigRho0 = 0.3301/pol0;
  double sigRho45 = 0.3217/pol45;
  double sigRho90 = 0.322/pol90;
  double sigRho135 = 0.3367/pol135;
  double errRho0 = sigRho0*(err0/pol0);
  double errRho45 = sigRho45*(err45/pol45);
  double errRho90 = sigRho90*(err90/pol90);
  double errRho135 = sigRho135*(err135/pol135);
  cout<<"SigRho0 = "<<sigRho0<<" +/- "<<errRho0<<endl;
  cout<<"SigRho45 = "<<sigRho45<<" +/- "<<errRho45<<endl;
  cout<<"SigRho90 = "<<sigRho90<<" +/- "<<errRho90<<endl;
  cout<<"SigRho135 = "<<sigRho135<<" +/- "<<errRho135<<endl;

  phiPars_t phiPars;
  phiPars.pol0 = pol0; 
  phiPars.pol45 = pol45; 
  phiPars.pol90 = pol90; 
  phiPars.pol135 = pol135; 

  phiPars.err0 = err0; 
  phiPars.err45 = err45; 
  phiPars.err90 = err90; 
  phiPars.err135 = err135; 
  
  phiPars.offset = fPhi0->GetParameter(2);
  phiPars.aCenter = fPhi0->GetParameter(3);
  phiPars.aOffset = fPhi0->GetParameter(4);


  return phiPars;
}
