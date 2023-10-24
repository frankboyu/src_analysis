#include <TLatex.h>

void plot(){
  TString infile  = "results/test2.root";
  TString outfile = "results/test2.pdf";

  TChain *EvTree = new TChain("EvTree");
  EvTree->Add(infile);

  TCanvas *c1 = new TCanvas("c1","c1", 1920,1080);

  //kinematic fit
  EvTree->Draw("CL>>h_CL(1000,0,1)","","");
  TH1F *h_CL = (TH1F*)gPad->GetPrimitive("h_CL");
  h_CL->SetTitle("confidence level;CL;counts");
  c1->SetLogy(1);
  c1->SaveAs(outfile+"(");

  EvTree->Draw("mass_phi>>h_mass_phi(400,0.8,1.2)","","ep");
  TH1F *h_mass_phi = (TH1F*)gPad->GetPrimitive("h_mass_phi");
  h_mass_phi->SetTitle("invariant mass of reconstructed #Phi meson;m(K^{+}K^{-})[GeV];counts/MeV");
  c1->SetLogy(0);
  c1->SaveAs(outfile);

  EvTree->Draw("Y_vtx:X_vtx>>h_vertex_XY(120,-6.0,6.0,120,-6.0,6.0)","","colz");
  TH2F *h_vertex_XY = (TH2F*)gPad->GetPrimitive("h_vertex_XY");
  h_vertex_XY->SetTitle("common vertex X and Y of kinematic fitting;vertex x[cm];vertex y[cm]");
  c1->SaveAs(outfile);

  EvTree->Draw("Z_vtx>>h_vertex_Z(120,40,100)","","");
  TH1F *h_vertex_Z = (TH1F*)gPad->GetPrimitive("h_vertex_Z");
  h_vertex_Z->SetTitle("common vertex Z of kinematic fitting;vertex z[cm];counts");
  c1->SaveAs(outfile);

  EvTree->Draw("theta_km*180/3.1415926:Pmag_km>>h_kminus(100,0,10,180,0,180)","","colz");
  TH2F *h_kminus = (TH2F*)gPad->GetPrimitive("h_kminus");
  h_kminus->SetTitle("#theta vs P distribution of K^{-};P[GeV];#theta[deg]");
  c1->SaveAs(outfile);

  EvTree->Draw("theta_kp*180/3.1415926:Pmag_kp>>h_kplus(100,0,10,180,0,180)","","colz");
  TH2F *h_kplus = (TH2F*)gPad->GetPrimitive("h_kplus");
  h_kplus->SetTitle("#theta vs P distribution of K^{+};P[GeV];#theta[deg]");
  c1->SaveAs(outfile);
  
  EvTree->Draw("theta_p*180/3.1415926:Pmag_p>>h_proton(100,0,10,180,0,180)","","colz");
  TH2F *h_proton = (TH2F*)gPad->GetPrimitive("h_proton");
  h_proton->SetTitle("#theta vs P distribution of proton;P[GeV];#theta[deg]");
  c1->SaveAs(outfile);

  EvTree->Draw("mass_kmp:mass_phi>>h_dalitz(120,0.8,2,400,1,5)","","colz");
  TH2F *h_dalitz = (TH2F*)gPad->GetPrimitive("h_dalitz");
  h_dalitz->SetTitle("Dalitz-like plot;m_{K^{-}K^{+}};m_{K^{-}p}");
  c1->SaveAs(outfile);

  EvTree->Draw("mass_rho>>h_mass_rho(2000,0,2)","","ep");
  TH1F *h_mass_rho = (TH1F*)gPad->GetPrimitive("h_mass_rho");
  h_mass_rho->SetTitle("invariant mass of reconstructed #rho meson assuming kaons as pions;m_{#rho};counts");
  c1->SaveAs(outfile);

  EvTree->Draw("mass_rho:mass_phi>>h_mass_2D(700,0.9,1.6,1000,0.2,1.2)","","colz");
  TH2F *h_mass_2D = (TH2F*)gPad->GetPrimitive("h_mass_2D");
  h_mass_2D->SetTitle("invariant mass of reconstructed #rho meson versus #Phi meson;m_{#Phi};m_{#rho}");
  c1->SaveAs(outfile);

  //photon info
  EvTree->Draw("dt>>h_dt(72,-18,18)","","");
  TH1F *h_dt = (TH1F*)gPad->GetPrimitive("h_dt");
  h_dt->SetTitle("#Delta t of beam photons;#Delta t;counts");
  c1->SaveAs(outfile);

  EvTree->Draw("dE>>h_dE_offtime(20,-1,1)","weight<0","");
  TH1F *h_dE_offtime = (TH1F*)gPad->GetPrimitive("h_dE_offtime");
  EvTree->Draw("dE>>h_dE_ontime(20,-1,1)","weight>0","");
  TH1F *h_dE_ontime = (TH1F*)gPad->GetPrimitive("h_dE_ontime");
  TH1F h_dE = *h_dE_ontime - (1.0/6.0)*(*h_dE_offtime);
  h_dE.SetTitle("energy balance;dE[GeV];counts");
  h_dE.Draw();
  c1->SaveAs(outfile);
  
  EvTree->Draw("t>>h_t_offtime(2500,-5,20)","weight<0","");
  TH1F *h_t_offtime = (TH1F*)gPad->GetPrimitive("h_t_offtime");
  EvTree->Draw("t>>h_t_ontime(2500,-5,20)","weight>0","");
  TH1F *h_t_ontime = (TH1F*)gPad->GetPrimitive("h_t_ontime");
  TH1F h_t = *h_t_ontime - (1.0/6.0)*(*h_t_offtime);
  h_t.SetTitle("-t;-t[GeV];counts");
  h_t.Draw();
  c1->SetLogy(1);
  c1->SaveAs(outfile);

  EvTree->Draw("u>>h_u_offtime(2500,-5,20)","weight<0","");
  TH1F *h_u_offtime = (TH1F*)gPad->GetPrimitive("h_u_offtime");
  EvTree->Draw("u>>h_u_ontime(2500,-5,20)","weight>0","");
  TH1F *h_u_ontime = (TH1F*)gPad->GetPrimitive("h_u_ontime");
  TH1F h_u = *h_u_ontime - (1.0/6.0)*(*h_u_offtime);
  h_u.SetTitle("-u;-u[GeV];counts");
  h_u.Draw();
  c1->SetLogy(0);
  c1->SaveAs(outfile);

  EvTree->Draw("kmiss>>h_kmiss_offtime(20,0,2)","weight<0","");
  TH1F *h_kmiss_offtime = (TH1F*)gPad->GetPrimitive("h_kmiss_offtime");
  EvTree->Draw("kmiss>>h_kmiss_ontime(20,0,2)","weight>0","");
  TH1F *h_kmiss_ontime = (TH1F*)gPad->GetPrimitive("h_kmiss_ontime");
  TH1F h_kmiss = *h_kmiss_ontime - (1.0/6.0)*(*h_kmiss_offtime);
  h_kmiss.SetTitle("k_{miss};k_{miss}[GeV];counts");
  h_kmiss.Draw();
  c1->SaveAs(outfile+")");
}
  //double weight = 6.59575e+00  //MF
  //double weight  = 2.03753e+01  //SRC
  /*
  TH1F *masskmp   = new TH1F( "masskmp", "",  bin, 0., 6);
  TH1F *masskpp   = new TH1F( "masskpp", "", bin, 0., 6);
 
  TH1F *protonpt  = new TH1F( "protonpt", "", bin, 0., 2.); 
  TH1F *kpluspt  = new TH1F( "kpluspt", "", bin, 0., 2.);
  TH1F *kminuspt = new TH1F( "kminuspt", "", bin, 0., 2.);
  TH1F *protonpz  = new TH1F( "protonpz", "", bin, 0., 10.);
  TH1F *kpluspz  = new TH1F( "kpluspz", "", bin, 0., 10.);
  TH1F *kminuspz = new TH1F( "kminuspz", "", bin, 0., 10.);

  TH1F *deltaphi   = new TH1F( "deltaphi", "", bin, 120, 230);
  TH1F *deltaphikmp     = new TH1F( "deltaphikmp", "", bin, 120, 230);
  TH1F *deltaphikpp     = new TH1F( "deltaphikpp", "", bin,  0, 360);

  TH1F *hpmissplus   = new TH1F( "hpmissplus", "", bin, -1, 3.5);
  TH1F *hpmissplusbgd   = new TH1F( "hpmissplusbgd", "", bin, -1, 3.5);
  TH1F *hs    = new TH1F( "hs",  "", bin, 0, 30.);
TH2F *massphi2D  = new TH2F( "massphi2D", "", 2*bin, 0., 5, 2*bin, 0., 20.);

TH1F *hkmiss   = new TH1F( "hkmiss", "", bin/2, 0, 1.);
TH2F *thetaCMt     = new TH2F("thetaCMt",    "", bin, -1, 1, bin, 0, 10);
TH1F *hemiss      = new TH1F( "hemiss", "", bin, -1.2, 5);
TH1F *hpmissmin   = new TH1F( "hpmissmin", "", bin, 0, 1.5);
TH1F *hpmissperp   = new TH1F( "hpmissperp", "", bin, 0, 1.5);

TH1F *hmmiss   = new TH1F( "hmmiss", "", bin, 0, 1.5);
TH1F *hminv2    = new TH1F( "hminv2", "", bin, 0, 2.);

TH1F *hpmiss   = new TH1F( "hpmiss", "", bin, 0, 2.);
TH2F *massphi2D2  = new TH2F( "massrh02D2", "", 2*bin, 0.,1, 2*bin, 0.5, 1.5);
TH2F *XY = new TH2F("XY", "", 100, xmin, xmax, 100, xmin, xmax);
TH2F *mass2D = new TH2F("mass2D", "", bins1, 0., 4., bins1, 0., 20.);
TH2F *massomega = new TH2F("massomega", "", bins1, 0.5, 2., bins1, -2, 2);
TH2F *massphiomega = new TH2F("massphiomega", "", bins1, 0., 20., bins1, -2, 2);
TH2F *massphiomega2 = new TH2F("massphiomega2", "", bins1, 0, 25., bins1, -2, 2);
TH2F *pminu    = new TH2F("pminu", "", bins1, 0, 30,bins1, 0,1.5);
	//if (mass_kmp*mass_kmp<mass_kmpcut) continue;
	//if (mass_kpp*mass_kpp<mass_kmpcut) continue;
	//if( theta_phi > 0.2) continue;  
	//if(mass_phi < 1.) continue; 
		
		if (massinv2[j]<0 ) continue;
		//if (t[j]<tcut) continue;
		//if (t[j]<1) continue;

		//if (u[j]<1) continue;
		//if (t[j]>10) continue;
		
		//if (omega[j] >  (mass_kmp*mass_kmp/10. - 0.3)) continue;
		//if (omega[j] <0) continue;
		//if (kmiss[j] > kmisscut )  continue;
		//if (Pmissperp[j] > pmissmincut) continue;
		//if (Mmiss[j]<0. ||  Mmiss[j]>1.5) continue;

TLine *lzmax=new TLine(zmax,0,zmax,Zvertex->GetMaximum());
lzmax->SetLineColor(kBlue);lzmax->SetLineWidth(3);
massphi1Dcut->SetFillColor(kRed-2);
protonpt->SetMinimum(1);
protonpt->Draw("HIST");
auto legend = new TLegend(0.5,0.6,0.85,0.85);
legend->AddEntry(hkmiss, "-t > 1 GeV", "epl");
legend->Draw();
gPad->SetLogy();
  */
