#import "models.h"
#include <iostream>
#include <stdio.h>
#include <string.h>

void filter_classic(string inFileName, string outFileName, string outHistName, string treeName = "gd_pi0pimprotprotinc__B4_F4_T1_S2")
{

  using namespace ROOT;
  using namespace std;
  EnableImplicitMT();

  TChain chain(treeName.c_str());
  chain.Add(inFileName.c_str());

  //RDataFrame rdf_raw("rhoMinus_p_p","tree_DSelector_deuterium_data.root");
  RDataFrame rdf_raw(chain);
  cout << "a\n";
  auto rdf_def = rdf_raw
    .Define("CL","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("N_p4","TLorentzVector(0,0,0,0.93892)")
    .Define("p1isLead","p1_p4_kin.T()>p2_p4_kin.T()")
    .Define("p_p4_meas","p1isLead ? p1_p4_meas : p2_p4_meas")
    .Define("p_p4_kin","p1isLead ? p1_p4_kin : p2_p4_kin")
    .Define("rec_p4_meas","!p1isLead ? p1_p4_meas : p2_p4_meas")
    .Define("rec_p4_kin","!p1isLead ? p1_p4_kin : p2_p4_kin")
    .Define("pi0_p4_meas","g1_p4_meas + g2_p4_meas")
    .Define("pi0_p4_kin","g1_p4_kin + g2_p4_kin")
    .Define("rho_p4_meas","pi0_p4_meas + pim_p4_meas")
    .Define("rho_p4_kin","pi0_p4_kin + pim_p4_kin")
    .Define("lead_p4_meas","rho_p4_meas + p_p4_meas")
    .Define("lead_p4_kin","rho_p4_kin + p_p4_kin")
    .Define("struck_p4_meas","lead_p4_meas - beam_p4_meas")
    .Define("struck_p4_kin","lead_p4_kin - beam_p4_kin")
    .Define("CM_p4_meas","struck_p4_meas + rec_p4_meas")
    .Define("CM_p4_kin","struck_p4_kin + rec_p4_kin")
    .Define("total_p4_meas","lead_p4_meas + rec_p4_meas")
    .Define("total_p4_kin","lead_p4_kin + rec_p4_kin")
    .Define("xVertex","pim_x4_kin.X()")
    .Define("yVertex","pim_x4_kin.Y()")
    .Define("zVertex","pim_x4_kin.Z()")
    .Define("pip1_p4_meas","TLorentzVector(p1_p4_meas.Vect(),sqrt(p1_p4_meas.Vect().Mag2() + 0.139570*0.139570))")
    .Define("pip2_p4_meas","TLorentzVector(p2_p4_meas.Vect(),sqrt(p2_p4_meas.Vect().Mag2() + 0.139570*0.139570))")
    .Define("omega1_p4_meas","pip1_p4_meas + pim_p4_meas + pi0_p4_meas")
    .Define("omega2_p4_meas","pip2_p4_meas + pim_p4_meas + pi0_p4_meas")
    .Define("omega1_m_meas","omega1_p4_meas.M()")
    .Define("omega2_m_meas","omega2_p4_meas.M()")
    .Define("pip1_p4_kin","TLorentzVector(p1_p4_kin.Vect(),sqrt(p1_p4_kin.Vect().Mag2() + 0.139570*0.139570))")
    .Define("pip2_p4_kin","TLorentzVector(p2_p4_kin.Vect(),sqrt(p2_p4_kin.Vect().Mag2() + 0.139570*0.139570))")
    .Define("omega1_p4_kin","pip1_p4_kin + pim_p4_kin + pi0_p4_kin")
    .Define("omega2_p4_kin","pip2_p4_kin + pim_p4_kin + pi0_p4_kin")
    .Define("omega1_m_kin","omega1_p4_kin.M()")
    .Define("omega2_m_kin","omega2_p4_kin.M()")
    .Define("pLead","p_p4_kin.P()")
    .Define("pRec","rec_p4_kin.P()")
    .Define("alphaCM_meas","CM_p4_meas.Minus()/N_p4.M()")
    .Define("proxy_plus_meas","2*N_p4.M() - (N_p4.M2() + struck_p4_meas.Perp2())/(2*N_p4.M() - struck_p4_meas.Minus())")
    .Define("proxy_E_meas","0.5*(proxy_plus_meas + struck_p4_meas.Minus())")
    .Define("proxy_pz_meas","0.5*(proxy_plus_meas - struck_p4_meas.Minus())")
    .Define("proxy_p4_meas","TLorentzVector(struck_p4_meas.X(),struck_p4_meas.Y(),proxy_pz_meas,proxy_E_meas)")
    .Define("alphaCM_kin","CM_p4_kin.Minus()/N_p4.M()")
    .Define("proxy_plus_kin","2*N_p4.M() - (N_p4.M2() + struck_p4_kin.Perp2())/(2*N_p4.M() - struck_p4_kin.Minus())")
    .Define("proxy_E_kin","0.5*(proxy_plus_kin + struck_p4_kin.Minus())")
    .Define("proxy_pz_kin","0.5*(proxy_plus_kin - struck_p4_kin.Minus())")
    .Define("proxy_p4_kin","TLorentzVector(struck_p4_kin.X(),struck_p4_kin.Y(),proxy_pz_kin,proxy_E_kin)")
    .Define("cosGamma_proxy","cos(rec_p4_kin.Angle(proxy_p4_kin.Vect()))")
    .Define("thetaP","p_p4_kin.Theta()*180./M_PI")
    .Define("thetaRec","rec_p4_kin.Theta()*180./M_PI")
    .Define("eLead","lead_p4_kin.T()")
    .Define("m2pi_meas","rho_p4_meas.M()")
    .Define("m2pi_kin","rho_p4_kin.M()")
    .Define("s","lead_p4_kin.M2()")
    .Define("t","(beam_p4_kin - rho_p4_kin).M2()")
    .Define("u","(beam_p4_kin - p_p4_kin).M2()")
    .Define("kmiss","N_p4.M()*sqrt((N_p4.M2() + struck_p4_kin.Perp2())/(struck_p4_kin.Minus()*(2*N_p4.M() - struck_p4_kin.Minus())) - 1)")
    .Define("Emiss2N","2*N_p4.M() - CM_p4_kin.T()");
  cout << "b\n";
  auto rdf_no_filter = rdf_def;//.Filter([](){ return true;},{});
  //auto rdf_vertex_filtered = rdf_no_filter.Filter("(zVertex > 51) && (zVertex < 79) && (xVertex*xVertex + yVertex*yVertex < 1*1)");
  auto rdf_vertex_filtered = rdf_no_filter.Filter([](double xVertex, double yVertex, double zVertex) {return (zVertex > 51) && (zVertex < 79) && (xVertex*xVertex + yVertex*yVertex < 1*1);},{"xVertex","yVertex","zVertex"});
  //auto rdf_energy_filtered = rdf_vertex_filtered.Filter("eLead > 7"); // Previously 7 GeV, now removing it as superfluous
  auto rdf_energy_filtered = rdf_vertex_filtered.Filter([](double eLead) {return eLead > 7;}, {"eLead"});
  //auto rdf_omega_filtered = rdf_energy_filtered.Filter("(omega1_m_kin > 0.9) && (omega2_m_kin > 0.9)");
  auto rdf_omega_filtered = rdf_energy_filtered.Filter([](double omega1_m_kin, double omega2_m_kin) {return (omega1_m_kin > 0.9) && (omega2_m_kin > 0.9);},{"omega1_m_kin","omega2_m_kin"});
  //auto rdf_bg_filtered = rdf_omega_filtered.Filter("alphaCM_kin/2.4 + thetaP/90. > 1");
  auto rdf_bg_filtered = rdf_omega_filtered.Filter([](double alphaCM_kin, double thetaP) {return alphaCM_kin/2.4 + thetaP/90. > 1.;},{"alphaCM_kin","thetaP"});
  //auto rdf_bg_filtered = rdf_omega_filtered.Filter("cosGamma_proxy < 4./3.*alphaCM_kin*N_p4.M() - 3.");
  //auto rdf_kmiss_filtered = rdf_bg_filtered.Filter("(kmiss > 0.4) && (kmiss < 1.2)");
  auto rdf_kmiss_filtered = rdf_bg_filtered.Filter([](double kmiss) {return (kmiss > 0.4) && (kmiss < 1.2);},{"kmiss"});
  //auto rdf_rec_filtered = rdf_kmiss_filtered.Filter("(pRec > 0.3) && (pRec < 1.0)");
  auto rdf_rec_filtered = rdf_kmiss_filtered.Filter([](double pRec) {return pRec > 0.3 && (pRec < 1.0);}, {"pRec"});
  //auto rdf_Emiss_filtered = rdf_rec_filtered.Filter("(fabs(Emiss2N) < 1)*(beam_p4_kin.T()>6.)");
  auto rdf_Emiss_filtered = rdf_rec_filtered.Filter([](double Emiss2N, TLorentzVector beam_p4_kin) {return (fabs(Emiss2N) < 1) && (beam_p4_kin.T()>6.);},{"Emiss2N","beam_p4_kin"});
  cout << "c\n";
  //auto rdf_with_diff = rdf_energy_filtered.Filter("(kmiss > 0.4) && (kmiss < 1.2) && (pRec > 0.3) && (pRec < 1.0) && (fabs(Emiss2N) < 1) && (beam_p4_kin.T()>6.) && (m2pi_kin > 0.65) && (m2pi_kin < 0.9) && (t < -1) && (u < -1)");

  // Take snapshot
  ColumnNames_t columns = {"kmiss","s","t","u"};
  rdf_Emiss_filtered.Snapshot("rhoMinus_p_p_filtered",outFileName);//,columns);
  cout << "d\n";
  int N_filters = 8;
  RNode rdfs [] = {rdf_no_filter , rdf_vertex_filtered , rdf_energy_filtered , rdf_omega_filtered , rdf_bg_filtered , rdf_kmiss_filtered , rdf_rec_filtered , rdf_Emiss_filtered};//, rdf_with_diff};
  string labels [] = {"no_cut" , "vertex_cut" , "energy_cut" , "omega_cut" , "diffractive_cut" , "kmiss_cut" , "recoil_cut" , "Emiss_cut"};// , "with_diff"};

  // Fill and write histograms
  TFile * histFile = new TFile(outHistName.c_str(),"RECREATE");

  histFile->cd();

  //vector<TH1*> hists;

  for (int i = 0; i < N_filters; i++) {

    cout << i << "\n";

    auto rdf = rdfs[i];
    string label = labels[i];

    TH2D h_classic = *rdf.Histo2D<double,double>(alphaCM_cosGammaProxy_model,"alphaCM_kin","cosGamma_proxy","accidweight");
    h_classic.SetName(("classic_" + label).c_str());
    //hists.push_back(&h_classic);
    //h_classic.Write();
    cout << "e\n";
    TH2D h_angular = *rdf.Histo2D<double,double>(alphaCM_thetaP_model,"alphaCM_kin","thetaP","accidweight");
    h_angular.SetName(("angular_" + label).c_str());
    //hists.push_back(&h_angular);
    //h_angular.Write();
    cout << "f\n";
    TH1D h_mass = *rdf.Histo1D<double>(mass_model,"m2pi_kin","accidweight");
    h_mass.SetName(("mass_" + label).c_str());
    //hists.push_back(&h_mass);
    //h_mass.Write();
    cout << "g\n";

    h_classic.Write();
    h_angular.Write();
    h_mass.Write();

  }

  // Before cuts
  TH1D h_z = *rdf_no_filter.Histo1D(zVertex_model,"zVertex","accidweight");
  hists.push_back(&h_z);
  //h_z.Write();
  TH2D h_xy = *rdf_no_filter.Histo2D(xyVertex_model,"xVertex","yVertex","accidweight");
  hists.push_back(&h_xy);
  //h_xy.Write();

  // After vertex cut
  TH1D h_energy = *rdf_vertex_filtered.Histo1D(E_model,"eLead","accidweight");
  hists.push_back(&h_energy);
  //h_energy.Write();

  // After energy cut
  TH1D h_omega1 = *rdf_energy_filtered.Histo1D(omega_model,"omega1_m_meas","accidweight");
  h_omega1.SetName("omega1_mass");
  hists.push_back(&h_omega1);
  //h_omega1.Write();
  TH1D h_omega2 = *rdf_energy_filtered.Histo1D(omega_model,"omega2_m_meas","accidweight");
  h_omega2.SetName("omega2_mass");
  hists.push_back(&h_omega2);
  //h_omega2.Write();

  // After diffractive background cut
  TH2D h_proton_momenta = *rdf_bg_filtered.Histo2D(proton_momenta_model,"pLead","pRec","accidweight");
  hists.push_back(&h_proton_momenta);
  //h_proton_momenta.Write();

  TH2D h_kmiss_rec = *rdf_bg_filtered.Histo2D(kmiss_rec_model,"kmiss","pRec","accidweight");
  hists.push_back(&h_kmiss_rec);
  //h_kmiss_rec.Write();

  // After momentum cuts
  TH1D h_Emiss = *rdf_rec_filtered.Histo1D(Emiss_model,"Emiss2N","accidweight");
  hists.push_back(&h_Emiss);
  //h_Emiss.Write();

  for (TH1* hist : hists)
    {
      hist->Write();
    }

  histFile->Close();


}
