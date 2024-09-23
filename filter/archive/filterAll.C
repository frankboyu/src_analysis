#include "analyzeData.h"
double returnNonZeroPoverE(double em_ebcal, double ep_ebcal, double em_p, double ep_p){
  if(em_ebcal>0){
    return em_p/em_ebcal;}
  else if(ep_ebcal>0){
    return ep_p/ep_ebcal;}
  else {
    return 0;}
}

double returnGreaterThetaVertex(double em_theta, double ep_theta, double vertex_em, double vertex_ep){
  if(em_theta>ep_theta){
    return vertex_em;}
  else {
    return vertex_ep;}
}

bool passEoverP(double p, double Efcal, double Ebcal, bool isFCAL, double sigma_minus, double sigma_plus)
{
  double fcal_mu=1.06319717;
  double fcal_sigma=0.085128472;
  double bcal_mu=1.028670304;
  double bcal_sigma=0.072119971;

  if (isFCAL)
  {
    return ((p/Efcal > fcal_mu + sigma_minus * fcal_sigma) && (p/Efcal < fcal_mu + sigma_plus * fcal_sigma));
  }
  else
    return ((p/Ebcal > bcal_mu + sigma_minus * bcal_sigma) && (p/Ebcal < bcal_mu + sigma_plus * bcal_sigma));
}

void filterAll(string inFileName, string outFileName, string treeName = "e_e_p_X", string track_shower_sel="1",
               string preBCut="0.030", string sigma_minus="-3", string sigma_plus="2", string EMissCut="1")
{

  using namespace ROOT;
  using namespace std;

  EnableImplicitMT();
  TChain chain(treeName.c_str());
  chain.Add(inFileName.c_str());

  RDataFrame rdf_raw(chain);

  auto rdf_def = rdf_raw
          .Define("mN","0.938")
          .Define("N2_p4","TLorentzVector(0,0,0,2*mN)")
          .Define("jpsi_p4_meas","ep_p4_meas + em_p4_meas")
          .Define("jpsi_p4_kin","ep_p4_kin + em_p4_kin")
          .Define("jpsi_m_meas","jpsi_p4_meas.M()")
          .Define("jpsi_m_kin","jpsi_p4_kin.M()")
          .Define("init_p4_kin","jpsi_p4_kin + p_p4_kin - beam_p4_kin")
          .Define("pMinus_init","init_p4_kin.Minus()")
          .Define("pPerp_init","init_p4_kin.Perp()")
          .Define("miss_p4_meas","N2_p4 + beam_p4_meas - jpsi_p4_meas - p_p4_meas")
          .Define("miss_p4_kin","N2_p4 + beam_p4_kin - jpsi_p4_kin - p_p4_kin")
          .Define("miss_m_meas","miss_p4_meas.M()")
          .Define("miss_m_kin","miss_p4_kin.M()")
          .Define("ep_EoverP","(ep_efcal+ep_ebcal)/ep_p4_kin.P()")
          .Define("em_EoverP","(em_efcal+em_ebcal)/em_p4_kin.P()")
          .Define("xVertex","p_x4_kin.X()")
          .Define("yVertex","p_x4_kin.Y()")
          .Define("zVertex","p_x4_kin.Z()")
          .Define("Delta_kin","beam_p4_kin - jpsi_p4_kin")
          .Define("t_kin","Delta_kin.M2()")
          .Define("minust_kin","-t_kin")
          .Define("pmiss","miss_p4_kin.P()")
          .Define("k","beam_p4_kin.P()")
          .Define("Delta_plus","1*(miss_p4_kin.Plus() - (mN*mN + miss_p4_kin.Perp2())/miss_p4_kin.Minus())")
          .Define("p_em","em_p4_kin.P()")
          .Define("em_theta","em_p4_kin.Theta()*180/3.14159265")
          .Define("p_ep","ep_p4_kin.P()")
          .Define("ep_theta","ep_p4_kin.Theta()*180/3.14159265")
          .Define("p_p","p_p4_kin.P()")
          .Define("p_theta","p_p4_kin.Theta()*180/3.14159265")
          .Define("jpsi_p","jpsi_p4_kin.P()")
          .Define("jpsi_theta","jpsi_p4_kin.Theta()*180/3.14159265")
          .Define("jpsi_m_pair","sqrt(jpsi_m_kin*jpsi_m_kin + 1*Delta_plus * jpsi_p4_kin.Minus())")
          .Define("jpsi_m_stand","sqrt((mN + beam_p4_kin.Plus() - p_p4_kin.Plus()) * jpsi_p4_kin.Minus() - jpsi_p4_kin.Perp2())")
          .Define("E_gamma","beam_p4_kin.T()")
          .Define("em_ep_angle", "(ep_p4_kin.Angle(em_p4_kin.Vect())*180/3.1415)")
          .Define("em_eprebcal_sinTheta", "em_eprebcal*sin(em_theta*3.1415/180)")
          .Define("ep_eprebcal_sinTheta", "ep_eprebcal*sin(ep_theta*3.1415/180)")
          .Define("em_pMinus","em_p4_kin.Minus()")
          .Define("ep_pMinus","ep_p4_kin.Minus()")
          .Define("p_pMinus","p_p4_kin.Minus()")
          .Define("jpsi_pMinus","jpsi_p4_kin.Minus()")
          .Define("em_pPlus","em_p4_kin.Plus()")
          .Define("ep_pPlus","ep_p4_kin.Plus()")
          .Define("p_pPlus","p_p4_kin.Plus()")
          .Define("jpsi_pPlus","jpsi_p4_kin.Plus()")
          .Define("em_perp","em_p4_kin.Perp()")
          .Define("ep_perp","ep_p4_kin.Perp()")
          .Define("p_perp","p_p4_kin.Perp()")
          .Define("jpsi_perp","jpsi_p4_kin.Perp()")
          .Define("jpsi_alpha","(jpsi_pMinus/mN)")
          .Define("p_alpha","(p_pMinus/mN)")
          .Define("n_alpha","(2-jpsi_alpha-p_alpha)")
          .Define("E_miss_stat_p","(mN + beam_p4_kin.T() - jpsi_p4_kin.T() - p_p4_kin.T())")
          .Define("E_prot_init","(jpsi_p4_kin.T() + p_p4_kin.T() - beam_p4_kin.T())")
          .Define("em_PoverE","em_p4_kin.P()/(em_efcal+em_ebcal)")
          .Define("ep_PoverE","ep_p4_kin.P()/(ep_efcal+ep_ebcal)")
          .Define("tot_p4_kin","(jpsi_p4_kin+p_p4_kin-beam_p4_kin)")
          .Define("tot_perp","(tot_p4_kin.Perp())")
          .Define("tot_pMinus","(tot_p4_kin.Minus())")
          .Define("alpha_miss","(tot_pMinus/mN)")
          .Define("kmiss","mN*sqrt((mN*mN + tot_p4_kin.Perp2())/(tot_p4_kin.Minus()*(2*mN - tot_p4_kin.Minus())) - 1)")
          .Define("mJPsi","(3.096)")
          .Define("t_kin_minus","(-(mJPsi*mJPsi-2*E_gamma*(jpsi_p4_kin.Minus())))")
          .Define("xVertexLepton","returnGreaterThetaVertex(em_theta, ep_theta, em_x4_meas.X(), ep_x4_meas.X())")
          .Define("yVertexLepton","returnGreaterThetaVertex(em_theta, ep_theta, em_x4_meas.Y(), ep_x4_meas.Y())")
          .Define("zVertexLepton","returnGreaterThetaVertex(em_theta, ep_theta, em_x4_meas.Z(), ep_x4_meas.Z())")
          .Define("xVertexProton","p_x4_meas.X()")
          .Define("yVertexProton","p_x4_meas.Y()")
          .Define("zVertexProton","p_x4_meas.Z()")
          .Define("deltaXYVertexLepton","sqrt(xVertexLepton*xVertexLepton+yVertexLepton*yVertexLepton)-sqrt(xVertexProton*xVertexProton+yVertexProton*yVertexProton)")
          .Define("deltaZVertex","zVertexLepton-zVertexProton")
          .Define("deltaTheta_p_ep","(p_p4_meas.Theta()-ep_p4_meas.Theta())*180/3.14159")
          .Define("deltaPhi_p_ep","(p_p4_meas.Phi()-ep_p4_meas.Phi())*180/3.14159")
          .Define("em_isFCAL","em_efcal > em_ebcal")
          .Define("ep_isFCAL","ep_efcal > ep_ebcal");

  auto rdf_fiducial=rdf_def.Filter("(jpsi_m_pair > 2) && (jpsi_m_pair < 3.5)")
          .Filter("(zVertex > 51) && (zVertex < 79) && (xVertex*xVertex + yVertex*yVertex < 1*1)")
          .Filter("(em_theta > 2) && (ep_theta > 2) && (p_theta > 2)")
          .Filter("(p_em > 0.4) && (p_ep > 0.4) && (p_p > 0.4)");

  string track_shower_cut="";
  switch (stoi(track_shower_sel)){
    case 0:
      track_shower_cut="(unusedShowers == 0) && (unusedTracks == 0)";
      break;
    case 1:
      track_shower_cut="(unusedShowers == 0) && (unusedTracks < 2)";
      break;
    case 2:
      track_shower_cut="(unusedShowers < 2) && (unusedTracks == 0)";
      break;
    case 3:
      track_shower_cut="(unusedShowers < 2) && (unusedTracks < 2)";
      break;
    default:
      track_shower_cut="(unusedShowers < 10) && (unusedTracks < 10)";
      break;
  }
  auto rdf_track_shower=rdf_fiducial.Filter(track_shower_cut.c_str());

  auto rdf_preShower=rdf_track_shower.Filter("(em_eprebcal_sinTheta > " +preBCut + ") || (em_efcal > 0 )")
          .Filter("(ep_eprebcal_sinTheta > " +preBCut + ") || (ep_efcal > 0 )");

  auto rdf_PoverE=rdf_preShower.Filter("(passEoverP(em_p4_kin.P(),em_efcal,em_ebcal,em_isFCAL," + sigma_minus + "," + sigma_plus + "))")
          .Filter("(passEoverP(ep_p4_kin.P(),ep_efcal,ep_ebcal,ep_isFCAL," + sigma_minus + "," + sigma_plus + "))");

  auto rdf_final=rdf_PoverE.Filter("(fabs(E_miss_stat_p) < " + EMissCut + ")");

  rdf_final.Snapshot(treeName.c_str(),outFileName.c_str());
  cout << "Finished!\n";
}
