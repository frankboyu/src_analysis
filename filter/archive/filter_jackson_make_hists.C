#import "hist_models.h"

void make_hists(string nucleus = "deuterium", string mode = "data")
{

  using namespace ROOT;
  using namespace std;
  EnableImplicitMT();
  
  string inFileName = "../files/trees/tree_filter_" + nucleus + "_" + mode + ".root";
  string histName = "../files/hists/hists_final_" + nucleus + "_" + mode + ".root";
  
  TChain chain("rhoMinus_p_p_filtered");
  chain.Add(inFileName.c_str());

  RDataFrame rdf_raw(chain);

  auto rdf_def = rdf_raw
    .Define("pCMx","CM_p4_kin.X()")
    .Define("pCMy","CM_p4_kin.Y()")
    .Define("perphat","beam_p4_kin.Vect().Cross(struck_p4_kin.Vect()).Unit()")
    .Define("pCMperp","CM_p4_kin.Vect().Dot(perphat)")
    .Define("cosGamma","cos(rec_p4_kin.Angle(struck_p4_kin.Vect()))")
    .Define("alpharec","rec_p4_kin.Minus()/N_p4.M()")
    .Define("alphamiss","struck_p4_kin.Minus()/N_p4.M()")
    .Define("thetarec","rec_p4_kin.Theta()*180/M_PI")
    .Define("thetamiss","struck_p4_kin.Theta()*180/M_PI")
    .Define("minust","-t")
    .Define("p_pim","pim_p4_kin.P()")
    .Define("theta_pim","pim_p4_kin.Theta()*180/M_PI")
    .Define("p_pi0","pi0_p4_kin.P()")
    .Define("theta_pi0","pi0_p4_kin.Theta()*180/M_PI")
    .Define("p_rho","rho_p4_kin.P()")
    .Define("theta_rho","rho_p4_kin.Theta()*180/M_PI")
    .Define("p_p","p_p4_kin.P()")
    .Define("theta_p","p_p4_kin.Theta()*180/M_PI")
    .Define("p_rec","rec_p4_kin.P()")
    .Define("theta_rec","rec_p4_kin.Theta()*180/M_PI");

  auto rdf_mass_filter = rdf_def.Filter("(m2pi_kin > 0.65) && (m2pi_kin < 0.9)").Filter("p_p>1.2");
  auto rdf_1 = rdf_mass_filter.Filter("(t < -1) && (u < -1)");
  auto rdf_15 = rdf_mass_filter.Filter("(t < -1.5) && (u < -1.5)");
  auto rdf_2 = rdf_mass_filter.Filter("(t < -2) && (u < -2)");
  auto rdf_25 = rdf_mass_filter.Filter("(t < -2.5) && (u < -2.5)");
  auto rdf_3 = rdf_mass_filter.Filter("(t < -3) && (u < -3)");
  auto rdf_35 = rdf_mass_filter.Filter("(t < -3.5) && (u < -3.5)");
  auto rdf_4 = rdf_mass_filter.Filter("(t < -4) && (u < -4)");

  //auto rdf = rdf_1;

  vector<decltype(rdf_1)> rdfs = {rdf_1,rdf_15,rdf_2,rdf_25,rdf_3,rdf_35,rdf_4};
  vector<string> tcuts = {"1","1.5","2","2.5","3","3.5","4"};
  
  // Fill and write histograms
  TFile * histFile = new TFile(histName.c_str(),"RECREATE");
  histFile->cd();
      
  TH1D h_t = *rdf_mass_filter.Histo1D(t_model,"minust","accidweight");
  h_t.Write();
      
  for (int i = 0; i < 7; i++)
    {
      auto rdf = rdfs[i];
      string tcut = tcuts[i];
      
      TH1D h_pCMx = *rdf.Histo1D(pCMx_model,"pCMx","accidweight");
      h_pCMx.SetName(("pCMx_" + tcut).c_str());
      h_pCMx.Write();
      
      TH1D h_pCMy = *rdf.Histo1D(pCMy_model,"pCMy","accidweight");
      h_pCMy.SetName(("pCMy_" + tcut).c_str());
      h_pCMy.Write();
      
      TH1D h_pCMperp = *rdf.Histo1D(pCMperp_model,"pCMperp","accidweight");
      h_pCMperp.SetName(("pCMperp_" + tcut).c_str());
      h_pCMperp.Write();
      
      TH1D h_rec = *rdf.Histo1D(rec_model,"pRec","accidweight");
      h_rec.SetName(("pRec_" + tcut).c_str());
      h_rec.Write();
      
      TH1D h_kmiss = *rdf.Histo1D(kmiss_model,"kmiss","accidweight");
      h_kmiss.SetName(("kmiss_" + tcut).c_str());
      h_kmiss.Write();
      
      TH1D h_cosGamma = *rdf.Histo1D(cosGamma_model,"cosGamma","accidweight");
      h_cosGamma.SetName(("cosGamma_" + tcut).c_str());
      h_cosGamma.Write();
      
      TH2D h_kmiss_rec = *rdf.Histo2D(kmiss_rec_model,"kmiss","pRec","accidweight");
      h_kmiss_rec.SetName(("kmiss_rec_" + tcut).c_str());
      h_kmiss_rec.Write();
      
      TH1D h_alpharec = *rdf.Histo1D(alpharec_model,"alpharec","accidweight");
      h_alpharec.SetName(("alpharec_" + tcut).c_str());
      h_alpharec.Write();
      
      TH1D h_alphamiss = *rdf.Histo1D(alphamiss_model,"alphamiss","accidweight");
      h_alphamiss.SetName(("alphamiss_" + tcut).c_str());
      h_alphamiss.Write();
      
      TH1D h_thetarec = *rdf.Histo1D(thetarec_model,"thetarec","accidweight");
      h_thetarec.SetName(("thetarec_" + tcut).c_str());
      h_thetarec.Write();
      
      TH1D h_thetamiss = *rdf.Histo1D(thetamiss_model,"thetamiss","accidweight");
      h_thetamiss.SetName(("thetamiss_" + tcut).c_str());
      h_thetamiss.Write();
      
      TH1D h_alphaCM = *rdf.Histo1D(alphaCM_model,"alphaCM_kin","accidweight");
      h_alphaCM.SetName(("alphaCM_" + tcut).c_str());
      h_alphaCM.Write();

      if (i==1)
	{
	  TH2D h_kin_pim = *rdf.Histo2D(kin_model,"p_pim","theta_pim");
	  h_kin_pim.SetName("kin_pim");
	  h_kin_pim.Write();
	    
	  TH2D h_kin_pi0 = *rdf.Histo2D(kin_model,"p_pi0","theta_pi0");
	  h_kin_pi0.SetName("kin_pi0");
	  h_kin_pi0.Write();
	    
	  TH2D h_kin_rho = *rdf.Histo2D(kin_model,"p_rho","theta_rho");
	  h_kin_rho.SetName("kin_rho");
	  h_kin_rho.Write();
	    
	  TH2D h_kin_p = *rdf.Histo2D(kin_model,"p_p","theta_p");
	  h_kin_p.SetName("kin_p");
	  h_kin_p.Write();
	    
	  TH2D h_kin_rec = *rdf.Histo2D(kin_model,"p_rec","theta_rec");
	  h_kin_rec.SetName("kin_rec");
	  h_kin_rec.Write();
	}

    }
  
  histFile->Close();
  
}
