#include <iostream>
#include <stdio.h>
#include <string.h>
#include <ROOT/RDataFrame.hxx>
#include <TH2D.h>
#include <TH1.h>
#include <ROOT/RDF/HistoModels.hxx>

using namespace ROOT;
using namespace std;
using namespace ROOT::Detail::RDF;



const double leftmass = 0.57;
const double rightmass = 0.97;

double rhoyield_anderror(TH1D *rhomass_hist, double &error) {
	const int leftbin = rhomass_hist->FindBin(leftmass);
	const int rightbin = rhomass_hist->FindBin(rightmass);	
	const int diffbin = rightbin - leftbin + 1;
	const double yleft = rhomass_hist->GetBinContent(leftbin);
	const double yright = rhomass_hist->GetBinContent(rightbin);
	const double yleft_error = rhomass_hist->GetBinError(leftbin);
	const double yright_error = rhomass_hist->GetBinError(rightbin);
	
	double integralerror;
	double integral = rhomass_hist->IntegralAndError(leftbin, rightbin, integralerror);
	error = sqrt(integralerror * integralerror + ((0.5 * 0.5) / (diffbin * diffbin) * (yleft_error * yleft_error + yright_error * yright_error)));
	return integral - 0.5 * (yright + yleft) * diffbin;
}


void filter_rho0_1p(string inputfilename, string treeName, string outfilename, string outhistname){
	TH1::SetDefaultSumw2();
	//using namespace ROOT;
	//using namespace std;
  
	TChain chain(treeName.c_str());
  chain.Add(inputfilename.c_str());

  RDataFrame rdf_raw(chain);

  cout << "Defining variables...\n";
	auto rdf_def = rdf_raw
		.Define("N_p4","TLorentzVector(0,0,0,0.93892)")
		.Define("prot_mom","p_p4_kin.P()")
		.Define("pipProt_p4_kin","TLorentzVector(p_p4_kin.Vect(),sqrt(p_p4_kin.Vect().Mag2() + 0.139570*0.139570))")
		.Define("pipProt_mom","pipProt_p4_kin.P()")
		.Define("pim_mom","pim_p4_kin.P()")
		.Define("pip_mom","pip_p4_kin.P()")
		.Define("prot_theta","180./M_PI * p_p4_kin.Theta()")
		.Define("pim_theta","180./M_PI *pim_p4_kin.Theta()")
		.Define("pip_theta","180./M_PI * pip_p4_kin.Theta()")
		//.Define("prot_z","prot_x4_kin.Z()")
		.Define("rho0_p4_kin","pip_p4_kin + pim_p4_kin")
		.Define("rho0_mass","rho0_p4_kin.M()")
		.Define("pmiss_p4","p_p4_kin + rho0_p4_kin - beam_p4_kin") 
		.Define("pmiss_mom","pmiss_p4.Vect().Mag()")
		.Define("pmiss_T","pmiss_p4.Pt()")
		.Define("t","(beam_p4_kin - rho0_p4_kin).M2()")
		//.iDefine("t_1","(other leg... 
		.Define("t_3M","(beam_p4_kin - rho0_p4_kin - pipProt_p4_kin).M2()")
		.Define("t_abs","abs(t)")
		//.Define("kmiss","N_p4.M()*sqrt((N_p4.M2() + p4_prot_kin.Perp2())/(p4_prot_kin.Minus()*(2*N_p4.M() - p4_prot_kin.Minus()))-1)")
		.Define("kmiss","N_p4.M()*sqrt((N_p4.M()*N_p4.M()+pmiss_T*pmiss_T)/(pmiss_p4.Minus()*(2*N_p4.M() - pmiss_p4.Minus()))-1)")
		.Define("y2pi","-t/(2* N_p4.M()*(beam_p4_kin.E() - rho0_p4_kin.E()))")
		.Define("y3pi","-t_3M/(2* N_p4.M()*(beam_p4_kin.E() - rho0_p4_kin.E() - pipProt_p4_kin.E()))")
		//.Define("W_2pi","(beam_p4_kin.P() + prot_mom - rho0_p4_kin.P())")
		.Define("W2_2pi","t +2*N_p4*(beam_p4_kin - rho0_p4_kin)+N_p4*N_p4")
		//.Define("W_3pi","(beam_p4_kin.P() + prot_mom - (rho0_p4_kin.P() + pipProt_p4_kin.P()))")
		//.Define("W2_3pi","W_3pi * W_3pi")
		.Define("deltapp", "(pip_p4_kin + p_p4_kin).M2()")
		.Define("delta0", "(pim_p4_kin + p_p4_kin).M2()")
		.Define("MX", "(beam_p4_kin - rho0_p4_kin - p_p4_kin + 2* N_p4).M2()")
		.Define("EX", "(beam_p4_kin - rho0_p4_kin - p_p4_kin).E()")
		.Define("piondiff","pim_p4_kin.P() - pip_p4_kin.P()")
		.Define("MXr", "double v=MX; if(MX<0) {v += 0.855;} return (const double)v;")
		.Define("RhomassX","double v = rho0_mass; if(MX<0) {v = sqrt(v*v + 0.855);} return v;")
		.Define("pmiss_theta","180./M_PI * pmiss_p4.Theta()");		


	cout << "Defining Filters...\n";
	auto rdf_no_filter = rdf_def;
	auto rdf_t1_filtered = rdf_no_filter.Filter("t_abs > 1");
	auto rdf_pmiss4_filtered = rdf_t1_filtered.Filter("pmiss_mom >0.4");
	auto rdf_protA10_filtered = rdf_pmiss4_filtered.Filter("prot_theta > 10");
	//auto rdf_mx3_filtered = rdf_protA10_filtered.Filter("MX < 3");
	//auto rdf_zvertex_filtered = rdf_no_filter.Filter("prot_z > 50 || prot_z < 80");	
	//auto rdf_MM_filter= rdf_no_filter.Filter("MX > 0"); 
	//auto rdf_y3pi_filter = rdf_no_filter.Filter("y3pi != 0");
	//auto rdf_Wypi_filter = rdf_protA10_filtered.Filter("W2_2pi < -3.2 * y3pi + 4.2");
	//auto rdf_Wypi_filter = rdf_y3pi_filter.Filter("W2_2pi < exp(2.3*(-1* y3pi + 0.85))");
	/*auto rdf_Wypi_filter = rdf_y3pi_filter.Filter([](const double &W2_2pi) {
		return W2_2pi < exp(2.3*(-1* y3pi + 0.85));
		},{"y3pi"});*/
	//auto rdf_rhocut = rdf_Wypi_filter.Filter("rho0_mass < 1");

	//auto rdf_mxW2_filtered = rdf_protA10_filtered.Filter("W2_2pi< -0.4* MX +4");

	auto rdf_mxy3pi_filtered = rdf_protA10_filtered.Filter("y3pi < -0.6* MX + 1.5");
	auto rdf_pipprot_filtered = rdf_mxy3pi_filtered.Filter([](const TLorentzVector &pip_p4_kin) {
   	return pip_p4_kin.P() < exp(-0.085 * 180.0 / TMath::Pi() * pip_p4_kin.Theta() + 2.85) - 0.1;
   	}, {"pip_p4_kin"});
	auto rdf_Deltapp_filtered = rdf_pipprot_filtered.Filter("deltapp > 1.3");
	auto rdf_Delta0_filtered = rdf_Deltapp_filtered.Filter("delta0 > 1.3");
	//auto rdf_diffbkgrnd_filtered = rdf_Delta0_filtered.Filter("yM");
	auto rdf_t15_filtered = rdf_Delta0_filtered.Filter("t_abs >  1.5");
	//auto rdf_kmiss_filtered = rdf_t15_filtered.Filter("kmiss > 0.4");
	auto rdf_rho0mass_filtered = rdf_t15_filtered.Filter("rho0_mass < 0.97 && rho0_mass > 0.57");

	TH2D* h_rho0mass;

	int N_filters = 10;
	ROOT::RDF::RNode rdfs [] = {rdf_no_filter ,rdf_t1_filtered , rdf_pmiss4_filtered, rdf_protA10_filtered, /*rdf_mx3_filtered }; rdf_mxW2_filtered,*/rdf_mxy3pi_filtered, /*rdf_zvertex_filtered,*/ rdf_pipprot_filtered, rdf_Deltapp_filtered,rdf_Delta0_filtered, rdf_t15_filtered, rdf_rho0mass_filtered};

	string labels[] = {"no_cuts", "t>1_filter","pmiss4_filter", "ProtAngle_10D_filter" /*,"Mx<3"};, "mxW2_filter"*/, "MXvsy3pi_filter",  /* "Z_vertex_filter",*/"Pi+Prot_confusion_filter" ,"DeltaPP_filter", "Delta0_filter", "t>1.5_filter",  "rho0Mass_filter"};

	cout << "Creating histogram file...\n";
	TFile * histfile = new TFile (outhistname.c_str(),"RECREATE");
	histfile->cd();

	cout << "Constructing histograms...\n";
	vector <TH1D*> histslist;
	vector<ROOT::RDF::RResultPtr<TH1D>> hists;
  vector<ROOT::RDF::RResultPtr<TH2D>> hists2d;
	TDirectory *rhomass = histfile->mkdir("Rho0Mass_Hists");
	/*TDirectory *fiducials = histfile->mkdir("Fiducial_Hists");
	TDirectory *deltas = histfile->mkdir("Delta_Hists");*/
	TDirectory *physics = histfile->mkdir("Physics_Hists");
 	
	for (int i = 0; i < N_filters; i++) {
			auto rdf = rdfs[i];
    	string label = labels[i];
    	cout << label << "\n";
		

			//fiducials->cd();
			//auto h_pip = rdf.Histo1D({("Pi+_" + label).c_str(), ";PiP Mom [GeV]",100,0,10}, "pip_mom","accidweight");
			//hists.push_back(h_pip);
			//auto h_pim = rdf.Histo1D({("Pi-_" + label).c_str(), ";PiM Mom [GeV]",100,0,10}, "pim_mom","accidweight");
			//hists.push_back(h_pim);
			//auto h_prot = rdf.Histo1D({("Prot_" + label).c_str(), ";Proton Mom [GeV]",100,0,10}, "prot_mom","accidweight");
			//hists.push_back(h_prot);
			//auto h_prot_z = rdf.Histo1D({("Prot_Z_" + label).c_str(), ";Proton z vertex [cm]; Counts",100,30,100},"prot_z","accidweight");
			//auto h_prot_z = rdf.Histo1D("prot_z");//,"accidweight");
			//hists.push_back(h_prot_z);
			auto h_prot_theta = rdf.Histo2D({("Prot_Theta_" + label).c_str(), "; Theta [Degrees];Proton Mom [GeV]",100,0,180,100,0,12}, "prot_theta","prot_mom","accidweight");
			hists2d.push_back(h_prot_theta);
			auto h_pim_theta = rdf.Histo2D({("Pi-_Theta_" + label).c_str(), "; Theta [Degrees];PiM Mom [GeV]",100,0,180,100,0,12}, "pim_theta","pim_mom","accidweight");
			hists2d.push_back(h_pim_theta);
			auto h_pip_theta = rdf.Histo2D({("Pi+_Theta_" + label).c_str(), "; Theta [Degrees];PiP Mom [GeV]",100,0,180,100,0,12},"pip_theta", "pip_mom","accidweight");
			hists2d.push_back(h_pip_theta);
			
			/*auto h_piptheta_t = rdf.Histo2D({("Pi+Theta_t_" + label).c_str(),";pi+Theta [degrees]; t [GeV]",100,0,180,100,0,10},"pip_theta","t_abs","accidweight");
			hists2d.push_back(h_piptheta_t);
			auto h_pimtheta_t = rdf.Histo2D({("Pi-Theta_t_" + label).c_str(),";pi-Theta [degrees]; t [GeV]",100,0,180,100,0,10},"pim_theta","t_abs","accidweight");
			hists2d.push_back(h_pimtheta_t);
			auto h_prottheta_t = rdf.Histo2D({("ProtTheta_t_" + label).c_str(),";protTheta [degrees]; t [GeV]",100,0,180,100,0,10},"prot_theta","t_abs","accidweight");
			hists2d.push_back(h_prottheta_t);*/
			//auto h_pipProt = rdf.Histo1D({("Pi+M_ProtVect_" + label).c_str(),";pipM_ProtVect [GeV];Counts",100,0,10},"pipProt_mom");
			//hists.push_back(h_pipProt);
			//fiducials->Close();

			//deltas->cd();
			auto h_deltapp = rdf.Histo1D({("DeltaPP_" + label).c_str(), ";DeltaPP [GeV]", 100,0,5},"deltapp","accidweight");
			hists.push_back(h_deltapp);
			auto h_delta0 = rdf.Histo1D({("Delta0_" + label).c_str(), ";Delta0 [GeV]", 100,0,5},"delta0","accidweight");
			hists.push_back(h_delta0);
			//TH2DModel delta_dalitz_model(("delta_dalitz_" + label).c_str(), ";Pip + Prot;Pim + Prot",100,0,3,100,0,3); exit (0);
			//auto h_Dalitz_delta = rdf.Histo2D<double,double>(delta_dalitz_model,"deltapp","delta0","accidweight");
			auto h_Dalitz_deltapp_rho0 = rdf.Histo2D({("DeltaPP_vs_Rho0_" + label).c_str(), ";rho0;DeltaPP;",100,0,5,100,0,5},"rho0_mass","deltapp","accidweight");
			hists2d.push_back(h_Dalitz_deltapp_rho0);
			auto h_Dalitz_delta0_rho0 = rdf.Histo2D({("Delta0_vs_Rho0_" + label).c_str(), ";rho0;Delta0",100,0,5,100,0,5},"rho0_mass","delta0","accidweight");
			hists2d.push_back(h_Dalitz_delta0_rho0);
			auto h_Dalitz_deltaS = rdf.Histo2D({("DeltaPP_vs_Delta0_" + label).c_str(), ";DeltaPP;Delta0",100,0,5,100,0,5},"deltapp","delta0","accidweight");
			hists2d.push_back(h_Dalitz_deltaS);

			//rhomass->cd();
			auto h_rho0_mass = rdf.Histo1D({("Rho0Mass_" + label).c_str(),";Rho0 Mass [GeV]",100,0,3},"rho0_mass","accidweight");
			hists.push_back(h_rho0_mass);
			auto h_rho0_pmiss = rdf.Histo2D({("Rho0Mass_vs_Pmiss_" + label).c_str(), ";rho0mass [GeV];pmiss [GeV]",100,0,3,100,0.4,1},"rho0_mass","pmiss_mom","accidweight");
			hists2d.push_back(h_rho0_pmiss);
			auto h_rho0_pmiss_binned = rdf.Histo2D({("Rho0Mass_vs_Pmiss_binned_" + label).c_str(), ";rho0mass [GeV];pmiss [GeV]",90,0,3,4,0.4,0.8},"rho0_mass","pmiss_mom","accidweight");
			hists2d.push_back(h_rho0_pmiss_binned);
			auto h_rho0_kmiss = rdf.Histo2D({("Rho0Mass_vs_Kmiss_" + label).c_str(), ";rho0mass [GeV]; kmiss [GeV]", 100,0,3,100,0,3},"rho0_mass","kmiss","accidweight");
			hists2d.push_back(h_rho0_kmiss);
			auto h_rho0_t = rdf.Histo2D({("Rho0_t_" + label).c_str(), ";rho0 mass [GeV]; t [GeV]", 100,0,3, 100,0,10},"rho0_mass","t_abs","accidweight");
			hists2d.push_back(h_rho0_t);

			//physics->cd();
			auto h_pmiss = rdf.Histo1D({("Pmiss_" + label).c_str(), ";Pmiss [GeV];Counts",100,0.4, 1},"pmiss_mom","accidweight");
			hists.push_back(h_pmiss);
			auto h_pmiss_kmiss = rdf.Histo2D({("Pmiss_Kmiss_" + label ).c_str(), ";pmiss [GeV];kmiss [GeV]", 100, 0,3,100,0,3},"pmiss_mom","kmiss","accidweight");
			hists2d.push_back(h_pmiss_kmiss);
			//hists2d.push_back(h_rho0_kmiss);
			auto h_kmiss = rdf.Histo1D({("Kmiss_" +label).c_str(), ";kmiss [GeV]; Counts", 100,0,3},"kmiss","accidweight");
			hists.push_back(h_kmiss);
			//auto h_yM = rdf.Histo1D("yM","accidweight");
			//auto h_y2pi = rdf.Histo1D({("y2Pi_" + label).c_str(), ";y2pi;Counts", 100,-2,4},"y2pi");//,"accidweight");
			//hists.push_back(h_y2pi);
			auto h_y3pi = rdf.Histo1D({("y3Pi_" + label).c_str(), ";y3pi;Counts", 100,-2,4},"y3pi");//,"accidweight");
			hists.push_back(h_y3pi);
			

			auto h_W2_2pi = rdf.Histo1D({("W2_2pi_" + label).c_str(), ";W2_2pi;Counts", 100,-5,20},"W2_2pi");//,"accidweight");
			hists.push_back(h_W2_2pi);
			//auto h_W2_3pi = rdf.Histo1D({("W2_3pi" + label).c_str(), ";W2_3pi;Counts", 100,-2,4},"W2+3pi");//,"accidweight");
			//hists.push_back(h_W2_3pi);
			auto h_y3pi_W2_2pi = rdf.Histo2D({("y3pi_W2_2pi_" + label).c_str(), ";y3pi;W2_2pi;Counts",100,-2,4,100,-5,20},"y3pi","W2_2pi");
			hists2d.push_back(h_y3pi_W2_2pi);

			//auto h_t = rdf.Histo1D({("t_" + label).c_str(), ";t [GeV]; Counts", 100, 0, 10}, "t", "accidweight");
			auto h_t = rdf.Histo1D({("t_"+ label).c_str(),";t [GeV];Counts",100,1.5,10},"t_abs","accidweight");
			hists.push_back(h_t);
			auto h_t_unweighted = rdf.Histo1D({("t_unweighted_"+ label).c_str(),";t [GeV];Counts",100,1.5,10},"t_abs");
			hists.push_back(h_t_unweighted);
			auto h_t_pmiss = rdf.Histo2D({("t_pmiss_"+ label).c_str(),"; t[GeV];pmiss [GeV]",100,1.5,10,100,0.4,1},"t_abs","pmiss_mom","accidweight");
			hists2d.push_back(h_t_pmiss);
			auto h_mx = rdf.Histo1D({("MX_"+ label).c_str(),"; MX",100,-2,5},"MX","accidweight");
			hists.push_back(h_mx);	
			auto h_rho_Mx = rdf.Histo2D({("RhoMass_MX_" + label).c_str(),";rhomass;MX",100,0,3,100,-5,5},"rho0_mass","MX","accidweight");
			hists2d.push_back(h_rho_Mx);
			auto h_ex = rdf.Histo1D({("EX_"+ label).c_str(),"; EX",100,-5,5},"EX","accidweight");
			hists.push_back(h_ex);	
			auto h_rho_Ex = rdf.Histo2D({("RhoMass_EX_" + label).c_str(),";rhomass;EX",100,0,3,100,-5,5},"rho0_mass","EX","accidweight");
			hists2d.push_back(h_rho_Ex);
			auto h_mx_w2 = rdf.Histo2D({("MX_W2_" + label).c_str(),";MX;W2_2pi",150,-2,25,100,-5,20},"MX","W2_2pi","accidweight");
			hists2d.push_back(h_mx_w2);
			auto h_mx_y3pi = rdf.Histo2D({("MX_y3pi_" + label).c_str(),";MX;Y3pi",100, -2,10, 100,-2,4},"MX","y3pi","accidweight");
			hists2d.push_back(h_mx_y3pi);
			/*auto h_piondiff = rdf.Histo1D({("piondiff_"+label).c_str(),";pim-pip",100,-10,10},"piondiff","accidweight");
			hists.push_back(h_piondiff);
			auto h_mxR = rdf.Histo1D({("MXR_"+ label).c_str(),"; MXr",100,-5,5},"MXr","accidweight");
			hists.push_back(h_mxR);	
			auto h_rho_MxR = rdf.Histo2D({("RhoMass_MXR_" + label).c_str(),";rhomass;MXr",100,0,3,100,-5,5},"rho0_mass","MXr","accidweight");
			hists2d.push_back(h_rho_MxR);
			auto h_RhomassX = rdf.Histo1D({("RhomassX_"+ label).c_str(),"; RhomassX",100,0,3},"RhomassX","accidweight");
			hists.push_back(h_RhomassX);	
			auto h_RhomassX_Mx = rdf.Histo2D({("RhoMass_MXR_" + label).c_str(),";RhomassX;MX",100,0,3,100,-5,5},"RhomassX","MX","accidweight");
			hists2d.push_back(h_RhomassX_Mx);
			auto h_pmiss_Theta = rdf.Histo1D({("Pmiss_theta_" + label).c_str(),";pmiss theta",100,0,180},"pmiss_theta","accidweight");
			hists.push_back(h_pmiss_Theta);*/


	/*		TF1 *f1 = new TF1("f1"," exp(-0.085 * x + 2.85) - 0.1",0,80);
			f1->SetLineWidth(2);
			f1->SetLineColor(kRed);
			auto c1=new TCanvas("","",1400,1400);c1->Divide(2,2);
			c1->cd(1);
      h_prot_theta->Draw("colz");
			f1 ->Draw("same");
      c1->cd(2);
			h_pim_theta->Draw("colz");
			f1->Draw("same");
			c1->cd(3);
			h_pip_theta->Draw("colz");
			f1->Draw("same");
			c1->SaveAs("mom_theta_f1_data.png");// + label).c_str()".png");

			auto c2=new TCanvas("","",1400,1400);c2->Divide(2,2);
			c2->cd(1);
      h_piptheta_t->Draw("colz");
      c2->cd(2);
			h_piptheta_t->Draw("colz");
			c2->cd(3);
			h_piptheta_t->Draw("colz");
			c2->SaveAs("theta_f1_data.png");// + label).c_str()".png"); */

	/*		TF1 *f2 = new TF1("f2", "exp(2.3*(-1* y3pi + 0.85))",-2,4;
			f2->SetLineWidth(2);
			f2->SetLineColor(kRed);
			auto c3 =new TCanvas("","",1400,1400);c3->Divide(2,2);
			c3->cd(1);
			h_y3pi_W2_2pi_no_cuts->Draw("colz");
			f2->Draw("same");
			c3->cd(2);
			h_y3pi_W2_2pi_Wypi->Draw("colz");
			f2->Draw("same");
			c3->cd(3);
			h_y3pi_W2_2pi_rhocut->Draw("colz");
			f2->Draw("same");
			c3->SaveAs("w2_y3pi_cut.png");

		exit(0);*/
			if (i == N_filters - 1 )  {
			h_rho0mass = (TH2D *) h_rho0_pmiss_binned-> Clone(); 
			};
			
}

  ////////////////////// Count the Rho0s//////////////////////

  ofstream outfile("output.txt");
  outfile << "#Yield vs PMiss\n"
          << "#[left edge] [bin center] [right edge] [counts] [yield] [error] \n";
  for (int i = 1; i <= h_rho0mass->GetNbinsY(); i++) {
      // projections
    TString temp = Form("proj%d",i);
    TH1D *thismass_spectrum = h_rho0mass->ProjectionX(temp.Data(), i, i);
    double error;
    double yield = rhoyield_anderror(thismass_spectrum, error);
    outfile << h_rho0mass->GetYaxis()->GetBinLowEdge(i) << " "
            << h_rho0mass->GetYaxis()->GetBinCenter(i) << " "
            << h_rho0mass->GetYaxis()->GetBinUpEdge(i) << " "
            /*<< h_rho0mass->GetYaxis()->GetBinContent(i) */<< " " << yield << " " << error << " \n"; 
  }
  outfile.close();

  // Writing out histograms
 	cout << "Writing out histograms...\n";

  for (ROOT::RDF::RResultPtr<TH1D> hist : hists)
    {
			histfile->cd("Physics_Hists");
			hist->Write();
    }
  for (ROOT::RDF::RResultPtr<TH2D> hist : hists2d)
    {
      histfile->cd("Rho0Mass_Hists");
			hist->Write();
    }

	
//f1->Write();
  histfile->Close();

/*	histfile->cd();
	histfile->Write();
	histfile->Close();*/

  // Take snapshot
  //cout << "Saving snapshot...\n";
  //rdf_Delta0_filtered.Snapshot("rho0_1p_filtered",outfilename);
}

