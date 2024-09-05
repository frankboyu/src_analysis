#include <iostream>
#include <cmath>
#include <stdio.h>
#include <string.h>

using namespace std;
using namespace ROOT;
using namespace RooFit;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_deuteron    = 1.875612;
double mass_proton      = 0.938272;
double mass_pion        = 0.139570;
double mass_phi         = 1.019461;
double RadToDeg         = 180.0 / 3.14159265;

void filter_phi_c_2H_test()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/flattree_phi_c_2H_data_ver02/*090213.root";
    string InputTree  = "flattree_phi_c_2H_data";
    string OutputFile = "output/filteredtree_phi_c_2H_data.root";
    string OutputTree = "filteredtree_phi_c_2H_data";

    EnableImplicitMT();

    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());
    RDataFrame rdf_raw(chain);

    vector<TH1*> hist_list;

    auto rdf_def = rdf_raw
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4","TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("phi_p4_meas","kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin","kp_p4_kin + km_p4_kin")
    .Define("missd_p4_meas","beam_p4_meas + target_p4 - phi_p4_meas")  // missd_p4_kin already defined in default branches
    // .Define("pip_p4_meas","TLorentzVector(kp_p4_meas.X(), kp_p4_meas.Y(), kp_p4_meas.Z(), sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pip_p4_kin","TLorentzVector(kp_p4_kin.X(), kp_p4_kin.Y(), kp_p4_kin.Z(), sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_meas","TLorentzVector(km_p4_meas.X(), km_p4_meas.Y(), km_p4_meas.Z(), sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_kin","TLorentzVector(km_p4_kin.X(), km_p4_kin.Y(), km_p4_kin.Z(), sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("beam_p4com_meas","beam_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("beam_p4com_kin","beam_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    // .Define("phi_p4com_meas","phi_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("phi_p4com_kin","phi_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    .Define("sqrt_s_meas", "(beam_p4_meas + target_p4).Mag()")
    .Define("sqrt_s_kin", "(beam_p4_kin + target_p4).Mag()")
    .Define("minus_t_meas", "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(beam_p4_kin - phi_p4_kin).Mag2()")
    .Define("minus_u_meas", "-(beam_p4_meas - missd_p4_meas).Mag2()")
    .Define("minus_u_kin", "-(beam_p4_kin - missd_p4_kin).Mag2()")
    .Define("phi_mass_meas","phi_p4_meas.M()")
    .Define("phi_mass_kin","phi_p4_kin.M()")
    .Define("phi_theta_meas","phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin","phi_p4_kin.Theta()*RadToDeg")
    // .Define("rho_mass_meas","pip_p4_meas + pim_p4_meas")
    // .Define("rho_mass_kin","pip_p4_kin + pim_p4_kin")
    .Define("y_phi_meas","minus_t_meas/(2*mass_deuteron*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin","minus_t_kin/(2*mass_deuteron*(beam_p4_kin.E()-phi_p4_kin.E()))")
    // .Define("DeltaE_meas","(s_meas - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_meas)) - phi_p4com_meas.E()")
    // .Define("DeltaE_kin","(s_kin - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_kin)) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    auto rdf_no_filter = rdf_def;
    auto rdf_cl_filtered = rdf_no_filter.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered = rdf_cl_filtered.Filter([](double kp_pidfom, double km_pidfom) {return (kp_pidfom > 0.01) && (km_pidfom > 0.01);}, {"kp_pidfom","km_pidfom"});
    auto rdf_y_phi_filtered = rdf_pidfom_filtered.Filter([](double y_phi_meas) {return y_phi_meas > 0.4;}, {"y_phi_meas"});
    auto rdf_phi_mass_filtered = rdf_y_phi_filtered.Filter([](double phi_mass_meas) {return phi_mass_meas > 1.01 && phi_mass_meas < 1.03;}, {"phi_mass_meas"});

    rdf_phi_mass_filtered.Snapshot("filteredtree_phi_c_2H_data",OutputFile);






// TH2DModel alphaCM_cosGammaProxy_model("alphaCM_cosGammaProxy", ";CM Lightcone Fraction;Cos Gamma Proxy",300,0,3,400,-1,1);
// TH2DModel alphaCM_thetaP_model("alphaCM_thetaP", ";CM Lightcone Fraction;Proton Angle [degrees]",300,0,3,360,0,180);
// TH1DModel zVertex_model("zVertex", ";z Vertex [cm];Counts",1000,0,200);
// TH2DModel xyVertex_model("xyVertex", ";x Vertex [cm];y Vertex [cm]",100,-5,5,100,-5,5);
// TH1DModel E_model("lead_energy", ";Lead Energy [GeV];Counts",1200,0,12);
// TH1DModel omega_model("omega_mass", ";Omega Mass [GeV];Counts",500,0,5);
// TH2DModel proton_momenta_model("proton_momenta", ";Lead Momentum [GeV];Recoil Momentum [GeV]",200,0,10,200,0,10);
// TH1DModel mass_model("rho_mass", ";2-Pion Mass [GeV];Counts",500,0,5);
// TH2DModel kmiss_rec_model("kmiss_rec", ";k_miss [GeV];Recoil Momentum [GeV]",200,0,2,200,0,2);
// TH1DModel Emiss_model("Emiss", ";2N Missing Energy [GeV];Counts",200,-10,10);

    int N_filters = 5;
    RNode rdfs [] = {rdf_no_filter, rdf_cl_filtered, rdf_pidfom_filtered, rdf_y_phi_filtered, rdf_phi_mass_filtered};
    string labels [] = {"no_cut", "cl_cut", "pidfom_cut", "y_phi_cut", "phi_mass_cut"};

    // Fill and write histograms
    TFile * histFile = new TFile(OutputFile.c_str(), "update");
    histFile->cd();

    TH1D h_z = *rdf_pidfom_filtered.Histo1D({"zVertex", ";z Vertex [cm];Counts",1000,0,2},"phi_mass_meas","accidweight");
    h_z.Write();


//   for (int i = 0; i < N_filters; i++) {

//     cout << i << "\n";

//     auto rdf = rdfs[i];
//     string label = labels[i];

//     TH2D h_classic = *rdf.Histo2D<double,double>(alphaCM_cosGammaProxy_model,"alphaCM_kin","cosGamma_proxy","accidweight");
//     h_classic.SetName(("classic_" + label).c_str());
//     //hists.push_back(&h_classic);
//     //h_classic.Write();
//     cout << "e\n";
//     TH2D h_angular = *rdf.Histo2D<double,double>(alphaCM_thetaP_model,"alphaCM_kin","thetaP","accidweight");
//     h_angular.SetName(("angular_" + label).c_str());
//     //hists.push_back(&h_angular);
//     //h_angular.Write();
//     cout << "f\n";
//     TH1D h_mass = *rdf.Histo1D<double>(mass_model,"m2pi_kin","accidweight");
//     h_mass.SetName(("mass_" + label).c_str());
//     //hists.push_back(&h_mass);
//     //h_mass.Write();
//     cout << "g\n";

//     h_classic.Write();
//     h_angular.Write();
//     h_mass.Write();

//   }

//   // Before cuts
//   TH1D h_z = *rdf_no_filter.Histo1D(zVertex_model,"zVertex","accidweight");
//   hists.push_back(&h_z);
//   //h_z.Write();
//   TH2D h_xy = *rdf_no_filter.Histo2D(xyVertex_model,"xVertex","yVertex","accidweight");
//   hists.push_back(&h_xy);
//   //h_xy.Write();

//   // After vertex cut
//   TH1D h_energy = *rdf_vertex_filtered.Histo1D(E_model,"eLead","accidweight");
//   hists.push_back(&h_energy);
//   //h_energy.Write();

//   // After energy cut
//   TH1D h_omega1 = *rdf_energy_filtered.Histo1D(omega_model,"omega1_m_meas","accidweight");
//   h_omega1.SetName("omega1_mass");
//   hists.push_back(&h_omega1);
//   //h_omega1.Write();
//   TH1D h_omega2 = *rdf_energy_filtered.Histo1D(omega_model,"omega2_m_meas","accidweight");
//   h_omega2.SetName("omega2_mass");
//   hists.push_back(&h_omega2);
//   //h_omega2.Write();

//   // After diffractive background cut
//   TH2D h_proton_momenta = *rdf_bg_filtered.Histo2D(proton_momenta_model,"pLead","pRec","accidweight");
//   hists.push_back(&h_proton_momenta);
//   //h_proton_momenta.Write();

//   TH2D h_kmiss_rec = *rdf_bg_filtered.Histo2D(kmiss_rec_model,"kmiss","pRec","accidweight");
//   hists.push_back(&h_kmiss_rec);
//   //h_kmiss_rec.Write();

//   // After momentum cuts
//   TH1D h_Emiss = *rdf_rec_filtered.Histo1D(Emiss_model,"Emiss2N","accidweight");
//   hists.push_back(&h_Emiss);
//   //h_Emiss.Write();

//   for (TH1* hist : hists)
//     {
//       hist->Write();
//     }

//   histFile->Close();


}
