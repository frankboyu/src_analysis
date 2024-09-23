#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstring>

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

void filter_phi_d_2H_data()
{
    string InputFile  = "/work/halld2/home/boyu/src_analysis/selection/output/test/flattree_phi_d_2H_data_1_initial.root";
    string InputTree  = "flattree_phi_d_2H_data";
    string OutputFile = "output/filteredtree_phi_d_2H_data.root";
    string OutputTree = "filteredtree_phi_d_2H_data";

    // Read input files
    cout << "Reading input files...\n";
    TChain chain(InputTree.c_str());
    chain.Add(InputFile.c_str());

    // Define data frame
    cout << "Defining data frame...\n";
    RDataFrame rdf_raw(chain);

    auto rdf_def = rdf_raw
    .Define("kin_cl","TMath::Prob(kin_chisq,kin_ndf)")
    .Define("target_p4","TLorentzVector(0, 0, 0, mass_deuteron)")
    .Define("phi_p4_meas","kp_p4_meas + km_p4_meas")
    .Define("phi_p4_kin","kp_p4_kin + km_p4_kin")
    // .Define("missd_p4_meas","beam_p4_meas + target_p4 - phi_p4_meas")  // missd_p4_kin already defined in default branches
    // .Define("pip_p4_meas","TLorentzVector(kp_p4_meas.X(), kp_p4_meas.Y(), kp_p4_meas.Z(), TMath::Sqrt(kp_p4_meas.P()*kp_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pip_p4_kin","TLorentzVector(kp_p4_kin.X(), kp_p4_kin.Y(), kp_p4_kin.Z(), TMath::Sqrt(kp_p4_kin.P()*kp_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_meas","TLorentzVector(km_p4_meas.X(), km_p4_meas.Y(), km_p4_meas.Z(), TMath::Sqrt(km_p4_meas.P()*km_p4_meas.P() + mass_pion*mass_pion))")
    // .Define("pim_p4_kin","TLorentzVector(km_p4_kin.X(), km_p4_kin.Y(), km_p4_kin.Z(), TMath::Sqrt(km_p4_kin.P()*km_p4_kin.P() + mass_pion*mass_pion))")
    // .Define("beam_p4com_meas","beam_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("beam_p4com_kin","beam_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    // .Define("phi_p4com_meas","phi_p4_meas.Boost(-(beam_p4_meas + target_p4).BoostVector())")
    // .Define("phi_p4com_kin","phi_p4_kin.Boost(-(beam_p4_kin + target_p4).BoostVector())")
    .Define("kp_p_meas","kp_p4_meas.P()")
    .Define("kp_p_kin","kp_p4_kin.P()")
    .Define("km_p_meas","km_p4_meas.P()")
    .Define("km_p_kin","km_p4_kin.P()")
    .Define("d_p_meas","d_p4_meas.P()")
    .Define("d_p_kin","d_p4_kin.P()")
    .Define("phi_mass_meas","phi_p4_meas.M()")
    .Define("phi_mass_kin","phi_p4_kin.M()")
    .Define("phi_theta_meas","phi_p4_meas.Theta()*RadToDeg")
    .Define("phi_theta_kin","phi_p4_kin.Theta()*RadToDeg")
    // .Define("missd_mass_meas","missd_p4_meas.M()")
    // .Define("missd_mass_kin","missd_p4_kin.M()")
    .Define("sqrt_s_meas", "(beam_p4_meas + target_p4).Mag()")
    .Define("sqrt_s_kin", "(beam_p4_kin + target_p4).Mag()")
    .Define("minus_t_meas", "-(beam_p4_meas - phi_p4_meas).Mag2()")
    .Define("minus_t_kin", "-(beam_p4_kin - phi_p4_kin).Mag2()")
    // .Define("minus_u_meas", "-(beam_p4_meas - missd_p4_meas).Mag2()")
    // .Define("minus_u_kin", "-(beam_p4_kin - missd_p4_kin).Mag2()")
    // .Define("rho_mass_meas","pip_p4_meas + pim_p4_meas")
    // .Define("rho_mass_kin","pip_p4_kin + pim_p4_kin")
    .Define("y_phi_meas","minus_t_meas/(2*mass_deuteron*(beam_p4_meas.E()-phi_p4_meas.E()))")
    .Define("y_phi_kin","minus_t_kin/(2*mass_deuteron*(beam_p4_kin.E()-phi_p4_kin.E()))")
    // .Define("DeltaE_meas","(s_meas - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_meas)) - phi_p4com_meas.E()")
    // .Define("DeltaE_kin","(s_kin - pow(mass_deuteron, 2) + pow(mass_phi, 2)) / (2.0 * sqrt(s_kin)) - phi_p4com_kin.E()")
    ;

    // Filter events and save to new tree
    cout << "Filtering events...\n";
    auto rdf_no_filter = rdf_def;
    auto rdf_cl_filtered = rdf_no_filter.Filter([](double kin_cl) {return kin_cl > 0.01 ;}, {"kin_cl"});
    auto rdf_pidfom_filtered = rdf_cl_filtered.Filter([](double kp_pidfom, double km_pidfom) {return (kp_pidfom > 0.01) && (km_pidfom > 0.01);}, {"kp_pidfom","km_pidfom"});
    auto rdf_y_phi_filtered = rdf_pidfom_filtered.Filter([](double y_phi_meas) {return y_phi_meas > 0.4;}, {"y_phi_meas"});
    auto rdf_phi_mass_filtered = rdf_y_phi_filtered.Filter([](double phi_mass_meas) {return phi_mass_meas > 1.01 && phi_mass_meas < 1.03;}, {"phi_mass_meas"});

    rdf_phi_mass_filtered.Snapshot("filteredtree_phi_d_2H_data",OutputFile);

    // Plot histograms
    cout << "Plotting histograms...\n";
    TFile * histFile = new TFile(OutputFile.c_str(), "update");
    histFile->cd();
    vector<TH1*> hist_list;

    int N_filters = 5;
    RNode rdfs [] = {rdf_no_filter, rdf_cl_filtered, rdf_pidfom_filtered, rdf_y_phi_filtered, rdf_phi_mass_filtered};
    string labels [] = {"no_cut", "cl_cut", "pidfom_cut", "y_phi_cut", "phi_mass_cut"};


    for (int i = 0; i < N_filters; i++)
    {
        auto rdf = rdfs[i];
        string label = labels[i];

        TDirectory * dir = histFile->mkdir(label.c_str());
        dir->cd();

        TH2D hist_dEdx_kp = *rdf.Histo2D({("dEdx_kp_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"kp_p_meas","kp_dedx_fdc","accidweight");
        hist_dEdx_kp.Write();

        TH2D hist_dEdx_km = *rdf.Histo2D({("dEdx_km_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 10.0, 200, 0.0, 2e-5},"km_p_meas","km_dedx_fdc","accidweight");
        hist_dEdx_km.Write();

        TH2D hist_dEdx_d = *rdf.Histo2D({("dEdx_d_"+ label).c_str(), ";p (GeV/c);dE/dx (keV/cm)", 100, 0.0, 3.0, 200, 0.0, 2e-5},"d_p_meas","d_dedx_cdc","accidweight");
        hist_dEdx_d.Write();

        TH1D hist_massKK = *rdf.Histo1D({("massKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);Counts", 400, 0.9, 1.3},"phi_mass_meas","accidweight");
        hist_massKK.Write();

        // TH1D hist_massmiss = *rdf.Histo1D({("massmiss_"+ label).c_str(), ";m_{miss} (GeV/c^{2});Counts",200, 1.0, 3.0},"missd_mass_meas","accidweight");
        // hist_massmiss.Write();

        // TH2D hist_massKK_massmiss = *rdf.Histo2D({("massKK_massmiss_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);m_{miss} (GeV/c^{2})", 1000, 0.0, 2.0, 200, 1.0, 3.0},"phi_mass_meas","missd_mass_meas","accidweight");
        // hist_massKK_massmiss.Write();

        TH2D hist_massKK_yphi = *rdf.Histo2D({("massKK_yphi_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);y_{#phi}", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","y_phi_meas","accidweight");
        hist_massKK_yphi.Write();

        TH2D hist_massKK_thetaKK = *rdf.Histo2D({("massKK_thetaKK_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);#theta_{K^{+}K^{-}} (deg)", 400, 0.9, 1.3, 200, 0.0, 20.0},"phi_mass_meas","phi_theta_meas","accidweight");
        hist_massKK_thetaKK.Write();

        TH2D hist_massKK_minust = *rdf.Histo2D({("massKK_minust_"+ label).c_str(), ";m_{K^{+}K^{-}} (GeV/c);-t (GeV^{2}/c^{2})", 400, 0.9, 1.3, 200, 0.0, 2.0},"phi_mass_meas","minus_t_meas","accidweight");
        hist_massKK_minust.Write();
    }

    histFile->Close();
}


    // TH2F *hist_MPhi_DeltaE          = new TH2F("hist_MPhi_DeltaE", "hist_MPhi_DeltaE", 400, 0.9, 1.3, 200, -2.0, 2.0);
    // TH2F *hist_MPhi_MRho            = new TH2F("hist_MPhi_MRho", "hist_MPhi_MRho", 400, 0.9, 1.3, 800, 0.3, 1.1);
    // TH2F *hist_ThetaPhi_yPhi        = new TH2F("hist_ThetaPhi_yPhi", "hist_ThetaPhi_yPhi", 200, 0.0, 20.0, 200, 0.0, 2.0);
    // TH2F *hist_MissingMass_yPhi     = new TH2F("hist_MissingMass_yPhi", "hist_MissingMass_yPhi", 200, 1.0, 3.0, 200, 0.0, 2.0);
    // TH2F *hist_DeltaE_yPhi          = new TH2F("hist_DeltaE_yPhi", "hist_DeltaE_yPhi", 200, -2.0, 2.0, 200, 0.0, 2.0);
    // TH2F *hist_DeltaE_MissingMass   = new TH2F("hist_DeltaE_MissingMass", "hist_DeltaE_MissingMass", 200, -2.0, 2.0, 200, 1.0, 3.0);

    //     hist_MPhi_DeltaE->Fill(phi_mass, phi_energy_expected - PhiP4COM->E(), WeightFactor);
    //     hist_MPhi_MRho->Fill(phi_mass, rho_mass, WeightFactor);
    //     hist_ThetaPhi_yPhi->Fill(phi_theta, y_phi, WeightFactor);
    //     hist_MissingMass_yPhi->Fill(MissingP4->M(), y_phi, WeightFactor);
    //     hist_DeltaE_yPhi->Fill(phi_energy_expected - PhiP4COM->E(), y_phi, WeightFactor);
    //     hist_DeltaE_MissingMass->Fill(phi_energy_expected - PhiP4COM->E(), MissingP4->M(), WeightFactor);












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

    //     TLorentzVector *PhiP4COM = new TLorentzVector(*KPlusP4 + *KMinusP4);
    //     PhiP4COM->Boost(-(*BeamP4+*TargetP4).BoostVector());
    //     TLorentzVector *KPlusP4AsPion = new TLorentzVector(KPlusP4->X(), KPlusP4->Y(), KPlusP4->Z(), sqrt(KPlusP4->P()*KPlusP4->P() + mass_pion*mass_pion));
    //     TLorentzVector *KMinusP4AsPion = new TLorentzVector(KMinusP4->X(), KMinusP4->Y(), KMinusP4->Z(), sqrt(KMinusP4->P()*KMinusP4->P() + mass_pion*mass_pion));

    //     double phi_energy_expected  = (pow(sqrt_s, 2) - pow(TargetP4->M(), 2) + pow(mass_phi, 2)) / (2. * sqrt_s);
    //     double rho_mass             = (*KPlusP4AsPion + *KMinusP4AsPion).M();


    // }

