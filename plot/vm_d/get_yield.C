#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;
using namespace ROOT;
using namespace ROOT::RDF;
using namespace ROOT::Detail::RDF;

double mass_kaon = 0.493677;

// Nominal fitting function
Double_t rel_bw_plus_linear             (Double_t *x, Double_t *par);
Double_t rel_bw                         (Double_t *x, Double_t *par);
Double_t linear                         (Double_t *x, Double_t *par);

// Alternative background functions
Double_t rel_bw_plus_fulllinear         (Double_t *x, Double_t *par);
Double_t rel_bw_plus_quadratic          (Double_t *x, Double_t *par);
Double_t rel_bw_plus_fullquadratic      (Double_t *x, Double_t *par);
Double_t rel_bw_plus_phenomenological   (Double_t *x, Double_t *par);
Double_t fulllinear                     (Double_t *x, Double_t *par);
Double_t quadratic                      (Double_t *x, Double_t *par);
Double_t fullquadratic                  (Double_t *x, Double_t *par);
Double_t phenomenological               (Double_t *x, Double_t *par);

// Alternative signal functions
Double_t rel_bw_noBL_plus_linear        (Double_t *x, Double_t *par);
Double_t nonrel_bw_plus_linear          (Double_t *x, Double_t *par);
Double_t rel_bw_noBL                    (Double_t *x, Double_t *par);
Double_t nonrel_bw                      (Double_t *x, Double_t *par);

// Simulation weight functions
double sim_weight_func_nominal      (double beam_energy_truth, double minust_truth);
double sim_weight_func_iterations   (double beam_energy_truth, double minust_truth, int iteration);
double sim_weight_func_systematic   (double beam_energy_truth, double minust_truth, int a1_variation, int b1_variation, int a2_variation, int b2_variation);

int get_yield(string channel, string reaction, string observable, string tag)
{
    // Read tree from input root file
    string input_treefile_name  = Form("/work/halld2/home/boyu/src_analysis/filter/output/filteredtree_%s_%s.root", channel.c_str(), reaction.c_str());
    string input_tree_name;
    if (reaction.find("recon") != string::npos)
        input_tree_name = Form("filteredtree_%s_recon", channel.c_str());
    else if (reaction.find("thrown") != string::npos)
        input_tree_name = Form("filteredtree_%s_thrown", channel.c_str());
    cout << "Input tree file: " << input_treefile_name << endl;
    cout << "Input tree: " << input_tree_name << endl;
    TChain chain(input_tree_name.c_str());
    chain.Add(input_treefile_name.c_str());
    RDataFrame rdf_raw(chain);
    auto rdf_input = RNode(rdf_raw);
    if (reaction.find("recon") != string::npos)
    {
        string dEdxCut          = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.65*d_momentum_meas+4.47) + 2.57)";
        string MissPMinusCut    = "miss_pminus_meas > -0.02";
        string KinFitChiSqCut   = "chisq_per_ndf_kin < 5.0";
        string MomentumCut      = "kp_momentum_meas > 0.400 && km_momentum_meas > 0.400 && d_momentum_meas > 0.400";
        string ThetaCut         = "kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        string VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 14.0";
        string VertexRCut       = "TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        string BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 18.0";

        if      (tag.find("dEdx_1.00") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-4.01*d_momentum_meas+4.88) + 3.26)";
        else if (tag.find("dEdx_1.50") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.85*d_momentum_meas+4.69) + 2.92)";
        else if (tag.find("dEdx_1.75") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.74*d_momentum_meas+4.58) + 2.73)";
        else if (tag.find("dEdx_2.50") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.41*d_momentum_meas+4.21) + 2.21)";
        else if (tag.find("dEdx_3.00") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.11*d_momentum_meas+3.90) + 1.83)";
        else if (tag.find("misspminus_0.0100") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.0100";
        else if (tag.find("misspminus_0.0150") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.0150";
        else if (tag.find("misspminus_0.0175") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.0175";
        else if (tag.find("misspminus_0.0250") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.0250";
        else if (tag.find("misspminus_0.0300") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.0300";
        else if (tag.find("chisquared_3.00") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 3.00";
        else if (tag.find("chisquared_3.50") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 3.50";
        else if (tag.find("chisquared_4.00") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 4.00";
        else if (tag.find("chisquared_4.50") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 4.50";
        else if (tag.find("chisquared_4.75") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 4.75";
        else if (tag.find("chisquared_5.50") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 5.50";
        else if (tag.find("chisquared_6.00") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 6.00";
        else if (tag.find("chisquared_7.00") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 7.00";
        else if (tag.find("momentum_0.350") != string::npos)
            MomentumCut    = "kp_momentum_meas > 0.350 && km_momentum_meas > 0.350 && d_momentum_meas > 0.350";
        else if (tag.find("momentum_0.375") != string::npos)
            MomentumCut    = "kp_momentum_meas > 0.375 && km_momentum_meas > 0.375 && d_momentum_meas > 0.375";
        else if (tag.find("momentum_0.425") != string::npos)
            MomentumCut    = "kp_momentum_meas > 0.425 && km_momentum_meas > 0.425 && d_momentum_meas > 0.425";
        else if (tag.find("momentum_0.450") != string::npos)
            MomentumCut    = "kp_momentum_meas > 0.450 && km_momentum_meas > 0.450 && d_momentum_meas > 0.450";
        else if (tag.find("theta_1.80") != string::npos)
            ThetaCut    = "kp_theta_meas > 1.80 && km_theta_meas > 1.80 && d_theta_meas > 1.80";
        else if (tag.find("theta_1.90") != string::npos)
            ThetaCut    = "kp_theta_meas > 1.90 && km_theta_meas > 1.90 && d_theta_meas > 1.90";
        else if (tag.find("theta_1.95") != string::npos)
            ThetaCut    = "kp_theta_meas > 1.95 && km_theta_meas > 1.95 && d_theta_meas > 1.95";
        else if (tag.find("theta_2.05") != string::npos)
            ThetaCut    = "kp_theta_meas > 2.05 && km_theta_meas > 2.05 && d_theta_meas > 2.05";
        else if (tag.find("theta_2.10") != string::npos)
            ThetaCut    = "kp_theta_meas > 2.10 && km_theta_meas > 2.10 && d_theta_meas > 2.10";
        else if (tag.find("theta_2.20") != string::npos)
            ThetaCut    = "kp_theta_meas > 2.20 && km_theta_meas > 2.20 && d_theta_meas > 2.20";
        else if (tag.find("vertexZ_13.00") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 13.00";
        else if (tag.find("vertexZ_13.50") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 13.50";
        else if (tag.find("vertexZ_13.75") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 13.75";
        else if (tag.find("vertexZ_14.25") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 14.25";
        else if (tag.find("vertexZ_14.50") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 14.50";
        else if (tag.find("vertexZ_15.00") != string::npos)
            VertexZCut       = "TMath::Abs(vertex_z_kin - 65.0) < 15.00";
        else if (tag.find("vertexR_0.50") != string::npos)
            VertexRCut       = "TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 0.50";
        else if (tag.find("vertexR_0.75") != string::npos)
            VertexRCut       = "TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 0.75";
        else if (tag.find("vertexR_1.25") != string::npos)
            VertexRCut       = "TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.25";
        else if (tag.find("vertexR_1.50") != string::npos)
            VertexRCut       = "TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.50";
        else if (tag.find("beamaccid_5") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 22.0";
        else if (tag.find("beamaccid_3") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 14.0";
        else if (tag.find("beamaccid_4out") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 2.0 || (TMath::Abs(beam_DeltaT_meas) > 6.0 && TMath::Abs(beam_DeltaT_meas) < 22.0)";

        rdf_input = rdf_input   .Filter(dEdxCut.c_str())
                                .Filter(MissPMinusCut.c_str())
                                .Filter(KinFitChiSqCut.c_str())
                                .Filter(MomentumCut.c_str())
                                .Filter(ThetaCut.c_str())
                                .Filter(VertexZCut.c_str())
                                .Filter(VertexRCut.c_str())
                                .Filter(BeamAccidCut.c_str());

        if      (tag.find("beamaccid_3") != string::npos)
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/6.0");
        else if (tag.find("beamaccid_5") != string::npos)
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/10.0");
        else
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/8.0");

        if      (tag.find("comboaccid_all") != string::npos)
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "NumCombos == 1.0 ? 1.0 : 1/NumCombos");
        else if (tag.find("comboaccid_none") != string::npos)
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "NumCombos == 1.0 ? 1.0 : 0.0");
        else
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "combo_accid_weight");

        if (reaction.find("model") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");
        else if (tag.find("simweight_iter0") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 0)");
        else if (tag.find("simweight_iter1") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 1)");
        else if (tag.find("simweight_iter2") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 2)");
        else if (tag.find("simweight_iter3") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 3)");
        else if (tag.find("simweight_iter4") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 4)");
        else if (tag.find("simweight_iter5") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 5)");
        else if (tag.find("simweight_iter6") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 6)");
        else if (tag.find("simweight_iter7") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 7)");
        else if (tag.find("simweight_syst") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_systematic(beam_energy_truth, minust_truth, tag)");
        else
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_nominal(beam_energy_truth, minust_truth)");

        rdf_input = rdf_input       .Define("yield_weight",             "beamaccid_weight_syst*combo_accid_weight_syst*sim_weight_syst")
                                    .Define("yield_weight_squared",     "yield_weight*yield_weight");
    }
    else if (reaction.find("thrown") != string::npos)
    {
        if (reaction.find("model") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");
        else if (tag.find("simweight_iter0") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_iterations(beam_energy_truth, minust_truth, 0)");
        else if (tag.find("simweight_syst") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_systematic(beam_energy_truth, minust_truth, tag)");

        else
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_nominal(beam_energy_truth, minust_truth)");

        rdf_input = rdf_input       .Define("yield_weight",             "sim_weight_syst")
                                    .Define("yield_weight_squared",     "yield_weight*yield_weight");
    }

    // Read bin edges from input text file
    string input_txt_name = Form("/work/halld2/home/boyu/src_analysis/plot/vm_d/configs/bins_%s_%s.txt", channel.c_str(), observable.c_str());
    cout << "Input text file: " << input_txt_name << endl;
    ifstream input_txt(input_txt_name.c_str());
    vector<vector<double>> bins;
    string this_line;
    while (getline(input_txt, this_line))
    {
        istringstream iss(this_line);
        vector<double> this_row;
        double this_value;
        while (iss >> this_value)
            this_row.push_back(this_value);
        bins.push_back(this_row);
    }
    input_txt.close();

    // Prepare output files
    string output_textfile_name = Form("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/yield_%s/yield_%s_%s_%s_%s.txt", channel.c_str(), channel.c_str(), reaction.c_str(), observable.c_str(), tag.c_str());
    string output_pdffile_name = Form("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/yield_%s/yield_%s_%s_%s_%s.pdf", channel.c_str(), channel.c_str(), reaction.c_str(), observable.c_str(), tag.c_str());
    cout << "Output text file: " << output_textfile_name << endl;
    cout << "Output PDF file: " << output_pdffile_name << endl;
    FILE *output_textfile = fopen(output_textfile_name.c_str(),"w");

    // Loop over bins and calculate yield
    cout << "Calculating yield..." << endl;
    TCanvas *canvas = new TCanvas("c1", "c1", 800, 600);
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    double yield, yield_err, energy_center, energy_width, t_center, t_width, angle_center, angle_width;
    string energy_cut, t_cut, angle_cut;
    string variable;

    double  fit_min     = 0.99;
    double  fit_max     = 1.08;
    double  bin_width   = 0.005;
    double  hist_min    = 0.9825;
    double  hist_max    = 1.1125;

    if (tag == "fitmax_1.06")
        fit_max = 1.06;
    else if (tag == "fitmax_1.07")
        fit_max = 1.07;
    else if (tag == "fitmax_1.09")
        fit_max = 1.09;
    else if (tag == "fitmax_1.10")
        fit_max = 1.10;

    if (tag == "fitwidth_0.0040")
        bin_width = 0.0040;
    else if (tag == "fitwidth_0.0048")
        bin_width = 0.0048;
    else if (tag == "fitwidth_0.0060")
        bin_width = 0.0060;
    else if (tag == "fitwidth_0.0075")
        bin_width = 0.0075;
    int num_bins    = (hist_max-hist_min)/bin_width;

    for (int i = 0; i < bins.size(); i++)
    {
        energy_center   = (bins[i][1] + bins[i][0]) / 2.0;
        energy_width    = (bins[i][1] - bins[i][0]) / 2.0;
        t_center        = (bins[i][3] + bins[i][2]) / 2.0;
        t_width         = (bins[i][3] - bins[i][2]) / 2.0;
        if (observable != "dsdt")
        {
            angle_center = (bins[i][5] + bins[i][4]) / 2.0;
            angle_width  = (bins[i][5] - bins[i][4]) / 2.0;
        }

        auto rdf_bin = RNode(rdf_input);
        TH1D hist_bin;

        if (reaction.find("recon") != string::npos)
        {
            energy_cut  = Form("beam_energy_kin>%.2f && beam_energy_kin<%.2f", bins[i][0], bins[i][1]);
            t_cut       = Form("minust_kin>%.3f && minust_kin<%.3f", bins[i][2], bins[i][3]);
            rdf_bin = rdf_input.Filter(energy_cut.c_str()).Filter(t_cut.c_str());
            if (reaction.find("data") != string::npos)
            {
                energy_center   = rdf_bin.Mean("beam_energy_kin").GetValue();
                energy_width    = rdf_bin.StdDev("beam_energy_kin").GetValue();
                t_center        = rdf_bin.Mean("minust_kin").GetValue();
                t_width         = rdf_bin.StdDev("minust_kin").GetValue();
            }

            if (observable == "dsdt")
            {
                cout << energy_cut << " && " << t_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]), ";m_{K^{+}K^{-}} (GeV/c);Counts", num_bins, hist_min, hist_max},"phi_mass_kin","yield_weight");
            }
            else
            {
                if (observable == "Wcostheta")
                    variable = "decay_costheta_helicity_kin";
                else if (observable == "Wdecayphi")
                    variable = "decay_phi_helicity_kin";
                else if (observable == "Wpolphi")
                    variable = "polarization_phi_com_kin";
                else if (observable == "Wpsi")
                    variable = "psi_helicity_kin";
                angle_cut = Form("%s>%.2f && %s<%.2f", variable.c_str(), bins[i][4], variable.c_str(), bins[i][5]);
                rdf_bin = rdf_bin.Filter(angle_cut.c_str());
                if (tag == "sideband")
                {
                    rdf_bin = rdf_bin.Filter("phi_mass_kin>1.05");
                }
                if (reaction.find("data") != string::npos)
                {
                    angle_center = rdf_bin.Mean(variable.c_str()).GetValue();
                    angle_width  = rdf_bin.StdDev(variable.c_str()).GetValue();
                }
                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f_%.2f_%.2f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]), ";m_{K^{+}K^{-}} (GeV/c);Counts", num_bins, hist_min, hist_max},"phi_mass_kin","yield_weight");
            }

            for (int i = 0; i < hist_bin.GetNbinsX(); i++)
            {
                if (hist_bin.GetBinContent(i) <= 0)
                {
                    // hist_bin.SetBinContent(i, 0);  hist_bin.SetBinError(i, 1);  // set non-positive bin content to zero, with error of 1
                    hist_bin.SetBinError(i, sqrt(hist_bin.GetBinError(i)*hist_bin.GetBinError(i) + 1.0));  // keep the negative bin content, but add error of 1 by quadrature
                }
            }

            if ((reaction.find("sim") != string::npos && tag != "fitsig_relBWsim") || tag == "sideband")
            {
                yield = rdf_bin.Sum("yield_weight").GetValue();
                yield_err = sqrt(rdf_bin.Sum("yield_weight_squared").GetValue());
                hist_bin.Draw();
            }
            else if (reaction.find("sim") != string::npos && tag == "fitsig_relBWsim")
            {
                TF1 fit_func("fit_func", rel_bw, 0.99, fit_max, 4);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitbkg_fulllinear")
            {
                TF1 fit_func("fit_func", rel_bw_plus_fulllinear, 0.99, fit_max, 6);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 0.0, 0.0);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.5, 200.0);
                fit_func.SetParLimits(5, -10.0, 10.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *fulllinear_function = new TF1("fulllinear_function", fulllinear, 0.99, fit_max, 2);
                fulllinear_function->SetParameter(0, fit_func.GetParameter(4));
                fulllinear_function->SetParameter(1, fit_func.GetParameter(5));
                fulllinear_function->SetLineColor(kBlue);
                fulllinear_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitbkg_quadratic")
            {
                TF1 fit_func("fit_func", rel_bw_plus_quadratic, 0.99, fit_max, 6);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20, 0.0);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.01, 100.0);
                fit_func.SetParLimits(5, 0.0, 2*mass_kaon);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *quadratic_function = new TF1("quadratic_function", quadratic, 0.99, fit_max, 2);
                quadratic_function->SetParameter(0, fit_func.GetParameter(4));
                quadratic_function->SetParameter(1, fit_func.GetParameter(5));
                quadratic_function->SetLineColor(kBlue);
                quadratic_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitbkg_fullquadratic")
            {
                TF1 fit_func("fit_func", rel_bw_plus_fullquadratic, 0.99, fit_max, 7);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20, 0.0, 0.0);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.01, 100.0);
                fit_func.SetParLimits(5, 0.0, 2*mass_kaon);
                fit_func.SetParLimits(6, -10.0, 10.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *fullquadratic_function = new TF1("fullquadratic_function", fullquadratic, 0.99, fit_max, 3);
                fullquadratic_function->SetParameter(0, fit_func.GetParameter(4));
                fullquadratic_function->SetParameter(1, fit_func.GetParameter(5));
                fullquadratic_function->SetParameter(2, fit_func.GetParameter(6));
                fullquadratic_function->SetLineColor(kBlue);
                fullquadratic_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitbkg_phenomenological")
            {
                TF1 fit_func("fit_func", rel_bw_plus_phenomenological, 0.99, fit_max, 6);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.01, 100.0);
                fit_func.SetParLimits(5, 0.01, 100.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *phenomenological_function = new TF1("phenomenological_function", phenomenological, 0.99, fit_max, 2);
                phenomenological_function->SetParameter(0, fit_func.GetParameter(4));
                phenomenological_function->SetParameter(1, fit_func.GetParameter(5));
                phenomenological_function->SetLineColor(kBlue);
                phenomenological_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitsig_noBL")
            {
                TF1 fit_func("fit_func", rel_bw_noBL_plus_linear, 0.99, fit_max, 5);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.5, 200.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_noBL_func = new TF1("rel_bw_noBL_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_noBL_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_noBL_func->SetLineColor(kRed);
                rel_bw_noBL_func->Draw("same");
                TF1 *linear_function = new TF1("linear_function", linear, 0.99, fit_max, 1);
                linear_function->SetParameter(0, fit_func.GetParameter(4));
                linear_function->SetLineColor(kBlue);
                linear_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else if (tag == "fitsig_nonrel")
            {
                TF1 fit_func("fit_func", nonrel_bw_plus_linear, 0.99, fit_max, 5);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.5, 200.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *nonrel_bw_func = new TF1("nonrel_bw_func", nonrel_bw, 0.99, fit_max, 4);
                nonrel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                nonrel_bw_func->SetLineColor(kRed);
                nonrel_bw_func->Draw("same");
                TF1 *linear_function = new TF1("linear_function", linear, 0.99, fit_max, 1);
                linear_function->SetParameter(0, fit_func.GetParameter(4));
                linear_function->SetLineColor(kBlue);
                linear_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            else
            {
                TF1 fit_func("fit_func", rel_bw_plus_linear, 0.99, fit_max, 5);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                fit_func.SetParLimits(4, 0.5, 200.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                // yield = 1;
                // yield_err = 1;
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *linear_function = new TF1("linear_function", linear, 0.99, fit_max, 1);
                linear_function->SetParameter(0, fit_func.GetParameter(4));
                linear_function->SetLineColor(kBlue);
                linear_function->Draw("same");
                yield = fit_func.GetParameter(0)/bin_width;
                yield_err = fit_func.GetParError(0)/bin_width;
            }
            canvas->Update();
            canvas->Print((output_pdffile_name+"(").c_str());
            canvas->Clear();
        }
        else if (reaction.find("thrown") != string::npos)
        {
            energy_cut  = Form("beam_energy_truth>%.2f && beam_energy_truth<%.2f", bins[i][0], bins[i][1]);
            t_cut       = Form("minust_truth>%.3f && minust_truth<%.3f", bins[i][2], bins[i][3]);
            rdf_bin     = rdf_bin.Filter(energy_cut.c_str()).Filter(t_cut.c_str());

            if (observable == "dsdt")
            {
                cout << energy_cut << " && " << t_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]), ";m_{K^{+}K^{-}} (GeV/c);Counts", num_bins, hist_min, hist_max},"phi_mass_truth","yield_weight");
            }
            else
            {
                if (observable == "Wcostheta")
                    variable = "decay_costheta_helicity_truth";
                else if (observable == "Wdecayphi")
                    variable = "decay_phi_helicity_truth";
                else if (observable == "Wpolphi")
                    variable = "polarization_phi_com_truth";
                else if (observable == "Wpsi")
                    variable = "psi_helicity_truth";
                angle_cut = Form("%s>%.2f && %s<%.2f", variable.c_str(), bins[i][4], variable.c_str(), bins[i][5]);
                rdf_bin = rdf_bin.Filter(angle_cut.c_str());
                if (tag == "sideband")
                {
                    rdf_bin = rdf_bin.Filter("phi_mass_kin>1.05");
                }

                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f_%.2f_%.2f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]), ";m_{K^{+}K^{-}} (GeV/c);Counts", num_bins, hist_min, hist_max},"phi_mass_truth","yield_weight");
            }

            yield           = rdf_bin.Sum("yield_weight").GetValue();
            yield_err       = sqrt(rdf_bin.Sum("yield_weight_squared").GetValue());
            hist_bin.Draw();
            canvas->Update();
            canvas->Print((output_pdffile_name+"(").c_str());
            canvas->Clear();
        }

        if (observable == "dsdt")
            fprintf(output_textfile, "%6.4f\t%6.4f\t%6.1f\t%6.1f\t%6.4f\t%6.4f\t%6.3f\t%6.3f\t%f\t%f\n", energy_center, energy_width, bins[i][0], bins[i][1], t_center, t_width, bins[i][2], bins[i][3], yield, yield_err);
            // fprintf(output_textfile, "%6.1f\t%6.1f\t%6.3f\t%6.3f\t%f\t%f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3], yield, yield_err);
            // fprintf(output_textfile, "%6.1f\t%6.1f\t%6.3f\t%6.3f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3]);
        else
            fprintf(output_textfile, "%6.3f\t%6.3f\t%6.1f\t%6.1f\t%6.3f\t%6.4f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.1f\t%6.1f\t%f\t%f\n", energy_center, energy_width, bins[i][0], bins[i][1], t_center, t_width, bins[i][2], bins[i][3], angle_center, angle_width, bins[i][4], bins[i][5], yield, yield_err);
            // fprintf(output_textfile, "%6.1f\t%6.1f\t%6.3f\t%6.3f\t%6.1f\t%6.1f\t%f\t%f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5], yield, yield_err);
            // fprintf(output_textfile, "%6.1f\t%6.1f\t%6.3f\t%6.3f\t%6.1f\t%6.1f\n", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]);
    }
    canvas->Print((output_pdffile_name+")").c_str());

    return 0;

}

Double_t rel_bw_plus_linear(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a linear background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Double_t slope = par[4];
    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
}

Double_t rel_bw(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum;
}

Double_t linear(Double_t *x, Double_t *par)
{
    // par[0] = linear slope
    return par[0] * (x[0] - 2*mass_kaon);
}

Double_t rel_bw_plus_fulllinear(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a linear background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Double_t slope = par[4];
    Double_t intercept = par[5];
    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon) + intercept;
}

Double_t rel_bw_plus_quadratic(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a quadratic background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope, par[5] = quadratic coefficient

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];

    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + par[4] * (x[0] - 2*mass_kaon) * (x[0] - par[5]);
}

Double_t rel_bw_plus_fullquadratic(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a quadratic background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope, par[5] = quadratic coefficient

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];

    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + par[4] * (x[0] - 2*mass_kaon) * (x[0] - par[5]) + par[6];
}

Double_t rel_bw_plus_phenomenological(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a quadratic background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = nonlinear coefficient, par[5] = quadratic coefficient

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];

    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Blatt-Weisskopf barrier factor squared
        Double_t R0 = 5.0677; // interaction radius in GeV^-1
        Double_t z = q * R0;
        Double_t z0 = q0 * R0;
        Double_t BL = (1 + z0*z0) / (1 + z*z);

        // Mass-dependent width
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime) * BL;

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + par[4] * sqrt(x[0]*x[0] - 4*mass_kaon*mass_kaon) + par[5] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
}

Double_t fulllinear(Double_t *x, Double_t *par)
{
    // par[0] = linear slope
    // par[1] = intercept at 2k mass
    return par[0] * (x[0] - 2*mass_kaon) + par[1];
}

Double_t quadratic(Double_t *x, Double_t *par)
{
    // par[0] = quadratic coefficient
    // par[1] = the other zero point
    return par[0] * (x[0] - 2*mass_kaon) * (x[0] - par[1]);
}

Double_t fullquadratic(Double_t *x, Double_t *par)
{
    // par[0] = quadratic coefficient
    // par[1] = the other zero point
    // par[2] = intercept at 2k mass
    return par[0] * (x[0] - 2*mass_kaon) * (x[0] - par[1]) + par[2];
}

Double_t phenomenological(Double_t *x, Double_t *par)
{
    // par[0] = nonlinear coefficient
    // par[1] = quadratic coefficient
    double mass_kaon = 0.493677;
    Double_t quadratic;
    if (x[0] < 2*mass_kaon)
        quadratic = 0.0;
    else
        quadratic = par[0] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
    return quadratic;
}

Double_t rel_bw_noBL_plus_linear(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner without Blatt-Weisskopf factor, plus a linear background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Double_t slope = par[4];
    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Mass-dependent width without Blatt-Weisskopf factor
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime);

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
}

Double_t rel_bw_noBL(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner without Blatt-Weisskopf factor
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Int_t L = 1;
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Kaon mass and momentum calculations
        Double_t mk = mass_kaon;
        Double_t q = sqrt((xprime*xprime - 4*mk*mk)/4); // momentum of kaon in phi rest frame
        Double_t q0 = sqrt((M0*M0 - 4*mk*mk)/4);     // momentum at pole mass

        // Mass-dependent width without Blatt-Weisskopf factor
        Double_t Gamma = Gamma0 * (q/q0) * (q/q0) * (q/q0) * (M0/xprime);

        // Relativistic Breit-Wigner
        Double_t denominator = (xprime*xprime - M0*M0)*(xprime*xprime - M0*M0) + M0*M0*Gamma*Gamma;

        convol_sum += (2 / TMath::Pi()) * xprime * M0 * Gamma / denominator * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] * convol_sum;
}

Double_t nonrel_bw_plus_linear(Double_t *x, Double_t *par)
{
    // Non-relativistic Breit-Wigner plus a linear background function
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width, par[4] = linear slope

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Double_t slope = par[4];
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Non-relativistic Breit-Wigner
        Double_t denominator = (xprime - M0)*(xprime - M0) + (Gamma0*Gamma0)/4;

        convol_sum += (1 / TMath::Pi()) * (Gamma0 / denominator) * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] / 2 * convol_sum + slope * (x[0] - 2*mass_kaon);
}

Double_t nonrel_bw(Double_t *x, Double_t *par)
{
    // Non-relativistic Breit-Wigner
    // par[0] = BW amplitude, par[1] = BW pole mass, par[2] = BW pole width, par[3] = Gaussian width

    Double_t M0 = par[1];
    Double_t Gamma0 = par[2];
    Double_t sigma = par[3];
    Double_t convol_sum = 0.0;
    Double_t convol_range = 10.0;
    Double_t convol_step = 0.001;

    for (double xprime = x[0]-convol_range; xprime < x[0]+convol_range; xprime += convol_step)
    {
        if (xprime < 2*mass_kaon)
            continue;

        // Non-relativistic Breit-Wigner
        Double_t denominator = (xprime - M0)*(xprime - M0) + (Gamma0*Gamma0)/4;

        convol_sum += (1 / TMath::Pi()) * (Gamma0 / denominator) * TMath::Gaus(x[0]-xprime, 0, sigma, true) * convol_step;
    }

    return par[0] / 2 * convol_sum;
}

// double sim_weight_func_pass1(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 3632.01;
//         b1 = 15.82;
//         a2 = 12.19;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 4842.04;
//         b1 = 17.22;
//         a2 = 17.43;
//         b2 = 3.13;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 2882.53;
//         b1 = 15.54;
//         a2 = 12.54;
//         b2 = 2.90;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass2(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6127.78;
//         b1 = 17.10;
//         a2 = 12.43;
//         b2 = 2.50;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 11328.73;
//         b1 = 19.75;
//         a2 = 23.50;
//         b2 = 3.54;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4602.99;
//         b1 = 16.69;
//         a2 = 13.21;
//         b2 = 2.95;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass3(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6417.46;
//         b1 = 17.21;
//         a2 = 12.42;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 12903.60;
//         b1 = 20.16;
//         a2 = 24.55;
//         b2 = 3.60;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4769.82;
//         b1 = 16.77;
//         a2 = 13.19;
//         b2 = 2.94;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass4(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6442.44;
//         b1 = 17.22;
//         a2 = 12.42;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 13181.89;
//         b1 = 20.22;
//         a2 = 24.74;
//         b2 = 3.61;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4781.31;
//         b1 = 16.78;
//         a2 = 13.18;
//         b2 = 2.94;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass5(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6444.93;
//         b1 = 17.22;
//         a2 = 12.42;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 13225.39;
//         b1 = 20.23;
//         a2 = 24.77;
//         b2 = 3.61;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4782.23;
//         b1 = 16.78;
//         a2 = 13.18;
//         b2 = 2.94;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass6(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6444.97;
//         b1 = 17.22;
//         a2 = 12.42;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 13232.23;
//         b1 = 20.23;
//         a2 = 24.78;
//         b2 = 3.61;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4782.18;
//         b1 = 16.78;
//         a2 = 13.18;
//         b2 = 2.94;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

// double sim_weight_func_pass7(double beam_energy_truth, double minust_truth)
// {
//     double a1, b1, a2, b2 = 0;
//     double normalization = 10;
//     if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
//         return 1.0;
//     else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
//     {
//         a1 = 6445.03;
//         b1 = 17.22;
//         a2 = 12.42;
//         b2 = 2.49;
//     }
//     else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
//     {
//         a1 = 13232.54;
//         b1 = 20.23;
//         a2 = 24.78;
//         b2 = 3.61;
//     }
//     else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
//     {
//         a1 = 4782.12;
//         b1 = 16.78;
//         a2 = 13.18;
//         b2 = 2.94;
//     }
//     return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
// }

double sim_weight_func_nominal(double beam_energy_truth, double minust_truth)
{
    // No energy dependence
    double a1, b1, a2, b2 = 0;
    double normalization = 10;
    if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
        return 1.0;
    else
    {
        a1 = 13232.54;
        b1 = 20.23;
        a2 = 24.78;
        b2 = 3.61;
    }
    return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;

    // With energy dependence
    // double a1, b1, a2, b2 = 0;
    // double normalization = 10;
    // if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
    //     return 1.0;
    // else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
    // {
    //     a1 = 6445.03;  b1 = 17.22; a2 = 12.42; b2 = 2.49;
    // }
    // else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
    // {
    //     a1 = 13232.54; b1 = 20.23; a2 = 24.78; b2 = 3.61;
    // }
    // else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
    // {
    //     a1 = 4782.12;  b1 = 16.78; a2 = 13.18; b2 = 2.94;
    // }
    // return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}

double sim_weight_func_iterations(double beam_energy_truth, double minust_truth, int iteration)
{
    // No energy dependence
    double a1, b1, a2, b2 = 0;
    double normalization = 10;
    if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
        return 1.0;
    else if (iteration == 0)
        return 1.0;
    else if (iteration == 1)
    {
        a1 = 6066;  b1 = 17.94; a2 = 18.72; b2 = 3.234;
    }



    else if (iteration == 2)
    {
        a1 = 11328.73; b1 = 19.75; a2 = 23.50; b2 = 3.54;
    }
    else if (iteration == 3)
    {
        a1 = 12903.60; b1 = 20.16; a2 = 24.55; b2 = 3.60;
    }
    else if (iteration == 4)
    {
        a1 = 13181.89; b1 = 20.22; a2 = 24.74; b2 = 3.61;
    }
    else if (iteration == 5)
    {
        a1 = 13225.39; b1 = 20.23; a2 = 24.77; b2 = 3.61;
    }
    else if (iteration == 6)
    {
        a1 = 13232.23; b1 = 20.23; a2 = 24.78; b2 = 3.61;
    }
    else if (iteration == 7)
    {
        a1 = 13232.54; b1 = 20.23; a2 = 24.78; b2 = 3.61;
    }
    return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;

    // With energy dependence
    // double a1, b1, a2, b2 = 0;
    // double normalization = 10;
    // if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
    //     return 1.0;
    // else if (iteration == "simweight_iter0")
    //     return 1.0;
    // else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
    // {
    //     if (iteration == "simweight_iter1")
    //     {
    //         a1 = 3632.01; b1 = 15.82; a2 = 12.19; b2 = 2.49;
    //     }
    //     else if (iteration == "simweight_iter2")
    //     {
    //         a1 = 6127.78; b1 = 17.10; a2 = 12.43; b2 = 2.50;
    //     }
    //     else if (iteration == "simweight_iter3")
    //     {
    //         a1 = 6417.46; b1 = 17.21; a2 = 12.42; b2 = 2.49;
    //     }
    //     else if (iteration == "simweight_iter4")
    //     {
    //         a1 = 6442.44; b1 = 17.22; a2 = 12.42; b2 = 2.49;
    //     }
    //     else if (iteration == "simweight_iter5")
    //     {
    //         a1 = 6444.93; b1 = 17.22; a2 = 12.42; b2 = 2.49;
    //     }
    //     else if (iteration == "simweight_iter6")
    //     {
    //         a1 = 6444.97; b1 = 17.22; a2 = 12.42; b2 = 2.49;
    //     }
    //     else if (iteration == "simweight_iter7")
    //     {
    //         a1 = 6445.03; b1 = 17.22; a2 = 12.42; b2 = 2.49;
    //     }
    // }
    // else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
    // {
    //     if (iteration == "simweight_iter1")
    //     {
    //         a1 = 4842.04;  b1 = 17.22; a2 = 17.43; b2 = 3.13;
    //     }
    //     else if (iteration == "simweight_iter2")
    //     {
    //         a1 = 11328.73; b1 = 19.75; a2 = 23.50; b2 = 3.54;
    //     }
    //     else if (iteration == "simweight_iter3")
    //     {
    //         a1 = 12903.60; b1 = 20.16; a2 = 24.55; b2 = 3.60;
    //     }
    //     else if (iteration == "simweight_iter4")
    //     {
    //         a1 = 13181.89; b1 = 20.22; a2 = 24.74; b2 = 3.61;
    //     }
    //     else if (iteration == "simweight_iter5")
    //     {
    //         a1 = 13225.39; b1 = 20.23; a2 = 24.77; b2 = 3.61;
    //     }
    //     else if (iteration == "simweight_iter6")
    //     {
    //         a1 = 13232.23; b1 = 20.23; a2 = 24.78; b2 = 3.61;
    //     }
    //     else if (iteration == "simweight_iter7")
    //     {
    //         a1 = 13232.54; b1 = 20.23; a2 = 24.78; b2 = 3.61;
    //     }
    // }
    // else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
    // {
    //     if (iteration == "simweight_iter1")
    //     {
    //         a1 = 2882.53; b1 = 15.54; a2 = 12.54; b2 = 2.90;
    //     }
    //     else if (iteration == "simweight_iter2")
    //     {
    //         a1 = 4602.99; b1 = 16.69; a2 = 13.21; b2 = 2.95;
    //     }
    //     else if (iteration == "simweight_iter3")
    //     {
    //         a1 = 4769.82; b1 = 16.77; a2 = 13.19; b2 = 2.94;
    //     }
    //     else if (iteration == "simweight_iter4")
    //     {
    //         a1 = 4781.31; b1 = 16.78; a2 = 13.18; b2 = 2.94;
    //     }
    //     else if (iteration == "simweight_iter5")
    //     {
    //         a1 = 4782.23; b1 = 16.78; a2 = 13.18; b2 = 2.94;
    //     }
    //     else if (iteration == "simweight_iter6")
    //     {
    //         a1 = 4782.18; b1 = 16.78; a2 = 13.18; b2 = 2.94;
    //     }
    //     else if (iteration == "simweight_iter7")
    //     {
    //         a1 = 4782.12; b1 = 16.78; a2 = 13.18; b2 = 2.94;
    //     }
    // }
    // return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}

double sim_weight_func_systematic(double beam_energy_truth, double minust_truth, int a1_variation, int b1_variation, int a2_variation, int b2_variation)
{
    // No energy dependence
    double a1_nominal, b1_nominal, a2_nominal, b2_nominal = 0;
    double a1_sigma, b1_sigma, a2_sigma, b2_sigma = 0;
    double a1, b1, a2, b2 = 0;
    double normalization = 10;
    if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
        return 1.0;
    else
    {
        a1_nominal = 13232.54; b1_nominal = 20.23; a2_nominal = 24.78; b2_nominal = 3.61;
        a1_sigma = 200.0; b1_sigma = 0.2; a2_sigma = 1.0; b2_sigma = 0.1;
    }

    // With energy dependence
    // double a1_nominal, b1_nominal, a2_nominal, b2_nominal = 0;
    // double a1_sigma, b1_sigma, a2_sigma, b2_sigma = 0;
    // double a1, b1, a2, b2 = 0;
    // double normalization = 10;
    // if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
    //     return 1.0;
    // else if (beam_energy_truth >= 6.0 && beam_energy_truth < 8.0)                            // simulation, weighted by the measured cross section
    // {
    //     a1_nominal = 6445.03;  b1_nominal = 17.22; a2_nominal = 12.42; b2_nominal = 2.49;
    //     a1_sigma = 100.0; b1_sigma = 0.1; a2_sigma = 0.5; b2_sigma = 0.05;
    // }
    // else if (beam_energy_truth >= 8.0 && beam_energy_truth < 9.0)
    // {
    //     a1_nominal = 13232.54; b1_nominal = 20.23; a2_nominal = 24.78; b2_nominal = 3.61;
    //     a1_sigma = 200.0; b1_sigma = 0.2; a2_sigma = 1.0; b2_sigma = 0.1;
    // }
    // else if (beam_energy_truth >= 9.0 && beam_energy_truth < 11.0)
    // {
    //     a1_nominal = 4782.12;  b1_nominal = 16.78; a2_nominal = 13.18; b2_nominal = 2.94;
    //     a1_sigma = 150.0; b1_sigma = 0.15; a2_sigma = 0.7; b2_sigma = 0.07;
    // }

        a1 = a1_nominal + a1_variation*a1_sigma;
        b1 = b1_nominal + b1_variation*b1_sigma;
        a2 = a2_nominal + a2_variation*a2_sigma;
        b2 = b2_nominal + b2_variation*b2_sigma;

    return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}
