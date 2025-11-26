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
double sim_weight_func_pass1(double beam_energy_truth, double minust_truth);

Double_t rel_bw(Double_t *x, Double_t *par);
Double_t rel_bw_plus_linear(Double_t *x, Double_t *par);
Double_t rel_bw_plus_fulllinear(Double_t *x, Double_t *par);
Double_t rel_bw_plus_quadratic(Double_t *x, Double_t *par);
Double_t rel_bw_plus_fullquadratic(Double_t *x, Double_t *par);
Double_t rel_bw_plus_phenomenological(Double_t *x, Double_t *par);
Double_t linear(Double_t *x, Double_t *par);
Double_t fulllinear(Double_t *x, Double_t *par);
Double_t quadratic(Double_t *x, Double_t *par);
Double_t fullquadratic(Double_t *x, Double_t *par);
Double_t phenomenological(Double_t *x, Double_t *par);

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
        string KinematicsCut    = "kp_momentum_meas > 0.45 && km_momentum_meas > 0.45 && d_momentum_meas > 0.45 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        string VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        string BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 18.0";

        if      (tag.find("dEdx_1.0") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-4.01*d_momentum_meas+4.88) + 3.26)";
        else if (tag.find("dEdx_1.5") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.85*d_momentum_meas+4.69) + 2.92)";
        else if (tag.find("dEdx_2.5") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.41*d_momentum_meas+4.21) + 2.21)";
        else if (tag.find("dEdx_3.0") != string::npos)
            dEdxCut         = "d_dedx_cdc_keV_per_cm_meas > (TMath::Exp(-3.11*d_momentum_meas+3.90) + 1.83)";
        else if (tag.find("misspminus_0.010") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.010";
        else if (tag.find("misspminus_0.015") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.015";
        else if (tag.find("misspminus_0.025") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.025";
        else if (tag.find("misspminus_0.030") != string::npos)
            MissPMinusCut   = "miss_pminus_meas > -0.030";
        else if (tag.find("chisquared_3.5") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 3.5";
        else if (tag.find("chisquared_4.0") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 4.0";
        else if (tag.find("chisquared_6.0") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 6.0";
        else if (tag.find("chisquared_7.0") != string::npos)
            KinFitChiSqCut   = "chisq_per_ndf_kin < 7.0";
        else if (tag.find("momentum_0.400") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.400 && km_momentum_meas > 0.400 && d_momentum_meas > 0.400 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        else if (tag.find("momentum_0.425") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.425 && km_momentum_meas > 0.425 && d_momentum_meas > 0.425 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        else if (tag.find("momentum_0.475") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.475 && km_momentum_meas > 0.475 && d_momentum_meas > 0.475 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        else if (tag.find("momentum_0.500") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.500 && km_momentum_meas > 0.500 && d_momentum_meas > 0.500 && kp_theta_meas > 2.0 && km_theta_meas > 2.0 && d_theta_meas > 2.0";
        else if (tag.find("theta_1.0") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.40 && km_momentum_meas > 0.40 && d_momentum_meas > 0.40 && kp_theta_meas > 1.0 && km_theta_meas > 1.0 && d_theta_meas > 1.0";
        else if (tag.find("theta_1.5") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.40 && km_momentum_meas > 0.40 && d_momentum_meas > 0.40 && kp_theta_meas > 1.5 && km_theta_meas > 1.5 && d_theta_meas > 1.5";
        else if (tag.find("theta_2.5") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.40 && km_momentum_meas > 0.40 && d_momentum_meas > 0.40 && kp_theta_meas > 2.5 && km_theta_meas > 2.5 && d_theta_meas > 2.5";
        else if (tag.find("theta_3.0") != string::npos)
            KinematicsCut    = "kp_momentum_meas > 0.40 && km_momentum_meas > 0.40 && d_momentum_meas > 0.40 && kp_theta_meas > 3.0 && km_theta_meas > 3.0 && d_theta_meas > 3.0";
        else if (tag.find("vertexZ_13.0") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 13.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        else if (tag.find("vertexZ_13.5") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 13.5 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        else if (tag.find("vertexZ_14.5") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.5 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        else if (tag.find("vertexZ_15.0") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 15.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.0";
        else if (tag.find("vertexR_0.50") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 0.50";
        else if (tag.find("vertexR_0.75") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 0.75";
        else if (tag.find("vertexR_1.25") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.25";
        else if (tag.find("vertexR_1.50") != string::npos)
            VertexCut        = "TMath::Abs(vertex_z_kin - 65.0) < 14.0 && TMath::Sqrt(vertex_x_kin*vertex_x_kin + vertex_y_kin*vertex_y_kin) < 1.50";
        else if (tag.find("beamaccid_5") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 22.0";
        else if (tag.find("beamaccid_3") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) < 14.0";
        else if (tag.find("beamaccid_4out") != string::npos)
            BeamAccidCut     = "TMath::Abs(beam_DeltaT_meas) > 2.0 && TMath::Abs(beam_DeltaT_meas) < 22.0";

        rdf_input = rdf_input   .Filter(dEdxCut.c_str())
                                .Filter(MissPMinusCut.c_str())
                                .Filter(KinFitChiSqCut.c_str())
                                .Filter(KinematicsCut.c_str())
                                .Filter(VertexCut.c_str())
                                .Filter(BeamAccidCut.c_str());

        if      (tag.find("beamaccid_3") != string::npos)
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/6.0");
        else if (tag.find("beamaccid_5") != string::npos)
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/10.0");
        else
            rdf_input = rdf_input   .Define("beamaccid_weight_syst",    "TMath::Abs(beam_DeltaT_meas) < 2.0 ? 1.0 : -1.0/8.0");

        if      (tag.find("comboaccid_1") != string::npos)
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "combo_accid_weight == 1.0 ? 1.0 : 1.0");
        else if (tag.find("comboaccid_-1") != string::npos)
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "combo_accid_weight == 1.0 ? 1.0 : -1.0");
        else
            rdf_input = rdf_input   .Define("combo_accid_weight_syst",  "combo_accid_weight");

        if      (tag.find("simweight_pass0") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");
        else if (tag.find("simweight_pass1") != string::npos)
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_pass1(beam_energy_truth, minust_truth)");
        else
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");

        rdf_input = rdf_input       .Define("yield_weight",             "beamaccid_weight_syst*combo_accid_weight_syst*sim_weight_syst")
                                    .Define("yield_weight_squared",     "yield_weight*yield_weight");
    }
    else if (reaction.find("thrown") != string::npos)
    {
        if      (tag == "simweight_pass0")
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");
        else if (tag == "simweight_pass1")
            rdf_input = rdf_input   .Define("sim_weight_syst",          "sim_weight_func_pass1(beam_energy_truth, minust_truth)");
        else
            rdf_input = rdf_input   .Define("sim_weight_syst",          "1.0");

        rdf_input = rdf_input       .Define("yield_weight",             "sim_weight_syst")
                                    .Define("yield_weight_squared",     "yield_weight*yield_weight");
    }

    // Read bin edges from input text file
    string input_txt_name = Form("configs/bins_%s_%s.txt", channel.c_str(), observable.c_str());
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
    string output_textfile_name = Form("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/yield_phi_d/yield_%s_%s_%s_%s.txt", channel.c_str(), reaction.c_str(), observable.c_str(), tag.c_str());
    string output_pdffile_name = Form("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/yield_phi_d/yield_%s_%s_%s_%s.pdf", channel.c_str(), reaction.c_str(), observable.c_str(), tag.c_str());
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

    double fit_max = 1.08;
    if (tag == "fitmax_1.06")
        fit_max = 1.06;
    else if (tag == "fitmax_1.07")
        fit_max = 1.07;
    else if (tag == "fitmax_1.09")
        fit_max = 1.09;
    else if (tag == "fitmax_1.10")
        fit_max = 1.10;

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
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]), ";m_{K^{+}K^{-}} (GeV/c);Counts", 24, 0.9825, 1.1025},"phi_mass_kin","yield_weight");
            }
            else
            {
                if (observable == "Wcostheta")
                    variable = "decay_costheta_helicity_kin";
                else if (observable == "Wphi")
                    variable = "decay_phi_helicity_kin";
                else if (observable == "WPhi")
                    variable = "polarization_phi_com_kin";
                else if (observable == "Wpsi")
                    variable = "psi_helicity_kin";
                angle_cut = Form("%s>%.2f && %s<%.2f", variable.c_str(), bins[i][4], variable.c_str(), bins[i][5]);
                rdf_bin = rdf_bin.Filter(angle_cut.c_str());
                if (reaction.find("data") != string::npos)
                {
                    angle_center = rdf_bin.Mean(variable.c_str()).GetValue();
                    angle_width  = rdf_bin.StdDev(variable.c_str()).GetValue();
                }
                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f_%.2f_%.2f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]), ";m_{K^{+}K^{-}} (GeV/c);Counts", 24, 0.9825, 1.1025},"phi_mass_kin","yield_weight");
            }

            for (int i = 0; i < hist_bin.GetNbinsX(); i++)
            {
                if (hist_bin.GetBinContent(i) < 0)
                {
                    hist_bin.SetBinContent(i, 0);
                    hist_bin.SetBinError(i, 1);
                }
            }

            if (tag == "fitfunc_none" || reaction.find("sim") != string::npos)
            {
                yield = rdf_bin.Sum("yield_weight").GetValue();
                yield_err = sqrt(rdf_bin.Sum("yield_weight_squared").GetValue());
                hist_bin.Draw();
            }
            else if (tag == "fitfunc_quadratic")
            {
                TF1 fit_func("fit_func", rel_bw_plus_quadratic, 0.99, fit_max, 5);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                // fit_func.FixParameter(3, 0.0035);
                fit_func.SetParLimits(4, 0.01, 100.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *quadratic_function = new TF1("quadratic_function", quadratic, 0.99, fit_max, 1);
                quadratic_function->SetParameter(0, fit_func.GetParameter(4));
                quadratic_function->SetLineColor(kBlue);
                quadratic_function->Draw("same");
                yield = fit_func.GetParameter(0)/0.005;
                yield_err = fit_func.GetParError(0)/0.005;
            }
            else
            {
                TF1 fit_func("fit_func", rel_bw_plus_linear, 0.99, fit_max, 5);
                fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                fit_func.SetParLimits(0, 0.01, 200);
                fit_func.FixParameter(1, 1.019456);
                fit_func.FixParameter(2, 0.00425);
                // fit_func.FixParameter(3, 0.0035);
                fit_func.SetParLimits(4, 0.5, 150.0);
                auto fit_results = hist_bin.Fit(&fit_func, "SLR");
                hist_bin.GetFunction("fit_func")->SetLineColor(kBlack);
                hist_bin.Draw("same");
                yield = 1;
                yield_err = 1;
                TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, fit_max, 4);
                rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                rel_bw_func->SetLineColor(kRed);
                rel_bw_func->Draw("same");
                TF1 *linear_function = new TF1("linear_function", linear, 0.99, fit_max, 1);
                linear_function->SetParameter(0, fit_func.GetParameter(4));
                linear_function->SetLineColor(kBlue);
                linear_function->Draw("same");
                yield = fit_func.GetParameter(0)/0.005;
                yield_err = fit_func.GetParError(0)/0.005;
            }
            canvas->Update();
            canvas->Print((output_pdffile_name+"(").c_str());
            canvas->Clear();
        }
        else if (reaction.find("thrown") != string::npos)
        {
            energy_cut  = Form("beam_energy_truth>%.2f && beam_energy_truth<%.2f", bins[i][0], bins[i][1]);
            t_cut       = Form("minust_truth>%.3f && minust_truth<%.3f", bins[i][2], bins[i][3]);
            rdf_bin = rdf_bin.Filter(energy_cut.c_str()).Filter(t_cut.c_str());

            if (observable == "dsdt")
            {
                cout << energy_cut << " && " << t_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]), ";m_{K^{+}K^{-}} (GeV/c);Counts", 24, 0.9825, 1.1025},"phi_mass_truth","yield_weight");
            }
            else
            {
                if (observable == "Wcostheta")
                    variable = "decay_costheta_helicity_truth";
                else if (observable == "Wphi")
                    variable = "decay_phi_helicity_truth";
                else if (observable == "WPhi")
                    variable = "polarization_phi_com_truth";
                else if (observable == "Wpsi")
                    variable = "psi_helicity_truth";
                angle_cut = Form("%s>%.2f && %s<%.2f", variable.c_str(), bins[i][4], variable.c_str(), bins[i][5]);
                rdf_bin = rdf_bin.Filter(angle_cut.c_str());

                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                hist_bin = *rdf_bin.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f_%.2f_%.2f", bins[i][0], bins[i][1], bins[i][2], bins[i][3], bins[i][4], bins[i][5]), ";m_{K^{+}K^{-}} (GeV/c);Counts", 24, 0.9825, 1.1025},"phi_mass_truth","yield_weight");
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
        else
            fprintf(output_textfile, "%6.4f\t%6.4f\t%6.1f\t%6.1f\t%6.4f\t%6.4f\t%6.3f\t%6.3f\t%6.4f\t%6.4f\t%6.1f\t%6.1f\t%f\t%f\n", energy_center, energy_width, bins[i][0], bins[i][1], t_center, t_width, bins[i][2], bins[i][3], angle_center, angle_width, bins[i][4], bins[i][5], yield, yield_err);
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

    // return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
    return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
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

    // return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
    return par[0] * convol_sum + slope*x[0] + intercept;
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

    // return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
    return par[0] * convol_sum + par[4] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
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

    // return par[0] * convol_sum + slope * (x[0] - 2*mass_kaon);
    return par[0] * convol_sum + par[4]*x[0]*x[0] + par[5]*x[0] + par[6];
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

Double_t fulllinear(Double_t *x, Double_t *par)
{
    // par[0] = linear slope
    // par[1] = linear intercept
    return par[0] * x[0] + par[1];
}

Double_t quadratic(Double_t *x, Double_t *par)
{
    // par[0] = quadratic coefficient
    double mass_kaon = 0.493677;
    Double_t quadratic;
    if (x[0] < 2*mass_kaon)
        quadratic = 0.0;
    else
        quadratic = par[0] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
    return quadratic;
}

Double_t fullquadratic(Double_t *x, Double_t *par)
{
    // par[0] = quadratic coefficient
    // par[1] = linear coefficient
    // par[2] = constant term
    double mass_kaon = 0.493677;
    Double_t quadratic;
    if (x[0] < 2*mass_kaon)
        quadratic = 0.0;
    else
        quadratic = par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
    return quadratic;
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

double sim_weight_func_pass1(double beam_energy_truth, double minust_truth)
{
    double a1 = 2813.72894997;
    double b1 = 15.13997936;
    double a2 = 17.88792021;
    double b2 = 2.98839991;
    double normalization = 10;
    if (beam_energy_truth < 0.01)   // data, with its truth variable set to zero as placeholder
        return 1.0;
    else                            // simulation, weighted by the measured cross section
        return (a1*TMath::Exp(-b1*minust_truth) + a2*TMath::Exp(-b2*minust_truth))/normalization;
}