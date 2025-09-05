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

Double_t rel_bw_plus_linear(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a linear background function
    // par[0] = BW amplitude
    // par[1] = BW pole mass
    // par[2] = BW pole width
    // par[3] = Gaussian width
    // par[4] = linear slope

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

Double_t rel_bw_plus_quadratic(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor, plus a quadratic background function
    // par[0] = BW amplitude
    // par[1] = BW pole mass
    // par[2] = BW pole width
    // par[3] = Gaussian width
    // par[4] = linear slope
    // par[5] = quadratic coefficient

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

Double_t rel_bw(Double_t *x, Double_t *par)
{
    // Relativistic Breit-Wigner with Blatt-Weisskopf factor
    // par[0] = BW amplitude
    // par[1] = BW pole mass
    // par[2] = BW pole width
    // par[3] = Gaussian width

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

Double_t voigt(Double_t *x, Double_t *par)
{
    // par[0] = Voigt amplitude
    // par[1] = Voigt mean
    // par[2] = Voigt sigma (Gaussian width)
    // par[3] = Voigt gamma (Lorentzian width)
    return par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
}

Double_t voigt_plus_linear(Double_t *x, Double_t *par)
{
    // par[0] = Voigt amplitude
    // par[1] = Voigt mean
    // par[2] = Voigt sigma (Gaussian width)
    // par[3] = Voigt gamma (Lorentzian width)
    // par[4] = linear slope
    // par[5] = linear intercept
    Double_t voigt = par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
    Double_t linear = par[4] * (x[0] - 2*mass_kaon);
    return voigt + linear;
}

Double_t voigt_plus_quadratic(Double_t *x, Double_t *par)
{
    // par[0] = Voigt amplitude
    // par[1] = Voigt mean
    // par[2] = Voigt sigma (Gaussian width)
    // par[3] = Voigt gamma (Lorentzian width)
    // par[4] = linear slope
    // par[5] = linear intercept
    Double_t voigt = par[0] * TMath::Voigt(x[0] - par[1], par[2], par[3]);
    Double_t quadratic;
    if (x[0] < 2*mass_kaon)
        quadratic = 0.0;
    else
        quadratic = par[4] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
    return voigt + quadratic;
}

Double_t linear(Double_t *x, Double_t *par)
{
    // par[0] = linear slope
    // return par[0] * (x[0] - 2*mass_kaon);
    return par[0] * (x[0] - 2*mass_kaon);
}

Double_t quadratic(Double_t *x, Double_t *par)
{
    // par[0] = linear slope
    // par[1] = linear intercept
    double mass_kaon = 0.493677;
    Double_t quadratic;
    if (x[0] < 2*mass_kaon)
        quadratic = 0.0;
    else
        quadratic = par[0] * (x[0]*x[0] - 4*mass_kaon*mass_kaon);
    return quadratic;
}

int get_yield(string channel, string reaction, string observable, string fitfunc, double masscut)
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
    RDataFrame rdf_all(chain);

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
    string output_textfile_name = Form("output/yield_%s_%s_%s_%s.txt", channel.c_str(), reaction.c_str(), observable.c_str(), fitfunc.c_str());
    string output_pdffile_name = Form("output/yield_%s_%s_%s_%s.pdf", channel.c_str(), reaction.c_str(), observable.c_str(), fitfunc.c_str());
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

        if (reaction.find("recon") != string::npos)
        {
            energy_cut  = Form("beam_energy_kin>%.2f && beam_energy_kin<%.2f", bins[i][0], bins[i][1]);
            t_cut       = Form("minust_kin>%.3f && minust_kin<%.3f", bins[i][2], bins[i][3]);
            auto rdf_energy_cut = rdf_all.Filter(energy_cut.c_str()).Define("event_weight_squared", "event_weight*event_weight");
            auto rdf_t_cut      = rdf_energy_cut.Filter(t_cut.c_str());
            if (reaction.find("data") != string::npos)
            {
                energy_center   = rdf_energy_cut.Mean("beam_energy_kin").GetValue();
                energy_width    = rdf_energy_cut.StdDev("beam_energy_kin").GetValue();
                t_center        = rdf_t_cut.Mean("minust_kin").GetValue();
                t_width         = rdf_t_cut.StdDev("minust_kin").GetValue();
            }

            if (observable == "dsdt")
            {
                cout << energy_cut << " && " << t_cut << endl;
                TH1D hist = *rdf_t_cut.Histo1D({Form("hist_%.1f_%.1f_%.3f_%.3f", bins[i][0], bins[i][1], bins[i][2], bins[i][3]), ";m_{K^{+}K^{-}} (GeV/c);Counts", 24, 0.9825, 1.1025},"phi_mass_kin","event_weight");
                for (int i = 0; i < hist.GetNbinsX(); i++)
                {
                    if (hist.GetBinContent(i) < 0)
                    {
                        hist.SetBinContent(i, 0);
                        hist.SetBinError(i, 0);
                    }
                }
                if (fitfunc == "none")
                {
                    yield = rdf_t_cut.Filter("phi_mass_kin>1.00 && phi_mass_kin<1.04").Sum("event_weight").GetValue();
                    yield_err = sqrt(rdf_t_cut.Filter("phi_mass_kin>1.00 && phi_mass_kin<1.04").Sum("event_weight_squared").GetValue());
                    hist.Draw();
                }
                else if (fitfunc == "quadratic")
                {
                    TF1 fit_func("fit_func", rel_bw_plus_quadratic, 0.99, 1.08, 5);
                    fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                    fit_func.SetParLimits(0, 0.01, 200);
                    fit_func.FixParameter(1, 1.019456);
                    fit_func.FixParameter(2, 0.00425);
                    // fit_func.FixParameter(3, 0.0035);
                    fit_func.SetParLimits(4, 0.01, 100.0);
                    auto fit_results = hist.Fit(&fit_func, "SLR");
                    hist.GetFunction("fit_func")->SetLineColor(kBlack);
                    hist.Draw("same");
                    TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, 1.08, 4);
                    rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                    rel_bw_func->SetLineColor(kRed);
                    rel_bw_func->Draw("same");
                    TF1 *quadratic_function = new TF1("quadratic_function", quadratic, 0.99, 1.08, 1);
                    quadratic_function->SetParameter(0, fit_func.GetParameter(4));
                    quadratic_function->SetLineColor(kBlue);
                    quadratic_function->Draw("same");
                    // yield = fit_func.GetParameter(0)/0.005;
                    // yield_err = fit_func.GetParError(0)/0.005;
                    yield = rel_bw_func->Integral(1.00, 1.04)/0.005;
                    yield_err = rel_bw_func->IntegralError(1.00, 1.04, fit_results->GetParams(), fit_results->GetCovarianceMatrix().GetMatrixArray())/0.005;
                }
                else if (fitfunc == "linear")
                {
                    TF1 fit_func("fit_func", rel_bw_plus_linear, 0.99, 1.08, 5);
                    fit_func.SetParameters(1.00, 1.02, 0.0035, 0.004, 20);
                    fit_func.SetParLimits(0, 0.01, 200);
                    fit_func.FixParameter(1, 1.019456);
                    fit_func.FixParameter(2, 0.00425);
                    // fit_func.FixParameter(3, 0.0035);
                    fit_func.SetParLimits(4, 0.5, 50.0);
                    auto fit_results = hist.Fit(&fit_func, "SLR");
                    hist.GetFunction("fit_func")->SetLineColor(kBlack);
                    hist.Draw("same");
                    yield = 1;
                    yield_err = 1;
                    TF1 *rel_bw_func = new TF1("rel_bw_func", rel_bw, 0.99, 1.08, 4);
                    rel_bw_func->SetParameters(fit_func.GetParameter(0), fit_func.GetParameter(1), fit_func.GetParameter(2), fit_func.GetParameter(3));
                    rel_bw_func->SetLineColor(kRed);
                    rel_bw_func->Draw("same");
                    TF1 *linear_function = new TF1("linear_function", linear, 0.99, 1.08, 1);
                    linear_function->SetParameter(0, fit_func.GetParameter(4));
                    linear_function->SetLineColor(kBlue);
                    linear_function->Draw("same");
                    // yield = fit_func.GetParameter(0)/0.005;
                    // yield_err = fit_func.GetParError(0)/0.005;
                    yield = rel_bw_func->Integral(1.00, 1.04)/0.005;
                    yield_err = rel_bw_func->IntegralError(1.00, 1.04, fit_results->GetParams(), fit_results->GetCovarianceMatrix().GetMatrixArray())/0.005;
                }
                canvas->Update();
                canvas->Print((output_pdffile_name+"(").c_str());
                canvas->Clear();
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
                auto rdf_angle_cut = rdf_t_cut.Filter(angle_cut.c_str());
                if (reaction.find("data") != string::npos)
                {
                    angle_center = rdf_angle_cut.Mean(variable.c_str()).GetValue();
                    angle_width  = rdf_angle_cut.StdDev(variable.c_str()).GetValue();
                }
                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                yield           = rdf_angle_cut.Sum("event_weight").GetValue();
                yield_err       = sqrt(rdf_angle_cut.Sum("event_weight_squared").GetValue());
            }
        }
        else if (reaction.find("thrown") != string::npos)
        {
            energy_cut  = Form("beam_energy_truth>%.2f && beam_energy_truth<%.2f", bins[i][0], bins[i][1]);
            t_cut       = Form("minust_truth>%.3f && minust_truth<%.3f", bins[i][2], bins[i][3]);
            // auto rdf_energy_cut = rdf_all.Filter(energy_cut.c_str());
            auto rdf_energy_cut = rdf_all.Filter(energy_cut.c_str()).Define("event_weight_squared", "event_weight*event_weight");
            auto rdf_t_cut      = rdf_energy_cut.Filter(t_cut.c_str());

            if (observable == "dsdt")
            {
                cout << energy_cut << " && " << t_cut << endl;
                // yield           = rdf_t_cut.Count().GetValue();
                // yield_err       = sqrt(yield);
                yield           = rdf_t_cut.Sum("event_weight").GetValue();
                yield_err       = sqrt(rdf_t_cut.Sum("event_weight_squared").GetValue());
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
                auto rdf_angle_cut = rdf_t_cut.Filter(angle_cut.c_str());
                cout << energy_cut << " && " << t_cut << " && " << angle_cut << endl;
                yield           = rdf_angle_cut.Count().GetValue();
                yield_err       = sqrt(yield);
            }
        }

        if (observable == "dsdt")
            fprintf(output_textfile, "%6.4f\t%6.4f\t%6.1f\t%6.1f\t%6.4f\t%6.4f\t%6.3f\t%6.3f\t%f\t%f\n", energy_center, energy_width, bins[i][0], bins[i][1], t_center, t_width, bins[i][2], bins[i][3], yield, yield_err);
        else
            fprintf(output_textfile, "%6.4f\t%6.4f\t%6.1f\t%6.1f\t%6.4f\t%6.4f\t%6.3f\t%6.3f\t%6.4f\t%6.4f\t%6.1f\t%6.1f\t%f\t%f\n", energy_center, energy_width, bins[i][0], bins[i][1], t_center, t_width, bins[i][2], bins[i][3], angle_center, angle_width, bins[i][4], bins[i][5], yield, yield_err);
    }
    canvas->Print((output_pdffile_name+")").c_str());

    return 0;

}