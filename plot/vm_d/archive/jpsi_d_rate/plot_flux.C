// ROOT macro to read genT, extract pBeam energies, and write a histogram.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>

double sigma_jpsi_p(double Ephoton, double t) {
    std::ifstream ifs("/work/halld2/home/boyu/src_analysis/selection/configs/dvmp_xsec_proton2025.csv");
    if (Ephoton < 8.2) return 0.0; // table only goes down to 8.2 GeV, so return 0 below that

    if (!ifs.is_open()) return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    double bestE=0, bestT=0, bestXS=std::numeric_limits<double>::quiet_NaN();
    double bestDist = std::numeric_limits<double>::infinity();
    bool any = false;

    while (std::getline(ifs, line)) {
        // trim leading spaces
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos==std::string::npos) continue;
        if (line[pos]=='#' || line[pos]=='%' || line[pos]=='q') continue; // skip comments

        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')) {
            // trim token
            size_t a = token.find_first_not_of(" \t\r\n");
            size_t b = token.find_last_not_of(" \t\r\n");
            if (a==std::string::npos) cols.push_back(""); else cols.push_back(token.substr(a, b-a+1));
        }

        try {
            double e = std::stod(cols[0]);
            double tt = std::stod(cols[2]);
            double xs = std::stod(cols[3]);

            any = true;

            // keep nearest (Euclidean) if no exact match
            double dx = e - Ephoton;
            double dt = tt + t;  // note: table has t as positive, but function argument is negative, so add to get distance
            double dist = std::sqrt(dx*dx + dt*dt);
            if (dist < bestDist) {
                bestDist = dist;
                bestE = e; bestT = tt; bestXS = xs;
            }
        } catch (...) {
            continue;
        }
    }

    if (!any) return std::numeric_limits<double>::quiet_NaN();
    return bestXS;
}

double sigma_jpsi_d(double Ephoton, double t) {
    std::ifstream ifs("/work/halld2/home/boyu/src_analysis/selection/configs/dvmp_xsec_deuteron2024.csv");
    if (Ephoton < 5.6) return 0.0; // table only goes down to 5.6 GeV, so return 0 below that

    if (!ifs.is_open()) return std::numeric_limits<double>::quiet_NaN();

    std::string line;
    double bestE=0, bestT=0, bestXS=std::numeric_limits<double>::quiet_NaN();
    double bestDist = std::numeric_limits<double>::infinity();
    bool any = false;

    while (std::getline(ifs, line)) {
        // trim leading spaces
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos==std::string::npos) continue;
        if (line[pos]=='#' || line[pos]=='%' || line[pos]=='q') continue; // skip comments

        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ',')) {
            // trim token
            size_t a = token.find_first_not_of(" \t\r\n");
            size_t b = token.find_last_not_of(" \t\r\n");
            if (a==std::string::npos) cols.push_back(""); else cols.push_back(token.substr(a, b-a+1));
        }

        try {
            double e = std::stod(cols[0]);
            double tt = std::stod(cols[1]);
            double xs = std::stod(cols[3]);

            any = true;

            // keep nearest (Euclidean) if no exact match
            double dx = e - Ephoton;
            double dt = tt + t;  // note: table has t as positive, but function argument is negative, so add to get distance
            double dist = std::sqrt(dx*dx + dt*dt);
            if (dist < bestDist) {
                bestDist = dist;
                bestE = e; bestT = tt; bestXS = xs;
            }
        } catch (...) {
            continue;
        }
    }

    if (!any) return std::numeric_limits<double>::quiet_NaN();
    return bestXS;
}

void plot_flux(const char* inputFile = "/work/halld2/home/boyu/src_analysis/sim/output/jpsi_p_1H_test/root/generator/genOut_gen_coherent_090213_000.root", const char* outputFile = "jpsi_p.root")
{
	TFile* in = TFile::Open(inputFile, "READ");
	if (!in || in->IsZombie()) {
		std::cerr << "Failed to open input file: " << inputFile << std::endl;
		return;
	}

	TTree* genT = nullptr;
	in->GetObject("genT", genT);
	if (!genT) {
		std::cerr << "Failed to find tree genT in file: " << inputFile << std::endl;
		in->Close();
		return;
	}

	TLorentzVector* pBeam = nullptr;
    TLorentzVector* pMeson = nullptr;
	genT->SetBranchAddress("pBeam", &pBeam);
    genT->SetBranchAddress("pMeson", &pMeson);

    double minust = 0.0; // placeholder value, replace with actual calculation if needed

	TH1D* hEBeam = new TH1D("hEBeam", "pBeam Energy;Energy;Entries", 90, 3.0, 12.0);

	const Long64_t nEntries = genT->GetEntries();
	for (Long64_t i = 0; i < nEntries; ++i) {
		genT->GetEntry(i);
		if (!pBeam) continue;

        minust = -(*pBeam - *pMeson).M2(); // calculate -t using the beam and meson 4-vectors

		hEBeam->Fill(pBeam->E(), sigma_jpsi_p(pBeam->E(), minust)); // weight by the cross section for proton case; use sigma_jpsi_d for deuteron case
	}

	TFile* out = TFile::Open(outputFile, "RECREATE");
	if (!out || out->IsZombie()) {
		std::cerr << "Failed to open output file: " << outputFile << std::endl;
		in->Close();
		return;
	}

	out->cd();
	hEBeam->Write();
	out->Close();
	in->Close();
}
