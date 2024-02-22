## src_analysis_flux

Scripts to calculate the photon flux

JLab: /work/halld2/home/boyu/src_analysis/flux

GitHub: https://github.com/frankboyu/src_analysis/tree/master/flux

### Usage

1.  Set up the environment

    `source env.sh`

2.  Run the script

    `python3.6 get_flux.py TARGET`

    TARGET: 2H, 4He, 12C, empty

### Output

flux_raw_RUN.txt: the raw, uncorrected flux of each run. Format: TAGH/TAGM bin number, scale_min, scale_max, flux, flux_err. The first row is E_endpoint_calib, E_endpoint, ps_acc_scale, ps_acc_emin, ps_acc_emax. The following rows start with TAGH, followed by TAGM.

flux_corr_RUN.txt: photon flux of each run. Format: combined bin number, TAGH/TAGM bin number, E_min (GeV), E_mid (GeV), E_max (GeV), flux (photons), flux_err_stats (photons)

flux_hist_RUN.txt: root file containing a flux histogram of each run. Used by the simulation as the input flux.

flux_total_TARGET.txt: the summed flux of each run