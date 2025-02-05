## src_analysis_flux

Scripts to calculate the photon flux

### Usage

To calculate the flux of each run

`sh run_flux.sh`

### Output

flux_raw_RUN.txt: the raw, uncorrected flux of each run. Format: TAGH/TAGM bin number, scale_min, scale_max, flux, flux_err. The first row is E_endpoint_calib, E_endpoint, ps_acc_scale, ps_acc_emin, ps_acc_emax. The following rows start with TAGH, followed by TAGM.

flux_corr_RUN.txt: photon flux of each run. Format: combined bin number, TAGH/TAGM bin number, E_min (GeV), E_mid (GeV), E_max (GeV), flux (photons), flux_err_stats (photons)

flux_hist_RUN.txt: root file containing a flux histogram of each run. Used by the simulation as the input flux.

flux_total_TARGET.txt: the summed flux of each run on the specified target.

lumi_summed_TARGET.txt: the summed luminosity of each tagger on the specified target. (same as flux for empty target)