#IMPORT PACKAGES AND FUNCTIONS
import os, sys, math, array, pprint, rcdb, ccdb, MySQLdb
import numpy as np
import ROOT as root
import matplotlib.pyplot as plt

#FUNCTION DEFINITIONS
def LoadCCDB():

    sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
    provider = ccdb.AlchemyProvider()
    provider.connect(sqlite_connect_str)
    provider.authentication.current_user_name = "psflux_user"
    return provider

rcdb_query = "@is_src_production and @status_approved and run_config == 'FCAL_BCAL_PS_SRC_m9.conf' and beam_on_current > 55.0 and target_type=='FULL & Ready Deuterium'"
db         = rcdb.RCDBProvider(os.environ.get('RCDB_CONNECTION'))
run_list   = db.select_runs(rcdb_query, 90001, 90662)
ccdb_variation = "default"
ccdb_conn = LoadCCDB()

fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 3)
axs = gs.subplots(sharey=True)
fig.subplots_adjust(wspace=0)

run_numbers = np.zeros(len(run_list))
run_numbers_ticks = np.zeros(len(run_list))
hodo_high = np.zeros(len(run_list))
hodo_low = np.zeros(len(run_list))
micro = np.zeros(len(run_list))

offset = 280

#LOOP OVER RUNS
for irun, run in enumerate(run_list):

    #INITIALIZATION
    run_number = run.number
    print("processing: {}".format(run_number))

    accidental_scaling_factor_assignment    = ccdb_conn.get_assignment("/ANALYSIS/accidental_scaling_factor", run_number, ccdb_variation)
    accidental_scaling_factor               = accidental_scaling_factor_assignment     .constant_set.data_table
    
    if (run_number < 90300):
        run_numbers[irun] = run_number
    else:
        run_numbers[irun] = run_number - offset
    run_numbers_ticks[irun] = run_numbers[irun]
    hodo_high[irun] = float(accidental_scaling_factor[0][0])
    micro[irun]     = float(accidental_scaling_factor[0][2])
    hodo_low[irun]  = float(accidental_scaling_factor[0][4])

axs[0].plot(run_numbers, hodo_high, 'o', label="Hodo High")
axs[1].plot(run_numbers, micro,     'o', label="Micro")
axs[2].plot(run_numbers, hodo_low,  'o', label="Hodo Low")
axs[0].text(0.02, 0.98, "Hodoscope high energy",    transform=axs[0].transAxes, verticalalignment="top", fontsize=12)
axs[1].text(0.02, 0.98, "Microscope",               transform=axs[1].transAxes, verticalalignment="top", fontsize=12)
axs[2].text(0.02, 0.98, "Hodoscope low energy",     transform=axs[2].transAxes, verticalalignment="top", fontsize=12)

axs[0].set_ylabel("Accidental Scaling Factor", fontsize=12)
axs[0].set_ylim(0.94, 1.06)
axs[0].set_xticks([90200, 90250, 90550-offset, 90600-offset])
axs[0].set_xticklabels(['90200', '90250', '90550', '90600'], fontsize=12)
axs[0].tick_params(axis='y', labelsize=12)
for i in range(1, 3):
    axs[i].set_xticks([90250, 90550-offset, 90600-offset])
    axs[i].set_xticklabels(['90250', '90550', '90600'], fontsize=12)

for ax in axs:
    ax.axhline(1.0, color='red', linestyle='--', linewidth=1.0)
    ax.set_xlabel("Run Number", fontsize=12)

plt.savefig("accidental_scaling_factor.png")
plt.close()
