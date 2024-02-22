#IMPORT PACKAGES AND FUNCTIONS
import os, sys, math, array, pprint, rcdb, ccdb, MySQLdb
import numpy as np
import ROOT as root

#FUNCTION DEFINITIONS
def LoadCCDB():  
    
    sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
#     sqlite_connect_str = "sqlite:///endpoint/ccdb.sqlite"
    provider = ccdb.AlchemyProvider()
    provider.connect(sqlite_connect_str)
    provider.authentication.current_user_name = "psflux_user"
    return provider

def PSAcceptance(x, par):
    
    emin = par[1]
    emax = par[2]

    if (x[0] >= 2*emin and x[0] < (emin + emax)):
        return par[0]*(1-2*emin/x[0])
    elif (x[0] >= (emin + emax) and x[0] < 2*emax):
        return par[0]*(2*emax/x[0] - 1)
    else:
        return 0

#CONSTANT DEFINITIONS
tagh_limits = [1, 125, 179, 274]  # all values are closed brackets
tagm_limits = [1, 102]            # all values are closed brackets
conv_length = 75e-6
Be_rl       = 35.28e-2      
conv_rl     = conv_length/Be_rl
ps_scale    = 1./((7/9.) * conv_rl)
ccdb_variation   = "default"

#FILE IO
target = sys.argv[1]
file_total = open("output/"+target+"/flux_total_"+target+".txt", "w")

#RCDB CONNECTION
if (target == 'deuterium'):
    rcdb_query = "@is_src_production and @status_approved and target_type=='FULL & Ready Deuterium'"
elif (target == 'helium'):
    rcdb_query = "@is_src_production and @status_approved and target_type=='FULL & Ready Helium'"
elif (target == 'carbon'):    
    rcdb_query = "@is_src_production and @status_approved and target_type=='FULL & Ready Carbon'"
elif (target == 'empty'):    
    rcdb_query = "@is_src_production and @status_approved and (target_type=='EMPTY & Ready' or target_type=='OFF')"
    
db         = rcdb.RCDBProvider("mysql://rcdb@hallddb/rcdb")
run_list   = db.select_runs(rcdb_query, 90001, 90662)

#LOOP OVER RUNS
for run in run_list:
    
    #INITIALIZATION
    run_number = run.number
    print("processing: {}".format(run_number))

    if (run_number == 90207 or run_number == 90620):
        print("run skipped")  # flux not processed yet
        continue
 
    #GET THE DATA TABLE FROM CCDB
    ccdb_conn = LoadCCDB()
    
    tagh_tagged_flux_assignment      = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh/tagged", run_number, ccdb_variation)
    tagm_tagged_flux_assignment      = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm/tagged", run_number, ccdb_variation)
    tagh_scaled_energy_assignment    = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/scaled_energy_range",      run_number, ccdb_variation)
    tagm_scaled_energy_assignment    = ccdb_conn.get_assignment("/PHOTON_BEAM/microscope/scaled_energy_range",     run_number, ccdb_variation)
    photon_endpoint_assignment       = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy",                    run_number, ccdb_variation)
    photon_endpoint_calib_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/endpoint_calib",           run_number, ccdb_variation)
    PS_accept_assignment             = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept",   run_number, ccdb_variation)
    
    tagh_tagged_flux      = tagh_tagged_flux_assignment     .constant_set.data_table
    tagm_tagged_flux      = tagm_tagged_flux_assignment     .constant_set.data_table
    tagh_scaled_energy    = tagh_scaled_energy_assignment   .constant_set.data_table
    tagm_scaled_energy    = tagm_scaled_energy_assignment   .constant_set.data_table
    photon_endpoint       = photon_endpoint_assignment      .constant_set.data_table
    photon_endpoint_calib = photon_endpoint_calib_assignment.constant_set.data_table 
    PS_accept             = PS_accept_assignment            .constant_set.data_table
    
    #PRINT THE RAW CONTENT IN THE DATABASE
    file_raw  = open("output/"+target+"/flux_raw_"+str(run_number)+".txt", "w")
    file_raw.write('{:>5.2f}    {:>6.3f}    {:>8.6f}    {:>7.5f}    {:>7.5f}\n\n'.format(float(photon_endpoint_calib[0][0]), float(photon_endpoint[0][0]), float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2])))
    
    for i in range(len(tagh_scaled_energy)):
        file_raw.write('{:>3.0f}    {:>8.6f}    {:>8.6f}    {:>9.3f}    {:>7.3f}\n'.format(float(tagh_scaled_energy[i][0]), float(tagh_scaled_energy[i][1]), float(tagh_scaled_energy[i][2]), float(tagh_tagged_flux[i][1]), float(tagh_tagged_flux[i][2])))
            
    file_raw.write('\n')
    
    for i in range(len(tagm_scaled_energy)):
        file_raw.write('{:>3.0f}    {:>8.6f}    {:>8.6f}    {:>9.3f}    {:>7.3f}\n'.format(float(tagm_scaled_energy[i][0]), float(tagm_scaled_energy[i][1]), float(tagm_scaled_energy[i][2]), float(tagm_tagged_flux[i][1]), float(tagm_tagged_flux[i][2])))
    
    file_raw.close()
    
    #PS ACCEPTANCE CORRECTION
    fPSAcceptance = root.TF1("PSAcceptance", PSAcceptance, 2.0, 12.0, 3);
    fPSAcceptance.SetParameters(float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
    
    #END POINT ENERGY CORRECTION
    delta_E = float(photon_endpoint[0][0]) - float(photon_endpoint_calib[0][0])
    
    #MIN AND MAX ENERGY THAT TAGM COVERS
    energy_tagm_min = float(photon_endpoint_calib[0][0])*float(tagm_scaled_energy[-1][1]) + delta_E
    energy_tagm_max = float(photon_endpoint_calib[0][0])*float(tagm_scaled_energy[0][2])  + delta_E
    
    #PRINT THE FLUX WITH CORRECTIONS AND IN TERMS OF ENERGY BINS
    file_corr = open("output/"+target+"/flux_corr_"+str(run_number)+".txt", "w")
    bin_count = 0
    total_flux = 0

    for i in range(tagh_limits[0]-1, tagh_limits[1]):
        
        energy_low    = float(photon_endpoint_calib[0][0])*float(tagh_scaled_energy[i][1]) + delta_E
        energy_high   = float(photon_endpoint_calib[0][0])*float(tagh_scaled_energy[i][2]) + delta_E
        energy_center = float(photon_endpoint_calib[0][0])*(float(tagh_scaled_energy[i][1]) + float(tagh_scaled_energy[i][2]))/2. + delta_E
            
        ps_acc = fPSAcceptance(energy_center)
        if (ps_acc == 0):
            continue
        
        counter_tagh  = float(tagh_tagged_flux[i][0])
        flux_tagh     = float(tagh_tagged_flux[i][1])*ps_scale/ps_acc
        flux_err_tagh = float(tagh_tagged_flux[i][2])*ps_scale/ps_acc
        
        bin_count += 1
        total_flux += flux_tagh
        file_corr.write('{:>3.0f}    {:>3.0f}    {:>13.10f}    {:>13.10f}    {:>13.10f}    {:>17.16e}    {:>17.16e}\n'.format(bin_count, counter_tagh, energy_low, energy_center, energy_high, flux_tagh, flux_err_tagh))  
    
    file_corr.write('\n')
            
    for i in range(tagm_limits[0]-1, tagm_limits[1]):
        
        energy_low    = float(photon_endpoint_calib[0][0])*float(tagm_scaled_energy[i][1]) + delta_E
        energy_high   = float(photon_endpoint_calib[0][0])*float(tagm_scaled_energy[i][2]) + delta_E
        energy_center = float(photon_endpoint_calib[0][0])*(float(tagm_scaled_energy[i][1]) + float(tagm_scaled_energy[i][2]))/2. + delta_E
            
        ps_acc = fPSAcceptance(energy_center)
        if (ps_acc == 0):
            continue
        
        counter_tagm  = float(tagm_tagged_flux[i][0])
        flux_tagm     = float(tagm_tagged_flux[i][1])*ps_scale/ps_acc
        flux_err_tagm = float(tagm_tagged_flux[i][2])*ps_scale/ps_acc
        
        bin_count += 1
        total_flux += flux_tagm
        file_corr.write('{:>3.0f}    {:>3.0f}    {:>13.10f}    {:>13.10f}    {:>13.10f}    {:>17.16e}    {:>17.16e}\n'.format(bin_count, counter_tagm, energy_low, energy_center, energy_high, flux_tagm, flux_err_tagm))  
    
    file_corr.write('\n')
    
    for i in range(tagh_limits[2]-1, tagh_limits[3]):
        
        energy_low    = float(photon_endpoint_calib[0][0])*float(tagh_scaled_energy[i][1]) + delta_E
        energy_high   = float(photon_endpoint_calib[0][0])*float(tagh_scaled_energy[i][2]) + delta_E
        energy_center = float(photon_endpoint_calib[0][0])*(float(tagh_scaled_energy[i][1]) + float(tagh_scaled_energy[i][2]))/2. + delta_E
            
        ps_acc = fPSAcceptance(energy_center)
        if (ps_acc == 0):
            continue
        
        counter_tagh  = float(tagh_tagged_flux[i][0])
        flux_tagh     = float(tagh_tagged_flux[i][1])*ps_scale/ps_acc
        flux_err_tagh = float(tagh_tagged_flux[i][2])*ps_scale/ps_acc
        
        bin_count += 1
        total_flux += flux_tagh
        file_corr.write('{:>3.0f}    {:>3.0f}    {:>13.10f}    {:>13.10f}    {:>13.10f}    {:>17.16e}    {:>17.16e}\n'.format(bin_count, counter_tagh, energy_low, energy_center, energy_high, flux_tagh, flux_err_tagh))  
    
    file_corr.close()
    file_total.write('{:>3.0f}    {:>17.16e}\n'.format(run_number, total_flux))
    
    #FILL THE FLUX INTO A HISTOGRAM
    flux_all = np.loadtxt("output/"+target+"/flux_corr_"+str(run_number)+".txt")
    hist_edges, hist_flux = array.array('d'), array.array('d')
    
    if (len(flux_all) > 0):
        
        hist_edges.append(flux_all[0][4])
        
        for i in range(len(flux_all)):
            
            hist_edges.append(flux_all[i][2])
            hist_flux. append(flux_all[i][5])
            
            if (i != len(flux_all)-1 and flux_all[i][2] > flux_all[i+1][4]):
                hist_edges.append(flux_all[i+1][4])
                hist_flux. append(0.0)

        hist_edges.reverse()
        hist_flux .reverse()

        hist_object = root.TH1D("TaggedFlux_"+str(run_number), "TaggedFlux_"+str(run_number), len(hist_flux), hist_edges)    
        for i in range(len(hist_flux)):
            hist_object.SetBinContent(i+1, hist_flux[i])

        hist_root = root.TFile("output/"+target+"/flux_hist_"+str(run_number)+".root", "recreate")  
        hist_object.Write()
        hist_root.Close()
    else:
        print("EMPTY FLUX! NO OUTPUT!")

file_total.close()