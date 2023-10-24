#!/usr/bin/env python

# Tool for creating flux histograms from CCDB (ver 0.2)
# Modified version of the flux originally used by GlueX

import os,sys
import rcdb
from ROOT import TFile,TGraph,TH1F,TF1,gRandom, FILE
from optparse import OptionParser
from array import array
from datetime import datetime
import pprint
import math
import MySQLdb

import ccdb
from ccdb import Directory, TypeTable, Assignment, ConstantSet

def LoadCCDB():
    sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
#    sqlite_connect_str = "sqlite:////home/somov/ccdb.sqlite"
    provider = ccdb.AlchemyProvider()                           # this class has all CCDB manipulation functions
    provider.connect(sqlite_connect_str)                        # use usual connection string to connect to database
    provider.authentication.current_user_name = "somov"   # to have a name in logs

    return provider

def loadCCDBContextList(runPeriod, restVer):
    dbhost = "hallddb.jlab.org"
    dbuser = 'datmon'
    dbpass = ''
    dbname = 'data_monitoring'

    conn=MySQLdb.connect(host=dbhost, user=dbuser, db=dbname)
    curs=conn.cursor()    

    cmd = "SELECT revision,ccdb_context FROM version_info WHERE run_period=%s AND data_type='recon' AND revision<=%s ORDER BY revision DESC"
    curs.execute(cmd, [runPeriod, restVer])
    rows=curs.fetchall()
    return rows

def PSAcceptance(x, par):

    min = par[1]
    max = par[2]

    if x[0] > 2*min and x[0] < min + max:
        return par[0]*(1-2*min/x[0])
    elif x[0] >= min + max:
        return par[0]*(2*max/x[0] - 1)
    
    return 0.

def main():

    RCDB_QUERY = "@is_production"
    RCDB_POLARIZATION = "" # AMO, PARA or PERP
    RCDB_POL_ANGLE = "" # 0, 45, 90, 135 (only 2017 and later)
    VARIATION = "default"
    CALIBTIME = datetime.now()
    RESTVERSION = 999


    run_type = 0      # 0 - empty; 1 - LHe;  2 - LD;  3 - Carbon 

        
    TARGETLENGTH  =  29.5 # length in CM
    Navagadro     =  6.02214e23 # atoms/mol
    units_cm2_b   =  1e-24 # 1e-24 cm^2 = 1 barn
    density       =  0;
    

    VARIATION   =  "default"

    conv_length   =  75e-6
    Be_rl         =  35.28e-2      
    conv_rl       =  conv_length/Be_rl;   # 2.1259 10-3   Used in the PS acceptance determination
        
    ps_scale = 1./((7/9.) * conv_rl)


    parser = OptionParser(usage = "plot_flux_ccdb.py --begin-run beginRun --end-run endRun")

    parser.add_option("-b","--begin-run", dest="begin_run",
                      help="Starting run for output")

    parser.add_option("-e","--end-run", dest="end_run",
                     help="Ending run for output")
    
    (options, args) = parser.parse_args(sys.argv)

    FIRST_RUN = int(options.begin_run)

    run_number = FIRST_RUN


    rcdb_conn = rcdb.RCDBProvider("mysql://rcdb@hallddb.jlab.org/rcdb")
    runs = rcdb_conn.select_runs("", run_number, run_number)

    n_ps_90253 = 8.96e6

    n_evt_run = 0

#    print n_ps_run


    for run in runs:
#        print run.get_condition('event_count').value
        n_evt_run  = int(run.get_condition('event_count').value)

    PS_frac   = 0

#   He, PS_frac  =  0.069 
#   n_evt_run    =  29.7e9

#   D, PS_frac   =  0.05 
#   n_evt_run    =  16.2e9

#   C, PS_frac   =  0.068 
#    n_evt_run    =  40.94e9

#    n_ps_run = n_evt_run*0.068


    atomic_mass = 0

    if run_number <  90206:   # He
        run_type     =  1
        density      =  0.1217 
        atomic_mass  =  4.003
        PS_frac      =  0.069

    if run_number == 90206 or run_number == 90252 or run_number == 90253:
        run_type = 0

    if run_number >= 90208  and run_number <= 90249:  #D
        run_type     =  2
        density      =  0.1638
        atomic_mass  =  2.014
        PS_frac      =  0.05

    if run_number > 90253:        
        run_type     =  3
        TARGETLENGTH =  1.848
        density      =  1.824
        atomic_mass  =  12.0108
        PS_frac      =  0.068


    n_ps_run    =  n_evt_run*PS_frac 

    scale_flux  =  n_ps_run/n_ps_90253


    OUTPUT_FILE_TAGH = "%d_tagh_ps_acc_cor.txt" %(FIRST_RUN)
    OUTPUT_FILE_TAGM = "%d_tagm_ps_acc_cor.txt" %(FIRST_RUN)
    OUTPUT_FILE_ROOT = "%d_flux.root" %(FIRST_RUN)


    ccdb_conn = LoadCCDB()
        

    htagh_flux      =  TH1F("TAGH Tagged Flux","TAGH_Tagged_Flux",300,0.5,299.5)
    htagh_flux_cor  =  TH1F("TAGH Tagged Flux Cor","TAGH_Tagged_Flux_Cor",300,0.5,299.5)

    htagm_flux      =  TH1F("TAGM Tagged Flux","TAGM_Tagged_Flux",102,0.5,102.5)
    htagm_flux_cor  =  TH1F("TAGM Tagged Flux Cor","TAGM_Tagged_Flux_Cor",102,0.5,102.5)


    photon_endpoint_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy",
                                                          run_number, VARIATION)
    
    photon_endpoint = photon_endpoint_assignment.constant_set.data_table
    
    tagh_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh_tagged", 
                                                           run_number, VARIATION)

#    tagh_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh/tagged", 
#                                                           run_number, VARIATION)
    
    tagh_tagged_flux = tagh_tagged_flux_assignment.constant_set.data_table
    
    tagm_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm_tagged", 
                                                           run_number, VARIATION)
 
#    tagm_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm/tagged", 
#                                                           run_number, VARIATION)
 

   
    tagm_tagged_flux = tagm_tagged_flux_assignment.constant_set.data_table


    tagh_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/scaled_energy_range",run_number, VARIATION)

    tagh_scaled_energy = tagh_scaled_energy_assignment.constant_set.data_table

    tagm_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/microscope/scaled_energy_range",run_number, VARIATION)

    tagm_scaled_energy = tagm_scaled_energy_assignment.constant_set.data_table

    # PS acceptance correction
    fPSAcceptance = TF1("PSAcceptance", PSAcceptance, 2.0, 12.0, 3);

    # Get parameters from CCDB 
    PS_accept_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept", run_number, VARIATION)
    
    PS_accept = PS_accept_assignment.constant_set.data_table
    

    fPSAcceptance.SetParameters(float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
    

    int_flux_7gev = 0;
    int_flux_8gev = 0;

#    print ps_scale

    ps_scale = ps_scale*scale_flux


    debug = 0

    low_en  = array('d')
    high_en = array('d')
    mean_en = array('d')

    flux_tmp = array('d')

    hist_low = array('d')
    hist_high = array('d')
    hist_flux = array('d')

    indx = 0

    for ii in range(125):
        
        if float(tagh_tagged_flux[ii][1]) > 0.:
      
            tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.            

            low_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][1]))
            high_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][2]))
            
            ps_acc = fPSAcceptance(tagh_energy)                        
      
            if ps_acc <= 0:
                ps_acc = 10
 

            flux_tmp.append(float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc)
            
            mean_en.append(tagh_energy)
            
            if debug == 1:
                print(ii,low_en[indx],high_en[indx],"%.0f" % flux_tmp[indx])

            indx = indx + 1

            
    for ii in range(102):
            
        if float(tagm_tagged_flux[ii][1]) > 0.:
        
            tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[ii][1]) + 
                                                        float(tagm_scaled_energy[ii][2]))/2.
            
            low_en.append( float(photon_endpoint[0][0])*float(tagm_scaled_energy[ii][1]))
            high_en.append( float(photon_endpoint[0][0])*float(tagm_scaled_energy[ii][2]))

            ps_acc = fPSAcceptance(tagm_energy)                        
      
            if ps_acc <= 0:
                ps_acc = 10                
     
            flux_tmp.append(float(tagm_tagged_flux[ii][1])*ps_scale/ps_acc )
            
            mean_en.append(tagm_energy)    
            
            if debug == 1:
                print(ii,low_en[indx],high_en[indx],"%.0f" % flux_tmp[indx])

            indx = indx + 1
            
    for ii in range(179,218):
            
        if float(tagh_tagged_flux[ii][1]) > 0.:
            
            tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
            
            low_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][1]))
            high_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][2]))
            
            ps_acc = fPSAcceptance(tagh_energy)                        
      
            if ps_acc <= 0:
                ps_acc = 10 

            flux_tmp.append(float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc)
            
            mean_en.append(tagh_energy)

            if debug == 1:
                print(ii,low_en[indx],high_en[indx],"%.0f" % flux_tmp[indx])

            indx = indx + 1


#  Add empty bins to the histogram

    nbin = 1
    
    hist_high.append(high_en[0])
    hist_flux.append(float(tagh_tagged_flux[0][1]))
            
    for ii in range(1,indx):
        
#        print('Add empty bins: ',low_en[ii],high_en[ii],low_en[ii-1],high_en[ii-1])

        if (high_en[ii] < low_en[ii-1]):
            hist_high.append(low_en[ii-1])
            hist_flux.append(0.)
            nbin = nbin + 1
            hist_high.append(high_en[ii])
            hist_flux.append(float(flux_tmp[ii]))
            nbin = nbin + 1
                
        elif (high_en[ii] >= low_en[ii-1]):
            hist_high.append(high_en[ii])
            hist_flux.append(float(flux_tmp[ii]))
            nbin = nbin + 1
                
    hist_high.append(low_en[indx-1])


    print('Number of bins and size of an array with edges',nbin,len(hist_high),len(hist_flux))

    hist_high1 = array('d')

    for ii in range(len(hist_high)):
        hist_high1.append(hist_high[len(hist_high) - ii -1])


#    for ii in range(nbin):
#        print('Test: ',hist_high[ii],hist_flux[ii],ii)


    htest  =  TH1F("TaggedFlux","TaggedFlux",nbin,hist_high1)

# Fill histogram

    for ii in range(nbin):
        htest.SetBinContent(nbin - ii,float(hist_flux[ii]))
 


    data_file_tagh = open(OUTPUT_FILE_TAGH,"w+")
    
    for ii in range(274):
        

#        htest.SetBinContent(169-int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][1]))

        htagh_flux.SetBinContent(int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][1]))
        htagh_flux.SetBinError(int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][2]))
        
        tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                    float(tagh_scaled_energy[ii][2]))/2.
        

        ps_acc = fPSAcceptance(tagh_energy)
  

#        print tagh_energy, ps_acc
      
        if ps_acc <= 0:
            ps_acc = 10
            
        htagh_flux_cor.SetBinContent(int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc)
        htagh_flux_cor.SetBinError(int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][2])*ps_scale/ps_acc)
                

        if tagh_energy > 7:
            int_flux_7gev = int_flux_7gev + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc

        if tagh_energy > 8:
            int_flux_8gev = int_flux_8gev + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc


        data_file_tagh.write("%4d  %10.3f  %10.3f \n" %(int(tagh_tagged_flux[ii][0]),float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc,float(tagh_tagged_flux[ii][2])*ps_scale/ps_acc))

        

    data_file_tagm = open(OUTPUT_FILE_TAGM,"w+")

    for jj in range(102):

        htagm_flux.SetBinContent(int(tagm_tagged_flux[jj][0]),float(tagm_tagged_flux[jj][1]))
        htagm_flux.SetBinError(int(tagm_tagged_flux[jj][0]),float(tagm_tagged_flux[jj][2]))

        tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[jj][1]) + 
                                                    float(tagm_scaled_energy[jj][2]))/2.

        ps_acc = fPSAcceptance(tagm_energy)


        if ps_acc <= 0:
            ps_acc = 10

        htagm_flux_cor.SetBinContent(int(tagm_tagged_flux[jj][0]),float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc)
        htagm_flux_cor.SetBinError(int(tagm_tagged_flux[jj][0]),float(tagm_tagged_flux[jj][2])*ps_scale/ps_acc)


        if tagm_energy > 7:
            int_flux_7gev = int_flux_7gev + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc

        if tagm_energy > 8:
            int_flux_8gev = int_flux_8gev + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc


        data_file_tagm.write("%4d  %10.3f  %10.3f \n" %(int(tagm_tagged_flux[jj][0]),float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc,float(tagm_tagged_flux[jj][2])*ps_scale/ps_acc))


    lumi_7gev = 0    
    lumi_8gev = 0


    if run_type != 0:
        
        lumi_7gev = int_flux_7gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12

        lumi_8gev = int_flux_8gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12  

#    print(" Total flux (E > 6 GeV)  (10e12):   %s    Flux  (pb-1)  %f "%(int_flux_7gev/1e12, lumi_7gev))

#    print(" Total flux (E > 7 GeV)  (10e12):   %s    Flux  (pb-1)  %f "%(int_flux_8gev/1e12, lumi_8gev))

#    density * TARGETLENGTH * Navagadro * units_cm2_b * units_g_mg


    froot_out = TFile(OUTPUT_FILE_ROOT, "recreate")

    htagh_flux.Write()
    htagh_flux_cor.Write()

    htagm_flux.Write()
    htagm_flux_cor.Write()

    htest.Write();


    froot_out.Close()

    print (" ")
    print(" Run number: %d"%(FIRST_RUN))
    print(" Beam energy %s "%(photon_endpoint[0][0]))
    if run_type == 0: 
        print(" Target     Empty  ")
    if run_type == 1: 
        print(" Target     LHe4  ")
    if run_type == 2: 
        print(" Target     LD  ")
    if run_type == 3: 
        print(" Target     Carbon  ")
    print("")
    print(" Number of events       %s  " %(n_evt_run))
    print(" Number of PS triggers  %s  " %(n_ps_run ))
    print (" ")
    print(" Total flux (E > 7 GeV)  (10e12):   %s    Lumi (per nucleus)  (pb-1)  %f "%(int_flux_7gev/1e12, lumi_7gev))
    print(" Total flux (E > 8 GeV)  (10e12):   %s    Lumi (per nucleus)  (pb-1)  %f "%(int_flux_8gev/1e12, lumi_8gev))

    print (" ")
    print(" Output file TAGH:  %s "%(OUTPUT_FILE_TAGH))
    print(" Output file TAGM:  %s "%(OUTPUT_FILE_TAGM))
    print(" Root file:         %s "%(OUTPUT_FILE_ROOT))

#    print run_type
    

## main function  
if __name__ == "__main__":
    main()
