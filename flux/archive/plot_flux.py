#!/usr/bin/env python

# Tool for creating flux histograms from CCDB (ver 0.2)
# Modified version of the flux originally used by GlueX

import os,sys
import rcdb
from ROOT import TFile,TGraph,TH1F,TH1D,TF1,gRandom, FILE
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

    int_flux_emin_tot = 0;
    int_flux_7gev_tot = 0;
    int_flux_8gev_tot = 0;

    lumi_emin_tot = 0    
    lumi_7gev_tot = 0    
    lumi_8gev_tot = 0


    low_en  = array('d')
    high_en = array('d')

    flux_tmp = array('d')
    
    hist_low = array('d')
    hist_high = array('d')
    hist_flux = array('d')
    
    indx = 0
    first_run = 0


    parser = OptionParser(usage = "plot_flux.py --run-number RunNumber --input-file InputFile --output-file OutputFile --emin-energy EminEnergy")

    parser.add_option("-r","--run-number", dest="run",
                      help="Run to analyze")
        
    parser.add_option("-f","--input-file", dest="input_file",
                      help="Input file")

    parser.add_option("-o","--output-file", dest="output_file",
                      help="Output file")


    parser.add_option("-e","--emin-energy", dest="emin_energy",
                      help="Minimum energy")
    


    (options, args) = parser.parse_args(sys.argv)
    
    input_file   =  options.input_file 

    output_file  =  options.output_file 

    run_number   =  options.run
        
    run_list     =  array('d')
    
    emin_energy  = 0 

    user_range   = 0

    if options.emin_energy is not  None:

        user_range   =  0
        emin_energy = float(options.emin_energy)
        
        if emin_energy > 0:
            user_range = 1



#    print(run_number)

    if run_number is None:

        with open(input_file) as f:

            print("\n Use list of runs from: %s\n" %(input_file))

            for line in f.readlines():
                
                aaa = int(line.split('_')[1])
                
                run_list.append(aaa)
                
                
    else: 
        run_list.append(int(run_number))
        print("\n Analyze run: %d\n" %(int(run_number)))
    

    print("Number of runs:  %d" %len(run_list)) 


#    print("\n Output file:  %s\n" %(output_file))


    for run in run_list:
        run_numb = int(run)
        print(run_numb)



    for run in run_list:

        run_number = int(run)

        atomic_mass = 0

        if run_number <  90206:   # He
            run_type     =  1
            density      =  0.1217 
            atomic_mass  =  4.003

        if run_number >= 90607  and run_number <= 90660:  # He
            run_type     =  1
            density      =  0.1217 
            atomic_mass  =  4.003  


        if run_number == 90206 or run_number == 90252 or run_number == 90253 or run_number == 90602 or run_number == 90661 or run_number == 90662:
            run_type = 0

        if run_number >= 90208  and run_number <= 90249:  #D
            run_type     =  2
            density      =  0.1638
            atomic_mass  =  2.014

        if run_number >= 90558  and run_number <= 90601:  #D
            run_type     =  2
            density      =  0.1638
            atomic_mass  =  2.014


        if run_number >= 90262 and run_number <= 90536:    #C
            run_type     =  3
            TARGETLENGTH =  1.848
            density      =  1.824
            atomic_mass  =  12.0108

        first_run = first_run + 1


        ccdb_conn = LoadCCDB()
        
                        
        photon_endpoint_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy",
                                                              run_number, VARIATION)
    
        photon_endpoint = photon_endpoint_assignment.constant_set.data_table
    

        tagh_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh/tagged", 
                                                               run_number, VARIATION)
    
        tagh_tagged_flux = tagh_tagged_flux_assignment.constant_set.data_table
    

        tagm_tagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm/tagged", 
                                                               run_number, VARIATION)
 
   
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
                
        if first_run == 1:
            for ii in range(125):     
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
                ps_acc = fPSAcceptance(tagh_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                low_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][1]))
                high_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][2]))
                flux_tmp.append(float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc)
                indx = indx + 1  
   
            for ii in range(102):           
                tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[ii][1]) + 
                                                            float(tagm_scaled_energy[ii][2]))/2.                
                ps_acc = fPSAcceptance(tagm_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                low_en.append( float(photon_endpoint[0][0])*float(tagm_scaled_energy[ii][1]))
                high_en.append( float(photon_endpoint[0][0])*float(tagm_scaled_energy[ii][2]))
                flux_tmp.append(float(tagm_tagged_flux[ii][1])*ps_scale/ps_acc )
                indx = indx + 1

            for ii in range(179,218):   
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
                ps_acc = fPSAcceptance(tagh_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                low_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][1]))
                high_en.append( float(photon_endpoint[0][0])*float(tagh_scaled_energy[ii][2]))
                flux_tmp.append(float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc)
                indx = indx + 1    
        else :
            indx_tmp = 0
            for ii in range(125):     
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
                ps_acc = fPSAcceptance(tagh_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                flux_tmp[indx_tmp] = flux_tmp[indx_tmp] + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc
                indx_tmp = indx_tmp + 1
            for ii in range(102):           
                tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[ii][1]) + 
                                                            float(tagm_scaled_energy[ii][2]))/2.                
                ps_acc = fPSAcceptance(tagm_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                flux_tmp[indx_tmp] = flux_tmp[indx_tmp] + float(tagm_tagged_flux[ii][1])*ps_scale/ps_acc 
                indx_tmp = indx_tmp + 1

            for ii in range(179,218):   
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
                ps_acc = fPSAcceptance(tagh_energy)
                if ps_acc <= 0:
                    ps_acc = 10  
                flux_tmp[indx_tmp] = flux_tmp[indx_tmp] + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc
                indx_tmp = indx_tmp + 1



        int_flux_7gev = 0;
        int_flux_8gev = 0;
        int_flux_emin = 0;

        for ii in range(274):    
        
            tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[ii][1]) + 
                                                        float(tagh_scaled_energy[ii][2]))/2.
        

            ps_acc = fPSAcceptance(tagh_energy)
  
      
            if ps_acc <= 0:
                ps_acc = 10

            if user_range == 1:
                if  tagh_energy > emin_energy:
                    int_flux_emin     = int_flux_emin     + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc
                    int_flux_emin_tot = int_flux_emin_tot + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc

            
            if tagh_energy > 7:
                int_flux_7gev     = int_flux_7gev     + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc
                int_flux_7gev_tot = int_flux_7gev_tot + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc

            if tagh_energy > 8:
                int_flux_8gev     = int_flux_8gev     + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc
                int_flux_8gev_tot = int_flux_8gev_tot + float(tagh_tagged_flux[ii][1])*ps_scale/ps_acc


        for jj in range(102):


            tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[jj][1]) + 
                                                        float(tagm_scaled_energy[jj][2]))/2.

            ps_acc = fPSAcceptance(tagm_energy)

            if ps_acc <= 0:
                ps_acc = 10

            if user_range == 1: 
                if tagm_energy > emin_energy:
                    int_flux_emin      =  int_flux_emin     + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc
                    int_flux_emin_tot  =  int_flux_emin_tot + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc



            if tagm_energy > 7:
                int_flux_7gev      =  int_flux_7gev     + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc
                int_flux_7gev_tot  =  int_flux_7gev_tot + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc

            if tagm_energy > 8:
                int_flux_8gev     = int_flux_8gev     + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc
                int_flux_8gev_tot = int_flux_8gev_tot + float(tagm_tagged_flux[jj][1])*ps_scale/ps_acc



        lumi_7gev = 0    
        lumi_8gev = 0
        lumi_emin = 0


        if run_type != 0:
            if user_range == 1: 
                lumi_emin     = int_flux_emin*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12
                lumi_emin_tot = lumi_emin_tot + int_flux_emin*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12    


            lumi_7gev     = int_flux_7gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12
            lumi_7gev_tot = lumi_7gev_tot + int_flux_7gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12

            lumi_8gev     = int_flux_8gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12  
            lumi_8gev_tot = lumi_8gev_tot + int_flux_8gev*density*TARGETLENGTH*Navagadro*units_cm2_b/atomic_mass/1e12  

            #    print(" Total flux (E > 6 GeV)  (10e12):   %s    Flux  (pb-1)  %f "%(int_flux_7gev/1e12, lumi_7gev))
            #    print(" Total flux (E > 7 GeV)  (10e12):   %s    Flux  (pb-1)  %f "%(int_flux_8gev/1e12, lumi_8gev))
            #    density * TARGETLENGTH * Navagadro * units_cm2_b * units_g_mg


            print(" ")
            print(" Run number: %d"%(run_number))
            print(" Beam energy %s "%(photon_endpoint[0][0]))
            if run_type == 0: 
                print(" Target     Empty  ")
            if run_type == 1: 
                print(" Target     LHe4  ")
            if run_type == 2: 
                print(" Target     LD  ")
            if run_type == 3: 
                print(" Target     Carbon  ")

            print(" Total flux (E > 7.0 GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(int_flux_7gev/1e12, lumi_7gev))
            print(" Total flux (E > 8.0 GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(int_flux_8gev/1e12, lumi_8gev))
            if user_range == 1:
                print(" Total flux (E > %.1f GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(emin_energy, int_flux_emin/1e12, lumi_emin))


    print(" \n")
    print(" =========================================== ")
    print(" ===========  ALL RUNS ===================== ")
    print(" =========================================== ")
    print(" \n")

    print(" Total flux (E > 7.0 GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(int_flux_7gev_tot/1e12, lumi_7gev_tot))
    print(" Total flux (E > 8.0 GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(int_flux_8gev_tot/1e12, lumi_8gev_tot))
    if user_range == 1:
        print(" Total flux (E > %.1f GeV)  (10e12):   %.7f    Lumi (per nucleus)  (pb-1)  %.7f "%(emin_energy,int_flux_emin_tot/1e12, lumi_emin_tot))




#  Add empty bins to the histogram

    nbin = 1
    
    hist_high.append(high_en[0])
#    hist_flux.append(float(tagh_tagged_flux[0][1]))
    hist_flux.append(float(flux_tmp[0]))
           
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


    htest  =  TH1D("TaggedFlux","TaggedFlux",nbin,hist_high1)

# Fill histogram

    for ii in range(nbin):
        htest.SetBinContent(nbin - ii,float(hist_flux[ii]))



    OUTPUT_FILE_ROOT = "flux.root"

    if options.output_file is None:        
        print("\n No output file specified, store histogram in flux.root \n")
    else:        
        OUTPUT_FILE_ROOT = output_file
        print("\n Store flux histogram in:  %s\n"%(OUTPUT_FILE_ROOT))

    froot_out = TFile(OUTPUT_FILE_ROOT, "recreate")


    htest.Write();

    froot_out.Close()

#    print run_type
    

## main function  
if __name__ == "__main__":
    main()
