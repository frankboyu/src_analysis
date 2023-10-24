from os import environ
from optparse import OptionParser
import os.path
import rcdb
import ccdb
from ccdb.cmd.console_context import ConsoleContext
import ccdb.path_utils
from ccdb import Directory, TypeTable, Assignment, ConstantSet
from array import array
from datetime import datetime
import mysql.connector
import time
import os
import getpass
import sys
import re
import subprocess
from subprocess import call
import glob
import hddm_s
import socket

def LoadCCDB():
        sqlite_connect_str = "mysql://ccdb_user@hallddb.jlab.org/ccdb"
        provider = ccdb.AlchemyProvider()                           # this class has all CCDB manipulation functions
        provider.connect(sqlite_connect_str)                        # use usual connection string to connect to database
        provider.authentication.current_user_name = "psflux_user"   # to have a name in logs

        return provider
    

def PSAcceptance(x, par0, par1, par2):
    min = par1
    max = par2

    if x > 2.*min and x < min + max:
        return par0*(1-2.*min/x)
    elif x >= min + max:
        return par0*(2.*max/x - 1)

    return 0.

def calcFluxCCDB_untagged(ccdb_conn, run, emin, emax):

        flux = 0.
        VARIATION = "default"
        CALIBTIME = datetime.now()
        CALIBTIME_USER = CALIBTIME
        CALIBTIME_ENERGY = CALIBTIME

        # Set livetime scale factor
        livetime_ratio = 0.0
        try:
                livetime_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/trig_live", run[0], VARIATION, CALIBTIME)
                livetime = livetime_assignment.constant_set.data_table
                if float(livetime[3][1]) > 0.0: # check that livetimes are non-zero
                       livetime_ratio = float(livetime[0][1])/float(livetime[3][1])
                else: # if bad livetime assume ratio is 1
                       livetime_ratio = 1.0
        except:
                livetime_ratio = 1.0 # default to unity if table doesn't exist

        # Conversion factors for total flux
        converterThickness = run[2]
        converterLength = 75e-6 # default is 75 um
        if converterThickness == "Be 750um":
                converterLength = 750e-6
        elif converterThickness != "Be 75um":
                print("Unknown converter thickness for run %s: %s, assuming Be 75um" % (run[0],run[2]))

        berilliumRL = 35.28e-2 # 35.28 cm
        radiationLength = converterLength/berilliumRL
        scale = livetime_ratio * 1./((7/9.) * radiationLength)

        photon_endpoint = array('d')
        tagm_untagged_flux = array('d')
        tagm_scaled_energy = array('d')
        tagh_untagged_flux = array('d')
        tagh_scaled_energy = array('d')

        try:
                photon_endpoint_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy", run[0], VARIATION, CALIBTIME_ENERGY)
                photon_endpoint = photon_endpoint_assignment.constant_set.data_table

                tagm_untagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm/untagged", run[0], VARIATION, CALIBTIME)
                tagm_untagged_flux = tagm_untagged_flux_assignment.constant_set.data_table
                tagm_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/microscope/scaled_energy_range", run[0], VARIATION, CALIBTIME_ENERGY)
                tagm_scaled_energy_table = tagm_scaled_energy_assignment.constant_set.data_table

                tagh_untagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh/untagged", run[0], VARIATION, CALIBTIME)
                tagh_untagged_flux = tagh_untagged_flux_assignment.constant_set.data_table
                tagh_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/scaled_energy_range", run[0], VARIATION, CALIBTIME_ENERGY)
                tagh_scaled_energy_table = tagh_scaled_energy_assignment.constant_set.data_table
                PS_accept_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept", run[0], VARIATION, CALIBTIME)
                PS_accept = PS_accept_assignment.constant_set.data_table
        except:
                print("Missing flux for run number = %d, skipping generation" % run[0])
                return -1.0


        # sum TAGM flux
        for tagm_flux, tagm_scaled_energy in zip(tagm_untagged_flux, tagm_scaled_energy_table):
                tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[1])+float(tagm_scaled_energy[2]))/2.

                if tagm_energy < emin or tagm_energy > emax:
                        continue

                psAccept = PSAcceptance(tagm_energy, float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
                if psAccept <= 0.0:
                        continue

                flux = flux + float(tagm_flux[1]) * scale / psAccept

	# sum TAGH flux
        for tagh_flux, tagh_scaled_energy in zip(tagh_untagged_flux, tagh_scaled_energy_table):
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[1])+float(tagh_scaled_energy[2]))/2.

                if tagh_energy < emin or tagh_energy > emax:
                        continue

                psAccept = PSAcceptance(tagh_energy, float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
                if psAccept <= 0.0:
                        continue

                flux = flux + float(tagh_flux[1]) * scale / psAccept

        return flux

def calcFluxCCDB_tagged(ccdb_conn, run, emin, emax):

        flux = 0.
        VARIATION = "default"
        CALIBTIME = datetime.now()
        CALIBTIME_USER = CALIBTIME
        CALIBTIME_ENERGY = CALIBTIME

        # Set livetime scale factor
        livetime_ratio = 0.0
        try:
                livetime_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/trig_live", run[0], VARIATION, CALIBTIME)
                livetime = livetime_assignment.constant_set.data_table
                if float(livetime[3][1]) > 0.0: # check that livetimes are non-zero
                       livetime_ratio = float(livetime[0][1])/float(livetime[3][1])
                else: # if bad livetime assume ratio is 1
                       livetime_ratio = 1.0
        except:
                livetime_ratio = 1.0 # default to unity if table doesn't exist

        # Conversion factors for total flux
        converterThickness = run[2]
        converterLength = 75e-6 # default is 75 um
        if converterThickness == "Be 750um":
                converterLength = 750e-6
        elif converterThickness != "Be 75um":
                print("Unknown converter thickness for run %s: %s, assuming Be 75um" % (run[0],run[2]))

        berilliumRL = 35.28e-2 # 35.28 cm
        radiationLength = converterLength/berilliumRL
        scale = livetime_ratio * 1./((7/9.) * radiationLength)

        photon_endpoint = array('d')
        tagm_untagged_flux = array('d')
        tagm_scaled_energy = array('d')
        tagh_untagged_flux = array('d')
        tagh_scaled_energy = array('d')

        try:
                photon_endpoint_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/endpoint_energy", run[0], VARIATION, CALIBTIME_ENERGY)
                photon_endpoint = photon_endpoint_assignment.constant_set.data_table

                tagm_untagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagm/tagged", run[0], VARIATION, CALIBTIME)
                tagm_untagged_flux = tagm_untagged_flux_assignment.constant_set.data_table
                tagm_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/microscope/scaled_energy_range", run[0], VARIATION, CALIBTIME_ENERGY)
                tagm_scaled_energy_table = tagm_scaled_energy_assignment.constant_set.data_table

                tagh_untagged_flux_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/tagh/tagged", run[0], VARIATION, CALIBTIME)
                tagh_untagged_flux = tagh_untagged_flux_assignment.constant_set.data_table
                tagh_scaled_energy_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/hodoscope/scaled_energy_range", run[0], VARIATION, CALIBTIME_ENERGY)
                tagh_scaled_energy_table = tagh_scaled_energy_assignment.constant_set.data_table
                PS_accept_assignment = ccdb_conn.get_assignment("/PHOTON_BEAM/pair_spectrometer/lumi/PS_accept", run[0], VARIATION, CALIBTIME)
                PS_accept = PS_accept_assignment.constant_set.data_table
        except:
                print("Missing flux for run number = %d, skipping generation" % run[0])
                return -1.0


        # sum TAGM flux
        for tagm_flux, tagm_scaled_energy in zip(tagm_untagged_flux, tagm_scaled_energy_table):
                tagm_energy = float(photon_endpoint[0][0])*(float(tagm_scaled_energy[1])+float(tagm_scaled_energy[2]))/2.

                if tagm_energy < emin or tagm_energy > emax:
                        continue

                psAccept = PSAcceptance(tagm_energy, float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
                if psAccept <= 0.0:
                        continue

                flux = flux + float(tagm_flux[1]) * scale / psAccept

	# sum TAGH flux
        for tagh_flux, tagh_scaled_energy in zip(tagh_untagged_flux, tagh_scaled_energy_table):
                tagh_energy = float(photon_endpoint[0][0])*(float(tagh_scaled_energy[1])+float(tagh_scaled_energy[2]))/2.

                if tagh_energy < emin or tagh_energy > emax:
                        continue

                psAccept = PSAcceptance(tagh_energy, float(PS_accept[0][0]), float(PS_accept[0][1]), float(PS_accept[0][2]))
                if psAccept <= 0.0:
                        continue

                flux = flux + float(tagh_flux[1]) * scale / psAccept

        return flux
    

db = rcdb.RCDBProvider("mysql://rcdb@hallddb.jlab.org/rcdb")
ccdb_conn = LoadCCDB()
table = db.select_runs("",73260,73263).get_values(['event_count','polarimeter_converter'],True)
print(table)

for runs in table:
    print(runs)
    tagged = calcFluxCCDB_tagged(ccdb_conn, runs, 6.0 , 12.0) / 1.e9
    untagged = calcFluxCCDB_untagged(ccdb_conn, runs, 6.0 , 12.0) / 1.e9
    
    print(tagged, untagged, tagged/untagged)