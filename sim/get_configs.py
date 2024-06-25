import sys
import numpy as np

reaction = sys.argv[1]
version  = sys.argv[2]

if (reaction[-2:] == '2H'):
    run_list = np.loadtxt("../flux/output/deuterium/flux_total_deuterium.txt")[:,0]
elif (reaction[-3:] == '4He'):
    run_list = np.loadtxt("../flux/output/helium/flux_total_helium.txt")[:,0]
elif (reaction[-3:] == '12C'):
    run_list = np.loadtxt("../flux/output/carbon/flux_total_carbon.txt")[:,0]

config_file = open("configs/wrapper/wrapper_"+reaction+"_ver"+version+".cfg", "r")
config_lines = config_file.readlines()
for (i, line) in enumerate(config_lines):
    if (line[0:11] == 'FLUX_TO_GEN'):
        line_flux = i
        for (j, character) in enumerate(line):
            if (line[j] == '9'):
                char_flux = j
                break

    if (line[0:3] == 'BKG'):
        line_bkg = i
        for (j, character) in enumerate(line):
            if (line[j] == 'R'):
                char_bkg = j
                break


for run_number in run_list:

    run_number = int(run_number)

    # Replace the flux histogram for every run
    line = config_lines[line_flux]
    new_line = line[:char_flux] + str(run_number) + line[char_flux+5:char_flux+22] + str(run_number) + line[char_flux+27:]
    config_lines[line_flux] = new_line

    # Replace the random background for runs without corresponding files
    if (run_number == 90140 or run_number == 90615 or run_number == 90616):
        line = config_lines[line_bkg]
        new_line = line[:char_bkg] + "None\n"
        config_lines[line_bkg] = new_line
    else:
        line = config_lines[line_bkg]
        new_line = line[:char_bkg] + "Random:recon-2021_11-ver01\n"
        config_lines[line_bkg] = new_line

    file = open("configs/wrapper/wrapper_"+reaction+"_ver"+version+"_"+str(run_number)+".cfg", "w")
    file.writelines(config_lines)