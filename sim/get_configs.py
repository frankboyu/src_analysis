import sys
import numpy as np

reaction = sys.argv[1]
version  = sys.argv[2]

config_file = open("configs/wrapper_"+reaction+"_ver"+version+".cfg", "r")
config_lines = config_file.readlines()
for (i, line) in enumerate(config_lines):
    if (line[0:11] == 'FLUX_TO_GEN'):
        index_line = i
        for (j, character) in enumerate(line):
            if (line[j] == '9'):
                index_charater = j
                break
        break
        
if (reaction[-2:] == '2H'):
    run_list = np.loadtxt("../flux/output/deuterium/flux_total_deuterium.txt")[:,0]
elif (reaction[-3:] == '4He'):
    run_list = np.loadtxt("../flux/output/helium/flux_total_helium.txt")[:,0]
elif (reaction[-3:] == '12C'):
    run_list = np.loadtxt("../flux/output/carbon/flux_total_carbon.txt")[:,0]
 
for run_number in run_list:

    run_number = int(run_number)
    
    line = config_lines[i]
    new_line = line[:j] + str(run_number) + line[j+5:j+22] + str(run_number) + line[j+27:]
    config_lines[i] = new_line  

    file = open("configs/wrapper_"+reaction+"_ver"+version+"_"+str(run_number)+".cfg", "w")
    file.writelines(config_lines) 