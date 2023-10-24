import sys
import numpy as np

channel = sys.argv[1]
config_file = open("configs/wrapper_"+channel+".cfg", "r")
config_lines = config_file.readlines()
for (i, line) in enumerate(config_lines):
    if (line[0:11] == 'FLUX_TO_GEN'):
        index_line = i
        for (j, character) in enumerate(line):
            if (line[j] == '9'):
                index_charater = j
                break
        break
        
if (channel[-11:-9] == '2H'):
    run_list = np.loadtxt("../flux/output/flux_total_deuterium.txt")[:,0]
elif (channel[-12:-9] == '4He'):
    run_list = np.loadtxt("../flux/output/flux_total_helium.txt")[:,0]
elif (channel[-12:-9] == '12C'):
    run_list = np.loadtxt("../flux/output/flux_total_carbon.txt")[:,0]
 
for run_number in run_list:

    run_number = int(run_number)
    
    line = config_lines[i]
    new_line = line[:j] + str(run_number) + line[j+5:j+22] + str(run_number) + line[j+27:]
    config_lines[i] = new_line  

    file = open("configs/wrapper_"+channel+"_"+str(run_number)+".cfg", "w")
    file.writelines(config_lines) 