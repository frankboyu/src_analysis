import sys
import numpy as np

reaction = sys.argv[1]
events   = sys.argv[2]
      
if (reaction[-2:] == '2H'):
    run_list = np.loadtxt("../flux/output/deuterium/flux_total_deuterium.txt")
elif (reaction[-3:] == '4He'):
    run_list = np.loadtxt("../flux/output/helium/flux_total_helium.txt")
elif (reaction[-3:] == '12C'):
    run_list = np.loadtxt("../flux/output/carbon/flux_total_carbon.txt")
 
run_list[:,1] = run_list[:,1]/run_list[:,1].sum()*int(events)

run_list = run_list.astype(int)

run_list[0,1] += int(events)-run_list[:,1].sum()

np.savetxt("list.txt", run_list, fmt='%i')