import sys
import numpy as np

channel = sys.argv[1]
events = sys.argv[2]
      
if (channel[-11:-9] == '2H'):
    run_list = np.loadtxt("../flux/output/flux_total_deuterium.txt")
elif (channel[-12:-9] == '4He'):
    run_list = np.loadtxt("../flux/output/flux_total_helium.txt")
elif (channel[-12:-9] == '12C'):
    run_list = np.loadtxt("../flux/output/flux_total_carbon.txt")
 
run_list[:,1] = run_list[:,1]/run_list[:,1].sum()*int(events)

run_list = run_list.astype(int)

run_list[0,1] += int(events)-run_list[:,1].sum()

np.savetxt("list.txt", run_list, fmt='%i')