import os, sys, math, array, pprint, rcdb, ccdb, MySQLdb
import numpy as np
import ROOT as root

db         = rcdb.RCDBProvider(os.environ.get('RCDB_CONNECTION'))
rcdb_query = "@is_src_production and @status_approved"
good_runs = db.select_runs(rcdb_query, 90001, 90662)

file = open("runs.txt", "w")

for run in good_runs:
    print(run)
    file.write(str(run)[-7:-2]+"\n")
