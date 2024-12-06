import rcdb
import sys

rcdb_query = "@is_src_production and @status_approved and target_type=='FULL & Ready Helium'"

db = rcdb.RCDBProvider("mysql://rcdb@hallddb/rcdb")

file = open("runs.txt", "w")

good_runs = db.select_runs(rcdb_query, 90001, 90662)

for run in good_runs:
    print(run)
    file.write(str(run)[-7:-2]+"\n")
