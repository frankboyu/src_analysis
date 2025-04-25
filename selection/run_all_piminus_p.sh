#!/bin/bash

# start
start=`date +%s`

# data, batch
sh run_batch.sh piminus_p_recon_data_2H_inc
sh run_batch.sh piminus_p_recon_data_2H_missprot
sh run_batch.sh piminus_p_recon_data_4He_inc
sh run_batch.sh piminus_p_recon_data_4He_misshe3
sh run_batch.sh piminus_p_recon_data_12C_inc
sh run_batch.sh piminus_p_recon_data_12C_missb11

sh run_merge.sh piminus_p_recon_data_2H_inc
sh run_merge.sh piminus_p_recon_data_2H_missprot
sh run_merge.sh piminus_p_recon_data_4He_inc
sh run_merge.sh piminus_p_recon_data_4He_misshe3
sh run_merge.sh piminus_p_recon_data_12C_inc
sh run_merge.sh piminus_p_recon_data_12C_missb11

# sim, batch
# sh run_batch.sh piminus_p_recon_sim_2H_inc_flat
# sh run_batch.sh piminus_p_recon_sim_2H_inc_model
# sh run_batch.sh piminus_p_recon_sim_2H_missprot_flat
# sh run_batch.sh piminus_p_recon_sim_2H_missprot_model
# sh run_batch.sh piminus_p_recon_sim_4He_inc_flat
# sh run_batch.sh piminus_p_recon_sim_4He_inc_model
# sh run_batch.sh piminus_p_recon_sim_4He_misshe3_flat
# sh run_batch.sh piminus_p_recon_sim_4He_misshe3_model
# sh run_batch.sh piminus_p_recon_sim_12C_inc_flat
# sh run_batch.sh piminus_p_recon_sim_12C_inc_model
# sh run_batch.sh piminus_p_recon_sim_12C_missb11_flat
# sh run_batch.sh piminus_p_recon_sim_12C_missb11_model

# sh run_merge.sh piminus_p_recon_sim_2H_inc_flat
# sh run_merge.sh piminus_p_recon_sim_2H_inc_model
# sh run_merge.sh piminus_p_recon_sim_2H_missprot_flat
# sh run_merge.sh piminus_p_recon_sim_2H_missprot_model
# sh run_merge.sh piminus_p_recon_sim_4He_inc_flat
# sh run_merge.sh piminus_p_recon_sim_4He_inc_model
# sh run_merge.sh piminus_p_recon_sim_4He_misshe3_flat
# sh run_merge.sh piminus_p_recon_sim_4He_misshe3_model
# sh run_merge.sh piminus_p_recon_sim_12C_inc_flat
# sh run_merge.sh piminus_p_recon_sim_12C_inc_model
# sh run_merge.sh piminus_p_recon_sim_12C_missb11_flat
# sh run_merge.sh piminus_p_recon_sim_12C_missb11_model

# sim, local
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_inc.root output/selectedhist_piminus_p_recon_sim_2H_inc_model.root
# mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_inc.root output/selectedhist_piminus_p_recon_sim_2H_inc_flat.root
# mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_missprot.root output/selectedhist_piminus_p_recon_sim_2H_missprot_model.root
# mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_missprot.root output/selectedhist_piminus_p_recon_sim_2H_missprot_flat.root
# mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_sim_4He_inc_model.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_sim_4He_inc_flat.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_misshe3.root output/selectedhist_piminus_p_recon_sim_4He_misshe3_model.root
# mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_misshe3.root output/selectedhist_piminus_p_recon_sim_4He_misshe3_flat.root
# mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_sim_12C_inc_model.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_sim_12C_inc_flat.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_missb11.root output/selectedhist_piminus_p_recon_sim_12C_missb11_model.root
# mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_missb11.root output/selectedhist_piminus_p_recon_sim_12C_missb11_flat.root
# mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_flat.root

# thrown, batch
# sh run_batch.sh piminus_p_thrown_tagged_2H_flat
# sh run_batch.sh piminus_p_thrown_tagged_2H_model
# sh run_batch.sh piminus_p_thrown_tagged_4He_flat
# sh run_batch.sh piminus_p_thrown_tagged_4He_model
# sh run_batch.sh piminus_p_thrown_tagged_12C_flat
# sh run_batch.sh piminus_p_thrown_tagged_12C_model

# sh run_merge.sh piminus_p_thrown_tagged_2H_flat
# sh run_merge.sh piminus_p_thrown_tagged_2H_model
# sh run_merge.sh piminus_p_thrown_tagged_4He_flat
# sh run_merge.sh piminus_p_thrown_tagged_4He_model
# sh run_merge.sh piminus_p_thrown_tagged_12C_flat
# sh run_merge.sh piminus_p_thrown_tagged_12C_model

# thrown, local
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_2H_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_2H_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_4He_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_4He_flat.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_12C_model.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/thrown/*.root' 'Thrown_Tree' 'piminus_p_thrown'
# mv output/selectedtree_piminus_p_thrown.root output/selectedtree_piminus_p_thrown_tagged_12C_flat.root

# bggen, batch
# sh run_batch.sh piminus_p_recon_bggen_4He_n_inc
# sh run_batch.sh piminus_p_recon_bggen_4He_n_misshe3
# sh run_batch.sh piminus_p_recon_bggen_4He_p_inc
# sh run_batch.sh piminus_p_recon_bggen_4He_p_misshe3
# sh run_batch.sh piminus_p_recon_bggen_12C_n_inc
# sh run_batch.sh piminus_p_recon_bggen_12C_n_missb11
# sh run_batch.sh piminus_p_recon_bggen_12C_p_inc
# sh run_batch.sh piminus_p_recon_bggen_12C_p_missb11

# sh run_merge.sh piminus_p_recon_bggen_4He_n_inc
# sh run_merge.sh piminus_p_recon_bggen_4He_n_misshe3
# sh run_merge.sh piminus_p_recon_bggen_4He_p_inc
# sh run_merge.sh piminus_p_recon_bggen_4He_p_misshe3
# sh run_merge.sh piminus_p_recon_bggen_12C_n_inc
# sh run_merge.sh piminus_p_recon_bggen_12C_n_missb11
# sh run_merge.sh piminus_p_recon_bggen_12C_p_inc
# sh run_merge.sh piminus_p_recon_bggen_12C_p_missb11

# bggen, local
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/4He_n_ver02/tree_ghe_pimprotinc__B4_F4_T1_S4_bggen/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_4He_inc.root output/selectedhist_piminus_p_recon_bggen_4He_inc_n.root
mv output/selectedtree_piminus_p_recon_bggen_4He_inc.root output/selectedtree_piminus_p_recon_bggen_4He_inc_n.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/4He_p_ver02/tree_ghe_pimprotinc__B4_F4_T1_S4_bggen/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_4He_inc.root output/selectedhist_piminus_p_recon_bggen_4He_inc_p.root
mv output/selectedtree_piminus_p_recon_bggen_4He_inc.root output/selectedtree_piminus_p_recon_bggen_4He_inc_p.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/4He_n_ver02/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_bggen/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_4He_misshe3.root output/selectedhist_piminus_p_recon_bggen_4He_misshe3_n.root
mv output/selectedtree_piminus_p_recon_bggen_4He_misshe3.root output/selectedtree_piminus_p_recon_bggen_4He_misshe3_n.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/4He_p_ver02/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_bggen/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_4He_misshe3.root output/selectedhist_piminus_p_recon_bggen_4He_misshe3_p.root
mv output/selectedtree_piminus_p_recon_bggen_4He_misshe3.root output/selectedtree_piminus_p_recon_bggen_4He_misshe3_p.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/12C_n_ver02/tree_gc12_pimprotinc__B4_F4_T1_S4_bggen/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_12C_inc.root output/selectedhist_piminus_p_recon_bggen_12C_inc_n.root
mv output/selectedtree_piminus_p_recon_bggen_12C_inc.root output/selectedtree_piminus_p_recon_bggen_12C_inc_n.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/12C_p_ver02/tree_gc12_pimprotinc__B4_F4_T1_S4_bggen/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_12C_inc.root output/selectedhist_piminus_p_recon_bggen_12C_inc_p.root
mv output/selectedtree_piminus_p_recon_bggen_12C_inc.root output/selectedtree_piminus_p_recon_bggen_12C_inc_p.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/12C_n_ver02/tree_gc12_pimprotmissb11__B4_F4_T1_S4_bggen/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_12C_missb11.root output/selectedhist_piminus_p_recon_bggen_12C_missb11_n.root
mv output/selectedtree_piminus_p_recon_bggen_12C_missb11.root output/selectedtree_piminus_p_recon_bggen_12C_missb11_n.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/bggen/output/12C_p_ver02/tree_gc12_pimprotmissb11__B4_F4_T1_S4_bggen/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_bggen_Tree' 'piminus_p_recon'
mv output/selectedhist_piminus_p_recon_bggen_12C_missb11.root output/selectedhist_piminus_p_recon_bggen_12C_missb11_p.root
mv output/selectedtree_piminus_p_recon_bggen_12C_missb11.root output/selectedtree_piminus_p_recon_bggen_12C_missb11_p.root

hadd output/selectedhist_piminus_p_recon_bggen_4He_inc.root output/selectedhist_piminus_p_recon_bggen_4He_inc_*.root
hadd output/selectedtree_piminus_p_recon_bggen_4He_inc.root output/selectedtree_piminus_p_recon_bggen_4He_inc_*.root
rm output/selectedhist_piminus_p_recon_bggen_4He_inc_*.root
rm output/selectedtree_piminus_p_recon_bggen_4He_inc_*.root
hadd output/selectedhist_piminus_p_recon_bggen_4He_misshe3.root output/selectedhist_piminus_p_recon_bggen_4He_misshe3_*.root
hadd output/selectedtree_piminus_p_recon_bggen_4He_misshe3.root output/selectedtree_piminus_p_recon_bggen_4He_misshe3_*.root
rm output/selectedhist_piminus_p_recon_bggen_4He_misshe3_*.root
rm output/selectedtree_piminus_p_recon_bggen_4He_misshe3_*.root
hadd output/selectedhist_piminus_p_recon_bggen_12C_inc.root output/selectedhist_piminus_p_recon_bggen_12C_inc_*.root
hadd output/selectedtree_piminus_p_recon_bggen_12C_inc.root output/selectedtree_piminus_p_recon_bggen_12C_inc_*.root
rm output/selectedhist_piminus_p_recon_bggen_12C_inc_*.root
rm output/selectedtree_piminus_p_recon_bggen_12C_inc_*.root
hadd output/selectedhist_piminus_p_recon_bggen_12C_missb11.root output/selectedhist_piminus_p_recon_bggen_12C_missb11_*.root
hadd output/selectedtree_piminus_p_recon_bggen_12C_missb11.root output/selectedtree_piminus_p_recon_bggen_12C_missb11_*.root
rm output/selectedhist_piminus_p_recon_bggen_12C_missb11_*.root
rm output/selectedtree_piminus_p_recon_bggen_12C_missb11_*.root

# end
end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"