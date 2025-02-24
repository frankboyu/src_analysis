#!/bin/bash

start=`date +%s`

# # data
# sh run_batch.sh piminus_p_recon_data_2H_inc
# sh run_batch.sh piminus_p_recon_data_2H_missprot
# sh run_batch.sh piminus_p_recon_data_4He_inc
# sh run_batch.sh piminus_p_recon_data_4He_misshe3
# sh run_batch.sh piminus_p_recon_data_12C_inc
# sh run_batch.sh piminus_p_recon_data_12C_missb11

# sim
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_inc.root output/selectedhist_piminus_p_recon_sim_2H_inc_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_inc.root output/selectedhist_piminus_p_recon_sim_2H_inc_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_flat_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_missprot.root output/selectedhist_piminus_p_recon_sim_2H_missprot_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_2H_missprot.root output/selectedhist_piminus_p_recon_sim_2H_missprot_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_flat_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_sim_4He_inc_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_sim_4He_inc_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_flat_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_misshe3.root output/selectedhist_piminus_p_recon_sim_4He_misshe3_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_misshe3.root output/selectedhist_piminus_p_recon_sim_4He_misshe3_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_flat_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_sim_12C_inc_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_sim_12C_inc_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_flat_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_missb11.root output/selectedhist_piminus_p_recon_sim_12C_missb11_model_all.root
# mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_model_all.root
# sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_missb11.root output/selectedhist_piminus_p_recon_sim_12C_missb11_flat_all.root
# mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_flat_all.root

sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*090213*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_2H_inc.root
mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotinc__B4_F4_T1_S4_gen_MF/*090213*.root' 'gd_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_2H_inc.root
mv output/selectedtree_piminus_p_recon_sim_2H_inc.root output/selectedtree_piminus_p_recon_sim_2H_inc_flat_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver03/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*090213*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_2H_missprot.root
mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_2H_ver04/root/tree_gd_pimprotmissprot__B4_F4_T1_S4_gen_MF/*090213*.root' 'gd_pimprotmissprot__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_2H_missprot.root
mv output/selectedtree_piminus_p_recon_sim_2H_missprot.root output/selectedtree_piminus_p_recon_sim_2H_missprot_flat_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*090061*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_4He_inc.root
mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotinc__B4_F4_T1_S4_gen_MF/*090061*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_4He_inc.root
mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_sim_4He_inc_flat_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver03/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*090061*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_4He_misshe3.root
mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_4He_ver04/root/tree_ghe_pimprotmisshe3__B4_F4_T1_S4_gen_MF/*090061*.root' 'ghe_pimprotmisshe3__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_4He_misshe3.root
mv output/selectedtree_piminus_p_recon_sim_4He_misshe3.root output/selectedtree_piminus_p_recon_sim_4He_misshe3_flat_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*090291*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_12C_inc.root
mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotinc__B4_F4_T1_S4_gen_MF/*090291*.root' 'gc12_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_12C_inc.root
mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_sim_12C_inc_flat_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver03/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*090291*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_12C_missb11.root
mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_model_one.root
sh run_local.sh '/work/halld2/home/boyu/src_analysis/sim/output/piminus_p_12C_ver04/root/tree_gc12_pimprotmissb11__B4_F4_T1_S4_gen_MF/*090291*.root' 'gc12_pimprotmissb11__B4_F4_T1_S4_Tree' 'piminus_p_recon'
rm output/selectedhist_piminus_p_recon_sim_12C_missb11.root
mv output/selectedtree_piminus_p_recon_sim_12C_missb11.root output/selectedtree_piminus_p_recon_sim_12C_missb11_flat_one.root

# # thrown
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

## bggen
# sh run_local.sh '/cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Helium-target-proton_3736/trees/tree_ghe_pimprotinc__B4_F4_T1_S4_bggen_upd/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_bggen_4He_inc_proton.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_bggen_4He_inc_proton.root
# sh run_local.sh '/cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Helium-target-neutron_3739/trees/tree_ghe_pimprotinc__B4_F4_T1_S4_bggen_upd/*.root' 'ghe_pimprotinc__B4_F4_T1_S4_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_4He_inc.root output/selectedhist_piminus_p_recon_bggen_4He_inc_neutron.root
# mv output/selectedtree_piminus_p_recon_sim_4He_inc.root output/selectedtree_piminus_p_recon_bggen_4He_inc_neutron.root
# sh run_local.sh '/cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Carbon-target-proton_3738/trees/tree_gc12_pimprotinc__B4_F4_T2_S5_bggen_upd/*.root' 'gc12_pimprotinc__B4_F4_T2_S5_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_bggen_12C_inc_proton.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_bggen_12C_inc_proton.root
# sh run_local.sh '/cache/halld/gluex_simulations/REQUESTED_MC/bggen_upd-2021-11-nucleus-Carbon-target-neutron_3742/trees/tree_gc12_pimprotinc__B4_F4_T2_S5_bggen_upd/*.root' 'gc12_pimprotinc__B4_F4_T2_S5_Tree' 'piminus_p_recon'
# mv output/selectedhist_piminus_p_recon_sim_12C_inc.root output/selectedhist_piminus_p_recon_bggen_12C_inc_neutron.root
# mv output/selectedtree_piminus_p_recon_sim_12C_inc.root output/selectedtree_piminus_p_recon_bggen_12C_inc_neutron.root

end=`date +%s`
echo "Time taken: $(echo "scale=2; ($end - $start) / 60" | bc -l) minutes"