echo 6.9   > input/theory_paras.txt
echo 11.0  >> input/theory_paras.txt
echo 4.00  >> input/theory_paras.txt
echo 40.0  >> input/theory_paras.txt
echo 10.00 >> input/theory_paras.txt

gfortran -ffixed-line-length-none -o run_theory.out get_edved_wkng_pol.f
./run_theory.out > output/table_phi_d_theory_intuition.txt
python get_theory_intuition.py