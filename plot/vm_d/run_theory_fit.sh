gfortran -ffixed-line-length-none -c get_edved_wkng_pol_callable.f -o exe_theory_callable.out
g++ get_sigma.C exe_theory_callable.out -o exe_theory_fit.out `root-config --cflags --libs` -lMinuit -lgfortran
./exe_theory_fit.out output/table_vm_d_dsdt.txt