import numpy as np
import matplotlib.pyplot as plt
import ROOT as root

mass_proton = 0.938272
mass_neutron = 0.939565
mass_piminus = 0.139570
Deg_To_Rad = np.pi/180

energy_list = np.arange(6.0, 11.0, 0.5)
theta_list = np.arange(1, 180, 1)

piminus_p = np.zeros((len(energy_list), len(theta_list)))
piminus_theta = np.zeros((len(energy_list), len(theta_list)))
proton_p = np.zeros((len(energy_list), len(theta_list)))
proton_theta = np.zeros((len(energy_list), len(theta_list)))
minus_t = np.zeros((len(energy_list), len(theta_list)))
minus_u = np.zeros((len(energy_list), len(theta_list)))

for i, energy in enumerate(energy_list):
    for j, theta in enumerate(theta_list):
        P4_gamma = root.TLorentzVector(0, 0, energy, energy)
        P4_neutron = root.TLorentzVector(0, 0, 0, mass_neutron)
        s = (P4_gamma + P4_neutron).Mag2()

        pf_cm = np.sqrt(((s-mass_piminus**2-mass_proton**2)**2-4*mass_piminus**2*mass_proton**2)/(4*s))

        P4_piminus = root.TLorentzVector()
        P4_proton = root.TLorentzVector()
        P4_piminus.SetXYZM(0, pf_cm*np.sin(theta*Deg_To_Rad), pf_cm*np.cos(theta*Deg_To_Rad), mass_piminus)
        P4_proton.SetXYZM(0, -pf_cm*np.sin(theta*Deg_To_Rad), -pf_cm*np.cos(theta*Deg_To_Rad), mass_proton)

        boost_cm = (P4_gamma+P4_neutron).BoostVector()

        P4_piminus.Boost(boost_cm)
        P4_proton.Boost(boost_cm)

        piminus_p[i][j] = P4_piminus.P()
        piminus_theta[i][j] = P4_piminus.Theta()/Deg_To_Rad
        proton_p[i][j] = P4_proton.P()
        proton_theta[i][j] = P4_proton.Theta()/Deg_To_Rad

        minus_t[i][j] = -(P4_gamma - P4_piminus).Mag2()
        minus_u[i][j] = -(P4_gamma - P4_proton).Mag2()


for i, energy in enumerate(energy_list):
    plt.plot(piminus_p[i], piminus_theta[i], label=format(energy, '.1f')+' GeV')
plt.xlabel('Pion momentum [GeV]')
plt.ylabel('Pion theta [deg]')
plt.title('Pion momentum vs theta')
plt.legend()
plt.grid()
plt.savefig('output/piminus_p_theta.png')
plt.close()

for i, energy in enumerate(energy_list):
    plt.plot(proton_p[i], proton_theta[i], label=format(energy, '.1f')+' GeV')
plt.xlabel('Proton momentum [GeV]')
plt.ylabel('Proton theta [deg]')
plt.title('Proton momentum vs theta')
plt.legend()
plt.grid()
plt.savefig('output/proton_p_theta.png')
plt.close()

for i, energy in enumerate(energy_list):
    plt.plot(piminus_theta[i], proton_theta[i], label=format(energy, '.1f')+' GeV')
plt.xlabel('Pion theta [deg]')
plt.ylabel('Proton theta [deg]')
plt.title('Pion theta vs Proton theta')
plt.legend()
plt.grid()
plt.savefig('output/piminus_theta_proton_theta.png')
plt.close()

for i, energy in enumerate(energy_list):
    plt.plot(piminus_p[i], proton_p[i], label=format(energy, '.1f')+' GeV')
plt.xlabel('Pion momentum [GeV]')
plt.ylabel('Proton momentum [GeV]')
plt.title('Pion momentum vs Proton momentum')
plt.legend()
plt.grid()
plt.savefig('output/piminus_p_proton_p.png')
plt.close()

for i, energy in enumerate(energy_list):
    plt.plot(minus_t[i], theta_list, label=format(energy, '.1f')+' GeV')
plt.xlabel('-t [GeV^2]')
plt.ylabel('ThetaCM [deg]')
plt.title('-t vs thetaCM')
plt.legend()
plt.grid()
plt.savefig('output/minus_t_thetaCM.png')