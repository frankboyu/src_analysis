import numpy as np
import matplotlib.pyplot as plt
import ROOT as root

mass_proton = 0.938272
mass_neutron = 0.939565
mass_piminus = 0.139570
Deg_To_Rad = np.pi/180

fig, axs = plt.subplots(1, 3, figsize=(15, 5), dpi=200)
axs[0].set_xlabel('Pion momentum [GeV]')
axs[0].set_ylabel('Pion theta [deg]')
axs[0].set_title('Pion momentum vs theta')
axs[0].grid()
axs[1].set_xlabel('Proton momentum [GeV]')
axs[1].set_ylabel('Proton theta [deg]')
axs[1].set_title('Proton momentum vs theta')
axs[1].grid()
axs[2].set_xlabel('Pion theta [deg]')
axs[2].set_ylabel('Proton theta [deg]')
axs[2].set_title('Pion theta vs Proton theta')
axs[2].grid()

for energy in np.arange(6.0, 11.0, 0.5):

    piminus_p = []
    piminus_theta = []
    proton_p = []
    proton_theta = []

    for theta in range(10, 171, 1):
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

        piminus_p.append(P4_piminus.P())
        piminus_theta.append(P4_piminus.Theta()/Deg_To_Rad)
        proton_p.append(P4_proton.P())
        proton_theta.append(P4_proton.Theta()/Deg_To_Rad)

    axs[0].plot(piminus_p, piminus_theta, label=format(energy, '.1f')+' GeV')
    axs[1].plot(proton_p, proton_theta, label=format(energy, '.1f')+' GeV')
    axs[2].plot(piminus_theta, proton_theta, label=format(energy, '.1f')+' GeV')


plt.legend()
plt.savefig('kinematics.png')