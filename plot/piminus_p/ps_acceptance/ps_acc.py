import numpy as np
import matplotlib.pyplot as plt

def PSAcceptance(x, par):
    emin = par[1]
    emax = par[2]

    if (x >= 2*emin and x < (emin + emax)):
        return par[0]*(1-2*emin/x)
    elif (x >= (emin + emax) and x < 2*emax):
        return par[0]*(2*emax/x - 1)
    else:
        return 0

par = [0.840087, 2.92615, 5.97122]
eg = np.linspace(0,13,130)
ps_acc = np.empty(130)
for i in range(130):
    ps_acc[i] = PSAcceptance(eg[i], par)

print("min energy:")
print(par[1]*2)
print("max energy:")
print(par[2]*2)

plt.figure(dpi=200)
plt.plot(eg, ps_acc)
plt.xlabel(r'$E_{\gamma}$ (GeV)')
plt.ylabel('PS Acceptance')
plt.savefig('ps_acc.png')