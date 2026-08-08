import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
from scipy.integrate import quad
import scipy
import ROOT as root

mass_kaon = 0.493677

def exponential(x, a, b, c):
    return np.exp(a*x+b)+c

def double_exponential(x, a, b, c, d, e):
    return np.exp(a*x+b)+np.exp(c*x+d)+e

def gaussian(x, a, b, c):
    return a*np.exp(-0.5*((x-b)/c)**2)

def double_gaussian(x, a1, b1, c1, a2, b2, c2):
    return a1*np.exp(-0.5*((x-b1)/c1)**2) + a2*np.exp(-0.5*((x-b2)/c2)**2)

def triple_gaussian(x, a1, b1, c1, a2, b2, c2, a3, c3):
    return a1*np.exp(-0.5*((x-b1)/c1)**2) + a2*np.exp(-0.5*((x-b2)/c2)**2) + a3*np.exp(-0.5*((x-b2)/c3)**2)

ROOT.gROOT.SetBatch(True)
nObj = 0

class File:
    def __init__(self,infile):
        if (isinstance(infile,ROOT.TFile)):
            self.TFile = infile
        else:
            self.TFile = ROOT.TFile(infile)

    def get(self,name,**kwargs):
        # if (not self.TFile.GetListOfKeys().Contains(name)):
            # raise ValueError("File does not contain specified object")
        h = self.TFile.Get(name)
        if (isinstance(h,ROOT.TH2)):
            return Hist2D(h,**kwargs)
        elif (isinstance(h,ROOT.TH1)):
            return Hist1D(h,**kwargs)
        else:
            raise ValueError("Object is not a supported type")

    def getNames(self):
        return [i.GetName() for i in self.TFile.GetListOfKeys()]

    def plotPoints(self,name,rebin=1,scale=1,**kwargs):
        h = self.get(name,rebin=rebin,scale=scale)
        if (isinstance(h,Hist1D)):
            return h.plotPoints(**kwargs)
        else:
            raise ValueError("This is not a 1D histogram and cannot be plotted with this method")

    def plotBand(self,name,rebin=1,scale=1,**kwargs):
        h = self.get(name,rebin=rebin,scale=scale)
        if (isinstance(h,Hist1D)):
            return h.plotBand(**kwargs)
        else:
            raise ValueError("This is not a 1D histogram and cannot be plotted with this method")

    def plotBar(self,name,rebin=1,scale=1,**kwargs):
        h = self.get(name,rebin=rebin,scale=scale)
        if (isinstance(h,Hist1D)):
            return h.plotBar(**kwargs)
        else:
            raise ValueError("This is not a 1D histogram and cannot be plotted with this method")

    def plotHeatmap(self,name,rebinx=1,rebiny=1,**kwargs):
        h = self.get(name,rebinx=rebinx,rebiny=rebiny)
        if (isinstance(h,Hist2D)):
            return h.plotHeatmap(**kwargs)
        else:
            raise ValueError("This is not a 2D histogram and cannot be plotted with this method")

class Hist1D:
    def __init__(self,hist,rebin=1,scale=1):
        global nObj
        self.TH1 = hist.Clone(str(nObj))
        nObj = nObj + 1

        if (rebin != 1):
            self.TH1.Rebin(rebin)

        g = ROOT.TGraphAsymmErrors(self.TH1)
        self.x = np.array(g.GetX())
        self.y = np.array(g.GetY())*scale
        xerr = []
        yerr = []
        for i in range(g.GetN()):
            xerr.append(g.GetErrorX(i))
            yerr.append(g.GetErrorY(i))
        self.xerr = np.array(xerr)
        self.yerr = np.array(yerr)*scale

    def scale(self,factor):
        self.y *= factor
        self.yerr *= factor

    def rebin(self,factor):
        if (factor != 1):
            self.TH1.Rebin(factor)
            g = ROOT.TGraphAsymmErrors(self.TH1)
            self.x = np.array(g.GetX())
            self.y = np.array(g.GetY())
            xerr = []
            yerr = []
            for i in range(g.GetN()):
                xerr.append(g.GetErrorX(i))
                yerr.append(g.GetErrorY(i))
            self.xerr = np.array(xerr)
            self.yerr = np.array(yerr)

    def areaNorm(self,reference):
        factor = np.sum(reference.y)/np.sum(self.y)
        self.scale(factor)
        return factor

    def plotPoints(self,**kwargs):
        return plt.errorbar(self.x,self.y,yerr=self.yerr,**kwargs)

    def plotBand(self,alpha=0.25,**kwargs):
        line = plt.plot(self.x,self.y,**kwargs)[0]
        band = plt.fill_between(self.x,self.y-self.yerr,self.y+self.yerr,
                                color=line.get_color(),zorder=line.zorder,
                                alpha=alpha)
        return line, band

    def plotBar(self,shift=0,**kwargs):
        bar = plt.bar(self.x+shift, self.y, **kwargs)

        return bar

class Hist2D:
    def __init__(self,hist,rebinx=1,rebiny=1):
        global nObj
        self.TH2 = hist.Clone(str(nObj))
        nObj = nObj + 1

        if (rebinx != 1):
            self.TH2.RebinX(rebinx)
        if (rebiny != 1):
            self.TH2.RebinY(rebiny)

        NX = self.TH2.GetNbinsX()
        NY = self.TH2.GetNbinsY()

        xedge = []
        yedge = []
        z = []
        zerr = []

        for j in range(NX):
            xedge.append(self.TH2.GetXaxis().GetBinLowEdge(j+1))
        xedge.append(self.TH2.GetXaxis().GetBinUpEdge(NX))

        for i in range(NY):
            yedge.append(self.TH2.GetYaxis().GetBinLowEdge(i+1))
            zcol = []
            zerrcol = []
            for j in range(NX):
                zval = self.TH2.GetBinContent(j+1,i+1)
                zcol.append(zval)
                zerrcol.append(self.TH2.GetBinError(j+1,i+1))
            z.append(zcol)
            zerr.append(zerrcol)
        yedge.append(self.TH2.GetYaxis().GetBinUpEdge(NY))

        self.xedge = np.array(xedge)
        self.yedge = np.array(yedge)
        self.z = np.array(z)
        self.zerr = np.array(zerr)

    def scale(self,factor):
        self.z *= factor
        self.zerr *= factor

    def areaNorm(self,reference):
        factor = np.sum(reference.z)/np.sum(self.z)
        self.scale(factor)
        return factor

    def plotHeatmap(self,kill_zeros=True,**kwargs):
        z = self.z
        if (kill_zeros):
            z[z==0] = np.nan
        return plt.pcolormesh(self.xedge,self.yedge,self.z,**kwargs)

    def projectionX(self,**kwargs):
        return Hist1D(self.TH2.ProjectionX(),**kwargs)

    def projectionY(self,**kwargs):
        return Hist1D(self.TH2.ProjectionY(),**kwargs)

print("################################################################# READ THE FILES #################################################################")

data_version = 'ver12'
sim_version  = '07'

print("Data version: " + data_version)
print("Sim version: " + sim_version)

file_pdf    = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_phi_d_hists.pdf")
file_data   = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_recon_data_"    + data_version                     + ".root")
file_sim    = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_recon_sim_"     + data_version + "_" + sim_version + ".root")
# file_tagged = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_thrown_tagged_" + data_version + "_" + sim_version + ".root")
# file_gen    = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_thrown_gen_"    + data_version + "_" + sim_version + ".root")

print("################################################################# EVENT SELECTION #################################################################")

print("KinFit Chi2")
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('KinFitChiSqCut/kinfit_cut_chisq_per_ndf_KinFitChiSqCut')
hist_sim  = file_sim. get('KinFitChiSqCut/kinfit_cut_chisq_per_ndf_KinFitChiSqCut')
hist_sim.scale(hist_data.y.max()/hist_sim.y.max())
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim. plotPoints(label="Sim",  marker='o', markersize=3, linestyle='None', color='blue')
plt.plot(np.linspace(0,10,1000), hist_data.y.max()*(scipy.stats.chi2.pdf(np.linspace(0,10,1000)*7, df=7)*7)/(scipy.stats.chi2.pdf(np.linspace(0,10,1000)*7, df=7)*7).max(), '-', color='green', label=r'$\chi^2$(ndf=7) PDF')
plt.plot([5, 5], [0, hist_data.y.max()*1.2], '-', color='red', label=r'Cut position')
plt.xlabel(r"$\chi^{2}$/NDF")
plt.ylabel(r"Counts")
plt.xlim(0, 10)
plt.ylim(0, round(hist_data.y.max()*1.2, -2))
plt.legend()
file_pdf.savefig()
plt.close()

print("Deuteron dE/dx: comparison between data and simulation")
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('dEdxCut/pid_cut_d_dEdx_cdc_dEdxCut')
hist_sim  = file_sim .get('dEdxCut/pid_cut_d_dEdx_cdc_dEdxCut')
# hist_sim.scale(0.002)
hist_sim.areaNorm(hist_data)
plt.axes(axs[0])
hist_data.plotHeatmap(vmin=0, vmax=50)
plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
plt.axes(axs[1])
hist_sim.plotHeatmap(vmin=0, vmax=50)
plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
for i in range(2):
    plt.axes(axs[i])
    plt.plot(np.arange(0.25, 3, 0.01), np.exp(-3.65*np.arange(0.25, 3, 0.01)+4.47)+2.57, color='red')
    plt.xlim(0, 2)
    plt.ylim(0, 40)
    plt.xlabel(r"$p_d$ (GeV/c)")
    plt.ylabel(r"$(dE/dx)^{\mathrm{CDC}}_d$ (keV/cm)")
file_pdf.savefig()
plt.close()

print("Deuteron dE/dx: combined fit of deuteron and proton bands")
plt.figure(figsize=(16,20), dpi=300)
fig, axs = plt.subplots(5, 3, figsize=(16,20), dpi=300)
axs = axs.flatten()
hist_data = file_data.get('dEdxCut/pid_cut_d_dEdx_cdc_dEdxCut')
rebin_factors  = np.array([4,  4,  4,    2,  1,  1,  1,  2,  2,  4,  4,   4,   8,   8], dtype=int)
dedx_p_edges   = np.array([45, 50, 52.5, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 200], dtype=int)
dedx_p_low     = dedx_p_edges[:-1]
dedx_p_high    = dedx_p_edges[1:]
dedx_p_centers = (dedx_p_low + dedx_p_high)/2
d_amplitude_value,  d_amplitude_err,  d_mean_value,   d_mean_err, d_sigma_value,  d_sigma_err  = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
p_amplitude1_value, p_amplitude1_err, p_mean_value,   p_mean_err, p_sigma1_value, p_sigma1_err = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
p_amplitude2_value, p_amplitude2_err, p_sigma2_value, p_sigma2_err                             = np.array([]), np.array([]), np.array([]), np.array([])
chisquared = np.array([])
contamination = np.array([])
for i in range(len(dedx_p_edges)-1):
    plt.sca(axs[i])
    hist_slice = Hist1D(hist_data.TH2.ProjectionY("_px",int(dedx_p_edges[i]),int(dedx_p_edges[i+1])))
    hist_slice.rebin(int(rebin_factors[i]))
    hist_slice.yerr = np.sqrt(hist_slice.yerr**2 + 1)  # Add a constant error of 1 to all bins to avoid zero error bars
    hist_slice.plotPoints(fmt='o', label='Data', markersize=3, color='black')
    plt.ylim(0, hist_slice.y.max()*1.2)
    plt.xlim(0, 40)

    fit_low = np.exp(-2.35*(dedx_p_centers[i]/100)+2.71)+0.24       # -2 sigma position for proton from the first trial fit, nominal position
    if (i == 0):
        fit_low = np.exp(-2.88*(dedx_p_centers[i]/100)+3.07)+0.66   # -1 sigma position for proton from the first trial fit, as the proton peak is not perfectly Gaussian at low momentum
    if (i == len(dedx_p_edges)-2):
        fit_low = np.exp(-3.42*(dedx_p_centers[i]/100)+4.33)+2.39   # -2 sigma position for deuteron from the first trial fit, last momentum bin where the proton peak is negligible
    fit_high    = np.exp(-4.64*(dedx_p_centers[i]/100)+5.70)+5.22   # +2 sigma position for deuteron from the first trial fit
    p0_d_mean   = np.exp(-4.22*(dedx_p_centers[i]/100)+5.18)+3.87
    p0_d_sigma  = np.exp(-4.46*(dedx_p_centers[i]/100)+5.47)+4.55 - p0_d_mean
    p0_p_mean   = np.exp(-3.62*(dedx_p_centers[i]/100)+3.70)+1.48
    p0_p_sigma  = np.exp(-4.00*(dedx_p_centers[i]/100)+4.06)+1.99 - p0_p_mean
    fit_mask = (hist_slice.x > fit_low) & (hist_slice.x < fit_high)
    hist_slice.x    = hist_slice.x[fit_mask]
    hist_slice.y    = hist_slice.y[fit_mask]
    hist_slice.yerr = hist_slice.yerr[fit_mask]
    if (i < len(dedx_p_edges)-2):
        popt, pcov = curve_fit(triple_gaussian, hist_slice.x, hist_slice.y, \
                                p0=[hist_slice.y.max()/2, p0_d_mean, p0_d_sigma, hist_slice.y.max()/2, p0_p_mean, p0_p_sigma, hist_slice.y.max()/2, p0_p_sigma], \
                                bounds=([0, 0, 0, 0, 0, 0, 0, 0], \
                                [hist_slice.y.max(), p0_d_mean+2*p0_d_sigma, 2*p0_d_sigma, hist_slice.y.max(), p0_p_mean+2*p0_p_sigma, 2*p0_p_sigma, hist_slice.y.max(), 2*p0_p_sigma]))
        residuals           = hist_slice.y - triple_gaussian(hist_slice.x, *popt)
        chisq               = np.sum((residuals**2) / triple_gaussian(hist_slice.x, *popt))
        d_amplitude_value   = np.append(d_amplitude_value,  popt[0])
        d_amplitude_err     = np.append(d_amplitude_err,    np.sqrt(pcov[0][0]))
        d_mean_value        = np.append(d_mean_value,       popt[1])
        d_mean_err          = np.append(d_mean_err,         np.sqrt(pcov[1][1]))
        d_sigma_value       = np.append(d_sigma_value,      abs(popt[2]))
        d_sigma_err         = np.append(d_sigma_err,        np.sqrt(pcov[2][2]))
        p_amplitude1_value  = np.append(p_amplitude1_value, popt[3])
        p_amplitude1_err    = np.append(p_amplitude1_err,   np.sqrt(pcov[3][3]))
        p_mean_value        = np.append(p_mean_value,       popt[4])
        p_mean_err          = np.append(p_mean_err,         np.sqrt(pcov[4][4]))
        p_sigma1_value      = np.append(p_sigma1_value,     abs(popt[5]))
        p_sigma1_err        = np.append(p_sigma1_err,       np.sqrt(pcov[5][5]))
        p_amplitude2_value  = np.append(p_amplitude2_value, popt[6])
        p_amplitude2_err    = np.append(p_amplitude2_err,   np.sqrt(pcov[6][6]))
        p_sigma2_value      = np.append(p_sigma2_value,     abs(popt[7]))
        p_sigma2_err        = np.append(p_sigma2_err,       np.sqrt(pcov[7][7]))
        chisquared          = np.append(chisquared,         chisq)
        contamination       = np.append(contamination, quad(lambda x: double_gaussian(x, popt[3], popt[4], popt[5], popt[6], popt[4], popt[7]), popt[1]-2.0*popt[2], np.inf)[0] / quad(lambda x: triple_gaussian(x, *popt), popt[1]-2.0*popt[2], np.inf)[0])
        plt.plot(np.arange(fit_low, fit_high, 0.01), triple_gaussian(np.arange(fit_low, fit_high, 0.01), *popt),                label="Triple Gaussian Fit",    color='red')
        plt.plot(np.arange(fit_low, fit_high, 0.01), gaussian(np.arange(fit_low, fit_high, 0.01), popt[0], popt[1], popt[2]),   label="Deuteron Component",     color='blue',   linestyle='--')
        plt.plot(np.arange(fit_low, fit_high, 0.01), double_gaussian(np.arange(fit_low, fit_high, 0.01), popt[3], popt[4], popt[5], popt[6], popt[4], popt[7]),   label="Proton Component",       color='green',  linestyle='--')
        if (i == 3):
            plt.text(0.5, 0.50, r'$A_1=%.2f \pm %.2f$' %        (popt[0], np.sqrt(pcov[0][0])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.45, r'$\mu_1=%.3f \pm %.3f$' %      (popt[1], np.sqrt(pcov[1][1])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.40, r'$\sigma_1=%.3f \pm %.3f$' %   (abs(popt[2]), np.sqrt(pcov[2][2])),    transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.35, r'$A_2=%.2f \pm %.2f$' %        (popt[3], np.sqrt(pcov[3][3])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.30, r'$\mu_2=%.3f \pm %.3f$' %      (popt[4], np.sqrt(pcov[4][4])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.25, r'$\sigma_2=%.3f \pm %.3f$' %   (abs(popt[5]), np.sqrt(pcov[5][5])),    transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.20, r'$\chi^2/ndf=%.2f$' %          (chisq/(len(hist_slice.x)-3)),          transform=plt.gca().transAxes, fontsize=10)
        else:
            plt.text(0.5, 0.70, r'$A_1=%.2f \pm %.2f$' %        (popt[0], np.sqrt(pcov[0][0])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.65, r'$\mu_1=%.3f \pm %.3f$' %      (popt[1], np.sqrt(pcov[1][1])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.60, r'$\sigma_1=%.3f \pm %.3f$' %   (abs(popt[2]), np.sqrt(pcov[2][2])),    transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.55, r'$A_2=%.2f \pm %.2f$' %        (popt[3], np.sqrt(pcov[3][3])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.50, r'$\mu_2=%.3f \pm %.3f$' %      (popt[4], np.sqrt(pcov[4][4])),         transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.45, r'$\sigma_2=%.3f \pm %.3f$' %   (abs(popt[5]), np.sqrt(pcov[5][5])),    transform=plt.gca().transAxes, fontsize=10)
            plt.text(0.5, 0.40, r'$\chi^2/ndf=%.2f$' %          (chisq/(len(hist_slice.x)-3)),          transform=plt.gca().transAxes, fontsize=10)
        plt.title(r'$p$ bin: [%.2f, %.2f] GeV/c' %  (dedx_p_low[i]/100, dedx_p_high[i]/100))
    else:
        popt, pcov          = curve_fit(gaussian, hist_slice.x, hist_slice.y, p0=[hist_slice.y.max(),p0_d_mean,p0_d_sigma])
        residuals           = hist_slice.y - gaussian(hist_slice.x, *popt)
        chisq               = np.sum((residuals**2) / gaussian(hist_slice.x, *popt))
        d_amplitude_value   = np.append(d_amplitude_value,  popt[0])
        d_amplitude_err     = np.append(d_amplitude_err,    np.sqrt(pcov[0][0]))
        d_mean_value        = np.append(d_mean_value,       popt[1])
        d_mean_err          = np.append(d_mean_err,         np.sqrt(pcov[1][1]))
        d_sigma_value       = np.append(d_sigma_value,      abs(popt[2]))
        d_sigma_err         = np.append(d_sigma_err,        np.sqrt(pcov[2][2]))
        chisquared          = np.append(chisquared,         chisq)
        contamination       = np.append(contamination,      0)
        plt.plot(np.arange(fit_low, fit_high, 0.01), gaussian(np.arange(fit_low, fit_high, 0.01), *popt), label="Gaussian Fit", color='red')
        plt.title(r'$p$ bin: [%.2f, %.2f] GeV/c' % (dedx_p_low[i]/100, dedx_p_high[i]/100))
        plt.text(0.5, 0.90, r'$A_2=%.2f \pm %.2f$' %            (popt[0], np.sqrt(pcov[0][0])),         transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.5, 0.85, r'$\mu_2=%.3f \pm %.3f$' %          (popt[1], np.sqrt(pcov[1][1])),         transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.5, 0.80, r'$\sigma_2=%.3f \pm %.3f$' %       (abs(popt[2]), np.sqrt(pcov[2][2])),    transform=plt.gca().transAxes, fontsize=10)
        plt.text(0.5, 0.75, r'$\chi^2/ndf=%.2f$' %              (chisq/(len(hist_slice.x)-3)),          transform=plt.gca().transAxes, fontsize=10)
fig.supxlabel(r'$dE/dx$ (keV/cm)', fontsize=20)
fig.supylabel('Counts', fontsize=20)
axs[12].legend()
file_pdf.savefig()
plt.close()

print("deuteron CDC dE/dx: average proton contamination")
plt.figure(figsize=(8,6))
plt.plot(dedx_p_centers/100, contamination, 'o-', color='magenta')
plt.xlabel('Deuteron Momentum (GeV/c)')
plt.ylabel('Average Proton Contamination')
plt.ylim(0, 0.005)
file_pdf.savefig()
plt.close()

print("deuteron CDC dE/dx: individual proton contaimations")
plt.figure(figsize=(8,6))
avg_purity, avg_efficiency, avg_rejection = np.array([]), np.array([]), np.array([])
for cut_values in np.arange(-3, -1, 0.25):
    total_signal, total_background, signal, background              = 0, 0, 0, 0
    individual_purity, individual_efficiency, individual_rejection  = np.array([]), np.array([]), np.array([])
    for i in range(len(dedx_p_edges)-2):
        total_signal        += quad(lambda x: gaussian(       x, d_amplitude_value[i],  d_mean_value[i], d_sigma_value[i]),                                                             -np.inf,                                     np.inf)[0]
        total_background    += quad(lambda x: double_gaussian(x, p_amplitude1_value[i], p_mean_value[i], p_sigma1_value[i], p_amplitude2_value[i], p_mean_value[i], p_sigma2_value[i]), -np.inf,                                     np.inf)[0]
        signal              += quad(lambda x: gaussian(       x, d_amplitude_value[i],  d_mean_value[i], d_sigma_value[i]),                                                             d_mean_value[i]+cut_values*d_sigma_value[i], np.inf)[0]
        background          += quad(lambda x: double_gaussian(x, p_amplitude1_value[i], p_mean_value[i], p_sigma1_value[i], p_amplitude2_value[i], p_mean_value[i], p_sigma2_value[i]), d_mean_value[i]+cut_values*d_sigma_value[i], np.inf)[0]
        individual_purity     = np.append(individual_purity,     signal/(signal+background))
        individual_efficiency = np.append(individual_efficiency, signal/total_signal)
        individual_rejection  = np.append(individual_rejection,  1-background/total_background)
    plt.plot(dedx_p_centers[:-1]/100, 1-individual_purity, 'o-', label='Cut at %.2f sigma' % cut_values)
    total_signal += quad(lambda x: gaussian(x, d_amplitude_value[-1], d_mean_value[-1], d_sigma_value[-1]), -np.inf, np.inf)[0]
    signal       += quad(lambda x: gaussian(x, d_amplitude_value[-1], d_mean_value[-1], d_sigma_value[-1]), d_mean_value[-1]+cut_values*d_sigma_value[-1], np.inf)[0]
    avg_purity     = np.append(avg_purity, signal/(signal+background))
    avg_efficiency = np.append(avg_efficiency, signal/total_signal)
    avg_rejection  = np.append(avg_rejection, 1-background/total_background)
plt.legend()
plt.xlabel('Deuteron Momentum (GeV/c)')
plt.ylabel('Average Proton Contamination')
plt.ylim(0, 0.005)
file_pdf.savefig()
plt.close()

print("deuteron CDC dE/dx: perfomance metrics")
plt.figure(figsize=(8,6))
plt.plot(np.arange(-3, -1, 0.25), avg_purity,     'o-', color='magenta', label='Purity')
plt.plot(np.arange(-3, -1, 0.25), avg_efficiency, 'o-', color='cyan',    label='Efficiency')
plt.plot(np.arange(-3, -1, 0.25), avg_rejection,  'o-', color='yellow',  label='Rejection')
plt.xlabel('Cut Values')
plt.ylabel('Performance Metrics')
plt.legend()
file_pdf.savefig()
plt.close()

# # hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
# # p_points = np.arange(0.25, 3, 0.01)
# # para_list = []
# # variation_list = np.array([-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, 0.0, 1.0, 2.0, 3.0])
# # color_list = ['violet', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'violet']
# # for N in variation_list:
# #     popt, pcov = curve_fit(exponential, dedx_p_centers/100, d_mean_value+N*d_sigma_value)
# #     para_list.append(popt)

# # plt.figure(figsize=(8,6))
# # plt.errorbar(dedx_p_centers/100, d_mean_value, xerr=(dedx_p_high-dedx_p_low)/200, yerr=d_mean_err, fmt='k.', label=r'$\mu$')
# # plt.errorbar(dedx_p_centers/100, d_sigma_value, xerr=(dedx_p_high-dedx_p_low)/200, yerr=d_sigma_err, fmt='b.', label=r'$\sigma$')
# # plt.xlim(0,2)
# # plt.ylim(0,40)
# # plt.xlabel(r'$p$ (GeV/c)')
# # plt.ylabel(r'$dE/dx$ (keV/cm)')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# # plt.figure(figsize=(8,6))
# # hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
# # dE_points = exponential(p_points, *para_list[4])
# # plt.scatter(dedx_p_centers/100, d_mean_value-2.0*d_sigma_value, color='r', marker='o', s=20, label=r'Data points of $\mu-2\sigma$')
# # plt.plot(p_points, dE_points, label=r'Exponential fit, $p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[4][0], para_list[4][1], para_list[4][2]), color='y')
# # plt.xlim(0,2)
# # plt.ylim(0,40)
# # plt.xlabel(r'$p$ (GeV/c)')
# # plt.ylabel(r'$dE/dx$ (keV/cm)')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# # plt.figure(figsize=(8,6))
# # hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
# # for i in range(len(variation_list)):
# #     N = variation_list[i]
# #     dE_points = exponential(p_points, *para_list[i])
# #     plt.plot(p_points, dE_points, label=r'$%.2f \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (N, para_list[i][0], para_list[i][1], para_list[i][2]), color=color_list[i])
# #     plt.scatter(dedx_p_centers/100, d_mean_value+N*d_sigma_value, color=color_list[i], marker='o', s=10)
# # plt.xlim(0,2)
# # plt.ylim(0,40)
# # plt.xlabel(r'$p$ (GeV/c)')
# # plt.ylabel(r'$dE/dx$ (keV/cm)')
# # plt.legend()
# # file_pdf.savefig()
# # plt.close()

# print("deuteron CDC dE/dx: fit the deuteron band only")
# plt.figure(figsize=(16,12))
# fig, axs = plt.subplots(4, 4, figsize=(20,16))
# axs = axs.flatten()
# hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
# dedx_p_edges = np.array([45, 50, 52.5, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 200], dtype=int)
# dedx_p_low = dedx_p_edges[:-1]
# dedx_p_high = dedx_p_edges[1:]
# dedx_p_centers = (dedx_p_low + dedx_p_high)/2
# d_mean_value, d_mean_err, d_sigma_value, d_sigma_err = np.array([]), np.array([]), np.array([]), np.array([])
# chisquared = np.array([])
# for i in range(len(dedx_p_edges)-1):
#     plt.sca(axs[i])
#     hist_slice = Hist1D(hist_data.TH2.ProjectionY("_px",int(dedx_p_edges[i]),int(dedx_p_edges[i+1])))
#     if (i < 1):
#         hist_slice.rebin(4)
#     elif (i < 8):
#         hist_slice.rebin(2)
#     # if (i < 2):
#         # hist_slice.rebin(2)
#     hist_slice.plotPoints(fmt='o', label='Data', markersize=3)
#     plt.ylim(0, hist_slice.y.max()*1.2)
#     plt.xlim(0, 40)
#     fit_low = np.exp(-3.67*(dedx_p_centers[i]/100)+4.48)+2.57
#     fit_high = np.exp(-4.58*(dedx_p_centers[i]/100)+5.66)+5.22
#     fit_mask = (hist_slice.x > fit_low) & (hist_slice.x < fit_high)
#     hist_slice.x = hist_slice.x[fit_mask]
#     hist_slice.y = hist_slice.y[fit_mask]
#     hist_slice.yerr = hist_slice.yerr[fit_mask]
#     popt, pcov = curve_fit(gaussian, hist_slice.x, hist_slice.y, p0=[hist_slice.y.max(),(fit_low+fit_high)/2,(fit_high-fit_low)/4])
#     d_mean_value = np.append(d_mean_value, popt[1])
#     d_mean_err = np.append(d_mean_err, np.sqrt(pcov[1][1]))
#     d_sigma_value = np.append(d_sigma_value, abs(popt[2]))
#     d_sigma_err = np.append(d_sigma_err, np.sqrt(pcov[2][2]))
#     residuals = hist_slice.y - gaussian(hist_slice.x, *popt)
#     chisq = np.sum((residuals**2) / gaussian(hist_slice.x, *popt))
#     chisquared = np.append(chisquared, chisq)
#     plt.plot(np.arange(fit_low, fit_high, 0.01), gaussian(np.arange(fit_low, fit_high, 0.01), *popt), \
#             label="Gaussian Fit", color='red', zorder=10)
#     plt.text(0.4, 0.75, r'$A=%.2f \pm %.2f$' % (popt[0], np.sqrt(pcov[0][0])), transform=plt.gca().transAxes, fontsize=12)
#     plt.text(0.4, 0.65, r'$\mu=%.3f \pm %.3f$' % (popt[1], np.sqrt(pcov[1][1])), transform=plt.gca().transAxes, fontsize=12)
#     plt.text(0.4, 0.55, r'$\sigma=%.3f \pm %.3f$' % (abs(popt[2]), np.sqrt(pcov[2][2])), transform=plt.gca().transAxes, fontsize=12)
#     plt.text(0.4, 0.45, r'$\chi^2/ndf=%.2f$' % (chisq/(len(hist_slice.x)-3)), transform=plt.gca().transAxes, fontsize=12)
#     plt.title(r'$p$ bin: [%.2f, %.2f] GeV/c' % (dedx_p_low[i]/100, dedx_p_high[i]/100))
# fig.supxlabel(r'$dE/dx$ (keV/cm)', fontsize=20)
# fig.supylabel('Counts', fontsize=20)
# axs[0].legend()
# file_pdf.savefig()
# plt.close()

# print("deuteron CDC dE/dx: chosen cut value")
# plt.figure(figsize=(8,6))
# p_points = np.arange(0.25, 3, 0.01)
# para_list = []
# variation_list = np.array([-3, -2.75, -2.5, -2.25, -2, -1.75, -1.5, -1.25, -1, 0.0, 1.0, 2.0, 3.0])
# color_list = ['violet', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'violet']
# for N in variation_list:
#     popt, pcov = curve_fit(exponential, dedx_p_centers/100, d_mean_value+N*d_sigma_value)
#     para_list.append(popt)
# hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
# for i in range(len(variation_list)):
#     N = variation_list[i]
#     dE_points = exponential(p_points, *para_list[i])
#     plt.plot(p_points, dE_points, \
#             label=r'$%.2f \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (N, para_list[i][0], para_list[i][1], para_list[i][2]), \
#             color=color_list[i])
#     plt.scatter(dedx_p_centers/100, d_mean_value+N*d_sigma_value, color=color_list[i], marker='o', s=10)
# plt.xlim(0,2)
# plt.ylim(0,40)
# plt.xlabel(r'$p$ (GeV/c)')
# plt.ylabel(r'$dE/dx$ (keV/cm)')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# print("deuteron CDC dE/dx: various cut options")
# plt.figure(figsize=(8,6))
# hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
# for i in range(len(variation_list)):
#     N = variation_list[i]
#     dE_points = exponential(p_points, *para_list[i])
#     plt.plot(p_points, dE_points, \
#             label=r'$%.2f \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (N, para_list[i][0], para_list[i][1], para_list[i][2]), \
#             color=color_list[i])
#     plt.scatter(dedx_p_centers/100, d_mean_value+N*d_sigma_value, color=color_list[i], marker='o', s=10)
# plt.xlim(0,2)
# plt.ylim(0,40)
# plt.xlabel(r'$p$ (GeV/c)')
# plt.ylabel(r'$dE/dx$ (keV/cm)')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# print("deuteron CDC dE/dx: chosen cut option")
# plt.figure(figsize=(8,6))
# hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
# dE_points = exponential(p_points, *para_list[4])
# plt.scatter(dedx_p_centers/100, d_mean_value-2.0*d_sigma_value, color='r', marker='o', s=20, label=r'Data points of $\mu-2\sigma$')
# plt.plot(p_points, dE_points, label=r'Exponential fit, $p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[4][0], para_list[4][1], para_list[4][2]), color='y')
# plt.xlim(0,2)
# plt.ylim(0,40)
# plt.xlabel(r'$p$ (GeV/c)')
# plt.ylabel(r'$dE/dx$ (keV/cm)')
# plt.legend()
# file_pdf.savefig()
# plt.close()

print("Exclusivity: comparison between data and simulation")
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('MissPMinusCut/exclusivity_cut_miss_pminus_MissPMinusCut')
hist_sim  = file_sim .get('MissPMinusCut/exclusivity_cut_miss_pminus_MissPMinusCut')
hist_sim.scale(hist_data.y.max()/hist_sim.y.max())
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
plt.plot([-0.02, -0.02], [0, 700], '-', color='red')
plt.xlim(-0.2, 0.2)
plt.ylim(0, 700)
plt.xlabel(r"$p^-_{\rm miss}$ (GeV/c)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

print("Exclusivity: comparison between data and simulation as pion")
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('MissPMinusCut/exclusivity_miss_pminus_as_pion_MissPMinusCut')
hist_sim  = file_sim .get('MissPMinusCut/exclusivity_miss_pminus_as_pion_MissPMinusCut')
hist_sim.scale(hist_data.y.max()/hist_sim.y.max())
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
plt.xlim(-0.2, 0.2)
plt.ylim(0, 700)
plt.xlabel(r"$p^-_{\rm miss}$ (GeV/c)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

print("Exclusivity: fit data with double Gaussian")
def double_gaussian(x, A1, mu1, sigma1, A2, mu2, sigma2):
    return A1 * np.exp(-(x - mu1)**2 / (2 * sigma1**2)) + A2 * np.exp(-(x - mu2)**2 / (2 * sigma2**2))
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('MissPMinusCut/exclusivity_cut_miss_pminus_MissPMinusCut')
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
fit_max = 0.02
fit_min = -0.1
fit_mask = (hist_data.x < fit_max) & (hist_data.x > fit_min)
popt, pcov = curve_fit(double_gaussian, hist_data.x[fit_mask], hist_data.y[fit_mask], p0=[150, -0.04, 0.01, 600, 0.0, 0.01], bounds=([0, -0.07, 0, 0, -0.01, 0], [np.inf, -0.04, 0.02, np.inf, 0.01, 0.02]))
chi_square = np.sum(((hist_data.y[fit_mask] - double_gaussian(hist_data.x[fit_mask], *popt)) ** 2) / hist_data.yerr[fit_mask]**2)
ndf = np.sum(fit_mask) - len(popt)
plt.plot(np.linspace(fit_min, fit_max, 1000), double_gaussian(np.linspace(fit_min, fit_max, 1000), *popt), '-', color='red', label='Total')
plt.plot(np.linspace(fit_min, fit_max, 1000), gaussian(np.linspace(fit_min, fit_max, 1000), popt[0], popt[1], popt[2]), '--', color='blue', label='Background')
plt.plot(np.linspace(fit_min, fit_max, 1000), gaussian(np.linspace(fit_min, fit_max, 1000), popt[3], popt[4], popt[5]), '--', color='green', label='Signal')
plt.text(-0.19, 600, f"$A_1$: {popt[0]:.4f} $\pm$ {np.sqrt(pcov[0, 0]):.4f}\n$\mu_1$: {1000*popt[1]:.4f} $\pm$ {1000*np.sqrt(pcov[1, 1]):.4f} MeV/c\n$\sigma_1$: {1000*popt[2]:.4f} $\pm$ {1000*np.sqrt(pcov[2, 2]):.4f} MeV/c", color='blue')
plt.text(-0.19, 500, f"$A_2$: {popt[3]:.4f} $\pm$ {np.sqrt(pcov[3, 3]):.4f}\n$\mu_2$: {1000*popt[4]:.4f} $\pm$ {1000*np.sqrt(pcov[4, 4]):.4f} MeV/c\n$\sigma_2$: {1000*popt[5]:.4f} $\pm$ {1000*np.sqrt(pcov[5, 5]):.4f} MeV/c", color='green')
plt.text(-0.19, 400, f"Chi-square/ndf: {chi_square/ndf:.2f}", color='red')
plt.xlim(-0.2, 0.2)
plt.ylim(0, 700)
plt.xlabel(r"$p^-_{\rm miss}$ (GeV/c)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

print("Exclusivity: signal significance")
fig = plt.figure(figsize=(8, 6), dpi=300)
cut_values      = np.linspace(-0.1, 0.0, 100)
significance    = np.zeros_like(cut_values)
efficiency      = np.zeros_like(cut_values)
purity          = np.zeros_like(cut_values)
rejection       = np.zeros_like(cut_values)
for i, cut in enumerate(cut_values):
    signal              = abs(quad(lambda x: gaussian(x, popt[3], popt[4], popt[5]), cut,     np.inf)[0])
    total_signal        = abs(quad(lambda x: gaussian(x, popt[3], popt[4], popt[5]), -np.inf, np.inf)[0])
    total_background    = abs(quad(lambda x: gaussian(x, popt[0], popt[1], popt[2]), -np.inf, np.inf)[0])
    background          = abs(quad(lambda x: gaussian(x, popt[0], popt[1], popt[2]), cut,     np.inf)[0])
    significance[i]     = signal/np.sqrt(signal+background)
    efficiency[i]       = signal/total_signal
    purity[i]           = signal/(signal+background)
    rejection[i]        = 1 - background/total_background
plt.plot(cut_values, significance, '-o', color='red')
plt.text(0.1, 0.2, f"Maximum significance: {np.max(significance):.2f} at cut value: {cut_values[np.argmax(significance)]:.4f} GeV/c", \
            color='black', fontsize=10, ha='left', va='bottom', transform=plt.gca().transAxes)
plt.text(0.1, 0.1, f"With efficiency of {efficiency[np.argmax(significance)]:.4f}, purity of {purity[np.argmax(significance)]:.4f}, and rejection of {rejection[np.argmax(significance)]:.4f}", \
            color='black', fontsize=10, ha='left', va='bottom', transform=plt.gca().transAxes)
plt.xlabel(r"$p^-_{\rm miss}$ Cut (GeV/c)")
plt.ylabel("Significance (S/sqrt(S+B))")
plt.title("Significance vs. Cut on Missing p minus")
plt.grid()
file_pdf.savefig()
plt.close()

print("Exclusivity: performance metrics")
fig = plt.figure(figsize=(8, 6), dpi=300)
plt.plot(cut_values, efficiency,    '-o', color='blue',  label='Efficiency')
plt.plot(cut_values, purity,        '-o', color='green', label='Purity')
plt.plot(cut_values, rejection,     '-o', color='red',   label='Rejection')
plt.xlabel(r"$p^-_{\rm miss}$ Cut (GeV/c)")
plt.ylabel("Efficiency / Purity / Rejection")
plt.title("Efficiency and Purity vs. Cut on Missing p minus")
plt.legend()
plt.grid()
file_pdf.savefig()
plt.close()

print("Kinematics: K+")
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('KinematicsCut/kinematics_cut_kp_KinematicsCut')
hist_sim  = file_sim .get('KinematicsCut/kinematics_cut_kp_KinematicsCut')
hist_sim.areaNorm(hist_data)
plt.axes(axs[0])
hist_data.plotHeatmap()
plt.text(0.95, 0.95, "data",       transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
plt.axes(axs[1])
hist_sim.plotHeatmap()
plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
for i in range(2):
    plt.axes(axs[i])
    plt.plot([0,10],    [2,2],  '-', color='red')
    plt.plot([0.4,0.4], [0,20], '-', color='red')
    plt.xlim(0, 10)
    plt.ylim(0, 20)
    plt.yticks(np.arange(0, 21, 2))
    plt.xlabel(r"$p_{K^+}$ (GeV/c)")
    plt.ylabel(r"$\theta_{K^+}$ (deg)")
file_pdf.savefig()
plt.close()

print("Kinematics: K-")
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('KinematicsCut/kinematics_cut_km_KinematicsCut')
hist_sim  = file_sim .get('KinematicsCut/kinematics_cut_km_KinematicsCut')
hist_sim.areaNorm(hist_data)
plt.axes(axs[0])
hist_data.plotHeatmap()
plt.text(0.95, 0.95, "data",       transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
plt.axes(axs[1])
hist_sim.plotHeatmap()
plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
for i in range(2):
    plt.axes(axs[i])
    plt.plot([0,10],    [2,2],  '-', color='red')
    plt.plot([0.4,0.4], [0,20], '-', color='red')
    plt.xlim(0, 10)
    plt.ylim(0, 20)
    plt.yticks(np.arange(0, 21, 2))
    plt.xlabel(r"$p_{K^-}$ (GeV/c)")
    plt.ylabel(r"$\theta_{K^-}$ (deg)")
file_pdf.savefig()
plt.close()

print("Kinematics: deuteron")
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('KinematicsCut/kinematics_cut_d_KinematicsCut')
hist_sim  = file_sim .get('KinematicsCut/kinematics_cut_d_KinematicsCut')
hist_sim.areaNorm(hist_data)
plt.axes(axs[0])
hist_data.plotHeatmap()
plt.text(0.95, 0.95, "data",       transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
plt.axes(axs[1])
hist_sim.plotHeatmap()
plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
for i in range(2):
    plt.axes(axs[i])
    plt.plot([0,10],    [2,2],   '-', color='red')
    plt.plot([0.4,0.4], [0,180], '-', color='red')
    plt.xlim(0, 2)
    plt.ylim(0, 180)
    plt.xlabel(r"$p_{d}$ (GeV/c)")
    plt.ylabel(r"$\theta_{d}$ (deg)")
file_pdf.savefig()
plt.close()

print("Vertex: x-y")
fig = plt.figure(figsize=(16, 7), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('VertexCut/vertex_cut_x_y_VertexCut')
hist_sim  = file_sim .get('VertexCut/vertex_cut_x_y_VertexCut')
plt.axes(axs[0])
hist_data.plotHeatmap(vmin=0, vmax=100)
plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), '-', color='red')
plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
plt.axes(axs[1])
hist_sim.plotHeatmap(vmin=0, vmax=100)
plt.plot(np.cos(np.linspace(0, 2*np.pi, 100)), np.sin(np.linspace(0, 2*np.pi, 100)), '-', color='red')
plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
for i in range(2):
    plt.axes(axs[i])
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.xlabel(r"Vertex $x$ (cm)")
    plt.ylabel(r"Vertex $y$ (cm)")
file_pdf.savefig()
plt.close()

print("Vertex: z")
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('VertexCut/vertex_cut_z_VertexCut')
hist_sim  = file_sim .get('VertexCut/vertex_cut_z_VertexCut')
# hist_sim.scale(0.002)
hist_sim.areaNorm(hist_data)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim .plotPoints(label="Sim",  marker='o', markersize=3, linestyle='None', color='blue')
plt.plot([51, 51], [0, 400], '-', color='red')
plt.plot([79, 79], [0, 400], '-', color='red')
plt.xlim(40, 90)
plt.ylim(0, hist_data.y.max()*1.2)
plt.xlabel(r"Vertex $z$ (cm)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

print("K+K- invariant mass")
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('NominalCut/observable_phi_mass_NominalCut')
hist_sim  = file_sim .get('NominalCut/observable_phi_mass_NominalCut')
hist_sim.scale(hist_data.y.max()/hist_sim.y.max())
# hist_sim.areaNorm(hist_data)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim .plotPoints(label="Sim",  marker='o', markersize=3, linestyle='None', color='blue')
# plt.plot(x_fit, y_fit, label="Rel BW + Linear Fit", color='red')
plt.plot([1.005, 1.005], [0, 3000], '-', color='red')
plt.plot([1.04, 1.04],   [0, 3000], '-', color='red')
plt.xlim(0.98, 1.1)
plt.ylim(0, hist_data.y.max()*1.2)
plt.xlabel(r"$M_{K^+K^-} (GeV/c^2)$", size=14)
plt.ylabel("Counts / 2 MeV", size=14)
plt.xticks(size=12)
plt.yticks(size=12)
plt.legend()
file_pdf.savefig()
plt.close()

print("################################################################# OBSERVABLE PLOTS #################################################################")

# # Beam energy
# print("Beam energy")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/beam_energy_PlotCut')
# hist_sim = file_sim.get('PlotCut/beam_energy_PlotCut')
# hist_sim.scale(0.0015)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=0.1, alpha=0.5, label='Data', color='green')
# plt.xlim(5, 11)
# plt.xlabel(r"Beam $E$ (GeV)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Beam Timing
# print("Beam timing")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/beam_DeltaT_PlotCut')
# hist_sim = file_sim.get('PlotCut/beam_DeltaT_PlotCut')
# hist_sim.scale(0.004)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(-20, 20)
# plt.xlabel(r"Beam $\Delta t$ (ns)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Beam Accidental Contamination
# print("Beam Accidental Contamination")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/beam_DeltaT_PlotCut')
# on_time_low = np.where(hist_data.x > -2)[0][0]
# on_time_high = np.where(hist_data.x > 2)[0][0]
# plt.plot(hist_data.x[on_time_low:on_time_high], hist_data.y[on_time_low:on_time_high], marker='o', markersize=3, linestyle='-', color='black', label="On-time")
# off_time_low = [np.where(hist_data.x > -18)[0][0], np.where(hist_data.x > -14)[0][0], np.where(hist_data.x > -10)[0][0], np.where(hist_data.x > -6)[0][0], np.where(hist_data.x > 2)[0][0], np.where(hist_data.x > 6)[0][0], np.where(hist_data.x > 10)[0][0], np.where(hist_data.x > 14)[0][0]]
# off_time_high = [np.where(hist_data.x > -14)[0][0], np.where(hist_data.x > -10)[0][0], np.where(hist_data.x > -6)[0][0], np.where(hist_data.x > -2)[0][0], np.where(hist_data.x > 6)[0][0], np.where(hist_data.x > 10)[0][0], np.where(hist_data.x > 14)[0][0], np.where(hist_data.x > 18)[0][0]]
# off_time_total = np.zeros_like(hist_data.y[on_time_low:on_time_high])
# for low, high in zip(off_time_low, off_time_high):
#     off_time_total += hist_data.y[low:high]
# plt.plot(hist_data.x[on_time_low:on_time_high], -off_time_total, marker='o', markersize=3, linestyle='-', color='red', label="Off-time")
# plt.xlim(-2, 2)
# plt.xlabel(r"Beam $\Delta t$ (ns)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# ###################################################################### K PLUS #####################################################################################

# # K+ timing
# print("K+ timing")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/kp_DeltaT_meas_PlotCut')
# hist_sim  = file_sim. get('PlotCut/kp_DeltaT_meas_PlotCut')
# # hist_sim.scale(0.002)
# hist_sim.areaNorm(hist_data)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim. plotPoints(label="Sim",  marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlabel(r"$K^+ \Delta t$ (ns)")
# plt.ylabel("Counts")
# plt.xlim(-2, 2)
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # K+ timing vs p
# print("K+ timing vs p")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_DeltaT_momentum_meas_PlotCut')
# hist_sim  = file_sim. get('PlotCut/kp_DeltaT_momentum_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(-2, 2)
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$\Delta t_{K^+}$ (ns)")
# file_pdf.savefig()
# plt.close()

# # K+ CDC dE/dx
# print("K+ CDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_dEdx_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_dEdx_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.plot(np.linspace(0, 10, 100), np.exp(-7.0*np.linspace(0, 10, 100)+3.0)+6.2, color='red')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.plot(np.linspace(0, 10, 100), np.exp(-7.0*np.linspace(0, 10, 100)+3.0)+6.2, color='red')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^+}^{\mathrm{CDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K+ FDC dE/dx
# print("K+ FDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_dEdx_fdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_dEdx_fdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^+}^{\mathrm{FDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K+ TOF dE/dx
# print("K+ TOF dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_dEdx_tof_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_dEdx_tof_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^+}^{\mathrm{TOF}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K+ ST dE/dx
# print("K+ ST dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_dEdx_st_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_dEdx_st_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^+}^{\mathrm{ST}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K+ kinematics
# print("K+ kinematics")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_kinematics_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_kinematics_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^+}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K+ kinematics in FDC only
# print("K+ kinematics in FDC only")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_kinematics_fdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_kinematics_fdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^+}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K+ kinematics in FDC and CDC
# print("K+ kinematics in FDC and CDC")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_kinematics_fdc_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_kinematics_fdc_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^+}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K+ kinematics in CDC only
# print("K+ kinematics in CDC only")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/kp_kinematics_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/kp_kinematics_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^+}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^+}$ (deg)")
# file_pdf.savefig()
# plt.close()

# ###################################################################### K MINUS #####################################################################################

# # K- timing
# print("K- timing")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/km_DeltaT_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_DeltaT_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(-2, 2)
# plt.xlabel(r"$K^-$ $\Delta t$ (ns)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # K- timing vs p
# print("K- timing vs p")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_DeltaT_momentum_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_DeltaT_momentum_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(-2, 2)
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$\Delta t_{K^-}$ (ns)")
# file_pdf.savefig()
# plt.close()

# # K- CDC dE/dx
# print("K- CDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_dEdx_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_dEdx_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.plot(np.linspace(0, 10, 100), np.exp(-7.0*np.linspace(0, 10, 100)+3.0)+6.2, color='red')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.plot(np.linspace(0, 10, 100), np.exp(-7.0*np.linspace(0, 10, 100)+3.0)+6.2, color='red')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^-}^{\mathrm{CDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K- FDC dE/dx
# print("K- FDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_dEdx_fdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_dEdx_fdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^-}^{\mathrm{FDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K- TOF dE/dx
# print("K- TOF dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_dEdx_tof_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_dEdx_tof_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^-}^{\mathrm{TOF}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K- ST dE/dx
# print("K- ST dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_dEdx_st_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_dEdx_st_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_{K^-}^{\mathrm{ST}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # K- kinematics
# print("K- kinematics")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_kinematics_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_kinematics_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^-}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K- kinematics in FDC only
# print("K- kinematics in FDC only")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_kinematics_fdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_kinematics_fdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^-}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K- kinematics in FDC and CDC
# print("K- kinematics in FDC and CDC")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_kinematics_fdc_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_kinematics_fdc_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^-}$ (deg)")
# file_pdf.savefig()
# plt.close()

# # K- kinematics in CDC only
# print("K- kinematics in CDC only")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/km_kinematics_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/km_kinematics_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.plot([0,10], [2,2], '-', color='red')
#     plt.plot([0.4,0.4], [0,20], '-', color='red')
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^-}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{K^-}$ (deg)")
# file_pdf.savefig()
# plt.close()

# ###################################################################### DEUTERON #####################################################################################

# # d timing
# print("d timing")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/d_DeltaT_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_DeltaT_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(-2, 2)
# plt.xlabel(r"$d \ \Delta t$ (ns)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # d timing vs p
# print("d timing vs p")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_DeltaT_momentum_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_DeltaT_momentum_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(-2, 2)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$d \ \Delta t$ (ns)")
# file_pdf.savefig()
# plt.close()

# # d CDC dE/dx
# print("d CDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_dEdx_cdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_dEdx_cdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_d^{\mathrm{CDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # d FDC dE/dx
# print("d FDC dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_dEdx_fdc_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_dEdx_fdc_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_d^{\mathrm{FDC}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # d TOF dE/dx
# print("d TOF dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_dEdx_tof_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_dEdx_tof_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_d^{\mathrm{TOF}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # d ST dE/dx
# print("d ST dE/dx")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_dEdx_st_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_dEdx_st_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$(dE/dx)_d^{\mathrm{ST}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # d dE/dx in CDC vs ST
# print("d dE/dx in CDC vs ST")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_dEdx_cdc_st_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_dEdx_cdc_st_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 40)
#     plt.ylim(0, 40)
#     plt.xlabel(r"$(dE/dx)_{d}^{\mathrm{CDC}}$ (keV/cm)")
#     plt.ylabel(r"$(dE/dx)_{d}^{\mathrm{ST}}$ (keV/cm)")
# file_pdf.savefig()
# plt.close()

# # d kinematics
# print("d kinematics")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/d_kinematics_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/d_kinematics_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 90)
#     plt.yticks(np.arange(0, 91, 10))
#     plt.xlabel(r"$p_{d}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{d}$ (deg)")
# file_pdf.savefig()
# plt.close()

# ###################################################################### PHI MESON #####################################################################################

# # phi mass
# print("Phi mass")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/phi_mass_kin_PlotCut')
# hist_sim  = file_sim. get('PlotCut/phi_mass_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim .plotPoints(label="Sim",  marker='o', markersize=3, linestyle='-',    color='blue')
# plt.xlim(0.9, 1.1)
# plt.xlabel(r"$m_{K^+ K^-}$ (GeV/$c$)" )
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # phi mass vs KinFit chi2
# print("Phi mass vs KinFit chi2")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/phi_mass_chisq_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/phi_mass_chisq_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0.9, 1.1)
#     plt.ylim(0, 10)
#     plt.xlabel(r"$m_{K^+ K^-}$ (GeV/$c$)" )
#     plt.ylabel(r"$\chi^2$/NDF")
# file_pdf.savefig()
# plt.close()

# # phi mass vs minus t
# print("Phi mass vs -t")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/phi_mass_minust_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/phi_mass_minust_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0.9, 1.1)
#     plt.ylim(0, 2)
#     plt.xlabel(r"$m_{K^+ K^-}$ (GeV/$c$)" )
#     plt.ylabel(r"$-t (\mathrm{GeV}^2/c^4)$")
# file_pdf.savefig()
# plt.close()

# # phi mass vs missing p minus
# print("Phi mass vs Missing p minus")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/phi_mass_miss_pminus_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/phi_mass_miss_pminus_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0.9, 1.1)
#     plt.ylim(-0.1, 0.1)
#     plt.xlabel(r"$m_{K^+ K^-}$ (GeV/$c$)" )
#     plt.ylabel(r"$p_{miss}^-$ (GeV/c)")
# file_pdf.savefig()
# plt.close()

# # phi meson kinematics
# print("Phi kinematics")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/phi_kinematics_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/phi_kinematics_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 90)
#     plt.yticks(np.arange(0, 91, 10))
#     plt.xlabel(r"$p_{\phi}$ (GeV/c)")
#     plt.ylabel(r"$\theta_{\phi}$ (deg)")
# file_pdf.savefig()
# plt.close()

# ###################################################################### MISSING P4 #####################################################################################

# # Missing energy
# print("Missing energy")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/miss_energy_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_energy_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
# plt.xlim(-2, 2)
# plt.xlabel(r"$E_{miss}$ (GeV)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Missing mass squared
# print("Missing mass squared")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/miss_masssquared_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_masssquared_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='-', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='blue')
# # plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
# plt.xlim(-0.1, 0.1)
# plt.xlabel(r"$m^2_{miss} (GeV^2/c^4)$")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Missing p minus
# print("Missing p minus")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/miss_pminus_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_pminus_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='-', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='blue')
# # plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
# plt.xlim(-0.1, 0.1)
# plt.xlabel(r"$p_{miss}^-$ (GeV/c)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Missing momentum
# print("Missing momentum")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/miss_momentum_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_momentum_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='-', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='blue')
# # plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
# plt.xlim(0, 2.0)
# plt.xlabel(r"$p_{miss} (GeV/c)$")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Missing momentum vs missing energy
# print("Missing momentum vs Missing energy")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/miss_momentum_energy_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_momentum_energy_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.ylim(-2, 2)
#     plt.xlim(0, 2)
#     plt.ylabel(r"$E_{miss}$ (GeV)")
#     plt.xlabel(r"$p_{miss} (GeV/c)$")
# file_pdf.savefig()
# plt.close()

# # Missing energy vs phi mass
# print("Missing energy vs Phi mass")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/miss_energy_phi_mass_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/miss_energy_phi_mass_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.05, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='bottom')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.05, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='bottom')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(-2, 2)
#     plt.ylim(0.9, 1.1)
#     plt.xlabel(r"$E_{miss}$ (GeV)")
#     plt.ylabel(r"$m_{K^+ K^-}$ (GeV/$c$)" )
# file_pdf.savefig()
# plt.close()

# ###################################################################### BACKGROUND #####################################################################################

# # Number of unused tracks
# print("Number of unused tracks")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/num_unused_tracks_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/num_unused_tracks_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
# plt.xlim(0, 5)
# plt.xlabel("Number of Unused Tracks")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Number of unused showers
# print("Number of unused showers")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/num_unused_showers_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/num_unused_showers_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
# plt.xlim(0, 10)
# plt.xlabel("Number of Unused Showers")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Beam accidental weight
# print("Beam accidental weight")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/beam_accid_weight_PlotCut')
# hist_sim = file_sim.get('PlotCut/beam_accid_weight_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=0.1, alpha=0.5, label='Data', color='green')
# plt.xlim(-0.5, 1.5)
# plt.xlabel("Beam Accidental Weight")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Combo accidental weight
# print("Combo accidental weight")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/combo_accid_weight_PlotCut')
# hist_sim = file_sim.get('PlotCut/combo_accid_weight_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
# plt.xlim(-0.5, 1.5)
# plt.xlabel("Combo Accidental Weight")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Vertex z
# print("Vertex z")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/vertex_z_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/vertex_z_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(40, 90)
# plt.xlabel(r"Vertex $z$ (cm)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Vertex x-y
# print("Vertex x-y")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/vertex_x_y_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/vertex_x_y_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(-1, 1)
#     plt.ylim(-1, 1)
#     plt.xlabel(r"Vertex $x$ (cm)")
#     plt.ylabel(r"Vertex $y$ (cm)")
# file_pdf.savefig()
# plt.close()

# # rho mass
# print("Rho mass")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/rho_mass_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/rho_mass_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(0.0, 1.0)
# plt.xlabel(r"$m_{\pi^+ \pi^-}$ (GeV/$c$)" )
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # rho missing p minus
# print("Rho missing p minus")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/rho_miss_pminus_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/rho_miss_pminus_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(-0.1, 0.1)
# plt.xlabel(r"$p_{miss, \pi^+ \pi^-}^-$ (GeV/c)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # KinFit FOM
# print("KinFit FOM")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/kinfit_fom_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/kinfit_fom_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(-5, 0)
# plt.yscale('log')
# plt.xlabel("KinFit FOM")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # KinFit Chi2
# print("KinFit Chi2")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/chisq_per_ndf_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/chisq_per_ndf_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(0, 10)
# plt.xlabel(r"$\chi^2$/NDF")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# ###################################################################### KINEMATICS #####################################################################################

# # minus t
# print("-t")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/minust_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/minust_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(0, 2)
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # minus t vs deuteron momentum
# print("-t vs Deuteron Momentum")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/minust_kin_d_momentum_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/minust_kin_d_momentum_meas_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(0, 2)
#     plt.xlabel(r"$p_d$ (GeV/c)")
#     plt.ylabel(r"$-t (GeV^2/c^4)$")
# file_pdf.savefig()
# plt.close()

# # minus t vs beam energy
# print("-t vs Beam Energy")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/beam_energy_minust_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/beam_energy_minust_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.xlim(5, 11)
#     plt.ylim(0, 2)
#     plt.xlabel(r"Beam $E$ (GeV)")
#     plt.ylabel(r"$-t (GeV^2/c^4)$")
# file_pdf.savefig()
# plt.close()

# # com scattering angle
# print("COM Scattering Angle")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/scatter_theta_com_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/scatter_theta_com_kin_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(0, 180)
# plt.xlabel(r"$\theta_{COM}$ (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # minus t vs com scattering angle
# print("-t vs COM Scattering Angle")
# fig = plt.figure(figsize=(16, 6), dpi=300)
# gs = fig.add_gridspec(1, 2)
# axs = gs.subplots()
# hist_data = file_data.get('PlotCut/minust_scatter_theta_com_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/minust_scatter_theta_com_kin_PlotCut')
# hist_sim.scale(0.002)
# plt.axes(axs[0])
# hist_data.plotHeatmap()
# plt.text(0.95, 0.95, "data", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "simulation", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(2):
#     plt.axes(axs[i])
#     plt.ylim(0, 40)
#     plt.xlim(0, 2)
#     plt.ylabel(r"$\theta_{COM}$ (deg)")
#     plt.xlabel(r"$-t (GeV^2/c^4)$")
# file_pdf.savefig()
# plt.close()

# # Coplanarity
# print("Coplanarity")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/coplanarity_meas_PlotCut')
# hist_sim = file_sim.get('PlotCut/coplanarity_meas_PlotCut')
# hist_sim.scale(0.002)
# hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='blue')
# plt.xlim(170, 190)
# plt.xlabel(r"Coplanarity (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# ###################################################################### DECAY ANGLES #####################################################################################

# # cos(vartheta)
# print("cos(vartheta)")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/decay_costheta_helicity_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/decay_costheta_helicity_kin_PlotCut')
# hist_sim.scale(0.018)
# hist_data.plotPoints(label="Data", marker='o', markersize=7, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=7, linestyle='-', color='blue')
# plt.xlim(-1, 1)
# plt.ylim(0, 2000)
# plt.xlabel(r"$\cos\vartheta$")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # Phi
# print("Phi")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/polarization_phi_com_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/polarization_phi_com_kin_PlotCut')
# hist_sim.scale(0.018)
# hist_data.plotPoints(label="Data", marker='o', markersize=7, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=7, linestyle='-', color='blue')
# plt.xlim(-180, 180)
# plt.ylim(0, 2000)
# plt.xlabel(r"$\Phi$ (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # varphi
# print("varphi")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/decay_phi_helicity_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/decay_phi_helicity_kin_PlotCut')
# hist_sim.scale(0.018)
# hist_data.plotPoints(label="Data", marker='o', markersize=7, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=7, linestyle='-', color='blue')
# plt.xlim(-180, 180)
# plt.ylim(0, 2000)
# plt.xlabel(r"$\varphi$ (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # psi
# print("psi")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_data = file_data.get('PlotCut/psi_helicity_kin_PlotCut')
# hist_sim = file_sim.get('PlotCut/psi_helicity_kin_PlotCut')
# hist_sim.scale(0.018)
# hist_data.plotPoints(label="Data", marker='o', markersize=7, linestyle='None', color='black')
# hist_sim.plotPoints(label="Sim", marker='o', markersize=7, linestyle='-', color='blue')
# plt.xlim(-180, 180)
# plt.ylim(0, 1000)
# plt.xlabel(r"$\psi$ (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# ##################################################################### TRUTH VARIABLES #####################################################################################

# mass_photon     = 0
# mass_deuteron   = 1.8756129
# mass_phi        = 1.019445
# Deg_To_Rad      = np.pi/180

# energy_list = np.arange(6.0, 11.0, 1.0)
# t_list = np.arange(-2.0, 0.0, 0.1)

# phi_p = np.zeros((len(energy_list), len(t_list)))
# phi_theta = np.zeros((len(energy_list), len(t_list)))
# deuteron_p = np.zeros((len(energy_list), len(t_list)))
# deuteron_theta = np.zeros((len(energy_list), len(t_list)))

# for i, energy in enumerate(energy_list):
#     for j, t in enumerate(t_list):
#         P4_gamma = root.TLorentzVector(0, 0, energy, energy)
#         P4_target = root.TLorentzVector(0, 0, 0, mass_deuteron)
#         s = (P4_gamma + P4_target).Mag2()

#         pi_cm = np.sqrt(((s-mass_phi**2-mass_photon**2)**2-4*mass_phi**2*mass_photon**2)/(4*s))
#         pf_cm = np.sqrt(((s-mass_phi**2-mass_deuteron**2)**2-4*mass_phi**2*mass_deuteron**2)/(4*s))
#         edi_cm = np.sqrt(pi_cm**2+mass_deuteron**2)
#         edf_cm = np.sqrt(pf_cm**2+mass_deuteron**2)
#         theta = np.arccos((t - 2*mass_deuteron**2 + 2*edi_cm*edf_cm)/(2*pi_cm*pf_cm))

#         P4_phi = root.TLorentzVector()
#         P4_deuteron = root.TLorentzVector()
#         P4_phi.SetXYZM(0, pf_cm*np.sin(theta), pf_cm*np.cos(theta), mass_phi)
#         P4_deuteron.SetXYZM(0, -pf_cm*np.sin(theta), -pf_cm*np.cos(theta), mass_deuteron)

#         boost_cm = (P4_gamma+P4_target).BoostVector()

#         P4_phi.Boost(boost_cm)
#         P4_deuteron.Boost(boost_cm)

#         phi_p[i][j] = P4_phi.P()
#         phi_theta[i][j] = P4_phi.Theta()/Deg_To_Rad
#         deuteron_p[i][j] = P4_deuteron.P()
#         deuteron_theta[i][j] = P4_deuteron.Theta()/Deg_To_Rad

# # Beam energy truth
# print("Beam energy truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim    = file_sim.get('PlotCut/beam_energy_truth_PlotCut')
# hist_tagged = file_tagged.get('NoCut/beam_energy_truth_NoCut')
# hist_gen    = file_gen.get('NoCut/beam_energy_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(5, 11)
# plt.xlabel(r"Beam $E$ Truth (GeV)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # phi mass truth
# print("Phi mass truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/phi_mass_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/phi_mass_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/phi_mass_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(0.9, 1.1)
# plt.xlabel(r"$m_{K^+ K^-}$ Truth $(GeV/c^2)$")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # K+ kinematics truth
# print("K+ kinematics truth")
# fig = plt.figure(figsize=(24, 6), dpi=300)
# gs = fig.add_gridspec(1, 3)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/kp_kinematics_truth_PlotCut')
# hist_tagged = file_tagged.get('NoCut/kp_kinematics_truth_NoCut')
# hist_gen = file_gen.get('NoCut/kp_kinematics_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(3):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^+}$ truth (GeV/c)")
#     plt.ylabel(r"$\theta_{K^+}$ truth (deg)")
# file_pdf.savefig()
# plt.close()

# # K- kinematics truth
# print("K- kinematics truth")
# fig = plt.figure(figsize=(24, 6), dpi=300)
# gs = fig.add_gridspec(1, 3)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/km_kinematics_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/km_kinematics_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/km_kinematics_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(3):
#     plt.axes(axs[i])
#     plt.xlim(0, 10)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{K^-}$ truth (GeV/c)")
#     plt.ylabel(r"$\theta_{K^-}$ truth (deg)")
# file_pdf.savefig()
# plt.close()

# # d kinematics truth
# print("d kinematics truth")
# fig = plt.figure(figsize=(32, 6), dpi=300)
# gs = fig.add_gridspec(1, 4)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/d_kinematics_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/d_kinematics_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/d_kinematics_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[3])
# for i in range(np.size(energy_list)):
#     plt.plot(deuteron_p[i], deuteron_theta[i], linestyle='-', label=f"E={energy_list[i]}")
#     plt.text(deuteron_p[i][-1], deuteron_theta[i][-1], f"{energy_list[i]:.1f}", fontsize=6, ha='right', va='top')
# for i in range(np.size(t_list)):
#     plt.plot(deuteron_p[:,i], deuteron_theta[:,i], linestyle='-')
#     plt.text(deuteron_p[-1][i], deuteron_theta[-1][i], f"{-t_list[i]:.1f}", fontsize=6, ha='left', va='bottom')
# plt.text(0.95, 0.95, "kinematics", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(4):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(50, 90)
#     plt.yticks(np.arange(50, 91, 10))
#     plt.xlabel(r"$p_{d}$ truth (GeV/c)")
#     plt.ylabel(r"$\theta_{d}$ truth(deg)")
# file_pdf.savefig()
# plt.close()

# # phi meson kinematics truth
# print("Phi kinematics truth")
# fig = plt.figure(figsize=(32, 6), dpi=300)
# gs = fig.add_gridspec(1, 4)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/phi_kinematics_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/phi_kinematics_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/phi_kinematics_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[3])
# for i in range(np.size(energy_list)):
#     plt.plot(phi_p[i], phi_theta[i], linestyle='-', label=f"E={energy_list[i]}")
#     plt.text(phi_p[i][-1], phi_theta[i][-1], f"{energy_list[i]:.1f}", fontsize=6, ha='right', va='top')
# for i in range(np.size(t_list)):
#     plt.plot(phi_p[:,i], phi_theta[:,i], linestyle='-')
#     plt.text(phi_p[-1][i], phi_theta[-1][i], f"{-t_list[i]:.1f}", fontsize=6, ha='left', va='bottom')
# plt.text(0.95, 0.95, "kinematics", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(4):
#     plt.axes(axs[i])
#     plt.xlim(5, 11)
#     plt.ylim(0, 20)
#     plt.yticks(np.arange(0, 21, 2))
#     plt.xlabel(r"$p_{\phi}$ truth (GeV/c)")
#     plt.ylabel(r"$\theta_{\phi}$ truth (deg)")
# file_pdf.savefig()
# plt.close()

# # phi meson and d theta truth
# print("Phi and d theta truth")
# fig = plt.figure(figsize=(32, 6), dpi=300)
# gs = fig.add_gridspec(1, 4)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/phi_d_theta_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/phi_d_theta_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/phi_d_theta_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[3])
# for i in range(np.size(energy_list)):
#     plt.plot(deuteron_theta[i], phi_theta[i], linestyle='-', label=f"E={energy_list[i]}")
#     plt.text(deuteron_theta[i][-1], phi_theta[i][-1], f"{energy_list[i]:.1f}", fontsize=6, ha='right', va='top')
# for i in range(np.size(t_list)):
#     plt.plot(deuteron_theta[:,i], phi_theta[:,i], linestyle='-')
#     plt.text(deuteron_theta[-1][i], phi_theta[-1][i], f"{-t_list[i]:.1f}", fontsize=6, ha='left', va='bottom')
# plt.text(0.95, 0.95, "kinematics", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(4):
#     plt.axes(axs[i])
#     plt.xlim(50, 90)
#     plt.ylim(0, 20)
#     plt.xlabel(r"$\theta_{d}$ truth (deg)")
#     plt.ylabel(r"$\theta_{\phi}$ truth (deg)")
# file_pdf.savefig()
# plt.close()

# # phi meson and d momentum truth
# print("Phi and d momentum truth")
# fig = plt.figure(figsize=(32, 6), dpi=300)
# gs = fig.add_gridspec(1, 4)
# axs = gs.subplots()
# hist_sim = file_sim.get('PlotCut/phi_d_momentum_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/phi_d_momentum_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/phi_d_momentum_truth_NoCut')
# hist_sim.scale(110)
# plt.axes(axs[0])
# hist_sim.plotHeatmap()
# plt.text(0.95, 0.95, "detected", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[1])
# hist_tagged.plotHeatmap()
# plt.text(0.95, 0.95, "tagged", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[2])
# hist_gen.plotHeatmap()
# plt.text(0.95, 0.95, "gen", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# plt.axes(axs[3])
# for i in range(np.size(energy_list)):
#     plt.plot(deuteron_p[i], phi_p[i], linestyle='-', label=f"E={energy_list[i]}")
#     plt.text(deuteron_p[i][-1], phi_p[i][-1], f"{energy_list[i]:.1f}", fontsize=6, ha='right', va='top')
# for i in range(np.size(t_list)):
#     plt.plot(deuteron_p[:,i], phi_p[:,i], linestyle='-')
#     plt.text(deuteron_p[-1][i], phi_p[-1][i], f"{-t_list[i]:.1f}", fontsize=6, ha='left', va='bottom')
# plt.text(0.95, 0.95, "kinematics", transform=plt.gca().transAxes, fontsize=16, fontstyle='italic', ha='right', va='top')
# for i in range(4):
#     plt.axes(axs[i])
#     plt.xlim(0, 2)
#     plt.ylim(5, 11)
#     plt.xlabel(r"$p_{d}$ truth (GeV/c)")
#     plt.ylabel(r"$p_{\phi}$ truth (GeV/c)")
# file_pdf.savefig()
# plt.close()

# # minus t truth
# print("-t truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/minust_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/minust_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/minust_truth_NoCut')
# hist_sim.scale(4)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(0, 2)
# plt.xlabel(r"$-t$ truth $(\rm{GeV}^2/c^4)$")
# plt.ylabel("Counts")
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # com scattering angle truth
# print("COM scattering angle truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/scatter_theta_com_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/scatter_theta_com_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/scatter_theta_com_truth_NoCut')
# hist_sim.scale(4)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(0, 180)
# plt.xlabel(r"$\theta_{COM}$ Truth (deg)")
# plt.ylabel("Counts")
# plt.yscale('log')
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # decay cos(vartheta) truth
# print("Decay cos(vartheta) truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/decay_costheta_helicity_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/decay_costheta_helicity_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/decay_costheta_helicity_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(-1, 1)
# plt.xlabel(r"Decay $\cos\vartheta$ Truth")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # decay phi truth
# print("Decay phi truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/decay_phi_helicity_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/decay_phi_helicity_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/decay_phi_helicity_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(-180, 180)
# plt.xlabel(r"Decay $\varphi$ Truth (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # polarization phi truth
# print("Polarization phi truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/polarization_phi_com_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/polarization_phi_com_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/polarization_phi_com_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(-180, 180)
# plt.xlabel(r"Polarization $\Phi$ Truth (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# # decay psi truth
# print("Decay psi truth")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/psi_helicity_truth_PlotCut')
# hist_gen = file_gen.get('NoCut/psi_helicity_truth_NoCut')
# hist_tagged = file_tagged.get('NoCut/psi_helicity_truth_NoCut')
# hist_sim.scale(110)
# hist_sim.plotPoints(label="detected", marker='o', markersize=3, linestyle='None', color='blue')
# hist_tagged.plotPoints(label="tagged", marker='o', markersize=3, linestyle='None', color='orange')
# hist_gen.plotPoints(label="gen", marker='o', markersize=3, linestyle='None', color='green')
# plt.xlim(-180, 180)
# plt.xlabel(r"Decay $\psi$ Truth (deg)")
# plt.ylabel("Counts")
# plt.legend()
# file_pdf.savefig()
# plt.close()

# ###################################################################### BIN MIGRATIONS #####################################################################################

# # beam energy bin migration
# print("beam energy bin migration")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/beam_energy_kin_truth_PlotCut')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"Beam $E$ detected")
# plt.xlabel(r"Beam $E$ thrown")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # minus t bin migration
# print("minus t bin migration")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/minust_kin_truth_PlotCut')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$-t$ reconstructed $(\mathrm{GeV}^2/c^4)$")
# plt.xlabel(r"$-t$ thrown $(\mathrm{GeV}^2/c^4)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

###################################################################### RESOLUTIONS #####################################################################################

def gaussian(x, a, b, c):
    return a*np.exp(-0.5*((x-b)/c)**2)

# # minus t resolution
# print("-t resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/minust_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta t (GeV^2/c^4)$")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.ylim(-0.1, 0.1)
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # costheta resolution
# print("costheta resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/decay_costheta_helicity_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta \cos\vartheta$")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.ylim(-0.5, 0.5)
# file_pdf.savefig()
# plt.close()

# # decayphi resolution
# print("decayphi resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/decay_phi_helicity_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta \varphi$ (deg)")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # polphi resolution
# print("polphi resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/polarization_phi_com_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta \Phi$ (deg)")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # psi resolution
# print("psi resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/psi_helicity_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta \psi$ (deg)")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # vertex z resolution
# print("Vertex z resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('VertexCut/resolution_vertex_z_VertexCut')
# for i in range(20, 80):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
#     print(f"Bin {i}: Center = {popt[1]:.2f} ± {popt[2]:.2f}")
#     print(f"Bin {i}: Width = {popt[2]:.2f}")
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta Z$ (cm)")
# plt.xlabel(r"$Z$ (cm)")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# # vertex z resolution w.r.t. t
# print("Vertex z resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('VertexCut/resolution_t_vertex_z_VertexCut')
# for i in range(4, 40):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
#     print(f"Bin {i}: Center = {popt[1]:.2f} ± {popt[2]:.2f}")
#     print(f"Bin {i}: Width = {popt[2]:.2f}")
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta Z$ (cm)")
# plt.xlabel(r"$-t (GeV^2)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

# phi mass resolution
# print("Phi mass resolution")
# fig = plt.figure(figsize=(8, 6), dpi=300)
# hist_sim = file_sim.get('PlotCut/phi_mass_diff_PlotCut')
# for i in range(4, np.size(hist_sim.xedge)-1):
#     this_ycenter = (hist_sim.yedge[:-1] + hist_sim.yedge[1:]) / 2
#     this_average = np.average(this_ycenter, weights=hist_sim.z[:, i])
#     this_std = np.sqrt(np.average((this_ycenter - this_average)**2, weights=hist_sim.z[:, i]))
#     popt, pcov = curve_fit(gaussian, this_ycenter, hist_sim.z[:, i], p0=[hist_sim.z[:, i].max(),this_average,this_std])
#     plt.errorbar(hist_sim.xedge[i]+(hist_sim.xedge[i+1]-hist_sim.xedge[i])/2, popt[1], yerr=popt[2], fmt='.', color='red')
# hist_sim.z[hist_sim.z<=0] = np.nan
# hist_sim.z = np.log10(hist_sim.z)
# hist_sim.plotHeatmap()
# plt.ylabel(r"$\Delta$ Phi Mass (GeV/c^2)")
# plt.xlabel(r"$-t (GeV^2/c^4)$")
# plt.colorbar(label=r"$\log_{10}$(Counts)")
# file_pdf.savefig()
# plt.close()

###################################################################### END #####################################################################################

file_pdf.close()