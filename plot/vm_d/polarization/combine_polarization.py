import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import ROOT as root
from scipy.integrate import quad
import scipy.stats

def lumi(energy_min, energy_max, length):
    lumi_table = np.loadtxt('/work/halld2/home/boyu/src_analysis/flux/output/2H/lumi_summed_2H.txt')
    length_total = 29.5

    integrated_lumi = np.zeros(energy_min.shape, dtype=float)
    for i in range(len(energy_min)):
        for j in range(len(lumi_table)):
            if (lumi_table[j,3] > energy_min[i]) and (lumi_table[j,3] < energy_max[i]):
                integrated_lumi[i] += lumi_table[j][5]

    return integrated_lumi/length_total * length

###################################################################### DEFINITIONS #####################################################################################

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

    def plotHeatmap(self,kill_zeros=True,kill_nonpositive=True,log_scale=False,**kwargs):
        z = self.z
        if (kill_zeros):
            z[z==0] = np.nan
        if (kill_nonpositive):
            z[z<=0] = np.nan
        if (log_scale):
            return plt.pcolormesh(self.xedge,self.yedge,np.log10(self.z),**kwargs)
        else:
            return plt.pcolormesh(self.xedge,self.yedge,self.z,**kwargs)

    def projectionX(self,**kwargs):
        return Hist1D(self.TH2.ProjectionX(),**kwargs)

    def projectionY(self,**kwargs):
        return Hist1D(self.TH2.ProjectionY(),**kwargs)

###################################################################### READ THE FILES #####################################################################################

file_d2 = File("outFiles/2021-11_TPol_D2.root")
file_he = File("outFiles/2021-11_TPol_He4.root")
file_c = File("outFiles/2021-11_TPol_C12.root")
file_all = File("outFiles/2021-11_TPol_all.root")

d2_0    = file_d2.get("hPol0")
d2_45   = file_d2.get("hPol45")
d2_90   = file_d2.get("hPol90")
d2_135  = file_d2.get("hPol135")

d2_0.y[0] = np.nan
d2_0.yerr[0] = np.nan
d2_45.y[0] = np.nan
d2_45.yerr[0] = np.nan
d2_90.y[0] = np.nan
d2_90.yerr[0] = np.nan
d2_135.y[0] = np.nan
d2_135.yerr[0] = np.nan

he_0    = file_he.get("hPol0")
he_45   = file_he.get("hPol45")
he_90   = file_he.get("hPol90")
he_135  = file_he.get("hPol135")

he_0.y[0] = np.nan
he_0.yerr[0] = np.nan
he_45.y[0] = np.nan
he_45.yerr[0] = np.nan
he_90.y[0] = np.nan
he_90.yerr[0] = np.nan
he_135.y[0] = np.nan
he_135.yerr[0] = np.nan

c_0     = file_c.get("hPol0")
c_45    = file_c.get("hPol45")
c_90    = file_c.get("hPol90")
c_135   = file_c.get("hPol135")

c_0.y[0] = np.nan
c_0.yerr[0] = np.nan
c_45.y[0] = np.nan
c_45.yerr[0] = np.nan
c_90.y[0] = np.nan
c_90.yerr[0] = np.nan
c_135.y[0] = np.nan
c_135.yerr[0] = np.nan

all_0   = file_all.get("hPol0")
all_45  = file_all.get("hPol45")
all_90  = file_all.get("hPol90")
all_135 = file_all.get("hPol135")

all_0.y[0] = np.nan
all_0.yerr[0] = np.nan
all_45.y[0] = np.nan
all_45.yerr[0] = np.nan
all_90.y[0] = np.nan
all_90.yerr[0] = np.nan
all_135.y[0] = np.nan
all_135.yerr[0] = np.nan

d2_combined_x               = d2_0.x
d2_combined_polarization    = (d2_0.y/d2_0.yerr**2+d2_45.y/d2_45.yerr**2+d2_90.y/d2_90.yerr**2+d2_135.y/d2_135.yerr**2)/(1/d2_0.yerr**2+1/d2_45.yerr**2+1/d2_90.yerr**2+1/d2_135.yerr**2)
d2_combined_error           = np.sqrt(1/(1/d2_0.yerr**2+1/d2_45.yerr**2+1/d2_90.yerr**2+1/d2_135.yerr**2))
he_combined_x               = he_0.x
he_combined_polarization    = (he_0.y/he_0.yerr**2+he_45.y/he_45.yerr**2+he_90.y/he_90.yerr**2+he_135.y/he_135.yerr**2)/(1/he_0.yerr**2+1/he_45.yerr**2+1/he_90.yerr**2+1/he_135.yerr**2)
he_combined_error           = np.sqrt(1/(1/he_0.yerr**2+1/he_45.yerr**2+1/he_90.yerr**2+1/he_135.yerr**2))
c_combined_x                = c_0.x
c_combined_polarization     = (c_0.y/c_0.yerr**2+c_45.y/c_45.yerr**2+c_90.y/c_90.yerr**2+c_135.y/c_135.yerr**2)/(1/c_0.yerr**2+1/c_45.yerr**2+1/c_90.yerr**2+1/c_135.yerr**2)
c_combined_error            = np.sqrt(1/(1/c_0.yerr**2+1/c_45.yerr**2+1/c_90.yerr**2+1/c_135.yerr**2))
all_combined_x              = all_0.x
all_combined_polarization   = (all_0.y/all_0.yerr**2+all_45.y/all_45.yerr**2+all_90.y/all_90.yerr**2+all_135.y/all_135.yerr**2)/(1/all_0.yerr**2+1/all_45.yerr**2+1/all_90.yerr**2+1/all_135.yerr**2)
all_combined_error          = np.sqrt(1/(1/all_0.yerr**2+1/all_45.yerr**2+1/all_90.yerr**2+1/all_135.yerr**2))

fig = plt.figure(figsize=(16,30), dpi=300)
gs = fig.add_gridspec(5,2)
axs = gs.subplots().flatten()
fig.supxlabel("Photon Energy (GeV)")
fig.supylabel("Polarization")

plt.axes(axs[0])
plt.errorbar(d2_0.x+0.02, d2_0.y, yerr=d2_0.yerr, label="D2 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_0.x+0.04, he_0.y, yerr=he_0.yerr, label="He 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_0.x+0.06, c_0.y, yerr=c_0.yerr, label="C 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_0.x, all_0.y, yerr=all_0.yerr, label="All 0 deg", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.xlim(5.8,10.8)
plt.legend()

plt.axes(axs[1])
plt.errorbar(all_0.x, (d2_0.y-all_0.y)/np.sqrt(np.abs(d2_0.yerr**2-all_0.yerr**2)), label="D2 - All 0 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_0.x, (he_0.y-all_0.y)/np.sqrt(np.abs(he_0.yerr**2-all_0.yerr**2)), label="He - All 0 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_0.x, (c_0.y-all_0.y)/np.sqrt(np.abs(c_0.yerr**2-all_0.yerr**2)), label="C - All 0 deg", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[2])
plt.errorbar(d2_45.x+0.02, d2_45.y, yerr=d2_45.yerr, label="D2 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_45.x+0.04, he_45.y, yerr=he_45.yerr, label="He 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_45.x+0.06, c_45.y, yerr=c_45.yerr, label="C 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_45.x, all_45.y, yerr=all_45.yerr, label="All 45 deg", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[3])
plt.errorbar(all_45.x, (d2_45.y-all_45.y)/np.sqrt(np.abs(d2_45.yerr**2-all_45.yerr**2)), label="D2 - All 45 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_45.x, (he_45.y-all_45.y)/np.sqrt(np.abs(he_45.yerr**2-all_45.yerr**2)), label="He - All 45 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_45.x, (c_45.y-all_45.y)/np.sqrt(np.abs(c_45.yerr**2-all_45.yerr**2)), label="C - All 45 deg", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[4])
plt.errorbar(d2_90.x+0.02, d2_90.y, yerr=d2_90.yerr, label="D2 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_90.x+0.04, he_90.y, yerr=he_90.yerr, label="He 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_90.x+0.06, c_90.y, yerr=c_90.yerr, label="C 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_90.x, all_90.y, yerr=all_90.yerr, label="All 90 deg", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[5])
plt.errorbar(all_90.x, (d2_90.y-all_90.y)/np.sqrt(np.abs(d2_90.yerr**2-all_90.yerr**2)), label="D2 - All Target 90 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_90.x, (he_90.y-all_90.y)/np.sqrt(np.abs(he_90.yerr**2-all_90.yerr**2)), label="He - All Target 90 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_90.x, (c_90.y-all_90.y)/np.sqrt(np.abs(c_90.yerr**2-all_90.yerr**2)), label="C - All Target 90 deg", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[6])
plt.errorbar(d2_135.x+0.02, d2_135.y, yerr=d2_135.yerr, label="D2 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_135.x+0.04, he_135.y, yerr=he_135.yerr, label="He 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_135.x+0.06, c_135.y, yerr=c_135.yerr, label="C 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_135.x, all_135.y, yerr=all_135.yerr, label="All 135 deg", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[7])
plt.errorbar(all_135.x, (d2_135.y-all_135.y)/np.sqrt(np.abs(d2_135.yerr**2-all_135.yerr**2)), label="D2 - All Target 135 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_135.x, (he_135.y-all_135.y)/np.sqrt(np.abs(he_135.yerr**2-all_135.yerr**2)), label="He - All Target 135 deg", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_135.x, (c_135.y-all_135.y)/np.sqrt(np.abs(c_135.yerr**2-all_135.yerr**2)), label="C - All Target 135 deg", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[8])
plt.errorbar(d2_combined_x, d2_combined_polarization, yerr=d2_combined_error, label="D2", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_combined_x+0.04, he_combined_polarization, yerr=he_combined_error, label="He4", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_combined_x+0.06, c_combined_polarization, yerr=c_combined_error, label="C12", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_combined_x+0.02, all_combined_polarization, yerr=all_combined_error, label="D2+He4+C12", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[9])
plt.errorbar(d2_combined_x, (d2_combined_polarization-all_combined_polarization)/np.sqrt(np.abs(d2_combined_error**2-all_combined_error**2)), label="D2 - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.errorbar(he_combined_x, (he_combined_polarization-all_combined_polarization)/np.sqrt(np.abs(he_combined_error**2-all_combined_error**2)), label="He - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.errorbar(c_combined_x, (c_combined_polarization-all_combined_polarization)/np.sqrt(np.abs(c_combined_error**2-all_combined_error**2)), label="C - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.savefig("polarization_angle.png")
plt.close()

fig = plt.figure(figsize=(12,24), dpi=300)
gs = fig.add_gridspec(4,2)
axs = gs.subplots().flatten()
fig.supxlabel("Photon Energy (GeV)")
fig.supylabel("Polarization")

plt.axes(axs[0])
plt.errorbar(d2_0.x+0.02, d2_0.y, yerr=d2_0.yerr, label="D2 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(d2_45.x+0.04, d2_45.y, yerr=d2_45.yerr, label="D2 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(d2_90.x+0.06, d2_90.y, yerr=d2_90.yerr, label="D2 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(d2_135.x+0.08, d2_135.y, yerr=d2_135.yerr, label="D2 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(d2_combined_x, d2_combined_polarization, yerr=d2_combined_error, label="D2 Combined", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[1])
plt.errorbar(d2_0.x, (d2_0.y-d2_combined_polarization)/np.sqrt(np.abs(d2_0.yerr**2-d2_combined_error**2)), label="0 deg - All angle D2", marker="o",markersize=5, linestyle="None")
plt.errorbar(d2_45.x, (d2_45.y-d2_combined_polarization)/np.sqrt(np.abs(d2_45.yerr**2-d2_combined_error**2)), label="45 deg - All angle D2", marker="o",markersize=5, linestyle="None")
plt.errorbar(d2_90.x, (d2_90.y-d2_combined_polarization)/np.sqrt(np.abs(d2_90.yerr**2-d2_combined_error**2)), label="90 deg - All angle D2", marker="o",markersize=5, linestyle="None")
plt.errorbar(d2_135.x, (d2_135.y-d2_combined_polarization)/np.sqrt(np.abs(d2_135.yerr**2-d2_combined_error**2)), label="135 deg - All angle D2", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[2])
plt.errorbar(he_0.x+0.02, he_0.y, yerr=he_0.yerr, label="He 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_45.x+0.04, he_45.y, yerr=he_45.yerr, label="He 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_90.x+0.06, he_90.y, yerr=he_90.yerr, label="He 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_135.x+0.08, he_135.y, yerr=he_135.yerr, label="He 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(he_combined_x, he_combined_polarization, yerr=he_combined_error, label="He Combined", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[3])
plt.errorbar(he_0.x, (he_0.y-he_combined_polarization)/np.sqrt(np.abs(he_0.yerr**2-he_combined_error**2)), label="0 deg - All angle He", marker="o",markersize=5, linestyle="None")
plt.errorbar(he_45.x, (he_45.y-he_combined_polarization)/np.sqrt(np.abs(he_45.yerr**2-he_combined_error**2)), label="45 deg - All angle He", marker="o",markersize=5, linestyle="None")
plt.errorbar(he_90.x, (he_90.y-he_combined_polarization)/np.sqrt(np.abs(he_90.yerr**2-he_combined_error**2)), label="90 deg - All angle He", marker="o",markersize=5, linestyle="None")
plt.errorbar(he_135.x, (he_135.y-he_combined_polarization)/np.sqrt(np.abs(he_135.yerr**2-he_combined_error**2)), label="135 deg - All angle He", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[4])
plt.errorbar(c_0.x+0.02, c_0.y, yerr=c_0.yerr, label="C 0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_45.x+0.04, c_45.y, yerr=c_45.yerr, label="C 45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_90.x+0.06, c_90.y, yerr=c_90.yerr, label="C 90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_135.x+0.08, c_135.y, yerr=c_135.yerr, label="C 135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(c_combined_x, c_combined_polarization, yerr=c_combined_error, label="C Combined", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[5])
plt.errorbar(c_0.x, (c_0.y-c_combined_polarization)/np.sqrt(np.abs(c_0.yerr**2-c_combined_error**2)), label="0 deg - All angle C", marker="o",markersize=5, linestyle="None")
plt.errorbar(c_45.x, (c_45.y-c_combined_polarization)/np.sqrt(np.abs(c_45.yerr**2-c_combined_error**2)), label="45 deg - All angle C", marker="o",markersize=5, linestyle="None")
plt.errorbar(c_90.x, (c_90.y-c_combined_polarization)/np.sqrt(np.abs(c_90.yerr**2-c_combined_error**2)), label="90 deg - All angle C", marker="o",markersize=5, linestyle="None")
plt.errorbar(c_135.x, (c_135.y-c_combined_polarization)/np.sqrt(np.abs(c_135.yerr**2-c_combined_error**2)), label="135 deg - All angle C", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.axes(axs[6])
plt.errorbar(all_0.x+0.02, all_0.y, yerr=all_0.yerr, label="0 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_45.x+0.04, all_45.y, yerr=all_45.yerr, label="45 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_90.x+0.06, all_90.y, yerr=all_90.yerr, label="90 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_135.x+0.08, all_135.y, yerr=all_135.yerr, label="135 deg", marker="o",markersize=2, linestyle="None")
plt.errorbar(all_combined_x, all_combined_polarization, yerr=all_combined_error, label="0+45+90+135", marker="o",markersize=2, linestyle="None")
plt.ylim(0,1)
plt.legend()

plt.axes(axs[7])
plt.errorbar(all_0.x, (all_0.y-all_combined_polarization)/np.sqrt(np.abs(all_0.yerr**2-all_combined_error**2)), label="0 deg - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_45.x, (all_45.y-all_combined_polarization)/np.sqrt(np.abs(all_45.yerr**2-all_combined_error**2)), label="45 deg - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_90.x, (all_90.y-all_combined_polarization)/np.sqrt(np.abs(all_90.yerr**2-all_combined_error**2)), label="90 deg - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.errorbar(all_135.x, (all_135.y-all_combined_polarization)/np.sqrt(np.abs(all_135.yerr**2-all_combined_error**2)), label="135 deg - All Target Combined", marker="o",markersize=5, linestyle="None")
plt.ylim(-5,5)
plt.xlim(5.8,10.8)
plt.fill_between([5.8,10.8], -1, 1, color="green", alpha=0.3)
plt.fill_between([5.8,10.8], 1, 4, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], -4, -1, color="yellow", alpha=0.3)
plt.fill_between([5.8,10.8], 4, 10, color="red", alpha=0.3)
plt.fill_between([5.8,10.8], -10, -4, color="red", alpha=0.3)

plt.savefig("polarization_target.png")
plt.close()

plt.figure(figsize=(10,6), dpi=300)
plt.errorbar(all_combined_x, all_combined_polarization, yerr=all_combined_error, marker="o",markersize=4, linestyle="None")

lumi_array = lumi(all_combined_x-0.1, all_combined_x+0.1, 28.0)

index_mask_one = np.where((all_combined_x - 6.0)*(all_combined_x - 7.8) < 0.0)
index_mask_two = np.where((all_combined_x - 7.8)*(all_combined_x - 8.8) < 0.0)
index_mask_three = np.where(all_combined_x > 8.8)

# avg_polarization_one        = np.mean(all_combined_polarization[np.where(all_combined_x < 7.8)])
# avg_polarization_one_err    = np.sqrt(np.sum(all_combined_error[np.where(all_combined_x < 7.8)]**2))/len(all_combined_error[np.where(all_combined_x < 7.8)])
# avg_polarization_two        = np.mean(all_combined_polarization[np.where((all_combined_x - 7.8)*(all_combined_x - 8.8) < 0.0)])
# avg_polarization_two_err    = np.sqrt(np.sum(all_combined_error[np.where((all_combined_x - 7.8)*(all_combined_x - 8.8) < 0.0)]**2))/len(all_combined_error[np.where((all_combined_x - 7.8)*(all_combined_x - 8.8) < 0.0)])
# avg_polarization_three      = np.mean(all_combined_polarization[np.where(all_combined_x > 8.8)])
# avg_polarization_three_err  = np.sqrt(np.sum(all_combined_error[np.where(all_combined_x > 8.8)]**2))/len(all_combined_error[np.where(all_combined_x > 8.8)])

avg_polarization_one        = np.average(all_combined_polarization[index_mask_one], weights=lumi_array[index_mask_one])
avg_polarization_one_err    = np.sqrt(np.sum((lumi_array[index_mask_one]*all_combined_error[index_mask_one])**2))/np.sum(lumi_array[index_mask_one])
avg_polarization_two        = np.average(all_combined_polarization[index_mask_two], weights=lumi_array[index_mask_two])
avg_polarization_two_err    = np.sqrt(np.sum((lumi_array[index_mask_two]*all_combined_error[index_mask_two])**2))/np.sum(lumi_array[index_mask_two])
avg_polarization_three      = np.average(all_combined_polarization[index_mask_three], weights=lumi_array[index_mask_three])
avg_polarization_three_err  = np.sqrt(np.sum((lumi_array[index_mask_three]*all_combined_error[index_mask_three])**2))/np.sum(lumi_array[index_mask_three])

plt.errorbar(6.8, avg_polarization_one, yerr=avg_polarization_one_err, xerr = 1.0, marker="o", markersize=5, color="red")
plt.errorbar(8.3, avg_polarization_two, yerr=avg_polarization_two_err, xerr = 0.5, marker="o", markersize=5, color="blue")
plt.errorbar(9.8, avg_polarization_three, yerr=avg_polarization_three_err, xerr = 1.0, marker="o", markersize=5, color="green")

plt.xlabel("Photon Energy (GeV)")
plt.ylabel("Polarization")
plt.savefig("polarization_average.png")
plt.close()