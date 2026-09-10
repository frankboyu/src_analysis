import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.integrate import quad
import scipy
import ROOT as root

mass_kaon = 0.493677

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

file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/plot/vm_d/output/plots_vm_d_efficiency_syst_test.pdf")


fig = plt.figure(figsize=(16, 16), dpi=300)
fig.suptitle("Pi+ data efficiency", fontsize=20, y=0.92)
gs = fig.add_gridspec(2, 2)
axs = gs.subplots().flatten()

plt.axes(axs[0])
file_eff_kp     = File("input/recon_efficiency_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2017-01 ver03", fontsize=15)

plt.axes(axs[1])
file_eff_kp     = File("input/recon_efficiency_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[2])
file_eff_kp     = File("input/recon_efficiency_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[3])
file_eff_kp     = File("input/recon_efficiency_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2020-01 ver01", fontsize=15)
plt.colorbar(ax=axs[3])

for i in range(4):
    axs[i].set_ylabel("Momentum [GeV]", fontsize=15)
    axs[i].set_xlabel("Theta [deg]", fontsize=15)
    axs[i].set_xticks(np.arange(0, 22, 2))
    axs[i].tick_params(axis='both', labelsize=15)

file_pdf.savefig()
plt.close()


fig = plt.figure(figsize=(16, 16), dpi=300)
fig.suptitle("Pi- data efficiency", fontsize=20, y=0.92)
gs = fig.add_gridspec(2, 2)
axs = gs.subplots().flatten()

plt.axes(axs[0])
file_eff_km     = File("input/recon_efficiency_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2017-01 ver03", fontsize=15)

plt.axes(axs[1])
file_eff_km     = File("input/recon_efficiency_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[2])
file_eff_km     = File("input/recon_efficiency_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[3])
file_eff_km     = File("input/recon_efficiency_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.plotHeatmap(vmin=0.8,vmax=1.0)
plt.title("2020-01 ver01", fontsize=15)
plt.colorbar(ax=axs[3])

for i in range(4):
    axs[i].set_ylabel("Momentum [GeV]", fontsize=15)
    axs[i].set_xlabel("Theta [deg]", fontsize=15)
    axs[i].set_xticks(np.arange(0, 22, 2))
    axs[i].tick_params(axis='both', labelsize=15)

file_pdf.savefig()
plt.close()

fig = plt.figure(figsize=(16, 16), dpi=300)
fig.suptitle("Pi+ / Pi- data efficiency ratio", fontsize=20, y=0.92)
gs = fig.add_gridspec(2, 2)
axs = gs.subplots().flatten()

plt.axes(axs[0])
file_eff_km     = File("input/recon_efficiency_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.z[efficiency_km.z==0] = np.nan
file_eff_kp     = File("input/recon_efficiency_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.z[efficiency_kp.z==0] = np.nan
efficiency_kp.z = efficiency_kp.z / efficiency_km.z
efficiency_kp.plotHeatmap(vmin=0.90,vmax=1.10)
plt.title("2017-01 ver03", fontsize=15)

plt.axes(axs[1])
file_eff_km     = File("input/recon_efficiency_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.z[efficiency_km.z==0] = np.nan
file_eff_kp     = File("input/recon_efficiency_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.z[efficiency_kp.z==0] = np.nan
efficiency_kp.z = efficiency_kp.z / efficiency_km.z
efficiency_kp.plotHeatmap(vmin=0.90,vmax=1.10)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[2])
file_eff_km     = File("input/recon_efficiency_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.z[efficiency_km.z==0] = np.nan
file_eff_kp     = File("input/recon_efficiency_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.z[efficiency_kp.z==0] = np.nan
efficiency_kp.z = efficiency_kp.z / efficiency_km.z
efficiency_kp.plotHeatmap(vmin=0.90,vmax=1.10)
plt.title("2018-01 ver02", fontsize=15)

plt.axes(axs[3])
file_eff_km     = File("input/recon_efficiency_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficMethod2_PVsTheta_MC")
efficiency_km.z[efficiency_km.z==0] = np.nan
file_eff_kp     = File("input/recon_efficiency_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficMethod2_PVsTheta_MC")
efficiency_kp.z[efficiency_kp.z==0] = np.nan
efficiency_kp.z = efficiency_kp.z / efficiency_km.z
efficiency_kp.plotHeatmap(vmin=0.90,vmax=1.10)
plt.title("2020-01 ver01", fontsize=15)
plt.colorbar(ax=axs[3])


for i in range(4):
    axs[i].set_ylabel("Momentum [GeV]", fontsize=15)
    axs[i].set_xlabel("Theta [deg]", fontsize=15)
    axs[i].set_xticks(np.arange(0, 22, 2))
    axs[i].tick_params(axis='both', labelsize=15)

file_pdf.savefig()
plt.close()

file_pdf.close()