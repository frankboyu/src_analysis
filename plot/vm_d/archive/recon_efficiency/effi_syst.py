import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
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


file_data   = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_recon_data_ver12.root")
file_sim    = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_exc_recon_sim_ver12_07.root")
kp_data = file_data.get("PlotCut/kinematics_cut_kp_PlotCut")
km_data = file_data.get("PlotCut/kinematics_cut_km_PlotCut")
kp_sim  = file_sim .get("NominalCut/kinematics_cut_kp_NominalCut")
km_sim  = file_sim .get("NominalCut/kinematics_cut_km_NominalCut")


file_eff_kp     = File("outEffic_F2018_ver02_misspip.root")
file_eff_km     = File("outEffic_F2018_ver02_misspim.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
weights_kp_data = np.zeros_like(efficiency_kp.z)
weights_km_data = np.zeros_like(efficiency_km.z)
weights_kp_sim  = np.zeros_like(efficiency_kp.z)
weights_km_sim  = np.zeros_like(efficiency_km.z)

for i in range(kp_data.z.shape[0]):
    for j in range(kp_data.z.shape[1]):
        flag = False
        for k in range(efficiency_kp.z.shape[0]):
            for l in range(efficiency_kp.z.shape[1]):
                if (kp_data.xedge[j] >= efficiency_kp.xedge[k] and kp_data.xedge[j+1] <= efficiency_kp.xedge[k+1] and
                    kp_data.yedge[i] >= efficiency_kp.yedge[l] and kp_data.yedge[i+1] <= efficiency_kp.yedge[l+1]):
                    weights_kp_data[k,l] += kp_data.z[i,j]
                    flag = True
                    break
            if (flag):
                break

for i in range(km_data.z.shape[0]):
    for j in range(km_data.z.shape[1]):
        flag = False
        for k in range(efficiency_km.z.shape[0]):
            for l in range(efficiency_km.z.shape[1]):
                if (km_data.xedge[j] >= efficiency_km.xedge[k] and km_data.xedge[j+1] <= efficiency_km.xedge[k+1] and
                    km_data.yedge[i] >= efficiency_km.yedge[l] and km_data.yedge[i+1] <= efficiency_km.yedge[l+1]):
                    weights_km_data[k,l] += km_data.z[i,j]
                    flag = True
                    break
            if (flag):
                break

for i in range(kp_sim.z.shape[0]):
    for j in range(kp_sim.z.shape[1]):
        flag = False
        for k in range(efficiency_kp.z.shape[0]):
            for l in range(efficiency_kp.z.shape[1]):
                if (kp_sim.xedge[j] >= efficiency_kp.xedge[k] and kp_sim.xedge[j+1] <= efficiency_kp.xedge[k+1] and
                    kp_sim.yedge[i] >= efficiency_kp.yedge[l] and kp_sim.yedge[i+1] <= efficiency_kp.yedge[l+1]):
                    weights_kp_sim[k,l] += kp_sim.z[i,j]
                    flag = True
                    break
            if (flag):
                break

for i in range(km_sim.z.shape[0]):
    for j in range(km_sim.z.shape[1]):
        flag = False
        for k in range(efficiency_km.z.shape[0]):
            for l in range(efficiency_km.z.shape[1]):
                if (km_sim.xedge[j] >= efficiency_km.xedge[k] and km_sim.xedge[j+1] <= efficiency_km.xedge[k+1] and
                    km_sim.yedge[i] >= efficiency_km.yedge[l] and km_sim.yedge[i+1] <= efficiency_km.yedge[l+1]):
                    weights_km_sim[k,l] += km_sim.z[i,j]
                    flag = True
                    break
            if (flag):
                break


fig = plt.figure(figsize=(16, 16), dpi=300)
gs = fig.add_gridspec(2, 2)
axs = gs.subplots().flatten()

plt.axes(axs[0])
plt.imshow(weights_kp_data,origin="lower",extent=[efficiency_kp.yedge[0],efficiency_kp.yedge[-1],efficiency_kp.xedge[0],efficiency_kp.xedge[-1]],aspect="auto")
plt.title("K+ Data")
plt.axes(axs[1])
plt.imshow(weights_km_data,origin="lower",extent=[efficiency_km.yedge[0],efficiency_km.yedge[-1],efficiency_km.xedge[0],efficiency_km.xedge[-1]],aspect="auto")
plt.title("K- Data")
plt.axes(axs[2])
plt.imshow(weights_kp_sim,origin="lower",extent=[efficiency_kp.yedge[0],efficiency_kp.yedge[-1],efficiency_kp.xedge[0],efficiency_kp.xedge[-1]],aspect="auto")
plt.title("K+ Sim")
plt.axes(axs[3])
plt.imshow(weights_km_sim,origin="lower",extent=[efficiency_km.yedge[0],efficiency_km.yedge[-1],efficiency_km.xedge[0],efficiency_km.xedge[-1]],aspect="auto")
plt.title("K- Sim")

for i in range(4):
    axs[i].set_xlabel("Momentum [GeV]")
    axs[i].set_ylabel("Theta [deg]")

plt.savefig("weights.png",dpi=300)
plt.close()

print("K+ Simulation weighted average efficiency: ")

file_eff_kp     = File("outEffic_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2017 ver03 Method 1: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2017 ver03 Method 2: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2018 Spring ver02 Method 1: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2018 Spring ver02 Method 2: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2018 Fall ver02 Method 1: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2018 Fall ver02 Method 2: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2020 ver01 Method 1: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

file_eff_kp     = File("outEffic_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2020 ver01 Method 2: ",np.sum(weights_kp_sim*efficiency_kp.z)/np.sum(weights_kp_sim))

print("\n")
print("K- Simulation weighted average efficiency: ")

file_eff_km     = File("outEffic_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2017 ver03 Method 1: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2017 ver03 Method 2: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2018 Spring ver02 Method 1: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2018 Spring ver02 Method 2: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2018 Fall ver02 Method 1: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2018 Fall ver02 Method 2: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2020 ver01 Method 1: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

file_eff_km     = File("outEffic_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2020 ver01 Method 2: ",np.sum(weights_km_sim*efficiency_km.z)/np.sum(weights_km_sim))

print("\n")
print("K+ Data weighted average efficiency: ")

file_eff_kp     = File("outEffic_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2017 ver03 Method 1: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_2017_ver03_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2017 ver03 Method 2: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2018 Spring ver02 Method 1: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_S2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2018 Spring ver02 Method 2: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2018 Fall ver02 Method 1: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_F2018_ver02_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2018 Fall ver02 Method 2: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod1_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
# print("2020 ver01 Method 1: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

file_eff_kp     = File("outEffic_2020_ver01_misspip.root")
efficiency_kp   = file_eff_kp.get("EfficRatioMethod2_PVsTheta")
efficiency_kp.z = np.where(efficiency_kp.z == 0, 1.03, efficiency_kp.z)
print("2020 ver01 Method 2: ",np.sum(weights_kp_data*efficiency_kp.z)/np.sum(weights_kp_data))

print("\n")
print("K- Data weighted average efficiency: ")

file_eff_km     = File("outEffic_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2017 ver03 Method 1: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_2017_ver03_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2017 ver03 Method 2: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2018 Spring ver02 Method 1: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_S2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2018 Spring ver02 Method 2: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2018 Fall ver02 Method 1: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_F2018_ver02_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2018 Fall ver02 Method 2: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod1_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
# print("2020 ver01 Method 1: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))

file_eff_km     = File("outEffic_2020_ver01_misspim.root")
efficiency_km   = file_eff_km.get("EfficRatioMethod2_PVsTheta")
efficiency_km.z = np.where(efficiency_km.z == 0, 1.03, efficiency_km.z)
print("2020 ver01 Method 2: ",np.sum(weights_km_data*efficiency_km.z)/np.sum(weights_km_data))