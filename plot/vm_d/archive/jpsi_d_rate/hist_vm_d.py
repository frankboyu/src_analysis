import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import ROOT as root
from scipy.integrate import quad
import scipy.stats

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


deuteron_file = File("/work/halld2/home/boyu/src_analysis/selection/output/test/selectedhist_jpsi_d_recon.root")
deuteron_hist = deuteron_file.get('PhotonEnergy')
proton_file = File("/work/halld2/home/boyu/src_analysis/selection/output/test/selectedhist_jpsi_p_recon.root")
proton_hist = proton_file.get('PhotonEnergy')

plt.figure(figsize=(8,6), dpi=300)
deuteron_hist.rebin(10)
proton_hist.rebin(10)
print(proton_hist.y.sum()/deuteron_hist.y.sum())
deuteron_hist.scale(2270/320*295/np.sum(proton_hist.y))
proton_hist.scale(2270/np.sum(proton_hist.y))
deuteron_hist.plotPoints(label=r'$J/\psi \ d$, GlueX LD2',color='blue', marker='.', linestyle='None', )
proton_hist.plotPoints(label=r'$J/\psi \ p$, GlueX LH2',color='orange', marker='.', linestyle='None')
plt.xlabel("Photon Energy (GeV)")
plt.ylabel("Counts / 0.1 GeV")
# plt.ylim(1, 1e4)
plt.legend()
plt.savefig("photon_energy_comparison.png")
plt.close()


deuteron_flux_file = File("jpsi_d.root")
deuteron_flux_hist = deuteron_flux_file.get('hEBeam')
proton_flux_file = File("jpsi_p.root")
proton_flux_hist = proton_flux_file.get('hEBeam')

deuteron_hist = deuteron_file.get('PhotonEnergy')
proton_hist = proton_file.get('PhotonEnergy')

plt.figure(figsize=(8,6), dpi=300)
deuteron_hist.rebin(10)
proton_hist.rebin(10)
deuteron_hist.y = deuteron_hist.y / deuteron_flux_hist.y[10:100]
proton_hist.y = proton_hist.y / proton_flux_hist.y[10:100]
deuteron_hist.yerr = deuteron_hist.yerr / deuteron_flux_hist.y[10:100]
proton_hist.yerr = proton_hist.yerr / proton_flux_hist.y[10:100]
deuteron_hist.plotPoints(label=r'$J/\psi \ d$',color='blue', marker='.', linestyle='None')
proton_hist.plotPoints(label=r'$J/\psi \ p$',color='orange', marker='.', linestyle='None')
plt.xlabel("Photon Energy (GeV)")
plt.ylabel("Efficiency")
plt.ylim(0,1.0)
plt.legend()
plt.savefig("photon_energy_comparison_flux_weighted.png")
plt.close()