import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import ROOT

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
    
    def plotHeatmap(self,kill_zeros=True,**kwargs):
        z = self.z
        if (kill_zeros):
            z[z==0] = np.nan
        return plt.pcolormesh(self.xedge,self.yedge,self.z,**kwargs)

    def projectionX(self,**kwargs):
        return Hist1D(self.TH2.ProjectionX(),**kwargs)

    def projectionY(self,**kwargs):
        return Hist1D(self.TH2.ProjectionY(),**kwargs)


# def exponential(x, a, b, c, d, e):
#     return np.exp(a*x+b)+c*x*x+d*x+e

def exponential(x, a, b, c):
    return np.exp(a*x+b)+c

points = np.loadtxt('points_cdc.txt', delimiter=',', skiprows=1)

popt, pcov = curve_fit(exponential, points[:,0], points[:,1])
print(popt)

file_data = File("../../../filter/output/filteredhist_phi_d_recon_exc_data_2H.root")

hist_data = file_data.get('NoCut/d_dEdx_cdc_meas_NoCut')
hist_data.plotHeatmap(label="CDC", vmin=0, vmax=100)
p_points = np.arange(0.25, 3, 0.01)

dE_points = exponential(p_points, *popt)
plt.plot(p_points, dE_points, label="Fit")

dE_points = np.exp(-3.3*p_points+4.1)+2.3
plt.plot(p_points, dE_points, label="-3.3, 4.1, 2.3")

plt.xlim(0,2)
plt.ylim(0,40)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
plt.legend(loc='upper right')
plt.savefig("fit_cdc.png")
plt.close()


points = np.loadtxt('points_st.txt', delimiter=',', skiprows=1)

popt, pcov = curve_fit(exponential, points[:,0], points[:,1])
print(popt)

file_data = File("../../../filter/output/filteredhist_phi_d_recon_exc_data_2H.root")

hist_data = file_data.get('NoCut/d_dEdx_st_meas_NoCut')
hist_data.plotHeatmap(label="ST", vmin=0, vmax=100)
p_points = np.arange(0.25, 3, 0.01)

dE_points = exponential(p_points, *popt)
plt.plot(p_points, dE_points, label="Fit")

dE_points = np.exp(-1.9*p_points+2.8)+0.6
plt.plot(p_points, dE_points, label="-1.9, 2.8, 0.6")

plt.xlim(0,2)
plt.ylim(0,20)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
plt.legend(loc='upper right')

plt.savefig("fit_st.png")