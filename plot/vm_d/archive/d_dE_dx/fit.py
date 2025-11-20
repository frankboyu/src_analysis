import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import ROOT
from matplotlib.backends.backend_pdf import PdfPages

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

    def plotHeatmap(self,kill_zeros=True,log_scale=False,**kwargs):
        z = self.z
        if (kill_zeros):
            z[z==0] = np.nan
        if (log_scale):
            return plt.pcolormesh(self.xedge,self.yedge,np.log10(self.z),**kwargs)
        else:
            return plt.pcolormesh(self.xedge,self.yedge,self.z,**kwargs)

    def projectionX(self,**kwargs):
        return Hist1D(self.TH2.ProjectionX(),**kwargs)

    def projectionY(self,**kwargs):
        return Hist1D(self.TH2.ProjectionY(),**kwargs)

def exponential(x, a, b, c):
    return np.exp(a*x+b)+c

def double_exponential(x, a, b, c, d, e):
    return np.exp(a*x+b)+np.exp(c*x+d)+e

def gaussian(x, a, b, c):
    return a*np.exp(-0.5*((x-b)/c)**2)


file_data = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_exc_data_2H_ver12.root")
hist_data = file_data.get('NominalCut/d_dEdx_cdc_meas_NominalCut')
pdf = PdfPages('fit_cdc.pdf')

dedx_points = np.arange(0, 40, 0.01)
dedx_p_edges = np.array([45, 50, 52.5, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 200], dtype=int)
dedx_p_low = dedx_p_edges[:-1]
dedx_p_high = dedx_p_edges[1:]
dedx_p_centers = (dedx_p_low + dedx_p_high)/2
dedx_p_mean = np.array([])
dedx_p_sigma = np.array([])

for i in range(len(dedx_p_edges)-1):
    plt.figure(figsize=(8,6))
    hist_slice = Hist1D(hist_data.TH2.ProjectionY("_px",int(dedx_p_edges[i]),int(dedx_p_edges[i+1])))
    if (i < 1):
        hist_slice.rebin(4)
    elif (i < 8):
        hist_slice.rebin(2)
    hist_slice.plotPoints()
    p0_count = hist_slice.y.max()
    p0_mean = np.average(hist_slice.x , weights=hist_slice.y)
    popt, pcov = curve_fit(gaussian, hist_slice.x, hist_slice.y, p0=[p0_count+4,p0_mean,4])
    dedx_p_mean = np.append(dedx_p_mean, popt[1])
    dedx_p_sigma = np.append(dedx_p_sigma, abs(popt[2]))
    plt.plot(dedx_points, gaussian(dedx_points, *popt), label="Gaussian Fit", color='red')
    plt.xlabel(r'$dE/dx$ (keV/cm)')
    plt.ylabel('Counts')
    plt.legend()
    pdf.savefig()
    plt.close()

hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')

plt.figure(figsize=(8,6))
hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-2.0*dedx_p_sigma)
p_points = np.arange(0.25, 3, 0.01)
dE_points = exponential(p_points, *popt)
# plt.plot(p_points, dE_points, label="Fit", color='red')
# plt.errorbar(dedx_p_centers/100, dedx_p_mean, yerr=2.0*dedx_p_sigma, xerr=(dedx_p_high-dedx_p_low)/200, fmt='ro', label='Mean', markersize=5)
plt.title('2.0 Sigma Band Fit')
plt.xlim(0,2)
plt.ylim(0,40)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
pdf.savefig()
plt.close()

para_list = []
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-1.0*dedx_p_sigma)
para_list.append(popt)
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-1.5*dedx_p_sigma)
para_list.append(popt)
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-2.0*dedx_p_sigma)
para_list.append(popt)
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-2.5*dedx_p_sigma)
para_list.append(popt)
popt, pcov = curve_fit(exponential, dedx_p_centers/100, dedx_p_mean-3.0*dedx_p_sigma)
para_list.append(popt)

plt.figure(figsize=(8,6))
hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
p_points = np.arange(0.25, 3, 0.01)
dE_points = exponential(p_points, *para_list[0])
plt.plot(p_points, dE_points, label=r'$1.0 \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[0][0], para_list[0][1], para_list[0][2]), color='green')
plt.scatter(dedx_p_centers/100, dedx_p_mean-dedx_p_sigma, color='green', marker='o', s=10)
dE_points = exponential(p_points, *para_list[1])
plt.plot(p_points, dE_points, label=r'$1.5 \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[1][0], para_list[1][1], para_list[1][2]), color='cyan')
plt.scatter(dedx_p_centers/100, dedx_p_mean-1.5*dedx_p_sigma, color='cyan', marker='o', s=10)
dE_points = exponential(p_points, *para_list[2])
plt.plot(p_points, dE_points, label=r'$2.0 \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[2][0], para_list[2][1], para_list[2][2]), color='red')
plt.scatter(dedx_p_centers/100, dedx_p_mean-2.0*dedx_p_sigma, color='red', marker='o', s=10)
dE_points = exponential(p_points, *para_list[3])
plt.plot(p_points, dE_points, label=r'$2.5 \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[3][0], para_list[3][1], para_list[3][2]), color='orange')
plt.scatter(dedx_p_centers/100, dedx_p_mean-2.5*dedx_p_sigma, color='orange', marker='o', s=10)
dE_points = exponential(p_points, *para_list[4])
plt.plot(p_points, dE_points, label=r'$3.0 \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[4][0], para_list[4][1], para_list[4][2]), color='yellow')
plt.scatter(dedx_p_centers/100, dedx_p_mean-3.0*dedx_p_sigma, color='yellow', marker='o', s=10)
plt.xlim(0,2)
plt.ylim(0,40)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
plt.legend()
pdf.savefig()
plt.close()



# points = np.loadtxt('points_st.txt', delimiter=',', skiprows=1)

# popt, pcov = curve_fit(exponential, points[:,0], points[:,1])
# print(popt)

# hist_data = file_data.get('DeuterondEdxST_After')
# hist_data.plotHeatmap(label="ST", vmin=0, vmax=500)
# p_points = np.arange(0.25, 3, 0.01)

# dE_points = exponential(p_points, *popt)
# plt.plot(p_points, dE_points, label="Fit")

# dE_points = np.exp(-1.9*p_points+2.8)+0.6
# plt.plot(p_points, dE_points, label="-1.9, 2.8, 0.6")

# plt.plot([0.4, 0.4], [0, 40], 'r-', label='p cut')

# plt.xlim(0,2)
# plt.ylim(0,20)
# plt.xlabel(r'$p$ (GeV/c)')
# plt.ylabel(r'$dE/dx$ (keV/cm)')
# plt.legend(loc='upper right')

# plt.savefig("fit_st.png", dpi=300)

pdf.close()