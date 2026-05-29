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
hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
pdf = PdfPages('fit_cdc.pdf')

dedx_points = np.arange(0, 40, 0.01)
dedx_p_edges = np.array([45, 50, 52.5, 55, 60, 65, 70, 75, 80, 90, 100, 120, 140, 200], dtype=int)
dedx_p_low = dedx_p_edges[:-1]
dedx_p_high = dedx_p_edges[1:]
dedx_p_centers = (dedx_p_low + dedx_p_high)/2
mean_value = np.array([])
mean_err = np.array([])
sigma_value = np.array([])
sigma_err = np.array([])
chisquared = np.array([])

for i in range(len(dedx_p_edges)-1):
    plt.figure(figsize=(8,6))
    hist_slice = Hist1D(hist_data.TH2.ProjectionY("_px",int(dedx_p_edges[i]),int(dedx_p_edges[i+1])))
    if (i < 1):
        hist_slice.rebin(4)
    elif (i < 8):
        hist_slice.rebin(2)
    hist_slice.plotPoints(fmt='o', label='Data')
    plt.ylim(0, hist_slice.y.max()*1.2)
    plt.xlim(0, 40)
    fit_low = np.exp(-3.67*(dedx_p_centers[i]/100)+4.48)+2.57
    fit_high = np.exp(-4.58*(dedx_p_centers[i]/100)+5.66)+5.22
    fit_mask = (hist_slice.x > fit_low) & (hist_slice.x < fit_high)
    hist_slice.x = hist_slice.x[fit_mask]
    hist_slice.y = hist_slice.y[fit_mask]
    hist_slice.yerr = hist_slice.yerr[fit_mask]
    popt, pcov = curve_fit(gaussian, hist_slice.x, hist_slice.y, p0=[hist_slice.y.max(),(fit_low+fit_high)/2,(fit_high-fit_low)/4])
    mean_value = np.append(mean_value, popt[1])
    mean_err = np.append(mean_err, np.sqrt(pcov[1][1]))
    sigma_value = np.append(sigma_value, abs(popt[2]))
    sigma_err = np.append(sigma_err, np.sqrt(pcov[2][2]))
    residuals = hist_slice.y - gaussian(hist_slice.x, *popt)
    chisq = np.sum((residuals**2) / gaussian(hist_slice.x, *popt))
    chisquared = np.append(chisquared, chisq)
    plt.plot(dedx_points[(dedx_points > fit_low) & (dedx_points < fit_high)], gaussian(dedx_points[(dedx_points > fit_low) & (dedx_points < fit_high)], *popt), label="Gaussian Fit", color='red')
    plt.text(0.02, 0.95, r'$\mu=%.3f \pm %.3f$ keV/cm' % (popt[1], np.sqrt(pcov[1][1])), transform=plt.gca().transAxes)
    plt.text(0.02, 0.90, r'$\sigma=%.3f \pm %.3f$ keV/cm' % (abs(popt[2]), np.sqrt(pcov[2][2])), transform=plt.gca().transAxes)
    plt.text(0.02, 0.85, r'$\chi^2/ndf=%.2f$' % (chisq/(len(hist_slice.x)-3)), transform=plt.gca().transAxes)
    plt.xlabel(r'$dE/dx$ (keV/cm)')
    plt.ylabel('Counts')
    plt.title(r'$p$ bin: [%.2f, %.2f] GeV/c' % (dedx_p_low[i]/100, dedx_p_high[i]/100))
    plt.legend()
    pdf.savefig()
    plt.close()

hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
p_points = np.arange(0.25, 3, 0.01)
para_list = []
variation_list = np.arange(-3,0.25,0.25)
color_list = ['violet', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red', 'orange', 'yellow', 'green', 'cyan', 'blue', 'violet']
for N in variation_list:
    popt, pcov = curve_fit(exponential, dedx_p_centers/100, mean_value+N*sigma_value)
    para_list.append(popt)

plt.figure(figsize=(8,6))
plt.errorbar(dedx_p_centers/100, mean_value, xerr=(dedx_p_high-dedx_p_low)/200, yerr=mean_err, fmt='k.', label=r'$\mu$')
plt.errorbar(dedx_p_centers/100, sigma_value, xerr=(dedx_p_high-dedx_p_low)/200, yerr=sigma_err, fmt='b.', label=r'$\sigma$')
plt.xlim(0,2)
plt.ylim(0,40)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
plt.legend()
pdf.savefig()
plt.close()

plt.figure(figsize=(8,6))
hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
dE_points = exponential(p_points, *para_list[4])
plt.scatter(dedx_p_centers/100, mean_value-2.0*sigma_value, color='r', marker='o', s=20, label=r'Data points of $\mu-2\sigma$')
plt.plot(p_points, dE_points, label=r'Exponential fit, $p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (para_list[4][0], para_list[4][1], para_list[4][2]), color='y')
plt.xlim(0,2)
plt.ylim(0,40)
plt.xlabel(r'$p$ (GeV/c)')
plt.ylabel(r'$dE/dx$ (keV/cm)')
plt.legend()
pdf.savefig()
plt.close()

plt.figure(figsize=(8,6))
hist_data.plotHeatmap(log_scale=True, vmin=0, vmax=5, cmap='jet')
for i in range(len(variation_list)):
    N = variation_list[i]
    dE_points = exponential(p_points, *para_list[i])
    plt.plot(p_points, dE_points, label=r'$%.2f \sigma, p_1=%.2f, p_2=%.2f, p_3=%.2f$' % (N, para_list[i][0], para_list[i][1], para_list[i][2]), color=color_list[i])
    plt.scatter(dedx_p_centers/100, mean_value+N*sigma_value, color=color_list[i], marker='o', s=10)
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