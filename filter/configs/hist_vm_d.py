import ROOT
import numpy as np
import matplotlib.pyplot as plt
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

    def scale(self,factor):
        self.z *= factor
        self.zerr *= factor

    def plotHeatmap(self,kill_zeros=True,**kwargs):
        z = self.z
        if (kill_zeros):
            z[z==0] = np.nan
        return plt.pcolormesh(self.xedge,self.yedge,self.z,**kwargs)

    def projectionX(self,**kwargs):
        return Hist1D(self.TH2.ProjectionX(),**kwargs)

    def projectionY(self,**kwargs):
        return Hist1D(self.TH2.ProjectionY(),**kwargs)


file_pdf = PdfPages("/work/halld2/home/boyu/src_analysis/filter/output/hist_vm_d.pdf")

#####################################################################################################################################################
file_data = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_exc_data_2H.root")
file_sim = File("/work/halld2/home/boyu/src_analysis/filter/output/filteredhist_phi_d_recon_exc_sim_2H.root")

# K+ kinematics
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('KinematicsCut/kp_kinematics_meas_KinematicsCut')
hist_sim = file_sim.get('KinematicsCut/kp_kinematics_meas_KinematicsCut')
hist_sim.scale(0.017)
plt.axes(axs[0])
hist_data.plotHeatmap()
plt.title("Data")
plt.axes(axs[1])
hist_sim.plotHeatmap()
plt.title("Simulation")
for i in range(2):
    plt.axes(axs[i])
    plt.plot([0,10], [2,2], '-', color='red')
    plt.plot([0.4,0.4], [0,20], '-', color='red')
    plt.xlim(0, 10)
    plt.ylim(0, 20)
    plt.yticks(np.arange(0, 21, 2))
    plt.xlabel(r"$p_{K^+}$ (GeV/c)")
    plt.ylabel(r"$\theta_{K^+}$ (deg)")
file_pdf.savefig()
plt.close()

# K- kinematics
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('KinematicsCut/km_kinematics_meas_KinematicsCut')
hist_sim = file_sim.get('KinematicsCut/km_kinematics_meas_KinematicsCut')
hist_sim.scale(0.017)
plt.axes(axs[0])
hist_data.plotHeatmap()
plt.title("Data")
plt.axes(axs[1])
hist_sim.plotHeatmap()
plt.title("Simulation")
for i in range(2):
    plt.axes(axs[i])
    plt.plot([0,10], [2,2], '-', color='red')
    plt.plot([0.4,0.4], [0,20], '-', color='red')
    plt.xlim(0, 10)
    plt.ylim(0, 20)
    plt.yticks(np.arange(0, 21, 2))
    plt.xlabel(r"$p_{K^-}$ (GeV/c)")
    plt.ylabel(r"$\theta_{K^-}$ (deg)")
file_pdf.savefig()
plt.close()

# CDC dE/dx
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
hist_sim = file_sim.get('dEdxCut/d_dEdx_cdc_meas_dEdxCut')
hist_sim.scale(0.017)
plt.axes(axs[0])
hist_data.plotHeatmap(vmin=0, vmax=50)
plt.title("Data")
plt.axes(axs[1])
hist_sim.plotHeatmap(vmin=0, vmax=50)
plt.title("Simulation")
for i in range(2):
    plt.axes(axs[i])
    plt.plot(np.arange(0.25, 3, 0.01), np.exp(-3.3*np.arange(0.25, 3, 0.01)+4.1)+2.3, color='red')
    plt.xlim(0, 2)
    plt.ylim(0, 40)
    plt.xlabel("$p$ (GeV/c)")
    plt.ylabel("CDC $dE/dx$ (keV/cm)")
file_pdf.savefig()
plt.close()

# ST dE/dx
fig = plt.figure(figsize=(16, 6), dpi=300)
gs = fig.add_gridspec(1, 2)
axs = gs.subplots()
hist_data = file_data.get('dEdxCut/d_dEdx_st_meas_dEdxCut')
hist_sim = file_sim.get('dEdxCut/d_dEdx_st_meas_dEdxCut')
hist_sim.scale(0.017)
plt.axes(axs[0])
hist_data.plotHeatmap(vmin=0, vmax=50)
plt.title("Data")
plt.axes(axs[1])
hist_sim.plotHeatmap(vmin=0, vmax=50)
plt.title("Simulation")
for i in range(2):
    plt.axes(axs[i])
    plt.plot(np.arange(0.25, 3, 0.01), np.exp(-1.9*np.arange(0.25, 3, 0.01)+2.8)+0.6, color='red')
    plt.xlim(0, 2)
    plt.ylim(0, 20)
    plt.xlabel("$p$ (GeV/c)")
    plt.ylabel("ST $dE/dx$ (keV/cm)")
file_pdf.savefig()
plt.close()

# K+K- invariant mass
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('NominalCut/phi_mass_kin_NominalCut')
hist_sim = file_sim.get('NominalCut/phi_mass_kin_NominalCut')
hist_sim.scale(0.017)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='-', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='red')
plt.xlim(0.98, 1.1)
plt.ylim(0, 1400)
plt.xlabel(r"$M_{K^+K^-} (GeV/c^2)$")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# -t
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/minust_kin_AllCut')
hist_sim = file_sim.get('AllCut/minust_kin_AllCut')
hist_sim.scale(0.017)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='red')
# plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
plt.xlim(0, 2)
plt.xlabel(r"$-t (GeV^2/c^4)$")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# Number of unused tracks
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/num_unused_tracks_meas_AllCut')
hist_sim = file_sim.get('AllCut/num_unused_tracks_meas_AllCut')
hist_sim.scale(0.017)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='red')
plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
plt.xlim(0, 5)
plt.xlabel("Number of Unused Tracks")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# Number of unused showers
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/num_unused_showers_meas_AllCut')
hist_sim = file_sim.get('AllCut/num_unused_showers_meas_AllCut')
hist_sim.scale(0.017)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='red')
plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
plt.xlim(0, 10)
plt.xlabel("Number of Unused Showers")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# KinFit Chi2
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/chisq_per_ndf_kin_AllCut')
hist_sim = file_sim.get('AllCut/chisq_per_ndf_kin_AllCut')
hist_sim.scale(0.017)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='None', color='red')
# plt.bar(hist_data.x, hist_data.y, width=1, alpha=0.5, label='Data', color='green')
plt.xlim(0, 10)
plt.xlabel(r"$\chi^{2}$/NDF")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# Missing energy
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/miss_energy_meas_AllCut')
hist_sim = file_sim.get('AllCut/miss_energy_meas_AllCut')
hist_sim.scale(0.02)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='None', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='red')
plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
plt.xlim(-2, 2)
plt.xlabel("Missing Energy (GeV)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

# Missing mass squared
fig = plt.figure(figsize=(8, 6), dpi=300)
hist_data = file_data.get('AllCut/miss_masssquared_meas_AllCut')
hist_sim = file_sim.get('AllCut/miss_masssquared_meas_AllCut')
hist_sim.scale(0.019)
hist_data.plotPoints(label="Data", marker='o', markersize=3, linestyle='-', color='black')
hist_sim.plotPoints(label="Sim", marker='o', markersize=3, linestyle='-', color='red')
# plt.bar(hist_data.x, hist_data.y, width=0.02, alpha=0.5, label='Data', color='green')
plt.xlim(-0.1, 0.1)
plt.xlabel(r"Missing Mass Squared ($GeV^2/c^4$)")
plt.ylabel("Counts")
plt.legend()
file_pdf.savefig()
plt.close()

file_pdf.close()