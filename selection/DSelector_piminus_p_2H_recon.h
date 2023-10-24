#ifndef DSelector_piminus_p_2H_recon_h
#define DSelector_piminus_p_2H_recon_h

#include <iostream>
#include <string>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"

class DSelector_piminus_p_2H_recon : public DSelector
{
	public:
    
		DSelector_piminus_p_2H_recon(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_piminus_p_2H_recon(){}

		void   Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:
    
		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool   dIsPolarizedFlag;
		bool   dIsPARAFlag;

		// MONTE CARLO INFORMATION
        bool dIsMC;

		// PARTICLE WRAPPERS
		DParticleComboStep*      dStep0Wrapper;
		DBeamParticle*           dComboBeamWrapper;
		DChargedTrackHypothesis* dPiMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;

		// CUSTOM HISTOGRAMS: BEFORE THE CUTS
        TH1F* dHist_NumUnusedTracks_Measured_Before;
        TH1F* dHist_NumUnusedShowers_Measured_Before;
        TH1F* dHist_PiMinusPIDFOM_Measured_Before;
        TH1F* dHist_ProtonPIDFOM_Measured_Before;
        TH1F* dHist_PhotonEnergy_Measured_Before;
        TH1F* dHist_ConfidenceLevel_KinFit_Before;
        TH1F* dHist_VertexZ_KinFit_Before;
        TH2F* dHist_VertexXY_KinFit_Before;
        TH1F* dHist_MissingPMinus_KinFit_Before;
        TH1F* dHist_MissingMomentum_KinFit_Before;
        
        // CUSTOM HISTOGRAMS: CUT EFFECTS
        TH1F* dHist_CutEffect_NumUnusedTracks;
        TH1F* dHist_CutEffect_NumUnusedShowers;
        TH1F* dHist_CutEffect_TrackPIDFOM;
        TH1F* dHist_CutEffect_PhotonEnergy;
        TH1F* dHist_CutEffect_ConfidenceLevel;
        TH1F* dHist_CutEffect_CommonVertex;
        TH1F* dHist_CutEffect_MissingMomentum;
         
        // CUSTOM HISTOGRAMS: AFTER THE CUTS
        TH1F* dHist_NumUnusedTracks_Measured_After;
        TH1F* dHist_NumUnusedShowers_Measured_After;
        TH1F* dHist_PiMinusPIDFOM_Measured_After;
        TH1F* dHist_ProtonPIDFOM_Measured_After;
        TH2F* dHist_PiMinusPVsdEdx_Measured_After;
        TH2F* dHist_ProtonPVsdEdx_Measured_After;
        TH1F* dHist_PhotonEnergy_Measured_After;
        TH1F* dHist_PhotonTiming_Measured_After;
        TH1F* dHist_ConfidenceLevel_KinFit_After;
        TH1F* dHist_ChiSquarePerNDF_KinFit_After;
        TH1F* dHist_VertexZ_KinFit_After;
        TH2F* dHist_VertexXY_KinFit_After; 
        TH1F* dHist_MissingPMinus_Measured_After;
        TH1F* dHist_MissingPMinus_KinFit_After;
        TH1F* dHist_MissingMassSquared_Measured_After;
        TH1F* dHist_MissingMassSquared_KinFit_After;
        TH1F* dHist_MissingMomentum_Measured_After;
        TH1F* dHist_MissingMomentum_KinFit_After;
        TH1F* dHist_SqrtS_Measured_After;
        TH1F* dHist_SqrtS_KinFit_After;
        TH1F* dHist_MinusT_Measured_After;
        TH1F* dHist_MinusT_KinFit_After;
        TH1F* dHist_MinusU_Measured_After;
        TH1F* dHist_MinusU_KinFit_After;
        TH1F* dHist_ThetaCM_Measured_After;
        TH1F* dHist_ThetaCM_KinFit_After;
        TH1F* dHist_Coplanarity_Measured_After;
        TH1F* dHist_Coplanarity_KinFit_After;
        TH2F* dHist_PiMinusPVsTheta_Measured_After;
        TH2F* dHist_PiMinusPVsTheta_KinFit_After;
        TH2F* dHist_ProtonPVsTheta_Measured_After;
        TH2F* dHist_ProtonPVsTheta_KinFit_After;

	ClassDef(DSelector_piminus_p_2H_recon, 0);
};

void DSelector_piminus_p_2H_recon::Get_ComboWrappers(void)
{
	dStep0Wrapper     = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dPiMinusWrapper   = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dProtonWrapper    = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
}

#endif
