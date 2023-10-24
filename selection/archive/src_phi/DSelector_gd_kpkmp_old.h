#ifndef DSelector_gd_kpkmp_h
#define DSelector_gd_kpkmp_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"

class DSelector_gd_kpkmp : public DSelector
{
	public:

		DSelector_gd_kpkmp(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_gd_kpkmp(){}

		void Init(TTree *tree);
		Bool_t Process(Long64_t entry);

	private:

		void Get_ComboWrappers(void);
		void Finalize(void);

		// BEAM POLARIZATION INFORMATION
		UInt_t dPreviousRunNumber;
		bool dIsPolarizedFlag; //else is AMO
		bool dIsPARAFlag; //else is PERP or AMO

		bool dIsMC;

		// ANALYZE CUT ACTIONS
		DHistogramAction_AnalyzeCutActions* dAnalyzeCutActions;

		//CREATE REACTION-SPECIFIC PARTICLE ARRAYS

		//Step 0
		DParticleComboStep* dStep0Wrapper;
		DBeamParticle* dComboBeamWrapper;
		DChargedTrackHypothesis* dKPlusWrapper;
		DChargedTrackHypothesis* dKMinusWrapper;
		DChargedTrackHypothesis* dProtonWrapper;
		DKinematicData* dMissingNeutronWrapper;

        TH1I* dHist_MultiCombo;
        //CUSTOM HISTOGRAMS: BEFORE THE CUTS
        //Combo
        TH1F* dHist_InvariantMass_Measured_Before;
        TH1F* dHist_InvariantMass_KinFit_Before;
        TH1F* dHist_ConfidenceLevel_KinFit_Before;
        TH1F* dHist_ChiSq_KinFit_Before;
        TH2F* dHist_ChiSqComparison_KinFit_Before;
        TH1F* dHist_MissingMassSquared_Measured_Before;
        TH1F* dHist_MissingMassSquared_KinFit_Before;
        TH1F* dHist_EnergyBalance_Measured_Before;
        TH1F* dHist_EnergyBalance_KinFit_Before;
        TH1F* dHist_VertexZ_KinFit_Before;
        TH2F* dHist_VertexXY_KinFit_Before;
        TH1F* dHist_MinusT_Measured_Before;
        TH1F* dHist_MinusT_KinFit_Before;
        TH1F* dHist_MinusU_Measured_Before;
        TH1F* dHist_MinusU_KinFit_Before;
        TH1F* dHist_Coplanarity_Measured_Before;
        TH1F* dHist_Coplanarity_KinFit_Before;
        //Beam photon
        TH1F* dHist_PhotonEnergy_Measured_Before;
        TH1F* dHist_PhotonTiming_Measured_Before;
        //Detected particles
        TH2F* dHist_KPlusPVsTheta_Measured_Before;
        TH2F* dHist_KPlusPVsTheta_KinFit_Before;
        TH2F* dHist_KMinusPVsTheta_Measured_Before;
        TH2F* dHist_KMinusPVsTheta_KinFit_Before;
        TH2F* dHist_ProtonPVsTheta_Measured_Before;
        TH2F* dHist_ProtonPVsTheta_KinFit_Before;
        //Undetected particles
        TH2F* dHist_MissingNeutronPVsTheta_Measured_Before;
        TH2F* dHist_MissingNeutronPVsTheta_KinFit_Before;
        TH2F* dHist_InitialProtonPVsTheta_Measured_Before;
        TH2F* dHist_InitialProtonPVsTheta_KinFit_Before;

        //CUSTOM HISTOGRAMS: DURING THE CUTS
        TH1F* dHist_InvariantMass_ConfidenceLevelCut;
        TH1F* dHist_InvariantMass_MissingMassSquaredCut;
        TH1F* dHist_InvariantMass_EnergyBalanceCut;
        TH1F* dHist_InvariantMass_CoplanarityCut;
        TH1F* dHist_InvariantMass_UnusedTracksCut;
        TH1F* dHist_InvariantMass_UnusedShowersCut;
        TH1F* dHist_InvariantMass_VertexCut;
        TH1F* dHist_InvariantMass_PhotonEnergyCut;
        TH1F* dHist_InvariantMass_InvariantMassCut;

        //CUSTOM HISTOGRAMS: AFTER THE CUTS
        //Combo
        TH1F* dHist_InvariantMass_Measured_After;
        TH1F* dHist_InvariantMass_KinFit_After;
        TH1F* dHist_ConfidenceLevel_KinFit_After;
        TH1F* dHist_ChiSq_KinFit_After;
        TH2F* dHist_ChiSqComparison_KinFit_After;
        TH1F* dHist_MissingMassSquared_Measured_After;
        TH1F* dHist_MissingMassSquared_KinFit_After;
        TH1F* dHist_EnergyBalance_Measured_After;
        TH1F* dHist_EnergyBalance_KinFit_After;       
        TH1F* dHist_VertexZ_KinFit_After;
        TH2F* dHist_VertexXY_KinFit_After;
        TH1F* dHist_MinusT_Measured_After;
        TH1F* dHist_MinusT_KinFit_After;
        TH1F* dHist_MinusU_Measured_After;
        TH1F* dHist_MinusU_KinFit_After;
        TH1F* dHist_Coplanarity_Measured_After;
        TH1F* dHist_Coplanarity_KinFit_After;
        //Beam photon
        TH1F* dHist_PhotonEnergy_Measured_After;
        TH1F* dHist_PhotonTiming_Measured_After;
        //Detected particles
        TH2F* dHist_KPlusPVsTheta_Measured_After;
        TH2F* dHist_KPlusPVsTheta_KinFit_After;
        TH2F* dHist_KMinusPVsTheta_Measured_After;
        TH2F* dHist_KMinusPVsTheta_KinFit_After;
        TH2F* dHist_ProtonPVsTheta_Measured_After;
        TH2F* dHist_ProtonPVsTheta_KinFit_After;
        //Undetected particles
        TH2F* dHist_MissingNeutronPVsTheta_Measured_After;
        TH2F* dHist_MissingNeutronPVsTheta_KinFit_After;
        TH2F* dHist_InitialProtonPVsTheta_Measured_After;
        TH2F* dHist_InitialProtonPVsTheta_KinFit_After;
        
        //CUSTOM HISTOGRAMS: WITH MULTI-COMBO WEIGHT
        //Combo
        TH1F* dHist_InvariantMass_Measured_Weighted;
        TH1F* dHist_InvariantMass_KinFit_Weighted;
        TH1F* dHist_ConfidenceLevel_KinFit_Weighted;
        TH1F* dHist_MissingMassSquared_Measured_Weighted;
        TH1F* dHist_MissingMassSquared_KinFit_Weighted;
        TH1F* dHist_EnergyBalance_Measured_Weighted;
        TH1F* dHist_EnergyBalance_KinFit_Weighted;        
        TH1F* dHist_VertexZ_KinFit_Weighted;
        TH2F* dHist_VertexXY_KinFit_Weighted;
        TH1F* dHist_MinusT_Measured_Weighted;
        TH1F* dHist_MinusT_KinFit_Weighted;
        TH1F* dHist_MinusU_Measured_Weighted;
        TH1F* dHist_MinusU_KinFit_Weighted;
        TH1F* dHist_Coplanarity_Measured_Weighted;
        TH1F* dHist_Coplanarity_KinFit_Weighted;
        //Beam photon
        TH1F* dHist_PhotonEnergy_Measured_Weighted;
        TH1F* dHist_PhotonTiming_Measured_Weighted;
        //Detected particles
        TH2F* dHist_KPlusPVsTheta_Measured_Weighted;
        TH2F* dHist_KPlusPVsTheta_KinFit_Weighted;
        TH2F* dHist_KMinusPVsTheta_Measured_Weighted;
        TH2F* dHist_KMinusPVsTheta_KinFit_Weighted;
        TH2F* dHist_ProtonPVsTheta_Measured_Weighted;
        TH2F* dHist_ProtonPVsTheta_KinFit_Weighted;
        //Undetected particles
        TH2F* dHist_MissingNeutronPVsTheta_Measured_Weighted;
        TH2F* dHist_MissingNeutronPVsTheta_KinFit_Weighted;
        TH2F* dHist_InitialProtonPVsTheta_Measured_Weighted;
        TH2F* dHist_InitialProtonPVsTheta_KinFit_Weighted;

	ClassDef(DSelector_gd_kpkmp, 0);
};

void DSelector_gd_kpkmp::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dProtonWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dMissingNeutronWrapper = dStep0Wrapper->Get_FinalParticle(3);
}

#endif // DSelector_gd_kpkmp_h
