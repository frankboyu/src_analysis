#ifndef DSelector_gd_kpkmnmissprot_h
#define DSelector_gd_kpkmnmissprot_h

#include <iostream>

#include "DSelector/DSelector.h"
#include "DSelector/DHistogramActions.h"
#include "DSelector/DCutActions.h"

#include "TH1I.h"
#include "TH2I.h"
#include "TH1F.h"
#include "TH2F.h"

class DSelector_gd_kpkmnmissprot : public DSelector
{
	public:

		DSelector_gd_kpkmnmissprot(TTree* locTree = NULL) : DSelector(locTree){}
		virtual ~DSelector_gd_kpkmnmissprot(){}

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
		DNeutralParticleHypothesis* dNeutronWrapper;
		DKinematicData* dMissingProtonWrapper;

        TH1I* dHist_MultiCombo;
        //CUSTOM HISTOGRAMS: BEFORE THE CUTS
        //Combo
        TH1F* dHist_InvariantMass_Measured_Before;
        TH1F* dHist_InvariantMass_KinFit_Before;
        TH1F* dHist_ConfidenceLevel_KinFit_Before;
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
        TH2F* dHist_NeutronPVsTheta_Measured_Before;
        TH2F* dHist_NeutronPVsTheta_KinFit_Before;
        //Undetected particles
        TH2F* dHist_MissingProtonPVsTheta_Measured_Before;
        TH2F* dHist_MissingProtonPVsTheta_KinFit_Before;
        TH2F* dHist_InitialNeutronPVsTheta_Measured_Before;
        TH2F* dHist_InitialNeutronPVsTheta_KinFit_Before;

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
        TH2F* dHist_NeutronPVsTheta_Measured_After;
        TH2F* dHist_NeutronPVsTheta_KinFit_After;
        //Undetected particles
        TH2F* dHist_MissingProtonPVsTheta_Measured_After;
        TH2F* dHist_MissingProtonPVsTheta_KinFit_After;
        TH2F* dHist_InitialNeutronPVsTheta_Measured_After;
        TH2F* dHist_InitialNeutronPVsTheta_KinFit_After;
        
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
        TH2F* dHist_NeutronPVsTheta_Measured_Weighted;
        TH2F* dHist_NeutronPVsTheta_KinFit_Weighted;
        //Undetected particles
        TH2F* dHist_MissingProtonPVsTheta_Measured_Weighted;
        TH2F* dHist_MissingProtonPVsTheta_KinFit_Weighted;
        TH2F* dHist_InitialNeutronPVsTheta_Measured_Weighted;
        TH2F* dHist_InitialNeutronPVsTheta_KinFit_Weighted;

	ClassDef(DSelector_gd_kpkmnmissprot, 0);
};

void DSelector_gd_kpkmnmissprot::Get_ComboWrappers(void)
{
	//Step 0
	dStep0Wrapper = dComboWrapper->Get_ParticleComboStep(0);
	dComboBeamWrapper = static_cast<DBeamParticle*>(dStep0Wrapper->Get_InitialParticle());
	dKPlusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(0));
	dKMinusWrapper = static_cast<DChargedTrackHypothesis*>(dStep0Wrapper->Get_FinalParticle(1));
	dNeutronWrapper = static_cast<DNeutralParticleHypothesis*>(dStep0Wrapper->Get_FinalParticle(2));
	dMissingProtonWrapper = dStep0Wrapper->Get_FinalParticle(3);
}

#endif // DSelector_gd_kpkmnmissprot_h
