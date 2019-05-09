// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

//#define FECORE_DLL
#define FECORE_API

#include "FECore/FECoreKernel.h"
#include "AngioFETask.h"
#include "FECore/FECoreFactory.h"
#include "FEAngioMaterial.h"
#include "AngioPlot.h"

#ifdef SVN
#include "svnrev.h"
#else
#define SVNREVISION 0
#endif

//-----------------------------------------------------------------------------
FECORE_EXPORT  unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	REGISTER_FECORE_CLASS(AngioFETask, FETASK_ID, "angio");
	REGISTER_FECORE_CLASS(FEAngioMaterial, FEMATERIAL_ID, "angio_mat");
	REGISTER_FECORE_CLASS(CommonAngioProperties, FEMATERIAL_ID, "angio_properties");

	// Distribution Classes
	REGISTER_FECORE_CLASS(FENormalDistribution, FEMATERIAL_ID, "normal_distribution");
	REGISTER_FECORE_CLASS(FEUniformDistribution, FEMATERIAL_ID, "uniform_distribution");
	REGISTER_FECORE_CLASS(FEExponentialDistribution, FEMATERIAL_ID, "exponential_distribution");
	REGISTER_FECORE_CLASS(FECauchyDistribution, FEMATERIAL_ID, "cauchy_distribution");
	REGISTER_FECORE_CLASS(FEChiSquaredDistribution, FEMATERIAL_ID, "chi_squared_distribution");
	REGISTER_FECORE_CLASS(FEWeibullDistribution, FEMATERIAL_ID, "weibull_distribution");
	REGISTER_FECORE_CLASS(FEGammaDistribution, FEMATERIAL_ID, "gamma_distribution");
	REGISTER_FECORE_CLASS(FEFixedDistribution, FEMATERIAL_ID, "fixed_distribution");

	// Stress Policy Classes
	REGISTER_FECORE_CLASS(SigmoidAngioStressPolicy, FEMATERIAL_ID, "sigmoid_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveVelAngioStressPolicy, FEMATERIAL_ID, "load_curve_vel_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveAngioStressPolicy, FEMATERIAL_ID, "load_curve_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveDenAngioStressPolicy, FEMATERIAL_ID, "load_curve_den_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveRefDenAngioStressPolicy, FEMATERIAL_ID, "load_curve_ref_den_angio_stress_policy");
	REGISTER_FECORE_CLASS(GrownSegmentsAngioStressPolicy, FEMATERIAL_ID, "grown_segments_angio_stress_policy");
	REGISTER_FECORE_CLASS(GrownSegmentsVelAngioStressPolicy, FEMATERIAL_ID, "grown_segments_vel_angio_stress_policy");

	// Seeder Classes
	REGISTER_FECORE_CLASS(ByElementFragmentSeeder, FEMATERIAL_ID, "by_element_fragment_seeder");
	REGISTER_FECORE_CLASS(ByElementFragmentSeederBiDirectional, FEMATERIAL_ID, "by_element_fragment_seeder_bidirectional");
	REGISTER_FECORE_CLASS(ByVolumeFragmentSeeder, FEMATERIAL_ID, "by_volume_fragment_seeder");
	REGISTER_FECORE_CLASS(ByVolumeFragmentSeederBiDirectional, FEMATERIAL_ID, "by_volume_fragment_seeder_bidirectional");

	//InitalModifiers any modifiers for angio_elements this will be run only once before much else happens
	REGISTER_FECORE_CLASS(InitialModifierManager, FEMATERIAL_ID, "im_manager");
	REGISTER_FECORE_CLASS(NodeDataInterpolationManager, FEMATERIAL_ID, "nodedata_interpolation_manager");
	REGISTER_FECORE_CLASS(FiberRandomizer, FEMATERIAL_ID, "fiber_randomizer");
	REGISTER_FECORE_CLASS(DensityInitializer, FEMATERIAL_ID, "density_initializer");
	REGISTER_FECORE_CLASS(DensityValuesNodeDataInterpolation, FEMATERIAL_ID, "ref_ecm_density");
	REGISTER_FECORE_CLASS(RepulseValuesNodeDataInterpolation, FEMATERIAL_ID, "repulse_value");


	// SegmentVelocityModifier CLasses
	REGISTER_FECORE_CLASS(SegmentGrowthVelocityManager, FEMATERIAL_ID, "segment_growth_velocity_manager");
	REGISTER_FECORE_CLASS(SegmentVelocityModifier, FEMATERIAL_ID, "segment_velocity_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityDensityScaleModifier, FEMATERIAL_ID, "segment_velocity_density_scale_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityRefDensityScaleModifier, FEMATERIAL_ID, "segment_velocity_ref_density_scale_modifier");
	REGISTER_FECORE_CLASS(SigmoidSegmentVelocity, FEMATERIAL_ID, "sigmoid_segment_velocity");
	REGISTER_FECORE_CLASS(GompertzSegmentVelocity, FEMATERIAL_ID, "gompertz_segment_velocity");

	// PSC Classes
	REGISTER_FECORE_CLASS(PreviousSegmentContributionManager, FEMATERIAL_ID, "previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS(PreviousSegmentPSC, FEMATERIAL_ID, "previous_segment_psc");

	// PDD Classes
	REGISTER_FECORE_CLASS(PositionDependentDirectionManager, FEMATERIAL_ID, "position_dependent_direction_manager");
	REGISTER_FECORE_CLASS(FiberPDD, FEMATERIAL_ID, "fiber_pdd");
	REGISTER_FECORE_CLASS(LaGrangePStrainPDD, FEMATERIAL_ID, "lagrange_principal_pdd");
	REGISTER_FECORE_CLASS(AnastamosisPDD, FEMATERIAL_ID, "anastamosis_pdd");
	REGISTER_FECORE_CLASS(ECMDensityGradientPDD, FEMATERIAL_ID, "ecm_density_gradient_pdd");
	REGISTER_FECORE_CLASS(RepulsePDD, FEMATERIAL_ID, "repulse_pdd");
	REGISTER_FECORE_CLASS(ConcentrationGradientPDD, FEMATERIAL_ID, "concentration_gradient_pdd");

	// ContributionMix Classes
	REGISTER_FECORE_CLASS(ContributionMixManager, FEMATERIAL_ID, "contribution_mix_manager");
	REGISTER_FECORE_CLASS(PSCPDDContributionMix, FEMATERIAL_ID, "psc_pdd_contribution_mix");

	// Species Manager Classes
	REGISTER_FECORE_CLASS(TipSpeciesManager, FEMATERIAL_ID, "SBM_Manager");
	REGISTER_FECORE_CLASS(TipSpecies, FEMATERIAL_ID, "SBM");

	// Mix Methods
	REGISTER_FECORE_CLASS(LinInterp, FEMATERIAL_ID, "LinInterp");
	REGISTER_FECORE_CLASS(LinRot, FEMATERIAL_ID, "LinRot");

	// Interpolation From Gauss Points to a Local Position
	REGISTER_FECORE_CLASS(PerElementVI, FEMATERIAL_ID, "per_element_vi");

	//BranchPolicy
	REGISTER_FECORE_CLASS(DelayedBranchingPolicy, FEMATERIAL_ID, "delayed_branching_policy");

	//branch related classes
	REGISTER_FECORE_CLASS(AzimuthAngleProbabilityDistribution, FEMATERIAL_ID, "azimuth_angle_probability_distribution");
	REGISTER_FECORE_CLASS(ZenithAngleProbabilityDistribution, FEMATERIAL_ID, "zenith_angle_probability_distribution");

	// Plot classes
	REGISTER_FECORE_CLASS(FEPlotAngioStress, FEPLOTDATA_ID, "angio stress");
	REGISTER_FECORE_CLASS(FEPlotVesselStress, FEPLOTDATA_ID, "vessel stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixStress, FEPLOTDATA_ID, "matrix stress");
	REGISTER_FECORE_CLASS(FEPlotVesselWeight, FEPLOTDATA_ID, "vessel weight");
	REGISTER_FECORE_CLASS(FEPlotMatrixWeight, FEPLOTDATA_ID, "matrix weight");
	REGISTER_FECORE_CLASS(FEPlotMatrixTangent, FEPLOTDATA_ID, "matrix tangent");
	REGISTER_FECORE_CLASS(FEPlotMatrixViscoStress, FEPLOTDATA_ID, "matrix visco stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixElasticStress, FEPLOTDATA_ID, "matrix elastic stress");

	REGISTER_FECORE_CLASS(FEPlotAngioECMDensity, FEPLOTDATA_ID, "angio ECM density");
	REGISTER_FECORE_CLASS(FEPlotBranches, FEPLOTDATA_ID, "branch_count");
	REGISTER_FECORE_CLASS(FEPlotAnastamoses, FEPLOTDATA_ID, "anastamoses");
	REGISTER_FECORE_CLASS(FEPlotSegmentLength, FEPLOTDATA_ID, "segment_length");
	REGISTER_FECORE_CLASS(FEPlotRefSegmentLength, FEPLOTDATA_ID, "reference_frame_segment_length");
	REGISTER_FECORE_CLASS(FEPlotPrimaryVesselDirection, FEPLOTDATA_ID, "primary_segment_direction");

	REGISTER_FECORE_CLASS(FEPlotMatrixElastic_m_Q, FEPLOTDATA_ID, "matrix elastic mQ");

	//REGISTER_FECORE_CLASS(TipDepositionBC, FEBC_ID, "tip_deposition_bc");

}

FECORE_EXPORT  void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 3;
	minor = 0;
	patch = SVNREVISION;
}

FECORE_EXPORT  void PluginCleanup()
{

}

