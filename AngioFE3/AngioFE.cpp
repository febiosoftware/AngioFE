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
	febio.SetActiveModule("solid");
	REGISTER_FECORE_CLASS_EXPLICIT(AngioFETask, FETASK_ID, "angio");
	REGISTER_FECORE_CLASS_EXPLICIT(FEAngioMaterial, FEMATERIAL_ID, "angio_mat");
	REGISTER_FECORE_CLASS_EXPLICIT(CommonAngioProperties, FEMATERIAL_ID, "angio_properties");

	// Distribution Classes
	REGISTER_FECORE_CLASS_EXPLICIT(FENormalDistribution, FEMATERIAL_ID, "normal_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEUniformDistribution, FEMATERIAL_ID, "uniform_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEExponentialDistribution, FEMATERIAL_ID, "exponential_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FECauchyDistribution, FEMATERIAL_ID, "cauchy_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEChiSquaredDistribution, FEMATERIAL_ID, "chi_squared_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEWeibullDistribution, FEMATERIAL_ID, "weibull_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEGammaDistribution, FEMATERIAL_ID, "gamma_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEFixedDistribution, FEMATERIAL_ID, "fixed_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEEllipticalDistribution, FEMATERIAL_ID, "elliptical_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEFisherDistribution, FEMATERIAL_ID, "fisher_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPrescribedDistribution, FEMATERIAL_ID, "prescribed_distribution");

	// Stress Policy Classes
	REGISTER_FECORE_CLASS_EXPLICIT(SigmoidAngioStressPolicy, FEMATERIAL_ID, "sigmoid_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(LoadCurveVelAngioStressPolicy, FEMATERIAL_ID, "load_curve_vel_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(LoadCurveAngioStressPolicy, FEMATERIAL_ID, "load_curve_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(LoadCurveDenAngioStressPolicy, FEMATERIAL_ID, "load_curve_den_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(LoadCurveRefDenAngioStressPolicy, FEMATERIAL_ID, "load_curve_ref_den_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(GrownSegmentsAngioStressPolicy, FEMATERIAL_ID, "grown_segments_angio_stress_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(GrownSegmentsVelAngioStressPolicy, FEMATERIAL_ID, "grown_segments_vel_angio_stress_policy");

	// Seeder Classes
	REGISTER_FECORE_CLASS_EXPLICIT(ByElementFragmentSeeder, FEMATERIAL_ID, "by_element_fragment_seeder");
	REGISTER_FECORE_CLASS_EXPLICIT(ByElementFragmentSeederBiDirectional, FEMATERIAL_ID, "by_element_fragment_seeder_bidirectional");
	REGISTER_FECORE_CLASS_EXPLICIT(ByVolumeFragmentSeeder, FEMATERIAL_ID, "by_volume_fragment_seeder");
	REGISTER_FECORE_CLASS_EXPLICIT(ByVolumeFragmentSeederBiDirectional, FEMATERIAL_ID, "by_volume_fragment_seeder_bidirectional");

	//InitalModifiers any modifiers for angio_elements this will be run only once before much else happens
	REGISTER_FECORE_CLASS_EXPLICIT(InitialModifierManager, FEMATERIAL_ID, "im_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(NodeDataInterpolationManager, FEMATERIAL_ID, "nodedata_interpolation_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(FiberRandomizer, FEMATERIAL_ID, "fiber_randomizer");
	REGISTER_FECORE_CLASS_EXPLICIT(DiscreteFiberEFDRandomizer, FEMATERIAL_ID, "discrete_fiber_efd_randomizer");
	REGISTER_FECORE_CLASS_EXPLICIT(EFDFiberInitializer, FEMATERIAL_ID, "efd_initializer");
	REGISTER_FECORE_CLASS_EXPLICIT(DensityInitializer, FEMATERIAL_ID, "density_initializer");
	REGISTER_FECORE_CLASS_EXPLICIT(RepulseInitializer, FEMATERIAL_ID, "repulse_value_initializer");
	REGISTER_FECORE_CLASS_EXPLICIT(DensityValuesNodeDataInterpolation, FEMATERIAL_ID, "ref_ecm_density");
	REGISTER_FECORE_CLASS_EXPLICIT(RepulseValuesNodeDataInterpolation, FEMATERIAL_ID, "repulse_value");


	// SegmentVelocityModifier CLasses
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentGrowthVelocityManager, FEMATERIAL_ID, "segment_growth_velocity_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentVelocityModifier, FEMATERIAL_ID, "segment_velocity_modifier");
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentVelocityDensityScaleModifier, FEMATERIAL_ID, "segment_velocity_density_scale_modifier");
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentVelocityRefDensityScaleModifier, FEMATERIAL_ID, "segment_velocity_ref_density_scale_modifier");
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentVelocity3PModifier, FEMATERIAL_ID, "segment_velocity_3P_modifier");
	REGISTER_FECORE_CLASS_EXPLICIT(SegmentVelocityFAModifier, FEMATERIAL_ID, "segment_velocity_fa_modifier");
	REGISTER_FECORE_CLASS_EXPLICIT(SigmoidSegmentVelocity, FEMATERIAL_ID, "sigmoid_segment_velocity");
	REGISTER_FECORE_CLASS_EXPLICIT(GompertzSegmentVelocity, FEMATERIAL_ID, "gompertz_segment_velocity");

	// PSC Classes
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoPreviousSegmentContributionManager, FEMATERIAL_ID, "proto_previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoPreviousSegmentPSC, FEMATERIAL_ID, "proto_previous_segment_psc");
	REGISTER_FECORE_CLASS_EXPLICIT(PreviousSegmentContributionManager, FEMATERIAL_ID, "previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(PreviousSegmentPSC, FEMATERIAL_ID, "previous_segment_psc");

	// PDD Classes
	REGISTER_FECORE_CLASS_EXPLICIT(PositionDependentDirectionManager, FEMATERIAL_ID, "position_dependent_direction_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(FiberPDD, FEMATERIAL_ID, "fiber_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(FractionalAnisotropyPDD, FEMATERIAL_ID, "fractional_anisotropy_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(LaGrangePStrainPDD, FEMATERIAL_ID, "lagrange_principal_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(AnastamosisPDD, FEMATERIAL_ID, "anastamosis_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(ECMDensityGradientPDD, FEMATERIAL_ID, "ecm_density_gradient_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(RepulsePDD, FEMATERIAL_ID, "repulse_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(ConcentrationGradientPDD, FEMATERIAL_ID, "concentration_gradient_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(FisherConcentrationGradientPDD, FEMATERIAL_ID, "fisher_concentration_gradient_pdd");

	// ProtoPDD Classes
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoPositionDependentDirectionManager, FEMATERIAL_ID, "proto_position_dependent_direction_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoFiberPDD, FEMATERIAL_ID, "proto_fiber_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoFractionalAnisotropyPDD, FEMATERIAL_ID, "proto_fractional_anisotropy_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoAnastamosisPDD, FEMATERIAL_ID, "proto_anastamosis_pdd");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoRepulsePDD, FEMATERIAL_ID, "proto_repulse_pdd");

	// ContributionMix Classes
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoContributionMixManager, FEMATERIAL_ID, "proto_contribution_mix_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(ContributionMixManager, FEMATERIAL_ID, "contribution_mix_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(ProtoPSCPDDContributionMix, FEMATERIAL_ID, "proto_psc_pdd_contribution_mix");
	REGISTER_FECORE_CLASS_EXPLICIT(PSCPDDContributionMix, FEMATERIAL_ID, "psc_pdd_contribution_mix");

	// Species Manager Classes
	REGISTER_FECORE_CLASS_EXPLICIT(CellSpeciesManager, FEMATERIAL_ID, "species_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(CellSBM, FEMATERIAL_ID, "SBM");
	REGISTER_FECORE_CLASS_EXPLICIT(CellSolute, FEMATERIAL_ID, "Solute");
	REGISTER_FECORE_CLASS_EXPLICIT(CellReactionManager, FEMATERIAL_ID, "cell_reaction_manager");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellReaction, FEMATERIAL_ID, "cell_reaction");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellReactionRateConst, FEMATERIAL_ID, "cell constant reaction rate");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellMassActionForward, FEMATERIAL_ID, "cell mass-action-forward");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellMassActionForwardEffective, FEMATERIAL_ID, "cell mass-action-forward-effective");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellMassActionReversible, FEMATERIAL_ID, "cell mass-action-reversible");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellMassActionReversibleEffective, FEMATERIAL_ID, "cell mass-action-reversible-effective");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellMichaelisMenten, FEMATERIAL_ID, "cell michaelis-menten");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellInternalization, FEMATERIAL_ID, "cell internalization");
	REGISTER_FECORE_CLASS_EXPLICIT(FECellSecretion, FEMATERIAL_ID, "cell secretion");


	// Mix Methods
	REGISTER_FECORE_CLASS_EXPLICIT(LinInterp, FEMATERIAL_ID, "LinInterp");
	REGISTER_FECORE_CLASS_EXPLICIT(LinRot, FEMATERIAL_ID, "LinRot");

	// Interpolation From Gauss Points to a Local Position
	REGISTER_FECORE_CLASS_EXPLICIT(PerElementVI, FEMATERIAL_ID, "per_element_vi");

	//BranchPolicy
	REGISTER_FECORE_CLASS_EXPLICIT(DelayedBranchingPolicy, FEMATERIAL_ID, "delayed_branching_policy");
	REGISTER_FECORE_CLASS_EXPLICIT(DelayedBranchingPolicyEFD, FEMATERIAL_ID, "delayed_branching_policy_efd");

	//branch related classes
	REGISTER_FECORE_CLASS_EXPLICIT(AzimuthAngleProbabilityDistribution, FEMATERIAL_ID, "azimuth_angle_probability_distribution");
	REGISTER_FECORE_CLASS_EXPLICIT(ZenithAngleProbabilityDistribution, FEMATERIAL_ID, "zenith_angle_probability_distribution");

	// Plot classes
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioStress, FEPLOTDATA_ID, "angio stress");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotVesselStress, FEPLOTDATA_ID, "vessel stress");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixStress, FEPLOTDATA_ID, "matrix stress");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotVesselWeight, FEPLOTDATA_ID, "vessel weight");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixWeight, FEPLOTDATA_ID, "matrix weight");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixTangent, FEPLOTDATA_ID, "matrix tangent");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixViscoStress, FEPLOTDATA_ID, "matrix visco stress");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixElasticStress, FEPLOTDATA_ID, "matrix elastic stress");

	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioECMDensity, FEPLOTDATA_ID, "angio ECM density");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioRepulseVal, FEPLOTDATA_ID, "angio repulse value");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioSPD, FEPLOTDATA_ID, "angio SPD");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioFractionalAnisotropy, FEPLOTDATA_ID, "angio fractional anisotropy");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotBranches, FEPLOTDATA_ID, "branch_count");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAnastamoses, FEPLOTDATA_ID, "anastamoses");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotSegmentLength, FEPLOTDATA_ID, "segment_length");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotRefSegmentLength, FEPLOTDATA_ID, "reference frame segment length");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotPrimaryVesselDirection, FEPLOTDATA_ID, "primary segment direction");
	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotAngioFiberDirection, FEPLOTDATA_ID, "angio fiber direction");

	REGISTER_FECORE_CLASS_EXPLICIT(FEPlotMatrixElastic_m_Q, FEPLOTDATA_ID, "matrix elastic mQ");

	//REGISTER_FECORE_CLASS_EXPLICIT(TipDepositionBC, FEBC_ID, "tip_deposition_bc");

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

