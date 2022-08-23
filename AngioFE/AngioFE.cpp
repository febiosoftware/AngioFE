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
	REGISTER_FECORE_CLASS(AngioFETask, "angio");
	REGISTER_FECORE_CLASS(FEAngioMaterial, "angio_mat");
	REGISTER_FECORE_CLASS(CommonAngioProperties, "angio_properties");

	// Distribution Classes
	REGISTER_FECORE_CLASS(FENormalDistribution, "normal_distribution");
	REGISTER_FECORE_CLASS(FEUniformDistribution, "uniform_distribution");
	REGISTER_FECORE_CLASS(FEExponentialDistribution, "exponential_distribution");
	REGISTER_FECORE_CLASS(FECauchyDistribution, "cauchy_distribution");
	REGISTER_FECORE_CLASS(FEChiSquaredDistribution, "chi_squared_distribution");
	REGISTER_FECORE_CLASS(FEWeibullDistribution, "weibull_distribution");
	REGISTER_FECORE_CLASS(FEGammaDistribution, "gamma_distribution");
	REGISTER_FECORE_CLASS(FEFixedDistribution, "fixed_distribution");
	REGISTER_FECORE_CLASS(FEEllipticalDistribution, "elliptical_distribution");
	REGISTER_FECORE_CLASS(FEFisherDistribution, "fisher_distribution");
	REGISTER_FECORE_CLASS(FEPrescribedDistribution, "prescribed_distribution");

	// Stress Policy Classes
	REGISTER_FECORE_CLASS(SigmoidAngioStressPolicy, "sigmoid_angio_stress_policy");
	REGISTER_FECORE_CLASS(SigmoidDensAngioStressPolicy, "sigmoid_dens_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveVelAngioStressPolicy, "load_curve_vel_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveAngioStressPolicy, "load_curve_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveDenAngioStressPolicy, "load_curve_den_angio_stress_policy");
	REGISTER_FECORE_CLASS(LoadCurveRefDenAngioStressPolicy, "load_curve_ref_den_angio_stress_policy");
	REGISTER_FECORE_CLASS(GrownSegmentsAngioStressPolicy, "grown_segments_angio_stress_policy");
	REGISTER_FECORE_CLASS(GrownSegmentsVelAngioStressPolicy, "grown_segments_vel_angio_stress_policy");

	// Seeder Classes
	REGISTER_FECORE_CLASS(ByElementFragmentSeeder, "by_element_fragment_seeder");
	REGISTER_FECORE_CLASS(ByElementFragmentSeederBiDirectional, "by_element_fragment_seeder_bidirectional");
	REGISTER_FECORE_CLASS(ByElementSetFragmentSeederBiDirectional, "by_element_set_fragment_seeder_bidirectional");
	REGISTER_FECORE_CLASS(ByVolumeFragmentSeeder, "by_volume_fragment_seeder");
	REGISTER_FECORE_CLASS(ByVolumeFragmentSeederBiDirectional, "by_volume_fragment_seeder_bidirectional");
	REGISTER_FECORE_CLASS(SingleCellSeeder, "single_cell_seeder");

	//InitalModifiers any modifiers for angio_elements this will be run only once before much else happens
	REGISTER_FECORE_CLASS(InitialModifierManager, "im_manager");
	REGISTER_FECORE_CLASS(NodeDataInterpolationManager, "nodedata_interpolation_manager");
	REGISTER_FECORE_CLASS(FiberRandomizer, "fiber_randomizer");
	REGISTER_FECORE_CLASS(DiscreteFiberEFDRandomizer, "discrete_fiber_efd_randomizer");
	REGISTER_FECORE_CLASS(EFDFiberInitializer, "efd_initializer");
	REGISTER_FECORE_CLASS(RotEFDFiberInitializer, "rot_efd_initializer");
	REGISTER_FECORE_CLASS(DensityInitializer, "density_initializer");
	REGISTER_FECORE_CLASS(RepulseInitializer, "repulse_value_initializer");
	REGISTER_FECORE_CLASS(DensityValuesNodeDataInterpolation, "ref_ecm_density");
	REGISTER_FECORE_CLASS(RepulseValuesNodeDataInterpolation, "repulse_value");


	// SegmentVelocityModifier CLasses
	REGISTER_FECORE_CLASS(SegmentGrowthVelocityManager, "segment_growth_velocity_manager");
	REGISTER_FECORE_CLASS(SegmentVelocityModifier, "segment_velocity_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityDensityScaleModifier, "segment_velocity_density_scale_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityRefDensityScaleModifier, "segment_velocity_ref_density_scale_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityDensityFAScaleModifier, "segment_velocity_density_fa_scale_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocity3PModifier, "segment_velocity_3P_modifier");
	REGISTER_FECORE_CLASS(SegmentVelocityFAModifier, "segment_velocity_fa_modifier");
	REGISTER_FECORE_CLASS(SigmoidSegmentVelocity, "sigmoid_segment_velocity");
	REGISTER_FECORE_CLASS(SigmoidAdjustedSegmentVelocity, "sigmoid_adjusted_segment_velocity");
	REGISTER_FECORE_CLASS(GompertzSegmentVelocity, "gompertz_segment_velocity");

	// PSC Classes
	REGISTER_FECORE_CLASS(ProtoPreviousSegmentContributionManager, "proto_previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS(ProtoPreviousSegmentPSC, "proto_previous_segment_psc");
	REGISTER_FECORE_CLASS(PreviousSegmentContributionManager, "previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS(PreviousSegmentPSC, "previous_segment_psc");

	// PDD Classes
	REGISTER_FECORE_CLASS(PositionDependentDirectionManager, "position_dependent_direction_manager");
	REGISTER_FECORE_CLASS(FiberPDD, "fiber_pdd");
	REGISTER_FECORE_CLASS(FractionalAnisotropyPDD, "fractional_anisotropy_pdd");
	REGISTER_FECORE_CLASS(LaGrangePStrainPDD, "lagrange_principal_pdd");
	REGISTER_FECORE_CLASS(AnastamosisPDD, "anastamosis_pdd");
	REGISTER_FECORE_CLASS(ECMDensityGradientPDD, "ecm_density_gradient_pdd");
	REGISTER_FECORE_CLASS(RepulsePDD, "repulse_pdd");
	REGISTER_FECORE_CLASS(ConcentrationGradientPDD, "concentration_gradient_pdd");
	REGISTER_FECORE_CLASS(FisherConcentrationGradientPDD, "fisher_concentration_gradient_pdd");

	// ProtoPDD Classes
	REGISTER_FECORE_CLASS(ProtoPositionDependentDirectionManager, "proto_position_dependent_direction_manager");
	REGISTER_FECORE_CLASS(ProtoFiberPDD, "proto_fiber_pdd");
	REGISTER_FECORE_CLASS(ProtoFractionalAnisotropyPDD, "proto_fractional_anisotropy_pdd");
	REGISTER_FECORE_CLASS(ProtoAnastamosisPDD, "proto_anastamosis_pdd");
	REGISTER_FECORE_CLASS(ProtoRepulsePDD, "proto_repulse_pdd");

	// ContributionMix Classes
	REGISTER_FECORE_CLASS(ProtoContributionMixManager, "proto_contribution_mix_manager");
	REGISTER_FECORE_CLASS(ContributionMixManager, "contribution_mix_manager");
	REGISTER_FECORE_CLASS(ProtoPSCPDDContributionMix, "proto_psc_pdd_contribution_mix");
	REGISTER_FECORE_CLASS(PSCPDDContributionMix, "psc_pdd_contribution_mix");
	REGISTER_FECORE_CLASS(DensFAContributionMix, "density_FA_contribution_mix");

	// Species Manager Classes
	REGISTER_FECORE_CLASS(CellSpeciesManager, "species_manager");
	REGISTER_FECORE_CLASS(CellSBM, "SBM");
	REGISTER_FECORE_CLASS(CellSolute, "Solute");
	REGISTER_FECORE_CLASS(CellReactionManager, "cell_reaction_manager");
	REGISTER_FECORE_CLASS(FECellReaction, "cell_reaction");
	REGISTER_FECORE_CLASS(FECellReactionRateConst, "cell constant reaction rate");
	REGISTER_FECORE_CLASS(FECellMassActionForward, "cell mass-action-forward");
	REGISTER_FECORE_CLASS(FECellMassActionForwardConstant, "cell mass-action-forward const");
	REGISTER_FECORE_CLASS(FECellMassActionForwardEffective, "cell mass-action-forward-effective");
	REGISTER_FECORE_CLASS(FECellMassActionReversible, "cell mass-action-reversible");
	REGISTER_FECORE_CLASS(FECellMassActionReversibleEffective, "cell mass-action-reversible-effective");
	REGISTER_FECORE_CLASS(FECellMichaelisMenten, "cell michaelis-menten");
	REGISTER_FECORE_CLASS(FECellInternalization, "cell internalization");
	REGISTER_FECORE_CLASS(FECellInternalizationConstant, "cell internalization constant");
	REGISTER_FECORE_CLASS(FECellSecretion, "cell secretion");
	REGISTER_FECORE_CLASS(FECellSecretionConstant, "cell secretion constant");


	// Mix Methods
	REGISTER_FECORE_CLASS(LinInterp, "LinInterp");
	REGISTER_FECORE_CLASS(LinRot, "LinRot");

	// Interpolation From Gauss Points to a Local Position
	REGISTER_FECORE_CLASS(PerElementVI, "per_element_vi");

	//BranchPolicy
	REGISTER_FECORE_CLASS(DelayedBranchingPolicyEFD, "delayed_branching_policy_efd");

	//branch related classes
	REGISTER_FECORE_CLASS(AzimuthAngleProbabilityDistribution, "azimuth_angle_probability_distribution");
	REGISTER_FECORE_CLASS(ZenithAngleProbabilityDistribution, "zenith_angle_probability_distribution");

	// Plot classes
	REGISTER_FECORE_CLASS(FEPlotAngioStress, "angio stress");
	REGISTER_FECORE_CLASS(FEPlotVesselStress, "vessel stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixStress, "matrix stress");
	REGISTER_FECORE_CLASS(FEPlotVesselWeight, "vessel weight");
	REGISTER_FECORE_CLASS(FEPlotMatrixWeight, "matrix weight");
	REGISTER_FECORE_CLASS(FEPlotVascularDensity, "vascular density");
	REGISTER_FECORE_CLASS(FEPlotMatrixTangent, "matrix tangent");
	REGISTER_FECORE_CLASS(FEPlotMatrixViscoStress, "matrix visco stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixElasticStress, "matrix elastic stress");

	REGISTER_FECORE_CLASS(FEPlotAngioECMDensity, "angio ECM density");
	REGISTER_FECORE_CLASS(FEPlotAngioRepulseVal, "angio repulse value");
	REGISTER_FECORE_CLASS(FEPlotAngioSPD, "angio SPD");
	REGISTER_FECORE_CLASS(FEPlotAngioFractionalAnisotropy, "angio fractional anisotropy");
	REGISTER_FECORE_CLASS(FEPlotBranches, "branch_count");
	REGISTER_FECORE_CLASS(FEPlotAnastamoses, "anastamoses");
	REGISTER_FECORE_CLASS(FEPlotSegmentLength, "segment_length");
	REGISTER_FECORE_CLASS(FEPlotRefSegmentLength, "reference frame segment length");
	REGISTER_FECORE_CLASS(FEPlotPrimaryVesselDirection, "primary segment direction");
	REGISTER_FECORE_CLASS(FEPlotAngioFiberDirection, "angio fiber direction");

	REGISTER_FECORE_CLASS(FEPlotMatrixElastic_m_Q, "matrix elastic mQ");

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

