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

	// PSC Classes
	REGISTER_FECORE_CLASS(PreviousSegmentContributionManager, FEMATERIAL_ID, "previous_segment_contribution_manager");
	REGISTER_FECORE_CLASS(PreviousSegmentPSC, FEMATERIAL_ID, "previous_segment_psc");

	// PDD Classes
	REGISTER_FECORE_CLASS(PositionDependentDirectionManager, FEMATERIAL_ID, "position_dependent_direction_manager");
	REGISTER_FECORE_CLASS(FiberPDD, FEMATERIAL_ID, "fiber_pdd");
	REGISTER_FECORE_CLASS(AnastamosisPDD, FEMATERIAL_ID, "anastamosis_pdd");
	REGISTER_FECORE_CLASS(ECMDensityGradientPDD, FEMATERIAL_ID, "ecm_density_gradient_pdd");
	REGISTER_FECORE_CLASS(RepulsePDD, FEMATERIAL_ID, "repulse_pdd");
	REGISTER_FECORE_CLASS(ConcentrationGradientPDD, FEMATERIAL_ID, "concentration_gradient_pdd");

	// ContributionMix Classes
	REGISTER_FECORE_CLASS(ContributionMixManager, FEMATERIAL_ID, "contribution_mix_manager");
	REGISTER_FECORE_CLASS(PSCPDDContributionMix, FEMATERIAL_ID, "psc_pdd_contribution_mix");

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

	REGISTER_FECORE_CLASS(TipDepositionBC, FEBC_ID, "tip_deposition_bc");

}

FECORE_EXPORT  void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 3;
	minor = 0;
	patch = SVNREVISION;
}

//-----------------------------------------------------------------------------

/*
FECORE_EXPORT  FECoreFactory * PluginGetFactory(int i)
{
	std::vector<FECoreFactory *> addon_classes{
		&angiofe_task_factory,
		&angio_mat_factory,
		//plot classes
		&plot_angio_stress,
		&plot_angio_ecm,
		&plot_branches,
		&plot_anastamoses,
		&plot_segment_length,
		&plot_reference_frame_segment_length,
		&plot_vessel_stress, &plot_matrix_stress,
		&plot_vessel_weight, &plot_matrix_weight, &plot_matrix_tangent, &plot_matrix_visco_stress,
		&plot_matrix_elastic_m_Q, &plot_matrix_elastic_stress, &plot_primary_vessel_direction,

		//angio stress policy classes
		&sigmoid_angio_stress_policy_factory, &load_curve_vel_angio_stress_policy_factory ,&load_curve_angio_stress_policy_factory,
		&grown_segments_angio_stress_policy_factory,
		&grown_segments_vel_angio_stress_policy_factory,

		/*
		//fiber initializers
		&null_fiber_initializer, &random_fiber_initializer,
		&random_fiber_initializer_non_mangling, &explicit_distribution_fiber_initializer,
		&random_fiber_initializer_pe, &ellipsoidal_fiber_initializer,
		*\/

		//branching factories
&delayed_branching_policy_factory,

//branching related classes
&azimuth_angle_probability_distribution_factory,
&zenith_angle_probability_distribution_factory,

//grow direction modifiers

//vessel contribution modifiers

//fragment seeders
&by_element_fragment_seeder_factory,
&by_element_fragment_seeder_bidirectional_factory,
&by_volume_fragment_seeder_factory,
&by_volume_fragment_seeder_bidirectional_factory,

//initial modifiers
&inital_modifier_manager_factory,
&nodedata_interpolation_manager_factory,
&fiber_randomizer_factory,
&density_initializer_factory,
&ecm_ref_density,
&repulse_value,

//boundary conditions

//random distribution
&cauchy_distribution_factory, &chi_squared_distribution_factory, &weibull_distribution_factory,
&gamma_distribution_factory, &normal_distribution_factory, &exponential_distribution_factory,
&uniform_distribution_factory, &fixed_distribution_factory,

//ggp's

// Segment Velocity Modifiers
&segment_growth_velocity_manager_factory,
&segment_velocity_modifier_factory,
&segment_velocity_density_scale_modifier_factory,
&segment_velocity_ref_density_scale_modifier_factory,

// PSC Modifier Classes
&previous_segment_contribution_manager_factory,
&previous_segment_psc_factory,

// PDD Modifier Classes
&position_dependent_direction_manager_factory,
&fiber_pdd_factory,
&anastamosis_pdd_factory,
&ecm_density_gradient_factory,
&repulse_factory,
&concentration_gradient_factory,

// ContributionMix Classes
&contribution_mix_manager_factory,
&psd_pdd_contribution_mix_factory,

// Interpolation Classes
&per_element_vi_factory,


/*
&plot2_ggp_factory, &gradient_plot2_ggp_factory,
&matrix_converter_ggp_factory, &forked_ggp_factory, &cross_ggp_factory,
&threshold_ggp_factory, &nodal_data_ggp_factory, &nodal_data_gradient_ggp_factory,
&arc_cos_ggp_factory, &arc_sin_ggp_factory, &cos_ggp_factory, &sin_ggp_factory,
&matrix_inverse_ggp_factory, &eigen_vectors_ggp_factory, &eigen_values_ggp_factory,
&setter_ggp_factory, &matrix_setter_ggp_factory,
&assert_ggp_factory, &unit_diagonal_factory, &direction_change_ggp_factory,
&matrix_mix_ggp_factory,
*\/
//other needed items
&common_angio_properties_factory

//boundary conditions
, &tip_deposition_bc
	};

	if (i < addon_classes.size())
	{
		return addon_classes[i];
	}
	return nullptr;

}
*/
//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginCleanup()
{

}

