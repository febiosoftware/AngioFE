// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

#define FECORE_DLL
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
FEPluginFactory_T<AngioFETask          , FETASK_ID    > angiofe_task_factory("angio"   );
FEPluginFactory_T<FEAngioMaterial      , FEMATERIAL_ID> angio_mat_factory   ("angio_mat"   );
FEPluginFactory_T<CommonAngioProperties, FEMATERIAL_ID> common_angio_properties_factory("angio_properties");

// Distribution Classes
FEPluginFactory_T<FENormalDistribution     , FEMATERIAL_ID> normal_distribution_factory     ("normal_distribution"     );
FEPluginFactory_T<FEUniformDistribution    , FEMATERIAL_ID> uniform_distribution_factory    ("uniform_distribution"    );
FEPluginFactory_T<FEExponentialDistribution, FEMATERIAL_ID> exponential_distribution_factory("exponential_distribution");
FEPluginFactory_T<FECauchyDistribution     , FEMATERIAL_ID> cauchy_distribution_factory     ("cauchy_distribution"     );
FEPluginFactory_T<FEChiSquaredDistribution , FEMATERIAL_ID> chi_squared_distribution_factory("chi_squared_distribution");
FEPluginFactory_T<FEWeibullDistribution    , FEMATERIAL_ID> weibull_distribution_factory    ("weibull_distribution"    );
FEPluginFactory_T<FEGammaDistribution      , FEMATERIAL_ID> gamma_distribution_factory      ("gamma_distribution"      );
FEPluginFactory_T<FEFixedDistribution, FEMATERIAL_ID> fixed_distribution_factory("fixed_distribution");

// Stress Policy Classes
FEPluginFactory_T<SigmoidAngioStressPolicy     , FEMATERIAL_ID> sigmoid_angio_stress_policy_factory("sigmoid_angio_stress_policy");
FEPluginFactory_T<LoadCurveVelAngioStressPolicy, FEMATERIAL_ID> load_curve_vel_angio_stress_policy_factory("load_curve_vel_angio_stress_policy");
FEPluginFactory_T<LoadCurveAngioStressPolicy   , FEMATERIAL_ID> load_curve_angio_stress_policy_factory("load_curve_angio_stress_policy");
FEPluginFactory_T<GrownSegmentsAngioStressPolicy, FEMATERIAL_ID> grown_segments_angio_stress_policy_factory("grown_segments_angio_stress_policy");
FEPluginFactory_T<GrownSegmentsVelAngioStressPolicy, FEMATERIAL_ID> grown_segments_vel_angio_stress_policy_factory("grown_segments_vel_angio_stress_policy");

// Seeder Classes
FEPluginFactory_T<ByElementFragmentSeeder             , FEMATERIAL_ID> by_element_fragment_seeder_factory("by_element_fragment_seeder");
FEPluginFactory_T<ByElementFragmentSeederBiDirectional, FEMATERIAL_ID> by_element_fragment_seeder_bidirectional_factory("by_element_fragment_seeder_bidirectional");
FEPluginFactory_T<ByVolumeFragmentSeeder              , FEMATERIAL_ID> by_volume_fragment_seeder_factory("by_volume_fragment_seeder");
FEPluginFactory_T<ByVolumeFragmentSeederBiDirectional , FEMATERIAL_ID> by_volume_fragment_seeder_bidirectional_factory("by_volume_fragment_seeder_bidirectional");

//InitalModifiers any modifiers for angio_elements this will be run only once before much else happens
FEPluginFactory_T<InitialModifierManager, FEMATERIAL_ID> inital_modifier_manager_factory("im_manager");
FEPluginFactory_T<NodeDataInterpolationManager, FEMATERIAL_ID> nodedata_interpolation_manager_factory("nodedata_interpolation_manager");
FEPluginFactory_T<FiberRandomizer       , FEMATERIAL_ID> fiber_randomizer_factory("fiber_randomizer");
FEPluginFactory_T<DensityInitializer    , FEMATERIAL_ID> density_initializer_factory("density_initializer");
FEPluginFactory_T<DensityValuesNodeDataInterpolation, FEMATERIAL_ID> ecm_ref_density("ref_ecm_density");
FEPluginFactory_T<RepulseValuesNodeDataInterpolation, FEMATERIAL_ID> repulse_value("repulse_value");


// SegmentVelocityModifier CLasses
FEPluginFactory_T<SegmentGrowthVelocityManager       , FEMATERIAL_ID> segment_growth_velocity_manager_factory("segment_growth_velocity_manager");
FEPluginFactory_T<SegmentVelocityModifier            , FEMATERIAL_ID> segment_velocity_modifier_factory("segment_velocity_modifier");
FEPluginFactory_T<SegmentVelocityDensityScaleModifier, FEMATERIAL_ID> segment_velocity_density_scale_modifier_factory("segment_velocity_density_scale_modifier");

// PSC Classes
FEPluginFactory_T<PreviousSegmentContributionManager, FEMATERIAL_ID> previous_segment_contribution_manager_factory("previous_segment_contribution_manager");
FEPluginFactory_T<PreviousSegmentPSC                , FEMATERIAL_ID> previous_segment_psc_factory("previous_segment_psc");

// PDD Classes
FEPluginFactory_T<PositionDependentDirectionManager, FEMATERIAL_ID> position_dependent_direction_manager_factory("position_dependent_direction_manager");
FEPluginFactory_T<FiberPDD                         , FEMATERIAL_ID> fiber_pdd_factory("fiber_pdd");
FEPluginFactory_T<AnastamosisPDD                   , FEMATERIAL_ID> anastamosis_pdd_factory("anastamosis_pdd");
FEPluginFactory_T<ECMDensityGradientPDD, FEMATERIAL_ID> ecm_density_gradient_factory("ecm_density_gradient_pdd");
FEPluginFactory_T<RepulsePDD, FEMATERIAL_ID> repulse_factory("repulse_pdd");
FEPluginFactory_T<ConcentrationGradientPDD, FEMATERIAL_ID> concentration_gradient_factory("concentration_gradient_pdd");

// ContributionMix Classes
FEPluginFactory_T<ContributionMixManager, FEMATERIAL_ID> contribution_mix_manager_factory("contribution_mix_manager");
FEPluginFactory_T<PSCPDDContributionMix , FEMATERIAL_ID> psd_pdd_contribution_mix_factory("psc_pdd_contribution_mix");

// Interpolation From Gauss Points to a Local Position
FEPluginFactory_T<PerElementVI, FEMATERIAL_ID> per_element_vi_factory("per_element_vi");

//BranchPolicy
FEPluginFactory_T<DelayedBranchingPolicy, FEMATERIAL_ID> delayed_branching_policy_factory("delayed_branching_policy");

//branch related classes
FEPluginFactory_T<AzimuthAngleProbabilityDistribution, FEMATERIAL_ID> azimuth_angle_probability_distribution_factory("azimuth_angle_probability_distribution");
FEPluginFactory_T<ZenithAngleProbabilityDistribution, FEMATERIAL_ID> zenith_angle_probability_distribution_factory("zenith_angle_probability_distribution");

/*
FEPluginFactory_T<GrowDirectionModifiers, FEMATERIAL_ID> grow_direction_modifiers_factory("grow_direction_modifiers");


FEPluginFactory_T<Plot2GGP, FEMATERIAL_ID> plot2_ggp_factory("plot2_ggp");
FEPluginFactory_T<GradientPlot2GGP, FEMATERIAL_ID> gradient_plot2_ggp_factory("gradient_plot2_ggp");
FEPluginFactory_T<MatrixConverterGGP, FEMATERIAL_ID> matrix_converter_ggp_factory("matrix_converter_ggp");
FEPluginFactory_T<ForkedGGP, FEMATERIAL_ID> forked_ggp_factory("forked_ggp");
FEPluginFactory_T<MatrixMixGGP, FEMATERIAL_ID> matrix_mix_ggp_factory("matrix_mix_ggp");
FEPluginFactory_T<EigenValuesGGP, FEMATERIAL_ID> eigen_values_ggp_factory("eigen_values_ggp");
FEPluginFactory_T<EigenVectorsGGP, FEMATERIAL_ID> eigen_vectors_ggp_factory("eigen_vectors_ggp");
FEPluginFactory_T<CrossGGP, FEMATERIAL_ID> cross_ggp_factory("cross_ggp");
FEPluginFactory_T<ThresholdGGP, FEMATERIAL_ID> threshold_ggp_factory("threshold_ggp");
FEPluginFactory_T<ArcCosGGP, FEMATERIAL_ID> arc_cos_ggp_factory("arccos_ggp");
FEPluginFactory_T<ArcSinGGP, FEMATERIAL_ID> arc_sin_ggp_factory("arcsin_ggp");
FEPluginFactory_T<CosGGP, FEMATERIAL_ID> cos_ggp_factory("cos_ggp");
FEPluginFactory_T<SinGGP, FEMATERIAL_ID> sin_ggp_factory("sin_ggp");
FEPluginFactory_T<UnitGGP, FEMATERIAL_ID> unit_diagonal_factory("unit_diagonal_ggp");
FEPluginFactory_T<SetterGGP, FEMATERIAL_ID> setter_ggp_factory("setter_ggp");
FEPluginFactory_T<MatrixSetterGGP, FEMATERIAL_ID> matrix_setter_ggp_factory("matrix_setter_ggp");
FEPluginFactory_T<MatrixInverseGGP, FEMATERIAL_ID> matrix_inverse_ggp_factory("matrix_inverse_ggp");
FEPluginFactory_T<AssertGGP, FEMATERIAL_ID> assert_ggp_factory("assert_ggp");
FEPluginFactory_T<NodalDataGGP, FEMATERIAL_ID> nodal_data_ggp_factory("nodal_data_ggp");
FEPluginFactory_T<NodalDataGradientGGP, FEMATERIAL_ID> nodal_data_gradient_ggp_factory("nodal_data_gradient_ggp");
FEPluginFactory_T<DirectionChangeGGP, FEMATERIAL_ID> direction_change_ggp_factory("direction_change_ggp");
*/


// Plot Classes
FEPluginFactory_T<FEPlotAngioStress          , FEPLOTDATA_ID> plot_angio_stress            ("angio stress"          );
FEPluginFactory_T<FEPlotVesselStress         , FEPLOTDATA_ID> plot_vessel_stress           ("vessel stress"         );
FEPluginFactory_T<FEPlotMatrixStress         , FEPLOTDATA_ID> plot_matrix_stress           ("matrix stress"         );
FEPluginFactory_T<FEPlotVesselWeight         , FEPLOTDATA_ID> plot_vessel_weight           ("vessel weight"         );
FEPluginFactory_T<FEPlotMatrixWeight         , FEPLOTDATA_ID> plot_matrix_weight           ("matrix weight"         );
FEPluginFactory_T<FEPlotMatrixTangent        , FEPLOTDATA_ID> plot_matrix_tangent          ("matrix tangent"        );
FEPluginFactory_T<FEPlotMatrixViscoStress    , FEPLOTDATA_ID> plot_matrix_visco_stress     ("matrix visco stress"   );
FEPluginFactory_T<FEPlotMatrixElasticStress  , FEPLOTDATA_ID> plot_matrix_elastic_stress   ("matrix elastic stress" );

FEPluginFactory_T<FEPlotAngioECMDensity      , FEPLOTDATA_ID> plot_angio_ecm               ("angio ECM density"     );
FEPluginFactory_T<FEPlotBranches, FEPLOTDATA_ID> plot_branches("branch_count");
FEPluginFactory_T<FEPlotAnastamoses, FEPLOTDATA_ID> plot_anastamoses("anastamoses");
FEPluginFactory_T<FEPlotSegmentLength, FEPLOTDATA_ID> plot_segment_length("segment_length");
FEPluginFactory_T<FEPlotRefSegmentLength, FEPLOTDATA_ID> plot_reference_frame_segment_length("reference_frame_segment_length");
FEPluginFactory_T<FEPlotPrimaryVesselDirection, FEPLOTDATA_ID> plot_primary_vessel_direction("primary_segment_direction");

FEPluginFactory_T<FEPlotMatrixElastic_m_Q, FEPLOTDATA_ID> plot_matrix_elastic_m_Q("matrix elastic mQ");

/*
FEPluginFactory_T<NullFiberInitializer, FEMATERIAL_ID> null_fiber_initializer("null_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializer, FEMATERIAL_ID> random_fiber_initializer("random_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializerNonMangling, FEMATERIAL_ID> random_fiber_initializer_non_mangling("random_fiber_initializer_non_mangling");
FEPluginFactory_T<ExplicitDistributionsFiberInitializer, FEMATERIAL_ID> explicit_distribution_fiber_initializer("explicit_distribution_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializerPE, FEMATERIAL_ID> random_fiber_initializer_pe("random_fiber_initializer_pe");
FEPluginFactory_T<EllipsoidalFiberInitializer, FEMATERIAL_ID> ellipsoidal_fiber_initializer("ellipsoidal_fiber_initializer");
*/


FEPluginFactory_T<TipDopingBC, FEBC_ID> tip_doping_bc("tip_doping_bc");

//-----------------------------------------------------------------------------
FECORE_EXPORT  unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
}

FECORE_EXPORT  void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 2;
	minor = 1;
	patch = SVNREVISION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  FECoreFactory * PluginGetFactory(int i)
{
	std::vector<FECoreFactory *> addon_classes{ 
		&angiofe_task_factory, &angio_mat_factory,
		&angio_mat_factory,
		//plot classes
		&plot_angio_stress, &plot_angio_stress,
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
		*/

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
		&gamma_distribution_factory,&normal_distribution_factory, &exponential_distribution_factory,
		&uniform_distribution_factory, &fixed_distribution_factory,
		
		//ggp's

		// Segment Velocity Modifiers
		&segment_growth_velocity_manager_factory,
		&segment_velocity_modifier_factory,
		&segment_velocity_density_scale_modifier_factory,

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
		*/
		//other needed items
		&common_angio_properties_factory
		
		//boundary conditions
		,&tip_doping_bc
	};

	if(i < addon_classes.size())
	{
		return addon_classes[i];
	}
	return nullptr;

}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginCleanup()
{

}

