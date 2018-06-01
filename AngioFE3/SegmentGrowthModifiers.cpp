#include "SegmentGrowthModifiers.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterialPoint.h"

vec3d PreviousSegmentPSC::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	return tip->GetDirection(mesh);
}

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{
	/*
	angio->
	// build the element data array
	for (int n = 0; n < 4; n++)
	{
		fiber_at_int_pts[n].clear();
		fiber_at_int_pts[n].resize(NE);

		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& e = sd->Element(i);
			int nint = e.GaussPoints();
			fiber_at_int_pts[n][i].assign(nint, 0.0);
		}
	}


	// this array will store the results
	FESPRProjection map;

	// loop over stress components

	// fill the ED array
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd->Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint * mp = el.GetMaterialPoint(j);
			FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
			vec3d fd = emp->m_F * emp->m_Q * vec3d(1, 0, 0);
			double lambda = fd.unit();

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
			fiber_at_int_pts[3][i][j] = lambda;
		}
	}

	for (int n = 0; n<4; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], fibers_at_nodes[n]);
	}
	*/
}

double PSCPDDContributionMix::ApplyModifiers(double prev, Tip* tip, FEMesh* mesh)
{
	return psc_weight;
}

void PSCPDDContributionMix::Update(FEMesh * mesh)
{
	// TODO: Implement.
}

BEGIN_PARAMETER_LIST(PSCPDDContributionMix, ContributionMix)
ADD_PARAMETER(psc_weight, FE_PARAM_DOUBLE, "psc_weight");
END_PARAMETER_LIST();

double SegmentVelocityModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	return prev *= segment_velocity_over_time;
}

bool SegmentVelocityModifier::Init()
{
	return true;
}

BEGIN_PARAMETER_LIST(SegmentVelocityModifier, SegmentGrowthVelocity)
ADD_PARAMETER(segment_velocity_over_time, FE_PARAM_DOUBLE, "segment_velocity_over_time");
END_PARAMETER_LIST();

double SegmentVelocityDensityScaleModifier::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_element, FEMesh* mesh)
{
	std::vector<double> density_at_integration_points;

	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		FEElasticMaterialPoint* elastic_mp = gauss_point->ExtractData<FEElasticMaterialPoint>();

		density_at_integration_points.push_back(angio_mp->ref_ecm_density * (1.0 / elastic_mp->m_J));
	}

	double density_at_point = interpolation_prop->Interpolate(angio_element->_elem, density_at_integration_points, natural_coords, mesh);
	double density_scale = m_density_scale_factor.x + m_density_scale_factor.y * exp(-m_density_scale_factor.z * density_at_point);

	return density_scale * prev;
}

bool SegmentVelocityDensityScaleModifier::Init()
{
	return true;
}

BEGIN_PARAMETER_LIST(SegmentVelocityDensityScaleModifier, SegmentGrowthVelocity)
ADD_PARAMETER(m_density_scale_factor, FE_PARAM_VEC3D, "density_scale_factor");
END_PARAMETER_LIST();

/*
vec3d FiberPDD::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	std::vector<double> x_data, y_data, z_data;
	std::vector<double> x_data_nodal, y_data_nodal, z_data_nodal;
	for(int i=0; i< tip->angio_element->_elem->GaussPoints();i++)
	{
		FEMaterialPoint * mp = tip->angio_element->_elem->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		vec3d axis(1, 0, 0);
		axis = emp->m_F * emp->m_Q * axis;
		x_data.push_back(axis.x);
		y_data.push_back(axis.y);
		z_data.push_back(axis.z);
	}
	FESolidDomain * se = dynamic_cast<FESolidDomain*>(tip->angio_element->_elem->GetDomain());
	//projection.Project(*se, x_data, x_data_nodal);
	//vec3d fiber_direction =
	return mix(prev, fiber_direction, contribution);
}
*/

// Managers
vec3d PreviousSegmentContributionManager::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	for(int i=0; i < psc_modifiers.size();i++)
	{
		prev = psc_modifiers[i]->ApplyModifiers(prev, tip, mesh);
	}
	return prev;
}


vec3d PositionDependentDirectionManager::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	for (int i = 0; i < pdd_modifiers.size(); i++)
	{
		prev = pdd_modifiers[i]->ApplyModifiers(prev, tip, mesh);
	}
	return prev;
}

double ContributionMixManager::ApplyModifiers(double prev, Tip* tip, FEMesh* mesh)
{
	for (int i = 0; i < cm_modifiers.size(); i++)
	{
		prev = cm_modifiers[i]->ApplyModifiers(prev, tip, mesh);
	}
	return prev;
}

double SegmentGrowthVelocityManager::ApplyModifiers(double prev, vec3d natural_coords, AngioElement* angio_elem, FEMesh* mesh)
{
	for (int i = 0; i < seg_vel_modifiers.size(); i++)
	{
		prev = seg_vel_modifiers[i]->ApplyModifiers(prev, natural_coords, angio_elem, mesh);
	}
	return prev;
}