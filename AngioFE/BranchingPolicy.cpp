#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"
#include <iostream>
#include "FEBioMech/FEElasticMaterialPoint.h"
#include "FECore/FEDomainMap.h"
#include <FECore/mathalg.h>
#include <unordered_map>

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(BranchPolicy, FEMaterialProperty)
	ADD_PROPERTY(azimuth_angle, "azimuth_angle");
	ADD_PROPERTY(zenith_angle, "zenith_angle");
	ADD_PROPERTY(interpolation_prop, "interpolation_prop");
	END_FECORE_CLASS()

BEGIN_FECORE_CLASS(DelayedBranchingPolicyEFD, BranchPolicy)
	ADD_PROPERTY(l2b, "length_to_branch");
	ADD_PROPERTY(t2e, "time_to_emerge");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(ZenithAngleProbabilityDistribution, ZenithAngle)
	ADD_PROPERTY(angle, "angle");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(AzimuthAngleProbabilityDistribution, AzimuthAngle)
	ADD_PROPERTY(angle, "angle");
END_FECORE_CLASS()

#pragma endregion FECoreClassDefs

void BranchPolicy::AddBranchTipEFD(AngioElement* angio_element, vec3d local_pos, Segment* parent_seg, double start_time, int vessel_id, int buffer_index, FEMesh* mesh)
{
	Tip* branch = new Tip();
	branch->time = start_time;
	branch->seg_time = start_time;
	branch->angio_element = angio_element;
	branch->face = angio_element;
	branch->SetLocalPosition(local_pos, mesh);
	branch->initial_fragment_id = angio_element->_angio_mat->m_pangio->AddFragment();
	branch->direction = GetBranchDirectionEFD(local_pos, parent_seg->Direction(mesh), angio_element, mesh);
	branch->use_direction = true;
	branch->is_branch = true;
	branch->parent = parent_seg;

	angio_element->active_tips[(buffer_index + 1 % 2)].at(angio_element).push_back(branch);
	angio_element->branch_count++;
}

vec3d BranchPolicy::GetBranchDirectionEFD(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh)
{
	auto se = angio_element->_elem;
	int nint = se->GaussPoints();
	/// Get the local SPD ///

	// vector containing the SPD for each gauss point in the element
	//! SL: Ask steve how to make array of this using nint rather than fixed value
	mat3ds SPDs_gausspts[8];
	double H[8];
	// get the local EFD
	for (int i = 0; i < nint; i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		// Get the SPD   
		angio_mp->UpdateSPD();
		mat3ds temp_SPD = angio_mp->angioSPD;
		SPDs_gausspts[i] = temp_SPD * (3.0 / temp_SPD.tr());
		//! TODO: calculate all distances from mp to nodes then normalize. 
		//! Currently assumes equidistance
		H[i] = 1.0 / double(nint);
	}
	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[8];
	//! Get the interpolated SPD from the shape function-weighted Average 
	//! Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_gausspts, H, nint);
	SPD_int = (3.0 / SPD_int.tr()) * SPD_int;

	FEEllipticalDistribution E(this->GetFEModel());
	E.spd = SPD_int;
	E.Init();
	vec3d fiber_direction = E.NextVec(angio_element->_rengine);

	vec3d normal = fiber_direction ^ parent_direction;
	//make the fiber direction be ortogonal to parent direction
	vec3d corrected_fiber_direction = normal ^ parent_direction;
	double z_ang = zenith_angle->GetZenithAngle(local_pos, parent_direction, angio_element);
	quatd zenith_rotation(z_ang,corrected_fiber_direction);
	double a_ang = azimuth_angle->GetAzimuthAngle(local_pos, parent_direction, angio_element);
	quatd azimuth_rotation(a_ang, parent_direction);
	mat3d a_R = azimuth_rotation.RotationMatrix();
	mat3d z_R = zenith_rotation.RotationMatrix();
	vec3d out_branch_dir =	a_R * z_R * parent_direction;
	return out_branch_dir;
}

void DelayedBranchingPolicyEFD::SetupBranchInfo(AngioElement* angio_elem)
{
	angio_elem->branch_info = new DelayBranchInfo();
	DelayBranchInfo* dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
	dbi->length_to_branch = l2b->NextValue(angio_elem->_rengine);
}

double ZenithAngleProbabilityDistribution::GetZenithAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element)
{
	return angle->NextValue(angio_element->_rengine);
}

double AzimuthAngleProbabilityDistribution::GetAzimuthAngle(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element)
{
	return angle->NextValue(angio_element->_rengine);
}