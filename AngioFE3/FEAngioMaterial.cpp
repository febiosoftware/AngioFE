#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEFiberMaterialPoint.h"
#include "FECore/FEElementTraits.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FEBioMech/FESPRProjection.h"
#include <iostream>
#include "angio3d.h"
#include "AngioElement.h"
#include "Segment.h"
#include "Tip.h"


//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	ADD_PARAMETER2(m_cultureParams.sprout_s_mag, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "a");
	ADD_PARAMETER2(m_cultureParams.sprout_s_range, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "b");
	ADD_PARAMETER2(m_cultureParams.sprout_s_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "N");

	ADD_PARAMETER2(m_cultureParams.m_length_adjustment, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "length_adjustment");
	ADD_PARAMETER2(m_cultureParams.m_vessel_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "vessel_width");
	ADD_PARAMETER(m_cultureParams.growth_length_over_time, FE_PARAM_DOUBLE, "growth_length_over_time");

	ADD_PARAMETER(m_cultureParams.ecm_control, FE_PARAM_INT, "ecm_seeder");
	ADD_PARAMETER2(m_cultureParams.m_matrix_density, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "matrix_density");

	ADD_PARAMETER(m_cultureParams.m_symmetry_plane, FE_PARAM_VEC3D, "symmetryplane");
	//uncategorized variables are incomplete
	ADD_PARAMETER(m_cultureParams.m_composite_material, FE_PARAM_INT, "composite_material");
	ADD_PARAMETER2(m_cultureParams.m_sprout_force, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "sprout_force");

	ADD_PARAMETER2(m_cultureParams.active_tip_threshold, FE_PARAM_INT, FE_RANGE_GREATER_OR_EQUAL(0), "active_tip_threshold");
	ADD_PARAMETER2(m_cultureParams.stress_radius, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "stress_radius");
	ADD_PARAMETER2(m_cultureParams.min_segment_length, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "min_segment_length");
	ADD_PARAMETER(m_cultureParams.m_weight_interpolation, FE_PARAM_DOUBLE, "direction_weight");
	
	ADD_PARAMETER(initial_segment_velocity, FE_PARAM_DOUBLE, "initial_segment_velocity");

END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	AddProperty(&common_properties, "common_properties");
	AddProperty(&matrix_material, "matrix");
	AddProperty(&angio_stress_policy, "angio_stress_policy");

	AddProperty(&pdd_manager, "pdd_manager");
	AddProperty(&psc_manager, "psc_manager");
	AddProperty(&cm_manager, "cm_manager");
	AddProperty(&velocity_manager, "velocity_manager");
	AddProperty(&im_manager, "im_manager");
	AddProperty(&branch_policy, "branch_policy");
	branch_policy.m_brequired = false;
}

FEAngioMaterial::~FEAngioMaterial()
{
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	// Create symmetry vectors


	if(!matrix_material->Init()) return false;

	if(!common_properties->vessel_material->Init()) return false;

	if (FEElasticMaterial::Init() == false) return false;

	//culture must be initialized here  so pangio is defined
	assert(m_pangio);



	// add the user sprouts
	std::vector<int> matls;
	
	FEMesh& mesh = GetFEModel()->GetMesh();
	
	FECoreBase * base = GetParent();
	if(base)
	{
		matls.emplace_back(base->GetID());
	}
	else
	{
		matls.emplace_back(this->GetID_ang());
	}



	return true;
}


void FEAngioMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	// get the material's coordinate system (if defined)
	FECoordSysMap* pmap = GetCoordinateSystemMap();
	//this allows the local coordinates to work correctly
	if (pmap)
	{
		FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);

		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		// compound the local map with the global material axes
		//mat3d Qlocal = pmap->LocalElementCoord(el, n);
		//pt.m_Q = pt.m_Q * Qlocal;

		vessel_elastic.m_Q = pt.m_Q;
		matrix_elastic.m_Q = pt.m_Q;

		FEElasticMaterial* vess_elastic = common_properties->vessel_material->GetElasticMaterial();
		FEElasticMaterial* mat_elastic = matrix_material->GetElasticMaterial();

		vess_elastic->SetLocalCoordinateSystem(el, n, *angioPt->vessPt);
		mat_elastic->SetLocalCoordinateSystem(el, n, *angioPt->matPt);
	}
	
}


//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::AngioStress(FEAngioMaterialPoint& angioPt)
{
	mat3ds val;
	val.zero();
	return val;
}

double FEAngioMaterial::FindDensityScale(FEAngioMaterialPoint * mp)
{
	return 1.0;
}



mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	s.zero();
	//should always be true but we should check
	assert(angioPt);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		vessel_elastic.m_rt = elastic_pt.m_rt;//spatial position
		vessel_elastic.m_r0 = elastic_pt.m_r0;//material position
		vessel_elastic.m_F = elastic_pt.m_F;//deformation gradient
		vessel_elastic.m_J = elastic_pt.m_J;//determinate
		
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;

		vessel_elastic.m_s = common_properties->vessel_material->Stress(*(angioPt->vessPt));
		matrix_elastic.m_s = matrix_material->Stress(*(angioPt->matPt));

		s = angioPt->m_as + angioPt->vessel_weight*vessel_elastic.m_s + angioPt->matrix_weight*matrix_elastic.m_s;
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEAngioMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	tens4ds s(0.0);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		s = angioPt->vessel_weight*common_properties->vessel_material->Tangent(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->Tangent(*angioPt->matPt);
	}
	return s;
}

double FEAngioMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
    
	// calculate strain energy density
	double sed = 0.0;
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		sed = angioPt->vessel_weight*common_properties->vessel_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->matPt);
	}
	return sed;
}

void FEAngioMaterial::SetSeeds(AngioElement* angio_elem)
{
	static int offset = 0;
	static int seed = m_pangio->m_fem->GetGlobalConstant("seed");
	angio_elem->_rengine.seed(offset + seed);
	offset++;
}

double FEAngioMaterial::GetSegmentVelocity(AngioElement * angio_element, vec3d local_pos, FEMesh* mesh)
{
	return velocity_manager->ApplyModifiers(1, local_pos, angio_element, mesh);
}

double FEAngioMaterial::GetInitialVelocity()
{
	return initial_segment_velocity;
}

double FEAngioMaterial::GetMin_dt(AngioElement* angio_elem, FEMesh* mesh)
{
	assert(angio_elem);
	assert(angio_elem->_elem);

	double V = m_pangio->GetMesh()->ElementVolume(*angio_elem->_elem);

	// Calculate the maximum growth velocity based on the potential velocities calculated at each of the Gauss points.
	double max_grow_velocity = 0.0;
	for (int i = 0; i < angio_elem->_elem->GaussPoints();i++)
	{
		vec3d pos(angio_elem->_elem->gr(i), angio_elem->_elem->gs(i), angio_elem->_elem->gt(i));
		double vel = angio_elem->_angio_mat->GetSegmentVelocity(angio_elem, pos, mesh);
		max_grow_velocity = std::max(max_grow_velocity, vel);
	}

	double min_side_length_so_far = std::numeric_limits<double>::max();

	for (int i = 0; i < angio_elem->_elem->GaussPoints(); i++)
	{
		double Gr[FEElement::MAX_NODES];
		double Gs[FEElement::MAX_NODES];
		double Gt[FEElement::MAX_NODES];
		angio_elem->_elem->shape_deriv(Gr, Gs, Gt, angio_elem->_elem->gr(i), angio_elem->_elem->gs(i), angio_elem->_elem->gt(i));

		// Calculate the length of each side of the element.
		vec3d er, es, et; // Basis vectors of the natural coordinates
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			er += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gr[j];
		}
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			es += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gs[j];
		}
		for (int j = 0; j < angio_elem->_elem->Nodes(); j++)
		{
			et += mesh->Node(angio_elem->_elem->m_node[j]).m_rt * Gt[j];
		}

		min_side_length_so_far = std::min(min_side_length_so_far, er.norm2()); // norm2() is magnitude without sqrt.
		min_side_length_so_far = std::min(min_side_length_so_far, es.norm2());
		min_side_length_so_far = std::min(min_side_length_so_far, et.norm2());
	}

	min_side_length_so_far = sqrt(min_side_length_so_far);

	// Basis vectors already represented correctly for HEX8 elements, need to take half of this value for TET4's.
	switch(angio_elem->_elem->Type())
	{
		case FE_Element_Type::FE_TET4G1:
		case FE_Element_Type::FE_TET4G4:
			min_side_length_so_far /= 2;
		default:
			break;
	}

	// Return dt = d / v.
	return (min_side_length_so_far * 0.25) / max_grow_velocity;
}

void FEAngioMaterial::GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index, double min_scale_factor)
{
	assert(angio_elem);
	for(int i=0; i < angio_elem->adjacency_list.size();i++)
	{
		if(angio_elem->adjacency_list.at(i))
		{
			std::vector<Tip*> & tips = angio_elem->adjacency_list.at(i)->active_tips[buffer_index].at(angio_elem);
			for(int j=0; j < tips.size();j++)
			{
				assert(tips.at(j)->angio_element == angio_elem);
				GrowthInElement(end_time, tips.at(j), i, buffer_index, min_scale_factor);
			}
			
		}
		
	}

	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index].at(angio_elem);
	for(int j=0; j < tips.size();j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		GrowthInElement(end_time, tips[j], -1, buffer_index, min_scale_factor);
	}
}

void FEAngioMaterial::GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index, double min_scale_factor)
{
	auto angio_element = active_tip->angio_element;
	
	assert(active_tip);
	auto mesh = m_pangio->GetMesh(); // DEBUG: active_tip->angio_element->_elem->m_nID == 528
	const double min_segm_len = 0.1;
	int next_buffer_index = (buffer_index + 1) % 2;
	double dt = end_time - active_tip->time;
	double grow_len;
	double grow_vel;
	if(active_tip->time < 0.0)
	{
		grow_len = angio_element->_angio_mat->GetInitialVelocity()  *dt;
	}
	else
	{
		grow_vel = this->GetSegmentVelocity(angio_element, active_tip->GetLocalPosition(), mesh);
		grow_len =  grow_vel*dt;
	}

	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];
	vec3d local_pos = active_tip->GetLocalPosition();
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, local_pos.x, local_pos.y, local_pos.z);
	vec3d er, es, et; // Basis vectors of the natural coordinates
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		er += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gr[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		es += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gs[j];
	}
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		et += mesh->Node(angio_element->_elem->m_node[j]).m_rt* Gt[j];
	}
	mat3d natc_to_global(er, es, et);
	mat3d global_to_natc = natc_to_global.inverse();
	vec3d psc_dir = psc_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip, mesh);
	vec3d pdd_dir = pdd_manager->ApplyModifiers(vec3d(1, 0, 0), active_tip, mesh);
	double alpha = cm_manager->ApplyModifiers(0, active_tip, mesh);
	vec3d global_dir = mix(psc_dir, pdd_dir , (alpha * dt));
	global_dir.unit();

	vec3d global_pos;
	double H[FEElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}

	// Just do the transformations
	vec3d nat_dir = global_to_natc * global_dir;
	double factor;
	bool proj_sucess = m_pangio->ScaleFactorToProjectToNaturalCoordinates(active_tip->angio_element->_elem, nat_dir, local_pos, factor, min_scale_factor); // TODO: HERE
	vec3d possible_natc = local_pos + (nat_dir*factor);
	double possible_grow_length = m_pangio->InElementLength(active_tip->angio_element->_elem, local_pos, possible_natc);
	if ((possible_grow_length >= grow_len) && proj_sucess)
	{
		//this is still in the same element
		vec3d real_natc = local_pos + (nat_dir * (factor * (grow_len / possible_grow_length)));
		//add the segment and add the new tip to the current angio element
		Tip * next = new Tip(active_tip, mesh);
		Segment * seg = new Segment();
		next->time = end_time;
		next->angio_element = angio_element;
		next->SetLocalPosition(real_natc);
		next->parent = seg;
		next->face = angio_element;
		next->growth_velocity = grow_vel;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G4 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		assert(next->angio_element->_elem->Type() == FE_Element_Type::FE_TET4G1 ? (local_pos.x + local_pos.y + local_pos.z) < 1.01 : true);
		angio_element->next_tips.at(angio_element).push_back(next);

		//move next to use the rest of the remaining dt

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}


		angio_element->recent_segments.push_back(seg);
		angio_element->final_active_tips.push_back(next);
		assert(next->angio_element);
	}
	else if(proj_sucess)
	{
		//the segment only grows for a portion of dt
		double tip_time_start = active_tip->time + (possible_grow_length / grow_len) * dt;
		assert(tip_time_start < end_time);
		//grow len should always be nonzero so this division should be okay
		//this segment is growing into a new element do the work to do this

		//next is just the end of the segment not the begining of the next segment
		Tip * next = new Tip();
		Segment * seg = new Segment();
		next->time = tip_time_start;

		next->angio_element = angio_element;

		//still need to update the local position of the tip
		next->SetLocalPosition(local_pos + (nat_dir * possible_grow_length));

		next->parent = seg;
		next->face = angio_element;
		next->initial_fragment_id = active_tip->initial_fragment_id;
		next->growth_velocity = grow_vel;

		seg->back = active_tip;
		seg->front = next;
		if (active_tip->parent)
		{
			seg->parent = active_tip->parent;
		}

		vec3d pos = next->GetPosition(mesh);
		angio_element->recent_segments.push_back(seg);

		for(int i=0; i < angio_element->adjacency_list.size();i++) // TODO: HERE
		{
			AngioElement * ang_elem = angio_element->adjacency_list[i];
			if(ang_elem)
			{
				FESolidDomain * sd = dynamic_cast<FESolidDomain*>(ang_elem->_elem->GetDomain());
				if(sd)
				{
					double r[3];//new natural coordinates
					sd->ProjectToElement(*ang_elem->_elem,pos,r);
					if(m_pangio->IsInBounds(ang_elem->_elem, r))
					{
						Tip * adj = new Tip(next, mesh);
						adj->angio_element = angio_element->adjacency_list[i];
						adj->face = angio_element;
						adj->SetLocalPosition(vec3d(r[0], r[1], r[2]));
						adj->use_direction = true;
						adj->growth_velocity = grow_vel;

						angio_element->active_tips[next_buffer_index].at(ang_elem).push_back(adj);

						// break; // Place the tip in exactly one element
					}
				}
			}
		}
	}
	else
	{
		//handle the difficult cases here
		//this occours when the tip is close to a face
		//this also might be when a segment is supposed to die
		assert(false);
	}
}

void FEAngioMaterial::PostGrowthUpdate(AngioElement* angio_elem, double end_time, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{
	//do the branching
	if(angio_elem->_angio_mat->branch_policy)
	{
		angio_elem->_angio_mat->branch_policy->AddBranches(angio_elem, buffer_index, mesh, feangio);
	}
	

	//do the buffer management
	for(auto iter = angio_elem->active_tips[buffer_index].begin(); iter != angio_elem->active_tips[buffer_index].end(); ++iter)
	{
		iter->second.clear();
	}
}

void FEAngioMaterial::Cleanup(AngioElement* angio_elem, double end_time, int buffer_index)
{
	for(int i=0;i < angio_elem->recent_segments.size();i++)
	{
		angio_elem->grown_segments.push_back(angio_elem->recent_segments[i]);
	}
	angio_elem->recent_segments.clear();
	angio_elem->processed_recent_segments = 0;
}

void FEAngioMaterial::PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index)
{
	angio_elem->next_tips.swap(angio_elem->active_tips[buffer_index]);
	for(auto iter=angio_elem->next_tips.begin(); iter != angio_elem->next_tips.end(); ++iter)
	{
		iter->second.clear();
	}
	angio_elem->final_active_tips.clear();
}

bool FEAngioMaterial::SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh)
{
	return common_properties->fseeder->SeedFragments(angio_elements, mesh, this, 0);
}

vec3d FEAngioMaterial::ApplySegmentDirectionModifiers(Tip * tip, double dt)
{
	
}

void FEAngioMaterial::UpdateGDMs()
{
	common_properties->UpdateGDMs();
}


//-----------------------------------------------------------------------------
FEMaterialPoint* FEAngioMaterial::CreateMaterialPointData()
{
	return new FEAngioMaterialPoint(new FEFiberMaterialPoint(FEElasticMaterial::CreateMaterialPointData()), common_properties->vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData());
}