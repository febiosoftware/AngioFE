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
	
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	AddProperty(&common_properties, "common_properties");
	AddProperty(&matrix_material, "matrix");
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


void FEAngioMaterial::UpdateAngioStresses()
{
	
	
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

double FEAngioMaterial::GetGrowthLengthOverUnitTime(AngioElement * angio_element, vec3d local_pos)
{
	return 140.0;
}

double FEAngioMaterial::GetMin_dt(AngioElement* angio_elem)
{
	assert(angio_elem);
	assert(angio_elem->_elem);

	double V = m_pangio->GetMesh()->ElementVolume(*angio_elem->_elem);

	double max_grow_length = 0.0;
	for (int i = 0; i < angio_elem->_elem->GaussPoints();i++)
	{
		vec3d pos(angio_elem->_elem->gr(i), angio_elem->_elem->gs(i), angio_elem->_elem->gt(i));
		double len = angio_elem->_angio_mat->GetGrowthLengthOverUnitTime(angio_elem, pos);
		max_grow_length = std::max(max_grow_length, len);

	}
	//get the smallest of the 12 possible edges
	//lut for the corners
	const vec3d corners[] = { {1,1,1}, {-1,1,1},{-1,-1,1}, {-1,1,-1}, {-1,-1,-1},{1,-1,1},{1,1,-1}, {1,-1,-1} };
	//lut for adjacency
	const std::pair<const size_t, const size_t> adjacency[] = { {0,1},{0,5},{0,6},
	{1,2},{1,3},
	{2,5},{2,4},
	{3,4},{3,6},
	{4,7},
	{5,7},
	{6,7}
	};
	double min_side_length = std::numeric_limits<double>::max();
	for(int i=0; i< sizeof(adjacency)/sizeof(adjacency[0]);i++)
	{
		vec3d p0 = m_pangio->Position(angio_elem->_elem, corners[adjacency[i].first]);
		vec3d p1 = m_pangio->Position(angio_elem->_elem, corners[adjacency[i].second]);
		vec3d res = p0 - p1;
		double len = res.norm();//consider using norm2 is possible
		assert(len != 0.0);
		min_side_length = std::min(min_side_length, len);
	}
	return (min_side_length*0.25)/max_grow_length;
}

void FEAngioMaterial::GrowSegments(AngioElement * angio_elem, double end_time, int buffer_index)
{
	assert(angio_elem);
	for(int i=0; i < angio_elem->adjacency_list.size();i++)
	{
		if(angio_elem->adjacency_list[i])
		{
			std::vector<Tip*> & tips = angio_elem->adjacency_list[i]->active_tips[buffer_index][angio_elem];
			for(int j=0; j < tips.size();j++)
			{
				assert(tips[j]->angio_element == angio_elem);
				GrowthInElement(end_time, tips[j], i, buffer_index);
			}
			
		}
		
	}
	std::vector<Tip*> & tips = angio_elem->active_tips[buffer_index][angio_elem];
	for(int j=0; j < tips.size();j++)
	{
		assert(tips[j]->angio_element == angio_elem);
		GrowthInElement(end_time, tips[j], -1, buffer_index);
	}
}

void FEAngioMaterial::GrowthInElement(double end_time, Tip * active_tip, int source_index, int buffer_index)
{
	auto angio_element = active_tip->angio_element;
	assert(active_tip);
	auto mesh = m_pangio->GetMesh();
	const double eps = 0.01;
	const double min_segm_len = 0.1;
	int next_buffer_index = (buffer_index + 1) % 2;
	double dt = end_time - active_tip->time;
	double grow_len = this->GetGrowthLengthOverUnitTime(angio_element, active_tip->local_pos) *dt;
	

	double Gr[FEElement::MAX_NODES];
	double Gs[FEElement::MAX_NODES];
	double Gt[FEElement::MAX_NODES];
	angio_element->_elem->shape_deriv(Gr, Gs, Gt, active_tip->local_pos.x, active_tip->local_pos.y, active_tip->local_pos.z);
	vec3d er, es, et;//basis vectors of the natural coordinates
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
	vec3d global_dir = active_tip->GetDirection(mesh);
	vec3d global_pos;
	double H[FEElement::MAX_NODES];
	angio_element->_elem->shape_fnc(H, active_tip->local_pos.x, active_tip->local_pos.y, active_tip->local_pos.z);
	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		global_pos += mesh->Node(angio_element->_elem->m_node[j]).m_rt* H[j];
	}
	//exclude the source index face from the raycasts, -1 for the same element
	for(int i=0; i < angio_element->inner_faces.Elements();i++)
	{
		if (i != source_index)
		{
			FESurfaceElement & fse = angio_element->inner_faces.Element(i);
			double len;
			double rs[2];
			
			//TODO: intersect only works for linear elements 
			if (angio_element->inner_faces.Intersect(fse, global_pos ,global_dir, rs, len, eps))
			{
				double tip_time_start = -1;
				//in element growth
				if (len >= grow_len)
				{
					len = grow_len;
					tip_time_start = end_time;

					if (len > min_segm_len)
					{
						vec3d natc_grow_vec = global_to_natc * (global_dir * len);
						//add the segment and add the new tip to the current angio element
						Tip * next = new Tip(active_tip, mesh);
						Segment * seg = new Segment();
						next->time = tip_time_start;
						next->local_pos = active_tip->local_pos + natc_grow_vec;
						next->parent = seg;
						next->face = angio_element;
						next->angio_element = angio_element;
						next->initial_fragment_id = active_tip->initial_fragment_id;
						angio_element->next_tips[angio_element].push_back(next);

						//move next to use the rest of the remaining dt
						
						seg->back = active_tip;
						seg->front = next;
						if(active_tip->parent)
						{
							seg->parent = active_tip->parent;
						}

						
						angio_element->recent_segments.push_back(seg);
						assert(next->angio_element);
						next->PrintTipInfo(mesh, "next tip(no collision with faces)");
					}
				}
				//hitting the face add the segment, then add the tip to the appropiate bucket
				//two tips will be created as the 
				else if(len > min_segm_len)
				{
					tip_time_start = active_tip->time + (len / grow_len) * dt;
					assert(tip_time_start < end_time);
					//grow len should always be nonzero so this division should be okay
					//this segment is growing into a new element do the work to do this

					vec3d natc_grow_vec = global_to_natc * (global_dir * len);
					//add the segment and add the new tip to the current angio element
					Tip * next = new Tip();
					Segment * seg = new Segment();
					next->time = tip_time_start;

					//need to change the angio element and update the local position z
					next->angio_element = angio_element->adjacency_list[i];
					//still need to update the local position of the tip
					next->local_pos = active_tip->local_pos + natc_grow_vec;
					

					next->parent = seg;
					next->face = angio_element;
					next->initial_fragment_id = active_tip->initial_fragment_id;
					//move next to use the rest of the remaining dt
					

					seg->back = active_tip;
					seg->front = next;
					if (active_tip->parent)
					{
						seg->parent = active_tip->parent;
					}
					

					angio_element->recent_segments.push_back(seg);
					

					if(angio_element->adjacency_list[i])
					{
						//map this tip into the adjacent element
						Tip * adj = new Tip(next, mesh);
						adj->angio_element = angio_element->adjacency_list[i];
						adj->face = angio_element;
						vec3d pos = next->GetPosition(mesh);
						FESolidDomain * sd = dynamic_cast<FESolidDomain*>(angio_element->_elem->GetDomain());
						assert(sd);
						double new_rst[3];
						sd->ProjectToElement(*angio_element->adjacency_list[i]->_elem, pos, new_rst);
						adj->local_pos = vec3d(new_rst[0], new_rst[1], new_rst[2]);
						adj->use_direction = true;

						angio_element->active_tips[next_buffer_index][angio_element->adjacency_list[i]].push_back(adj);

						adj->PrintTipInfo(mesh, "adj tip");
					}
					else
					{
						next->angio_element = angio_element;
					}
					next->PrintTipInfo(mesh, "next tip(collision with face)");
					assert(next->angio_element);
				}
				if ((len >= 0) && (len < min_segm_len))
				{
					//assert(false);
					//make a copy of the tip and remap it to all valid elements

				}
				//if len is negative do nothing
			}
		}
	}
}

void FEAngioMaterial::PostGrowthUpdate(AngioElement* angio_elem, double end_time, int buffer_index)
{
	angio_elem->active_tips[buffer_index].clear();
}

void FEAngioMaterial::Cleanup(AngioElement* angio_elem, double end_time, int buffer_index)
{
	for(int i=0;i < angio_elem->recent_segments.size();i++)
	{
		angio_elem->grown_segments.push_back(angio_elem->recent_segments[i]);
	}
	angio_elem->recent_segments.clear();
}

void FEAngioMaterial::PrepBuffers(AngioElement* angio_elem, double end_time, int buffer_index)
{
	angio_elem->next_tips.swap(angio_elem->active_tips[buffer_index]);
	angio_elem->next_tips.clear();
}

bool FEAngioMaterial::SeedFragments(std::vector<AngioElement *>& angio_elements)
{
	return common_properties->fseeder->SeedFragments(angio_elements, this,0);
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