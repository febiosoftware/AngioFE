///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FELoadCurve.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMix/FESolute.h"
#include "FEBioMix/FEMultiphasic.h"
#include "FEBioLib/FEBioModel.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FECore/FESolidDomain.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FECore/FEElemElemList.h"
#include "angio3d.h"
#include <ctime>
#include <future>
#include <algorithm>
#include <unordered_set>
#include <cfloat> 
#include <FEBioMix/FETriphasic.h>
#include <FECore/FENodeDataMap.h>
#include "NodeDataInterpolation.h"
#include <omp.h>
#include <iostream>
#include <math.h>
#include "FECore/FEDomainMap.h"

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat);

// need to rename from vegfID to solute id
bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID);

void FEAngio::SetSeeds()
{
	// gets seed for random engine
	rengine.seed(int (m_fem->GetGlobalConstant("seed")));
	// for each angio element
	for(int i=0; i <angio_elements.size();i++)
	{
		// call set seeds for the element
		angio_elements[i]->_angio_mat->SetSeeds(angio_elements[i]);
	}
}

bool FEAngio::SeedFragments()
{
	bool rv = true;
	// for each angio element across materials
	for(auto iter = elements_by_material.begin(); iter != elements_by_material.end(); ++iter)
	{
		// Call seed fragments for the element
		rv &= iter->first->SeedFragments(iter->second, GetMesh());
	}
	return rv;
}

// ensures that the FE step matches the angio step
void FEAngio::ApplydtToTimestepper(double dt, bool initial)
{
	FEAnalysis * fea = m_fem->GetCurrentStep();
	// if the FE step is bigger than the angio step
	if(fea->m_dt > dt)
	{
		//set the fe step to the angio step
		fea->m_dt = dt;
	}
}

// 
void FEAngio::GrowSegments(double min_scale_factor, double bounds_tolerance, double min_angle, int growth_substeps)
{
	// get the current time
	FETimeInfo time_info = m_fem->GetTime();
	// get the current FE step
	FEAnalysis * cs = m_fem->GetCurrentStep();
	auto   mesh = GetMesh();
	static double min_dt;
	
	size_t angio_element_count = angio_elements.size();
	

	// Is it valid for this to run more than once? only 
	// if the current time is greater than the next time
	if(time_info.currentTime >= next_time)
	{	// 
		min_dt = std::numeric_limits<double>::max();
		// parallel for with value min_dt shared
		#pragma omp parallel for shared(min_dt)
		// for each angio element
		for (int i = 0; i < angio_element_count; ++i)
		{
			// determines the minimum timestep for each angio element to allow growth safely
			double temp_dt = angio_elements[i]->_angio_mat->GetMin_dt(angio_elements[i], mesh);
			// this part must be executed one thread at a time to prevent overwriting
			#pragma omp critical
			{
				// change min_dt to the lesser of min_dt and temp dt
				min_dt = std::min(min_dt, temp_dt);
			}
		}
		// check that the max_angio time is not greater than the angio dt max.
		if (min_dt > max_angio_dt) { min_dt = max_angio_dt; }
		// check that the step is not too small
		if ((min_angio_dt != -1) && (min_dt < min_angio_dt)) 
		{ 
			min_dt = min_angio_dt; 
		}
				
		// current time is next time plus min_dt
		double ctime = next_time + min_dt;
		// allow 16 iterations to be run on each thread as they become available
		#pragma omp parallel for schedule(dynamic, 16)
		// for each angio element
		for (int j = 0; j <angio_element_count; j++)
		{
			// prep the buffers
			angio_elements[j]->_angio_mat->PrepBuffers(angio_elements[j], next_time, buffer_index);
		}
//worse performace than separate declarations
//#pragma omp parallel
		{
			// for each growth substep
			for (int i = 0; i <growth_substeps; i++)
			{
				// for 16 iteration chunks
				#pragma omp parallel for schedule(dynamic, 16)
				// for each angio element
				for (int j = 0; j <angio_element_count; j++)
				{
					// Grow segments
					angio_elements.at(j)->_angio_mat->GrowSegments(angio_elements.at(j), ctime, buffer_index, min_scale_factor,bounds_tolerance, min_angle);
				}

				#pragma omp parallel for schedule(dynamic, 16)
				// for each angio element 
				for (int j = 0; j <angio_element_count; j++)
				{
					// Post growth update for each element
					angio_elements[j]->_angio_mat->PostGrowthUpdate(angio_elements[j], ctime, cs->m_tend , min_scale_factor, buffer_index, mesh, this);
				}
				// update the buffer index
				buffer_index = (buffer_index + 1) % 2;
			}
		}
	}
	ApplydtToTimestepper(min_dt, true);

	//do the output
	fileout->save_vessel_state(*this);
	fileout->save_final_cells_txt(*this);


	//do the cleanup if needed
	if(time_info.currentTime >= next_time)
	{
		next_time += min_dt;
		printf("\nangio dt chosen is: %lg next angio time\n", min_dt);
		#pragma omp parallel for schedule(dynamic, 16)
		for (int j = 0; j <angio_element_count; j++)
		{
			angio_elements[j]->_angio_mat->Cleanup(angio_elements[j], next_time, buffer_index);
		}
		buffer_index = (buffer_index + 1) % 2;
	}

}

void FEAngio::ProtoGrowSegments(double min_scale_factor, double bounds_tolerance, double min_angle, int growth_substeps)
{
	//grow from time -1 to time 0
	auto   time_info = m_fem->GetTime();
	auto   mesh = GetMesh();
	double min_dt = std::numeric_limits<double>::max();
	size_t angio_element_count = angio_elements.size();

	// while mechanical analysis still hasn't begun
	while(next_time < 0)
	{
		// Is it valid for this to run more than once? only 
		if (time_info.currentTime >= next_time)
		{
			#pragma omp parallel for shared(min_dt)
			// for each angio element
			for (int i = 0; i < angio_element_count; ++i)
			{
				// get min dt for all element
				double temp_dt = angio_elements[i]->_angio_mat->GetMin_dt(angio_elements[i], mesh);
				// only one thread can perform this at a time
				#pragma omp critical
				{
					// set min_dt to the lower of min_dt and temp_dt
					min_dt = std::min(min_dt, temp_dt);
				}
			}

			min_dt = 1.0 / (growth_substeps*6.0);

			//min_dt = 0;

			// update current time
			double ctime = next_time + min_dt;
			//clamp the dt to hit zero
			if(ctime > 0)
			{
				min_dt = 0.0 - next_time;
				ctime = 0;
			}

			// prepare buffers for each angio element
			#pragma omp parallel for schedule(dynamic, 32)
			for (int j = 0; j <angio_element_count; j++)
			{
				angio_elements[j]->_angio_mat->PrepBuffers(angio_elements[j], next_time, buffer_index);
			}
			//worse performace than separate declarations
			//#pragma omp parallel
			{
				int n = growth_substeps;
				//				int n = 3;
				for (int i = 0; i <n; i++)
				{
					#pragma omp parallel for schedule(dynamic, 32)
					// for each element
					for (int j = 0; j <angio_element_count; j++)
					{
						// grow the segments
						angio_elements.at(j)->_angio_mat->ProtoGrowSegments(angio_elements.at(j), ctime, buffer_index, min_scale_factor, bounds_tolerance, min_angle);
					}

					#pragma omp parallel for schedule(dynamic, 32)
					// for each angio element
					for (int j = 0; j <angio_element_count; j++)
					{
						// update the elements
						angio_elements[j]->_angio_mat->ProtoPostGrowthUpdate(angio_elements[j], ctime, min_scale_factor, buffer_index, mesh, this);
					}
					buffer_index = (buffer_index + 1) % 2;
				}
			}

		}

		//do the cleanup if needed
		if (time_info.currentTime >= next_time)
		{
			next_time += min_dt;
			//printf("\nproto angio dt chosen is: %lg next angio time\n", min_dt);
			#pragma omp parallel for schedule(dynamic, 16)
			for (int j = 0; j <angio_element_count; j++)
			{
				angio_elements[j]->_angio_mat->Cleanup(angio_elements[j], next_time, buffer_index);
			}
			buffer_index = (buffer_index + 1) % 2;
		}
	}
	//this may be a bit conservative for a first step but should produce good results
	ApplydtToTimestepper(min_dt);
	//do the output, this will output all of the segments
	fileout->bulk_save_vessel_state(*this);
	fileout->save_initial_cells_txt(*this);
}

// takes a given element then first counts the tips in that element then adds all adjacent elements to an inspection list and keep adding tips and adjacent elements until
// they are out of bounds. 
void FEAngio::GetActiveFinalTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip *> & tips)
{
	std::unordered_set<AngioElement *> visited;
	std::set<AngioElement *> next;
	std::vector<vec3d> element_bounds;
	visited.reserve(1000);
	tips.reserve(500);
	pangio->ExtremaInElement(angio_element->_elem, element_bounds);
	// add element to be inspected
	next.insert(angio_element);
	// while there are still elements in next
	while (next.size())
	{
		// get angio element
		AngioElement * cur = *next.begin();
		// remove it from the buffer
		next.erase(next.begin());
		// add it to the visited buffer
		visited.insert(cur);
		std::vector<vec3d> cur_element_bounds;
		// Get the extrema for the current element
		pangio->ExtremaInElement(cur->_elem, cur_element_bounds);
		// get the minimum distance from the bounds of the element that we're calculating tips in the radius of to the iteration's element
		double cdist = FEAngio::MinDistance(element_bounds, cur_element_bounds);
		// if the distance is in the desired range
		if (cdist <= radius)
		{
			//add the tips and the add all unvisited adjacent elements to next
			for (int i = 0; i < cur->final_active_tips.size(); i++)
			{
				tips.push_back(cur->final_active_tips[i]);
			}

			for (int i = 0; i < cur->face_adjacency_list.size(); i++)
			{
				if (!visited.count(cur->face_adjacency_list[i]))
				{
					next.insert(cur->face_adjacency_list[i]);
				}
			}
		}
	}
}

void FEAngio::GetActiveTipsInRadius(AngioElement* angio_element, double radius, int buffer, FEAngio* pangio, std::vector<Tip *> & tips, int excliude)
{
	std::unordered_set<AngioElement *> visited;
	std::set<AngioElement *> next;
	std::vector<vec3d> element_bounds;
	visited.reserve(1000);
	tips.reserve(500);
	pangio->ExtremaInElement(angio_element->_elem, element_bounds);

	for (int i = 0; i < angio_element->face_adjacency_list.size(); i++)
	{
		next.insert(angio_element->face_adjacency_list[i]);
	}
	visited.insert(angio_element);
	while (next.size())
	{
		AngioElement * cur = *next.begin();
		next.erase(next.begin());
		visited.insert(cur);
		std::vector<vec3d> cur_element_bounds;
		pangio->ExtremaInElement(cur->_elem, cur_element_bounds);
		double cdist = FEAngio::MinDistance(element_bounds, cur_element_bounds); 
		//std::cout << cdist << endl;
		if (cdist <= radius)
		{
			for(auto iter = cur->active_tips[buffer].begin(); iter != cur->active_tips[buffer].end(); ++iter)
			{
				//add the tips and the add all unvisited adjacent elements to next
				for (int i = 0; i < iter->second.size(); i++)
				{
					Tip * tip = iter->second[i];
					if(tip->initial_fragment_id != excliude)
					{
						tips.push_back(tip);
					}
				}
			}
			for (int i = 0; i < cur->face_adjacency_list.size(); i++)
			{
				if (!visited.count(cur->face_adjacency_list[i]))
				{
					next.insert(cur->face_adjacency_list[i]);
				}
			}
		}
	}
}

void FEAngio::GetGrownTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip *> & tips)
{
	std::unordered_set<AngioElement *> visited;
	std::set<AngioElement *> next;
	std::vector<vec3d> element_bounds;
	visited.reserve(1000);
	tips.reserve(500);
	pangio->ExtremaInElement(angio_element->_elem, element_bounds);

	for (int i = 0; i < angio_element->face_adjacency_list.size(); i++)
	{
		next.insert(angio_element->face_adjacency_list[i]);
	}
	visited.insert(angio_element);
	while (next.size())
	{
		AngioElement * cur = *next.begin();
		next.erase(next.begin());
		visited.insert(cur);
		std::vector<vec3d> cur_element_bounds;
		pangio->ExtremaInElement(cur->_elem, cur_element_bounds);
		double cdist = FEAngio::MinDistance(element_bounds, cur_element_bounds);
		if (cdist <= radius)
		{
			for (auto iter = cur->grown_segments.begin(); iter != cur->grown_segments.end(); ++iter)
			{
				
				tips.push_back((*iter)->front);
				tips.push_back((*iter)->back);
			}

			for (int i = 0; i < cur->face_adjacency_list.size(); i++)
			{
				if (!visited.count(cur->face_adjacency_list[i]))
				{
					next.insert(cur->face_adjacency_list[i]);
				}
			}
		}
	}
}


vec3d FEAngio::ReferenceCoordinates(Tip * tip) const
{
	vec3d r(0, 0, 0);
	FEMesh & mesh = m_fem->GetMesh();
	//Point has already been positioned
	FESolidElement * se= tip->angio_element->_elem;

	double arr[FESolidElement::MAX_NODES];
	// get tip local coordinates
	vec3d local_pos = tip->GetLocalPosition();
	se->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	// for each node in the element
	for (int j = 0; j < se->Nodes(); j++)
	{
		// add the initial global position
		r += mesh.Node(se->m_node[j]).m_r0* arr[j];
	}

	return r;
}

double FEAngio::NaturalCoordinatesUpperBound_r(int et)
{
	switch(et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return 1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 1.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return 1.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

double FEAngio::NaturalCoordinatesUpperBound_s(int et)
{
	switch (et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return 1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 1.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return 1.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

double FEAngio::NaturalCoordinatesUpperBound_t(int et)
{
	switch (et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return 1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 1.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return 1.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

double FEAngio::NaturalCoordinatesLowerBound_r(int et)
{
	switch (et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return -1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 0.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return 0.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

double FEAngio::NaturalCoordinatesLowerBound_s(int et)
{
	switch (et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return -1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 0.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return 0.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

double FEAngio::NaturalCoordinatesLowerBound_t(int et)
{
	switch (et)
	{
	case FE_HEX8G1:
	case FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
	case FE_HEX8RI:
		return -1.0;
	case FE_TET4G1:
	case FE_TET4G4:
	case FE_TET10G1:
	case FE_TET10G4:
	case FE_TET10G4RI1:
	case FE_TET10G8:
	case FE_TET10G8RI4:
	case FE_TET10GL11:
	case FE_TET15G4:
	case FE_TET15G8:
	case FE_TET15G11:
	case FE_TET15G15:
	case FE_TET15G15RI4:
	case FE_TET20G15:
		return 0.0;
	case FE_PENTA15G21:
	case FE_PENTA15G8:
	case FE_PENTA6G6:
		return -1.0;
	default:
		assert(false);
	}
	return std::numeric_limits<double>::max();
}

// clamping checks if an element is in the prescribed range and if not moves it to the closest acceptable value in the range
// This clamps a position (natc) to the range given by the natural coordinates for an element type et
vec3d FEAngio::clamp_natc(int et, vec3d natc)
{
	return vec3d(std::max(std::min(NaturalCoordinatesUpperBound_r(et),natc.x ), NaturalCoordinatesLowerBound_r(et)),
		std::max(std::min(NaturalCoordinatesUpperBound_s(et), natc.y), NaturalCoordinatesLowerBound_s(et)),
		std::max(std::min(NaturalCoordinatesUpperBound_t(et), natc.z), NaturalCoordinatesLowerBound_t(et))
		);
}

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : ztopi(std::uniform_real_distribution<double>(0, PI)), 
	zto2pi(std::uniform_real_distribution<double>(0, 2 * PI)), n1to1(std::uniform_real_distribution<double>(-1, 1))
{	
	// Input random seed number
	m_irseed = 0;

	m_fem = dynamic_cast<FEBioModel *>(&fem);
	assert(m_fem);
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
}

//-----------------------------------------------------------------------------
FEBioModel* FEAngio::GetFEModel() const
{
	return m_fem;
}

FEMesh * FEAngio::GetMesh() const
{
	return &m_fem->GetMesh();
}

// consider moving this to FEBio
// this returns the range of elements (why is this necessary?)
std::pair<int, int> FEAngio::GetMinMaxElementIDs() const
{
	FEMesh * mesh = GetMesh();
	assert(mesh);

	// get the ID ranges
	int m_minID = -1;
	int m_maxID = -1;
	int NDOM = mesh->Domains();
	// for each material domain
	for (int i = 0; i<NDOM; ++i)
	{
		// get the mesh of the domain
		FEDomain& dom = mesh->Domain(i);
		// get the number of elements in the domain
		int NE = dom.Elements();
		// for each element
		for (int j = 0; j<NE; ++j)
		{
			//  get the element number
			FEElement& el = dom.ElementRef(j);
			// get the element id
			int eid = el.GetID();

			if ((eid < m_minID) || (m_minID == -1)) m_minID = eid;
			if ((eid > m_maxID) || (m_maxID == -1)) m_maxID = eid;
		}
	}
	return { m_minID, m_maxID };
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	//create any classes which have nontrivial destructors 

	// Init all the FE stuff
	//must be done first initializes material
	if (InitFEM() == false)
	{
		return false;
	}
	if (m_fem->Init() == false) return false;

	// start timer
	time(&m_start);

	return true;
}


//-----------------------------------------------------------------------------
// Initialize FE model.
bool FEAngio::InitFEM()
{
	// for each material
	for (int i = 0; i < m_fem->Materials(); i++)
	{
		FEAngioMaterial * cmat = nullptr;
		// get the material id
		FEMaterial * mat = m_fem->GetMaterial(i);
		int id = -1;
		// get the angio component of the material
		cmat = GetAngioComponent(mat);
		
		// if an angio material was successfully found
		if (cmat)
		{
			id = mat->GetID();
			assert(id != -1);
			// add the current angio material to the material pointer
			m_pmat.emplace_back(cmat);
			// add the current angio material's id to the ids
			m_pmat_ids.emplace_back(id);
			//TODO: check that material parameters are set here
			//cmat->ApplySym();
			cmat->SetFEAngio(this);
		}
	}
	//assert(m_pmat.size());
	feLog("%d Angio materials found. Stress approach will be used.", m_pmat.size());

	// register the angio callback

	m_fem->AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED | CB_STEP_ACTIVE | CB_INIT, this, CallbackHandler::CB_ADD_FRONT);

	return true;
}
//fills in the adjacnecy information 
void FEAngio::FillInAdjacencyInfo(FEMesh * mesh,FEElemElemList * eel, AngioElement *angio_element, int elem_index)
{
	//one element of adjaceny per face
	FESolidElement * se = angio_element->_elem;
	assert(se);
	// for each face in the solid element
	for(int i=0; i < se->Faces();i++)
	{
		// get the neighboring element
		FEElement * elem = eel->Neighbor(elem_index, i);
		FESolidElement * nse = dynamic_cast<FESolidElement*>(elem);
		if(nse)
		{
			//next look to see if there is an angio element that uses this element
			auto ae_iter = se_to_angio_elem.find(nse);
			if(ae_iter != se_to_angio_elem.end())
			{
				angio_element->face_adjacency_list.push_back(ae_iter->second.first);
			}
			else
			{
				assert(false);
			}
		}
	}
}

void FEAngio::FillInFaces(FEMesh * mesh, AngioElement * angio_element)
{
	angio_element->adjacency_list = angio_elements_to_all_adjacent_elements[angio_element];
}

void FEAngio::SetupAngioElements()
{
	//setup the element access data structure
	auto min_max = GetMinMaxElementIDs();
	auto min_elementID = std::get<0>(min_max);
	//initiialize the FEAngioElementData Structures

	FEMesh * mesh = GetMesh();
	FEModel * model = GetFEModel();
	//this stores the element adjacency information
	FEElemElemList eel;
	eel.Create(mesh);

	SetupNodesToElement(min_elementID);

	//we are using the lut in FEMesh so the speed should be good
	//the table now contains only angio elements
	// for each element 
	for (int i = 0; i < mesh->Elements(); i++)
	{
		// get the solid element
		FESolidElement * elem = dynamic_cast<FESolidElement*>(mesh->FindElementFromID(min_elementID + i));
		if (elem)
		{
			// get the material and angio part of the material
			FEMaterial * mat = model->GetMaterial(elem->GetMatID());
			FEAngioMaterial*angio_mat = GetAngioComponent(mat);
			AngioElement * angio_element = nullptr;
			if (angio_mat)
			{
				// create a new angio element and add it to the container
				angio_element = new AngioElement(elem, angio_mat, mat, mesh);
				angio_elements.push_back(angio_element);
				//if this is the last element create a pointer for the angio element vector and an array for the elements by material
				if (elements_by_material.find(angio_mat) == elements_by_material.end())
				{
					std::vector<AngioElement *> ang_elem;
					elements_by_material[angio_mat] = ang_elem;
				}
				// if this is the last angio material append it to the angio material container
				if(std::find(angio_materials.begin(), angio_materials.end() ,angio_mat) == angio_materials.end())
				{
					angio_materials.push_back(angio_mat);
				}
				// add the element to the elements by material
				elements_by_material[angio_mat].push_back(angio_element);
				// copy the solid element information
				se_to_angio_elem[elem] = { angio_element,i };
			}
			//now fill in the version with holes
			angio_elements_with_holes.push_back(angio_element);
		}
	}
	// for each angio element fill in the adjacency info
	for(int i=0; i < angio_elements.size();i++)
	{
		FillInAdjacencyInfo(mesh, &eel,  angio_elements[i], se_to_angio_elem[angio_elements[i]->_elem].second);
	}
	// for each angio element copy the solid element information
	for(int i=0; i < angio_elements.size();i++)
	{
		se_to_angio_element[angio_elements[i]->_elem] = angio_elements[i];
	}
	// for each angio element copy the nodal solid element data.
	for(int i=0; i < angio_elements.size();i++)
	{
		for(int j=0; j < angio_elements[i]->_elem->Nodes();j++)
		{
			FENode & node = mesh->Node(angio_elements[i]->_elem->m_node[j]);
			std::vector<FESolidElement *> & adj_elements = nodes_to_elements[&node];
			for(int k=0;k< adj_elements.size();k++)
			{
				auto iter = se_to_angio_element.find(adj_elements[k]);
				if(iter != se_to_angio_element.end() && (iter->second != angio_elements[i]))
				{
					std::vector<AngioElement*> & ang_elems = angio_elements_to_all_adjacent_elements[angio_elements[i]];
					if(std::find(ang_elems.begin(), ang_elems.end(), iter->second) == ang_elems.end())
					{
						ang_elems.push_back(iter->second);
					}
				}	
			}
		}
	}
	// for each angio element fill in face info
	for (int i = 0; i < angio_elements.size(); i++)
	{
		FillInFaces(mesh, angio_elements[i]);
	}
	//make sure that no maps are inserted into during multithreaded code
	// for each angio element fill in the tip information
	for (int i = 0; i < angio_elements.size(); i++)
	{
		std::vector<AngioElement*> & adj_list = angio_elements_to_all_adjacent_elements[angio_elements[i]];
		std::vector<Tip*> tips;
		for(int j=0; j < adj_list.size();j++)
		{
			
			angio_elements[i]->active_tips[0][adj_list[j]] = tips;
			angio_elements[i]->active_tips[1][adj_list[j]] = tips;
			angio_elements[i]->next_tips[adj_list[j]] = tips;
		}
		angio_elements[i]->active_tips[0][angio_elements[i]] = tips;
		angio_elements[i]->active_tips[1][angio_elements[i]] = tips;
		angio_elements[i]->next_tips[angio_elements[i]] = tips;
	}
}


double FEAngio::GetDoubleFromDataStore(int record, int elem_id, int item)
{
	DataStore & ds = m_fem->GetDataStore();
	return ds.GetDataRecord(record)->Evaluate(elem_id, item);
}
#ifndef __linux__
mat3d FEAngio::unifromRandomRotationMatrix(angiofe_random_engine & rengine) const
{
	//collagen fibers are right handed so the following transformation is legal
	//the following will only produce right handed bases for the collagen fibers which is molecularly accurate
	double alpha = zto2pi(rengine);
	double beta = zto2pi(rengine);
	double gamma = zto2pi(rengine);


	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, - s_alpha*s_gamma+c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
		);
	//rv = rv.transinv();
#ifndef NDEBUG
	//do some testing that bases run tthrough this are still orthogonal
	vec3d x(1, 0, 0);
	vec3d y(0, 1, 0);
	vec3d z(0, 0, 1);
	vec3d xt = rv * x;
	vec3d yt = rv * y;
	vec3d zt = rv * z;
	double tol = 0.01;
	assert(xt * yt < tol && xt * yt > -tol);
	assert(xt * zt < tol && xt * zt > -tol);
	assert(zt * yt < tol && zt * yt > -tol);

	double noq = xt.norm();
	double nom1 = yt.norm();
	double nom2 = zt.norm();
	tol = 0.01;
	assert(noq < (1 + tol) && noq >(1 - tol));
	assert(nom1 < (1 + tol) && nom1 >(1 - tol));
	assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif
	return rv;
}
#else
mat3d FEAngio::unifromRandomRotationMatrix(angiofe_random_engine & rengine)
{
	//collagen fibers are right handed so the following transformation is legal
	//the following will only produce right handed bases for the collagen fibers which is molecularly accurate
	double alpha = zto2pi(rengine);
	double beta = zto2pi(rengine);
	double gamma = zto2pi(rengine);

	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, -s_alpha*s_gamma + c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
	);
	//rv = rv.transinv();
#ifndef NDEBUG
	//do some testing that bases run tthrough this are still orthogonal
	vec3d x(1, 0, 0);
	vec3d y(0, 1, 0);
	vec3d z(0, 0, 1);
	vec3d xt = rv * x;
	vec3d yt = rv * y;
	vec3d zt = rv * z;
	double tol = 0.01;
	assert(xt * yt < tol && xt * yt > -tol);
	assert(xt * zt < tol && xt * zt > -tol);
	assert(zt * yt < tol && zt * yt > -tol);

	double noq = xt.rnorm();
	double nom1 = yt.rnorm();
	double nom2 = zt.rnorm();
	tol = 0.01;
	assert(noq < (1 + tol) && noq >(1 - tol));
	assert(nom1 < (1 + tol) && nom1 >(1 - tol));
	assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif
	return rv;
}
#endif

mat3d FEAngio::rotationMatrix(double alpha, double beta, double gamma) const
{
	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, -s_alpha*s_gamma + c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
	);
	return rv;
}

vec3d FEAngio::uniformRandomDirection(angiofe_random_engine& rengine)
{
	//to revert this set this to return vrand
	double theta = zto2pi(rengine);
	double phi = zto2pi(rengine);
	double sintheta = sin(theta);
	vec3d dir(sintheta * cos(phi), sintheta * sin(phi), cos(theta));
	return dir;
}
vec3d FEAngio::uniformInUnitCube()
{
	//force the order of initialization should allow results on different platforms to converge 
	double x = n1to1(rengine);
	double y = n1to1(rengine);
	double z = n1to1(rengine);
	return vec3d(x,y,z);
}

vec3d FEAngio::gradient(FESolidElement * se, std::vector<double> & fn, vec3d pt)
{
	assert(se);
	FESolidDomain* domain = dynamic_cast<FESolidDomain*>(se->GetMeshPartition());
	double Ji[3][3];
	domain->invjact(*se, Ji, pt.x, pt.y, pt.z);
	double Gr[FESolidElement::MAX_NODES], Gs[FESolidElement::MAX_NODES], Gt[FESolidElement::MAX_NODES];
	se->shape_deriv(Gr, Gs, Gt, pt.x, pt.y, pt.z);
	double out_values[FESolidElement::MAX_NODES];
	se->project_to_nodes(&fn[0], out_values);


	vec3d gradf;
	size_t N = se->Nodes();
	for (size_t i = 0; i<N; ++i)
	{
		double Gx, Gy, Gz;
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
		Gy = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
		Gz = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];

		// calculate pressure gradient
		gradf.x += Gx*out_values[i];
		gradf.y += Gy*out_values[i];
		gradf.z += Gz*out_values[i];
	}

	return gradf;
}
//-----------------------------------------------------------------------------
static void solve_3x3(double A[3][3], double b[3], double x[3])
{
	// get the determinant
	double D = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[2][1] * A[0][2] \
		- A[1][1] * A[2][0] * A[0][2] - A[2][2] * A[1][0] * A[0][1] - A[0][0] * A[2][1] * A[1][2];

	assert(D != 0);
	// get the minors of A
	double Ai[3][3];
	Ai[0][0] = A[1][1] * A[2][2] - A[2][1] * A[1][2];
	Ai[0][1] = A[2][1] * A[0][2] - A[0][1] * A[2][2];
	Ai[0][2] = A[0][1] * A[1][2] - A[1][1] * A[0][2];

	Ai[1][0] = A[2][0] * A[1][2] - A[1][0] * A[2][2];
	Ai[1][1] = A[0][0] * A[2][2] - A[2][0] * A[0][2];
	Ai[1][2] = A[1][0] * A[0][2] - A[0][0] * A[1][2];

	Ai[2][0] = A[1][0] * A[2][1] - A[2][0] * A[1][1];
	Ai[2][1] = A[2][0] * A[0][1] - A[0][0] * A[2][1];
	Ai[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	// solve Ax=b for x as x = A^-1 b
	x[0] = (Ai[0][0] * b[0] + Ai[0][1] * b[1] + Ai[0][2] * b[2]) / D;
	x[1] = (Ai[1][0] * b[0] + Ai[1][1] * b[1] + Ai[1][2] * b[2]) / D;
	x[2] = (Ai[2][0] * b[0] + Ai[2][1] * b[1] + Ai[2][2] * b[2]) / D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]);
	r[1] = b[1] - (A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]);
	r[2] = b[2] - (A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]);

	double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
#endif
}

FEAngioMaterial * FEAngio::GetAngioComponent(FEMaterial * mat)
{
	// try to cast an angio material
	auto angm = dynamic_cast<FEAngioMaterial*>(mat);
	// if it is not an angio material check if it is a multiphasic material
	if (!angm)
	{
		// try to cast a multiphasic material
		FEMultiphasic * mmat = dynamic_cast<FEMultiphasic*>(mat);
		// if it is a multiphasic material then get the solid component of the multiphasic material and cast it to the angio material
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	// if it wasn't an angio material or multiphasic material see if it was triphasic
	if (!angm)
	{
		FETriphasic * mmat = dynamic_cast<FETriphasic*>(mat);
		// if a triphasic material is found cast the solid to an angio material.
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	return angm;
}

vec3d FEAngio::Position(FESolidElement * se, vec3d local) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(se);
	se->shape_fnc(arr, local.x, local.y, local.z);
	vec3d rc(0, 0, 0);

	auto mesh = GetMesh();
	for (int j = 0; j < se->Nodes(); j++)
	{
		rc += mesh->Node(se->m_node[j]).m_rt* arr[j];
	}
	return rc;
}

void FEAngio::CalculateSegmentLengths(FEMesh* mesh)
{
	int angio_elements_size = int (angio_elements.size());
#pragma omp parallel for schedule(dynamic, 16)
	for(int i=0; i < angio_elements_size;i++)
	{
		AngioElement * angio_element = angio_elements[i];
		angio_element->global_segment_length = 0.0;//reset this
		for(int j=0; j < angio_element->grown_segments.size();j++)
		{
			angio_element->global_segment_length += angio_element->grown_segments[j]->Length(mesh);
		}
	}
}

void FEAngio::AdjustMatrixVesselWeights(class FEMesh* mesh)
{
	std::cout << "Adjusting weights" << endl;
	double time = GetFEModel()->GetTime().currentTime;
	int angio_elements_size = int (angio_elements.size());
#pragma omp parallel for schedule(dynamic, 16)
	for (int i = 0; i < angio_elements_size; i++)
	{
		AngioElement * angio_element = angio_elements[i];
		double vessel_volume = PI * angio_element->_angio_mat->vessel_radius *angio_element->_angio_mat->vessel_radius * angio_element->refernce_frame_segment_length;
		double element_volume = mesh->CurrentElementVolume(*angio_element->_elem);
		/*double vessel_weight = vessel_volume/element_volume;*/
		angio_element->vessel_weight = vessel_volume / element_volume;
		double matrix_weight = 1 - angio_element->vessel_weight;
		// get vascular density in mm/mm3
		double vascular_density = (angio_element->refernce_frame_segment_length * 1e-3) / (element_volume * 1e-9);
		double FA = 0;
		for (int j = 0; j < angio_element->_elem->GaussPoints(); j++)
		{
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(angio_element->_elem->GetMaterialPoint(j));
			FA += angioPt->angioFA;
		}
		FA = FA / angio_element->_elem->GaussPoints();
		if ((angio_element->vessel_weight*(1.0) > angio_element->_angio_mat->thresh_vess_weight) & (angio_element->vasc_thresh_time < 0.0)) {
			angio_element->vasc_thresh_time = time;
		}

		for(int j=0; j < angio_element->_elem->GaussPoints();j++)
		{
			FEAngioMaterialPoint * angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(angio_element->_elem->GetMaterialPoint(j));
			angioPt->vessel_weight = angio_element->vessel_weight;
			angioPt->matrix_weight = matrix_weight;
			angioPt->vascular_density = vascular_density;
		}
	}
}



//-----------------------------------------------------------------------------
bool FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;
	FEMesh * mesh = GetMesh();
	if(nwhen == CB_INIT)
	{
		SetupAngioElements();

		SetSeeds();

		// only output to the logfile (not to the screen)
		//felog.SetMode(Logfile::LOG_FILE);

		//currently the destructors are not called for classes created by FEBio this allows destructors to be called
		fileout = new Fileout(*this);

		//do the seeding of the tips/segments
		if (!SeedFragments())
		{
			printf("fragment seeding failed\n");
			throw;
		}
		// --- Output initial state of model ---
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			// Output initial microvessel state

			// save active tips
			fileout->save_active_tips(*this);
		}

		//start the protogrowth phase
		min_scale_factor = m_fem->GetGlobalConstant("min_scale_factor");
		bounds_tolerance = m_fem->GetGlobalConstant("bounds_tolerance");
		min_angle = m_fem->GetGlobalConstant("min_angle"); min_angle = cos(((PI/180)*min_angle));
		max_angio_dt = m_fem->GetGlobalConstant("max_angio_dt"); if (max_angio_dt == 0) { max_angio_dt = 0.25; }
		min_angio_dt = m_fem->GetGlobalConstant("min_angio_dt"); 
		growth_substeps = int (m_fem->GetGlobalConstant("growth_substeps"));
		bounce = m_fem->GetGlobalConstant("bounce"); 
		size_t angio_element_count = angio_elements.size();
		FEMesh * mesh = GetMesh();
		// Do angio material initialization for the mat_axis
		// Need to loop over all angio domains assigned
		for (int i = 0; i < mesh->Domains(); i++)
		{
			FEDomain &test_dom = mesh->Domain(i);
			// see if it is an angio domain
			FEMaterial* mat = test_dom.GetMaterial();
			FEAngioMaterial* test_angmat = dynamic_cast<FEAngioMaterial*>(mat->ExtractProperty<FEElasticMaterial>());
			if (test_angmat)
			{
				FEElementSet* elset = mesh->FindElementSet(test_dom.GetName());
				FEDomainMap* map = new FEDomainMap(FE_MAT3D, FMT_MATPOINTS);
				map->Create(elset);
				map->fillValue(mat3d::identity());
				FEParam* matax = test_angmat->FindParameter("mat_axis");
				// create parameter
				FEParamMat3d& p = matax->value<FEParamMat3d>();
				// create evaluator
				FEMappedValueMat3d* val = fecore_alloc(FEMappedValueMat3d, GetFEModel());
				val->setDataMap(map);
				//set the valuator to the model parameter
				p.setValuator(val);
			}
		}
#pragma omp parallel for schedule(dynamic, 16)
		// apply the initial modifiers to each angio element
		for (int i = 0; i < angio_element_count; i++)
		{
			if(angio_elements[i]->_angio_mat->im_manager)
			{
				angio_elements[i]->_angio_mat->im_manager->ApplyModifier(angio_elements[i], mesh, this);
			}
		}
		//now override any specific information that needs it
		// for each element by material do the node interpolations
		for(auto iter = elements_by_material.begin(); iter != elements_by_material.end(); ++iter)
		{
			if(iter->first->nodedata_interpolation_manager)
			{
				iter->first->nodedata_interpolation_manager->DoInterpolations(this, mesh, iter->first);
			}
		}

		//setup the branch info before doing any growth
#pragma omp parallel for schedule(dynamic, 16)
		for (int i = 0; i < angio_element_count; i++)
		{
			if (angio_elements[i]->_angio_mat->proto_branch_policy)
			{
				angio_elements[i]->_angio_mat->proto_branch_policy->SetupBranchInfo(angio_elements[i]);
			}
			else if (angio_elements[i]->_angio_mat->branch_policy)
			{
				angio_elements[i]->_angio_mat->branch_policy->SetupBranchInfo(angio_elements[i]);
			}
		}
		// start the growth timer
		grow_timer.start();
		// grow segments
		ProtoGrowSegments(min_scale_factor, bounds_tolerance, min_angle, growth_substeps);
		grow_timer.stop();

		CalculateSegmentLengths(mesh);
		AdjustMatrixVesselWeights(mesh);

		update_angio_stress_timer.start();
		for (int i = 0; i < angio_materials.size(); i++)
		{
			angio_materials[i]->angio_stress_policy->UpdateScale();
		}
#pragma omp parallel for schedule(dynamic, 16)
		for (int i = 0; i < angio_element_count; i++)
		{
			angio_elements[i]->_angio_mat->angio_stress_policy->AngioStress(angio_elements[i], this, mesh);
		}
		update_angio_stress_timer.stop();
	}

	if (nwhen == CB_UPDATE_TIME)
	{
		static int index = 0;
		// grab the time information
		
		auto & ti = fem.GetTime();

		update_branch_policy_timestep_timer.start();
		for (size_t i = 0; i < m_pmat.size(); i++)
		{	
			if(m_pmat[i]->branch_policy)
			{
				m_pmat[i]->branch_policy->TimeStepUpdate(ti.currentTime);
			}
		}
		update_branch_policy_timestep_timer.stop();

		// update velocity scales
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			for (int j = 0; j < m_pmat[i]->velocity_manager->seg_vel_modifiers.size(); j++)
			{
				m_pmat[i]->velocity_manager->seg_vel_modifiers[j]->UpdateScale();
			}
		}

		//new function to find the start time grow time and if this is the final iteration this timestep
		grow_timer.start();
		GrowSegments(min_scale_factor, bounds_tolerance,min_angle,growth_substeps);
		grow_timer.stop();

		update_angio_stress_timer.start();
		for (int i = 0; i < angio_materials.size(); i++)
		{
			angio_materials[i]->angio_stress_policy->UpdateScale();
		}

		AdjustMatrixVesselWeights(mesh);
		CalculateSegmentLengths(mesh);

		size_t angio_element_count = angio_elements.size();
#pragma omp parallel for schedule(static, 16)
		for (int i = 0; i < angio_element_count; i++)
		{
			angio_elements[i]->_angio_mat->angio_stress_policy->AngioStress(angio_elements[i], this,mesh);
		}
		update_angio_stress_timer.stop();
	}
	else if (nwhen == CB_MAJOR_ITERS)
	{

		fileout->save_feangio_stats(*this);
		ResetTimers();
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			// save active tips
			CalculateSegmentLengths(mesh);
			// Print the status of angio3d to the user    
			fileout->printStatus(*this, fem.GetTime().currentTime);
		}
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			Output();
		}
		//force any destructors to be called that need it
		delete fileout;
		fileout = nullptr;
	}
	return true;
}

void FEAngio::ResetTimers()
{
	grow_timer.reset();
	update_branch_policy_timestep_timer.reset();
	update_angio_stress_timer.reset();
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{	
	//write out the timeline of branchpoints
	fileout->save_timeline(*this);
	//fileout->save_winfiber(*this);
	fileout->save_final_vessel_csv(*this);
}

//returns the scale factor needed to scale the ray to grow to the boundary of the natural coordinates
bool FEAngio::ScaleFactorToProjectToNaturalCoordinates(FESolidElement* se, vec3d & dir, vec3d & pt, double & sf, double min_sf) const
{
	// TODO: Expose min_sf to the user.

	std::vector<double> possible_values;

	// Determine the natural coordinate bounds based on the element type.
	vec3d upper_bounds(FEAngio::NaturalCoordinatesUpperBound_r(se->Type()), FEAngio::NaturalCoordinatesUpperBound_s(se->Type()), FEAngio::NaturalCoordinatesUpperBound_t(se->Type()));
	vec3d lower_bounds(FEAngio::NaturalCoordinatesLowerBound_r(se->Type()), FEAngio::NaturalCoordinatesLowerBound_s(se->Type()), FEAngio::NaturalCoordinatesLowerBound_t(se->Type()));

	// Determine the limiting direction from the current position to a face based on the element type. Store it in the memory allocated as the factor/sf.
	switch (se->Type())
	{
	case FE_Element_Type::FE_HEX8G1:
	case FE_Element_Type::FE_HEX8G8:
	case FE_Element_Type::FE_HEX8RI:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
		if (dir.x > 0)
		{
			double temp = (upper_bounds.x - pt.x) / dir.x;
			if (temp >= 0)
				possible_values.push_back(temp);
		}
		else if (dir.x < 0)
		{
			double temp = (lower_bounds.x - pt.x) / dir.x;
			if (temp >= 0)
				possible_values.push_back(temp);
		}
		if (dir.y > 0)
		{
			double temp = (upper_bounds.y - pt.y) / dir.y;
			if (temp >= 0)
				possible_values.push_back(temp);
		}
		else if (dir.y < 0)
		{
			double temp = (lower_bounds.y - pt.y) / dir.y;
			if (temp >= 0)
				possible_values.push_back(temp);
		}

		if (dir.z > 0)
		{
			double temp = (upper_bounds.z - pt.z) / dir.z;
			if (temp >= 0)
				possible_values.push_back(temp);
		}
		else if (dir.z < 0)
		{
			double temp = (lower_bounds.z - pt.z) / dir.z;
			if (temp >= 0)
				possible_values.push_back(temp);
		}
		if (possible_values.size())
		{
			auto double_it = std::min_element(possible_values.begin(), possible_values.end());
			sf = *double_it;

			// if the scale factor is under the min scale factor remove it from the possible locations
			while (sf < min_sf)
			{
				possible_values.erase(double_it);

				// if all possible factors are too small return false
				if (possible_values.size() == 0)
				{
					return false;
				}

				// update sf to the next smallest value
				double_it = std::min_element(possible_values.begin(), possible_values.end());
				sf = *double_it;
			}
			return true;
		}
		break;
	case FE_Element_Type::FE_TET4G1:
	case FE_Element_Type::FE_TET4G4:
	case FE_Element_Type::FE_TET10G1:
	case FE_Element_Type::FE_TET10G4:
	case FE_Element_Type::FE_TET10G4RI1:
	case FE_Element_Type::FE_TET10G8:
	case FE_Element_Type::FE_TET10G8RI4:
	case FE_Element_Type::FE_TET10GL11:
		{
			assert((pt.x + pt.y + pt.z) < 1.01);

			// double sol0 = (1 - (pt.x + pt.y + pt.z)) / (dir.x + dir.y + dir.z);

			/*
			 * For each face of TET4G4, check to see what scale factor is needed to have the given vector (pt, dir)
			 * reach the plane shared by the face.
			 */

			 /*
			  * FACE 0
			  *    v0: (0, 1, 0)
			  *    v1: (0, 0, 1)
			  *    v2: (1, 0, 0)
			  */
			vec3d relative_v1 = vec3d(0, 0, 1) - vec3d(0, 1, 0); // v1 - v0
			vec3d relative_v2 = vec3d(1, 0, 0) - vec3d(0, 1, 0); // v2 - v0
			mat3d temp0(-dir, relative_v1, relative_v2);
			vec3d sol0;
			double det0 = temp0.det();
			if (det0 != 0.0)
			{
				temp0 = temp0.inverse();
				sol0 = temp0 * (pt - vec3d(0, 1, 0)); // pt - v0
			}

			/*
			 * FACE 1
			 *    v0: (0, 0, 0)
			 *    v1: (0, 0, 1)
			 *    v2: (1, 0, 0)
			 */
			 // relative_v1 = vec3d(0, 0, 1) - vec3d(0, 0, 0); // v1 - v0
			 // relative_v2 = vec3d(1, 0, 0) - vec3d(0, 0, 0); // v2 - v0
			mat3d temp1(-dir, vec3d(0, 0, 1), vec3d(1, 0, 0));
			vec3d sol1;
			double det1 = temp1.det();
			if (det1 != 0.0)
			{
				temp1 = temp1.inverse();
				sol1 = temp1 * (pt); // pt - v0
			}

			/*
			 * FACE 2
			 *    v0: (0, 0, 0)
			 *    v1: (1, 0, 0)
			 *    v2: (0, 1, 0)
			 */
			 // relative_v1 = vec3d(1, 0, 0) - vec3d(0, 0, 0); // v1 - v0
			 // relative_v2 = vec3d(0, 1, 0) - vec3d(0, 0, 0); // v2 - v0
			mat3d temp2(-dir, vec3d(1, 0, 0), vec3d(0, 1, 0));
			vec3d sol2;
			double det2 = temp2.det();
			if (det2 != 0.0)
			{
				temp2 = temp2.inverse();
				sol2 = temp2 * (pt); // pt - v0
			}

			/*
			* FACE 3
			*    v0: (0, 0, 0)
			*    v1: (0, 0, 1)
			*    v2: (0, 1, 0)
			*/
			// relative_v1 = vec3d(0, 0, 1) - vec3d(0, 0, 0); // v1 - v0
			// relative_v2 = vec3d(0, 1, 0) - vec3d(0, 0, 0); // v2 - v0
			mat3d temp3(-dir, vec3d(0, 0, 1), vec3d(0, 1, 0));
			vec3d sol3;
			double det3 = temp3.det();
			if (det3 != 0.0)
			{
				temp3 = temp3.inverse();
				sol3 = temp3 * (pt); // pt - v0
			}

			/*
			 *  Check to make sure that the scale factor is greater than 0 to ensure growth is simulated. Then ensure
			 * that the v1 and v2 components are within the bounds of the face's natural coordinate system. Last,
			 * to ensure that the solution is within the face, the sum of v1 and v2 components must be less than
			 * 1.01 (scaled to account for floating point calculation inaccuracies).
			 *
			 * sol.x: scale factor
			 * sol.y: v1 component of solution
			 * sol.z: v2 component of solution
			 */
			if (sol0.x > min_sf && (sol0.y <= 1.0) && (sol0.y >= 0) && (sol0.z <= 1.0) && (sol0.z >= 0) && ((sol0.y + sol0.z) < 1.01) && det0 != 0.0)
			{
				sf = sol0.x;
				// assert((sol0.y + sol0.z) < 1.01);

				return true;
			}
			if (sol1.x > min_sf && (sol1.y <= 1.0) && (sol1.y >= 0) && (sol1.z <= 1.0) && (sol1.z >= 0) && ((sol1.y + sol1.z) < 1.01) && det1 != 0.0)
			{
				sf = sol1.x;
				// assert((sol1.y + sol1.z) < 1.01);

				return true;
			}
			if (sol2.x > min_sf && (sol2.y <= 1.0) && (sol2.y >= 0) && (sol2.z <= 1.0) && (sol2.z >= 0) && ((sol2.y + sol2.z) < 1.01) && det2 != 0.0)
			{
				sf = sol2.x;
				// assert((sol2.y + sol2.z) < 1.01);

				return true;
			}
			if (sol3.x > min_sf && (sol3.y <= 1.0) && (sol3.y >= 0) && (sol3.z <= 1.0) && (sol3.z >= 0) && ((sol3.y + sol3.z) < 1.01) && det3 != 0.0)
			{
				sf = sol3.x;
				// assert((sol3.y + sol3.z) < 1.01);

				return true;
			}

			return false;
		}
		break;
		default:
			assert(false);
	}
	
	return false;
}

double FEAngio::InElementLength(FESolidElement * se, vec3d pt0, vec3d pt1) const
{
	// Get the mesh
	auto mesh = GetMesh();
	// Get the number of nodes for the element's type
	double H[FESolidElement::MAX_NODES];
	vec3d g_pt0, g_pt1;
	// Get the shape function values at the first point
	se->shape_fnc(H, pt0.x, pt0.y, pt0.z);
	// For each node get the current position of tip 1's positional shape function scale
	for (int j = 0; j < se->Nodes(); j++)
	{
		g_pt0 += mesh->Node(se->m_node[j]).m_rt* H[j];
	}
	// For each node get the current position of tip 2's positional shape function scale
	se->shape_fnc(H, pt1.x, pt1.y, pt1.z);
	for (int j = 0; j < se->Nodes(); j++)
	{
		g_pt1 += mesh->Node(se->m_node[j]).m_rt* H[j];
	}
	// Return the L2 norm (length) of this within the element itself.
	return (g_pt0 - g_pt1).norm();
}

double FEAngio::GetConcentration(FEMaterial* mat, FEMaterialPoint * mp, int sol_id)
{
	FEMultiphasic * multiph = dynamic_cast<FEMultiphasic*>(mat);
	if(multiph)
	{
		 return multiph->Concentration(*mp, sol_id);
	}
	else
	{
		FETriphasic * tri = dynamic_cast<FETriphasic*>(mat);
		if(tri)
		{
			return tri->Concentration(*mp, sol_id);
		}
	}
	assert(false);
	return -1.0;
}

void FEAngio::SetupNodesToElement(int min_element_id)
{
	auto mesh = GetMesh();
	for (int i = 0; i < mesh->Elements(); i++)
	{
		FESolidElement * elem = dynamic_cast<FESolidElement*>(mesh->FindElementFromID(min_element_id + i));
		if(elem)
		{
			for(int j=0;j < elem->Nodes();j++)
			{
				FENode & node = mesh->Node(elem->m_node[j]);
				nodes_to_elements[&node].push_back(elem);
			}
		}
	}
}

bool FEAngio::ProjectToElement(FESolidElement& el, const vec3d& p, FEMesh* mesh, double r[3])
{
	const int MN = FESolidElement::MAX_NODES;
	vec3d rt[MN];

	// get the element nodal coordinates
	int ne = el.Nodes();
	for (int i = 0; i<ne; ++i) rt[i] = mesh->Node(el.m_node[i]).m_rt;

	r[0] = r[1] = r[2] = 0;
	const double tol = 1e-5;
	double dr[3], norm;
	double H[MN], Gr[MN], Gs[MN], Gt[MN];
	int max_iter = 1000000;
	int iter = 0;
	do
	{
		// evaluate shape functions
		el.shape_fnc(H, r[0], r[1], r[2]);

		// evaluate shape function derivatives
		el.shape_deriv(Gr, Gs, Gt, r[0], r[1], r[2]);

		// solve for coordinate increment
		double R[3] = { 0 }, A[3][3] = { 0 };
		for (int i = 0; i<ne; ++i)
		{
			R[0] += rt[i].x*H[i];
			R[1] += rt[i].y*H[i];
			R[2] += rt[i].z*H[i];

			A[0][0] -= rt[i].x*Gr[i]; A[0][1] -= rt[i].x*Gs[i]; A[0][2] -= rt[i].x*Gt[i];
			A[1][0] -= rt[i].y*Gr[i]; A[1][1] -= rt[i].y*Gs[i]; A[1][2] -= rt[i].y*Gt[i];
			A[2][0] -= rt[i].z*Gr[i]; A[2][1] -= rt[i].z*Gs[i]; A[2][2] -= rt[i].z*Gt[i];
		}
		R[0] = p.x - R[0];
		R[1] = p.y - R[1];
		R[2] = p.z - R[2];

		solve_3x3(A, R, dr);
		r[0] -= dr[0];
		r[1] -= dr[1];
		r[2] -= dr[2];

		norm = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
		iter++;
	} while (norm > tol && iter < max_iter);
	return norm <= tol;
}

void FEAngio::GetElementsContainingNode(FENode * node, std::vector<FESolidElement*> & elements)
{
	elements = nodes_to_elements[node];
}

bool FEAngio::IsInBounds(FESolidElement* se, double r[3], double eps)
{
	vec3d nat_pos(r[0], r[1], r[2]);
	switch(se->Type())
	{
	case FE_Element_Type::FE_TET4G1:
	case FE_Element_Type::FE_TET4G4:
	case FE_Element_Type::FE_TET10G1:
	case FE_Element_Type::FE_TET10G4:
	case FE_Element_Type::FE_TET10G4RI1:
	case FE_Element_Type::FE_TET10G8:
	case FE_Element_Type::FE_TET10G8RI4:
	case FE_Element_Type::FE_TET10GL11:

		return (r[0] <= NaturalCoordinatesUpperBound_r(se->Type()) + eps) && (r[0] >= NaturalCoordinatesLowerBound_r(se->Type())-eps) &&
			(r[1] <= NaturalCoordinatesUpperBound_s(se->Type())+ eps) && (r[1] >= NaturalCoordinatesLowerBound_s(se->Type())-eps) &&
			(r[2] <= NaturalCoordinatesUpperBound_t(se->Type())+ eps) && (r[2] >= NaturalCoordinatesLowerBound_t(se->Type())- eps) && ((nat_pos.x + nat_pos.y + nat_pos.z) <= 1 + eps);

	}
	//consider doing this with tolerances
	return (r[0] <= NaturalCoordinatesUpperBound_r(se->Type())+ eps) && (r[0] >= NaturalCoordinatesLowerBound_r(se->Type())-eps) &&
		(r[1] <= NaturalCoordinatesUpperBound_s(se->Type()) + eps) && (r[1] >= NaturalCoordinatesLowerBound_s(se->Type())-eps) &&
		(r[2] <= NaturalCoordinatesUpperBound_t(se->Type()) + eps) && (r[2] >= NaturalCoordinatesLowerBound_t(se->Type())- eps);
}

// 
void FEAngio::ExtremaInElement(FESolidElement * se, std::vector<vec3d> & extrema) const
{
	auto mesh = GetMesh();
	// make extrema the size of the number of nodes in the element type
	extrema.reserve(se->Nodes());
	switch(se->Type())
	{
	case FE_Element_Type::FE_TET4G1:
	case FE_Element_Type::FE_TET4G4:
	case FE_Element_Type::FE_TET10G1:
	case FE_Element_Type::FE_TET10G4:
	case FE_Element_Type::FE_TET10G4RI1:
	case FE_Element_Type::FE_TET10G8:
	case FE_Element_Type::FE_TET10G8RI4:
	case FE_Element_Type::FE_TET10GL11:
	case FE_Element_Type::FE_HEX8G1:
	case FE_Element_Type::FE_HEX8G8:
	case FE_HEX20G27:
	case FE_HEX20G8:
	case FE_HEX27G27:
		// for each node in the element
		for (int j = 0; j < se->Nodes(); j++)
		{
			// add the current global position to extrema
			extrema.push_back(mesh->Node(se->m_node[j]).m_rt);
		}

		break;

	default:
		assert(false);
	}
}

double FEAngio::MinDistance(std::vector<vec3d> & element_bounds0, std::vector<vec3d>  & element_bounds1)
{
	assert(element_bounds0.size());
	assert(element_bounds1.size());
	//just do the sqrt at the end to not have to do it a lot
	double min_dist = (element_bounds0[0] - element_bounds1[0]).norm2();
	for(int i=0;i < element_bounds0.size();i++)
	{
		for(int j=0; j < element_bounds1.size();j++)
		{
			double dist = (element_bounds0[i] - element_bounds1[j]).norm2();
			if (dist < min_dist)
			{
				min_dist = dist;
			}
		}
	}
	return sqrt(min_dist);
}

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the fiber array
	int N = mesh.Nodes();
	fiber.assign(N, vec3d(0,0,0));

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local fiber orientation at integration points
			vector<double> fx(nint), fy(nint), fz(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mpoint->ExtractData<FEElasticMaterialPoint>();

				//Ask Steve about this
				FEMesh* mesh = dom.GetMesh();
				FEElementSet* elset = mesh->FindElementSet(dom.GetName());
				int local_index = elset->GetLocalIndex(el);

				FEMaterial* Mat_a = dom.GetMaterial()->ExtractProperty<FEElasticMaterial>();
				// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
				FEParam* matax = Mat_a->FindParameter("mat_axis");
				FEParamMat3d& p = matax->value<FEParamMat3d>();
				FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
				FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());
				mat3d m_Q = map->valueMat3d(pt);

				mat3d m = m_Q;
				//mat3d m = pt.m_Q;
				
				// grab the first column as the fiber orientation
				fx[n] = m[0][0];
				fy[n] = m[1][0];
				fz[n] = m[2][0];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln), gy(neln), gz(neln);
			el.project_to_nodes(&fx[0], &gx[0]);
			el.project_to_nodes(&fy[0], &gy[0]);
			el.project_to_nodes(&fz[0], &gz[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				vec3d f(gx[i], gy[i], gz[i]);
				fiber[ni] += f;
			}
		}
	}

	// normalize the fibers
	for (int i=0; i<N; ++i) fiber[i].unit();

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	density.resize(N, 0.0);
	anisotropy.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> den(nint);
			vector<double> anis(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			vector<double> fx(neln);
			el.project_to_nodes(&den[0], &gx[0]);
			el.project_to_nodes(&anis[0], &fx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				density[ni] = gx[i];
				if(fx[i] > anisotropy[ni])
					anisotropy[ni] = fx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	concentration.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> con(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FESolutesMaterialPoint& spt = *mpoint->ExtractData<FESolutesMaterialPoint>();
				con[n]=spt.m_c[vegfID];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			el.project_to_nodes(&con[0], &gx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				concentration[ni] = gx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}