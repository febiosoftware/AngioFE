#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"
#include <iostream>

void BranchPolicy::AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index, FEMesh* mesh)
{
	Tip * branch = new Tip();
	branch->time = start_time;
	branch->angio_element = angio_element;
	branch->face = angio_element;
	branch->SetLocalPosition(local_pos, mesh);
	branch->initial_fragment_id = vessel_id;
	branch->direction = GetBranchDirection(local_pos, parent_direction,angio_element,mesh);
	branch->use_direction = true;
	branch->is_branch = true;

	angio_element->active_tips[(buffer_index + 1 % 2)].at(angio_element).push_back(branch);
	angio_element->branch_count++;
}

vec3d BranchPolicy::GetBranchDirection(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh)
{
	//get the fiber direction to form a basis to determine the branch directions
	std::vector<quatd> gauss_data;
	std::vector<quatd> nodal_data;
	for (int i = 0; i< angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		vec3d axis(1, 0, 0);
		axis = emp->m_F * emp->m_Q * axis;
		gauss_data.push_back({ axis });
	}

	quatd rv = interpolation_prop->Interpolate(angio_element->_elem, gauss_data, local_pos, mesh);
	vec3d fiber_direction = rv.GetVector();
	vec3d normal = fiber_direction ^ parent_direction;
	vec3d corrected_fiber_direction = normal ^ parent_direction;//make the fiber direction be ortogonal to parent direction

	quatd zenith_rotation(zenith_angle->GetZenithAngle(local_pos, parent_direction, angio_element), corrected_fiber_direction);
	quatd azimuth_rotation(azimuth_angle->GetAzimuthAngle(local_pos, parent_direction, angio_element), parent_direction);

	return azimuth_rotation.RotationMatrix() * zenith_rotation.RotationMatrix() * parent_direction;
}



void DelayedBranchingPolicy::AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, double final_time, double min_scale_factor, FEMesh* mesh, FEAngio* feangio)
{
	//this algorithm is n^2 there propably exists an nlogn algorithm
	DelayBranchInfo * dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
	//calculate where and when the branches shoud occour
	if(angio_elem->processed_recent_segments == angio_elem->recent_segments.size())
	{
		return;
	}
	//get points of interest 
	std::set<double> poi;
	for(int i= angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size();i++)
	{
		Segment * seg = angio_elem->recent_segments[i];
		poi.insert(seg->front->time);
		poi.insert(seg->back->time);
	}
	//order the points of interst
	std::vector<double> ordered_poi;
	ordered_poi.reserve(poi.size());
	for(auto iter = poi.begin(); iter != poi.end(); ++iter)
	{
		ordered_poi.push_back(*iter);
	}
	//may do nothing due to using set
	std::sort(ordered_poi.begin(), ordered_poi.end());

	//construct branch points
	std::list<BranchPoint> bps;
	for(int i=0; i < ordered_poi.size() -1;i++)
	{
		BranchPoint bpt(ordered_poi[i], ordered_poi[i+1]);
		bps.push_back(bpt);
	}
	//add segments to all branch points that contain them
	for (int i = angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size(); i++)
	{
		Segment * seg = angio_elem->recent_segments[i];
		for(auto iter = bps.begin(); iter != bps.end(); ++iter)
		{
			if(seg->front->time >= iter->_end_time && seg->back->time <= iter->_start_time)
			{
				iter->current_segments.push_back(seg);
			}
		}
	}


	//verify the container with the branch points has been constructed correctly
#ifndef NDEBUG
	for(auto iter = bps.begin(); iter != bps.end();++iter)
	{
		auto next = iter;
		++next;
		if(next != bps.end())
		{
			assert(iter->_end_time <= next->_start_time);
		}
		assert(iter->_start_time < iter->_end_time);
		for (int i = 0; i < iter->current_segments.size(); i++)
		{
			//check the segment is within the time bounds
			assert(iter->current_segments[i]->front->time >= iter->_end_time);
			assert(iter->current_segments[i]->back->time <= iter->_start_time);
		}
	}
#endif

	//compute dependent portion of the BranchPoint
	for(auto iter= bps.begin(); iter != bps.end();++iter)
	{
		double length = 0.0;
		for(int i =0; i < iter->current_segments.size();i++)
		{
			Segment * seg = iter->current_segments[i];
			
			length += feangio->InElementLength(angio_elem->_elem, seg->back->GetLocalPosition(), seg->front->GetLocalPosition());
		}
		iter->length = length;
		iter->discrete_sections = int (floor(length/discretization_length));
	}
	angio_elem->processed_recent_segments = int (angio_elem->recent_segments.size());
	assert(dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info));
	double & remaing_l2b = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info)->length_to_branch;

	//process all of the branch points, creates the futre points at which branches will occour
#ifndef NDEBUG
	if(bps.size()>25 || angio_elem->recent_segments.size() > 100)
	{
		std::cout << "large number of branch points" << std::endl;
	}
#endif

	while(bps.size())
	{
		BranchPoint & cur = bps.front();
		if(cur.processed < cur.length)
		{
			double remaing_length = cur.length - cur.processed;
			if(remaing_length >= remaing_l2b)
			{
				double distance_in = cur.processed + remaing_l2b;
				int actual_section = int (std::floor(distance_in/discretization_length));
				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size() ];
				double start_time = mix(cur._end_time, cur._start_time, (distance_in / cur.length));
				vec3d tip_pos = seg->NatcAtTime(start_time);
				//add a tip at the appropriate localation
				//AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index);
				double time_to_wait = t2e->NextValue(angio_elem->_rengine);
				time_to_wait = start_time + time_to_wait;
				//need to test this with multiple steps
				if(time_to_wait >= final_time)
				{
					cur.processed += remaing_l2b;
					remaing_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				else
				{
					FutureBranch fb = FutureBranch(tip_pos, seg, time_to_wait);
					dbi->future_branches.push_back(fb);
					cur.processed += remaing_l2b;
					remaing_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				cur.processed += min_scale_factor;
				
			}
			else
			{
				remaing_l2b -= remaing_length;
				cur.processed = cur.length;
			}
		}
		else
		{
			bps.pop_front();
		}
	}

	for(auto iter = dbi->future_branches.begin(); iter != dbi->future_branches.end(); ++iter)
	{
		while((iter != dbi->future_branches.end()) && (iter->_start_time < end_time))
		{
			//grow this segment
			AddBranchTip(angio_elem, iter->_local_pos, iter->_parent->Direction(mesh), iter->_start_time, iter->_parent->GetInitialFragmentID(), buffer_index,mesh);
			iter = dbi->future_branches.erase(iter);
		}
		if (iter == dbi->future_branches.end())
			break;
	}
}
void DelayedBranchingPolicy::SetupBranchInfo(AngioElement * angio_elem)
{
	angio_elem->branch_info = new DelayBranchInfo();
	DelayBranchInfo *dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
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