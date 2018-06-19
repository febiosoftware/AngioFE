#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"

void BranchPolicy::AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index)
{
	Tip * branch = new Tip();
	branch->time = start_time;
	branch->angio_element = angio_element;
	branch->SetLocalPosition(local_pos);
	branch->initial_fragment_id = vessel_id;
	branch->direction = GetBranchDirection(local_pos, parent_direction);
	branch->use_direction = true;
	branch->is_branch = true;

	angio_element->active_tips[(buffer_index + 1 % 2)].at(angio_element).push_back(branch);
}

vec3d BranchPolicy::GetBranchDirection(vec3d local_pos, vec3d parent_direction)
{
	return vec3d(1, 0, 0);
}



void DelayedBranchingPolicy::AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, FEMesh* mesh, FEAngio* feangio)
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
		iter->discrete_sections = floor(length/discretization_length);
	}
	angio_elem->processed_recent_segments = angio_elem->recent_segments.size();
	assert(dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info));
	double & remaing_l2b = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info)->length_to_branch;

	//process all of the branch points, creates the futre points at which branches will occour
	while(bps.size())
	{
		BranchPoint & cur = bps.front();
		if(cur.processed < cur.length)
		{
			double remaing_length = cur.length - cur.processed;
			if(remaing_length >= remaing_l2b)
			{
				double distance_in = cur.processed + remaing_l2b;
				int actual_section = std::floor(distance_in/discretization_length);
				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size() ];
				double start_time = mix(cur._end_time, cur._start_time, (distance_in / cur.length));
				vec3d tip_pos = seg->NatcAtTime(start_time);
				//add a tip at the appropriate localation
				//AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index);
				double time_to_wait = t2e->NextValue(angio_elem->_rengine);
				FutureBranch fb = FutureBranch(tip_pos, seg, start_time + time_to_wait);
				dbi->future_branches.push_back(fb);
				cur.processed += remaing_l2b;
				remaing_l2b = l2b->NextValue(angio_elem->_rengine);
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
			AddBranchTip(angio_elem, iter->_local_pos, iter->_parent->Direction(mesh), iter->_start_time, iter->_parent->GetInitialFragmentID(), buffer_index);
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