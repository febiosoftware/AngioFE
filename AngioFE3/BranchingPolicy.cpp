#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"

void BranchPolicy::AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index)
{
	Tip * branch = new Tip();
	branch->time = start_time;
	branch->SetLocalPosition(local_pos);
	branch->initial_fragment_id = vessel_id;
	branch->direction = GetBranchDirection(local_pos, parent_direction);
	branch->use_direction = true;
	branch->angio_element = angio_element;
	branch->is_branch = true;

	angio_element->active_tips[buffer_index].at(angio_element).push_back(branch);
}

vec3d BranchPolicy::GetBranchDirection(vec3d local_pos, vec3d parent_direction)
{
	return vec3d(1, 0, 0);
}



void DelayedBranchingPolicy::AddBranches(AngioElement * angio_elem, int buffer_index, FEMesh* mesh, FEAngio* feangio)
{
	//assemble a list of when the number of segments in an element at a time change
	//ordered by time from earliest to latest
	std::list<BranchPoints> bps;
	for(int i=angio_elem->processed_recent_segments; i< angio_elem->recent_segments.size();i++)
	{
		if(bps.size() == 1)
		{
			BranchPoints bpt;
			Segment * seg = angio_elem->recent_segments[i];
			bpt.current_segments.push_back(seg);
			bpt.start_time = seg->back->time;
			bpt.end_time = seg->front->time;
			bps.push_back(bpt);
		}

		//TODO: finish this by processing any number of segments

	}
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

	//process all of the branch points
	while(bps.size())
	{
		auto cur = bps.front();
		if(cur.processed < cur.length)
		{
			double remaing_length = cur.length - cur.processed;
			if(remaing_length >= remaing_l2b)
			{
				double distance_in = cur.processed + remaing_l2b;
				int actual_section = std::floor(distance_in/discretization_length);
				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size() ];
				double start_time = mix(cur.end_time, cur.start_time, (distance_in / cur.length));
				vec3d tip_pos = seg->NatcAtTime(start_time);
				//add a tip at the appropriate localation
				AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index);
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
}
void DelayedBranchingPolicy::SetupBranchInfo(AngioElement * angio_elem)
{
	angio_elem->branch_info = new DelayBranchInfo();
	DelayBranchInfo *dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
	dbi->length_to_branch = l2b->NextValue(angio_elem->_rengine);
}