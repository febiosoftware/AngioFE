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
	DelayBranchInfo * dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
	//calculate where and when the branches shoud occour
	if(angio_elem->processed_recent_segments == angio_elem->recent_segments.size())
	{
		return;
	}
	//assemble a list of when the number of segments in an element at a time change
	//ordered by time from earliest to latest
	std::list<BranchPoint> bps;

	{
		Segment * seg = angio_elem->recent_segments[angio_elem->processed_recent_segments];
		BranchPoint bpt(seg->back->time, seg->front->time);
		bpt.current_segments.push_back(seg);
		bps.push_back(bpt);
	}

	
	for(int i=angio_elem->processed_recent_segments + 1; i< angio_elem->recent_segments.size();i++)
	{
		bool last_placed = false;
		bool first_placed = true;
		auto iter = bps.begin();
		for(; iter != bps.end();++iter)
		{
			Segment * seg = angio_elem->recent_segments[i];
			//insert leading branch point
			if(iter->_start_time > seg->back->time)
			{
				//break up the bpt
				BranchPoint bpt(seg->back->time, iter->_start_time);
				bpt.current_segments.push_back(seg);
				bps.insert(iter, bpt);
				++iter;
				//does this work in high density situations?
				first_placed = true;
				break;
			}
			//do a middle split
			else if(iter->_start_time < seg->back->time && iter->_end_time > seg->back->time)
			{
				BranchPoint bpt(seg->back->time, iter->_end_time);
				iter->_end_time = seg->back->time;
				for(int j=0; j < iter->current_segments.size();j++)
				{
					bpt.current_segments.push_back(iter->current_segments[j]);
				}
				++iter;
				bps.insert(iter, bpt);
				//does this work in high density situations?
				first_placed = true;
				break;
			}
			//well aligned segments
			else if(iter->_start_time == seg->back->time)
			{
				//iter now points to the correct bpt 
				first_placed = true;
				break;
			}
		}
		//now advance to the end
		for (; iter != bps.end() && first_placed; ++iter)
		{
			//still need to place
			Segment * seg = angio_elem->recent_segments[i];
			if (iter->_end_time > seg->front->time)
			{
				//break up the bpt
				
				BranchPoint bpt(seg->front->time, iter->_end_time);
				
				iter->_end_time = seg->front->time;
				for(int j=0; j < iter->current_segments.size();j++)
				{
					bpt.current_segments.push_back(iter->current_segments[j]);
				}
				iter->current_segments.push_back(seg);
				auto next = iter;
				++next;
				bps.insert(next, bpt);
				last_placed = true;
				break;
			}
			else if (iter->_end_time == seg->front->time)
			{
				//add the segment to the current bucket
				iter->current_segments.push_back(seg);
				//iter now points to the correct bpt 
				last_placed = true;
				break;
			}
			else
			{
				//add the segment to the current bucket
				iter->current_segments.push_back(seg);
				last_placed = true;
			}
		}
		if (!last_placed)
		{
			/*
			auto prev = iter;
			--prev;
			Segment * seg = angio_elem->recent_segments[i];
			if(seg->back->time == prev->_start_time && seg->front->time == prev->_end_time)
			{
				prev->current_segments.push_back(seg);
			}
			//do a split

			//trailing disconnected segment
			else
			{
				BranchPoint bpt(seg->back->time, seg->front->time);
				bps.insert(iter, bpt);
			}
			*/
			
		}
	}

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

	//process all of the branch points
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