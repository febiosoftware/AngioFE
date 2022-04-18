#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"
#include <iostream>
#include "FEBioMech/FEElasticMaterialPoint.h"
#include "FECore/FEDomainMap.h"
#include <FECore/mathalg.h>

void BranchPolicy::AddBranchTipEFD(AngioElement * angio_element, vec3d local_pos, Segment* parent_seg, double start_time, int vessel_id, int buffer_index, FEMesh* mesh)
{
	Tip * branch = new Tip();
	branch->TipCell = new FECell();
	branch->time = start_time;
	branch->TipCell->time = start_time;
	branch->seg_time = start_time;
	branch->angio_element = angio_element;
	branch->TipCell->angio_element = angio_element;
	branch->face = angio_element;
	branch->SetLocalPosition(local_pos, mesh);
	branch->initial_fragment_id = vessel_id;
	// create a completely new cell id
	branch->TipCell->initial_cell_id = angio_element->_angio_mat->GetCommonAngioProperties()->fseeder->IncrementCellCounter();
	branch->direction = GetBranchDirectionEFD(local_pos, parent_seg->Direction(mesh), angio_element, mesh);
	branch->use_direction = true;
	branch->is_branch = true;
	branch->parent = parent_seg;

	angio_element->active_tips[(buffer_index + 1 % 2)].at(angio_element).push_back(branch);
	angio_element->branch_count++;
}

vec3d BranchPolicy::GetBranchDirectionEFD(vec3d local_pos, vec3d parent_direction, AngioElement* angio_element, FEMesh* mesh)
{
	/// Get the local SPD ///

	// vector containing the SPD for each gauss point in the element
	mat3ds SPDs_gausspts[8];
	double H[8];
	// get the local EFD
	for (int i = 0; i< angio_element->_elem->GaussPoints(); i++)
	{
		// get the angio point
		FEMaterialPoint* gauss_point = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(gauss_point);
		// Get the SPD   
		angio_mp->UpdateSPD();
		mat3ds temp_SPD = angio_mp->angioSPD;
		SPDs_gausspts[i] = temp_SPD*(3.0 / temp_SPD.tr());
		// TODO: calculate all distances from mp to nodes then normalize. Currently assumes equidistance
		H[i] = 1.0 / angio_element->_elem->GaussPoints();
	}
	// array containing the SPD for each node in the element
	mat3ds SPDs_nodes[8];
	// Get the interpolated SPD from the shape function-weighted Average Structure Tensor
	mat3ds SPD_int = weightedAverageStructureTensor(SPDs_gausspts, H, angio_element->_elem->GaussPoints());
	SPD_int = (3.0 / SPD_int.tr())*SPD_int;

	FEEllipticalDistribution E(this->GetFEModel());
	E.spd = SPD_int;
	E.Init();
	vec3d fiber_direction = E.NextVec(angio_element->_rengine);

	vec3d normal = fiber_direction ^ parent_direction;
	vec3d corrected_fiber_direction = normal ^ parent_direction;//make the fiber direction be ortogonal to parent direction

	quatd zenith_rotation(zenith_angle->GetZenithAngle(local_pos, parent_direction, angio_element), corrected_fiber_direction);
	quatd azimuth_rotation(azimuth_angle->GetAzimuthAngle(local_pos, parent_direction, angio_element), parent_direction);

	return azimuth_rotation.RotationMatrix() * zenith_rotation.RotationMatrix() * parent_direction;
}

//void DelayedBranchingPolicyEFD::AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, double final_time, double min_scale_factor, FEMesh* mesh, FEAngio* feangio)
//{
//	//this algorithm is n^2 there propably exists an nlogn algorithm
//	DelayBranchInfo * dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
//	//calculate where and when the branches shoud occur
//	// check to see if we have processed all segments completely
//	if (angio_elem->processed_recent_segments == angio_elem->recent_segments.size())
//	{
//		return;
//	}
//	//get points of interest. These are the segment endpoints
//	std::set<double> poi;
//	// for each recent segment in the element get the recent segments and get time for the front and back
//	for (int i = angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size(); i++)
//	{
//		Segment* seg = angio_elem->recent_segments[i];
//		poi.insert(seg->front->time);
//		poi.insert(seg->back->time);
//	}
//	//order the points of interest (segments) by time
//	std::vector<double> ordered_poi;
//	ordered_poi.reserve(poi.size());
//	for (auto iter = poi.begin(); iter != poi.end(); ++iter)
//	{
//		ordered_poi.push_back(*iter);
//	}
//	//may do nothing due to using set
//	// sort by the initial time?
//	std::sort(ordered_poi.begin(), ordered_poi.end());
//
//	//construct branch points
//	std::list<BranchPoint> bps;
//	// for each poi (segment) create a branchpoint
//	for (int i = 0; i < ordered_poi.size() - 1; i++)
//	{
//		BranchPoint bpt(ordered_poi[i], ordered_poi[i + 1]);
//		bps.push_back(bpt);
//	}
//	//if the branch belongs to a given segment (checks by comparing each branch start/end time to those of the segments. Seems slow/roundabout and assumes uniqueness)
//	for (int i = angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size(); i++)
//	{
//		Segment * seg = angio_elem->recent_segments[i];
//		// compare each segment to all branch points
//		for (auto iter = bps.begin(); iter != bps.end(); ++iter)
//		{
//			// determine if and which segment the branch belongs to
//			// assign the segments to the branch if they match. Compares start/end times of segments and points of interest
//			if ((seg->front->time >= iter->_end_time && seg->back->time <= iter->_start_time))
//			{
//				iter->current_segments.push_back(seg);
//			}
//		}
//	}
//
//	//compute dependent portion of the BranchPoint
//	// unclear if this only calculates in element length when looking at length to branch.
//
//	// for each branch point 
//	for (auto iter = bps.begin(); iter != bps.end(); ++iter)
//	{
//		double length = 0.0;
//		// get the total length from the segments
//		for (int i = 0; i < iter->current_segments.size(); i++)
//		{
//			Segment * seg = iter->current_segments[i];
//			// Get the length of the segment in this element. 
//			length += feangio->InElementLength(angio_elem->_elem, seg->back->GetLocalPosition(), seg->front->GetLocalPosition());
//		}
//		// assign the length of the segment
//		iter->length = length;
//		// appears to be a way to discretize a vessel in an element into multiple smaller segments to improve curvature. Right now the default discretization length is 1 so this is doing nothing. 
//		iter->discrete_sections = int(floor(length / discretization_length));
//	}
//	angio_elem->processed_recent_segments = int(angio_elem->recent_segments.size());
//	assert(dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info));
//	double & remaining_l2b = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info)->length_to_branch;
//
//	//process all of the branch points, creates the future points at which branches will occur
//#ifndef NDEBUG
//	if (bps.size()>25 || angio_elem->recent_segments.size() > 100)
//	{
//		std::cout << "large number of branch points" << std::endl;
//	}
//#endif
//
//	while (bps.size())
//	{
//		BranchPoint & cur = bps.front();
//		// if there's still stuff to process
//		if (cur.processed < cur.length)
//		{
//			// remaining_length is the difference in the current length and the processed length
//			double remaining_length = cur.length - cur.processed;
//			// if the remaining length is >= the remaining l2b introduce new branch.
//			if (remaining_length >= remaining_l2b)
//			{
//				// determine the distance at which the branch should be inserted.
//				double distance_in = cur.processed + remaining_l2b;
//				// again not sure what this floor operation is for.
//				int actual_section = int(std::floor(distance_in / discretization_length));
//				// Make a new segment from the current segment of interest
//				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size()];
//				// assign the start time based on a weighted mix of the start and end time of the parent segment. The weight is the ratio of the distance until the branch and the total current length.
//				double start_time = mix(cur._end_time, cur._start_time, (distance_in / cur.length));
//				// assign the position of the branch back tip
//				vec3d tip_pos = seg->NatcAtTime(start_time);
//				//add a tip at the appropriate location
//				// This next line was commented out? Get an access violation error with it.
//				//AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index, mesh);
//				double time_to_wait = t2e->NextValue(angio_elem->_rengine);
//				// assign the time that the branch will emerge.
//				time_to_wait = start_time + time_to_wait;
//				//need to test this with multiple steps
//				// if the branch should have emerged but hasn't...
//				if (time_to_wait >= final_time)
//				{
//					// enforce that the branch emerges
//					cur.processed += remaining_l2b;
//					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
//				}
//				// else assign it as a future branch and let that class handle it
//				else
//				{
//					FutureBranch fb = FutureBranch(tip_pos, seg, time_to_wait);
//					dbi->future_branches.push_back(fb);
//					cur.processed += remaining_l2b;
//					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
//				}
//				cur.processed += min_scale_factor;
//			}
//			else
//			{
//				remaining_l2b -= remaining_length;
//				cur.processed = cur.length;
//			}
//		}
//		else
//		{
//			bps.pop_front();
//		}
//	}
//
//	for (auto iter = dbi->future_branches.begin(); iter != dbi->future_branches.end(); ++iter)
//	{
//		while ((iter != dbi->future_branches.end()) && (iter->_start_time < end_time))
//		{
//			//grow this segment
//			//AddBranchTipEFD(angio_elem, iter->_local_pos, iter->_parent->Direction(mesh), iter->_start_time, iter->_parent->GetInitialFragmentID(), buffer_index, mesh);
//			iter = dbi->future_branches.erase(iter);
//		}
//		if (iter == dbi->future_branches.end())
//			break;
//	}
//}

void DelayedBranchingPolicyEFD::SetupBranchInfo(AngioElement * angio_elem)
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