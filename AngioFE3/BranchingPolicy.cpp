#include "BranchingPolicy.h"
#include "angio3d.h"
#include "Segment.h"
#include "FEAngio.h"
#include <iostream>
#include "FEBioMech/FEElasticMaterialPoint.h"
#include "FECore/FEDomainMap.h"
#include <FECore/mathalg.h>

void BranchPolicy::AddBranchTip(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index, FEMesh* mesh)
{
	Tip * branch = new Tip();
	branch->time = start_time;
	branch->angio_element = angio_element;
	branch->face = angio_element;
	branch->SetLocalPosition(local_pos, mesh);
	branch->initial_fragment_id = vessel_id;
	branch->direction = GetBranchDirection(local_pos, parent_direction, angio_element, mesh);
	branch->use_direction = true;
	branch->is_branch = true;

	angio_element->active_tips[(buffer_index + 1 % 2)].at(angio_element).push_back(branch);
	angio_element->branch_count++;
}

void BranchPolicy::AddBranchTipEFD(AngioElement * angio_element, vec3d local_pos, vec3d parent_direction, double start_time, int vessel_id, int buffer_index, FEMesh* mesh)
{
	Tip * branch = new Tip();
	branch->time = start_time;
	branch->angio_element = angio_element;
	branch->face = angio_element;
	branch->SetLocalPosition(local_pos, mesh);
	branch->initial_fragment_id = vessel_id;
	branch->direction = GetBranchDirectionEFD(local_pos, parent_direction, angio_element, mesh);
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

		//Ask Steve about this
		//get the FE domain
		FEDomain* Dom = dynamic_cast<FEDomain*>(angio_element->_elem->GetMeshPartition());
		//
		FEElementSet* elset = mesh->FindElementSet(Dom->GetName());
		int local_index = elset->GetLocalIndex(*angio_element->_elem);

		FEAngioMaterial* Mat_ang = Dom->GetMaterial()->ExtractProperty<FEAngioMaterial>();
		FEMaterial * Mat_a = Mat_ang->GetMatrixMaterial();
		// assumes that materials mat_axis is already mapped which we'll need to do somewhere else.
		/*FEParam* matax = Mat_a->FindParameter("mat_axis");
		FEParamMat3d& p = matax->value<FEParamMat3d>();
		FEMappedValueMat3d* val = dynamic_cast<FEMappedValueMat3d*>(p.valuator());
		FEDomainMap* map = dynamic_cast<FEDomainMap*>(val->dataMap());*/
		mat3d m_Q = Mat_a->GetLocalCS(*mp);

		axis = emp->m_F * m_Q * axis;
		gauss_data.push_back({ axis });

		/*axis = emp->m_F * emp->m_Q * axis;
		gauss_data.push_back({ axis });*/
	}

	quatd rv = interpolation_prop->Interpolate(angio_element->_elem, gauss_data, local_pos, mesh);
	vec3d fiber_direction = rv.GetVector();
	vec3d normal = fiber_direction ^ parent_direction;
	vec3d corrected_fiber_direction = normal ^ parent_direction;//make the fiber direction be ortogonal to parent direction

	quatd zenith_rotation(zenith_angle->GetZenithAngle(local_pos, parent_direction, angio_element), corrected_fiber_direction);
	quatd azimuth_rotation(azimuth_angle->GetAzimuthAngle(local_pos, parent_direction, angio_element), parent_direction);

	return azimuth_rotation.RotationMatrix() * zenith_rotation.RotationMatrix() * parent_direction;
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

	/// Get the sampled fiber direction ///

	std::vector<pair<double, int>> v;
	mat3d ax;
	ax.setCol(0, vec3d(SPD_int.xx(), SPD_int.xy(), SPD_int.xz()));
	ax.setCol(1, vec3d(SPD_int.xy(), SPD_int.yy(), SPD_int.yz()));
	ax.setCol(2, vec3d(SPD_int.xz(), SPD_int.yz(), SPD_int.zz()));
	v.push_back(pair<double, int>(ax.col(0).norm(), 0));
	v.push_back(pair<double, int>(ax.col(1).norm(), 1));
	v.push_back(pair<double, int>(ax.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	vec3d axis_0 = ax.col(i); axis_0.unit();
	vec3d axis_1 = ax.col(j); axis_1.unit();
	vec3d axis_2 = ax.col(k); axis_2.unit();
	double r0 = (ax.col(i).norm());
	double r1 = (ax.col(j).norm());
	double r2 = (ax.col(k).norm());
	FEEllipticalDistribution E0(this->GetFEModel());
	FEEllipticalDistribution E1(this->GetFEModel());
	E0.a = r0; E0.b = r1; E0.Init();
	E1.a = r0; E1.b = r2; E1.Init();
	double theta_12 = E0.NextValue(angio_element->_rengine);
	double theta_13 = E1.NextValue(angio_element->_rengine);
	
	// rotate the primary direction by theta_12 about the normal between them
	vec3d axis = mix3d_t(axis_0, axis_1, theta_12); axis.unit();
	vec3d fiber_direction = mix3d_t(axis, axis_2, theta_13); fiber_direction.unit();

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
	//calculate where and when the branches shoud occur
	// check if all segments were processed and return if so.
	if(angio_elem->processed_recent_segments == angio_elem->recent_segments.size())
	{
		return;
	}
	//get points of interest
	std::set<double> poi;
	// for each processed segment in the element get the recent segments and get time for the front and back
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
	// for each poi in order create the branch point
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
			// if the branch is supposed to occur after this segment starts but before it finishes then assign the branch to the segment.
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
	// unclear if this only calculates in element length when looking at length to branch.

	// for each branch point 
	for(auto iter= bps.begin(); iter != bps.end();++iter)
	{
		double length = 0.0;
		// for each segment along the branch point calculate the length
		for(int i =0; i < iter->current_segments.size();i++)
		{
			Segment * seg = iter->current_segments[i];
			// Get the length of the segment in this element. Unclear if this means it only considers element length as opposed to the total length.
			length += feangio->InElementLength(angio_elem->_elem, seg->back->GetLocalPosition(), seg->front->GetLocalPosition());
		}
		iter->length = length;
		// appears to be a way to discretize a vessel in an element into multiple smaller segments to improve curvature. Right now the default discretization length is 1 so this is doing nothing. 
		iter->discrete_sections = int (floor(length/discretization_length));
	}
	angio_elem->processed_recent_segments = int (angio_elem->recent_segments.size());
	assert(dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info));
	double & remaining_l2b = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info)->length_to_branch;

	//process all of the branch points, creates the futre points at which branches will occur
#ifndef NDEBUG
	if(bps.size()>25 || angio_elem->recent_segments.size() > 100)
	{
		std::cout << "large number of branch points" << std::endl;
	}
#endif

	while(bps.size())
	{
		BranchPoint & cur = bps.front();
		// what is processed. Not sure what this is comparing here
		// cur.processed is the processed length for the branch while cur.length is the remaining length in the element?
		//std::cout << "cur.processed is " << cur.processed << " cur.length is " << cur.length << endl;
		if(cur.processed < cur.length)
		{
			// remaining_length is the difference in the current length and the processed length
			double remaining_length = cur.length - cur.processed;
			// if the remaining length is >= the remaining l2b introduce new branch.
			if(remaining_length >= remaining_l2b)
			{
				// determine the distance at which the branch should be inserted.
				double distance_in = cur.processed + remaining_l2b;
				// again not sure what this floor operation is for.
				int actual_section = int (std::floor(distance_in/discretization_length));
				// Make a new segment from the current segment of interest
				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size() ];
				// assign the start time based on a weighted mix of the start and end time of the parent segment. The weight is the ratio of the distance until the branch and the total current length.
				double start_time = mix(cur._end_time, cur._start_time, (distance_in / cur.length));
				// assign the position of the branch back tip
				vec3d tip_pos = seg->NatcAtTime(start_time);
				//add a tip at the appropriate location
				//AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index);
				double time_to_wait = t2e->NextValue(angio_elem->_rengine);
				// assign the time that the branch will emerge.
				time_to_wait = start_time + time_to_wait;
				//need to test this with multiple steps
				// if the branch should have emerged but hasn't...
				if(time_to_wait >= final_time)
				{
					// enforce that the branch emerges
					cur.processed += remaining_l2b;
					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				// else assign it as a future branch and let that class handle it
				else
				{
					FutureBranch fb = FutureBranch(tip_pos, seg, time_to_wait);
					dbi->future_branches.push_back(fb);
					cur.processed += remaining_l2b;
					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				cur.processed += min_scale_factor;
				
			}
			else
			{
				remaining_l2b -= remaining_length;
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

void DelayedBranchingPolicyEFD::AddBranches(AngioElement * angio_elem, int buffer_index, double end_time, double final_time, double min_scale_factor, FEMesh* mesh, FEAngio* feangio)
{
	//this algorithm is n^2 there propably exists an nlogn algorithm
	DelayBranchInfo * dbi = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info);
	//calculate where and when the branches shoud occur
	// check if all segments were processed and return if so.
	// if there's no new segments to process return
	if (angio_elem->processed_recent_segments == angio_elem->recent_segments.size())
	{
		return;
	}
	//get points of interest. These are the segment endpoints
	std::set<double> poi;
	// for each processed segment in the element get the recent segments and get time for the front and back
	for (int i = angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size(); i++)
	{
		Segment * seg = angio_elem->recent_segments[i];
		poi.insert(seg->front->time);
		poi.insert(seg->back->time);
	}
	//order the points of interest (segments) by time
	std::vector<double> ordered_poi;
	ordered_poi.reserve(poi.size());
	for (auto iter = poi.begin(); iter != poi.end(); ++iter)
	{
		ordered_poi.push_back(*iter);
	}
	//may do nothing due to using set
	// sort by the initial time?
	std::sort(ordered_poi.begin(), ordered_poi.end());

	//construct branch points
	std::list<BranchPoint> bps;
	// for each poi (segment) create a branchpoint
	for (int i = 0; i < ordered_poi.size() - 1; i++)
	{
		BranchPoint bpt(ordered_poi[i], ordered_poi[i + 1]);
		bps.push_back(bpt);
	}
	//if the branch belongs to a given segment (checks by comparing each branch start/end time to those of the segments. Seems slow/roundabout and assumes uniqueness)
	for (int i = angio_elem->processed_recent_segments; i < angio_elem->recent_segments.size(); i++)
	{
		Segment * seg = angio_elem->recent_segments[i];
		// compare each segment to all branch points
		for (auto iter = bps.begin(); iter != bps.end(); ++iter)
		{
			// determine if and which segment the branch belongs to
			// assign the segments to the branch if they match. Compares start/end times of segments and points of interest
			if (seg->front->time >= iter->_end_time && seg->back->time <= iter->_start_time)
			{
				iter->current_segments.push_back(seg);
			}
		}
	}


	//verify the container with the branch points has been constructed correctly
#ifndef NDEBUG
	for (auto iter = bps.begin(); iter != bps.end(); ++iter)
	{
		auto next = iter;
		++next;
		if (next != bps.end())
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
	// unclear if this only calculates in element length when looking at length to branch.

	// for each branch point 
	for (auto iter = bps.begin(); iter != bps.end(); ++iter)
	{
		double length = 0.0;
		// for each segment to process
		for (int i = 0; i < iter->current_segments.size(); i++)
		{
			Segment * seg = iter->current_segments[i];
			// Get the length of the segment in this element. 
			length += feangio->InElementLength(angio_elem->_elem, seg->back->GetLocalPosition(), seg->front->GetLocalPosition());
		}
		// assign the length of the segment
		iter->length = length;
		// appears to be a way to discretize a vessel in an element into multiple smaller segments to improve curvature. Right now the default discretization length is 1 so this is doing nothing. 
		iter->discrete_sections = int(floor(length / discretization_length));
	}
	angio_elem->processed_recent_segments = int(angio_elem->recent_segments.size());
	assert(dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info));
	double & remaining_l2b = dynamic_cast<DelayBranchInfo*>(angio_elem->branch_info)->length_to_branch;

	//process all of the branch points, creates the future points at which branches will occur
#ifndef NDEBUG
	if (bps.size()>25 || angio_elem->recent_segments.size() > 100)
	{
		std::cout << "large number of branch points" << std::endl;
	}
#endif

	while (bps.size())
	{
		BranchPoint & cur = bps.front();
		// if there's still stuff to process
		if (cur.processed < cur.length)
		{
			// remaining_length is the difference in the current length and the processed length
			double remaining_length = cur.length - cur.processed;
			// if the remaining length is >= the remaining l2b introduce new branch.
			if (remaining_length >= remaining_l2b)
			{
				// determine the distance at which the branch should be inserted.
				double distance_in = cur.processed + remaining_l2b;
				// again not sure what this floor operation is for.
				int actual_section = int(std::floor(distance_in / discretization_length));
				// Make a new segment from the current segment of interest
				Segment * seg = cur.current_segments[actual_section % cur.current_segments.size()];
				// assign the start time based on a weighted mix of the start and end time of the parent segment. The weight is the ratio of the distance until the branch and the total current length.
				double start_time = mix(cur._end_time, cur._start_time, (distance_in / cur.length));
				// assign the position of the branch back tip
				vec3d tip_pos = seg->NatcAtTime(start_time);
				//add a tip at the appropriate location
				// This next line was commented out? Get an access violation error with it.
				//AddBranchTip(angio_elem, tip_pos, seg->Direction(mesh),start_time, seg->GetInitialFragmentID(), buffer_index, mesh);
				double time_to_wait = t2e->NextValue(angio_elem->_rengine);
				// assign the time that the branch will emerge.
				time_to_wait = start_time + time_to_wait;
				//need to test this with multiple steps
				// if the branch should have emerged but hasn't...
				if (time_to_wait >= final_time)
				{
					// enforce that the branch emerges
					cur.processed += remaining_l2b;
					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				// else assign it as a future branch and let that class handle it
				else
				{
					FutureBranch fb = FutureBranch(tip_pos, seg, time_to_wait);
					dbi->future_branches.push_back(fb);
					cur.processed += remaining_l2b;
					remaining_l2b = l2b->NextValue(angio_elem->_rengine);
				}
				cur.processed += min_scale_factor;

			}
			else
			{
				remaining_l2b -= remaining_length;
				cur.processed = cur.length;
			}
		}
		else
		{
			bps.pop_front();
		}
	}

	for (auto iter = dbi->future_branches.begin(); iter != dbi->future_branches.end(); ++iter)
	{
		while ((iter != dbi->future_branches.end()) && (iter->_start_time < end_time))
		{
			//grow this segment
			AddBranchTipEFD(angio_elem, iter->_local_pos, iter->_parent->Direction(mesh), iter->_start_time, iter->_parent->GetInitialFragmentID(), buffer_index, mesh);
			iter = dbi->future_branches.erase(iter);
		}
		if (iter == dbi->future_branches.end())
			break;
	}
}
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