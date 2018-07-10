#pragma once
#include <FECore/FEElement.h>
#include <random>
#include <FECore/FEMaterial.h>
#include <unordered_map>
#include <FECore/FESurface.h>
#include <unordered_set>


class FEAngioMaterial;
class Segment;
class Tip;
class BranchInfo;

class AngioElement
{
public:
	AngioElement(FESolidElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat, FEMesh * mesh) : _elem(elem), _angio_mat(angio_mat), _mat(mat)
	{
		//active_tips[0][this]
	};
	//  
	FESolidElement * _elem=nullptr;
	FEAngioMaterial * _angio_mat=nullptr;
	FEMaterial * _mat= nullptr;
	std::mt19937_64 _rengine;

	//begin buffers
	std::vector<AngioElement *> adjacency_list;//elements that share a node with the element
	std::vector<AngioElement *> face_adjacency_list;//elements that share a face with the element
	std::vector<Segment *> grown_segments;
	std::vector<Segment *> recent_segments;
	int processed_recent_segments = 0;
	//this might be further optimized to a lookup into a constant lookup table given the element type that is doing the looking up ... its possible to do this at compile time
	//this should reduce this to a jump based on the element type
	std::unordered_map<AngioElement*, std::vector<Tip *>> active_tips[2];
	std::unordered_map<AngioElement*, std::vector<Tip *>> next_tips;
	std::vector<Tip *> current_tips;//to be used in stress calculations

	std::vector<Tip*> final_active_tips;

	BranchInfo * branch_info = nullptr;
	int branch_count = 0;
	double global_segment_length = 0.0;
};
