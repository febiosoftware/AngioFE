#pragma once
#include <FECore/FEElement.h>
#include <random>
#include <FECore/FEMaterial.h>
#include <unordered_map>


class FEAngioMaterial;
class Segment;
class Tip;

class AngioElement
{
public:
	AngioElement() {};
	AngioElement(FESolidElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat) : _elem(elem), _angio_mat(angio_mat), _mat(mat)
	{
		//active_tips[0][this]
	};
	//  
	FESolidElement * _elem=nullptr;
	FEAngioMaterial * _angio_mat=nullptr;
	FEMaterial * _mat= nullptr;
	std::mt19937_64 _rengine;

	//begin buffers
	std::vector<AngioElement *> adjacency_list;
	std::vector<Segment *> grown_segments;
	std::vector<Segment *> recent_segments;
	//this might be further optimized to a lookup into a constant lookup table given the element type that is doing the looking up ... its possible to do this at compile time
	//this should reduce this to a jump based on the element type
	std::unordered_map<AngioElement*, Tip *> active_tips[2];
	int active_tips_index = 0;

	std::vector<FESurfaceElement*>  inner_faces;

};
