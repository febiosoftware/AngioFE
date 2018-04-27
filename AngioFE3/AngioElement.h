#pragma once
#include <FECore/FEElement.h>
#include <random>
#include <FECore/FEMaterial.h>
#include <unordered_map>
#include <FECore/FESurface.h>


class FEAngioMaterial;
class Segment;
class Tip;

class AngioElement
{
public:
	AngioElement(FESolidElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat, FEMesh * mesh) : _elem(elem), _angio_mat(angio_mat), _mat(mat), inner_faces(mesh)
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
	std::unordered_map<AngioElement *, int> angio_element_to_adjacency_index;
	std::vector<Segment *> grown_segments;
	std::vector<Segment *> recent_segments;
	//this might be further optimized to a lookup into a constant lookup table given the element type that is doing the looking up ... its possible to do this at compile time
	//this should reduce this to a jump based on the element type
	std::unordered_map<AngioElement*, std::vector<Tip *>> active_tips[2];
	std::unordered_map<AngioElement*, std::vector<Tip *>> next_tips;
	std::vector<Tip *> current_tips;//to be used in stress calculations

	//std::vector<FESurfaceElement*>  inner_faces;
	FESurface inner_faces;

	int padding[16];
};
