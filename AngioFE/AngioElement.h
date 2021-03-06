///////////////////////////////////////////////////////////////////////
// AngioElement.h
// Used for operations based at the elemental level. Includes growth 
// handling. This class mirrors FESolidElements (but does not inherit 
// from FESolidElement).
///////////////////////////////////////////////////////////////////////

#pragma once
#include <FECore/FESolidElement.h>
#include <random>
#include <FECore/FEMaterial.h>
#include <unordered_map>
#include <FECore/FESurface.h>
#include <unordered_set>

class FEAngioMaterial;
class Segment;
class Tip;
class BranchInfo;
class FECell;
 
// Base class
class AngioElement
{
public:

	//! constructor
	AngioElement(FESolidElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat, FEMesh * mesh) : _elem(elem), _angio_mat(angio_mat), _mat(mat)
	{

	};

	//! get the length of the segments within the element at a given time
	double GetLengthAtTime(FEMesh* mesh, double time) const;

	//! pointer to element  
	FESolidElement * _elem = nullptr;
	//! pointer to angio material
	FEAngioMaterial * _angio_mat = nullptr;
	//! pointer to top level material(may be angio material or multiphasic material)
	FEMaterial * _mat = nullptr;
	//! random engine for the elment (sequence all random rolls within the element sequentially wrt time)
	std::mt19937_64 _rengine;

	//begin buffers
	//! elements that share a node with the element
	std::vector<AngioElement *> adjacency_list;
	//! elements that share a face with the element
	std::vector<AngioElement *> face_adjacency_list;
	//! segments that are contained within this element 
	std::vector<Segment *> grown_segments;
	//! segments that have grown within this element since the previous growth step
	std::vector<Segment *> recent_segments;
	//! count to make sure all recent segments are handled correctly
	int processed_recent_segments = 0;
	// MP4: This might be further optimized to a constant lookup table given the element type that is doing the looking up ... its possible to do this at compile time this should reduce this to a jump based on the element type
	//! map of tips to the element which they should try to grow into on the next growth substep
	std::unordered_map<AngioElement*, std::vector<Tip *>> active_tips[2];
	//! map of tips to the element which they should try to grow into on the next growth step
	std::unordered_map<AngioElement*, std::vector<Tip *>> next_tips;
	//! to be used in stress calculations
	std::vector<Tip *> current_tips;
	//! to be used in stress calculations
	std::vector<Tip*> final_active_tips;
	//! to be used for I/O and looking at just the tip cells
	//std::vector<FECell *> tip_cells;
	std::unordered_map<int, FECell*> tip_cells;
	//! pointer to any additional information for branching
	BranchInfo * branch_info = nullptr;
	//! count of branches
	int branch_count = 0;
	//! segments containing branch points that have not been added
	std::vector<Segment*> active_branch_segs;
	//! segment length within the global frame within this element
	double global_segment_length = 0.0;
	//! segment length within the reference frame within this element
	double reference_frame_segment_length = 0.0;
	//! number of anastamoses that have occured
	int anastamoses = 0;
	//! time that the element becomes overly vascularized
	double vasc_thresh_time = -1.0;
	//! element vessel weight
	double vessel_weight;
};
