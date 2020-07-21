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


 
//! Contains the data that is needed by the plugin on a per element basis. This class mirrors FESolidElements but does not inhertit from FESolidElement. If this inherited from FESolidElement a large number of custom domain classes would need to be created.
class AngioElement
{
public:

	//! constructor
	AngioElement(FESolidElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat, FEMesh * mesh) : _elem(elem), _angio_mat(angio_mat), _mat(mat)
	{
		//active_tips[0][this]
	};

	//! get the length of the segments within the element at a given time
	double GetLengthAtTime(FEMesh* mesh, double time) const;
	//! get the angio fractional anisotropy
	void UpdateAngioFractionalAnisotropy();
	//! Update the angioSPA based on deformation/rotation
	void UpdateSPA();
	//! Construct the elliptical distribution between 2 orthogonal SPA and sample it
	double GetEllipseAngle(const double a, const double b, const double dist_min, const double dist_max, const int n);
	//double GetEllipseAngle2(const double a, const double b);

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
	//this might be further optimized to a constant lookup table given the element type that is doing the looking up ... its possible to do this at compile time
	//this should reduce this to a jump based on the element type
	//! map of tips to the element which they should try to grow into on the next growth substep
	std::unordered_map<AngioElement*, std::vector<Tip *>> active_tips[2];
	//! map of tips to the element which they should try to grow into on the next growth step
	std::unordered_map<AngioElement*, std::vector<Tip *>> next_tips;
	//! to be used in stress calculations
	std::vector<Tip *> current_tips;
	//! to be used in stress calculations
	std::vector<Tip*> final_active_tips;
	//! pointer to any additional information for branching
	BranchInfo * branch_info = nullptr;
	//! count of branches
	int branch_count = 0;
	//! segment length within the global frame within this element
	double global_segment_length = 0.0;
	//! segment length within the reference frame within this element
	double refernce_frame_segment_length = 0.0;
	//! number of anastamoses that have occured
	int anastamoses = 0;
	//! initial orientation of spa
	mat3d initial_angioSPA;
	//! updated angioSPA
	mat3d angioSPA;
	//! Angio fractional anisotropy
	double angioFA;
};
