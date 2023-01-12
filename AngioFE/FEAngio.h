#pragma once
#include "StdAfx.h"
#include "Fileout.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "FECore/FESolidElement.h" //is this include correct or should I just forward declare the class
#include <FEBioLib/FEBioModel.h>
#include <future>
#include "AngioElement.h"
#include <random>
#include <unordered_map>

class FEElemElemList;
//-----------------------------------------------------------------------------
class FEModel;
class FEAngioMaterial;
class FEAngioMaterialBase;

//! The FEAngio class contains all the functionality of the growth model.
class FEAngio
{
public:
	//! constructor
	explicit FEAngio(FEModel& fem);
	//! destructor
	~FEAngio();

	//! initialize the FEAnio class
	bool Init();

	//! returns the min and max element id's within the model
	std::pair<int, int> GetMinMaxElementIDs() const;

	//! Get the FE model
	FEBioModel* GetFEModel() const;

	//! return the mesh associated with this model
	FEMesh* GetMesh() const;

	//! reset all perfomance timers
	void ResetTimers();

	//! get the angio material component of a material if it exists, returns nullptr if there is no angio component
	static FEAngioMaterial* GetAngioComponent(FEMaterial* mat);

	//! gets the global position of the local coordinates given an element
	vec3d Position(FESolidElement* se, vec3d local) const;

	//! updates all of the segment lengths in angio elements
	void CalculateSegmentLengths(FEMesh* mesh);

	//! updates the weighting between the matrix and vessel submaterials
	void AdjustMatrixVesselWeights(FEMesh* mesh);

	//! updates the weighting between the matrix and vessel submaterials
	bool CheckSpecies(FEMesh* mesh);

#ifdef WIN32
	//! generate a rotation matrix in which all rotations are equally probable
	mat3d unifromRandomRotationMatrix(angiofe_random_engine& rengine) const;
#else
	//! needed to compile on gcc see documentation for comments on random
	mat3d unifromRandomRotationMatrix(angiofe_random_engine& rengine);
#endif


	//! generate a rotation matrix from euler angles
	mat3d rotationMatrix(double alpha, double beta, double gamma) const;
	//! returns a random direction(equal probability of each area on the surface of the unit sphere)
	vec3d uniformRandomDirection(angiofe_random_engine& rengine);
	//! returns a random number within the unit cube
	vec3d uniformInUnitCube();

	//! accessors for the DataStore
	double GetDoubleFromDataStore(int record, int elem_id, int item = 0);

	//! calcualtes the gradient at the given natural coordinates, the function variables are stored at the integration points
	static vec3d gradient(FESolidElement* se, std::vector<double>& fn, vec3d pt);

	//these freindships are for displaying/reading the data and are okay
	friend class Fileout;
	friend class FEPlotAngioECMDensity;
	friend class FEPlotAngioRepulseVal;
	friend class FEPlotAngioECMAlpha;
	friend class FEPlotAngioGradient;
	friend class FEPlotBranches;
	friend class FEPlotAnastamoses;
	friend class FEPlotSegmentLength;
	friend class FEPlotRefSegmentLength;
	friend class FEPlotPrimaryVesselDirection;
	//the following friendships are bad and need removed eventually
	//TODO: remove the freindship, creation in the old way requires this
	//or consider making the node and element data public
	friend class FEAngioMaterial;
	friend class NodeDataInterpolationManager;
	friend class NodeDataInterpolation;

	//! Natural coordinate bounds based on the element type.
	static double NaturalCoordinatesUpperBound_r(int et);
	//! Natural coordinates bounds
	static double NaturalCoordinatesUpperBound_s(int et);
	//! Natural coordinates bounds
	static double NaturalCoordinatesUpperBound_t(int et);
	//! Natural coordinates bounds
	static double NaturalCoordinatesLowerBound_r(int et);
	//! Natural coordinates bounds
	static double NaturalCoordinatesLowerBound_s(int et);
	//! Natural coordinates bounds
	static double NaturalCoordinatesLowerBound_t(int et);

	//! clamp the vector to the natural coordinate system
	static vec3d clamp_natc(int et, vec3d natc);

	//! returns the extrema in global coordinates of the given elem
	void ExtremaInElement(FESolidElement* se, std::vector<vec3d>& extrema) const;

	//! returns the minimun distance between the two sets of points
	static double MinDistance(std::vector<vec3d>& element_bounds0, std::vector<vec3d>& element_bounds1);

	//! grow the segments
	void GrowSegments();

	//! grow the segments before t=0
	void ProtoGrowSegments();

	//! get the final active tips in a radius, used in stress calculations
	static void GetActiveFinalTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip*>& tips);
	//! get the tips in a radius, used in stress calculations
	static void GetActiveTipsInRadius(AngioElement* angio_element, double radius, int buffer, FEAngio* pangio, std::vector<Tip*>& tips, int exclude);
	//! get the tips in a radius, used in stress calculations
	static void GetGrownTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip*>& tips);

	//! returns the length between the two points as if they are conencted by a line segment in the natrual coordinates of the element only the difference between the points if the elements are linear
	double InElementLength(FESolidElement* se, vec3d pt0, vec3d pt1) const;

	//! returns the concentration of a given solute at a material point
	static double GetConcentration(FEMaterial* mat, FEMaterialPoint* mp, int sol_id);
private:

	//! Initialize the nodal ECM values
	bool InitECMDensity();

	//! Init FE stuff
	bool InitFEM();

	//! finalize Finite element initialization
	void FinalizeFEM();

	//! populate the adjacency information for the mesh
	void FillInAdjacencyInfo(FEMesh* mesh, FEElemElemList* eel, AngioElement* angio_element, int elem_index);

	//! setup initialization of angio elements(only run once)
	void SetupAngioElements();

	//! set the seeds in angio elements
	void SetSeeds();

	//! seed the tips within the angio elements
	bool SeedFragments();

	//! modify the step size taken by febio
	void ApplydtToTimestepper(double dt);

	//! return the tip position in the reference frame
	vec3d ReferenceCoordinates(Tip* tip) const;

	//! setup the map from nodes to elements
	void SetupNodesToElement(int min_element_id);

	//! do the final output
	void Output();

	//! returns the scale factor that would project this ray onto the unit cube, pt must be within the unit cube
	bool ScaleFactorToProjectToNaturalCoordinates(FESolidElement* se, vec3d& dir, vec3d& pt, double& sf) const;

	//! returns whether or not a ray can be projected onto the surface of an element
	bool ProjectToElement(FESolidElement& el, const vec3d& p, FEMesh* mesh, double r[3]);

	//! populates the vector with elements that contain the node
	void GetElementsContainingNode(FENode* node, std::vector<FESolidElement*>& elements);

	//! returns whether or not the natural coordinates are valid
	static bool IsInBounds(FESolidElement* se, double r[3], double eps = 0.000001);

	//! fill in the adjacency information for an angio element
	void FillInFaces(FEMesh* mesh, AngioElement* angio_element);

	//! the callback that controls the plugin's execution
	static bool feangio_callback(FEModel* pfem, unsigned int nwhen, void* pd)
	{
		FEAngio* pfa = reinterpret_cast<FEAngio*>(pd);
		return pfa->OnCallback(pfem, nwhen);
	}

	//! callback 
	bool OnCallback(FEModel* pfem, unsigned int nwhen);

public:	// parameters read directly from file

	// miscellaneous
	//! Seed number for the random number generator
	unsigned int	m_irseed;

	std::vector<FEAngioMaterial*>	m_pmat;	//!< the angio-material pointer
	std::vector<int>                m_pmat_ids;//!< the material id's

	//! global random engine
	angiofe_random_engine rengine;

	//! pointer to the fbio model
	FEBioModel* m_fem = nullptr;
	//! solid element to angio element
	std::unordered_map<FESolidElement*, AngioElement*> se_to_angio_element;
	//! nodes to elements
	std::unordered_map<FENode*, std::vector<FESolidElement*>> nodes_to_elements;
	//! all node adjacent elements
	std::unordered_map<AngioElement*, std::vector<AngioElement*>> angio_elements_to_all_adjacent_elements;//adjacent is shares a node with the element
	//! elements by material
	std::unordered_map < FEAngioMaterial*, std::vector<AngioElement*>> elements_by_material;

	//! increment fragment number
	int AddFragment();

	//! increment cell number
	int AddCell();

	std::unordered_map<int, FECell*> cells;
private:
	//both nodes and elements id's go from 1 to n+1 for n items
	//first element is padding so the id can be used to lookup the data for that node

	time_t m_start = 0;			// time of start
	Fileout* fileout = nullptr;		// output manager

	std::uniform_real_distribution<double> ztopi;
	std::uniform_real_distribution<double> zto2pi;
	std::uniform_real_distribution<double> n1to1;
	Timer grow_timer;
	Timer update_angio_stress_timer;
	Timer update_branch_policy_timestep_timer;

	std::vector<AngioElement*> angio_elements;//the dense list of angio elements
	std::unordered_map<FESolidElement*, std::pair<AngioElement*, int>> se_to_angio_elem;//the int is the index which is used in neighbor lookups
	std::vector<AngioElement*> angio_elements_with_holes;//the possibly sparse list of elements .. used to serialize data
	std::vector<FEAngioMaterial*> angio_materials;
	int buffer_index = 0;

	double next_time = -1.0;

	const double eps = 0.001;
	static int fragment_id_counter;
	static int cell_id_counter;
	double min_scale_factor = 0.01;
	double bounds_tolerance = 1e-2;
	int growth_substeps = 3;
	double max_angio_dt = 0.25;
	double min_angio_dt = -1.0;
	double min_segment_length = 5.0;
	int auto_stepper_key = -1;
	int bounce = 1;
	// determines whether vessels bounce or grow along a wall. Bounce condition is a symmetry condition.
};
