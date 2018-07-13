#pragma once
#include "StdAfx.h"
#include "Fileout.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "FECore/FESolidDomain.h" //isd this include correct or should i just forward declare the class
#include <FEBioLib/FEBioModel.h>
#include <future>
#include "AngioElement.h"
#include <random>


class FEElemElemList;
//-----------------------------------------------------------------------------
class FEModel;
class FEAngioMaterial;
class FEAngioMaterialBase;

//-----------------------------------------------------------------------------
// The FEAngio class contains all the functionality of the growth model.
class FEAngio
{
public:
	explicit FEAngio(FEModel& fem);
	~FEAngio();

	// initialize the FEAnio class
	bool Init();

	std::pair<int, int> GetMinMaxElementIDs() const;

	// Get the FE model
	FEBioModel* GetFEModel() const;

	//check the get FEModel above it may not be useful in any way
	FEMesh * GetMesh() const;

	void ResetTimers();

	static FEAngioMaterial * GetAngioComponent(FEMaterial * mat);
	
	//gets the global position of the local coordinates given an element
	vec3d Position(FESolidElement * se, vec3d local) const;

	void CalculateSegmentLengths(FEMesh* mesh);
	void AdjustMatrixVesselWeights(FEMesh* mesh);


	mat3d unifromRandomRotationMatrix(angiofe_random_engine & rengine) const;
	mat3d rotationMatrix(double alpha, double beta, double gamma);
	vec3d uniformRandomDirection(angiofe_random_engine& rengine);
	vec3d uniformInUnitCube();

	//accessors for the DataStore
	double GetDoubleFromDataStore(int record, int elem_id, int item = 0);

	//calcualtes the gradient at the given natural coordinates
	static vec3d gradient(FESolidElement * se, std::vector<double> & fn, vec3d pt, int size =1,int offset=0);

	//these freindships are for displaying/reading the data and are okay
	friend class Fileout;
	friend class FEPlotAngioECMDensity;
	friend class FEPlotAngioECMAlpha;
	friend class FEPlotAngioGradient;
	friend class FEPlotBranches;
	friend class FEPlotSegmentLength;
	friend class FEPlotRefSegmentLength;
	//the following friendships are bad and need removed eventually
	//TODO: remove the freindship, creation in the old way requires this
	//or consider making the node and element data public
	friend class FEAngioMaterial;
	friend class NodeDataInterpolationManager;
	friend class NodeDataInterpolation;

	// Natural coordinate bounds based on the element type.
	static double NaturalCoordinatesUpperBound_r(int et);
	static double NaturalCoordinatesUpperBound_s(int et);
	static double NaturalCoordinatesUpperBound_t(int et);
	static double NaturalCoordinatesLowerBound_r(int et);
	static double NaturalCoordinatesLowerBound_s(int et);
	static double NaturalCoordinatesLowerBound_t(int et);

	static vec3d clamp_natc(int et, vec3d natc);

	//returns the extrema in global coordinates of the given elem
	void ExtremaInElement(FESolidElement * se, std::vector<vec3d> & extrema) const;

	//returns the minimun distance between the two sets of points
	static double MinDistance(std::vector<vec3d> & element_bounds0, std::vector<vec3d>  & element_bounds1);

	void GrowSegments(double min_scale_factor, double bounds_tolerance, double min_angle, int growth_substeps);

	void ProtoGrowSegments(double min_scale_factor, double bounds_tolerance, double min_angle, int growth_substeps);

	static void GetActiveFinalTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip *> & tips);

	//returns the length between the two points as if they are conencted by a line segment in the natrual coordinates of the element
	//only the difference between the points if the elements are linear
	double InElementLength(FESolidElement * se, vec3d pt0, vec3d pt1) const;

	double GetConcentration(FEMaterial* mat, FEMaterialPoint * mp, int sol_id) const;
private:
	
	// Initialize the nodal ECM values
	bool InitECMDensity();

	// Init FE stuff
	bool InitFEM();

	void FinalizeFEM();

	void FillInAdjacencyInfo(FEMesh * mesh, FEElemElemList * eel, AngioElement *angio_element, int elem_index);

	void SetupAngioElements();

	void SetSeeds();

	bool SeedFragments();

	void ApplydtToTimestepper(double dt, bool initial = false);
	

	vec3d ReferenceCoordinates(Tip * tip) const;

	void SetupNodesToElement(int min_element_id);


	// do the final output
	void Output();

	//returns the scale factor that would project this ray onto the unit cube, pt must be within the unit cube
	bool ScaleFactorToProjectToNaturalCoordinates(FESolidElement* se, vec3d & dir, vec3d & pt, double & sf, double min_sf = 0.00001) const;

	
	

	void GetElementsContainingNode(FENode * node, std::vector<FESolidElement*> & elements);

	static bool IsInBounds(FESolidElement* se, double r[3], double eps= 0.000001);

	

	void FillInFaces(FEMesh * mesh, AngioElement * angio_element);

	static bool feangio_callback(FEModel* pfem, unsigned int nwhen, void* pd)
	{
		FEAngio* pfa = reinterpret_cast<FEAngio*>(pd);
		return pfa->OnCallback(pfem, nwhen);
	}

	bool OnCallback(FEModel* pfem, unsigned int nwhen);

public:	// parameters read directly from file

	// miscellaneous
	unsigned int	m_irseed;			// Seed number for the random number generator
    
	int		total_bdyf;
	int		FE_state;			// State counter to count the number of solved FE states

	std::vector<FEAngioMaterial*>	m_pmat;	//!< the angio-material pointer
	std::vector<int>                m_pmat_ids;

	angiofe_random_engine rengine;

	FEBioModel * m_fem;//just do the cast once
private:
	//both nodes and elements id's go from 1 to n+1 for n items
	//first element is padding so the id can be used to lookup the data for that node

    time_t m_start = 0;			// time of start
	Fileout * fileout = nullptr;		// output manager
	
	std::uniform_real_distribution<double> ztopi;
	std::uniform_real_distribution<double> zto2pi;
	std::uniform_real_distribution<double> n1to1;
	Timer grow_timer;
	Timer mesh_stiffness_timer;
	Timer update_sprout_stress_scaling_timer;
	Timer update_angio_stress_timer;
	Timer update_gdms_timer;
	Timer update_ecm_timer;
	Timer material_update_timer;

	std::vector<AngioElement *> angio_elements;//the dense list of angio elements
	std::unordered_map<FESolidElement*,std::pair<AngioElement *,int>> se_to_angio_elem;//the int is the index which is used in neighbor lookups
	std::vector<AngioElement *> angio_elements_with_holes;//the possibly sparse list of elements .. used to serialize data
	std::unordered_map < FEAngioMaterial *, std::vector<AngioElement *>> elements_by_material;
	std::vector<FEAngioMaterial*> angio_materials;
	int buffer_index = 0;

	std::unordered_map<FESolidElement *, AngioElement *> se_to_angio_element;
	std::unordered_map<FENode *, std::vector<FESolidElement*>> nodes_to_elements;
	std::unordered_map<AngioElement *, std::vector<AngioElement *>> angio_elements_to_all_adjacent_elements;//adjacent is shares a node with the element
	double next_time = -1;


	const double eps = 0.001;

	double min_scale_factor;
	double bounds_tolerance;
	double min_angle;
	int growth_substeps;
	int auto_stepper_key = -1;
};
