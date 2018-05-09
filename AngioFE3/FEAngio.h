#pragma once
#include "StdAfx.h"
#include "Fileout.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "FECore/FESolidDomain.h" //isd this include correct or should i just forward declare the class
#include <FEBioLib/FEBioModel.h>
#include <future>
#include "AngioElement.h"


class FEElemElemList;
//-----------------------------------------------------------------------------
class FEModel;
class FEAngioMaterial;
class FEAngioMaterialBase;

//-----------------------------------------------------------------------------
// This class represents the time parameters
class SimulationTime
{
public:
	double	t;		// current time value
	double	dt;		// time increment from last time
	double	maxt;	// end time of simulation
	SimulationTime() { t = dt = maxt = 0.0; }
};

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


	mat3d unifromRandomRotationMatrix();
	mat3d rotationMatrix(double alpha, double beta, double gamma);
	vec3d uniformRandomDirection();
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
	//the following friendships are bad and need removed eventually
	//TODO: remove the freindship, creation in the old way requires this
	//or consider making the node and element data public
	friend class FEAngioMaterial;
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

	void ApplydtToTimestepper(double dt);

	void GrowSegments();

	vec3d ReferenceCoordinates(Tip * tip) const;
	// do the final output
	void Output();

	//returns the scale factor that would project this direction onto the unit cube if the direction was a ray that started at pt, pt must be within the unit cube
	double ScaleFactorToProjectToUnitCube(vec3d & dir, vec3d & pt) const;

	//returns the length between the two points as if they are conencted by a line segment in the natrual coordinates of the element
	//only the difference between the points if the elements are linear
	double InElementLength(FESolidElement * se, vec3d pt0, vec3d pt1) const;

	bool IsInBounds(double r[3]);

	static bool feangio_callback(FEModel* pfem, unsigned int nwhen, void* pd)
	{
		FEAngio* pfa = reinterpret_cast<FEAngio*>(pd);
		pfa->OnCallback(pfem, nwhen);
		return true;
	}

	void OnCallback(FEModel* pfem, unsigned int nwhen);

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

	SimulationTime	m_time;		// simulation time

    time_t m_start = 0;			// time of start
	Fileout * fileout = nullptr;		// output manager
	
	std::uniform_real_distribution<double> ztopi;
	std::uniform_real_distribution<double> zto2pi;
	std::uniform_real_distribution<double> n1to1;
	Timer grow_timer;
	Timer mesh_stiffness_timer;
	Timer update_sprout_stress_scaling_timer;
	Timer update_gdms_timer;
	Timer update_ecm_timer;
	Timer material_update_timer;

	std::vector<AngioElement *> angio_elements;//the dense list of angio elements
	std::unordered_map<FESolidElement*,std::pair<AngioElement *,int>> se_to_angio_elem;//the int is the index which is used in neighbor lookups
	std::vector<AngioElement *> angio_elements_with_holes;//the possibly sparse list of elements .. used to serialize data
	std::unordered_map < FEAngioMaterial *, std::vector<AngioElement *>> elements_by_material;
	int buffer_index = 0;

	double next_time = 0.0;

};
