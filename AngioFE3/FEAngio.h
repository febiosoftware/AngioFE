#pragma once
#include "StdAfx.h"
#include "Fileout.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "FECore/FESolidDomain.h" //isd this include correct or should i just forward declare the class
#include <FEBioLib/FEBioModel.h>
#include <future>
#include "AngioElement.h"


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


	void SetSeeds();

	void GrowSegments();
	// do the final output
	void Output();

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

	std::vector<AngioElement> angio_elements;
	std::unordered_map < FEAngioMaterial *, std::vector<AngioElement *>> elements_by_material;
};
