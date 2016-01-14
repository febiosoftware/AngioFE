///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include "FECore\FESolidDomain.h"
#include "FEBioMech\FEElasticMaterial.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "Elem.h"
#include "angio3d.h"
#include "Culture.h"
#include <iostream>


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, FEMaterial* pmat);

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : m_fem(fem), m_grid(fem.GetMesh())
{
	// Create the culture
	m_pCult = new Culture(*this);

	// Body force counter
	total_bdyf = 0;
	
	FE_state = 0;

	comp_mat = 0;
	m_pmat = 0;

	phi_stiff_factor = 1.0;
    m_ntime = 1;

	// flag for generating fibers (0 = random, 3 = element orientation)
	m_matrix_cond = 0;

	// flatten fiber option (default to false)
	m_bzfibflat = 0;

	// initialize time stepping parameters
	m_time.dt = 0.25;

	// TODO: What are these and make these user parameters
    m_dtA = 0.0637;
    m_dtB = 9.0957;
    m_dtC = 2.6073;

	// vessel_width - Diameter of microvessels (Default: 7 um)
	m_vessel_width = 7;

	m_bsprout_verify = 0;				// Sprout verification problem flag

	// boundary conditions
	m_cgelbc = 'u';					// Gel boundary conditions ('u' unconstrained)
	m_Sx = 0.;						// Location of the x-symmetry plane
	m_Sy = 0.;						// Location of the y-symmetry plane
	m_Sz = 0.;						// Location of the z-symmetry plane

	// Input random seed number
	m_irseed = (unsigned int) time(0);

	// default sprout force parameters
   	m_sproutf   = 1.0;		// Sprout force magnitude	
	m_tip_range = 250.;		// Sprout force range
	m_spfactor = 0.;		// Sprout force directional factor
	m_bsp_sphere = 0;		// Switch between local directional (0), local isotropic (1), and global isotropic (2) sprout froce representations
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
	delete m_pCult;
}

//-----------------------------------------------------------------------------
FEModel& FEAngio::GetFEModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
// find the angio material component
FEAngioMaterial* FindAngioMaterial(FEMaterial* pm)
{
	FEMaterial* pmat = pm->FindComponentByType("angio");
	if (pmat)
	{
		FEAngioMaterial* pma = dynamic_cast<FEAngioMaterial*>(pmat);
		return pma;
	}
	return 0;
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	// Seed the random number generator
	srand(m_irseed);

	// Print out the seed number for the random generator
	fileout.printrandseed(m_irseed);					

	// Initialize Culture class
	if (m_pCult->Init() == false) return false;

	// create the grid based on the FEBio mesh
	if (m_grid.Init() == false) return false;

	// assign ECM densities to grid nodes
	if (InitECMDensity() == false) return false;

	// assign collagen fibers to grid nodes
	if (InitCollagenFibers() == false) return false;

	// Seed initial fragments 
	// NOTE: must be done after InitECMDensity() and InitCollagenFibers().
	m_pCult->SeedFragments(m_time);

	// Init all the FE stuff
	if (InitFEM() == false) return false;
	
	// start timer
	time(&m_start);

	return true;
}

//-----------------------------------------------------------------------------
double FEAngio::RunTime()
{
	time_t stop;
    time(&stop);
	return (double) difftime(stop, m_start);
}

//-----------------------------------------------------------------------------
// Initialize FE model.
bool FEAngio::InitFEM()
{
	// See if an "angio" material is defined.
	bool bmat = true;
	FEAngioMaterial* pma = FindAngioMaterial(m_fem.GetMaterial(0));
	if (pma == 0) bmat = false;
	else 
	{
		if (m_Sx != 0.){ pma->Sx = m_Sx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (m_Sy != 0.){ pma->Sy = m_Sy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (m_Sz != 0.){ pma->Sz = m_Sz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pma->ApplySym();

		m_pmat = pma;
		update_sprout_stress_scaling();
	}

	// If the angio material is not defined we apply the "old" body force approach
	if (bmat == false)
	{
		FESproutBodyForce* pbf = new FESproutBodyForce(&m_fem);			// Define the one-and-only bodyforce	
		m_pbf = pbf; // --- SAM ---
	
		m_fem.AddBodyLoad(pbf);										// Add the body force to the FEmodel
		FEParameterList& pl = pbf->GetParameterList();					// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);							// Get the magnitude parameter
		FEParam* pb = pl.Find("b"); assert(pb);							// Get the range parameter
	
		pa->value<double>() = (1.0/4.0)*0.001*m_sproutf;	// Set the mangnitude parameter using the input file
		pb->value<double>() = 1.0/m_tip_range;				// Set the range parameter using the input file
		pbf->m_factor = m_spfactor;							// Set the directional factor parameter using the input file
	
		if (m_bsp_sphere == 1) pbf->m_factor = 0.;				// If using local isotropic sprout force, set directional factor to 0

		if (m_bsp_sphere == 2){									// If using global isotropic sprout froce, set directional factor and range to 0 
			pbf->m_factor = 0.; 
			pb->value<double>() = 0;}					

		//cout << pbf->m_factor << endl;									// Print out the directional factor parameter

		if (m_Sx != 0.){ pbf->Sx = m_Sx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (m_Sy != 0.){ pbf->Sy = m_Sy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (m_Sz != 0.){ pbf->Sz = m_Sz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pbf->ApplySym();												// Apply any symmetry to the bodyforce class
	}

	//// FEBIO - Initialize the FE model
	m_pmat->SetFEAngio(this);

	// report if the stress or body force approach will be used
	if (bmat)
		felog.printf("Angio material found. Stress approach will be used.");
	else
		felog.printf("Angio materia NOT found. Body-force appraoch will be used.");

	// register the callback
	m_fem.AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED, this);

	// Do the model initialization
	if (m_fem.Init() == false) return false;

	// apply the intial sprout forces
	apply_sprout_forces(1, 0.5);

	// Adjust the stiffness of the mesh based on microvessel volume
	adjust_mesh_stiffness();

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::FILE_ONLY);

	// --- Output initial state of model ---

	// Output initial microvessel state
	fileout.save_vessel_state(*this);

	// Output time information
	fileout.save_time(*this);

	// Output initial collagen fiber orientation
	fileout.writeCollFib(GetGrid(), true);

	return true;
}

//-----------------------------------------------------------------------------
// Initialize the nodal ECM values
bool FEAngio::InitECMDensity()
{
	int NN = m_grid.Nodes();
	vector<double> density(NN, 0.0);

	if (m_grid.m_coll_den == 0.0)
	{
		// get the material
		FEMaterial* pm = m_fem.GetMaterial(0);
		FEMaterial* pmat = pm->FindComponentByType("angio");
		if (pmat == 0) return false;

		if (CreateDensityMap(density, pmat) == false) return false;
	}
	else
	{
		for (int i=0; i<NN; ++i) density[i] = m_grid.m_coll_den;
	}

	// assign ECM density
	for (int i = 0; i < NN; ++i)								
	{
		Node& node = m_grid.GetNode(i);

		node.m_ecm_den0 = density[i];	
		node.m_ecm_den = node.m_ecm_den0;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Initialize nodal collagen fiber directions
bool FEAngio::InitCollagenFibers()
{
	int NN = m_grid.Nodes();
	vector<vec3d> fiber;
	fiber.resize(NN, vec3d(0,0,0));

	switch (m_matrix_cond)
	{
	case 0: // random orientation
		{
			for (int i=0; i<NN; ++i)
			{
				fiber[i] = vrand();
			}
		}
		break;
	case 3:	// from element's local coordinates
		{
			// get the material
			FEMaterial* pm = m_fem.GetMaterial(0);
			FEMaterial* efd = pm->FindComponentByType("EFD neo-Hookean");
			if (efd == 0) return false;
			else
			{
				if (CreateFiberMap(fiber, efd) == false) return false;
			}
		}
	}

	// assign collagen fibers
	for (int i=0; i<NN; ++i)
	{
		Node& node = m_grid.GetNode(i);

		vec3d v = fiber[i];

		// flatten if requested
		if (m_bzfibflat == 1) v.z *= 0.25;

		// normalize the vector
		v.unit();

		// assign the node
		node.m_collfib0 = v;
		node.m_collfib = node.m_collfib0;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;

	if (nwhen == CB_UPDATE_TIME)
	{
		// grab the time information
		m_time.t = fem.m_ftime;
		m_time.dt = fem.GetCurrentStep()->m_dt;

		// do a growth step
		m_pCult->Grow(m_time);

		// update sprout stress scaling
		update_sprout_stress_scaling();

		// Update the positions of the body forces
		update_body_forces(1.0);
	}
	else if (nwhen == CB_MAJOR_ITERS)
	{
		// update the grid data
		m_grid.Update();

		// update the culture
		m_pCult->Update();

		++FE_state;

		// Save the current vessel state
		fileout.save_vessel_state(*this);

		// Output time information	
		fileout.save_time(*this);
		
		// Print the status of angio3d to the user    
		fileout.printStatus(*this);
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		Output();
	}
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{
	// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	fileout.output_params(*this);
	
	// Output data file
	fileout.dataout(*this);
		
	// Output final collagen fiber orientation
	fileout.writeCollFib(GetGrid(), false);

	// Output final matrix density
	fileout.writeECMDen(GetGrid());
}

//-----------------------------------------------------------------------------
// Apply sprout forces to the mesh for each active vessel tip
void FEAngio::apply_sprout_forces(int load_curve, double scale)
{
	double magnitude = scale*m_sproutf;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

	//#pragma omp parallel for
	for (list<SegIter>::iterator tip_it = m_pCult->m_active_tips.begin(); tip_it != m_pCult->m_active_tips.end(); ++tip_it)		// For each active growth tip...
	{
		Segment& seg = (*(*tip_it));												// Obtain the growth tip

		if (seg.tip(0).bactive)
		{
			vec3d tip = seg.tip(0).pos();												// Obtain the position of the active tip
			
			// Calculate the directional unit vector of the sprout (notice negative sign)
			vec3d sprout_vect = -seg.uvect();

			(*tip_it)->tip(0).bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, load_curve);				// Create a new body force, set the tips body force ID
		}
		
		if (seg.tip(1).bactive)
		{	
			vec3d tip = seg.tip(1).pos();												// Obtain the position of the active tip
			
			// Calculate the directional unit vector of the sprout
			vec3d sprout_vect = seg.uvect();

			(*tip_it)->tip(1).bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, load_curve);				// Create a new body force, set the tips body force ID
		}
	}

	return;
}

//-----------------------------------------------------------------------------
// Add a new body force entry into the body force field applyied to the mesh
int FEAngio::create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, int load_curve)
{
	total_bdyf++;							// Iterate the total body force counter							

	if (m_pmat)
	{
		m_pmat->AddSprout(vec3d(xpt, ypt, zpt), sprout_vect);
		return m_pmat->Sprouts() - 1;
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));					// Obtain the body force class
		pbf->AddSprout(vec3d(xpt, ypt, zpt),sprout_vect);												// Add a new component to the body force for this active sprout tip
		return pbf->Sprouts() - 1;																		// Return the ID number for the body force
	}
}


//-----------------------------------------------------------------------------
// Update the sprout forces after a deformation
void FEAngio::update_body_forces(double scale)
{
	vec3d sprout_vect;												// Sprout direction vector

	vec3d tip(0,0,0);
	double magnitude = scale*m_sproutf;								// Magnitude of the sprout force

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

	if (m_pmat)
	{
		//m_pmat->scale = magnitude;
		int NSP = m_pmat->Sprouts();
		for (int i=0; i<NSP; ++i)
		{
			FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(i);
			sp.bactive = false;
		}
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));			// Obtain the sprout body force field
		int NSP = pbf->Sprouts();										// Obtain the number of sprouts
		for (int i = 0; i < NSP; i++)									// Deactivate all sprout force components
		{
			FESproutBodyForce::SPROUT& sp = pbf->GetSprout(i);				// Obtain the sprout force component
			sp.active = false;												// Deactive the sprout force component
		}
		FEParameterList& pl = pbf->GetParameterList();										// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);												// Get the sprout force magnitude parameter
		pa->value<double>() = magnitude*m_sproutf;													// Set the sprout force magnitude parameter
	}

	FEMesh& mesh = m_fem.GetMesh();									// Obtain the FE mesh

	//#pragma omp parallel for
	for (SegIter frag_it = m_pCult->SegmentBegin(); frag_it != m_pCult->SegmentEnd(); ++frag_it)		// Iterate through each segment in the model...
	{
		const Segment& seg = (*frag_it);								// Obtain the segment, keep it constant to prevent changes

		if (((seg.tip(0).bactive) || (seg.tip(0).BC == 1)) && (seg.tip(0).bdyf_id >= 0)){		  // Turn on the body force for any active -1 segment OR -1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[0] == -1) && (seg.bdyf_id[0] >= 0)){									// Turn on the body force for any active -1 segment
			tip = seg.tip(0).pos();																	// Obtain the tip position

			sprout_vect = seg.tip(0).pos() - seg.tip(1).pos();												// Calculate the sprout directional vector
			sprout_vect = sprout_vect/sprout_vect.norm();			

			update_angio_sprout(seg.tip(0).bdyf_id, true, tip, sprout_vect);
			}
		
		if (((seg.tip(1).bactive) || (seg.tip(1).BC == 1)) && (seg.tip(1).bdyf_id >= 0)){		  // Turn on the body force for any active +1 segment OR +1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[1] == 1) && (seg.bdyf_id[1] >= 0)){									// Turn on the body force for any active +1 segment
			tip = seg.tip(1).pos();																	// Obtain the tip position
			
			sprout_vect = seg.tip(1).pos() - seg.tip(0).pos();												// Calculate the sprout directional vector
			sprout_vect.unit();
			update_angio_sprout(seg.tip(1).bdyf_id, true, tip, sprout_vect);
			}
	}
	
	return;
}

//-----------------------------------------------------------------------------
void FEAngio::update_angio_sprout(int id, bool bactive, const vec3d& rc, const vec3d& sprout_vect)
{
	if (m_pmat)
	{
		FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(id);
		sp.bactive = true;
		sp.sprout = sprout_vect;
		m_pmat->UpdateSprout(sp, rc);
	}
	else
	{
		FESproutBodyForce::SPROUT& sp = m_pbf->GetSprout(id);		// Obtain the sprout component 
		sp.active = true;											// Set the sprout component to active
		sp.rc = rc;													// Set the tip position
		sp.sprout = sprout_vect;									// Set the sprout force directional vector
	}
}

///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////
// TODO: vessel lengths are always positive now, so we need to fix the logic here.
void FEAngio::adjust_mesh_stiffness()
{
	if (comp_mat == 0)													// If a composite consitutive model isn't being used, exit
		return;

	Grid& grid = m_grid;
		
	int elem_num = 0;													// Element number
	vec3d vess_vect;													// Vessel vector

	int NE = m_grid.Elems();
	for (int i = 0; i < NE; ++i)
	{							
		Elem& el = grid.GetElement(i);
		el.alpha = 0.;						// Set the vessel volume fraction, alpha, to zero
		el.fiber_orient = vec3d(0,0,0);		// Set the vessel orientation vector to 0 (this is the element fiber_orient vector, which does not contain collagen fiber orientation information but rather the material fiber direction for a transversely isotropic material model)
	}
	
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1/(double)Nsub;									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

	for (SegIter frag_it = m_pCult->SegmentBegin(); frag_it != m_pCult->SegmentEnd(); ++frag_it)		// For each segment...
	{
		Segment subunit;												// Segment subdivision placeholder
	
		Segment& seg = (*frag_it);												// Obtain the segment
		
		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length() > 0.){											// If it's a +1 segment...
					subunit.tip(0).pt.r = seg.tip(0).pos();
				}
				else{															// If it's a -1 segment...
					subunit.tip(1).pt.r = seg.tip(1).pos();
//					subunit.m_length = -1.;
				}
		}
			
			// Calculate the subdivision
			if (seg.length() > 0.){										// If it's a +1 segment...			
				subunit.tip(1).pt.r = subunit.tip(0).pos() + seg.uvect()*(sub_scale*seg.length());     
			}
			else{														// If it's a -1 segment...
				subunit.tip(0).pt.r = subunit.tip(1).pos() + seg.uvect()*(sub_scale*seg.length());     
			}

			subunit.Update();										// Find the length of the subdivision
			
			mid = (subunit.tip(1).pos() + subunit.tip(0).pos())*0.5;

			elem_num = m_grid.findelem(mid);				// Find the element that the midpoint is within

			// Calculate the orientation of the subdivision
			if (seg.length() > 0.){										// If it's a +1 segment...
				vess_vect = subunit.tip(1).pos() - subunit.tip(0).pos();
			}
			else{														// If it's a -1 segment...
				vess_vect = subunit.tip(0).pos() - subunit.tip(1).pos();
			}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
				vess_vect = vess_vect/vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				Elem& el = m_grid.GetElement(elem_num);

				elem_volume = el.m_volume;					// Calculate the volume of the element

				subunit_volume = pi*(m_vessel_width/2.)*(m_vessel_width/2.)*fabs(subunit.length());		// Find the volume of the subdivision
				volume_fraction = subunit_volume/elem_volume;				// Calculate the volume fraction

				el.alpha = el.alpha + volume_fraction;	// Add the volume fraction for each subdivision to alpha
			
				// Calculate the vessel orientation vector 
				if ((el.fiber_orient.x == 0) && (el.fiber_orient.y == 0) && (el.fiber_orient.z == 0)){	// If the vessel orientation vector hasn't been assigned yet...
					el.fiber_orient = vess_vect;			// Set the vessel orientation vector					
					el.fiber_orient.z = vess_vect.z;
				}
				else{														// If it has been...	
					el.fiber_orient.x = (el.fiber_orient.x + vess_vect.x)/2;	// Average together the vessel orientation vector
					el.fiber_orient.y = (el.fiber_orient.y + vess_vect.y)/2;
					el.fiber_orient.z = (el.fiber_orient.z + vess_vect.z)/2;}
			}
			
			// Set the origin of the next subdivision to the end of the current one
			if (seg.length() > 0.){
				subunit.tip(0).pt.r = subunit.tip(1).pos();
			}
			else{
				subunit.tip(1).pt.r = subunit.tip(0).pos();
			}
		}
	}
	
	double alpha = 0.;													// Volume fraction for the composite material model
	vec3d e1; vec3d e2; vec3d e3;										// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
	
	FEMesh& mesh = m_fem.GetMesh();										// Get the FE mesh
	int J = mesh.Domains();												// Find the number of domains within the mesh
    int num_elem = 0;

	for (int k = 0; k < J; ++k){
		FEDomain& d = mesh.Domain(k);										// Obtain the domain
		
		for (int j = 0; j < d.Elements(); ++j)								// For each element within the domain...
		{
			FEElement& e = d.ElementRef(j);										// Obtain the element from the domain
			int nint = e.GaussPoints();											// Obtain the number of gauss points
		
			Elem& eg = grid.GetElement(num_elem);
			alpha = eg.alpha;											// Obtain alpha from the grid element

			// Set e1 to the vessel orientation vector
			e1 = eg.fiber_orient;

			if ((e1.x == 0) && (e1.y == 0) && (e1.z == 0)){						// If there is not vessels in the element, set the material basis to the global coordinate basis
				e1 = vec3d(1,0,0);
				e2 = vec3d(0,1,0);
				e3 = vec3d(0,0,1);}
			else{																// Else, set the other two directions to be orthogonal to the vessel orientation
				e2.y = 1;
				e2 = e1^e2;
				e3 = e1^e2;}
		
			for (int n = 0; n < nint; ++n)										// For each gauss point...
			{
				//FEMaterialPoint& mp = *e.m_State[n];								// Obtain the material point
				//FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();	// Obtain the mixture material point
				//pt.m_w[0] = alpha;													// Set the first weight factor to alpha
				//pt.m_w[1] = 1.0 - alpha;											// Set the second weight factor to 1 - alpha
			
				FEMaterialPoint& mp = *(e.GetMaterialPoint(n)->GetPointData(1));    // returns the second component of the mixture
				FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>(); // get the mixture material point
				pt.m_w[0] = alpha;
				pt.m_w[1] = 1.0 - alpha;

				if (comp_mat == 2){													// If the transversely isotropic material is being used...
					FEElasticMaterialPoint& pt2 = *mp.ExtractData<FEElasticMaterialPoint>();
					pt2.m_Q[0][0] = e1.x;												// Set the first column of Q to e1
					pt2.m_Q[1][0] = e1.y;
					pt2.m_Q[2][0] = e1.z;
					pt2.m_Q[0][1] = e2.x;												// Set the second column of Q to e2
					pt2.m_Q[1][1] = e2.y;
					pt2.m_Q[2][1] = e2.z;
					pt2.m_Q[0][2] = e3.x;												// Set the third column of Q to e3
					pt2.m_Q[1][2] = e3.y;
					pt2.m_Q[2][2] = e3.z;}
				}

			num_elem++;
		}
	}

	return;
}

//-----------------------------------------------------------------------------
// Calculate the density gradient for each element (this function may not work properly)
void FEAngio::update_ecm_den_grad()
{
	assert(false);
/*	vec3d elem_den_grad;
	double ex = 0.; double ey = 0.; double ez = 0.;	
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;
	
	dN1 = m_grid.shapefun_d1(ex, ey, ez, 1);
    dN2 = m_grid.shapefun_d1(ex, ey, ez, 2);
	dN3 = m_grid.shapefun_d1(ex, ey, ez, 3);
	dN4 = m_grid.shapefun_d1(ex, ey, ez, 4);
	dN5 = m_grid.shapefun_d1(ex, ey, ez, 5);
	dN6 = m_grid.shapefun_d1(ex, ey, ez, 6);
	dN7 = m_grid.shapefun_d1(ex, ey, ez, 7);
	dN8 = m_grid.shapefun_d1(ex, ey, ez, 8);
	
	Grid& grid = m_grid;
	int NN = grid.Nodes();
	for (int j = 0; j < NN; ++j){
		grid.nodes[j].updated = false;}


	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){
		Elem& elem = grid.ebin[i];
		
		elem_den_grad.x = ((*elem.n1).ecm_den)*dN1.x + ((*elem.n2).ecm_den)*dN2.x + ((*elem.n3).ecm_den)*dN3.x + ((*elem.n4).ecm_den)*dN4.x + ((*elem.n5).ecm_den)*dN5.x + ((*elem.n6).ecm_den)*dN6.x + ((*elem.n7).ecm_den)*dN7.x + ((*elem.n8).ecm_den)*dN8.x;
		elem_den_grad.y = ((*elem.n1).ecm_den)*dN1.y + ((*elem.n2).ecm_den)*dN2.y + ((*elem.n3).ecm_den)*dN3.y + ((*elem.n4).ecm_den)*dN4.y + ((*elem.n5).ecm_den)*dN5.y + ((*elem.n6).ecm_den)*dN6.y + ((*elem.n7).ecm_den)*dN7.y + ((*elem.n8).ecm_den)*dN8.y;
		elem_den_grad.z = ((*elem.n1).ecm_den)*dN1.z + ((*elem.n2).ecm_den)*dN2.z + ((*elem.n3).ecm_den)*dN3.z + ((*elem.n4).ecm_den)*dN4.z + ((*elem.n5).ecm_den)*dN5.z + ((*elem.n6).ecm_den)*dN6.z + ((*elem.n7).ecm_den)*dN7.z + ((*elem.n8).ecm_den)*dN8.z;
	
		if ((*elem.n1).updated == false){
			(*elem.n1).ecm_den_grad = elem_den_grad;
			(*elem.n1).updated = true;}
		else if ((*elem.n1).updated == true){
			(*elem.n1).ecm_den_grad = (elem_den_grad + (*elem.n1).ecm_den_grad)*0.5;}

		if ((*elem.n2).updated == false){
			(*elem.n2).ecm_den_grad = elem_den_grad;
			(*elem.n2).updated = true;}
		else if ((*elem.n2).updated == true){
			(*elem.n2).ecm_den_grad = (elem_den_grad + (*elem.n2).ecm_den_grad)*0.5;}
		
		if ((*elem.n3).updated == false){
			(*elem.n3).ecm_den_grad = elem_den_grad;
			(*elem.n3).updated = true;}
		else if ((*elem.n3).updated == true){
			(*elem.n3).ecm_den_grad = (elem_den_grad + (*elem.n3).ecm_den_grad)*0.5;}

		if ((*elem.n4).updated == false){
			(*elem.n4).ecm_den_grad = elem_den_grad;
			(*elem.n4).updated = true;}
		else if ((*elem.n4).updated == true){
			(*elem.n4).ecm_den_grad = (elem_den_grad + (*elem.n4).ecm_den_grad)*0.5;}

		if ((*elem.n5).updated == false){
			(*elem.n5).ecm_den_grad = elem_den_grad;
			(*elem.n5).updated = true;}
		else if ((*elem.n5).updated == true){
			(*elem.n5).ecm_den_grad = (elem_den_grad + (*elem.n5).ecm_den_grad)*0.5;}

		if ((*elem.n6).updated == false){
			(*elem.n6).ecm_den_grad = elem_den_grad;
			(*elem.n6).updated = true;}
		else if ((*elem.n6).updated == true){
			(*elem.n6).ecm_den_grad = (elem_den_grad + (*elem.n6).ecm_den_grad)*0.5;}

		if ((*elem.n7).updated == false){
			(*elem.n7).ecm_den_grad = elem_den_grad;
			(*elem.n7).updated = true;}
		else if ((*elem.n7).updated == true){
			(*elem.n7).ecm_den_grad = (elem_den_grad + (*elem.n7).ecm_den_grad)*0.5;}

		if ((*elem.n8).updated == false){
			(*elem.n8).ecm_den_grad = elem_den_grad;
			(*elem.n8).updated = true;}
		else if ((*elem.n8).updated == true){
			(*elem.n8).ecm_den_grad = (elem_den_grad + (*elem.n8).ecm_den_grad)*0.5;}
	}
*/
}


//-----------------------------------------------------------------------------
void FEAngio::update_sprout_stress_scaling()
{
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	if (m_pmat)
		m_pmat->scale = y0 + a/(1 + exp(-(m_time.t - x0)/b));
	
	return;
}

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
	// get the material's coordinate system
	FECoordSysMap* pmap = pmat->GetCoordinateSystemMap();

	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the fiber array
	int N = mesh.Nodes();
	fiber.assign(N, vec3d(0,0,0));

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local fiber orientation at integration points
			vector<double> fx(nint), fy(nint), fz(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mpoint->ExtractData<FEElasticMaterialPoint>();
				mat3d m = pt.m_Q;
				
				// grab the first column as the fiber orientation
				fx[n] = m[0][0];
				fy[n] = m[1][0];
				fz[n] = m[2][0];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln), gy(neln), gz(neln);
			el.project_to_nodes(&fx[0], &gx[0]);
			el.project_to_nodes(&fy[0], &gy[0]);
			el.project_to_nodes(&fz[0], &gz[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				vec3d f(gx[i], gy[i], gz[i]);
				fiber[ni] += f;
			}
		}
	}

	// normalize the fibers
	for (int i=0; i<N; ++i) fiber[i].unit();

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateDensityMap(vector<double>& density, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	density.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> den(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMixtureMaterialPoint* mPt = dynamic_cast<FEElasticMixtureMaterialPoint*>(mpoint);

				if(mPt)
				{
					vector<FEMaterialPoint*> mPtV = mPt->m_mp;
					for (int i=0; i<(int)mPtV.size(); ++i)
					{
						FEAngioMaterialPoint* angioPt = dynamic_cast<FEAngioMaterialPoint*>(mPtV[i]);
						if(angioPt)
						{
							den[n] = angioPt->m_D;
							break;
						}
					}
				}
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			el.project_to_nodes(&den[0], &gx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				density[ni] = gx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}
