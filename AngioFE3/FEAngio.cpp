///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEDataLoadCurve.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMix/FESolute.h"
#include "FEBioMix/FEMultiphasic.h"
#include "FEBioLib/FEBioModel.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FECore/FESolidDomain.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FECore/FEElemElemList.h"
#include "angio3d.h"
#include <ctime>
#include <future>
#include <algorithm>
#include <cfloat> 
#include <FEBioMix/FETriphasic.h>
#include "GrowDirectionModifier.h"
#include <omp.h>


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat);

bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID);

void FEAngio::SetSeeds()
{
	for(int i=0; i <angio_elements.size();i++)
	{
		angio_elements[i]->_angio_mat->SetSeeds(angio_elements[i]);
	}
}

bool FEAngio::SeedFragments()
{
	bool rv = true;
	for(auto iter = elements_by_material.begin(); iter != elements_by_material.end(); ++iter)
	{
		rv &= iter->first->SeedFragments(iter->second);
	}
	return rv;
}

void FEAngio::ApplydtToTimestepper(double dt)
{
	FETimeStepController &tsc = m_fem->GetCurrentStep()->m_timeController;
	if(tsc.m_nmplc == -1)
	{
		//there is no load curve being used, just set dtmax
		tsc.m_dtmax = dt;
	}
	else
	{
		//not implemented yet
		//use the lower of the loadcurve vs dt
		assert(false);
	}
}

void FEAngio::GrowSegments()
{
	auto time_info = m_fem->GetTime();
	double min_dt = std::numeric_limits<double>::max();
	size_t angio_element_count = angio_elements.size();
	//is it valid for this to run more than once?
	if(time_info.currentTime >= next_time)
	{
		
		
		#pragma omp parallel for shared(min_dt)
		for (int i = 0; i < angio_element_count; ++i)
		{
			double temp_dt = angio_elements[i]->_angio_mat->GetMin_dt(angio_elements[i]);
			#pragma omp critical
			{
				min_dt = std::min(min_dt, temp_dt);
			}
		}
		
		double ctime = time_info.currentTime + min_dt;

		#pragma omp parallel for schedule(dynamic, 16)
		for (int j = 0; j <angio_element_count; j++)
		{
			angio_elements[j]->_angio_mat->PrepBuffers(angio_elements[j], next_time, buffer_index);
		}

//#pragma omp parallel
		{
			int n = 3;
			for (int i = 0; i <n; i++)
			{
//#pragma omp  parallel for schedule(dynamic, 16)
				for (int j = 0; j <angio_element_count; j++)
				{
					angio_elements[j]->_angio_mat->GrowSegments(angio_elements[j], ctime, buffer_index);
				}

#pragma omp parallel for schedule(dynamic, 16)
				for (int j = 0; j <angio_element_count; j++)
				{
					angio_elements[j]->_angio_mat->PostGrowthUpdate(angio_elements[j], ctime, buffer_index);
				}
				buffer_index = (buffer_index + 1) % 2;
			}
		}

		ApplydtToTimestepper(min_dt);
	}
	//do the output
	fileout->save_vessel_state(*this);

	//do the cleanup if needed
	if(time_info.currentTime >= next_time)
	{
		next_time += min_dt;
		printf("\nangio dt chosen is: %lg next angio time\n", min_dt);
		#pragma omp parallel for schedule(dynamic, 16)
		for (int j = 0; j <angio_element_count; j++)
		{
			angio_elements[j]->_angio_mat->Cleanup(angio_elements[j], next_time, buffer_index);
		}
		buffer_index = (buffer_index + 1) % 2;
	}

}

vec3d FEAngio::ReferenceCoordinates(Tip * tip) const
{
	vec3d r(0, 0, 0);
	FEMesh & mesh = m_fem->GetMesh();
	//Point has already been positioned
	FESolidElement * se= tip->angio_element->_elem;

	double arr[FEElement::MAX_NODES];
	se->shape_fnc(arr, tip->local_pos.x, tip->local_pos.y, tip->local_pos.z);
	for (int j = 0; j < se->Nodes(); j++)
	{
		r += mesh.Node(se->m_node[j]).m_r0* arr[j];
	}

	return r;
}

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : ztopi(std::uniform_real_distribution<double>(0, PI)), 
	zto2pi(std::uniform_real_distribution<double>(0, 2 * PI)), n1to1(std::uniform_real_distribution<double>(-1, 1))
{
	// Body force counter
	total_bdyf = 0;
	
	FE_state = 0;

	// initialize time stepping parameters
	m_time.dt = 1.0;
	// Input random seed number
	m_irseed = 0;

	m_fem = dynamic_cast<FEBioModel *>(&fem);
	assert(m_fem);
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
}

//-----------------------------------------------------------------------------
FEBioModel* FEAngio::GetFEModel() const
{
	return m_fem;
}

FEMesh * FEAngio::GetMesh() const
{
	return &m_fem->GetMesh();
}

//consider moving this to FEBio
std::pair<int, int> FEAngio::GetMinMaxElementIDs() const
{
	FEMesh * mesh = GetMesh();
	assert(mesh);

	// get the ID ranges
	int m_minID = -1;
	int m_maxID = -1;
	int NDOM = mesh->Domains();
	for (int i = 0; i<NDOM; ++i)
	{
		FEDomain& dom = mesh->Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j<NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			int eid = el.GetID();
			if ((eid < m_minID) || (m_minID == -1)) m_minID = eid;
			if ((eid > m_maxID) || (m_maxID == -1)) m_maxID = eid;
		}
	}
	return { m_minID, m_maxID };
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	//create any classes which have nontrivial destructors 

	// Init all the FE stuff
	//must be done first initializes material
	if (InitFEM() == false) return false;
	
	SetupAngioElements();

	SetSeeds();


	FinalizeFEM();

	//do the seeding of the tips/segments
	if(!SeedFragments())
	{
		printf("fragment seeding failed\n");
		return false;
	}

	// start timer
	time(&m_start);

	return true;
}


//-----------------------------------------------------------------------------
// Initialize FE model.
bool FEAngio::InitFEM()
{
	for (int i = 0; i < m_fem->Materials(); i++)
	{
		FEAngioMaterial * cmat = nullptr;
		FEMaterial * mat = m_fem->GetMaterial(i);
		int id = -1;
		cmat = GetAngioComponent(mat);
		

		if (cmat)
		{
			id = mat->GetID();
			assert(id != -1);
			m_pmat.emplace_back(cmat);
			m_pmat_ids.emplace_back(id);
			//TODO: check that material parameters are set here
			//cmat->ApplySym();
			cmat->SetFEAngio(this);
		}
	}
	//assert(m_pmat.size());

	felog.printf("%d Angio materials found. Stress approach will be used.", m_pmat.size());

	// register the callback
	m_fem->AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED | CB_STEP_ACTIVE, this);

	// Do the model initialization
	if (m_fem->Init() == false) return false;

	return true;
}
void FEAngio::FinalizeFEM()
{

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::LOG_FILE);

	//currently the destructors are not called for classes created by FEBio this allows destructors to be called
	fileout = new Fileout(*this);

	// --- Output initial state of model ---
	if (!m_fem->GetGlobalConstant("no_io"))
	{
		// Output initial microvessel state
		fileout->save_vessel_state(*this);

		// save active tips
		fileout->save_active_tips(*this);
	}
}
//fills in the adjacnecy information 
void FEAngio::FillInAdjacencyInfo(FEMesh * mesh,FEElemElemList * eel, AngioElement *angio_element, int elem_index)
{
	//one element of adjaceny per face
	FESolidElement * se = angio_element->_elem;
	assert(se);
	for(int i=0; i < mesh->Faces(*se);i++)
	{
		FEElement * elem = eel->Neighbor(elem_index, i);
		FESolidElement * nse = dynamic_cast<FESolidElement*>(elem);
		if(nse)
		{
			//next look to see if there is an angio element that uses this element
			auto ae_iter = se_to_angio_elem.find(nse);
			if(ae_iter != se_to_angio_elem.end())
			{
				angio_element->adjacency_list.push_back(ae_iter->second.first);
			}
			else
			{
				assert(false);
			}
		}
		else
		{
			angio_element->adjacency_list.push_back(nullptr);
		}
	}
}

void FillInFaces(FEMesh * mesh, AngioElement * angio_element)
{
	//for hex elements the faces produced have the normals pointing outwards with a CCW winding ... verify this is true for other elements
	int nodes[FEElement::MAX_NODES];
	int faces = mesh->Faces(*angio_element->_elem);
	FEFacetSet fs(mesh);
	fs.Create(faces);
	for(int i=0; i < faces;i++)
	{
		FEFacetSet::FACET& face = fs.Face(i);
		int face_nodes = mesh->GetFace(*angio_element->_elem, i, nodes);
		face.ntype = face_nodes;
		//do any remapping of nodes here
		//might do most of the remappings .. will fail if there is a central node in a face
		assert(face_nodes>=3);
		std::vector<int> v_nodes;
		for(int j=0; j < face_nodes;j++)
		{
			v_nodes.push_back(nodes[j]);
		}
		std::reverse(v_nodes.begin() + 1, v_nodes.end());
		for(int j=0; j < v_nodes.size();j++)
		{
			face.node[j] = v_nodes[j];
		}

	}
	angio_element->inner_faces.BuildFromSet(fs);
}

void FEAngio::SetupAngioElements()
{
	//setup the element access data structure
	auto min_max = GetMinMaxElementIDs();
	auto min_elementID = std::get<0>(min_max);
	//initiialize the FEAngioElementData Structures

	FEMesh * mesh = GetMesh();
	FEModel * model = GetFEModel();
	//this stores the element adjacency information
	FEElemElemList eel;
	eel.Create(mesh);

	//we are using the lut in FEMesh so the speed should be good
	//the table now contains only angio elements
	for (int i = 0; i < mesh->Elements(); i++)
	{
		FESolidElement * elem = dynamic_cast<FESolidElement*>(mesh->FindElementFromID(min_elementID + i));
		if (elem)
		{
			FEMaterial * mat = model->GetMaterial(elem->GetMatID());
			FEAngioMaterial*angio_mat = GetAngioComponent(mat);
			AngioElement * angio_element = nullptr;
			if (angio_mat)
			{
				angio_element = new AngioElement(elem, angio_mat, mat, mesh);
				angio_elements.push_back(angio_element);
				if (elements_by_material.find(angio_mat) == elements_by_material.end())
				{
					std::vector<AngioElement *> ang_elem;
					elements_by_material[angio_mat] = ang_elem;
				}
				elements_by_material[angio_mat].push_back(angio_element);
				se_to_angio_elem[elem] = { angio_element,i };
			}
			//now fill in the version with holes
			angio_elements_with_holes.push_back(angio_element);
		}
	}
	for(int i=0; i < angio_elements.size();i++)
	{
		FillInAdjacencyInfo(mesh, &eel,  angio_elements[i], se_to_angio_elem[angio_elements[i]->_elem].second);
	}
	for (int i = 0; i < angio_elements.size(); i++)
	{
		FillInFaces(mesh, angio_elements[i]);
	}
}


double FEAngio::GetDoubleFromDataStore(int record, int elem_id, int item)
{
	DataStore & ds = m_fem->GetDataStore();
	return ds.GetDataRecord(record)->Evaluate(elem_id, item);
}

mat3d FEAngio::unifromRandomRotationMatrix()
{
	//collagen fibers are right handed so the following transformation is legal
	//the following will only produce right handed bases for the collagen fibers which is molecularly accurate
	double alpha = zto2pi(rengine);
	double beta = zto2pi(rengine);
	double gamma = zto2pi(rengine);


	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, - s_alpha*s_gamma+c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
		);
	//rv = rv.transinv();
#ifndef NDEBUG
	//do some testing that bases run tthrough this are still orthogonal
	vec3d x(1, 0, 0);
	vec3d y(0, 1, 0);
	vec3d z(0, 0, 1);
	vec3d xt = rv * x;
	vec3d yt = rv * y;
	vec3d zt = rv * z;
	double tol = 0.01;
	assert(xt * yt < tol && xt * yt > -tol);
	assert(xt * zt < tol && xt * zt > -tol);
	assert(zt * yt < tol && zt * yt > -tol);

	double noq = xt.norm();
	double nom1 = yt.norm();
	double nom2 = zt.norm();
	tol = 0.01;
	assert(noq < (1 + tol) && noq >(1 - tol));
	assert(nom1 < (1 + tol) && nom1 >(1 - tol));
	assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif
	return rv;
}
mat3d FEAngio::rotationMatrix(double alpha, double beta, double gamma)
{
	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, -s_alpha*s_gamma + c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
	);
	return rv;
}

vec3d FEAngio::uniformRandomDirection()
{
	//to revert this set this to return vrand
	double theta = zto2pi(rengine);
	double phi = zto2pi(rengine);
	double sintheta = sin(theta);
	vec3d dir(sintheta * cos(phi), sintheta * sin(phi), cos(theta));
	return dir;
}
vec3d FEAngio::uniformInUnitCube()
{
	//force the order of initialization should allow results on different platforms to converge 
	double x = n1to1(rengine);
	double y = n1to1(rengine);
	double z = n1to1(rengine);
	return vec3d(x,y,z);
}

vec3d FEAngio::gradient(FESolidElement * se, std::vector<double> & fn, vec3d pt, int size, int offset)
{
	assert(se);
	FESolidDomain* domain = dynamic_cast<FESolidDomain*>(se->GetDomain());
	double Ji[3][3];
	domain->invjact(*se, Ji, pt.x, pt.y, pt.z);
	double Gr[FEElement::MAX_NODES], Gs[FEElement::MAX_NODES], Gt[FEElement::MAX_NODES];
	se->shape_deriv(Gr, Gs, Gt, pt.x, pt.y, pt.z);

	vec3d gradf;
	size_t N = se->Nodes();
	assert(N == size*fn.size());
	for (size_t i = 0; i<N; ++i)
	{
		double Gx, Gy, Gz;
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
		Gy = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
		Gz = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];

		// calculate pressure gradient
		gradf.x += Gx*fn[i*size + offset];
		gradf.y += Gy*fn[i*size + offset];
		gradf.z += Gz*fn[i*size + offset];
	}

	return gradf;
}
//-----------------------------------------------------------------------------
static void solve_3x3(double A[3][3], double b[3], double x[3])
{
	double D = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[2][1] * A[0][2] \
		- A[1][1] * A[2][0] * A[0][2] - A[2][2] * A[1][0] * A[0][1] - A[0][0] * A[2][1] * A[1][2];

	assert(D != 0);

	double Ai[3][3];
	Ai[0][0] = A[1][1] * A[2][2] - A[2][1] * A[1][2];
	Ai[0][1] = A[2][1] * A[0][2] - A[0][1] * A[2][2];
	Ai[0][2] = A[0][1] * A[1][2] - A[1][1] * A[0][2];

	Ai[1][0] = A[2][0] * A[1][2] - A[1][0] * A[2][2];
	Ai[1][1] = A[0][0] * A[2][2] - A[2][0] * A[0][2];
	Ai[1][2] = A[1][0] * A[0][2] - A[0][0] * A[1][2];

	Ai[2][0] = A[1][0] * A[2][1] - A[2][0] * A[1][1];
	Ai[2][1] = A[2][0] * A[0][1] - A[0][0] * A[2][1];
	Ai[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	x[0] = (Ai[0][0] * b[0] + Ai[0][1] * b[1] + Ai[0][2] * b[2]) / D;
	x[1] = (Ai[1][0] * b[0] + Ai[1][1] * b[1] + Ai[1][2] * b[2]) / D;
	x[2] = (Ai[2][0] * b[0] + Ai[2][1] * b[1] + Ai[2][2] * b[2]) / D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]);
	r[1] = b[1] - (A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]);
	r[2] = b[2] - (A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]);

	double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
#endif
}





FEAngioMaterial * FEAngio::GetAngioComponent(FEMaterial * mat)
{
	auto angm = dynamic_cast<FEAngioMaterial*>(mat);
	if (!angm)
	{
		FEMultiphasic * mmat = dynamic_cast<FEMultiphasic*>(mat);
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	if (!angm)
	{
		FETriphasic * mmat = dynamic_cast<FETriphasic*>(mat);
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	/* //Biphasic Material doesn't have a get solid function 
	 * //this fucntion should be moved to an interface in febio
	if (!angm)
	{
		FEBiphasic * mmat = dynamic_cast<FEBiphasic*>(mat);
		angm = dynamic_cast<FEAngioMaterialBase*>(mmat->GetSolid());
	}
	*/
	return angm;
}

vec3d FEAngio::Position(FESolidElement * se, vec3d local) const
{
	double arr[FEElement::MAX_NODES];
	assert(se);
	se->shape_fnc(arr, local.x, local.y, local.z);
	vec3d rc(0, 0, 0);

	auto mesh = GetMesh();
	for (int j = 0; j < se->Nodes(); j++)
	{
		rc += mesh->Node(se->m_node[j]).m_rt* arr[j];
	}
	return rc;
}

//-----------------------------------------------------------------------------
void FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;
	FETimeInfo fti = fem.GetTime();
	m_time.t = fti.currentTime;
	m_time.dt = fti.timeIncrement;
	
	static bool start = false;
	if (!start)
	{
		grow_timer.start();
		GrowSegments();
		grow_timer.stop();
		
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			//expermentally calculate the stress from the vessels here
			m_pmat[i]->UpdateAngioStresses();
		}
		start = true;
	}

	if (nwhen == CB_UPDATE_TIME)
	{
		static int index = 0;
		// grab the time information
		
		

		update_gdms_timer.start();
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->UpdateGDMs();
		}
		update_gdms_timer.stop();


		//new function to find the start time grow time and if this is the final iteration this timestep
		grow_timer.start();
		GrowSegments();
		grow_timer.stop();
		
		mesh_stiffness_timer.start();
		mesh_stiffness_timer.stop();

		update_sprout_stress_scaling_timer.start();

		update_sprout_stress_scaling_timer.stop();

		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			//expermentally calculate the stress from the vessels here
			m_pmat[i]->UpdateAngioStresses();
		}
	}
	else if (nwhen == CB_MAJOR_ITERS)
	{
		// update the grid data
		update_ecm_timer.start();
		update_ecm_timer.stop();

		material_update_timer.start();
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			
		}
		material_update_timer.stop();
		fileout->save_feangio_stats(*this);
		ResetTimers();
		++FE_state;
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			// Save the current vessel state
			fileout->save_vessel_state(*this);

			// save active tips
			fileout->save_active_tips(*this);

			// Print the status of angio3d to the user    
			fileout->printStatus(*this);
		}
		
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			Output();
		}
		//force any destructors to be called that need it
		delete fileout;
		fileout = nullptr;
			
	}
}

void FEAngio::ResetTimers()
{
	grow_timer.reset();
	mesh_stiffness_timer.reset();
	update_sprout_stress_scaling_timer.reset();
	update_gdms_timer.reset();
	update_ecm_timer.reset();
	material_update_timer.reset();
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{	
	//write out the timeline of branchpoints
	fileout->save_timeline(*this);
	//fileout->save_winfiber(*this);
	fileout->save_final_vessel_csv(*this);
}


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
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
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	density.resize(N, 0.0);
	anisotropy.resize(N, 0.0);

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
			vector<double> anis(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			vector<double> fx(neln);
			el.project_to_nodes(&den[0], &gx[0]);
			el.project_to_nodes(&anis[0], &fx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				density[ni] = gx[i];
				if(fx[i] > anisotropy[ni])
					anisotropy[ni] = fx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	concentration.resize(N, 0.0);

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
			vector<double> con(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FESolutesMaterialPoint& spt = *mpoint->ExtractData<FESolutesMaterialPoint>();
				con[n]=spt.m_c[vegfID];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			el.project_to_nodes(&con[0], &gx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				concentration[ni] = gx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}