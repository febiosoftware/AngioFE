#pragma once
#include "StdAfx.h"
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include <FECore/FESolidDomain.h>
#include "AngioElement.h"
#include "FEProbabilityDistribution.h"
#include <FECore/FEOctreeSearch.h>

//! Base class for all fragment seeders
class FragmentSeeder : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FragmentSeeder)
public:
	//! seed the fragments within the elements
	virtual bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) = 0;
	//! constructor for class
	explicit FragmentSeeder(FEModel* model);
	virtual ~FragmentSeeder() {}
	bool Init() override {
		return initial_segment_length->Init();
	}
	//! Get a random set of valid natural coordinates
	vec3d GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(FEMesh* mesh, AngioElement* angio_element, angiofe_random_engine random_engine);
	//! return an incremented cell id
	int IncrementCellCounter();
	bool proto_mat_cross = true;
protected:
	//! number of fragments to seed
	int number_fragments = 0;
	//! used to initialize fragment ids
	static int initial_fragment_id_counter;
	//! number of cells
	int number_cells = 0;
	//! used to initialize cell ids
	static int initial_cell_id_counter;
	double cell_radius = 10e-6;
	//! Probability Distribution
	FEProbabilityDistribution* initial_segment_length = nullptr;
	DECLARE_FECORE_CLASS()
};

//! seed fragments with all elements having equal probabilities of being choosen
class ByElementFragmentSeeder : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByElementFragmentSeeder(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

//! seed fragments with all elements having equal probabilities of being choosen, seeds 2 oppositly directed tips at each loaction that is choosen
class ByElementFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByElementFragmentSeederBiDirectional(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

//! seed fragments with all elements having equal probabilities of being choosen, seeds 2 oppositly directed tips at each loaction that is choosen
class ByElementSetFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByElementSetFragmentSeederBiDirectional(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

//! seed fragments with each unit of volume having the same probability of being choosen
class ByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByVolumeFragmentSeeder(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

//! seed fragments with each unit of volume having the same probability of being choosen, seeds 2 oppositly directed tips at each loaction that is choosen
class ByVolumeFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByVolumeFragmentSeederBiDirectional(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

//! seed fragments with all elements having equal probabilities of being choosen
class SingleCellSeeder : public FragmentSeeder
{
public:
	//! constructor for class
	explicit SingleCellSeeder(FEModel* model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
	vec3d initial_position = vec3d(0, 0, 0);
protected:
	DECLARE_FECORE_CLASS()
private:
	FEOctreeSearch		m_search;
	FESolidElement* m_el;
	double				m_q[3];
};