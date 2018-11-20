#pragma once
#include "StdAfx.h"
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include <FECore/FESolidDomain.h>
#include "AngioElement.h"

class FragmentSeeder : public FEMaterial
{
public:
	//! seed the fragments within the elements
	virtual bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) = 0;
	//! constructor for class
	explicit FragmentSeeder(FEModel * model);
	virtual ~FragmentSeeder() {}
	//! Get a random set of valid natural coordinates
	vec3d GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(FEMesh* mesh, AngioElement* angio_element, angiofe_random_engine random_engine);
protected:
	//! number of fragments to seed
	int number_fragments = 0;
	//! used to initialize fragment ids
	static int initial_fragment_id_counter;
	//! initial vessel length
	double initial_vessel_length = 20.0;
	//! parameter list
	DECLARE_PARAMETER_LIST();
};

class ByElementFragmentSeeder : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByElementFragmentSeeder(FEModel * model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};
class ByElementFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByElementFragmentSeederBiDirectional(FEModel * model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

class ByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByVolumeFragmentSeeder(FEModel * model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

class ByVolumeFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	//! constructor for class
	explicit ByVolumeFragmentSeederBiDirectional(FEModel * model);
	//! Seed the fragments
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};