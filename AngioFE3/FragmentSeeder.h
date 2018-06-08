#pragma once
#include "StdAfx.h"
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include "CultureParameters.h"
#include <FECore/FESolidDomain.h>
#include "AngioElement.h"

class FragmentSeeder : public FEMaterial
{
public:
	virtual bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) = 0;
	explicit FragmentSeeder(FEModel * model);
	virtual ~FragmentSeeder() {}
	vec3d GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(FEMesh* mesh, AngioElement* angio_element, angiofe_random_engine random_engine);
protected:
	int number_fragments = 0;
	static int initial_fragment_id_counter;
	double initial_vessel_length = 20.0;
	DECLARE_PARAMETER_LIST();
};

class ByElementFragmentSeeder : public FragmentSeeder
{
public:
	explicit ByElementFragmentSeeder(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};
class ByElementFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	explicit ByElementFragmentSeederBiDirectional(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

class ByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	explicit ByVolumeFragmentSeeder(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};

class ByVolumeFragmentSeederBiDirectional : public FragmentSeeder
{
public:
	explicit ByVolumeFragmentSeederBiDirectional(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index) override;
};