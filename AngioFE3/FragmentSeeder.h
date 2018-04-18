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
	virtual bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEAngioMaterial* angio_mat) = 0;
	explicit FragmentSeeder(FEModel * model);
	virtual ~FragmentSeeder() {}
protected:
	int number_fragments = 0;
	double initial_vessel_length = 20.0;
	DECLARE_PARAMETER_LIST();
};

class ByElementFragmentSeeder :public FragmentSeeder
{
public:
	explicit ByElementFragmentSeeder(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEAngioMaterial* angio_mat)override;
};

class ByVolumeFragmentSeeder :public FragmentSeeder
{
public:
	explicit ByVolumeFragmentSeeder(FEModel * model);
	bool SeedFragments(std::vector<AngioElement *> &angio_elements, FEAngioMaterial* angio_mat)override;
};