#include "StdAfx.h"
#include "FragmentSeeder.h"
#include "FEAngioMaterial.h"

FragmentSeeder::FragmentSeeder(FEModel * model) : FEMaterial(model)
{

}

BEGIN_PARAMETER_LIST(FragmentSeeder, FEMaterial)
ADD_PARAMETER(number_fragments, FE_PARAM_INT, "number_fragments");
ADD_PARAMETER(initial_vessel_length, FE_PARAM_DOUBLE, "initial_vessel_length");
END_PARAMETER_LIST();

bool ByElementFragmentSeeder::SeedFragments(std::vector<AngioElement *> &angio_elements, FEAngioMaterial* angio_mat)
{
	FEMesh * mesh = angio_mat->m_pangio->GetMesh();
	std::uniform_int_distribution<int> edist(0, angio_elements.size());

	if(angio_elements.size() == 0)
	{
		return false;
	}

	for (int i = 0; i < number_fragments; ++i)
	{
		Tip t;
		t.local_pos = angio_mat->m_pangio->uniformInUnitCube();
		t.angio_element = angio_elements[edist(angio_mat->m_pangio->rengine)];
		t.time = 0;
		t.is_branch = true;
		t.direction = angio_mat->m_pangio->uniformRandomDirection();
	}

	return true;
}