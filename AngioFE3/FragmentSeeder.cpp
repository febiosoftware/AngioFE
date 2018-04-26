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

ByElementFragmentSeeder::ByElementFragmentSeeder(FEModel * model) :FragmentSeeder(model)
{
	
}

bool ByElementFragmentSeeder::SeedFragments(std::vector<AngioElement *> &angio_elements, FEAngioMaterial* angio_mat, int buffer_index)
{
	FEMesh * mesh = angio_mat->m_pangio->GetMesh();
	std::uniform_int_distribution<int> edist(0, angio_elements.size()-1);

	if(angio_elements.size() == 0)
	{
		return false;
	}

	for (int i = 0; i < number_fragments; ++i)
	{
		Tip* t = new Tip();
		t->local_pos = angio_mat->m_pangio->uniformInUnitCube();
		int elem_index = edist(angio_mat->m_pangio->rengine);
		t->angio_element = angio_elements[elem_index];
		t->time = 0;
		t->is_branch = true;
		t->direction = angio_mat->m_pangio->uniformRandomDirection();
		t->face = t->angio_element;
		//finally add this to the AngioElement
		t->angio_element->next_tips[t->angio_element].push_back(t);
	}

	return true;
}