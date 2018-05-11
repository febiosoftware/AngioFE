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
		Tip* r0 = new Tip();
		int elem_index = edist(angio_mat->m_pangio->rengine);
		r0->angio_element = angio_elements[elem_index];
		std::uniform_real_distribution<double> dist(FEAngio::NaturalCoordinatesLowerBound(r0->angio_element->_elem->Type()), FEAngio::NaturalCoordinatesUpperBound(r0->angio_element->_elem->Type()));
		r0->local_pos =  vec3d(dist(angio_mat->m_pangio->rengine), dist(angio_mat->m_pangio->rengine), dist(angio_mat->m_pangio->rengine));
		r0->time = 0;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection();
		r0->face = r0->angio_element;
		//finally add this to the AngioElement

		//repeat the step if not placed correctly
		if(r0->angio_element->_angio_mat->ProtoGrowthInElement(r0, -1, initial_vessel_length))
		{
			i--;
			continue;
		}
	}

	return true;
}