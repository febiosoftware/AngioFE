#include "StdAfx.h"
#include "FragmentSeeder.h"
#include "FEAngioMaterial.h"
#include "angio3d.h"
#include "FEProbabilityDistribution.h"
#include <iostream>
#include "Tip.h"
#include "FECell.h"

// Variable used to ensure that all initial Tips that are part of the same fragment
// share an ID so that they do not anastamose with each other.
int FragmentSeeder::initial_fragment_id_counter = 0;
int FragmentSeeder::initial_cell_id_counter = 0;

FragmentSeeder::FragmentSeeder(FEModel * model) : FEMaterial(model)
{
	AddClassProperty(this, &initial_segment_length, "initial_segment_length");
}

vec3d FragmentSeeder::GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(FEMesh* mesh, AngioElement* angio_element, angiofe_random_engine random_engine)
{
	int element_type = angio_element->_elem->Type();

	switch(element_type)
	{	
	case FE_Element_Type::FE_TET4G4:
	case FE_Element_Type::FE_TET4G1:
	case FE_Element_Type::FE_TET10G1:
	case FE_Element_Type::FE_TET10G4:
	case FE_Element_Type::FE_TET10G4RI1:
	case FE_Element_Type::FE_TET10G8:
	case FE_Element_Type::FE_TET10G8RI4:
	case FE_Element_Type::FE_TET10GL11:
		{
			std::uniform_real_distribution<double> udv(0, 1);

			// Get a random point within the "far face" (defined in FragmentSeeder.h) of the element.
			std::uniform_real_distribution<double> side_dist(0, 1); // Natural coordinate side length is 1.

			// Sampled scale factor.
			double sf_1 = side_dist(random_engine);
			double sf_2 = side_dist(random_engine);

			// The two sides of the face that form the parllelogram that will be sampled from.
			vec3d s1(vec3d(1, 0, 0) - vec3d(0, 0, 1));
			vec3d s2(vec3d(0, 1, 0) - vec3d(0, 0, 1));

			// Unitize and scale the sides by the scale factor.
			s1.unit();
			s2.unit();
			vec3d s1_s = (s1 * sf_1);
			vec3d s2_s = (s2 * sf_2);

			// Calculate the sample point.
			vec3d sample_pt = vec3d(0, 0, 1) + s1_s + s2_s;

			// The point is in the parallelogram, but not the half that we want.
			// Use the same scale factors, but apply them from a point on the opposite
			// side of the parellelogram so that the result ends up in the correct face.
			if (sample_pt.z < 0)
			{
				sample_pt = vec3d(1, 1, -1) - s1_s - s2_s;
			}

			double volume_sample = udv(random_engine);
			double scale_factor = cbrt(volume_sample);

			return sample_pt * scale_factor;
		}
	default:
		// HEX and PENTA element Types
		std::uniform_real_distribution<double> r_dist(FEAngio::NaturalCoordinatesLowerBound_r(element_type), FEAngio::NaturalCoordinatesUpperBound_r(element_type));
		std::uniform_real_distribution<double> s_dist(FEAngio::NaturalCoordinatesLowerBound_s(element_type), FEAngio::NaturalCoordinatesUpperBound_s(element_type));
		std::uniform_real_distribution<double> t_dist(FEAngio::NaturalCoordinatesLowerBound_t(element_type), FEAngio::NaturalCoordinatesUpperBound_t(element_type));

		return vec3d(r_dist(random_engine), s_dist(random_engine), t_dist(random_engine));
	}
}

int FragmentSeeder::IncrementCellCounter() {
	return initial_cell_id_counter++;
}
BEGIN_FECORE_CLASS(FragmentSeeder, FEMaterial)
	ADD_PARAMETER(number_fragments, "number_fragments");
	ADD_PARAMETER(proto_mat_cross, "proto_mat_cross");
	ADD_PARAMETER(cell_radius, "cell_radius");
END_FECORE_CLASS();

ByElementFragmentSeeder::ByElementFragmentSeeder(FEModel * model) : FragmentSeeder(model)
{
	
}

bool ByElementFragmentSeeder::SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index)
{
	std::uniform_int_distribution<int> edist(0, int(angio_elements.size()-1));

	if(angio_elements.size() == 0)
	{
		return false;
	}

	for (int i = 0; i < number_fragments; ++i)
	{
		Tip* r0 = new Tip();
		r0->TipCell = new FECell();
		r0->initial_fragment_id = initial_fragment_id_counter++;
		r0->TipCell->initial_cell_id = initial_cell_id_counter++;
		int elem_index = edist(angio_mat->m_pangio->rengine);
		r0->angio_element = angio_elements[elem_index];
		vec3d local_pos = GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(mesh, r0->angio_element, r0->angio_element->_rengine);
		r0->SetLocalPosition(local_pos, mesh);
		r0->TipCell->SetLocalPosition(local_pos, mesh);
		r0->time = -1;
		r0->TipCell->time = -1;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection(r0->angio_element->_rengine);
		r0->face = r0->angio_element;
		r0->TipCell->cell_radius = cell_radius;
		r0->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		// Finally add this to the AngioElement.
		r0->SetProtoGrowthLength(initial_segment_length);
		r0->angio_element->next_tips.at(r0->angio_element).push_back(r0);
	}
	return true;
}

ByElementFragmentSeederBiDirectional::ByElementFragmentSeederBiDirectional(FEModel * model) : FragmentSeeder(model)
{

}

bool ByElementFragmentSeederBiDirectional::SeedFragments(std::vector<AngioElement *> &angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index)
{
	std::uniform_int_distribution<int> edist(0, int (angio_elements.size()-1));

	if (angio_elements.size() == 0)
	{
		return false;
	}

	for (int i = 0; i < number_fragments; ++i)
	{
		Tip* r0 = new Tip();
		r0->TipCell = new FECell();
		r0->initial_fragment_id = initial_fragment_id_counter++;
		r0->TipCell->initial_cell_id = initial_cell_id_counter++;
		int elem_index = edist(angio_mat->m_pangio->rengine);
		r0->angio_element = angio_elements[elem_index];
		r0->TipCell->angio_element = r0->angio_element;
		r0->TipCell->ParentTip = r0;
		vec3d local_pos = GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(mesh, r0->angio_element, r0->angio_element->_rengine);
		r0->SetLocalPosition(local_pos, mesh);
		r0->TipCell->SetLocalPosition(local_pos, mesh);
		r0->time = -1;
		r0->TipCell->time = -1;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection(r0->angio_element->_rengine);
		r0->face = r0->angio_element;
		r0->TipCell->cell_radius = cell_radius;
		r0->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r0->SetProtoGrowthLength(initial_segment_length);
		FEModel* fem = GetFEModel();
		r0->TipCell->InitSpecies(mesh);
		// Finally add this to the AngioElement.
		r0->angio_element->next_tips.at(r0->angio_element).push_back(r0);
		// Now add an oppositely directed tip.
		Tip * r1 = new Tip();
		r1->TipCell = new FECell();
		r1->angio_element = r0->angio_element;
		r1->TipCell->angio_element = r1->angio_element;
		r1->TipCell->ParentTip = r1;
		r1->face = r0->face;
		r1->time = r0->time;
		r1->TipCell->time = r1->time;
		r1->growth_velocity = r0->growth_velocity;
		r1->SetLocalPosition(r0->GetLocalPosition(), mesh);
		r1->TipCell->SetLocalPosition(r1->GetLocalPosition(), mesh);
		r1->SetProtoGrowthLength(r0);
		r1->initial_fragment_id = initial_fragment_id_counter++;
		r1->TipCell->initial_cell_id = initial_cell_id_counter++;
		r1->use_direction = true;
		r1->direction = -r0->direction; r1->direction.unit();
		r1->TipCell->cell_radius = cell_radius;
		r1->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r1->TipCell->InitSpecies(mesh);
		r1->angio_element->next_tips.at(r1->angio_element).push_back(r1);
		Segment* seg = new Segment();
		seg->front = r0;
		seg->back = r1;
		r0->parent = seg;
		r1->parent = seg;
		seg->parent = seg;
	}
	return true;
}

ByVolumeFragmentSeeder::ByVolumeFragmentSeeder(FEModel * model) : FragmentSeeder(model)
{
	
}

bool ByVolumeFragmentSeeder::SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index)
{
	if (angio_elements.size() == 0)
	{
		return false;
	}

	// Random engine used to sample from distributions.
	angiofe_random_engine & random_engine = angio_mat->m_pangio->rengine;

	// A uniform distribution must be used, so take a sample from (0, total_volume) and find which element that the sample corresponds to by ordering 
	// the element's "starting" and "ending" volumes in two separate parallel arrays.
	double* element_beginning_volume = new double[angio_elements.size()];
	double* element_ending_volume    = new double[angio_elements.size()];

	// In order to have the sample be from a uniform distribution, accumulate the total volume of all the elements and take a sample from (0, total_volume).
	double total_volume = 0;

	for (size_t i = 0; i < angio_elements.size(); i++)
	{
		element_beginning_volume[i] = total_volume;
		total_volume += mesh->ElementVolume(*angio_elements.at(i)->_elem);
		element_ending_volume[i] = total_volume;
	}

	// Form the distribution.
	std::uniform_real_distribution<double> vol_dist(0.0, total_volume);

	// Based on the sample taken from the distribution, find the element that corresponds to the sample from the parallel arrays.
	for (int i = 0; i < number_fragments; i++)
	{
		double volume_sample = vol_dist(random_engine);

		// Perform a binary search of the pre-sorted parallel arrays (only requires one of them) to find the index of the element that will get the tip.
		size_t element_index = findElement(volume_sample, 0, int (angio_elements.size()), element_beginning_volume, element_ending_volume);

		// Build the fragment using a tip and build it in the proper element.
		Tip* r0 = new Tip();
		r0->TipCell = new FECell();
		r0->initial_fragment_id = initial_fragment_id_counter++;
		r0->TipCell->initial_cell_id = initial_cell_id_counter++;
		r0->angio_element = angio_elements[element_index];
		r0->TipCell->angio_element = r0->angio_element;
		vec3d local_pos = GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(mesh, r0->angio_element, r0->angio_element->_rengine);
		r0->SetLocalPosition(local_pos, mesh);
		r0->TipCell->SetLocalPosition(local_pos, mesh);
		r0->SetProtoGrowthLength(initial_segment_length);
		r0->time = -1;
		r0->TipCell->time = -1;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection(r0->angio_element->_rengine);
		r0->face = r0->angio_element;
		r0->TipCell->cell_radius = cell_radius;
		r0->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;

		// Finally add this to the AngioElement.
		r0->angio_element->next_tips.at(r0->angio_element).push_back(r0);
	}

	delete[] element_beginning_volume;
	delete[] element_ending_volume;

	return true;
}

ByVolumeFragmentSeederBiDirectional::ByVolumeFragmentSeederBiDirectional(FEModel * model) : FragmentSeeder(model)
{

}

bool ByVolumeFragmentSeederBiDirectional::SeedFragments(std::vector<AngioElement *>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index)
{
	if (angio_elements.size() == 0)
	{
		return false;
	}

	// Random engine used to sample from distributions.
	angiofe_random_engine & random_engine = angio_mat->m_pangio->rengine;

	// A uniform distribution must be used, so take a sample from (0, total_volume) and find which element that the sample corresponds to by ordering 
	// the element's "starting" and "ending" volumes in two separate parallel arrays.
	double* element_beginning_volume = new double[angio_elements.size()];
	double* element_ending_volume = new double[angio_elements.size()];

	// In order to have the sample be from a uniform distribution, accumulate the total volume of all the elements and take a sample from (0, total_volume).
	double total_volume = 0;

	for (size_t i = 0; i < angio_elements.size(); i++)
	{
		element_beginning_volume[i] = total_volume;
		total_volume += mesh->ElementVolume(*angio_elements.at(i)->_elem);
		element_ending_volume[i] = total_volume;
	}

	// Form the distribution.
	std::uniform_real_distribution<double> vol_dist(0.0, total_volume);

	// Based on the sample taken from the distribution, find the element that corresponds to the sample from the parallel arrays.
	for (int i = 0; i < number_fragments; i++)
	{
		double volume_sample = vol_dist(random_engine);

		// Perform a binary search of the pre-sorted parallel arrays (only requires one of them) to find the index of the element that will get the tip.
		size_t element_index = findElement(volume_sample, 0, int(angio_elements.size()-1), element_beginning_volume, element_ending_volume);

		// Build the fragment using a tip and build it in the element.
		Tip* r0 = new Tip();
		r0->TipCell = new FECell();
		r0->initial_fragment_id = initial_fragment_id_counter++;
		r0->TipCell->initial_cell_id = initial_cell_id_counter++;
		r0->angio_element = angio_elements[element_index];
		r0->TipCell->angio_element = r0->angio_element;
		vec3d local_pos = GetRandomVectorPositionWithinNaturalCoordinateBoundsByElementType(mesh, r0->angio_element, angio_mat->m_pangio->rengine);
		r0->SetLocalPosition(local_pos, mesh);
		r0->TipCell->SetLocalPosition(local_pos, mesh);
		r0->time = -1;
		r0->TipCell->time = -1;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection(angio_mat->m_pangio->rengine);
		r0->face = r0->angio_element;
		r0->TipCell->cell_radius = cell_radius;
		r0->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r0->SetProtoGrowthLength(initial_segment_length);

		// Finally add this to the AngioElement.
		r0->angio_element->next_tips.at(r0->angio_element).push_back(r0);

		// Now add an oppositely directed tip.
		Tip * r1 = new Tip(r0, mesh);
		r1->TipCell = new FECell();
		r1->initial_fragment_id = initial_fragment_id_counter++;
		r1->TipCell->initial_cell_id = initial_cell_id_counter++;
		r1->TipCell->angio_element = r1->angio_element;
		r1->use_direction = true;
		r1->direction = -r0->direction;
		r1->TipCell->cell_radius = cell_radius;
		r1->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r1->angio_element->next_tips.at(r1->angio_element).push_back(r1);
		r1->SetProtoGrowthLength(r0);
		Segment* seg = new Segment();
		seg->front = r0;
		seg->back = r1;
		r0->parent = seg;
		r1->parent = seg;
		seg->parent = seg;
	}

	delete[] element_beginning_volume;
	delete[] element_ending_volume;

	return true;
}

SingleCellSeeder::SingleCellSeeder(FEModel* model) : FragmentSeeder(model)
{

}

bool SingleCellSeeder::SeedFragments(std::vector<AngioElement*>& angio_elements, FEMesh* mesh, FEAngioMaterial* angio_mat, int buffer_index)
{
	std::uniform_int_distribution<int> edist(0, int(angio_elements.size() - 1));

	if (angio_elements.size() == 0)
	{
		return false;
	}

	for (int i = 0; i < number_fragments; ++i)
	{
		Tip* r0 = new Tip();
		r0->TipCell = new FECell();
		r0->initial_fragment_id = initial_fragment_id_counter++;
		r0->TipCell->initial_cell_id = initial_cell_id_counter++;
		int elem_index = edist(angio_mat->m_pangio->rengine);
		r0->angio_element = angio_elements[elem_index];
		r0->TipCell->angio_element = r0->angio_element;
		//r0->TipCell->Init();
		r0->TipCell->ParentTip = r0;
		vec3d local_pos = initial_position;
		r0->SetLocalPosition(local_pos, mesh);
		r0->TipCell->SetLocalPosition(local_pos, mesh);
		r0->time = -1;
		r0->TipCell->time = -1;
		r0->use_direction = true;
		r0->direction = angio_mat->m_pangio->uniformRandomDirection(r0->angio_element->_rengine);
		r0->face = r0->angio_element;
		r0->TipCell->cell_radius = cell_radius;
		r0->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r0->SetProtoGrowthLength(initial_segment_length);
		FEModel* fem = GetFEModel();
		r0->TipCell->InitSpecies(mesh);
		// Finally add this to the AngioElement.
		r0->angio_element->next_tips.at(r0->angio_element).push_back(r0);
		// Now add an oppositely directed tip.
		Tip* r1 = new Tip();
		r1->TipCell = new FECell();
		r1->angio_element = r0->angio_element;
		r1->TipCell->angio_element = r1->angio_element;
		r1->TipCell->ParentTip = r1;
		r1->face = r0->face;
		r1->time = r0->time;
		r1->TipCell->time = r1->time;
		r1->growth_velocity = r0->growth_velocity;
		r1->SetLocalPosition(r0->GetLocalPosition(), mesh);
		r1->TipCell->SetLocalPosition(r1->GetLocalPosition(), mesh);
		r1->SetProtoGrowthLength(r0);
		r1->initial_fragment_id = initial_fragment_id_counter++;
		r1->TipCell->initial_cell_id = initial_cell_id_counter++;
		r1->use_direction = true;
		r1->direction = -r0->direction; r1->direction.unit();
		r1->TipCell->cell_radius = cell_radius;
		r1->TipCell->cell_volume = pow(cell_radius, 3) * 4 / 3 * PI;
		r1->TipCell->InitSpecies(mesh);
		r1->angio_element->next_tips.at(r1->angio_element).push_back(r1);
	}
	return true;
}
BEGIN_FECORE_CLASS(SingleCellSeeder,FragmentSeeder)
ADD_PARAMETER(initial_position, "initial_position");
END_FECORE_CLASS();