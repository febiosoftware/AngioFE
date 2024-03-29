#include "Tip.h"
#include "Segment.h"
#include <iostream>
#include <algorithm>
#include "FEAngio.h"
#include "FEProbabilityDistribution.h"
#include <unordered_map>

vec3d Tip::GetDirection(FEMesh* mesh) const
{
	if (is_branch || use_direction)
		return direction;
	else if (parent)
		return parent->Direction(mesh);
	else
	{
		assert(false);
		vec3d d;
		return d;
	}
}

vec3d Tip::GetPosition(FEMesh* mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	auto se = angio_element->_elem;
	se->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < se->Nodes(); j++)
		rc += mesh->Node(se->m_node[j]).m_rt * arr[j];

	return rc;
}

vec3d Tip::GetRefPosition(FEMesh* mesh) const
{
	double arr[FESolidElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	auto se = angio_element->_elem;
	se->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < se->Nodes(); j++)
		rc += mesh->Node(se->m_node[j]).m_r0 * arr[j];

	return rc;
}

void Tip::PrintTipInfo(FEMesh* mesh, std::string title) const
{
#ifndef NDEBUG
	std::cout << title << std::endl;
	PrintTipInfo(mesh);
#endif
}

void Tip::PrintTipInfo(FEMesh* mesh) const
{
#ifndef NDEBUG
	std::cout << "local position: " << local_pos.x << " , " << local_pos.y << " , " << local_pos.z << std::endl;
	vec3d pos = GetPosition(mesh);
	std::cout << "global position: " << pos.x << " , " << pos.y << " , " << pos.z << std::endl;
	std::cout << "element id: " << angio_element->_elem->GetID() << std::endl;
	if (face)
		std::cout << "face id: " << face->_elem->GetID() << std::endl;
	else
		std::cout << "no next face" << std::endl;

	vec3d dir = GetDirection(mesh);
	std::cout << "direction: " << dir.x << " , " << dir.y << " , " << dir.z << std::endl;

	std::cout << "is branch: " << is_branch << std::endl;
	std::cout << "use direction: " << use_direction << std::endl;
	std::cout << std::endl;
#endif
}

// create a new tip based on the pre-existing tip.
Tip::Tip(Tip* other, FEMesh* mesh)
{
	// get the other tip's parameters.
	angio_element = other->angio_element;
	face = other->face;
	time = other->time;
	seg_time = other->seg_time;
	growth_velocity = other->growth_velocity;
	SetProtoGrowthLength(other);
	local_pos = other->local_pos;
	// deparent the new tip
	initial_fragment_id = other->initial_fragment_id;
	direction = other->GetDirection(mesh);
	direction.unit();
	FEModel* fem = angio_element->_mat->GetFEModel();
}

void Tip::SetLocalPosition(vec3d pos, FEMesh* mesh)
{
	assert(angio_element);
	auto se = angio_element->_elem;
	assert(se->Type() == FE_Element_Type::FE_TET4G4 
						? (pos.x + pos.y + pos.z) < 1.01 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G4 
						? (pos.x < 1.01) && (pos.x >= -0.01) 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G4 
						? (pos.y < 1.01) && (pos.y >= -0.01) 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G4 
						? (pos.z < 1.01) && (pos.z >= -0.01) 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G1 
						? (pos.x + pos.y + pos.z) < 1.01 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G1 
						? (pos.x < 1.01) && (pos.x >= -0.01) 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G1 
						? (pos.y < 1.01) && (pos.y >= -0.01) 
						: true);
	assert(se->Type() == FE_Element_Type::FE_TET4G1 
						? (pos.z < 1.01) && (pos.z >= -0.01) 
						: true);

	local_pos = pos;
	local_pos.x 
		= std::max(
			std::min(	FEAngio::NaturalCoordinatesUpperBound_r(se->Type()), 
						local_pos.x),
			FEAngio::NaturalCoordinatesLowerBound_r(se->Type()));
	local_pos.y 
		= std::max(
			std::min(	FEAngio::NaturalCoordinatesUpperBound_s(se->Type()), 
						local_pos.y), 
			FEAngio::NaturalCoordinatesLowerBound_s(se->Type()));
	local_pos.z 
		= std::max(
			std::min(	FEAngio::NaturalCoordinatesUpperBound_t(se->Type()), 
						local_pos.z), 
			FEAngio::NaturalCoordinatesLowerBound_t(se->Type()));
}

vec3d Tip::GetLocalPosition() const
{
	return local_pos;
}

void Tip::SetProtoGrowthLength(FEProbabilityDistribution* dist)
{
	proto_growth_length = 0.5 * dist->NextValue(angio_element->_rengine);
}

void Tip::SetProtoGrowthLength(Tip* tip)
{
	proto_growth_length = tip->GetProtoGrowthLength();
}

void Tip::SetProtoGrowthLength(double initial_length)
{
	proto_growth_length = initial_length;
}

double Tip::GetProtoGrowthLength() 
{
	return proto_growth_length;
}