#include "Tip.h"
#include "Segment.h"
#include <iostream>

vec3d Tip::GetDirection(FEMesh* mesh) const
{
	if(is_branch || use_direction)
	{
		return direction;
	}
	else if(parent)
	{
		return parent->Direction(mesh);
	}
	else
	{
		assert(false);
		vec3d d;
		return d;
	}
}

vec3d Tip::GetPosition(FEMesh * mesh) const
{
	double arr[FEElement::MAX_NODES];
	assert(angio_element);
	assert(angio_element->_elem);
	angio_element->_elem->shape_fnc(arr, local_pos.x, local_pos.y, local_pos.z);
	vec3d rc(0, 0, 0);

	for (int j = 0; j < angio_element->_elem->Nodes(); j++)
	{
		rc += mesh->Node(angio_element->_elem->m_node[j]).m_rt* arr[j];
	}
	return rc;
}
void Tip::PrintTipInfo(FEMesh *mesh, std::string title) const
{
	std::cout << title << std::endl;
	PrintTipInfo(mesh);
}

void Tip::PrintTipInfo(FEMesh *mesh) const
{
	std::cout << "local position: " << local_pos.x << " , " << local_pos.y << " , " << local_pos.z << std::endl;
	vec3d pos = GetPosition(mesh);
	std::cout << "global position: " << pos.x << " , " << pos.y << " , " << pos.z << std::endl;
	std::cout << "element id: " << angio_element->_elem->GetID() << std::endl;
	if(face)
	{
		std::cout << "face id: " << face->_elem->GetID() << std::endl;
	}
	else
	{
		std::cout << "no next face" << std::endl;
	}

	vec3d dir = GetDirection(mesh);
	std::cout << "direction: " << dir.x << " , " << dir.y << " , " << dir.z << std::endl;
	
	std::cout << "is branch: " << is_branch << std::endl;
	std::cout << "use direction: " << use_direction << std::endl;
	std::cout << std::endl;
}

Tip::Tip(Tip * other, FEMesh * mesh)
{
	time = other->time;
	growth_velocity = other->growth_velocity;
	//deparent the new tip
	initial_fragment_id = other->initial_fragment_id;
	direction = other->GetDirection(mesh);
}