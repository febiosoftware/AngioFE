#include "Tip.h"
#include "Segment.h"

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

Tip::Tip(Tip * other, FEMesh * mesh)
{
	local_pos = other->local_pos;
	angio_element = other->angio_element;
	time = other->time;
	growth_velocity = other->growth_velocity;
	face = other->face;
	//deparent the new tip
	initial_fragment_id = other->initial_fragment_id;
	direction = other->GetDirection(mesh);
}