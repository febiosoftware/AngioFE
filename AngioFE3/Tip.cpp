#include "Tip.h"
#include "Segment.h"

vec3d Tip::GetDirection(FEMesh* mesh) const
{
	if(is_branch)
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