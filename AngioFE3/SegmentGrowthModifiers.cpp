#include "SegmentGrowthModifiers.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>

vec3d PreviousSegmentPSC::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	return tip->GetDirection(mesh);
}

void FiberPDD::Update(FEMesh * mesh)
{
	
}

/*
vec3d FiberPDD::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	std::vector<double> x_data, y_data, z_data;
	std::vector<double> x_data_nodal, y_data_nodal, z_data_nodal;
	for(int i=0; i< tip->angio_element->_elem->GaussPoints();i++)
	{
		FEMaterialPoint * mp = tip->angio_element->_elem->GetMaterialPoint(i);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		vec3d axis(1, 0, 0);
		axis = emp->m_Q * axis;
		x_data.push_back(axis.x);
		y_data.push_back(axis.y);
		z_data.push_back(axis.z);
	}
	FESolidDomain * se = dynamic_cast<FESolidDomain*>(tip->angio_element->_elem->GetDomain());
	//projection.Project(*se, x_data, x_data_nodal);
	//vec3d fiber_direction =
	return mix(prev, fiber_direction, contribution);
}
*/