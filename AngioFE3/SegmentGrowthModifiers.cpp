#include "SegmentGrowthModifiers.h"
#include "Tip.h"
#include "angio3d.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FESolidDomain.h>

vec3d PreviousSegmentPSC::ApplyModifiers(vec3d prev, Tip* tip, FEMesh* mesh)
{
	return tip->GetDirection(mesh);
}

void FiberPDD::Update(FEMesh * mesh, FEAngio* angio)
{
	/*
	angio->
	// build the element data array
	for (int n = 0; n < 4; n++)
	{
		fiber_at_int_pts[n].clear();
		fiber_at_int_pts[n].resize(NE);

		for (int i = 0; i<NE; ++i)
		{
			FESolidElement& e = sd->Element(i);
			int nint = e.GaussPoints();
			fiber_at_int_pts[n][i].assign(nint, 0.0);
		}
	}


	// this array will store the results
	FESPRProjection map;

	// loop over stress components

	// fill the ED array
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = sd->Element(i);
		int nint = el.GaussPoints();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint * mp = el.GetMaterialPoint(j);
			FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
			vec3d fd = emp->m_F * emp->m_Q * vec3d(1, 0, 0);
			double lambda = fd.unit();

			//store the data
			fiber_at_int_pts[0][i][j] = fd.x;
			fiber_at_int_pts[1][i][j] = fd.y;
			fiber_at_int_pts[2][i][j] = fd.z;
			fiber_at_int_pts[3][i][j] = lambda;
		}
	}

	for (int n = 0; n<4; ++n)
	{
		// project to nodes
		map.Project(*sd, fiber_at_int_pts[n], fibers_at_nodes[n]);
	}
	*/
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
		axis = emp->m_F * emp->m_Q * axis;
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
