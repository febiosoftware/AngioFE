#include "VariableInterpolation.h"
#include "FECore/FEElement.h"
#include "FECore/FEMesh.h"
#include "angio3d.h"

double PerElementVI::Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh)
{
	double ao[FEElement::MAX_NODES];
	double H[FEElement::MAX_NODES];
	assert(values_at_gauss_points.size() == se->GaussPoints());
	//this hack might will only work with vectors
	se->project_to_nodes(&values_at_gauss_points[0], ao);

	
	se->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	double val = 0.0;
	for (int j = 0; j < se->Nodes(); j++)
	{
		val += ao[j]* H[j];
	}
	return val;
}

quatd PerElementVI::Interpolate(FESolidElement *se, std::vector<quatd> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh)
{
	double nw[FEElement::MAX_NODES], nx[FEElement::MAX_NODES], ny[FEElement::MAX_NODES], nz[FEElement::MAX_NODES];
	double gw[FEElement::MAX_NODES], gx[FEElement::MAX_NODES], gy[FEElement::MAX_NODES], gz[FEElement::MAX_NODES];
	double H[FEElement::MAX_NODES];
	assert(values_at_gauss_points.size() == se->GaussPoints());
	//this hack might will only work with vectors
	//se->project_to_nodes(&values_at_gauss_points[0], ao);
	for(int i=0; i < se->GaussPoints();i++)
	{
		gw[i] = values_at_gauss_points[i].w;
		gx[i] = values_at_gauss_points[i].x;
		gy[i] = values_at_gauss_points[i].y;
		gz[i] = values_at_gauss_points[i].z;
	}
	se->project_to_nodes(gw, nw);
	se->project_to_nodes(gx, nx);
	se->project_to_nodes(gy, ny);
	se->project_to_nodes(gz, nz);

	se->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	quatd val;
	for (int j = 0; j < se->Nodes(); j++)
	{
		val.w += nw[j] * H[j];
		val.x += nx[j] * H[j];
		val.y += ny[j] * H[j];
		val.z += nz[j] * H[j];
	}
	return val;
}

vec3d LinInterp::ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) {
	return mix(psc_dir, pdd_dir, contribution);
}

vec3d LinRot::ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) {
	return mix3d(psc_dir, pdd_dir, contribution);
}