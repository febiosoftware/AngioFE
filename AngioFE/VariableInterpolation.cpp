#include "VariableInterpolation.h"
#include "FECore/FESolidElement.h"
#include "FECore/FEMesh.h"
#include "angio3d.h"
#include <iostream>

double PerElementVI::Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh)
{
	double ao[FESolidElement::MAX_NODES];
	double H[FESolidElement::MAX_NODES];
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
	double nw[FESolidElement::MAX_NODES], nx[FESolidElement::MAX_NODES], ny[FESolidElement::MAX_NODES], nz[FESolidElement::MAX_NODES];
	double gw[FESolidElement::MAX_NODES], gx[FESolidElement::MAX_NODES], gy[FESolidElement::MAX_NODES], gz[FESolidElement::MAX_NODES];
	double H[FESolidElement::MAX_NODES];
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

vec3d LinInterp::ApplyMixAxis(vec3d psc_dir, vec3d pdd_dir, double contribution) {
	if (psc_dir * pdd_dir < 0) {pdd_dir = -pdd_dir; }
	return mix(psc_dir, pdd_dir, contribution);
}

vec3d LinRot::ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) {
	return mix3d(psc_dir, pdd_dir, contribution);
}

vec3d LinRot::ApplyMixAxis(vec3d psc_dir, vec3d pdd_dir, double contribution) {
	if (psc_dir * pdd_dir < 0) 	{pdd_dir = -pdd_dir;}
	return mix3d(psc_dir, pdd_dir, contribution);
}