#pragma once
#include <FECore/FEMaterial.h>

class FESolidElement;


//this may need to be updated to do more quadrature related interpolation of values
//this might do the calculations more than is needed if multiple tips are in the same element at the same time
//consider rewritting this once data has been collected on actual simulations
class FEVariableInterpolation : public FEMaterial
{
public:
	explicit FEVariableInterpolation(FEModel * pfem) : FEMaterial(pfem){}
	virtual double Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) = 0;
};

//discontinuous interp[olation of values at the gauss points
class PerElementVI : public FEVariableInterpolation
{
public:
	explicit PerElementVI(FEModel * pfem) : FEVariableInterpolation(pfem) {}
	double Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) override;
};