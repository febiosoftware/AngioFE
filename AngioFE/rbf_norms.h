#pragma once
#include <FECore/FEMaterial.h>
/*

class rbf_norm : public FEMaterialProperty
{
public:
	explicit rbf_norm(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual double rbfnorm(vec3d pos, vec3d center) = 0;
protected:
	double max_intpt_contribution = 100000;
	double cut_off_point = 1e-6;
};

class gaussian_rbf_norm : public rbf_norm
{
public:
	explicit gaussian_rbf_norm(FEModel* pfem) : rbf_norm(pfem) {}
	double rbfnorm(vec3d pos, vec3d center) override;
};

*/