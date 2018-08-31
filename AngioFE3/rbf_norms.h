#pragma once
#include <FECore/FEMaterial.h>

class rbf_norm : public FEMaterial
{
public:
	explicit rbf_norm(FEModel* pfem) : FEMaterial(pfem) {}
	virtual double norm(vec3d pos, vec3d center) = 0;
protected:
	double max_intpt_contribution = 100000;
	double cut_off_point = 1e-6;
};

class gaussian_rbf_norm : public rbf_norm
{
public:
	explicit gaussian_rbf_norm(FEModel* pfem) : rbf_norm(pfem) {}
	double norm(vec3d pos, vec3d center) override;
};
