#pragma once
#include <FECore/FEMaterial.h>

class FESolidElement;


//! this may need to be updated to do more quadrature related interpolation of values
//! this might do the calculations more than is needed if multiple tips are in the same element at the same time
//! consider rewritting this once data has been collected on actual simulations

//! base class for doing interpolation to all locations within an element
class FEVariableInterpolation : public FEMaterial {
public:
	//! constructor
	explicit FEVariableInterpolation(FEModel * pfem) : FEMaterial(pfem){}
	//! interpolate doubles
	virtual double Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) = 0;
	//! interpolate directions
	virtual quatd Interpolate(FESolidElement *se, std::vector<quatd> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) = 0;
};

//! discontinuous interpolation of values at the gauss points
// issue
class PerElementVI : public FEVariableInterpolation {
public:
	//! constructor
	explicit PerElementVI(FEModel * pfem) : FEVariableInterpolation(pfem){}
	//! interpolate values on a per element basis
	double Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) override;
	//! interpolate values on a per element basis
	quatd Interpolate(FESolidElement *se, std::vector<quatd> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh) override;
};

// FEMixMethod, class for interpolating the contribution of two vectors.
// issue
class FEMixMethod : public FEMaterial {
public:
	// constructor
	explicit FEMixMethod(FEModel * pfem) : FEMaterial(pfem){}
	// interpolate doubles
	virtual vec3d ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) = 0;
	// interpolate along an axis (selects direction of growth)
	virtual vec3d ApplyMixAxis(vec3d psc_dir, vec3d pdd_dir, double contribution) = 0;
};

// Legacy method. This mixes the two vectors using a linear interpolation of their components
class LinInterp : public FEMixMethod {
public:
	explicit LinInterp(FEModel * pfem) : FEMixMethod(pfem) {}
	vec3d ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) override;
	vec3d ApplyMixAxis(vec3d psc_dir, vec3d pdd_dir, double contribution) override;
};

// New method. This determines a rotation matrix for the rotation in the plane spdnned by the two vectors then the first vector is rotated towards the second by a scale between 0-1 where 0 returns the original vector and 1 returns the second.
class LinRot : public FEMixMethod {
public:
	explicit LinRot(FEModel * pfem) : FEMixMethod(pfem) {}
	vec3d ApplyMix(vec3d psc_dir, vec3d pdd_dir, double contribution) override;
	vec3d ApplyMixAxis(vec3d psc_dir, vec3d pdd_dir, double contribution) override;
};