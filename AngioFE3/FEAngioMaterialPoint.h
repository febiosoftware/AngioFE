#pragma once
#include "StdAfx.h"
#include "AngioElement.h"
#include <FECore/FEMaterialPoint.h>

//! per integration point values for angiofe
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	explicit FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt);

	//! The init function is used to intialize data
	void Init() override;

	//! update the material point to a given time
	void Update(const FETimeInfo& timeInfo) override;

	void UpdateAngioFractionalAnisotropy();
	//! Update the angioSPD based on deformation/rotation
	void UpdateSPD();
	//! Construct the elliptical distribution between 2 orthogonal SPD and sample it
	double GetEllipseAngle(const double a, const double b, const double dist_min, const double dist_max, const int n, AngioElement* angio_elem);
	//double GetEllipseAngle2(const double a, const double b);
	//! Use Gerard's suggested approach for getting an ellipse angle
	double GetEllipseAngleAteshian(const double a, const double b, const double dist_min, const double dist_max, const int n, AngioElement* angio_elem);

	//! copy material point data (for running restarts) todo Is this still used?
	FEMaterialPoint* Copy() override;

	//! copy material point data (for running restarts) todo Is this still used?
	void Serialize(DumpStream& dmp) override;

	//! ecm density in the reference configuration
	double ref_ecm_density;

	//! the value which may be used to repulse vessels from material points
	double repulse_value = 0.0;

	//! the temporary to store angio stress
	mat3ds m_as;

	//! the weight of the vessel material
	double vessel_weight;

	//! the weight of the matrix material
	double matrix_weight;

	//! previous deformation gradient
	mat3d Fp = mat3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
	mat3d Fi = mat3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
	bool nhit = 0;

	//! pointer to material point of the vessel
	FEMaterialPoint* vessPt = nullptr;
	//! pointer to matrix material point
	FEMaterialPoint* matPt = nullptr;

	//! initial orientation of spd
	mat3ds initial_angioSPD;
	//! updated angioSPD
	mat3ds angioSPD;
	//! Angio fractional anisotropy
	double angioFA;
	vec3d angio_fiber_dir = vec3d(1, 0, 0);

	//! parameter list
	//DECLARE_FECORE_CLASS();
	// only works with classes that derive from FECoreBase

public:
	//! return the angio material point from a material point if it exists, return nullptr otherwise
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};