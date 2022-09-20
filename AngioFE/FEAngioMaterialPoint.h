#pragma once
#include "StdAfx.h"
#include "AngioElement.h"
#include <FECore/FEMaterialPoint.h>

//! per integration point values for angiofe
class FEAngioMaterialPoint : public FEMaterialPointData
{
public:
	//! constructor
	explicit FEAngioMaterialPoint(FEMaterialPointData* pt, FEMaterialPointData* vesselPt, FEMaterialPointData* matrixPt);

	//! The init function is used to intialize data
	void Init() override;

	//! update the material point to a given time
	void Update(const FETimeInfo& timeInfo) override;

	void UpdateAngioFractionalAnisotropy();
	//! Update the angioSPD based on deformation/rotation
	void UpdateSPD();
	
	//! copy material point data (for running restarts) todo Is this still used?
	FEMaterialPointData* Copy() override;

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

	double vascular_density;

	//! pointer to material point of the vessel
	FEMaterialPoint* vessPt = nullptr;

	//! pointer to matrix material point
	FEMaterialPoint* matPt = nullptr;

	//! initial orientation of efd/spd
	mat3ds initial_angioSPD;
	//! updated angioSPD
	mat3ds angioSPD;
	//! Angio fractional anisotropy
	double angioFA;
	//! Angio fiber direction
	vec3d angio_fiber_dir = vec3d(1, 0, 0);

public:
	//! return the angio material point from a material point if it exists, return nullptr otherwise
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};