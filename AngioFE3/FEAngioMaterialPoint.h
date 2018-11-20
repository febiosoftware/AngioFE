#pragma once
#include "StdAfx.h"
#include <FECore/FEMaterialPoint.h>

//-----------------------------------------------------------------------------
// A new material point class is defined to store the elastic parameters for 
// each integration point.
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	explicit FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt);

	//! The init function is used to intialize data
	void Init() override;

	//! update the material point to a given time
	void Update(const FETimeInfo& timeInfo) override;

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

	//! pointer to material point of the vessel
	FEMaterialPoint* vessPt;
	//! pointer to matrix material point
	FEMaterialPoint* matPt;

	//! parameter list
	DECLARE_PARAMETER_LIST();

public:
	//! return the angio material point from a material point if it exists, return nullptr otherwise
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};