#include "StdAfx.h"
#include <algorithm>
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMixture.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include "angio3d.h"

//-----------------------------------------------------------------------------
FEAngioMaterialPoint::FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt) : FEMaterialPoint(pt)
{
	vessPt = vesselPt;
	matPt = matrixPt;
	matrix_weight = 1.0;
	vessel_weight = 0.0;
	ref_ecm_density = 3.0;
	vessPt->SetPrev(this);
	matPt->SetPrev(this);
	m_as.zero();
}

//-----------------------------------------------------------------------------
//! The init function is used to intialize data
void FEAngioMaterialPoint::Init()
{
	FEMaterialPoint::Init();
	vessPt->Init();
	matPt->Init();
}

void FEAngioMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPoint::Update(timeInfo);
	matPt->Update(timeInfo);
	vessPt->Update(timeInfo);
}

//-----------------------------------------------------------------------------
// define the material parameters
//BEGIN_FECORE_CLASS(FEAngioMaterialPoint, FEMaterialPoint)
//ADD_PARAMETER(ref_ecm_density, "ref_ecm_density");
//END_FECORE_CLASS();






//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
FEMaterialPoint* FEAngioMaterialPoint::Copy()
{
	FEAngioMaterialPoint* pt = new FEAngioMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
void FEAngioMaterialPoint::Serialize(DumpStream& dmp)
{

	FEMaterialPoint::Serialize(dmp);
}

FEAngioMaterialPoint* FEAngioMaterialPoint::FindAngioMaterialPoint(FEMaterialPoint* mp)
{
	FEAngioMaterialPoint* angioPt = dynamic_cast<FEAngioMaterialPoint*>(mp);
	if (angioPt)
		return angioPt;

	FEMaterialPoint* pt = mp;
	while (pt)
	{
		angioPt = dynamic_cast<FEAngioMaterialPoint*>(pt);
		if (angioPt)
			return angioPt;

		FEElasticMixtureMaterialPoint* mixtureP = dynamic_cast<FEElasticMixtureMaterialPoint*>(pt);
		if (mixtureP)
		{
			for (int i = 0; i<mixtureP->Components(); i++)
			{
				//TODO: is the recursion needed or not(is this search too deep?)
				angioPt = FindAngioMaterialPoint(mixtureP->GetPointData(i));
				if (angioPt)
				{
					return angioPt;
				}
			}
		}

		pt = pt->Next();
	}

	return nullptr;
}

void FEAngioMaterialPoint::UpdateAngioFractionalAnisotropy()
{
	//sort the spd
	std::vector<pair<double, int>> v;
	double d[3]; vec3d r[3];
	angioSPD.eigen2(d, r);
	v.push_back(pair<double, int>(d[0], 0));
	v.push_back(pair<double, int>(d[1], 1));
	v.push_back(pair<double, int>(d[2], 2));
	sort(v.begin(), v.end(), sortinrev);
	// determine the "eigenvalues"
	double a = v[0].first;
	double b = v[1].first;
	double c = v[2].first;
	// calculate the fractional anisotropy
	angioFA = sqrt(0.5)*(sqrt(pow(a - b, 2) + pow(b - c, 2) + pow(c - a, 2)) / (sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))));
}

void FEAngioMaterialPoint::UpdateSPD()
{
	// get the material point data and polar decomposition of F
	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	mat3d F = emp->m_F;
	mat3d Uinv = emp->RightStretchInverse();
	mat3d R = F*Uinv;
	// get the initial semiprincipal axes magnitudes and directions
	double B[3]; vec3d N[3];
	initial_angioSPD = (3.0 / initial_angioSPD.tr())*initial_angioSPD;
	initial_angioSPD.eigen2(B, N);
	// get the stretch ratios
	double b[3]; vec3d n[3]; double lam[3];
	lam[0] = (F*N[0]).norm(); b[0] = B[0] * lam[0];
	lam[1] = (F*N[1]).norm(); b[1] = B[1] * lam[1];
	lam[2] = (F*N[2]).norm(); b[2] = B[2] * lam[2];
	// Assemble the diagonal matrix
	mat3dd d = mat3dd(b[0], b[1], b[2]);
	// Assemble the matrix containing the rotated directions
	mat3d q; q.setCol(0, R*N[0]); q.setCol(1, R*N[1]); q.setCol(2, R*N[2]);
	// Construct the SPD from the eigen components
	mat3d p = q*d*q.transpose();
	// input to SPD
	angioSPD = mat3ds(p[0][0], p[1][1], p[2][2], p[0][1], p[1][2], p[0][2]);
	angioSPD = (3.0 / angioSPD.tr())*angioSPD;	
}