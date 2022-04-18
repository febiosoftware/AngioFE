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
	// the FA can be calculated as std/rms of the ODF. We can assume each direction is one sample and perform the calculation on the eigenvalues.
	double sum = d[0]+d[1]+d[2];
	double mean = sum / 3.0;
	double sq_sum = (d[0] - mean) * (d[0] - mean) + (d[1] - mean) * (d[1] - mean) + (d[2] - mean) * (d[2] - mean);
	double stdev = sqrt(sq_sum / 2.0);
	double rms = sqrt((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / 3.0);
	angioFA = stdev/rms;
}

void FEAngioMaterialPoint::UpdateSPD()
{
	// get the material point data and polar decomposition of F
	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	// get the net deformation gradient which includes the mapping of a sphere to the initial angio SPD and the applied deformation.
	mat3d Fnet = emp->m_F*initial_angioSPD; 
	// get the Left cauchy-green deformation tensor
	mat3d Bh = Fnet * Fnet.transpose();
	mat3ds B = mat3ds(Bh[0][0], Bh[1][1], Bh[2][2], Bh[0][1], Bh[1][2], Bh[0][2]);
	// get the ellipsoid associated with B
	double l2[3];
	vec3d v[3];
	B.eigen2(l2, v);
	// get V from B
	mat3ds V = dyad(v[0]) * sqrt(l2[0]) + dyad(v[1]) * sqrt(l2[1]) + dyad(v[2]) * sqrt(l2[2]);
	double b[3]; vec3d n[3];
	// get the ellipsoid associated with V
	V.eigen2(b, n);
	//// Construct the SPD from components of V's eigensolution.
	mat3dd d = mat3dd(b[0], b[1], b[2]);
	mat3d q; q.setCol(0, n[0]); q.setCol(1, n[1]); q.setCol(2, n[2]);
	mat3d p = q*d*q.transpose();
	// input to SPD
	angioSPD = mat3ds(p[0][0], p[1][1], p[2][2], p[0][1], p[1][2], p[0][2]);
	angioSPD = (3.0 / angioSPD.tr())*angioSPD;	
}