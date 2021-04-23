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

double FEAngioMaterialPoint::GetEllipseAngle(const double a, const double b, const double dist_min, const double dist_max, const int n, AngioElement* angio_elem)
{
	// initialize vals
	double p_i = dist_min;
	double p_f = dist_max;
	double div = (PI / 180)*((n + 1) / n);

	// theta
	std::vector<double> t(n);
	// radius
	std::vector<double> r_t(n);
	// radius cumulative sum
	std::vector<double> rc_t(n);
	// area vector
	std::vector<double> a_t(n);
	// cumulative area vector
	std::vector<double> ac_t(n);

	t.at(0) = dist_min;
	a_t.at(0) = 0;
	r_t.at(0) = (a*b) / sqrt(pow(b*cos(t.at(0)), 2) + pow(a*sin(t.at(0)), 2));
	for (int i = 1; i < n; i++)
	{
		// assign the angle
		t.at(i) = (p_i + i*div);
		// get the radius
		r_t.at(i) = (a*b) / sqrt(pow(b*cos(t.at(i)), 2) + pow(a*sin(t.at(i)), 2));
		// get the area
		a_t.at(i) = (r_t.at(i)*r_t.at(i - 1))*sin(div);
	}
	// get the cumulative sum
	std::partial_sum(a_t.begin(), a_t.end(), ac_t.begin());
	// divide cumulative sum by sum
	std::transform(ac_t.begin(), ac_t.end(), ac_t.begin(),
		std::bind(std::divides<double>(), std::placeholders::_1, ac_t.at(n - 1)));
	// construct uniform distribution
	std::uniform_real_distribution<double> ud = std::uniform_real_distribution<double>(0.0, 1.0);
	// get a random number
	
	double rn = ud(angio_elem->_rengine);
	
	// find the rc_t value closest to the random number and get the position
	int fi = std::distance(ac_t.begin(), std::lower_bound(ac_t.begin(), ac_t.end(), rn));
	return t.at(fi);
}

double FEAngioMaterialPoint::GetEllipseAngleAteshian(const double a, const double b, const double dist_min, const double dist_max, const int n, AngioElement* angio_elem)
{
	// initialize vals
	double p_i = dist_min;
	double p_f = dist_max;
	double div = (PI / 180)*((n + 1) / n);

	// theta
	std::vector<double> t(n);
	// radius
	std::vector<double> r_t(n);
	// radius cumulative sum
	std::vector<double> rc_t(n);
	// area vector
	std::vector<double> a_t(n);
	// cumulative area vector
	std::vector<double> ac_t(n);

	t.at(0) = dist_min;
	a_t.at(0) = 0;
	r_t.at(0) = (a*b) / sqrt(pow(b*cos(t.at(0)), 2) + pow(a*sin(t.at(0)), 2));
	for (int i = 1; i < n; i++)
	{
		// assign the angle
		t.at(i) = (p_i + i*div);
		// get the radius
		r_t.at(i) = (a*b) / sqrt(pow(b*cos(t.at(i)), 2) + pow(a*sin(t.at(i)), 2));
		// get the area
		a_t.at(i) = (r_t.at(i)*r_t.at(i - 1))*sin(div);
	}
	// get the cumulative sum
	std::partial_sum(a_t.begin(), a_t.end(), ac_t.begin());
	// divide cumulative sum by sum
	std::transform(ac_t.begin(), ac_t.end(), ac_t.begin(),
		std::bind(std::divides<double>(), std::placeholders::_1, ac_t.at(n - 1)));
	// construct uniform distribution
	std::uniform_real_distribution<double> ud = std::uniform_real_distribution<double>(0.0, 1.0);
	// get a random number

	double rn = ud(angio_elem->_rengine);

	// find the rc_t value closest to the random number and get the position
	int fi = std::distance(ac_t.begin(), std::lower_bound(ac_t.begin(), ac_t.end(), rn));
	return t.at(fi);
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
	// get the material point data
	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	mat3d F = emp->m_F;
	mat3d R = emp->LeftStretchInverse()*F;
	// get the normed SPD
	mat3d A = (3.0 / initial_angioSPD.tr())*initial_angioSPD;
	// get the spd magnitudes b
	double b[3];
	b[0] = A.col(0).norm(); b[1] = A.col(1).norm(); b[2] = A.col(2).norm();
	// get the rotated directions in Q
	mat3d Q;
	Q.setCol(0, R*A.col(0)/b[0]);
	Q.setCol(1, R*A.col(1)/b[1]);
	Q.setCol(2, R*A.col(2)/b[2]);
	// get the stretch of each axis for d
	double d[3];
	d[0] = (F*Q.col(0)).norm();
	d[1] = (F*Q.col(1)).norm();
	d[2] = (F*Q.col(2)).norm();
	// Assemble the diagonal matrix
	mat3dd D = mat3dd(b[0]*d[0],b[1]*d[1],b[2]*d[2]);
	// Construct the SPD from the eigen components
	mat3d P = Q*D*Q.inverse();
	// input to SPD
	angioSPD = mat3ds(P[0][0], P[1][1], P[2][2], P[0][1], P[1][2], P[0][2]);
	angioSPD = (3.0 / angioSPD.tr())*angioSPD;	
}