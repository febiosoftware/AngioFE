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
	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	//create a matrix with the spd on the diagonal
	mat3d F = emp->m_F;
	
	mat3d spd_d;
	// Apply the deformation gradient to each principal axis (column)
	double d[3]; vec3d r[3];
	this->initial_angioSPD.eigen(d, r);
	spd_d.setCol(0, F*d[0] * r[0]);
	spd_d.setCol(1, F*d[1] * r[1]);
	spd_d.setCol(2, F*d[2] * r[2]);
	// get the vectors of the principal directions and sort in descending order
	std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(spd_d.col(0).norm(), 0));
	v.push_back(pair<double, int>(spd_d.col(1).norm(), 1));
	v.push_back(pair<double, int>(spd_d.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	//transform the 2nd and 3rd PDs so that they are orthogonal to the first PD.
	// set the largest direction to be equal to the stretched spd
	mat3d spd;
	spd.setCol(i, spd_d.col(i));

	// find the perpendicular projection of the 2nd PD onto the 1st PD.
	vec3d u = spd_d.col(j) - spd_d.col(i)*((spd_d.col(j)*spd_d.col(i)) / spd_d.col(i).norm2());
	spd.setCol(j, u);

	// get the unit vector along the new 3rd direction
	// find the 3rd orthogonal direction by getting the cross product of the 1st and 2nd PDs
	vec3d n3 = spd.col(i) ^ spd.col(j); n3.unit();

	// project the 3rd PD onto the cross product direction.
	vec3d w = n3*((spd_d.col(k) * n3) / n3.norm2());
	spd.setCol(k, w);
	angioSPD = mat3ds(spd[0][0], spd[1][1], spd[2][2], spd[0][1], spd[1][2], spd[0][2]);
}