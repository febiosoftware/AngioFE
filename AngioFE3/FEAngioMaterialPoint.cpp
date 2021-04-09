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

	

	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	mat3d F = emp->m_F;
	mat3d V = emp->LeftStretch();
	mat3d R = emp->LeftStretchInverse()*F;
	mat3d AT = initial_angioSPD;
	//AT = R*AT*R.transpose();
	AT.setCol(0, R.transpose()*AT.col(0));
	AT.setCol(1, R.transpose()*AT.col(1));
	AT.setCol(2, R.transpose()*AT.col(2));
	//AT = AT*(3.0 / AT.trace());
	//AT = V*AT;
	double v0 = V.col(0).norm(); double v1 = V.col(1).norm(); double v2 = V.col(2).norm();
	AT.setCol(0, AT.col(0)*v0);
	AT.setCol(1, AT.col(1)*v1);
	AT.setCol(2, AT.col(2)*v2);
	angioSPD = mat3ds(AT[0][0], AT[1][1], AT[2][2], AT[0][1], AT[1][2], AT[0][2]);
	angioSPD = (3.0 / angioSPD.tr())*angioSPD;

	//FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	//mat3d F = emp->m_F;
	//mat3d spd_d;
	//// Apply the deformation gradient to each principal axis (column)
	//double d[3]; vec3d r[3];
	//initial_angioSPD.eigen(d, r);
	//spd_d.setCol(0, F*d[0] * r[0]);
	//spd_d.setCol(1, F*d[1] * r[1]);
	//spd_d.setCol(2, F*d[2] * r[2]);
	//spd_d = spd_d.sym();
	////spd_d = F.sym()*initial_angioSPD;
	//// get the vectors of the principal directions and sort in descending order
	//std::vector<pair<double, int>> v;
	//v.push_back(pair<double, int>(spd_d.col(0).norm(), 0));
	//v.push_back(pair<double, int>(spd_d.col(1).norm(), 1));
	//v.push_back(pair<double, int>(spd_d.col(2).norm(), 2));
	//sort(v.begin(), v.end(), sortinrev);

	//// store the indices
	//int i = v[0].second;
	//int j = v[1].second;
	//int k = v[2].second;

	////transform the 2nd and 3rd PDs so that they are orthogonal to the first PD.
	//// set the largest direction to be equal to the stretched spd
	//mat3d spd;
	//spd.setCol(i, spd_d.col(i));

	//// find the perpendicular projection of the 2nd PD onto the 1st PD.
	//vec3d u = spd_d.col(j) - spd_d.col(i)*((spd_d.col(j)*spd_d.col(i)) / spd_d.col(i).norm2());
	//spd.setCol(j, u);

	//// get the unit vector along the new 3rd direction
	//// find the 3rd orthogonal direction by getting the cross product of the 1st and 2nd PDs
	//vec3d n3 = spd.col(i) ^ spd.col(j); n3.unit();

	//// project the 3rd PD onto the cross product direction.
	//vec3d w = n3*((spd_d.col(k) * n3) / n3.norm2());
	//spd.setCol(k, w);
	//angioSPD = mat3ds(spd[0][0], spd[1][1], spd[2][2], spd[0][1], spd[1][2], spd[0][2]);
	//angioSPD = (3.0 / angioSPD.tr())*angioSPD;

	////  end old part ////


	/*mat3d spd_d;
	spd_d = F.sym()*this->angioSPD;
	double d[3]; vec3d r[3];
	initial_angioSPD.eigen(d, r);
	vec3d p0[3];
	p0[0] = r[0] * d[0];
	p0[1] = r[1] * d[1];
	p0[2] = r[2] * d[2];
	double lambda[3];
	lambda[0] = (F*r[0]).norm();
	lambda[1] = (F*r[1]).norm();
	lambda[2] = (F*r[2]).norm();
	mat3d R; mat3ds U;
	F.right_polar(R, U);
	vec3d p[3];
	p[0] = R*p0[0] * lambda[0];
	p[1] = R*p0[1] * lambda[1];
	p[2] = R*p0[2] * lambda[2];
	double dp[3]; vec3d rp[3];*/



	//spd_d.setCol(0, F * r[0] * d[0]);
	//spd_d.setCol(1, F * r[1] * d[1]);
	//spd_d.setCol(2, F * r[2] * d[2]);

	//spd_d.setCol(0, R * r[0] * d[0]);
	//spd_d.setCol(1, R * r[1] * d[1]);
	//spd_d.setCol(2, R * r[2] * d[2]);

	/*angioSPD = spd_d.sym();*/

	//double d[3]; vec3d r[3];
	//initial_angioSPD.eigen(d, r);
	//spd_d.setCol(0, F * r[0] * d[0]);
	//spd_d.setCol(1, F * r[1] * d[1]);
	//spd_d.setCol(2, F * r[2] * d[2]);
	//// Apply the deformation gradient to each principal axis (column)
	//
	///*spd_d.setCol(0, F*initial_angioSPD.col(0));
	//spd_d.setCol(1, F*initial_angioSPD.col(1));
	//spd_d.setCol(2, F*initial_angioSPD.col(2));*/
	//// get the vectors of the principal directions and sort in descending order
	//std::vector<pair<double, int>> v;
	//v.push_back(pair<double, int>(spd_d.col(0).norm(), 0));
	//v.push_back(pair<double, int>(spd_d.col(1).norm(), 1));
	//v.push_back(pair<double, int>(spd_d.col(2).norm(), 2));
	//sort(v.begin(), v.end(), sortinrev);

	//// store the indices
	//int i = v[0].second;
	//int j = v[1].second;
	//int k = v[2].second;

	////transform the 2nd and 3rd PDs so that they are orthogonal to the first PD.
	//// set the largest direction to be equal to the stretched spa
	//mat3d spd;
	//spd.setCol(i, spd_d.col(i));

	//// find the perpendicular projection of the 2nd PD onto the 1st PD.
	//vec3d u = spd_d.col(j) - spd_d.col(i)*((spd_d.col(j)*spd_d.col(i)) / spd_d.col(i).norm2());
	//spd.setCol(j, u);

	//// get the unit vector along the new 3rd direction
	//// find the 3rd orthogonal direction by getting the cross product of the 1st and 2nd PDs
	//vec3d n3 = spd.col(i) ^ spd.col(j); n3.unit();

	//// project the 3rd PD onto the cross product direction.
	//vec3d w = n3*((spd_d.col(k) * n3) / n3.norm2());
	//spd.setCol(k, w);
	//angioSPD = spd.sym();

	////if (nhit)
	////{
	//	FEElasticMaterialPoint* emp = this->ExtractData<FEElasticMaterialPoint>();
	//	//create a matrix with the spd on the diagonal
	//	mat3d F = emp->m_F;
	//	/*mat3d I = mat3d(1, 0, 0, 0, 1, 0, 0, 0, 1);*/
	//	//mat3d E = (F.transpose()*F - I) / 2;
	//	//mat3d Eh = (F.transpose()*F) / 2;
	//	mat3d spd_d;
	//	/*vec3d r[3];
	//	r[0] = vec3d(this->initial_angioSPD.xx(), this->initial_angioSPD.xy(), this->initial_angioSPD.xz());
	//	r[1] = vec3d(this->initial_angioSPD.xy(), this->initial_angioSPD.yy(), this->initial_angioSPD.yz());
	//	r[2] = vec3d(this->initial_angioSPD.xz(), this->initial_angioSPD.yz(), this->initial_angioSPD.zz());*/
	//	/*r[0] = F.sym()*r[0];
	//	r[1] = F.sym()*r[1];
	//	r[2] = F.sym()*r[2];*/
	//	/*mat3d Fd = F*Fp.inverse();
	//	Fi = Fd.inverse();*/
	//	//Fp = Fd*Fi;
	//	/*r[0] = vec3d(this->angioSPD.xx(), this->angioSPD.xy(), this->angioSPD.xz());
	//	r[1] = vec3d(this->angioSPD.xy(), this->angioSPD.yy(), this->angioSPD.yz());
	//	r[2] = vec3d(this->angioSPD.xz(), this->angioSPD.yz(), this->angioSPD.zz());*/
	//	/*r[0] = Fd.sym()*r[0];
	//	r[1] = Fd.sym()*r[1];
	//	r[2] = Fd.sym()*r[2];*/
	//	/*r[0] = (E + I)*r[0];
	//	r[1] = (E + I)*r[1];
	//	r[2] = (E + I)*r[2];*/
	//	/*r[0] = (Eh)*r[0];
	//	r[1] = (Eh)*r[1];
	//	r[2] = (Eh)*r[2];*/
	//	/*spd_d.setCol(0, r[0]);
	//	spd_d.setCol(1, r[1]);
	//	spd_d.setCol(2, r[2]);
	//	angioSPD = spd_d.sym();*/
	//	//nhit = 0;
	////}
	// Apply the deformation gradient to each principal axis (column)
	//double d[3]; vec3d r[3];
	//this->angioSPD.eigen(d, r);
	/*vec3d r0 = vec3d(this->initial_angioSPD.xx(), this->initial_angioSPD.xy(), this->initial_angioSPD.xz());
	vec3d r1= vec3d(this->initial_angioSPD.xy(), this->initial_angioSPD.yy(), this->initial_angioSPD.yz());
	vec3d r2 = vec3d(this->initial_angioSPD.xz(), this->initial_angioSPD.yz(), this->initial_angioSPD.zz());
	spd_d.setCol(0, F*r0);
	spd_d.setCol(1, F*r1);
	spd_d.setCol(2, F*r2);*/
	// get the vectors of the principal directions and sort in descending order
	//Apply the deformation gradient to each principal axis (column)
	/*double d[3]; vec3d r[3];
	initial_angioSPD.eigen(d, r);
	spd_d.setCol(0, F * r[0] * d[0]);
	spd_d.setCol(1, F * r[1] * d[1]);
	spd_d.setCol(2, F * r[2] * d[2]);*/
	/*vec3d r[3]; 
	r[0] = vec3d(this->angioSPD.xx(), this->angioSPD.xy(), this->angioSPD.xz());
	r[1] = vec3d(this->angioSPD.xy(), this->angioSPD.yy(), this->angioSPD.yz());
	r[2] = vec3d(this->angioSPD.xz(), this->angioSPD.yz(), this->angioSPD.zz());*/
	/*spd_d.setCol(0, r[0]);
	spd_d.setCol(1, r[1]);
	spd_d.setCol(2, r[2]);*/
	// Multiply the spd by e+I 
	//spd_d = (F+F.transpose())*0.5*spd_d;
	//spd_d = F.sym()*spd_d;
	/*std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(spd_d.col(0).norm(), 0));
	v.push_back(pair<double, int>(spd_d.col(1).norm(), 1));
	v.push_back(pair<double, int>(spd_d.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);*/

	//// store the indices
	/*int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;*/

	////transform the 2nd and 3rd PDs so that they are orthogonal to the first PD.
	//// set the largest direction to be equal to the stretched spd
	/*mat3d spd;
	spd.setCol(i, spd_d.col(i));*/

	//// find the perpendicular projection of the 2nd PD onto the 1st PD.
	/*-vec3d u = spd_d.col(j) - spd_d.col(i)*((spd_d.col(j)*spd_d.col(i)) / spd_d.col(i).norm2());
	spd.setCol(j, u);

	//// get the unit vector along the new 3rd direction
	//// find the 3rd orthogonal direction by getting the cross product of the 1st and 2nd PDs
	vec3d n3 = spd.col(i) ^ spd.col(j); n3.unit();

	//// project the 3rd PD onto the cross product direction.
	vec3d w = n3*((spd_d.col(k) * n3) / n3.norm2());
	spd.setCol(k, w);
	angioSPD = mat3ds(spd[0][0], spd[1][1], spd[2][2], spd[0][1], spd[1][2], spd[0][2]);
	//angioSPD = (3.0 / angioSPD.tr())*angioSPD;*/
}