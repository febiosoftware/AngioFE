#include "AngioElement.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "Tip.h"
#include <iostream>
#include <math.h>
#include <numeric>
#include <algorithm>

#ifndef PI
#define PI 3.14159265358979
#endif

double AngioElement::GetLengthAtTime(FEMesh* mesh, double time) const
{
	double len = 0.0;
	for(int i=0; i < grown_segments.size(); i++)
	{
		len += grown_segments[i]->LengthAtTime(mesh, time);
	}
	return len;
}


void AngioElement::UpdateAngioFractionalAnisotropy()
{
	//sort the spd
	std::vector<pair<double, int>> v;
	double d[3]; vec3d r[3];
	angioSPD.eigen2(d,r);
	v.push_back(pair<double, int>(d[0], 0));
	v.push_back(pair<double, int>(d[1], 1));
	v.push_back(pair<double, int>(d[2], 2));
	sort(v.begin(), v.end(), sortinrev);
	// determine the "eigenvalues"
	double a = v[0].first;
	double b = v[1].first;
	double c = v[2].first;

	// calculate the fractional anisotropy
	angioFA = sqrt(0.5)*(sqrt(pow(a - b, 2) + pow(b - c, 2) + pow(c - a, 2))/(sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))));
}

void AngioElement::UpdateSPD()
{
	mat3d F; F.zero();
	for (int j = 0; j < _elem->GaussPoints(); j++)
	{
		FEMaterialPoint * mp = _elem->GetMaterialPoint(j);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		//create a matrix with the spd on the diagonal
		F += emp->m_F;
	}
	F /= _elem->GaussPoints();
	mat3d spd_d;
	// Apply the deformation gradient to each principal axis (column)
	double d[3]; vec3d r[3];
 	initial_angioSPD.eigen(d, r);
	spd_d.setCol(0, F*d[0]*r[0]);
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
	angioSPD = mat3ds(spd[0][0],spd[1][1],spd[2][2],spd[0][1],spd[1][2],spd[0][2]);
	angioSPD = angioSPD*(3 / angioSPD.tr());
}