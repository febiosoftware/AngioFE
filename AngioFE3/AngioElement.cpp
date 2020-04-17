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

/*bool sortinrev(const pair<double, int> &a, const pair<double, int> &b)
{
	return (a.first > b.first);
}*/

double AngioElement::GetEllipseAngle(const double a, const double b)
{
	// initialize vals
	double p_i = -90;
	double p_f = PI / 2;
	int n = 180;
	double div = (PI/180)*(181/180);

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
		std::bind(std::divides<double>(), std::placeholders::_1, ac_t.at(n-1)));
	// construct uniform distribution
	std::uniform_real_distribution<double> ud = std::uniform_real_distribution<double>(0.0, 1.0);
	// get a random number
	double rn = ud(_rengine);
	// find the rc_t value closest to the random number and get the position
	int i = std::distance(ac_t.begin(), std::lower_bound(ac_t.begin(), ac_t.end(), rn));
	return t.at(i);
}

void AngioElement::UpdateAngioFractionalAnisotropy()
{
	//sort the spa
	std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(angioSPA.col(0).norm(), 0));
	v.push_back(pair<double, int>(angioSPA.col(1).norm(), 1));
	v.push_back(pair<double, int>(angioSPA.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// determine the "eigenvalues"
	double a = v[0].first;
	double b = v[1].first;
	double c = v[2].first;

	// calculate the fractional anisotropy
	angioFA = sqrt(0.5)*(sqrt(pow(a - b, 2) + pow(b - c, 2) + pow(c - a, 2))/(sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2))));
}

void AngioElement::UpdateSPA()
{
	mat3d F; F.zero();
	for (int j = 0; j < _elem->GaussPoints(); j++)
	{
		FEMaterialPoint * mp = _elem->GetMaterialPoint(j);
		FEElasticMaterialPoint* emp = mp->ExtractData<FEElasticMaterialPoint>();
		//create a matrix with the spa on the diagonal
		F += emp->m_F;
	}
	F /= _elem->GaussPoints();
	mat3d spa_d;
	// Apply the deformation gradient to each principal axis (column)
	spa_d.setCol(0, F*initial_angioSPA.col(0));
	spa_d.setCol(1, F*initial_angioSPA.col(1));
	spa_d.setCol(2, F*initial_angioSPA.col(2));
	// get the vectors of the principal directions and sort in descending order
	std::vector<pair<double, int>> v;
	v.push_back(pair<double, int>(spa_d.col(0).norm(), 0));
	v.push_back(pair<double, int>(spa_d.col(1).norm(), 1));
	v.push_back(pair<double, int>(spa_d.col(2).norm(), 2));
	sort(v.begin(), v.end(), sortinrev);

	// store the indices
	int i = v[0].second;
	int j = v[1].second;
	int k = v[2].second;

	//transform the 2nd and 3rd PDs so that they are orthogonal to the first PD.
	// set the largest direction to be equal to the stretched spa
	mat3d spa;
	spa.setCol(i, spa_d.col(i));

	// find the perpendicular projection of the 2nd PD onto the 1st PD.
	vec3d u = spa_d.col(j) - spa_d.col(i)*((spa_d.col(j)*spa_d.col(i)) / spa_d.col(i).norm2());
	spa.setCol(j, u);

	// get the unit vector along the new 3rd direction
	// find the 3rd orthogonal direction by getting the cross product of the 1st and 2nd PDs
	vec3d n3 = spa.col(i) ^ spa.col(j); n3.unit();

	// project the 3rd PD onto the cross product direction.
	vec3d w = n3*((spa_d.col(k) * n3) / n3.norm2());
	spa.setCol(k, w);
	angioSPA = spa;
}