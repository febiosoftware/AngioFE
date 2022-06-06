#include "FEProbabilityDistribution.h"
#include "FECore/FEModel.h"
#include "FECore/FELoadCurve.h"
#include <FECore/FELoadCurve.h>
#include <numeric>
#include <FECore/FECoreBase.h>
#include <iostream>
#include <algorithm>

BEGIN_FECORE_CLASS(FEProbabilityDistribution, FEMaterialProperty)
ADD_PARAMETER(max_retries, "max_retries");
END_FECORE_CLASS();

void FEProbabilityDistribution::SetLoadCurveToStep(const char * param)
{
	//if load curves are used they must use step interpolation
	FEParam * m = FindParameter(ParamString(param));
	assert(m);

	FEModel * model = GetFEModel();
	FELoadCurve* mlc = dynamic_cast<FELoadCurve*>(model->GetLoadController(m));
	//assert(mlc);
	if (mlc)
	{
		mlc->SetInterpolation(PointCurve::STEP);
	}
}

bool FEProbabilityDistribution::ChangeInParam(const char * param, double time, double & prev, double & new_p)
{
	FEParam * m = FindParameter(ParamString(param));
	assert(m);

	FEModel * model = GetFEModel();
	FELoadCurve* mlc = dynamic_cast<FELoadCurve*>(model->GetLoadController(m));
	assert(mlc);
	if (mlc)
	{
		new_p = mlc->GetValue(time);
		if (new_p != prev)
		{
			return true;
		}
	}
	return false;
}

//implemenations of FENormalDistribution
double FENormalDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = nd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FENormalDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FENormalDistribution::Init()
{
	nd = std::normal_distribution<double>(mean, stddev);

	return true;
}

void FENormalDistribution::TimeStepUpdate(double current_time)
{
	nd = std::normal_distribution<double>(mean, stddev);
}

BEGIN_FECORE_CLASS(FENormalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(mean, "mean");
ADD_PARAMETER(stddev, "stddev");
END_FECORE_CLASS();

double FEUniformDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = rd(re);
		if (val >= 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FEUniformDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FEUniformDistribution::Init()
{
	rd = std::uniform_real_distribution<double>(a,b);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");

	return true;
}

void FEUniformDistribution::TimeStepUpdate(double current_time)
{
	if(time_clamped)
	{
		rd = std::uniform_real_distribution<double>(a, b - current_time);
	}
	else
	{
		rd = std::uniform_real_distribution<double>(a, b);
	}
}

BEGIN_FECORE_CLASS(FEUniformDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
ADD_PARAMETER(time_clamped, "time_clamped");
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FEExponentialDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = ed(re);
		if (val > 0.0)
			return mult*val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FEExponentialDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FEExponentialDistribution::Init()
{
	ed = std::exponential_distribution<double>(lambda);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("lambda");
	SetLoadCurveToStep("mult");
	return true;
}

void FEExponentialDistribution::TimeStepUpdate(double current_time)
{
	ed = std::exponential_distribution<double>(lambda);
}

BEGIN_FECORE_CLASS(FEExponentialDistribution, FEProbabilityDistribution)
ADD_PARAMETER(lambda, "lambda");
ADD_PARAMETER(mult, "mult");
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FECauchyDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = cd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FECauchyDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FECauchyDistribution::Init()
{
	cd = std::cauchy_distribution<double>(a, b);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");

	return true;
}
void FECauchyDistribution::TimeStepUpdate(double current_time)
{
	cd = std::cauchy_distribution<double>(a, b);
}

BEGIN_FECORE_CLASS(FECauchyDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FEChiSquaredDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = cd(re);
		if (val > 0.0)
			return mult*val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FEChiSquaredDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FEChiSquaredDistribution::Init()
{
	cd = std::chi_squared_distribution<double>(dof);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("dof");
	SetLoadCurveToStep("mult");

	return true;
}

void FEChiSquaredDistribution::TimeStepUpdate(double current_time)
{
	cd = std::chi_squared_distribution<double>(dof);
}

BEGIN_FECORE_CLASS(FEChiSquaredDistribution, FEProbabilityDistribution)
ADD_PARAMETER(dof, "dof");
ADD_PARAMETER(mult, "mult");
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FEWeibullDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = wd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FEWeibullDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FEWeibullDistribution::Init()
{
	wd = std::weibull_distribution<double>(a, b);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");
	return true;
}

void FEWeibullDistribution::TimeStepUpdate(double current_time)
{
	wd = std::weibull_distribution<double>(a, b);
}

BEGIN_FECORE_CLASS(FEWeibullDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
END_FECORE_CLASS();

//implemenations of FEGammaDistribution
double FEGammaDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = gd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vec3d FEGammaDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

bool FEGammaDistribution::Init()
{
	gd = std::gamma_distribution<double>(alpha, beta);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("alpha");
	SetLoadCurveToStep("beta");

	return true;
}

void FEGammaDistribution::TimeStepUpdate(double current_time)
{
	gd = std::gamma_distribution<double>(alpha, beta);
}

BEGIN_FECORE_CLASS(FEGammaDistribution, FEProbabilityDistribution)
ADD_PARAMETER(alpha, "alpha");
ADD_PARAMETER(beta, "beta");
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FEFixedDistribution::NextValue(angiofe_random_engine & re)
{
	return value;
}

vec3d FEFixedDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1, 0, 0);
}

void FEFixedDistribution::TimeStepUpdate(double current_time)
{
	
}

bool FEFixedDistribution::Init()
{
	return true;
}

BEGIN_FECORE_CLASS(FEFixedDistribution, FEProbabilityDistribution)
ADD_PARAMETER(value, "value");
END_FECORE_CLASS();

//implemenations of FEGammaDistribution
vec3d FEEllipticalDistribution::NextVec(angiofe_random_engine & re)
{
	bool found = false;
	while(!found)
	{
		std::uniform_real_distribution<double> rd = std::uniform_real_distribution<double> (-1,1);
		vec3d rv;
		// get the random direction with magnitude up to the value of r^efd_exp
		double w = std::max(std::max(d[0], d[1]), d[2]);
		rv.x = pow(w,efd_exp)*rd(re);
		rv.y = pow(w,efd_exp)*rd(re);
		rv.z = pow(w,efd_exp)*rd(re);
		// Find the value of the standard ellipse equation
		double ellipse_par_eq = (rv.x * rv.x) / (d[0] * d[0]) + (rv.y * rv.y) / (d[1] * d[1]) + (rv.z * rv.z) / (d[2] * d[2]);
		// Find the radius to the random vector point
		double r = rv.norm();
		// Find the max value which is equal to r^efd_exp for the given orientation.
		double r_max = sqrt(1 / ellipse_par_eq) * r;
		// If the sampled point is within the power elliptical distribution
		if (r < pow(r_max,efd_exp)) {
			// rotate the random direction from the global basis into the EFD basis and normalize
			rv = Q*rv;
			rv = rv.normalized();
			found = true;
			return rv;
		}
	}	
}

double FEEllipticalDistribution::NextValue(angiofe_random_engine & re)
{
	return 0;
}

bool FEEllipticalDistribution::Init()
{
	spd.eigen(d, v);
	// set the eigenvector matrix
	Q.setCol(0, v[0]); Q.setCol(1, v[1]); Q.setCol(2, v[2]);
	return true;
}

void FEEllipticalDistribution::TimeStepUpdate(double current_time)
{

}

BEGIN_FECORE_CLASS(FEEllipticalDistribution, FEProbabilityDistribution)
END_FECORE_CLASS();

//implemenations of FEGammaDistribution
vec3d FEFisherDistribution::NextVec(angiofe_random_engine & re)
{
	//cdf.resize(resolution);
	ODF.resize(resolution);
	//bins.resize(resolution);
	double prior_val = 0.0;
	// Create the cdf
	for (int i = 0; i < resolution; i++) {
		//bins.at(i) = i;
		ODF.at(i) = AREAL[i]*(k / (2 * PI*(exp(k) - exp(-k))) * exp(k*(mu*dir[0])));
		//cdf.at(i) = ODF.at(i) + prior_val;
		//prior_val = cdf.at(i);
	}
	vec3d rand_pt;
	double maxODF = *std::max_element(ODF.begin(), ODF.end());
	bool found = false;
	while (!found) {
		rand_pt.x = ud(re)*maxODF;
		rand_pt.y = ud(re)*maxODF;
		rand_pt.z = ud(re)*maxODF;
		std::vector<double> dist;
		for (int i = 0; i < resolution; i++) {
			dist.emplace_back(rand_pt*dir[i]);
		}
		//double rand_dir = *std::min_element(dist.begin(),dist.end());
		int rand_pt_indx = std::min_element(dist.begin(), dist.end()) - dist.begin();
		double check = ODF[rand_pt_indx];
		//if (rand_pt.norm() <= check) {
		//	return rand_pt.unit();
		//}
	}

	// get the cumulative sum
	std::partial_sum(ODF.begin(), ODF.end(), cdf.begin());
	// divide cumulative sum by sum
	std::transform(cdf.begin(), cdf.end(), cdf.begin(),
		std::bind(std::divides<double>(), std::placeholders::_1, cdf.at(resolution - 1)));
	double rn = ud(re);
	// find the value closest to the random number and get the position
	int fi = std::distance(cdf.begin(), std::lower_bound(cdf.begin(), cdf.end(), rn));
	return dir[fi];

}

double FEFisherDistribution::NextValue(angiofe_random_engine & re)
{
	return 0;
}

bool FEFisherDistribution::Init()
{
	ODF.reserve(resolution);
	for (int i = 0; i < resolution; i++) {
		dir[i].x = cos(PHIL[i])*sin(THETAL[i]);	dir[i].y = sin(PHIL[i])*sin(THETAL[i]);	dir[i].z = cos(THETAL[i]);
		ODF[i] = 1.0/resolution;
	}
	return true;
}

void FEFisherDistribution::TimeStepUpdate(double current_time)
{

}

BEGIN_FECORE_CLASS(FEFisherDistribution, FEProbabilityDistribution)
END_FECORE_CLASS();

//implemenations of FENormalDistribution
double FEPrescribedDistribution::NextValue(angiofe_random_engine & re)
{
	// get a random number
	double rn = ud(re);
	// find the rc_t value closest to the random number and get the position
	int fi = std::distance(cdf.begin(), std::lower_bound(cdf.begin(), cdf.end(), rn));
	return bins.at(fi);
}

vec3d FEPrescribedDistribution::NextVec(angiofe_random_engine & re)
{
	return vec3d(1,0,0);
}

void FEPrescribedDistribution::TimeStepUpdate(double current_time)
{

}

// gets the load curve, solves the cdf
bool FEPrescribedDistribution::Init()
{
	FEParam * m = FindParameter(ParamString("distribution"));
	assert(m);

	FEModel * model = GetFEModel();
	FELoadCurve* mlc = dynamic_cast<FELoadCurve*>(model->GetLoadController(m));
	prescribed_distribution = mlc->GetFunction().GetPoints();

	//Create the cdf
	n = prescribed_distribution.size();
	bins.resize(n);
	pdf.resize(n);
	cdf; cdf.resize(n);
	double prior_val = 0.0;
	for (int i = 0; i < cdf.size(); i++)
	{
		bins.at(i) = prescribed_distribution.at(i).x();
		pdf.at(i) = prescribed_distribution.at(i).y();
		cdf.at(i) =  prescribed_distribution.at(i).y() + prior_val;
		prior_val = cdf.at(i);
	}
	// get the cumulative sum
	std::partial_sum(pdf.begin(), pdf.end(), cdf.begin());
	// divide cumulative sum by sum
	std::transform(cdf.begin(), cdf.end(), cdf.begin(),
		std::bind(std::divides<double>(), std::placeholders::_1, cdf.at(n - 1)));
	return true;
}

BEGIN_FECORE_CLASS(FEPrescribedDistribution, FEProbabilityDistribution)
ADD_PARAMETER(distribution, "distribution");
END_FECORE_CLASS();