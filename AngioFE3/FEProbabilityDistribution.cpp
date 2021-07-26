#include "FEProbabilityDistribution.h"
#include "FECore/FEModel.h"
#include "FECore/FELoadCurve.h"
#include <FECore/FELoadCurve.h>
#include <numeric>
#include <FECore/FECoreBase.h>

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


BEGIN_FECORE_CLASS(FEProbabilityDistribution, FEMaterial)
ADD_PARAMETER(max_retries, "max_retries");
END_FECORE_CLASS();

void FEProbabilityDistribution::SetLoadCurveToStep(const char * param)
{
	//if load curves are used they must use step interpolation
	FEParam * m = FindParameter(ParamString(param));
	assert(m);
	
	FEModel * model = GetFEModel();
	FELoadCurve* mlc = dynamic_cast<FELoadCurve*>(model->GetLoadController(m));
	assert(mlc);
	if (mlc)
	{
		mlc->SetInterpolation(FEPointFunction::STEP);
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
double FEEllipticalDistribution::NextValue(angiofe_random_engine & re)
{
	// get a random number
	double rn = ud(re);
	// find the rc_t value closest to the random number and get the position
	int fi = std::distance(lc_t.begin(), std::lower_bound(lc_t.begin(), lc_t.end(), rn));
	double val = t.at(fi);
	return val;
}

bool FEEllipticalDistribution::Init()
{
	double div = (PI / 180)*((n + 1) / n);
	// theta
	t.resize(n);
	// radius
	r_t.resize(n);
	// radius cumulative sum
	rc_t.resize(n);
	// area vector
	l_t.resize(n);
	// cumulative area vector
	lc_t.resize(n);

	t.at(0) = p_i;
	l_t.at(0) = 0;
	r_t.at(0) = (a*b) / sqrt(pow(b*cos(t.at(0)), 2) + pow(a*sin(t.at(0)), 2));
	for (int i = 1; i < n; i++)
	{
		// assign the angle
		t.at(i) = (p_i + i*div);
		// get the radius
		r_t.at(i) = (a*b) / sqrt(pow(b*cos(t.at(i)), 2) + pow(a*sin(t.at(i)), 2));
		// get the length
		l_t.at(i) = sqrt(pow(r_t.at(i), 2.0) + pow(r_t.at(i - 1.0), 2.0) - 2.0 * r_t.at(i)*r_t.at(i - 1.0)*cos(div));
	}
	// get the cumulative sum
	std::partial_sum(l_t.begin(), l_t.end(), lc_t.begin());
	// divide cumulative sum by sum
	std::transform(lc_t.begin(), lc_t.end(), lc_t.begin(),
		std::bind(std::divides<double>(), std::placeholders::_1, lc_t.at(n - 1)));


	//gd = std::gamma_distribution<double>(a, b);
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");

	return true;
}

void FEEllipticalDistribution::TimeStepUpdate(double current_time)
{

}



BEGIN_FECORE_CLASS(FEEllipticalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, "a");
ADD_PARAMETER(b, "b");
END_FECORE_CLASS();