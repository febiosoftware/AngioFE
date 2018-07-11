#include "FEProbabilityDistribution.h"
#include "FECore/FEModel.h"
#include "FECore/LoadCurve.h"
#include <FECore/FEDataLoadCurve.h>

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


BEGIN_PARAMETER_LIST(FENormalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(mean, FE_PARAM_DOUBLE, "mean");
ADD_PARAMETER(stddev, FE_PARAM_DOUBLE, "stddev");
END_PARAMETER_LIST();


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

BEGIN_PARAMETER_LIST(FEUniformDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
ADD_PARAMETER(time_clamped, FE_PARAM_BOOL, "time_clamped");
END_PARAMETER_LIST();


BEGIN_PARAMETER_LIST(FEProbabilityDistribution, FEMaterial)
ADD_PARAMETER(max_retries, FE_PARAM_INT, "max_retries");
END_PARAMETER_LIST();

void FEProbabilityDistribution::SetLoadCurveToStep(const char * param)
{
	//if load curves are used they must use step interpolation
	FEParam * m = FindParameter(ParamString(param));

	int mlci = m->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci > 0)
	{
		FEDataLoadCurve * mlc = dynamic_cast<FEDataLoadCurve*>(model->GetLoadCurve(mlci));
		assert(mlc);
		mlc->SetInterpolation(FEDataLoadCurve::STEP);
	}
}

bool FEProbabilityDistribution::ChangeInParam(const char * param, double time, double & prev, double & new_p)
{
	FEParam * m = FindParameter(ParamString(param));

	int mlci = m->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci > 0)
	{
		FELoadCurve * mlc = model->GetLoadCurve(mlci);
		new_p = mlc->Value(time);
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


BEGIN_PARAMETER_LIST(FEExponentialDistribution, FEProbabilityDistribution)
ADD_PARAMETER(lambda, FE_PARAM_DOUBLE, "lambda");
ADD_PARAMETER(mult, FE_PARAM_DOUBLE, "mult");
END_PARAMETER_LIST();

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



BEGIN_PARAMETER_LIST(FECauchyDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();


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



BEGIN_PARAMETER_LIST(FEChiSquaredDistribution, FEProbabilityDistribution)
ADD_PARAMETER(dof, FE_PARAM_DOUBLE, "dof");
ADD_PARAMETER(mult, FE_PARAM_DOUBLE, "mult");
END_PARAMETER_LIST();


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



BEGIN_PARAMETER_LIST(FEWeibullDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();


//implemenations of FENormalDistribution
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



BEGIN_PARAMETER_LIST(FEGammaDistribution, FEProbabilityDistribution)
ADD_PARAMETER(alpha, FE_PARAM_DOUBLE, "alpha");
ADD_PARAMETER(beta, FE_PARAM_DOUBLE, "beta");
END_PARAMETER_LIST();

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


BEGIN_PARAMETER_LIST(FEFixedDistribution, FEProbabilityDistribution)
ADD_PARAMETER(value, FE_PARAM_DOUBLE, "value");
END_PARAMETER_LIST();