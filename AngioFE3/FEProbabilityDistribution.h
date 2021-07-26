#pragma once
#include "StdAfx.h"
#include <FECore/FEMaterial.h>
#include <numeric>

#include "FECore/FEModel.h"
#include "FECore/FELoadCurve.h"
#include <FECore/FECoreBase.h>

#ifndef PI
#define PI 3.14159265358979
#endif

//! pure virtual base class for probability distributions
class FEProbabilityDistribution : public FEMaterial
{
public:
	//! constructor
	explicit FEProbabilityDistribution(FEModel* pfem) : FEMaterial(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	virtual double NextValue(angiofe_random_engine & re) = 0;

	//! updates the distribution to a given time
	virtual void TimeStepUpdate(double current_time) = 0;
	
protected:
	//! number of times the distribution is sampled to attempt to get a valid value before the nextvalue returns nan
	int max_retries = 10;

	//! load curves that are used for distributions must use a step functions
	void SetLoadCurveToStep(const char * param);

	//! returns whether a parameter has changed over time
	bool ChangeInParam(const char * param, double time, double & prev, double & new_p);
private:
	DECLARE_FECORE_CLASS();
};

//! get random numbers from a normal distribution
class FENormalDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FENormalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

private:
	double mean = 1.0;//distribution's mean
	double stddev = 1.0;//distribution's standard deviation

	std::normal_distribution<double> nd;
	DECLARE_FECORE_CLASS();
};

//! get random numbers from a uniform distribution
class FEUniformDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEUniformDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

private:
	double a = 0.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation

	std::uniform_real_distribution<double> rd;
	bool time_clamped = true;


	DECLARE_FECORE_CLASS();
};

//! get random numbers from a exponential distribution
//! http://www.cplusplus.com/reference/random/exponential_distribution/
class FEExponentialDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEExponentialDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

private:
	double lambda = 1.0;//distribution's lambda
	double mult = 1.0;//multiplier

	std::exponential_distribution<double> ed;

	double prev_lambda = 1.0;
	double prev_mult = mult;

	DECLARE_FECORE_CLASS();
};

//! get random numbers from a cauchy distribution
//! https://www.cplusplus.com/reference/random/cauchy_distribution/
class FECauchyDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FECauchyDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
private:
	double a = 1.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation

	std::cauchy_distribution<double> cd;

	double prev_a = a;
	double prev_b = b;


	DECLARE_FECORE_CLASS();
};

//! get random numbers from a chisquared distribution
//! https://www.cplusplus.com/reference/random/chi_squared_distribution/
class FEChiSquaredDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEChiSquaredDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
private:
	double dof = 1.0;//distribution's x^2
	double mult = 1.0;

	std::chi_squared_distribution<double> cd;

	double prev_dof = dof;
	double prev_mult = mult;
	

	DECLARE_FECORE_CLASS();
};

//! get random numbers from a weibull distribution
//! https://www.cplusplus.com/reference/random/weibull_distribution/
class FEWeibullDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEWeibullDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

private:
	double a = 1.0;
	double b = 1.0;

	std::weibull_distribution<double> wd;

	double prev_a = a;
	double prev_b = b;

	DECLARE_FECORE_CLASS();
};

//! get random numbers from a gamma distribution
//! http://www.cplusplus.com/reference/random/gamma_distribution/
class FEGammaDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEGammaDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

private:
	double alpha = 1.0;
	double beta = 1.0;

	std::gamma_distribution<double> gd;

	double prev_alpha = alpha;
	double prev_beta = beta;

	DECLARE_FECORE_CLASS();
};

//! return a given value
class FEFixedDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEFixedDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	
	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
private:
	double value = 1.0;

	DECLARE_FECORE_CLASS();
};

//! get random numbers from an elliptical distribution
class FEEllipticalDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEEllipticalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//! generates the next value in the given sequence which fits a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	//! performs initialization
	bool Init() override;

	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;

	double a = 1.0;
	double b = 1.0;

private:


	// construct uniform distribution
	std::uniform_real_distribution<double> ud = std::uniform_real_distribution<double>(0.0, 1.0);

	//std::gamma_distribution<double> gd;

	double prev_alpha = a;
	double prev_beta = b;
	double p_i = -PI / 2;
	double p_f = PI / 2;
	const int n = 180;

	// theta
	std::vector<double> t;
	// radius
	std::vector<double> r_t;
	// radius cumulative sum
	std::vector<double> rc_t;
	// area vector
	std::vector<double> l_t;
	// cumulative area vector
	std::vector<double> lc_t;

	DECLARE_FECORE_CLASS();
};