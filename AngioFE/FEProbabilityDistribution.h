#pragma once
#include "StdAfx.h"
#include <FECore/FEMaterial.h>
#include <numeric>
#include<FEBioMech/geodesic.h>
#include "FECore/FEModel.h"
#include "FECore/FELoadCurve.h"
#include <FECore/FECoreBase.h>

#ifndef PI
#define PI 3.14159265358979
#endif

//! pure virtual base class for probability distributions
class FEProbabilityDistribution : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FEProbabilityDistribution)
public:
	//! constructor
	explicit FEProbabilityDistribution(FEModel* pfem) : FEMaterialProperty(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	virtual double NextValue(angiofe_random_engine & re) = 0;
	//! generates the next value in the given sequence for a given distribution
	virtual vec3d NextVec(angiofe_random_engine & re) = 0;
	//! updates the distribution to a given time
	virtual void TimeStepUpdate(double current_time) = 0;
protected:
	//! number of times the distribution is sampled to attempt to get a valid 
	//! value before the nextvalue returns nan
	int max_retries = 10;
	//! load curves that are used for distributions must use a step functions
	void SetLoadCurveToStep(const char * param);
	//! returns whether a parameter has changed over time
	bool ChangeInParam(const char * param, double time, double & prev, double & new_p);
	DECLARE_FECORE_CLASS()
};

//! get random numbers from a normal distribution
class FENormalDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FENormalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries 
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double mean = 1.0;//distribution's mean
	double stddev = 1.0;//distribution's standard deviation
	std::normal_distribution<double> nd;
};

//! get random numbers from a uniform distribution
class FEUniformDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEUniformDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double a = 0.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation
	std::uniform_real_distribution<double> rd;
	bool time_clamped = true;
};

//! get random numbers from a exponential distribution
//! http://www.cplusplus.com/reference/random/exponential_distribution/
class FEExponentialDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEExponentialDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double lambda = 1.0;//distribution's lambda
	double mult = 1.0;//multiplier
	std::exponential_distribution<double> ed;
	double prev_lambda = 1.0;
	double prev_mult = mult;
};

//! get random numbers from a cauchy distribution
//! https://www.cplusplus.com/reference/random/cauchy_distribution/
class FECauchyDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FECauchyDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double a = 1.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation
	std::cauchy_distribution<double> cd;
	double prev_a = a;
	double prev_b = b;
};

//! get random numbers from a chisquared distribution
//! https://www.cplusplus.com/reference/random/chi_squared_distribution/
class FEChiSquaredDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEChiSquaredDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double dof = 1.0;//distribution's x^2
	double mult = 1.0;
	std::chi_squared_distribution<double> cd;
	double prev_dof = dof;
	double prev_mult = mult;
};

//! get random numbers from a weibull distribution
//! https://www.cplusplus.com/reference/random/weibull_distribution/
class FEWeibullDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEWeibullDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double a = 1.0;
	double b = 1.0;
	std::weibull_distribution<double> wd;
	double prev_a = a;
	double prev_b = b;
};

//! get random numbers from a gamma distribution
//! http://www.cplusplus.com/reference/random/gamma_distribution/
class FEGammaDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEGammaDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double alpha = 1.0;
	double beta = 1.0;
	std::gamma_distribution<double> gd;
	double prev_alpha = alpha;
	double prev_beta = beta;
};

//! return a given value
class FEFixedDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEFixedDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double value = 1.0;
};

//! get random numbers from an elliptical distribution
class FEEllipticalDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEEllipticalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	vec3d NextVec(angiofe_random_engine & re) override;
	double NextValue(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
	mat3ds spd = mat3ds();
	double efd_exp = 1;
protected:
	DECLARE_FECORE_CLASS()
private:
	vec3d rv;
	mat3d Q = mat3d(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
	double d[3];
	vec3d v[3];
};

//! get random numbers from a Fisherdistribution
class FEFisherDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEFisherDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	vec3d NextVec(angiofe_random_engine & re) override;
	double NextValue(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
	vec3d mu;
	double k;
	static const int resolution = 320;
	double theta[resolution];
	double phi[resolution];
	double area[resolution];
	vec3d dir[resolution];
protected:
	DECLARE_FECORE_CLASS()
private:
	std::vector<double> ODF;
	std::vector<double> cdf;
	std::vector<double> bins;
	std::uniform_real_distribution<double> ud 
		= std::uniform_real_distribution<double>(-1.0, 1.0);
};

//! return a given value. Can be used to set data to a load curve as well?
class FEPrescribedDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FEPrescribedDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;
	vec3d NextVec(angiofe_random_engine & re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double distribution = 1.0;
	// construct uniform distribution
	std::uniform_real_distribution<double> ud 
		= std::uniform_real_distribution<double>(0.0, 1.0);
	std::vector<vec2d> prescribed_distribution;
	int n = -1;
	std::vector<double> pdf; 
	std::vector<double> cdf; 
	std::vector<double> bins; 
};

//! get random numbers from an rational distribution
class FERationalDistribution : public FEProbabilityDistribution
{
public:
	//! constructor
	explicit FERationalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}
	//! generates the next value in the given sequence for a given distribution
	//! this value cannot be zero or less if the value is zero or less the result 
	//! will be redrawn up to max_retries
	//! nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine& re) override;
	vec3d NextVec(angiofe_random_engine& re) override;
	//! performs initialization
	bool Init() override;
	//! updates the distribution to a given time
	void TimeStepUpdate(double current_time) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	double scale = 1.0;
	double p0 = 0.106;
	double p1 = 0.6905;
	double q0 = -89.14;
	double q1 = 2369.0;
	double max_initial_length = 800.0;
	double min_initial_length = 50.0;
	int points = 100;
	double distribution = 1.0;
	// construct uniform distribution
	std::uniform_real_distribution<double> ud 
		= std::uniform_real_distribution<double>(0.0, 1.0);
	std::vector<vec2d> rational_distribution;
	int n = -1;
	std::vector<double> pdf;
	std::vector<double> cdf;
	std::vector<double> bins;
};