#pragma once
#include <FECore/FEMaterial.h>

//! forward declarations
class FEAngio;
class AngioElement;
class Tip;

//! base class for calculating the stress from the vascular network
class AngioStressPolicy : public FEMaterialProperty
{
	FECORE_BASE_CLASS(AngioStressPolicy)
public:
	//! constructor
	explicit AngioStressPolicy(FEModel* pfem) : FEMaterialProperty(pfem) {}
	virtual ~AngioStressPolicy() {}
	//! calculates the stress at the gauss point for a given element
	virtual void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) = 0;
	//! update scale factor on a per timestep basis
	virtual void UpdateScale() = 0;
	//! helper function for getting values out of a load curve at a given time
	void UpdateToLoadCurve(const char* param_name, double & value);
	//! get the density scale
	double GetDensScale(AngioElement* angio_element, Tip* tip, FEMesh* mesh, const int is_ref);
protected:
	//! maximum distance to apply a single stress field (in um)
	double radius = 1000.0;
	//! density scale factor
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

//! basic sigmoid stress policies
//! legacy stress calculations
class SigmoidAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit SigmoidAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~SigmoidAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sigmoid curve parameters
	double scale=0.027;
	double y0 = -0.004, x0 = 2.0, a = 1.0081, b=0.5436;
	//! sprout force field parameters
	double sprout_mag = 0.02;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! legacy stress calculations
class SigmoidDensAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit SigmoidDensAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~SigmoidDensAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sigmoid parameters
	double scale = 1.0;
	double y0 = 1.0e-8, x0 = 7.0, a = 1.0, b = 0.5435;
	//! sprout body force parameters
	double sprout_mag = 0.0252;		
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! Loadcurve based stress policies
//! a stress calculation at the tips which includes velocity
class LoadCurveVelAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveVelAngioStressPolicy(FEModel* pfem) 
		: AngioStressPolicy(pfem) {}
	virtual ~LoadCurveVelAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! stress at the tips. Most trusted stress policy
class LoadCurveAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

class LoadCurveDenAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveDenAngioStressPolicy(FEModel* pfem) 
		: AngioStressPolicy(pfem) {}
	virtual ~LoadCurveDenAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

class LoadCurveRefDenAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveRefDenAngioStressPolicy(FEModel* pfem) 
		: AngioStressPolicy(pfem) {}
	virtual ~LoadCurveRefDenAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! Grown segments policies. These calculate stress at all segments, not just tips.
//! stress calculated from the segments, includes velocity
class GrownSegmentsVelAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit GrownSegmentsVelAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~GrownSegmentsVelAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! stress calculations from the segments.
class GrownSegmentsAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit GrownSegmentsAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~GrownSegmentsAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	DECLARE_FECORE_CLASS()
private:
	//! sprout body force parameters
	double sprout_mag = 0.027;
	double fan_exponential = 2.0;
	//! maximum distance to apply a single stress field (in um)
	double sprout_range = 200.0;
	double sprout_radius_multiplier = 3.0;
};

//! a stress calculation that does nothing
class NullAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit NullAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~NullAngioStressPolicy() {}
	//! performs initialization
	bool Init() override { return true; }
	//! update scale factor on a per timestep basis
	void UpdateScale() override {}
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override {}
};
