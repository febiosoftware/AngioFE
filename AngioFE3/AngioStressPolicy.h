#pragma once
#include <FECore/FEMaterial.h>



class FEAngio;
class AngioElement;
class Tip;

//! base class for calculating the stress from the vascular network
class AngioStressPolicy : public FEMaterial
{
public:
	//! constructor
	explicit AngioStressPolicy(FEModel* pfem) : FEMaterial(pfem) {}
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
	//! how far the stress policy lets vessels effect the stress within angio elements
	double radius = 1000;
	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);
};

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
	//! calculates the stress at the gauss point for a given element, this is scaled by a sigmoid which is used for legacy reasons
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double scale=0.027;
	double y0= -0.004, x0 =2, a= 1.0081, b=0.5436;
	double sprout_mag = 0.02;
	double fan_exponential = 2;
	double sprout_range= 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
	//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
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
	//! calculates the stress at the gauss point for a given element, this is scaled by a sigmoid which is used for legacy reasons
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double scale = 1.0;
	double y0 = 0.0, x0 = 7, a = 1.0, b = 0.5435;
	double sprout_mag = 0.0252;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
	//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

//! a stress calculation at the tips which includes velocity
class LoadCurveVelAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveVelAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveVelAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve and velocity
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below ~5% of a tip on top of a gauss point
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
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

class LoadCurveDenAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveDenAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveDenAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

class LoadCurveRefDenAngioStressPolicy : public AngioStressPolicy
{
public:
	//! constructor
	explicit LoadCurveRefDenAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveRefDenAngioStressPolicy() {}
	//! performs initialization
	bool Init() override;
	//! update scale factor on a per timestep basis
	void UpdateScale() override;
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

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
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve and velocity, this is based on all tips not just active tips
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
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
	//! calculates the stress at the gauss point for a given element, this is scaled by a load curve, this is based on all tips not just active tips
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
protected:
	//! parameter list
	DECLARE_FECORE_CLASS();
private:
	double sprout_mag = 0.027;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
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
	//! calculates the stress at the gauss point for a given element,  this does nothing
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override {}
};
