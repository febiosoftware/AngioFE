#pragma once
#include <FECore/FEMaterial.h>



class FEAngio;
class AngioElement;
class Tip;

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
protected:
	//! how far the stress policy lets vessels effect the stress within angio elements
	double radius = 1000;
};

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
	DECLARE_PARAMETER_LIST();
private:
	double scale=0.0;
	double y0= -0.004, x0 =2, a= 1.0081, b=0.5436;
	double sprout_mag = 3.72e-12;
	double fan_exponential = 2;
	double sprout_range= 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
	//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

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
	DECLARE_PARAMETER_LIST();
private:
	double sprout_mag = 3.72e-12;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below ~5% of a tip on top of a gauss point
};
//used to update the old code
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
	DECLARE_PARAMETER_LIST();
private:
	double sprout_mag = 3.72e-12;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

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
	DECLARE_PARAMETER_LIST();
private:
	double sprout_mag = 3.72e-12;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

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
	DECLARE_PARAMETER_LIST();
private:
	double sprout_mag = 3.72e-12;
	double fan_exponential = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

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
