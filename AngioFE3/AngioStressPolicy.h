#pragma once
#include <FECore/FEMaterial.h>



class FEAngio;
class AngioElement;
class Tip;

class AngioStressPolicy : public FEMaterial
{
public:
	explicit AngioStressPolicy(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~AngioStressPolicy() {}
	//will set m_as for all gauss points of the element
	virtual void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) = 0;
	virtual void UpdateScale(double time) = 0;
	static void GetActiveFinalTipsInRadius(AngioElement* angio_element, double radius, FEAngio* pangio, std::vector<Tip *> & tips);
protected:
	double radius = 1000;
};

class SigmoidAngioStressPolicy : public AngioStressPolicy
{
public:
	explicit SigmoidAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~SigmoidAngioStressPolicy() {}
	bool Init() override;
	void UpdateScale(double time) override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
	
private:
	double scale=0.0;
	double y0= -0.004, x0 =2, a= 1.0081, b=0.5436;
	double sprout_mag = 3.72e-12;
	double sprout_width = 2;
	double sprout_range= 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
	//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};
//probably the right way to do this calculation
class LoadCurveVelAngioStressPolicy : public AngioStressPolicy
{
public:
	explicit LoadCurveVelAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveVelAngioStressPolicy() {}
	bool Init() override;
	void UpdateScale(double time) override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
	double sprout_mag = 3.72e-12;
	double sprout_width = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};
//used to update the old code
class LoadCurveAngioStressPolicy : public AngioStressPolicy
{
public:
	explicit LoadCurveAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~LoadCurveAngioStressPolicy() {}
	bool Init() override;
	void UpdateScale(double time) override;
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override;
	double sprout_mag = 3.72e-12;
	double sprout_width = 2;
	double sprout_range = 200;//used to calculate the falloff of stress
	double sprout_radius_multiplier = 3;//multiplied by sprout range implicitly gives the cutoff for the tips that are included
										//the default value cuts of tips stresses that are below 5% of a tip on top of a gauss point
};

class NullAngioStressPolicy : public AngioStressPolicy
{
public:
	explicit NullAngioStressPolicy(FEModel* pfem) : AngioStressPolicy(pfem) {}
	virtual ~NullAngioStressPolicy() {}
	bool Init() override { return true; }
	void UpdateScale(double time) override {}
	void AngioStress(AngioElement* angio_element, FEAngio* pangio, FEMesh* mesh) override {}
	double scale = 0.0;
};
