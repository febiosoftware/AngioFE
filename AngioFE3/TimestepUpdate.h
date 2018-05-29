#pragma once
#include <FECore/FEMaterial.h>
#include "AngioElement.h"

class TimeStepUpdate : public FEMaterial
{
public:
	virtual void FETimeStepUpdate(){}
	virtual void AngioTimeStepUpdate(AngioElement* angio_element){}
	virtual void Initialization(AngioElement* angio_element) { }
private:

};

class RandomFiberManager : public TimeStepUpdate
{
public:
	void FETimeStepUpdate()override;
	void AngioTimeStepUpdate(AngioElement* angio_element)override;
	void Initialization(AngioElement* angio_element)override;
	bool Init() override;
private:
};