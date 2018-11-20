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

//! initializes the fibers to random direcitons
class RandomFiberManager : public TimeStepUpdate
{
public:
	//! constructor
	void FETimeStepUpdate()override;
	//! Update the the fibers in an element to a given time
	void AngioTimeStepUpdate(AngioElement* angio_element)override;
	//! Initializes the fibers within an element
	void Initialization(AngioElement* angio_element)override;
	//! performs initialization
	bool Init() override;
private:
};