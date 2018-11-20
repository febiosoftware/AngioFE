#pragma once
#include <FECore/BC.h>


class TipDopingBC: public FEPrescribedDOF
{
public:
	//! constructor
	TipDopingBC(FEModel* pfem): FEPrescribedDOF(pfem){}
	//! idealy changes how the concentrations are updated to be based on the location of tips
	void Update() override;
};

