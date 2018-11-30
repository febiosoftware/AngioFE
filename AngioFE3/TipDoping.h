#pragma once
#include <FECore/BC.h>


//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class TipDopingBC: public FEPrescribedDOF
{
public:
	//! constructor
	TipDopingBC(FEModel* pfem): FEPrescribedDOF(pfem){}
	//! idealy changes how the concentrations are updated to be based on the location of tips
	void Update() override;
};

