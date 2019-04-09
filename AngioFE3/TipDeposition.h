#pragma once
#include <FECore/BC.h>


//! (Incomplete) Prescribed dof to allow tips to deposit chemicals within a multiphasic material
class TipDepositionBC: public FEPrescribedDOF
{
public:
	//! constructor
	TipDepositionBC(FEModel* pfem): FEPrescribedDOF(pfem){}
	//! idealy changes how the concentrations are updated to be based on the location of tips
	void Update() override;
};

