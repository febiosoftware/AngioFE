#pragma once
#include <FECore/BC.h>
#include "ChemicalDepositInfo.h"


class TipDopingBC: public FEPrescribedDOF
{
public:
	TipDopingBC(FEModel* pfem): FEPrescribedDOF(pfem){}
	void Update() override;
};

