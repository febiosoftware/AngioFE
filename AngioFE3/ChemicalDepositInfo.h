#pragma once
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FEBioMix/FEMultiphasic.h>

//provides the descriptions for how much chemical property is released by a tip
//


class ChemicalDepositInfo : public FEMaterial
{
public:
	explicit ChemicalDepositInfo(FEModel* pfem) : FEMaterial(pfem) {}
	virtual void DepositChemicals(class FESolutesMaterialPoint * smp, double amount)=0;
};

class SoluteDepositInfo : public ChemicalDepositInfo
{
public:
	explicit SoluteDepositInfo(FEModel* pfem) : ChemicalDepositInfo(pfem) {}
	void DepositChemicals(class FESolutesMaterialPoint * smp, double amount) override;
	int solute_id=0;
};
