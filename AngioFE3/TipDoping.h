#pragma once
#include <FECore/FEMaterial.h>
#include "Tip.h"
#include "rbf_norms.h"
#include "IntptSelector.h"
#include "ChemicalDepositInfo.h"

class TipDoping :public FEMaterial
{
public:
	explicit TipDoping(FEModel* pfem);
	void DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh);
	FEPropertyT<rbf_norm> rnorm;
	FEPropertyT<IntptSelector> selector;
	FEPropertyT<ChemicalDepositInfo> chemical_deposit_info;
	double chemical_rate = 1.0;
};

class TipDopingManager : public FEMaterial
{
public:
	explicit TipDopingManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&tip_effects, "tip_effect", false); }
	void DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh);
	FEVecPropertyT<TipDoping> tip_effects;
};
