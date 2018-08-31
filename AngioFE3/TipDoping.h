#pragma once
#include <FECore/FEMaterial.h>
#include "Tip.h"
#include "rbf_norms.h"
#include "IntptSelector.h"

class TipDoping :public FEMaterial
{
public:
	explicit TipDoping(FEModel* pfem) : FEMaterial(pfem) {}
	void DopeAtTip(Tip * tip, double dt, int solute_id, double solute_amount, FEAngio* feangio, FEMesh* mesh);
	FEPropertyT<rbf_norm> norm;
	FEPropertyT<IntptSelector> selector;
};

class TipDopingManager : public FEMaterial
{
public:
	explicit TipDopingManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&tip_effects, "tip_effect", false); }
	void DopeAtTip(Tip * tip, double dt, int solute_id, double solute_amount, FEAngio* feangio, FEMesh* mesh);
	FEVecPropertyT<TipDoping> tip_effects;
};
