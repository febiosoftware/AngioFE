#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/BC.h>
#include "Tip.h"
#include "rbf_norms.h"
#include "IntptSelector.h"
#include "ChemicalDepositInfo.h"

class TipDoping :public FEMaterial
{
public:
	explicit TipDoping(FEModel* pfem);
	bool Init()override;
	void DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh);
	void SetupPDOF(FEAngio * feangio, FEAngioMaterial* angio_material);
	void Update(FEAngio* pangio, FEAngioMaterial* angio_material);
	FEPropertyT<rbf_norm> rnorm;
	FEPropertyT<IntptSelector> selector;
	FEPropertyT<ChemicalDepositInfo> chemical_deposit_info;
	double chemical_rate = 1.0;
private:
	FEPrescribedDOF * fe_P_dof = nullptr;
	FENodeSet * fenode_set = nullptr;
	std::unordered_set<int> nodes_set;
	std::vector<int> node_vec;
	std::vector<double> scale_values;
	std::string dof_name = "c1";
};

class TipDopingBC: public FEPrescribedDOF
{
public:
	TipDopingBC(FEModel* pfem): FEPrescribedDOF(pfem){}
	void Update() override;
};

class TipDopingManager : public FEMaterial
{
public:
	explicit TipDopingManager(FEModel* pfem) : FEMaterial(pfem) { AddProperty(&tip_effects, "tip_effect", false); }
	bool Init()override;
	void DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh);
	FEVecPropertyT<TipDoping> tip_effects;
};
