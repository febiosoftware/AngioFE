#include "TipDoping.h"

TipDoping::TipDoping(FEModel* pfem): FEMaterial(pfem)
{
	AddProperty(&rnorm, "rnorm");
	AddProperty(&selector, "selector");
	AddProperty(&chemical_deposit_info, "chemical_deposit_info");
}

void TipDopingManager::DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh)
{
	for(int i=0; i < tip_effects.size();i++)
	{
		tip_effects[i]->DopeAtTip(tip, dt, feangio, mesh);
	}
}

void TipDoping::DopeAtTip(Tip * tip, double dt, FEAngio* feangio, FEMesh* mesh)
{
	//get the relevant intergation points
	std::vector<FEMaterialPoint*> int_pts;
	selector->SelectIntPtr(int_pts, tip->angio_element->_elem, tip->GetLocalPosition(), feangio, mesh);

	//relative weights of the relevant integration points
	std::vector<double> weights;
	assert(int_pts.size());
	double total_weights = 0.0;
	for(int i=0; i < int_pts.size();i++)
	{
		FEElasticMaterialPoint* emp = int_pts[i]->ExtractData<FEElasticMaterialPoint>();
		assert(emp);
		double weight = rnorm->rbfnorm(emp->m_rt, tip->GetPosition(mesh));
		weights.push_back(weight);
		total_weights += weight;
	}

	double amount = dt * chemical_rate;
	assert(total_weights > 0.0);

	for (int i = 0; i < int_pts.size(); i++)
	{
		FESolutesMaterialPoint * smp = int_pts[i]->ExtractData<FESolutesMaterialPoint>();
		assert(smp);
		chemical_deposit_info->DepositChemicals(smp, (amount* weights[i])/total_weights);
	}
}