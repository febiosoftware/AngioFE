#include "TipDoping.h"

TipDoping::TipDoping(FEModel* pfem): FEMaterial(pfem)
{
	AddProperty(&norm, "norm");
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
	std::vector<FEMaterialPoint*> int_pts;
	selector->SelectIntPtr(int_pts, tip->angio_element->_elem, tip->GetLocalPosition(), feangio, mesh);
	std::vector<double> weights;
	assert(int_pts.size());
	for(int i=0; i < int_pts.size();i++)
	{
		
	}
}