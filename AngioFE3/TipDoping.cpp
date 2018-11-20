#include "TipDoping.h"
#include "FEAngio.h"

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

void TipDopingBC::Update()
{
	FEPrescribedDOF::Update();
}

bool TipDopingManager::Init()
{
	return FEMaterial::Init();
}

bool TipDoping::Init()
{
	return FEMaterial::Init();
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

void TipDoping::SetupPDOF(FEAngio * feangio, FEAngioMaterial* angio_material)
{
	//get the collection of nodes to modify
	std::vector<AngioElement*> & angio_elements = feangio->elements_by_material[angio_material];
	FEMesh * mesh = feangio->GetMesh();

	for(size_t i=0; i < angio_elements.size();i++)
	{
		for(int j=0; j < angio_elements[i]->_elem->Nodes(); j++)
		{
			nodes_set.insert(angio_elements[i]->_elem->m_lnode[j]);
		}
	}
	
	for(auto iter= nodes_set.begin(); iter != nodes_set.end();++iter)
	{
		node_vec.push_back(*iter);
	}

	fenode_set = new FENodeSet();
	fenode_set->add(node_vec);

	fe_P_dof = new FEPrescribedDOF(GetFEModel());
	fe_P_dof->AddNodes(*fenode_set, 0);
	fe_P_dof->SetRelativeFlag(true);

	//get the dof that will be modified
	int dof = GetFEModel()->GetDOFS().GetDOF(dof_name.c_str());
	fe_P_dof->SetDOF(dof);
	scale_values.resize(node_vec.size(), 0.0);
}

void TipDoping::Update(FEAngio* pangio, FEAngioMaterial* angio_material)
{
	//zero before accumulating values
	scale_values.assign(scale_values.size(), 0.0);


	//set the scale values
	for(int i=0; i < node_vec.size();i++)
	{
		fe_P_dof->SetNodeScale(node_vec[i], scale_values[i]);
	}
}