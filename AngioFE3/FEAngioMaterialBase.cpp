#include "StdAfx.h"
#include "FEAngioMaterialBase.h"
#include "KDTree/kdtree.h"

#include <FECore/FEElement.h>
#include "FEAngioMaterial.h"


std::vector<double> units3d(3, 1.0);

//-----------------------------------------------------------------------------
vec3d FEAngioMaterialBase::CurrentPosition(FESolidElement * pe, double r, double s, double t) const
{
	double arr[FEElement::MAX_NODES];
	FEMesh * mesh = m_pangio->GetMesh();
	vec3d rc(0, 0, 0);

	assert(pe);
	pe->shape_fnc(arr, r, s, t);
	for (int j = 0; j < pe->Nodes(); j++)
	{
		rc += mesh->Node(pe->m_node[j]).m_rt* arr[j];
	}
	return rc;
}


//begin implementation of FEAngioBase
FEAngioMaterialBase::FEAngioMaterialBase()
{
	scale = 1.0;

	m_pangio = nullptr;
}

void FEAngioMaterialBase::UpdateFiberManager()
{
	fiber_manager->Update();
}

void FEAngioMaterialBase::AdjustMeshStiffness(FEMaterial* mat)
{
	FEMesh & mesh = m_pangio->m_fem->GetMesh();
	if (m_cultureParams.m_composite_material == 0)													// If a composite consitutive model isn't being used, exit
		return;

	int elem_num = 0;													// Element number
	vec3d vess_vect;													// Vessel vector
	std::vector<int> matls;
	matls.emplace_back(this->GetID_ang());
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1 / static_cast<double>(Nsub);									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

																		//Zero the element items needed
																		//break even on core in field model
	m_pangio->ForEachElement([&](FESolidElement & se, FESolidDomain & d)
	{
		int elemnum = se.GetID();

	}, matls);

	return;
}

void FEAngioMaterialBase::UpdateSproutStressScaling()
{
	//TODO: make these user parameters and get better names for these
	double y0 = -0.004; double x0 = 2.0; double b = 0.5436; double a = 1.0081;

	auto time = m_pangio->CurrentSimTime();

	scale = y0 + a / (1 + exp(-(time.t - x0) / b));
}


void FEAngioMaterialBase::Update()
{
	fiber_manager->Update();
}

bool FEAngioMaterialBase::Overwrite() const
{
	return ecm_initializer->overwrite();
}

double FEAngioMaterialBase::GetAnisotropy() const
{
	return m_cultureParams.GetWeightInterpolation(1.0);
}