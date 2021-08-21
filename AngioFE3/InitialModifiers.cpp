#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"
#include <iostream>
#include "angio3d.h"
#include "FECore/FEDomainMap.h"

void InitialModifierManager::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for(int i=0; i< initial_modifiers.size();i++)
	{
		initial_modifiers[i]->ApplyModifier(angio_element, mesh, feangio);
	}
}

// Rename this to AlignedFiberRandomizer or DiscreteFiberRandomizer
// take a given angio element and give it a randomized discrete fiber direction
void FiberRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	// for each integration point in the element
	for(int i=0; i < angio_element->_elem->GaussPoints();i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		// get the elastic material point
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		// multiply the elastic material point's global orientation by a random matrix
		mat3d RandMat = feangio->unifromRandomRotationMatrix(angio_element->_rengine);
		// get local domain index of element
		angio_pt->angio_fiber_dir = RandMat*angio_pt->angio_fiber_dir;
	}
}

// take a given angio element and give it a randomized discrete fiber direction based on user input spd
void DiscreteFiberEFDRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	// for each integration point in the element
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->initial_angioSPD = m_SPD;
		FEEllipticalDistribution E(this->GetFEModel());
		E.spd = angio_pt->initial_angioSPD;
		E.Init();
		angio_pt->angio_fiber_dir = E.NextVec(angio_element->_rengine);
	}
}

BEGIN_FECORE_CLASS(DiscreteFiberEFDRandomizer, InitialModifier)
ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS();

void EFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	FESolidElement* se = angio_element->_elem;
	for (auto i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->initial_angioSPD = m_SPD;
		angio_pt->angioSPD = m_SPD;
		angio_pt->UpdateAngioFractionalAnisotropy();
	}
}

BEGIN_FECORE_CLASS(EFDFiberInitializer, InitialModifier)
ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS();

void DensityInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		angio_pt->ref_ecm_density = initial_density;
	}
}

BEGIN_FECORE_CLASS(DensityInitializer, InitialModifier)
ADD_PARAMETER(initial_density, "initial_density");
END_FECORE_CLASS();

void RepulseInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		angio_pt->repulse_value = initial_repulse_value;
	}
}

BEGIN_FECORE_CLASS(RepulseInitializer, InitialModifier)
ADD_PARAMETER(initial_repulse_value, "initial_repulse_value");
END_FECORE_CLASS();