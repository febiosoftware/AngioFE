#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"
#include <iostream>
#include "angio3d.h"
#include "FECore/FEDomainMap.h"

#pragma region FECoreClassDefs
BEGIN_FECORE_CLASS(InitialModifierManager,FEMaterialProperty)
	ADD_PROPERTY(initial_modifiers, "initial_modifier", FEProperty::Optional);
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(DiscreteFiberEFDRandomizer, InitialModifier)
	ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(EFDFiberInitializer, InitialModifier)
	ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(RotEFDFiberInitializer, InitialModifier)
	ADD_PARAMETER(m_SPD, "spd");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(DensityInitializer, InitialModifier)
	ADD_PARAMETER(initial_density, "initial_density");
END_FECORE_CLASS()

BEGIN_FECORE_CLASS(RepulseInitializer, InitialModifier)
	ADD_PARAMETER(initial_repulse_value, "initial_repulse_value");
END_FECORE_CLASS()
#pragma endregion FECoreClassDefs

void InitialModifierManager::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i< initial_modifiers.size(); i++)
		initial_modifiers[i]->ApplyModifier(angio_element, mesh, feangio);
}

//! Rename this to AlignedFiberRandomizer or DiscreteFiberRandomizer
//! Take a given angio element and give it a randomized discrete fiber direction
void FiberRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	auto se = angio_element->_elem;
	// for each integration point in the element
	for (int i = 0; i < se->GaussPoints(); i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = se->GetMaterialPoint(i);
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
	auto se = angio_element->_elem;
	// for each integration point in the element
	for (int i = 0; i < se->GaussPoints(); i++)
	{
		// get the material point of the integration point
		FEMaterialPoint * mp = se->GetMaterialPoint(i);
		// get the angio material point
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->initial_angioSPD = m_SPD(*mp);
		FEEllipticalDistribution E(this->GetFEModel());
		E.spd = angio_pt->initial_angioSPD;
		E.Init();
		angio_pt->angio_fiber_dir = E.NextVec(angio_element->_rengine);
	}
}

void EFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	FESolidElement* se = angio_element->_elem;
	for (auto i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->initial_angioSPD = m_SPD(*mp);
		angio_pt->angioSPD = m_SPD(*mp);
		angio_pt->UpdateAngioFractionalAnisotropy();
	}
}

void RotEFDFiberInitializer::ApplyModifier(AngioElement* angio_element, FEMesh* mesh, FEAngio* feangio)
{
	FESolidElement* se = angio_element->_elem;
	auto zto2pi = std::uniform_real_distribution<double>(0.0, 2.0 * PI);
	double alpha = zto2pi(angio_element->_rengine);
	double beta = 0.0;
	double gamma = zto2pi(angio_element->_rengine);

	mat3d RandMat = feangio->rotationMatrix(alpha,beta,gamma);
	
	for (auto i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint* mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint* angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		
		mat3d P = RandMat*m_SPD(*mp)*RandMat.transpose();
		angio_pt->initial_angioSPD = mat3ds(P(0,0), P(1,1), P(2,2), P(0,1), P(1,2), P(0,2));
		angio_pt->angioSPD = m_SPD(*mp);
		angio_pt->UpdateAngioFractionalAnisotropy();
	}
}

void DensityInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	auto se = angio_element->_elem;
	for (int i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->ref_ecm_density = initial_density(*mp);
	}
}

void RepulseInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	auto se = angio_element->_elem;
	for (int i = 0; i < se->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = se->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		angio_pt->repulse_value = initial_repulse_value;
	}
}