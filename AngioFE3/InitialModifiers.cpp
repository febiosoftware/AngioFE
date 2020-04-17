#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"
#include <iostream>
//#include <FEBioMech/FEEllipsoidalFiberDistribution.h>
//#include <FEBioMech/FEFiberDensityDistribution.h>

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
		emp->m_Q = emp->m_Q * feangio->unifromRandomRotationMatrix(angio_element->_rengine);
		//then propogate the rotation to the angio submaterials
		FEElasticMaterialPoint * mat_emp = angio_pt->matPt->ExtractData<FEElasticMaterialPoint>();
		mat_emp->m_Q = emp->m_Q;
		// create a vessel material point and copy the elastic material point to it
		FEElasticMaterialPoint * ves_emp = angio_pt->vessPt->ExtractData<FEElasticMaterialPoint>();
		ves_emp->m_Q = emp->m_Q;
	}
}

void EFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	// norm the initial axes
	initial_axes_a.unit();
	initial_axes_b.unit();
	initial_axes_c.unit();
	
	// scale initial axes by values	
	angio_element->initial_angioSPA.setCol(0, initial_axes_a*initial_spa.x);
	angio_element->initial_angioSPA.setCol(1, initial_axes_b*initial_spa.y);
	angio_element->initial_angioSPA.setCol(2, initial_axes_c*initial_spa.z);
	
	// store the deformed SPA
	angio_element->angioSPA = angio_element->initial_angioSPA;
	angio_element->UpdateAngioFractionalAnisotropy();
}

BEGIN_PARAMETER_LIST(EFDFiberInitializer, InitialModifier)
ADD_PARAMETER(initial_axes_a, FE_PARAM_VEC3D, "initial_axes_a");
ADD_PARAMETER(initial_axes_b, FE_PARAM_VEC3D, "initial_axes_b");
ADD_PARAMETER(initial_axes_c, FE_PARAM_VEC3D, "initial_axes_c");
ADD_PARAMETER(initial_spa, FE_PARAM_VEC3D, "initial_spa");
END_PARAMETER_LIST();

void DensityInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++)
	{
		FEMaterialPoint *mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint *angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);

		angio_pt->ref_ecm_density = initial_density;
	}
}

BEGIN_PARAMETER_LIST(DensityInitializer, InitialModifier)
ADD_PARAMETER(initial_density, FE_PARAM_DOUBLE, "initial_density");
END_PARAMETER_LIST();