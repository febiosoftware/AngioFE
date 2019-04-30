#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"
//#include <FEBioMech/FEEllipsoidalFiberDistribution.h>
//#include <FEBioMech/FEFiberDensityDistribution.h>

void InitialModifierManager::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for(int i=0; i< initial_modifiers.size();i++)
	{
		initial_modifiers[i]->ApplyModifier(angio_element, mesh, feangio);
	}
}
/*
void IsotropicEFDFiberInitializer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for (int i = 0; i < angio_element->_elem->GaussPoints(); i++) {
		//FEParam* mf = pfem->FindParameter("m_spa");
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		
		//FEEllipsodialFiberDensityDistribution * angio_FDD = mp->ExtractData<FEEllipsodialFiberDensityDistribution>();
	}
}
*/
// Rename this to AlignedFiberRandomizer or DiscreteFiberRandomizer
void FiberRandomizer::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for(int i=0; i < angio_element->_elem->GaussPoints();i++)
	{
		FEMaterialPoint * mp = angio_element->_elem->GetMaterialPoint(i);
		FEAngioMaterialPoint * angio_pt = FEAngioMaterialPoint::FindAngioMaterialPoint(mp);
		FEElasticMaterialPoint * emp = mp->ExtractData<FEElasticMaterialPoint>();
		emp->m_Q = emp->m_Q * feangio->unifromRandomRotationMatrix(angio_element->_rengine);
		//then propogate the rotation to the angio submaterials
		FEElasticMaterialPoint * mat_emp = angio_pt->matPt->ExtractData<FEElasticMaterialPoint>();
		mat_emp->m_Q = emp->m_Q;
		FEElasticMaterialPoint * ves_emp = angio_pt->vessPt->ExtractData<FEElasticMaterialPoint>();
		ves_emp->m_Q = emp->m_Q;
	}
}

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