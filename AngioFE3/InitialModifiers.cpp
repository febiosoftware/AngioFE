#include "InitialModifiers.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"

void InitialModifierManager::ApplyModifier(AngioElement * angio_element, FEMesh * mesh, FEAngio* feangio)
{
	for(int i=0; i< initial_modifiers.size();i++)
	{
		initial_modifiers[i]->ApplyModifier(angio_element, mesh, feangio);
	}
}



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
