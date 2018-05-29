#include "TimestepUpdate.h"
#include "FEAngioMaterial.h"

void RandomFiberManager::FETimeStepUpdate()
{
	
}

void RandomFiberManager::AngioTimeStepUpdate(AngioElement* angio_element)
{
	
}

void RandomFiberManager::Initialization(AngioElement* angio_element)
{
	FESolidElement *se = angio_element->_elem;
	assert(se);
	for (int k = 0; k < se->GaussPoints(); k++)
	{
		FEMaterialPoint * mp = se->GetMaterialPoint(k);
		FEAngioMaterialPoint * angiopt = mp->ExtractData<FEAngioMaterialPoint>();
		FEElasticMaterialPoint *  emp = mp->ExtractData<FEElasticMaterialPoint>();

		FEElasticMaterialPoint *  emp_matrix = angiopt->matPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint *  emp_vessel = angiopt->vessPt->ExtractData<FEElasticMaterialPoint>();

		mat3d rm = angio_element->_angio_mat->m_pangio->unifromRandomRotationMatrix(angio_element->_rengine);
		emp->m_Q = emp->m_Q * rm;
		emp_matrix->m_Q = emp->m_Q;
		emp_vessel->m_Q = emp->m_Q;
	}
	angio_element->_angio_mat->m_pangio->uniformRandomDirection(angio_element->_rengine);
}

bool RandomFiberManager::Init()
{
	return true;
}