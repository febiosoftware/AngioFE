#include "ChemicalDepositInfo.h"

void SoluteDepositInfo::DepositChemicals(class FESolutesMaterialPoint * smp, double amount)
{	
	int bc = GetFEModel()->GetDOFS().GetDOF("c1");
	//modify the effective concentration, actual concentration is constant wrt time/deformation
	assert(solute_id < smp->m_c.size());
	smp->m_c[solute_id] += amount;
}
