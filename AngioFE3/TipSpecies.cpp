#include "TipSpecies.h"
#include "FEAngio.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"
#include <iostream>

BEGIN_FECORE_CLASS(TipSpecies, FEMaterial)
ADD_PARAMETER(SBM_ID, "SBM_ID");
ADD_PARAMETER(prod_rate, "prod_rate");
END_FECORE_CLASS();