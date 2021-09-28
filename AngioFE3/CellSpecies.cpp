#include "CellSpecies.h"
#include "FEAngio.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"
#include <iostream>

BEGIN_FECORE_CLASS(CellSpecies, FEMaterial)
ADD_PARAMETER(species_ID, "species_ID");
ADD_PARAMETER(prod_rate, "prod_rate");
END_FECORE_CLASS();