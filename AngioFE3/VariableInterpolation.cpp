#include "VariableInterpolation.h"
#include "FECore/FEElement.h"
#include "FECore/FEMesh.h"

double PerElementVI::Interpolate(FESolidElement *se, std::vector<double> & values_at_gauss_points, vec3d local_pos, FEMesh* mesh)
{
	double ao[FEElement::MAX_NODES];
	double H[FEElement::MAX_NODES];
	assert(values_at_gauss_points.size() == se->GaussPoints());
	//this hack might will only work with vectors
	se->project_to_nodes(&values_at_gauss_points[0], ao);

	
	se->shape_fnc(H, local_pos.x, local_pos.y, local_pos.z);
	double val = 0.0;
	for (int j = 0; j < se->Nodes(); j++)
	{
		val += ao[j]* H[j];
	}
	return val;
}