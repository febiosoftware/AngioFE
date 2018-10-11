#include "AngioElement.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "Tip.h"

double AngioElement::GetLengthAtTime(FEMesh* mesh, double time) const
{
	double len = 0.0;
	for(int i=0; i < grown_segments.size(); i++)
	{
		len += grown_segments[i]->LengthAtTime(mesh, time);
	}
	return len;
}