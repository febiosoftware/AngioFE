#include "AngioElement.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "Tip.h"
#include <iostream>
#include <math.h>
#include <numeric>
#include <algorithm>

#ifndef PI
#define PI 3.14159265358979
#endif

double AngioElement::GetLengthAtTime(FEMesh* mesh, double time) const
{
	double len = 0.0;
	for(int i=0; i < grown_segments.size(); i++)
	{
		len += grown_segments[i]->LengthAtTime(mesh, time);
	}
	return len;
}