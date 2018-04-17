#pragma once
#include <FECore/FEElement.h>
#include <random>

class FEAngioMaterial;

class AngioElement
{
public:
	AngioElement() {};
	AngioElement(FEElement * elem, FEAngioMaterial * angio_mat, FEMaterial * mat) : _elem(elem), _angio_mat(angio_mat), _mat(mat) {};
	//  
	FEElement * _elem=nullptr;
	FEAngioMaterial * _angio_mat=nullptr;
	FEMaterial * _mat= nullptr;
	std::mt19937_64 _rengine;
};
