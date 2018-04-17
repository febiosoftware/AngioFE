#include "GrowDirectionModifier.h"
#include "FECore/ElementDataRecord.h"
#include "FEBioPlot/FEBioPlotFile2.h"
#include "FEAngio.h"
#include "angio3d.h"

/*
GrowDirectionModifier::GrowDirectionModifier(FEModel * model) : FEMaterial(model)
{

}

void GrowDirectionModifier::SetCulture(Culture * cp)
{
	culture = cp;
}

GrowDirectionModifiers::GrowDirectionModifiers(FEModel* model) : FEMaterial(model)
{
	AddProperty(&grow_direction_modifiers, "gdm");
}


vec3d GrowDirectionModifiers::ApplyModifiers(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterialBase* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		previous_dir = grow_direction_modifiers[i]->GrowModifyGrowDirection(previous_dir, tip, mat, branch, start_time, grow_time, seg_length);
	}
	return previous_dir;
}

void GrowDirectionModifiers::SetCulture(Culture * c)
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		grow_direction_modifiers[i]->SetCulture(c);
	}
}

void GrowDirectionModifiers::Update()
{
	for (int i = 0; i < grow_direction_modifiers.size(); i++)
	{
		grow_direction_modifiers[i]->Update();
	}
}



void GDMArchive::reset()
{
	fpdata.clear();
}

void GDMArchive::WriteData(int nid, std::vector<float>& data)
{
	if(fpdata.size() <= nid)
	{
		vector<float> temp;
		fpdata.push_back(temp);
	}
		
	fpdata[nid-1].resize(data.size());
	std::copy(data.begin(), data.end(), fpdata[nid-1].begin());
	//consider preallocating unrolled data
	if(gradient_defined)
	{
		//only reset at once per step
		//may be a hack
		if(nid == 1)
		{
			node_hit_count.clear();
			projected_nodal_data.clear();
			node_hit_count.resize(angio->GetMesh()->Nodes(), 0);
			projected_nodal_data.resize(angio->GetMesh()->Nodes(), 0);
		}
		
		std::vector<int> dl;
		angio->GetMesh()->DomainListFromMaterial(angio->m_pmat_ids, dl);

			FESolidDomain & d = reinterpret_cast<FESolidDomain&>(angio->GetMesh()->Domain(dl[nid-1]));
			for (int j = 0; j < d.Elements(); j++)
			{
				FESolidElement & e = reinterpret_cast<FESolidElement&>(d.ElementRef(j));
				for (int k = 0; k< e.Nodes(); k++)
				{
					node_hit_count[e.m_node[k]]++;
					projected_nodal_data[e.m_node[k]] += fpdata[nid-1][j];
				}
			}
		}
}
mat3dd GDMArchive::GetDataMat3dd(int domain, int element_index)
{
	const int index = 3 * element_index;
	return mat3dd(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2]);
}
mat3ds GDMArchive::GetDataMat3ds(int domain, int element_index)
{
	const int index = 6 * element_index;
	return mat3ds(fpdata[domain][index], fpdata[domain][index+1], fpdata[domain][index+2],
		fpdata[domain][index+3], fpdata[domain][index+4], fpdata[domain][index+5]);

}
mat3d  GDMArchive::GetDataMat3d(int domain, int element_index)
{
	const int index = 9 * element_index;
	return mat3d(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2],
		fpdata[domain][index + 3], fpdata[domain][index + 4], fpdata[domain][index + 5],
		fpdata[domain][index + 6], fpdata[domain][index + 7], fpdata[domain][index + 8]);
}
float  GDMArchive::GetDataFloat(int domain, int element_index)
{
	const int index = element_index;
	return fpdata[domain][index];
}
vec3d  GDMArchive::GetDataVec3d(int domain, int element_index)
{
	const int index = 3 * element_index;
	return vec3d(fpdata[domain][index], fpdata[domain][index + 1], fpdata[domain][index + 2]);
}
//gradient version of element functions
vec3d  GDMArchive::GetDataGradientFloat(int domain, int element_index, Segment::TIP& tip,int size,int offset)
{
	FESolidElement * tse = &tip.pt.ndomain->Element(element_index);
	std::vector<double> v2g;
	v2g.resize(tse->Nodes());
	for(int i =0; i < tse->Nodes();i++)
	{
		v2g[i] = projected_nodal_data[tse->m_node[i]] / node_hit_count[tse->m_node[i]];
	}
	vec3d d1 = FEAngio::gradient(tse,v2g, tip.pt.q, size, offset);
	return d1;
}


//node versions of fucntions
mat3dd GDMArchive::GetDataMat3dd(int domain, int element_index, Segment::TIP& tip)
{
	mat3dd r(0,0,0);
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = 3 * se->m_node[j];
			r +=  mat3dd(fpdata[domain][ni], fpdata[domain][ni+1], fpdata[domain][ni+2]) * arr[j];
		}
	}
	return r;
}
mat3ds GDMArchive::GetDataMat3ds(int domain, int element_index, Segment::TIP& tip)
{
	mat3ds r(0, 0, 0, 0, 0, 0);
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = 6 * se->m_node[j];
			r += mat3ds(fpdata[domain][ni], fpdata[domain][ni + 1], fpdata[domain][ni + 2],
				fpdata[domain][ni + 3], fpdata[domain][ni + 4], fpdata[domain][ni + 5]) * arr[j];
		}
	}
	return r;
}
mat3d  GDMArchive::GetDataMat3d(int domain, int element_index, Segment::TIP& tip)
{
	mat3d r(0, 0, 0, 0, 0, 0, 0, 0, 0);
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = 9 * se->m_node[j];
			r += mat3d(fpdata[domain][ni], fpdata[domain][ni + 1], fpdata[domain][ni + 2],
				fpdata[domain][ni + 3], fpdata[domain][ni + 4], fpdata[domain][ni + 5],
				fpdata[domain][ni + 6], fpdata[domain][ni + 7], fpdata[domain][ni + 8]) * arr[j];
		}
	}
	return r;
}
float  GDMArchive::GetDataFloat(int domain, int element_index, Segment::TIP& tip)
{
	float r = 0.0;
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = se->m_node[j];
			r += fpdata[domain][ni] * arr[j];
		}
	}
	return r;
}
vec3d  GDMArchive::GetDataVec3d(int domain, int element_index, Segment::TIP& tip)
{
	vec3d r = vec3d(0,0,0);
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = 3 * se->m_node[j];
			r += vec3d(fpdata[domain][ni], fpdata[domain][ni+1], fpdata[domain][ni+2]) * arr[j];
		}
	}
	return r;
}

bool Plot2GGP::Init()
{
	FEModel * model = GetFEModel();
	FEBioModel * bm = dynamic_cast<FEBioModel*>(model);
	assert(bm);
	FEBioPlotFile2 * pf2 = dynamic_cast<FEBioPlotFile2 *>(bm->GetPlotFile());
	if (!pf2)
	{
		printf("plot file not in plot file 2 format\n");
		return false;
	}
	const list<FEBioPlotFile2::DICTIONARY_ITEM>& domain_var_list = pf2->GetDictionary().DomainVariableList();


	for (auto i = domain_var_list.begin(); i != domain_var_list.end(); ++i)
	{
		if (!strcmp(i->m_szname, field_name))
		{
			record_index = i;
			return GGP::Init();
		}
	}
	printf("\nfield %s not found\n", field_name);
	//data record name not found
	return false;
}

mat3d Plot2GGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	switch (record_index->m_nfmt)
	{
	case Storage_Fmt::FMT_ITEM:
		{
			switch(record_index->m_ntype)
			{
			case Var_Type::PLT_FLOAT:
				in(0,0) = archive.GetDataFloat(mat->domains[0], tip.pt.elemindex);
				break;
			case Var_Type::PLT_MAT3F:
				{
					//make sure this order is correct
					mat3d temp = archive.GetDataMat3d(mat->domains[0], tip.pt.elemindex);
					in = in * temp;
				}
				break;
			case Var_Type::PLT_MAT3FD:
				{
					mat3dd temp = archive.GetDataMat3dd(mat->domains[0], tip.pt.elemindex);
					in = in *temp;
				}
				break;
			case Var_Type::PLT_MAT3FS:
				{
					mat3ds temp = archive.GetDataMat3ds(mat->domains[0], tip.pt.elemindex);
					in = in *temp;
				}
				break;
			case Var_Type::PLT_VEC3F:
				{
					vec3d col = archive.GetDataVec3d(mat->domains[0], tip.pt.elemindex);
					in(0, 0) = col.x;
					in(1, 1) = col.y;
					in(2, 2) = col.z;
				}
				break;
			default:
				assert(false);
			}
		}
		break;
	case Storage_Fmt::FMT_NODE:
		{
			switch (record_index->m_ntype)
			{
			case Var_Type::PLT_FLOAT:
				in(0, 0) = archive.GetDataFloat(mat->domains[0], tip.pt.elemindex, tip);
				break;
			case Var_Type::PLT_MAT3F:
			{
				//make sure this order is correct
				mat3d temp = archive.GetDataMat3d(mat->domains[0], tip.pt.elemindex, tip);
				in = in * temp;
			}
			break;
			case Var_Type::PLT_MAT3FD:
			{
				mat3dd temp = archive.GetDataMat3dd(mat->domains[0], tip.pt.elemindex, tip);
				in = in *temp;
			}
			break;
			case Var_Type::PLT_MAT3FS:
			{
				mat3ds temp = archive.GetDataMat3ds(mat->domains[0], tip.pt.elemindex, tip);
				in = in *temp;
			}
			break;
			case Var_Type::PLT_VEC3F:
			{
				vec3d col = archive.GetDataVec3d(mat->domains[0], tip.pt.elemindex, tip);
				in(0, 0) = col.x;
				in(1, 0) = col.y;
				in(2, 0) = col.z;
			}
			break;
			default:
				assert(false);
			}
		}
		break;
	case Storage_Fmt::FMT_REGION:
	default:
		assert(false);
	}

	return GGP::Operation(in, fin, mat, tip);
}

mat3d GradientPlot2GGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	switch (record_index->m_nfmt)
	{
	case Storage_Fmt::FMT_ITEM:
	{
		switch (record_index->m_ntype)
		{
		case Var_Type::PLT_FLOAT:
		{
			vec3d grad = archive.GetDataGradientFloat(mat->domains[0], tip.pt.elemindex, tip, size, offset);
			in[0][0] = grad.x;
			in[1][1] = grad.y;
			in[2][2] = grad.z;
		}
			break;
		default:
			assert(false);
		}
	}
	break;
	case Storage_Fmt::FMT_NODE:
	{
		switch (record_index->m_ntype)
		{
		default:
			assert(false);
		}
	}
	break;
	case Storage_Fmt::FMT_REGION:
	default:
		assert(false);
	}

	return GGP::Operation(in, fin, mat, tip);
}


BEGIN_PARAMETER_LIST(Plot2GGP, GGP)
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");
END_PARAMETER_LIST();

MatrixConverterGGP::MatrixConverterGGP(FEModel * model) : GGP(model)
{
	//initialize the convrsion matrix as an identity matrix
	for(int i=0; i < 9;i++)
	{
		for(int j=0; j < 9; j++)
		{
			m[i][j] = 0.0;
		}
	}

	for(int i=0;i < 9;i++)
	{
		m[i][i] = 1.0;
	}
}

mat3d MatrixConverterGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	double unroll[9];
	double results[9];
	int index = 0;
	for(int i =0;i<3;i++)
	{
		for(int j=0; j < 3;j++)
		{
			unroll[index] = in[i][j];
			index++;
		}
	}

	for(int i =0; i < 9;i++)
	{
		double sum = 0.0;
		for(int j=0; j< 9;j++)
		{
			sum += unroll[j] * m[i][j];
		}
		results[i] = sum;
	}


	index = 0;
	for (int i = 0; i<3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			in[i][j] = results[index];
			index++;
		}
	}
	return GGP::Operation(in, fin, mat, tip);
}

BEGIN_PARAMETER_LIST(MatrixConverterGGP, GGP)
ADD_PARAMETER(m[0][0], FE_PARAM_DOUBLE, "m11");
ADD_PARAMETER(m[0][1], FE_PARAM_DOUBLE, "m12");
ADD_PARAMETER(m[0][2], FE_PARAM_DOUBLE, "m13");
ADD_PARAMETER(m[0][3], FE_PARAM_DOUBLE, "m14");
ADD_PARAMETER(m[0][4], FE_PARAM_DOUBLE, "m15");
ADD_PARAMETER(m[0][5], FE_PARAM_DOUBLE, "m16");
ADD_PARAMETER(m[0][6], FE_PARAM_DOUBLE, "m17");
ADD_PARAMETER(m[0][7], FE_PARAM_DOUBLE, "m18");
ADD_PARAMETER(m[0][8], FE_PARAM_DOUBLE, "m19");

ADD_PARAMETER(m[1][0], FE_PARAM_DOUBLE, "m21");
ADD_PARAMETER(m[1][1], FE_PARAM_DOUBLE, "m22");
ADD_PARAMETER(m[1][2], FE_PARAM_DOUBLE, "m23");
ADD_PARAMETER(m[1][3], FE_PARAM_DOUBLE, "m24");
ADD_PARAMETER(m[1][4], FE_PARAM_DOUBLE, "m25");
ADD_PARAMETER(m[1][5], FE_PARAM_DOUBLE, "m26");
ADD_PARAMETER(m[1][6], FE_PARAM_DOUBLE, "m27");
ADD_PARAMETER(m[1][7], FE_PARAM_DOUBLE, "m28");
ADD_PARAMETER(m[1][8], FE_PARAM_DOUBLE, "m29");

ADD_PARAMETER(m[2][0], FE_PARAM_DOUBLE, "m31");
ADD_PARAMETER(m[2][1], FE_PARAM_DOUBLE, "m32");
ADD_PARAMETER(m[2][2], FE_PARAM_DOUBLE, "m33");
ADD_PARAMETER(m[2][3], FE_PARAM_DOUBLE, "m34");
ADD_PARAMETER(m[2][4], FE_PARAM_DOUBLE, "m35");
ADD_PARAMETER(m[2][5], FE_PARAM_DOUBLE, "m36");
ADD_PARAMETER(m[2][6], FE_PARAM_DOUBLE, "m37");
ADD_PARAMETER(m[2][7], FE_PARAM_DOUBLE, "m38");
ADD_PARAMETER(m[2][8], FE_PARAM_DOUBLE, "m39");

ADD_PARAMETER(m[3][0], FE_PARAM_DOUBLE, "m41");
ADD_PARAMETER(m[3][1], FE_PARAM_DOUBLE, "m42");
ADD_PARAMETER(m[3][2], FE_PARAM_DOUBLE, "m43");
ADD_PARAMETER(m[3][3], FE_PARAM_DOUBLE, "m44");
ADD_PARAMETER(m[3][4], FE_PARAM_DOUBLE, "m45");
ADD_PARAMETER(m[3][5], FE_PARAM_DOUBLE, "m46");
ADD_PARAMETER(m[3][6], FE_PARAM_DOUBLE, "m47");
ADD_PARAMETER(m[3][7], FE_PARAM_DOUBLE, "m48");
ADD_PARAMETER(m[3][8], FE_PARAM_DOUBLE, "m49");

ADD_PARAMETER(m[4][0], FE_PARAM_DOUBLE, "m51");
ADD_PARAMETER(m[4][1], FE_PARAM_DOUBLE, "m52");
ADD_PARAMETER(m[4][2], FE_PARAM_DOUBLE, "m53");
ADD_PARAMETER(m[4][3], FE_PARAM_DOUBLE, "m54");
ADD_PARAMETER(m[4][4], FE_PARAM_DOUBLE, "m55");
ADD_PARAMETER(m[4][5], FE_PARAM_DOUBLE, "m56");
ADD_PARAMETER(m[4][6], FE_PARAM_DOUBLE, "m57");
ADD_PARAMETER(m[4][7], FE_PARAM_DOUBLE, "m58");
ADD_PARAMETER(m[4][8], FE_PARAM_DOUBLE, "m59");

ADD_PARAMETER(m[5][0], FE_PARAM_DOUBLE, "m61");
ADD_PARAMETER(m[5][1], FE_PARAM_DOUBLE, "m62");
ADD_PARAMETER(m[5][2], FE_PARAM_DOUBLE, "m63");
ADD_PARAMETER(m[5][3], FE_PARAM_DOUBLE, "m64");
ADD_PARAMETER(m[5][4], FE_PARAM_DOUBLE, "m65");
ADD_PARAMETER(m[5][5], FE_PARAM_DOUBLE, "m66");
ADD_PARAMETER(m[5][6], FE_PARAM_DOUBLE, "m67");
ADD_PARAMETER(m[5][7], FE_PARAM_DOUBLE, "m68");
ADD_PARAMETER(m[5][8], FE_PARAM_DOUBLE, "m69");

ADD_PARAMETER(m[6][0], FE_PARAM_DOUBLE, "m71");
ADD_PARAMETER(m[6][1], FE_PARAM_DOUBLE, "m72");
ADD_PARAMETER(m[6][2], FE_PARAM_DOUBLE, "m73");
ADD_PARAMETER(m[6][3], FE_PARAM_DOUBLE, "m74");
ADD_PARAMETER(m[6][4], FE_PARAM_DOUBLE, "m75");
ADD_PARAMETER(m[6][5], FE_PARAM_DOUBLE, "m76");
ADD_PARAMETER(m[6][6], FE_PARAM_DOUBLE, "m77");
ADD_PARAMETER(m[6][7], FE_PARAM_DOUBLE, "m78");
ADD_PARAMETER(m[6][8], FE_PARAM_DOUBLE, "m79");

ADD_PARAMETER(m[7][0], FE_PARAM_DOUBLE, "m81");
ADD_PARAMETER(m[7][1], FE_PARAM_DOUBLE, "m82");
ADD_PARAMETER(m[7][2], FE_PARAM_DOUBLE, "m83");
ADD_PARAMETER(m[7][3], FE_PARAM_DOUBLE, "m84");
ADD_PARAMETER(m[7][4], FE_PARAM_DOUBLE, "m85");
ADD_PARAMETER(m[7][5], FE_PARAM_DOUBLE, "m86");
ADD_PARAMETER(m[7][6], FE_PARAM_DOUBLE, "m87");
ADD_PARAMETER(m[7][7], FE_PARAM_DOUBLE, "m88");
ADD_PARAMETER(m[7][8], FE_PARAM_DOUBLE, "m89");

ADD_PARAMETER(m[8][0], FE_PARAM_DOUBLE, "m91");
ADD_PARAMETER(m[8][1], FE_PARAM_DOUBLE, "m92");
ADD_PARAMETER(m[8][2], FE_PARAM_DOUBLE, "m93");
ADD_PARAMETER(m[8][3], FE_PARAM_DOUBLE, "m94");
ADD_PARAMETER(m[8][4], FE_PARAM_DOUBLE, "m95");
ADD_PARAMETER(m[8][5], FE_PARAM_DOUBLE, "m96");
ADD_PARAMETER(m[8][6], FE_PARAM_DOUBLE, "m97");
ADD_PARAMETER(m[8][7], FE_PARAM_DOUBLE, "m98");
ADD_PARAMETER(m[8][8], FE_PARAM_DOUBLE, "m99");

END_PARAMETER_LIST();


mat3d ForkedGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d temp = nest->Operation(in, fin, mat, tip);
	return GGP::Operation(temp, fin, mat, tip);
}

BEGIN_PARAMETER_LIST(MatrixMixGGP, GGP)
ADD_PARAMETER(alpha, FE_PARAM_DOUBLE, "alpha");
END_PARAMETER_LIST();

mat3d MatrixMixGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d o = other->Operation(in, fin, mat, tip);
	mat3d c = GGP::Operation(in, fin, mat, tip);
	double omalpha = 1 - alpha;
	return mat3d(
		o[0][0] * omalpha + c[0][0] * alpha, o[0][1] * omalpha + c[0][1] * alpha, o[0][2] * omalpha + c[0][2] * alpha,
		o[1][0] * omalpha + c[1][0] * alpha, o[1][1] * omalpha + c[1][1] * alpha, o[1][2] * omalpha + c[1][2] * alpha, 
		o[2][0] * omalpha + c[2][0] * alpha, o[2][1] * omalpha + c[2][1] * alpha, o[2][2] * omalpha + c[2][2] * alpha);
}

mat3d EigenValuesGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3ds temp(in[0][0], in[1][1], in[2][2], in[0][1], in[1][2], in[0][2]);
	double eig_val[3];
	vec3d eig_vec[3];
	temp.eigen(eig_val, eig_vec);

	mat3d rv(eig_val[0], 0.0,0.0,
		0.0, eig_val[1], 0.0,
		0.0, 0.0, eig_val[2]);
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d EigenVectorsGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3ds temp(in[0][0], in[1][1], in[2][2], in[0][1], in[1][2], in[0][2]);
	double eig_val[3];
	vec3d eig_vec[3];
	temp.eigen(eig_val, eig_vec);

	mat3d rv(eig_vec[0].x, eig_vec[0].y, eig_vec[0].z,
		eig_vec[1].x, eig_vec[1].y, eig_vec[1].z,
		eig_vec[2].x, eig_vec[2].y, eig_vec[2].z);
	return GGP::Operation(rv, fin, mat, tip);
}


mat3d CrossGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d v1_rez = v1->Operation(in, fin, mat, tip);
	mat3d v2_rez = v2->Operation(in, fin, mat, tip);
	vec3d a(v1_rez[0][0], v1_rez[1][1], v1_rez[2][2]);
	vec3d b(v2_rez[0][0], v2_rez[1][1], v2_rez[2][2]);
	vec3d res = a ^ b;

	mat3d rv(res.x, 0.0, 0.0,
		0.0, res.y, 0.0,
		0.0, 0.0, res.z);
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d ThresholdGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d tr_rez = threshold->Operation(in, fin, mat, tip);
	vec3d x = vec;
	x = tr_rez * x;

	if(x.x > condition)
	{
		mat3d rv = statement->Operation(in, fin, mat, tip);
		return GGP::Operation(rv, fin, mat, tip);
	}
	return GGP::Operation(in, fin, mat, tip);
}

BEGIN_PARAMETER_LIST(ThresholdGGP, GGP)
ADD_PARAMETER(vec.x, FE_PARAM_DOUBLE, "vec_x");
ADD_PARAMETER(vec.y, FE_PARAM_DOUBLE, "vec_y");
ADD_PARAMETER(vec.z, FE_PARAM_DOUBLE, "vec_z");
ADD_PARAMETER(condition, FE_PARAM_DOUBLE, "condition");
END_PARAMETER_LIST();

mat3d ArcCosGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d rv(acos(in[0][0]), acos(in[0][1]), acos(in[0][2]),
		acos(in[1][0]), acos(in[1][1]), acos(in[1][2]),
		acos(in[2][0]), acos(in[2][1]), acos(in[2][2]));
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d ArcSinGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d rv(asin(in[0][0]), asin(in[0][1]), asin(in[0][2]),
		asin(in[1][0]), asin(in[1][1]), asin(in[1][2]),
		asin(in[2][0]), asin(in[2][1]), asin(in[2][2]));
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d CosGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d rv(cos(in[0][0]), cos(in[0][1]), cos(in[0][2]),
		cos(in[1][0]), cos(in[1][1]), cos(in[1][2]),
		cos(in[2][0]), cos(in[2][1]), cos(in[2][2]));
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d SinGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	mat3d rv(sin(in[0][0]), sin(in[0][1]), sin(in[0][2]),
		sin(in[1][0]), sin(in[1][1]), sin(in[1][2]),
		sin(in[2][0]), sin(in[2][1]), sin(in[2][2]));
	return GGP::Operation(rv, fin, mat, tip);
}

mat3d SetterGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	return GGP::Operation(m, invec, mat, tip);
}

BEGIN_PARAMETER_LIST(SetterGGP, GGP)
ADD_PARAMETER(m[0][0], FE_PARAM_DOUBLE, "m11");
ADD_PARAMETER(m[0][1], FE_PARAM_DOUBLE, "m12");
ADD_PARAMETER(m[0][2], FE_PARAM_DOUBLE, "m13");

ADD_PARAMETER(m[1][0], FE_PARAM_DOUBLE, "m21");
ADD_PARAMETER(m[1][1], FE_PARAM_DOUBLE, "m22");
ADD_PARAMETER(m[1][2], FE_PARAM_DOUBLE, "m23");

ADD_PARAMETER(m[2][0], FE_PARAM_DOUBLE, "m31");
ADD_PARAMETER(m[2][1], FE_PARAM_DOUBLE, "m32");
ADD_PARAMETER(m[2][2], FE_PARAM_DOUBLE, "m33");


ADD_PARAMETER(invec.x, FE_PARAM_DOUBLE, "in_x");
ADD_PARAMETER(invec.y, FE_PARAM_DOUBLE, "in_y");
ADD_PARAMETER(invec.z, FE_PARAM_DOUBLE, "in_z");
END_PARAMETER_LIST();


mat3d MatrixSetterGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	return GGP::Operation(m, fin, mat, tip);
}

BEGIN_PARAMETER_LIST(MatrixSetterGGP, GGP)
ADD_PARAMETER(m[0][0], FE_PARAM_DOUBLE, "m11");
ADD_PARAMETER(m[0][1], FE_PARAM_DOUBLE, "m12");
ADD_PARAMETER(m[0][2], FE_PARAM_DOUBLE, "m13");

ADD_PARAMETER(m[1][0], FE_PARAM_DOUBLE, "m21");
ADD_PARAMETER(m[1][1], FE_PARAM_DOUBLE, "m22");
ADD_PARAMETER(m[1][2], FE_PARAM_DOUBLE, "m23");

ADD_PARAMETER(m[2][0], FE_PARAM_DOUBLE, "m31");
ADD_PARAMETER(m[2][1], FE_PARAM_DOUBLE, "m32");
ADD_PARAMETER(m[2][2], FE_PARAM_DOUBLE, "m33");
END_PARAMETER_LIST();

mat3d MatrixInverseGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	return GGP::Operation(in.inverse(), fin, mat, tip);
}

mat3d UnitGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	vec3d temp(in[0][0], in[1][1], in[2][2]);
	temp.unit();
	in[0][0] = temp.x;
	in[1][1] = temp.y;
	in[2][2] = temp.z;
	return GGP::Operation(in, fin, mat, tip);
}

mat3d AssertGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	assert(in[0][0] <= (m[0][0] + tolerance));
	assert(in[0][0] >= (m[0][0] - tolerance));
	assert(in[0][1] <= (m[0][1] + tolerance));
	assert(in[0][1] >= (m[0][1] - tolerance));
	assert(in[0][2] <= (m[0][2] + tolerance));
	assert(in[0][2] >= (m[0][2] - tolerance));

	assert(in[1][0] <= (m[1][0] + tolerance));
	assert(in[1][0] >= (m[1][0] - tolerance));
	assert(in[1][1] <= (m[1][1] + tolerance));
	assert(in[1][1] >= (m[1][1] - tolerance));
	assert(in[1][2] <= (m[1][2] + tolerance));
	assert(in[1][2] >= (m[1][2] - tolerance));

	assert(in[2][0] <= (m[2][0] + tolerance));
	assert(in[2][0] >= (m[2][0] - tolerance));
	assert(in[2][1] <= (m[2][1] + tolerance));
	assert(in[2][1] >= (m[2][1] - tolerance));
	assert(in[2][2] <= (m[2][2] + tolerance));
	assert(in[2][2] >= (m[2][2] - tolerance));

	return GGP::Operation(in, fin, mat, tip);
}

BEGIN_PARAMETER_LIST(AssertGGP, GGP)
ADD_PARAMETER(m[0][0], FE_PARAM_DOUBLE, "m11");
ADD_PARAMETER(m[0][1], FE_PARAM_DOUBLE, "m12");
ADD_PARAMETER(m[0][2], FE_PARAM_DOUBLE, "m13");

ADD_PARAMETER(m[1][0], FE_PARAM_DOUBLE, "m21");
ADD_PARAMETER(m[1][1], FE_PARAM_DOUBLE, "m22");
ADD_PARAMETER(m[1][2], FE_PARAM_DOUBLE, "m23");

ADD_PARAMETER(m[2][0], FE_PARAM_DOUBLE, "m31");
ADD_PARAMETER(m[2][1], FE_PARAM_DOUBLE, "m32");
ADD_PARAMETER(m[2][2], FE_PARAM_DOUBLE, "m33");


ADD_PARAMETER(tolerance, FE_PARAM_DOUBLE, "tolerance");
END_PARAMETER_LIST();



BEGIN_PARAMETER_LIST(NodalDataGGP, GGP)
ADD_PARAMETER(offset, FE_PARAM_INT, "offset");
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");

END_PARAMETER_LIST();



mat3d NodalDataGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	double val = 0.0; 
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			int  ni = se->m_node[j];
			val += data[ni] * arr[j];
		}
	}
	in[0][0] = val;
	return GGP::Operation(in, fin, mat, tip);
}



BEGIN_PARAMETER_LIST(NodalDataGradientGGP, GGP)
ADD_PARAMETER(offset, FE_PARAM_INT, "offset");
ADD_PARAMETER(field_name, FE_PARAM_STRING, "field_name");

END_PARAMETER_LIST();



mat3d NodalDataGradientGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	double val[FEElement::MAX_NODES];
	FESolidDomain * d = tip.pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(tip.pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			
			int  ni = se->m_node[j];
			val[j] = data[ni] * arr[j];
		}
		vec3d grad = d->gradient(*se, val, se->GaussPoints());
		in[0][0] = grad.x;
		in[1][1] = grad.y;
		in[2][2] = grad.z;
	}
	
	return GGP::Operation(in, fin, mat, tip);
}


mat3d DirectionChangeGGP::Operation(mat3d in, vec3d fin, FEAngioMaterialBase* mat, Segment::TIP& tip)
{
	vec3d temp = fin;
	vec3d temp1(1, 1, 1);
	temp = in * temp;
	temp.unit();

	mat3d mi = vector->Operation(in, fin, mat, tip);
	temp1 = mi * temp1;
	temp1.unit();
	if(acos(temp*temp1) > angle)
	{
		in *= mat3d(-1, 0, 0,
			0, -1, 0,
			0, 0, -1);
	}

	return GGP::Operation(in, fin, mat, tip);
}

void DirectionChangeGGP::Update()
{
	if (vector)
	{
		vector->Update();
	}
	GGP::Update();
}
void DirectionChangeGGP::SetCulture(Culture * cp)
{
	if (vector)
		vector->SetCulture(cp);

	GGP::SetCulture(cp);
}

BEGIN_PARAMETER_LIST(DirectionChangeGGP, GGP)
ADD_PARAMETER(angle, FE_PARAM_DOUBLE, "angle");
END_PARAMETER_LIST();

*/