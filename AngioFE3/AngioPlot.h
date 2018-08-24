#pragma once
#include <FECore/FEPlotData.h>

//base class for computationally intensive plot variables that cannot be cached
//-----------------------------------------------------------------------------
class CIFEDomainData : public FEDomainData
{
public:
	explicit CIFEDomainData(Var_Type t, Storage_Fmt s) : FEDomainData(t, s) {}
	virtual void Update(FEModel* model) = 0;//update any cached data as this fucntion

};

//-----------------------------------------------------------------------------
class FEPlotAngioStress : public FEDomainData
{
public:
	explicit FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixStress : public FEDomainData
{
public:
	explicit FEPlotMatrixStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixViscoStress : public FEDomainData
{
public:
	explicit FEPlotMatrixViscoStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixElasticStress : public FEDomainData
{
public:
	explicit FEPlotMatrixElasticStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixElastic_m_Q : public FEDomainData
{
public:
	explicit FEPlotMatrixElastic_m_Q(FEModel* pfem) : FEDomainData(PLT_MAT3F, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotVesselStress : public FEDomainData
{
public:
	explicit FEPlotVesselStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotVesselWeight : public FEDomainData
{
public:
	explicit FEPlotVesselWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixWeight : public FEDomainData
{
public:
	explicit FEPlotMatrixWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixTangent : public FEDomainData
{
public:
	explicit FEPlotMatrixTangent(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotBranches : public FEDomainData
{
public:
	explicit FEPlotBranches(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAnastamoses : public FEDomainData
{
public:
	explicit FEPlotAnastamoses(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotSegmentLength : public FEDomainData
{
public:
	explicit FEPlotSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotRefSegmentLength : public FEDomainData
{
public:
	explicit FEPlotRefSegmentLength(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMDensity : public FEDomainData
{
public:
	explicit FEPlotAngioECMDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str)  override;
};

//-----------------------------------------------------------------------------
class FEPlotPrimaryVesselDirection : public FEDomainData
{
public:
	explicit FEPlotPrimaryVesselDirection(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str)  override;
private:
	//maps element id to segment direction in the reference configuration
	std::unordered_map<int, vec3d> accumulated_seg_direction;
};
