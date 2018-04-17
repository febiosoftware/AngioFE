#pragma once
#include <FECore/vec3d.h>
#include "FECore/FEElement.h"
#include "FECore/FESolidDomain.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "Segment.h"
#include "CultureParameters.h"
#include "ECMInitializer.h"
#include <FECore/FENormalProjection.h>
#include "FiberManager.h"
#include "CommonAngioProperites.h"


class FEAngioMaterialPoint;

//contains the shared functionality, material pointers may need to be passed in to get portions of the functionality
class FEAngioMaterialBase
{
public:
	FEAngioMaterialBase();
	virtual ~FEAngioMaterialBase(){}

	// calculate the current spatial position, given an element and local coordinates
	vec3d CurrentPosition(FESolidElement * pe, double r, double s, double t) const;

	void AdjustMeshStiffness(FEMaterial* mat);

	void UpdateFiberManager();

	//! Assign a grid
	void SetFEAngio(FEAngio* pangio) { m_pangio = pangio; }

	void UpdateSproutStressScaling();

	void Update();

	bool Overwrite() const;

	double GetAnisotropy() const;

	//begin virtual functions
	virtual void InitializeFibers()=0;

	virtual mat3ds AngioStress(FEAngioMaterialPoint& mp)=0;

	virtual void FinalizeInit()=0;

	virtual void UpdateECM()=0;

	virtual void UpdateGDMs()=0;

	virtual bool InitECMDensity(FEAngio * angio)=0;

	virtual void ApplySym()=0;

	virtual void SetupSurface() = 0;

	virtual void UpdateAngioStresses() = 0;

	virtual FEMaterial * GetMatrixMaterial()=0;

	virtual CommonAngioProperties * GetCommonAngioProperties()=0;

	virtual int GetID_ang() = 0;

	virtual FEMaterial * GetMaterial() = 0;

	FEAngio * m_pangio;
	CultureParameters m_cultureParams;

	ECMInitializer * ecm_initializer;

	double scale;

	int sym_planes[7];

	bool sym_on;

	vec3d sym;
	double sym_vects[7][3];

	std::vector<int> domains;
	std::vector<FEDomain*> domainptrs;
	std::unordered_map<int, int> node_map;//maps global node number to domain local node numbers
	std::unordered_map<FEDomain *, int> meshOffsets;
	FESurface * exterior_surface;
	FENormalProjection * normal_proj;

	FiberManager * fiber_manager;

	int mat_id;

protected:
};