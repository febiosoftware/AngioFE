#pragma once
#include <FECore/FECoreKernel.h>
#include <FECore/FECoreTask.h>

//
class FEAngio;

//! Allows angiofe to setup all the items it needs for the simulations
class AngioFETask : public FECoreTask
{
public:
	//! constructor
	explicit AngioFETask(FEModel* pfem);
	~AngioFETask(void);

	//! initalize the task. 
	bool Init(const char* szfile) override;

	//! run the task
	bool Run() override;

private:
	FEAngio*	m_pangio = nullptr;
};
