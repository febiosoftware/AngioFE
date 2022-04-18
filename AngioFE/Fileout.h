///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "StdAfx.h"
//#include "zlib.h"

class FEAngio;

//! class that contains the streams for doing output
class Fileout
{
public:
	//! constructor
	Fileout(FEAngio& angio);
	virtual ~Fileout();
	//! print the status to the console
	void printStatus(FEAngio& angio, double time);
	//! save the vessel state
	void save_vessel_state(FEAngio& angio);
	//! saves all segments not just recent ones, should only be used in initialization
	void bulk_save_vessel_state(FEAngio& angio);
	//! saves all cells not just recent ones, should only be used in initialization
	void save_initial_cells_txt(FEAngio& angio);
	//! save the active tips
	void save_active_tips(FEAngio& angio) const;
	//! save the timeline of branch points
	void save_timeline(FEAngio& angio);

	//! timing statistics
	void save_feangio_stats(FEAngio& angio);
	//! save the final vascular network as a csv
	static void save_final_vessel_csv(FEAngio & angio);
	//! save the cell positions as a txt
	static void save_final_cells_txt(FEAngio & angio);

private:
	int getBranchCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getBranchCount_pm(FEAngio& angio);
	int getSegmentCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getSegmentCount_pm(FEAngio& angio);
	double getSegmentLength(FEAngio& angio, double time);
	// get per angio material
	std::vector <double> getSegmentLength_pm(FEAngio& angio, double time);
	int getTipCount(FEAngio& angio);
	// get per angio material
	std::vector <int> getTipCount_pm(FEAngio& angio);
	std::ofstream logstream;

	FILE*  vessel_state_stream = nullptr;
	FILE*  feangio_state_stream = nullptr;
	FILE*  cell_state_stream = nullptr;
};
