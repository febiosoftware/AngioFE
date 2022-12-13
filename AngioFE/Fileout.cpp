///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include <unordered_set>
#include <FECore/FEAnalysis.h>
#include <FECore/FSPath.h>
#include <FEBioLib/FEBioModel.h>
#include <FEBioLib/FEBioConfig.h>
#include <FEBioLib/febiolib_api.h>

// needed for GetInputFileName()
#pragma comment(lib, "psapi.lib")

using namespace std;

string Fileout::m_sfile;

//-----------------------------------------------------------------------------
Fileout::Fileout(FEAngio& angio)
{
	m_sfile = angio.m_fem->GetInputFileName();
	size_t dot = m_sfile.rfind('.');
	m_sfile = m_sfile.substr(0, dot);

	logstream.open(m_sfile + "_log.csv");
	//write the headers
	logstream << "Time,Material,Segments,Total Length,Vessels,Branches,Anastamoses,Active Tips" << endl;

	// write the line file
	vessel_state_stream = fopen((m_sfile + ".ang2").c_str(), "wb");//check the parameters consider setting the compression level
	// initialize version and line numbers
	unsigned int magic = 0xfdb97531;
	unsigned int version = 1;
	unsigned int num_bitmasks = angio.m_fem->Materials() / 32 + 1;
	// write the magic number, version, and number of materials
	fwrite(&magic, sizeof(unsigned int), 1, vessel_state_stream);
	fwrite(&version, sizeof(unsigned int), 1, vessel_state_stream);
	fwrite(&num_bitmasks, sizeof(unsigned int), 1, vessel_state_stream);
	// for each material
	for (unsigned int j = 0; j < num_bitmasks; j++)
	{
		unsigned int c_bitmask = 0;
		unsigned int place_holder = 1;
		for (int i = 0; i < 32; i++)
		{
			int index = i + j * 32;
			if (index == angio.m_fem->Materials())
				break;
			FEMaterial* mat = angio.m_fem->GetMaterial(index);
			FEAngioMaterial* angio_mat = angio.GetAngioComponent(mat);
			if (angio_mat)
			{
				c_bitmask |= place_holder;
			}
			place_holder = place_holder << 1;
		}
		fwrite(&c_bitmask, sizeof(unsigned int), 1, vessel_state_stream);
	}

	feangio_state_stream = fopen((m_sfile + "_time_stats.csv").c_str(), "wt");
	fprintf(feangio_state_stream, "%-64s,%-64s,%-64s,%-64s\n",
		"Timestep", "Growth Process", "Branch Policy Update", "Update Stress");
	cell_state_stream = fopen((m_sfile + "_cells.txt").c_str(), "wt"); // cells
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
	logstream.close();
	fclose(vessel_state_stream);
	fclose(cell_state_stream);
	fclose(feangio_state_stream);
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio, double time)
{
	//just store the status in a csv
	//logstream << "Time,Segments,Total Length,Vessels,Branchs,Anastamoses,Active Tips" << endl;
	std::vector <int> SC = getSegmentCount_pm(angio);
	std::vector <double> SL = getSegmentLength_pm(angio, time);
	std::vector <int> BC = getBranchCount_pm(angio);
	std::vector <int> TC = getTipCount_pm(angio);
	for (size_t i = 0; i < angio.m_fem->Materials(); i++) {
		logstream << angio.GetFEModel()->GetTime().currentTime << "," << angio.m_fem->GetMaterial(i)->GetName() << "," << SC[i] << "," << SL[i] << ",," << BC[i] << ",," << TC[i] << std::endl;
	}
}

void PrintSegment(vec3d r0, vec3d r1)
{
#ifndef NDEBUG
	//std::cout << "r0: " << r0.x << "," << r0.y << "," << r0.z << std::endl;
	//std::cout << "r1: " << r1.x << "," << r1.y << "," << r1.z << std::endl;
#endif
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	unsigned int new_seg_count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		new_seg_count += angio.angio_elements[i]->recent_segments.size();
	}

	//write segcount and time
	fwrite(&new_seg_count, sizeof(unsigned int), 1, vessel_state_stream);
	auto ti = angio.GetFEModel()->GetTime();
	float rtime = static_cast<float>(ti.currentTime);
	fwrite(&rtime, sizeof(float), 1, vessel_state_stream);// == sizeof(float); SL: Not sure what the equivalency was doing

	unsigned int crc_segcount = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		for (int j = 0; j < angio.angio_elements[i]->recent_segments.size(); j++)
		{
			crc_segcount++;
			Tip* back_tip = angio.angio_elements[i]->recent_segments[j]->back;
			Tip* front_tip = angio.angio_elements[i]->recent_segments[j]->front;
			vec3d r0 = angio.ReferenceCoordinates(back_tip);
			vec3d r1 = angio.ReferenceCoordinates(front_tip);

			PrintSegment(r0, r1);
			//consider checking fwrites return values
			float cur = static_cast<float>(r0.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			cur = static_cast<float>(r1.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			r1 = angio.ReferenceCoordinates(front_tip);
		}
	}
	assert(crc_segcount == new_seg_count);
}

void Fileout::bulk_save_vessel_state(FEAngio& angio)
{
	unsigned int new_seg_count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		new_seg_count += angio.angio_elements[i]->grown_segments.size();
	}

	//write segcount and time
	fwrite(&new_seg_count, sizeof(unsigned int), 1, vessel_state_stream);
	auto ti = angio.GetFEModel()->GetTime();
	float rtime = static_cast<float>(ti.currentTime);
	fwrite(&rtime, sizeof(float), 1, vessel_state_stream);// == sizeof(float); SL: Not sure what this equivalency is doing

	unsigned int crc_segcount = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		for (int j = 0; j < angio.angio_elements[i]->grown_segments.size(); j++)
		{
			crc_segcount++;
			Tip* back_tip = angio.angio_elements[i]->grown_segments[j]->back;
			Tip* front_tip = angio.angio_elements[i]->grown_segments[j]->front;
			vec3d r0 = angio.ReferenceCoordinates(back_tip);
			vec3d r1 = angio.ReferenceCoordinates(front_tip);

			PrintSegment(r0, r1);
			//consider checking fwrites return values
			float cur = static_cast<float>(r0.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			cur = static_cast<float>(r1.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			r1 = angio.ReferenceCoordinates(front_tip);
		}
	}
	assert(crc_segcount == new_seg_count);
}

//-----------------------------------------------------------------------------
// Save active points
void Fileout::save_active_tips(FEAngio& angio) const
{

}
void Fileout::save_timeline(FEAngio& angio)
{
#ifndef NDEBUG

#endif
}

void Fileout::save_final_vessel_csv(FEAngio& angio)
{
	FILE* final_vessel_file = fopen((m_sfile + "_vessels.csv").c_str(), "wt");
	assert(final_vessel_file);
	fprintf(final_vessel_file, "x0,y0,z0,x1,y1,z1,start time\n");
	FEMesh* mesh = angio.GetMesh();
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		for (int j = 0; j < angio.angio_elements[i]->grown_segments.size(); j++) {
			auto seg = *(angio.angio_elements[i]->grown_segments[j]);
			vec3d p0 = seg.front->GetPosition(mesh);
			vec3d p1 = seg.back->GetPosition(mesh);
			fprintf(final_vessel_file, "%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f\n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, seg.back->time);
		}

	}
	fclose(final_vessel_file);
}

void Fileout::save_final_cells_txt(FEAngio& angio)
{
	FILE* final_cell_file = fopen((m_sfile + "_cells.txt").c_str(), "at");
	assert(final_cell_file);
	FEMesh* mesh = angio.GetMesh();
	for (auto iter = angio.cells.begin(); iter != angio.cells.end(); iter++) {
		auto cell = iter->second;
		vec3d p = cell->GetPosition(mesh);
		int time = angio.GetFEModel()->GetCurrentStep()->m_ntimesteps;
		fprintf(final_cell_file, "%d,%d,%-12.5e,%-12.5e,%-12.5e", time, cell->initial_cell_id, p.x, p.y, p.z);
		//! Print solute values
		for (int isol = 0; isol < cell->Solutes.size(); isol++) {
			fprintf(final_cell_file, ",%-12.5e", cell->Solutes[isol]->GetInt());
		}
		//! Print SBM values
		for (int isbm = 0; isbm < cell->SBMs.size(); isbm++) {
			fprintf(final_cell_file, ",%-12.5e", cell->SBMs[isbm]->GetInt());
		}
		fprintf(final_cell_file, "\n");
	}
	fclose(final_cell_file);
}

void Fileout::save_initial_cells_txt(FEAngio& angio)
{
	FILE* final_cell_file = fopen((m_sfile + "_cells.txt").c_str(), "at");
	assert(final_cell_file);
	FEMesh* mesh = angio.GetMesh();
	fprintf(final_cell_file, "*timestep,cell id,X Pos,Y Pos,Z Pos");
	CellSpeciesManager* m_species = angio.angio_elements[0]->_angio_mat->cell_species_manager;
	if (m_species) {
		for (int isol = 0; isol < m_species->cell_solute_prop.size(); isol++) {
			int sol_id = m_species->cell_solute_prop[isol]->GetSoluteID();
			fprintf(final_cell_file, ",Sol %d", sol_id);
		}
		// assign the properties to each SBM. Initialize body loads and add to the SBMs container.
		for (int isbm = 0; isbm < m_species->cell_SBM_prop.size(); isbm++)
		{
			int sbm_id = m_species->cell_SBM_prop[isbm]->GetSBMID();
			fprintf(final_cell_file, ",SBM %d", sbm_id);
		}
	}
	fprintf(final_cell_file, "\n");
	for (auto iter = angio.cells.begin(); iter != angio.cells.end(); iter++)
	{
		auto cell = iter->second;
		vec3d p = cell->GetPosition(mesh);
		fprintf(final_cell_file, "%d,%d,%-12.5e,%-12.5e,%-12.5e", 0, cell->initial_cell_id, p.x, p.y, p.z);
		//! Print solute values
		for (int isol = 0; isol < cell->Solutes.size(); isol++) {
			fprintf(final_cell_file, ",%-12.5e", cell->Solutes[isol]->GetInt());
		}
		//! Print SBM values
		for (int isbm = 0; isbm < cell->SBMs.size(); isbm++) {
			fprintf(final_cell_file, ",%-12.5e", cell->SBMs[isbm]->GetInt());
		}
		fprintf(final_cell_file, "\n");
	}
	fclose(final_cell_file);
}

void Fileout::save_feangio_stats(FEAngio& angio)
{
	FETimeInfo& ti = angio.m_fem->GetTime();
	const int STR_SIZE = 64;
	char grow_time[STR_SIZE];
	char update_branch_policy[STR_SIZE];
	char update_stress[STR_SIZE];
	angio.grow_timer.time_str(grow_time);
	angio.update_branch_policy_timestep_timer.time_str(update_branch_policy);
	angio.update_angio_stress_timer.time_str(update_stress);
	fprintf(feangio_state_stream, "%-12.7f,%-64s,%-64s,%-64s\n",
		ti.currentTime, grow_time, update_branch_policy, update_stress);
	fflush(feangio_state_stream);
}

int Fileout::getBranchCount(FEAngio& angio)
{
	int count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		count += angio.angio_elements[i]->branch_count;
	}
	return count;
}

std::vector <int> Fileout::getBranchCount_pm(FEAngio& angio)
{
	std::vector <int> count;
	//	for (int i = 0; i < angio.m_pmat_ids.size(); i++) {
	for (int i = 0; i < angio.m_fem->Materials(); i++) {
		count.emplace_back(0);
	}
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		FESolidElement* se = angio.angio_elements[i]->_elem;
		int mat_id = se->GetMatID();
		count[mat_id] += angio.angio_elements[i]->branch_count;
	}
	return count;
}

int Fileout::getSegmentCount(FEAngio& angio)
{
	int count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		count += int(angio.angio_elements[i]->grown_segments.size());
	}
	return count;
}

std::vector <int> Fileout::getSegmentCount_pm(FEAngio& angio)
{
	std::vector <int> count;
//	for (int i = 0; i < angio.m_pmat_ids.size(); i++) {
	for (int i = 0; i < angio.m_fem->Materials(); i++) {
	count.emplace_back(0);
	}
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		FESolidElement* se = angio.angio_elements[i]->_elem;
		int mat_id = se->GetMatID();
		count[mat_id] += angio.angio_elements[i]->grown_segments.size();
	}
	return count;
}

double Fileout::getSegmentLength(FEAngio& angio, double time)
{
	FEMesh* mesh = angio.GetMesh();
	double len = 0;
#pragma omp parallel for
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		double temp = angio.angio_elements[i]->GetLengthAtTime(mesh, time);
#pragma omp critical
		len += temp;
	}
	return len;
}

std::vector <double> Fileout::getSegmentLength_pm(FEAngio& angio, double time)
{
	std::vector <double> len;
	//for (int i = 0; i < angio.m_pmat_ids.size(); i++) {
		for (int i = 0; i < angio.m_fem->Materials(); i++) {
		len.emplace_back(0);
	}
	FEMesh* mesh = angio.GetMesh();
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		FESolidElement* se = angio.angio_elements[i]->_elem;
		int mat_id = se->GetMatID();
		double temp = angio.angio_elements[i]->GetLengthAtTime(mesh, time);
#pragma omp critical
		len[mat_id] += temp;
	}
	return len;
}

int Fileout::getTipCount(FEAngio& angio)
{
	int count = 0;
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		for (auto iter = angio.angio_elements[i]->next_tips.begin(); iter != angio.angio_elements[i]->next_tips.end(); ++iter)
		{
			count += int(iter->second.size());
		}
	}
	return count;
}

std::vector <int> Fileout::getTipCount_pm(FEAngio& angio)
{
	std::vector <int> count;
	//	for (int i = 0; i < angio.m_pmat_ids.size(); i++) {
	for (int i = 0; i < angio.m_fem->Materials(); i++) {
		count.emplace_back(0);
	}
	for (int i = 0; i < angio.angio_elements.size(); i++)
	{
		FESolidElement* se = angio.angio_elements[i]->_elem;
		int mat_id = se->GetMatID();
		for (auto iter = angio.angio_elements[i]->next_tips.begin(); iter != angio.angio_elements[i]->next_tips.end(); ++iter)
			//for (auto iter = angio.angio_elements[i]->active_tips.begin(); iter != angio.angio_elements[i]->active_tips.end(); ++iter)
		{
			count[mat_id] += int(iter->second.size());
			//count[mat_id] += iter->second.size();
		}
	}
	return count;
}