///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "Culture.h"
#include "angio3d.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include <iostream>

#include "FEAngio.h"

using namespace std;

//-----------------------------------------------------------------------------
Fileout::Fileout()
{
    logstream.open("out_log.ang");
    //stream3 = fopen("tracking.ang","wt");   // tracking.ang: time step, model time, total length in culture, number of branches in culture
	m_stream = fopen("out_data.ang","wt");                                        // data.ang: Store 3D coordinates of begining and end of each vessel segment
																			// as well as total length of the segment}
	m_stream2 = fopen("out_vess_state.ang","wt");						// Open the stream for the vessel state data file		
	bf_stream = fopen("out_bf_state.ang","wt");						// Open the stream for the body force state data file

	time_stream = fopen("out_time.ang","wt");						// Open the stream for the time and state data file
	time_write_headers = true;										// Set the time and state data file to write the headers on its first output
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	fclose(m_stream2);
	fclose(bf_stream);
	fclose(time_stream);
}

//-----------------------------------------------------------------------------
void Fileout::writeTracking(FEAngio& angio)
{
    fprintf(m_stream3,"%-12.7f %-12.7f %-12.7f %-5i\n",angio.m_dt,angio.m_t,angio.m_total_length,angio.m_num_branches);   // Write to tracking.ang
    
    return;
}

//-----------------------------------------------------------------------------
void Fileout::closeTracking()
{
    fclose(m_stream3);                                                        // Close stream to 'tracking.ang' (stream3) 
    
    return;
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio)
{
    cout << endl << "Time: " << angio.m_t << endl;                             // Print out current time to user
	//cout << "dt: " << data.dt << endl;
    cout << "Segments: " << angio.m_nsegs << endl;                             // Print out current number of segments to user
	cout << "Total Length: " << angio.m_total_length << endl;                  // Print out the current total length to user
	cout << "Branch Points: " << angio.m_num_branches << endl;                 // Print out the current number of branches to user
	cout << "Anastomoses: " << angio.cult.m_num_anastom << endl << endl;            // Print out the current number of anastomoses to user
    
    logstream << endl << "Time: " << angio.m_t << endl;                        // Print out current time to log file
	//logstream << "dt: " << data.dt << endl;
    logstream << "Segments: " << angio.m_nsegs << endl;                        // Print out current number of segments to log file
	logstream << "Total Length: " << angio.m_total_length << endl;             // Print out the current total length to log file
	logstream << "Branch Points: " << angio.m_num_branches << endl;            // Print out the current number of branches to log file
	logstream << "Anastomoses: " << angio.cult.m_num_anastom << endl << endl;       // Print out the current number of anastomoses to log file
        
    return;
}

//-----------------------------------------------------------------------------
void Fileout::dataout(FEAngio &feangio)
{
    writeData(feangio);                                    // Create and write to 'data.ang'
	
    //writeGrid(data, grid);                              // Create and write to 'grid.ang'
   
    //writeNodes(angio.data, angio.grid);
    
    //writeEconn(angio.data, angio.grid);
    
    //writeBC(angio.grid);
    
    //writeGrad(data, grid);                              // Create and write to 'grad.ang'
    
    //writeAngle(frag);                                   // Create and write to 'angle.ang'

    printtime(feangio);                                        // Display the run-time to the user

    return;
}

//-----------------------------------------------------------------------------
void Fileout::writeData(FEAngio &feangio)
{
	list<Segment>::iterator it;
	
	fprintf(m_stream,"%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-5s %-5s\n","Time","X1","Y1","Z1","X2","Y2","Z2","Length","Vess","Label");  // Write column labels to data.ang
	
	for (it = feangio.cult.m_frag.begin(); it != feangio.cult.m_frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		vec3d& r0 = it->m_tip[0].rt;
		vec3d& r1 = it->m_tip[1].rt;
		fprintf(m_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->seg_num,it->label);
	}
	fclose(m_stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeNodes(FEAngio& angio)
{
    Grid& grid = angio.grid;
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_nodes.ang","wt");                                       
	Node node;

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i){
	    node = grid.nodes[i];
	    fprintf(stream2, "%-5.2i %-12.7f %-12.7f %-12.7f\n", node.id, node.rt.x, node.rt.y, node.rt.z);
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeEconn(FEAngio& angio)
{
/*    /// File output: 'econn.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_econn.ang","wt");                                       
	Elem elem;
	
	for (int i = 0; i < grid.Ne; ++i){
	    elem = grid.ebin[i];
	    fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, elem.n1.id, elem.n2.id, elem.n3.id, elem.n4.id, elem.n5.id, elem.n6.id,  elem.n7.id,  elem.n8.id);
	}  
	                                                                      	
	fclose(stream2); */                                                       
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeCollFib(Grid &grid, bool initial)
{
	FILE *node_stream;
	
	if (initial == true)
		node_stream = fopen("out_coll_fib_init.ang","wt");
	else
		node_stream = fopen("out_coll_fib.ang","wt");

	//fprintf(node_stream,"%-5s %-12s %-12s %-12s %-12s %-12s\n","ID"," X "," Y "," Z ","THETA","ETA");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].collfib.x, grid.nodes[i].collfib.y, grid.nodes[i].collfib.z);
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeECMDen(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den.ang","wt");

	//fprintf(node_stream,"%-5s %-12s %-12s %-12s %-12s %-12s\n","ID"," X "," Y "," Z ","ECM_DEN","ECM_DEN0");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].ecm_den, grid.nodes[i].ecm_den0);
	}

	fclose(node_stream);

	return;
}


//-----------------------------------------------------------------------------
void Fileout::writeECMDenGrad(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den_grad.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].ecm_den_grad.x, grid.nodes[i].ecm_den_grad.y, grid.nodes[i].ecm_den_grad.z);
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeECMDenStore(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_density_store.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i ", grid.nodes[i].id+1);

		for (int j = 0; j < grid.nodes[i].ecm_den_store.size(); j++)
			fprintf(node_stream,"%-12.7f ",grid.nodes[i].ecm_den_store[j]);
		
		fprintf(node_stream,"\n");
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeECMFibrilStore(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_fibril_store.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i ", grid.nodes[i].id+1);

		for (int j = 0; j < grid.nodes[i].ecm_den_store.size(); j++)
			fprintf(node_stream,"%-12.7f %-12.7f %-12.7f ",grid.nodes[i].ecm_fibril_store[j].x,grid.nodes[i].ecm_fibril_store[j].y,grid.nodes[i].ecm_fibril_store[j].z);
		
		fprintf(node_stream,"\n");
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeBC(Grid &grid)
{
    /// File output: 'eBC.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_eBC.ang","wt");                                       
	
	int BC_violate[6] = {0};

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){
	    Elem& elem = grid.ebin[i];
	    
	    for (int j = 0; j < 6; ++j)
	        BC_violate[j] = 0;
	    
	    if ((elem.f1.BC == true) || (elem.f2.BC == true) || (elem.f3.BC == true) || (elem.f4.BC == true) || (elem.f5.BC == true) || (elem.f6.BC == true)){
	        if (elem.f1.BC == true)
	            BC_violate[0] = 1;
	        if (elem.f2.BC == true)
	            BC_violate[1] = 1;
	        if (elem.f3.BC == true)
	            BC_violate[2] = 1;
	        if (elem.f4.BC == true)
	            BC_violate[3] = 1;
	        if (elem.f5.BC == true)
	            BC_violate[4] = 1;
	        if (elem.f6.BC == true)
	            BC_violate[5] = 1;   
	        
	        fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, BC_violate[0], BC_violate[1], BC_violate[2], BC_violate[3], BC_violate[4], BC_violate[5]);}
	}  
	                                                                      	
	fclose(stream2);                                                        
}


//-----------------------------------------------------------------------------
void Fileout::printtime(FEAngio& angio)
{
	time_t stop;
    time(&stop);                                                // Stop the timer
	
	double t_seconds = (double) difftime(stop, angio.m_start);                 // Calculate the simulation time in seconds
	
	cout << endl << "Simulation time: " << t_seconds << " seconds (" 
	    << floor(t_seconds/60) << " minutes)." << endl << endl;                // Show the user how long the simulation took (in seconds)
    
    logstream << endl << "Simulation time: " << t_seconds << " seconds (" << floor(t_seconds/60) << " minutes)." << endl << endl;  
}

//-----------------------------------------------------------------------------
void Fileout::printrandseed(int randseed)
{
	logstream << endl << "Rand seed:" << randseed << endl << endl;
}

//-----------------------------------------------------------------------------
void Fileout::writeSegConn(list<Segment> &frag)
{
	list<Segment>::iterator it;
	
	///// File output: 'out_seg_conn.ang' /////
	
	FILE *stream;                                                           // Open stream to 'out_seg_conn.ang' (stream)
	stream = fopen("out_seg_conn.ang","wt");                                        
	
	fprintf(stream,"%-5s %-5s %-5s %-5s %-5s\n","SNum","Node1","     ","Node2","     ");  // Write column labels to data.ang

	for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		fprintf(stream,"%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i\n",it->seg_num,it->seg_conn[0][0],it->seg_conn[0][1],it->seg_conn[1][0],it->seg_conn[1][1]);  // Write to seg_conn.ang
	
	}
	
	fclose(stream);                                                         // Close stream to 'data.ang' (stream)
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	list<Segment>::iterator it;													// Iterator for the segment container FRAG
		
	fprintf(m_stream2,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length");  // Write column labels to out_vess_state.ang
	
	for (it = angio.cult.m_frag.begin(); it != angio.cult.m_frag.end(); ++it)								// Iterate through all segments in frag list container (it)
	{
		vec3d& r0 = it->m_tip[0].rt;
		vec3d& r1 = it->m_tip[1].rt;
		fprintf(m_stream2,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",angio.FE_state,it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length);  // Write to out_vess_state.ang
	}
	
	return;
}

//-----------------------------------------------------------------------------
// Save positions of the body forces at the current time step (This function needs to be re-written)
void Fileout::save_bdy_forces(FEAngio& angio)
{
	FEModel& fem = angio.GetFEModel();
	int NBF = fem.BodyLoads();

	fprintf(bf_stream,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length"); 

	for (int i = 0; i < NBF; ++i)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
		FEParameterList& pl = pbf->GetParameterList();
		FEParam* pa = pl.Find("a");
		FEParam* prc = pl.Find("rc");

		if (pa && prc)
		{
			if (pa->value<double>() != 0.0)
				fprintf(bf_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",angio.FE_state, angio.m_t, prc->value<vec3d>().x, prc->value<vec3d>().y, prc->value<vec3d>().z, prc->value<vec3d>().x + 1.0, prc->value<vec3d>().y + 1.0, prc->value<vec3d>().z + 1.0, 1.73); 
		}
	}

	return;
}

//-----------------------------------------------------------------------------
// Save the current time information			
void Fileout::save_time(FEAngio& angio)
{
	if (time_write_headers == true){												// If this is the first time writing to out_time.ang
		fprintf(time_stream,"%-5s %-12s %-12s\n","State","Vess Time","FE Time");		// Print the column labels
		time_write_headers = false;}													// Turn off the headers flag
	
	fprintf(time_stream,"%-5.2i %-12.7f %-12.7f\n",angio.FE_state, angio.m_t, ((double)angio.FE_state - 1.0)*angio.FE_time_step);	// Print out the FE state, the vessel growth model time, and the FE time
	
	return;
}
