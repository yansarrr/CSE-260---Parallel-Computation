#include "stimulus.h"
#include "obstacle.h"
#include "compute.h"
#include "plotter.h"
#include <list>
#include <random>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <assert.h>
using namespace std;

#ifdef _MPI_
#include <mpi.h>
#include <malloc.h>

#endif
#include <math.h>

///
/// Compute object
///
/// upon construction, runs a simulation.
///
Compute::Compute(Buffers &u, Plotter *plt, ControlBlock &cb, const int _myRank,
		 const int seedv):
    u(u), plt(plt), cb(cb), myRank(_myRank), seedv(seedv), M(u.M), N(u.N), RANDSTIM(0)
{


#ifdef _MPI_
    // detect global edge
    topGlobalEdge = myRank / cb.px == 0;
    botGlobalEdge = myRank / cb.px == (cb.py - 1);
    leftGlobalEdge = myRank % cb.px == 0;
    rightGlobalEdge = myRank % cb.px == (cb.px - 1);

#else    
    topGlobalEdge = true;
    botGlobalEdge = true;
    leftGlobalEdge = true;
    rightGlobalEdge = true;
#endif
}

#ifdef _MPI_
// possible problems:
// initial conditions are only applied on inerior cells excluding ghost cells, thus we need update u.cur

void Compute::ghostExchange(Buffers& u) {
    // consts
    enum Direction {TOP, BOTTOM, LEFT, RIGHT, N_DIR};
    enum Role {SEND, RECV, N_ROLE};
    array<int, N_DIR> neighborAt = {myRank - cb.px, myRank + cb.px, myRank - 1, myRank + 1};
    array<bool, N_DIR> hasEdge = {topGlobalEdge, botGlobalEdge, leftGlobalEdge, rightGlobalEdge};

MPI_Barrier(MPI_COMM_WORLD);

    int nRequestDir = 4;

    for (bool edge : hasEdge) {
        if (edge) nRequestDir--;
    }

    MPI_Request request[nRequestDir * 2];

    // cout << logHead << "sync ghostCells (Barrier: avoid update data override unfinished calculation of partner)!" << endl;

    int iRequest = 0;
    // send to all neighbors
    // north/top neighbor
    // cout << logHead << "request to top neighbor..." << endl;
    vector<double *> pE_prev_tosend; // only consider LEFT/RIGHT

    array<double *, 4> ghostCells = {nullptr};

    for (int dir = 0; dir < N_DIR; dir++) {
        if (hasEdge[dir]) continue;
        int partner = neighborAt[dir];

        double *E_prev_tosend;
        int length_tosend;

        if (dir == TOP || dir == BOTTOM) {
            int send_i = (dir == TOP) ? 1 : u.gridM - 2;
            length_tosend = u.gridN - 2;

            E_prev_tosend = u.cur(send_i, 1);
        } else {
            int send_j = (dir == LEFT) ? 1 : u.gridN - 2;
            length_tosend = u.gridM - 2;

            E_prev_tosend = (double *) memalign(16, sizeof(double) * length_tosend);
            assert(E_prev_tosend != nullptr);
            for (int i = 1; i <= u.gridM - 2; i++) {
                E_prev_tosend[i - 1] = u.curV(i, send_j);
            } 
            pE_prev_tosend.push_back(E_prev_tosend);
        }

        ghostCells[dir] = (double *) memalign(16, sizeof(double) * length_tosend);
        assert(ghostCells[dir] != nullptr);
        MPI_Isend(E_prev_tosend, length_tosend, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &request[iRequest * 2 + SEND]);
        MPI_Irecv(ghostCells[dir], length_tosend, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &request[iRequest * 2 + RECV]);
        iRequest++;
    }

// sync
// cout << logHead << "waiting for immediate send and recv finish..." << endl;
    MPI_Waitall(nRequestDir * 2, request, MPI_STATUSES_IGNORE);

    // refill ghost cells
    for (int j = 1; j < u.gridN - 1; j++) {
        if (!topGlobalEdge) *u.cur(0, j) = ghostCells[TOP][j - 1];
        if (!botGlobalEdge) *u.cur(u.gridM - 1, j) = ghostCells[BOTTOM][j - 1];
    }

    for (int i = 1; i < u.gridM - 1; i++) {
        if (!leftGlobalEdge) *u.cur(i, 0) = ghostCells[LEFT][i - 1];
        if (!rightGlobalEdge) *u.cur(i, u.gridN - 1) = ghostCells[RIGHT][i - 1];
    }

    // free packing
    for (double * pE : pE_prev_tosend) {
        free(pE);
    }

    for (double* pG : ghostCells) {
    if (pG != nullptr) free(pG);
    }

MPI_Barrier(MPI_COMM_WORLD);
 }
#endif
    
///
/// Simulate
///
/// calls class specific calcU and calcEdgeU
///
void Compute::Simulate()
{
    const unsigned int t = 1;	 // timestep
    const unsigned int h = 1;    // grid space
    const double c = 0.29;   // velocity
    const double kappa = c*t/h;

    mt19937 generator(seedv);

    uniform_int_distribution<int> randRow(1, cb.m-10);
    uniform_int_distribution<int> randCol(1, cb.n-10);
    uniform_int_distribution<int> randEvent(0, 100);


    u.setAlpha((c*t/h) * (c*t/h));

    list<Stimulus *> sList;
    list<Obstacle *> oList;
    int iter = 0;

    ifstream f(cb.configFileName);
    if (cb.config.count("objects")){
	for (int i=0; i<cb.config["objects"].size(); i++){
	    auto const &ob = cb.config["objects"][i];
	    if (ob["type"] == "sine"){
		sList.push_back(new StimSine(u, ob["row"], ob["col"],
					     ob["start"], ob["duration"],
					     ob["period"]));
	    }else if (ob["type"] == "rectobstacle"){
		oList.push_back(new Rectangle(u, ob["row"], ob["col"],
					      ob["height"], ob["width"]));
	    }
	}
    }else{
	fprintf(stderr, "Using hardcoded stimulus\n");
	Rectangle obstacleA(u, cb.m/2+5, cb.n/2, 45, 5);
	Rectangle obstacleB(u, cb.m/2-50, cb.n/2, 45, 5);
	sList.push_back(new StimSine(u, cb.m/2, cb.n/3, 0 /*start*/, 500/*duration*/, 10 /*period*/));
    }

    ///
    /// generate stimulus
    ///
    /// once quiet (non-deterministic),
    /// we exit this loop and go into a loop that
    /// continues until iterations is exhausted
    ///
    while (!sList.empty() && iter < cb.niters){
	for (auto it = begin(sList); it!= end(sList);){
	    if (!(*it)->doit(iter)){
		delete *it;
		it = sList.erase(it);
	    }else{
		it++;
	    }
	}

#ifdef _MPI_
    // calc interior, then update boundary or ghost cells on next, and pass to cur and prv
    ghostExchange(u);
#endif

	calcU(u);


	calcEdgeU(u, kappa);

	if (cb.plot_freq && iter % cb.plot_freq == 0)
	    plt->updatePlot(iter, u.gridM, u.gridN);

	// DEBUG start
	//	u.printActive(iter);
	// DEBUG end

	u.AdvBuffers();


	iter++;
    }

    ///
    /// all stimulus done
    /// keep simulating till end
    ///
    for (;iter < cb.niters; iter++){

#ifdef _MPI_
    // calc interior, then update boundary or ghost cells on next, and pass to cur and prv
    ghostExchange(u);
#endif

	calcU(u);
	//	if (cb.plot_freq && iter % cb.plot_freq == 0)
	//	    plt->updatePlot(iter, u.gridM, u.gridN);

	calcEdgeU(u, kappa);


	if ((cb.plot_freq!=0) && (iter % cb.plot_freq == 0))
	    plt->updatePlot(iter, u.gridM, u.gridN);


	// DEBUG
	// u.printActive(iter);
	u.AdvBuffers();
    }
}

TwoDWave::TwoDWave(Buffers &u, Plotter *plt, ControlBlock &cb, const int _myRank,
		 const int seedv):
    Compute(u, plt, cb, _myRank, seedv){};

///
/// compute the interior cells
///
///
void TwoDWave::calcU(Buffers &u)
{

    // interior always starts at 2,2, ends at gridN
    for (int i=1; i<u.gridM-1; i++){
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(i,j) =
		u.alpV(i,j) *
		(u.curV(i-1,j) + u.curV(i+1,j) +
		 u.curV(i,j-1) + u.curV(i,j+1) - 4 * u.curV(i,j)) +
		2 * u.curV(i,j) - u.preV(i,j);
	}
    }
}

///
/// compute edges
///
/// compute interior edges. These are not ghost cells but cells that rely
/// on either ghost cell values or boundary cell values.
///
void TwoDWave::calcEdgeU(Buffers &u, const double kappa)
{
    // set the boundary conditions to absorbing boundary conditions (ABC)
    // du/dx = -1/c du/dt   x=0
    // du/dx = 1/c du/dt    x=N-1
    // conditions for an internal boundary (ie.g. ghost cells)
    // top edge


    // top global edge (instead of ghost cells)
    if (topGlobalEdge){
	// top row absorbing boundary condition
	int i = 0;
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(i,j) = u.curV(i+1,j) +
		((kappa-1)/(kappa+1)) * (u.nxtV(i+1,j) - u.curV(i,j));
	}
    }

    // bottom edge (instead of ghost cells)
    if (botGlobalEdge){
	int i = u.gridM-1;
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(i,j) = u.curV(i-1,j) +
		((kappa-1)/(kappa+1)) * (u.nxtV(i-1,j) - u.curV(i,j));
	}
    }

    // left edge
    if (leftGlobalEdge){
	int j = 0;
	for (int i=1; i<u.gridM-1; i++){
	    *u.nxt(i,j) = u.curV(i,j+1) +
		((kappa-1)/(kappa+1)) * (u.nxtV(i,j+1) - u.curV(i,j));
	}
    }
    // right edge
    if (rightGlobalEdge){
	int j = u.gridN-1;
	for (int i=1; i<u.gridM-1; i++){
	    *u.nxt(i,j) = u.curV(i,j-1) +
		((kappa-1)/(kappa+1)) * (u.nxtV(i,j-1) - u.curV(i,j));
	}
    }
}

//!
//! Use a different propgation model
//! This model shifts values in the horizontal direction
//!
DebugPropagate::DebugPropagate(Buffers &u, Plotter *plt, ControlBlock &cb, const int _myRank,
		 const int seedv):
    Compute(u, plt, cb, _myRank, seedv){};

//!
//! compute the interior cells
//!
void DebugPropagate::calcU(Buffers &u)
{

    // interior always starts at 2,2, ends at gridN-3
    for (int i=2; i<u.gridM-2; i++){
	for (int j=2; j<u.gridN-2; j++){
	    *u.nxt(i,j) = u.curV(i, j-1);
	}
    }
}

//!
//! compute edges
//! (either interior edges or global edges)
//!
void DebugPropagate::calcEdgeU(Buffers &u, const double kappa)
{
    if (topGlobalEdge){
	// top row absorbing boundary condition
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(1,j) = 0;
	}
    }else{
	int i = 1;
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(i,j) = u.curV(i, j-1);
	}
    }

    // bottom edge
    if (botGlobalEdge){
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(u.gridM-2,j) = 0;
	}
    }else{
	int i=u.gridM-2;
	for (int j=1; j<u.gridN-1; j++){
	    *u.nxt(i,j) = u.curV(i, j-1);
	}
    }

    // left edge
    if (leftGlobalEdge){
	for (int i=1; i<u.gridM-1; i++){
	    *u.nxt(i,1) = 0.0;
	}
    }else{
	int j=1;
	for (int i=1; i<u.gridM-1; i++){
	    *u.nxt(i,j) = u.curV(i, j-1);
	}
    }
    // right edge
    if (rightGlobalEdge){
	for (int i=1; i<u.gridM-1; i++){
	    // right column
	    *u.nxt(i,u.gridN-2) = 0.0;
	}
    }else{
	int j=u.gridN-2;
	for (int i=1; i<u.gridM-1; i++){
	    *u.nxt(i,j) = u.curV(i, j-1);
	}
    }
}


