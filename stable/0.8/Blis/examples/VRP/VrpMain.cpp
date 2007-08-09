/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <iostream>

#include "CoinError.hpp"
#include "CoinTime.hpp"

#include "BlisConfig.h"

#include "OsiSolverInterface.hpp"
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif

#include "VrpModel.h"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

//#############################################################################

int main(int argc, char *argv[]) 
{
    try{
	// Set up lp solver
#ifdef  COIN_HAS_CLP
        OsiClpSolverInterface lpSolver;
	lpSolver.getModelPtr()->setDualBound(1.0e10);
	lpSolver.messageHandler()->setLogLevel(0);
#endif
	
	// Create VRP model 
	VrpModel model;
	model.setSolver(&lpSolver);
	
#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, model);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, model); 
#endif

	// Search for best solution
	broker.search(&model);
	
	// Report the best solution found and its ojective value
	broker.printBestSolution();
    }
    catch(CoinError& er) {
	std::cerr << "\nBLIS ERROR: \"" << er.message() 
		  << "\""<< std::endl
		  << "             from function \"" << er.methodName()
		  << "\""<< std::endl
		  << "             from class \"" << er.className()
		  << "\"" << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }
    
    
    return 0;
}

//#############################################################################

