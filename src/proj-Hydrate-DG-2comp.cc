// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <time.h>

#include "dune/proj-Hydrate-DG-2comp/IncludesDUNE.hh"
#include "dune/proj-Hydrate-DG-2comp/IncludesProblem.hh"


//#define ALUGRID		already defined in IncludeDUNE.hh
//#define PARALLEL		already defined in IncludeDUNE.hh

int main(int argc, char **argv)
{
	try
	{
		// Maybe initialize MPI
		Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
		std::cout << "Hello World! This is Test-TwoPhaseFlow." << std::endl;
		if (Dune::MPIHelper::isFake)
			std::cout << "This is a sequential program." << std::endl;
		else
		{
			if (helper.rank() == 0)
			{
				std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
				std::cout << "parallel run on " << helper.rank() << " processor" << std::endl;
			}
		}

		if (argc != 5)
		{
			if (helper.rank() == 0)
			{
				std::cout << "usage: ./proj_Hydrate_DG_2comp <level> <dt> <t_END> <t_OP>" << std::endl;
			}
			return 1;
		}

		int level;
		sscanf(argv[1], "%d", &level); //argument 'refinement level'
		std::cout << "level=" << argv[1] << "\t";

		double dt;
		sscanf(argv[2], "%lg", &dt); //argument 'time step size'
		std::cout << "dt=" << argv[2] << "\t";

		double t_END;
		sscanf(argv[3], "%lg", &t_END); //argument 'end time'
		std::cout << "t_END=" << argv[3] << "\t";

		double t_OP;
		sscanf(argv[4], "%lg", &t_OP); //argument 'op Interval'
		std::cout << "t_OP=" << argv[4] << std::endl;

		/*____________________________________________*/

		IncludeProblemSpecifications::ProblemSpecifications parameters;
		const int dim = parameters.dimension;
		const double X = parameters.X_length;
		const double Y = parameters.Y_length;
		const double Z = parameters.Z_length;
		/*____________________________________________*/

		
		Dune::FieldVector<double, dim> L(0.0); // L represents the right top node of the rectangular/cuboidal domain
		L[0] = X;
		if (dim == 2)
		{
			L[1] = Z;
		}
		else if (dim == 3)
		{
			L[1] = Y;
			L[2] = Z;
		}
		std::array<int, dim> N(Dune::filledArray<dim, int>(level));
		N[0] = parameters.X_cells * level; //gridParams.gridScaleRatio;
		if (dim == 2)
			N[1] = parameters.Z_cells * level;
		else if (dim == 3)
		{
			N[1] = parameters.Y_cells; // * level;
			N[2] = parameters.Z_cells; // * level;
		}
#ifdef YASP
		typedef Dune::YaspGrid<dim> Grid;
		std::bitset<dim> periodic(false);

		int overlap = 2;
		std::shared_ptr<Grid> grid = std::shared_ptr<Grid>(new Grid(L, N, periodic, overlap, Dune::MPIHelper::getCollectiveCommunication()));
		grid->refineOptions(false); // keep overlap in cells

		typedef Grid::LeafGridView GV;
		const GV &gv = grid->leafGridView();
		
#elif UG

		typedef Dune::UGGrid<dim> Grid;
		//Grid(UGCollectiveCommunication comm =CollectiveCommunication<MPI_Comm>);
		auto ll = Dune::FieldVector<Grid::ctype, dim>{{0, 0}};
		auto ur = Dune::FieldVector<Grid::ctype, dim>{{L[0], L[1]}};
		std::array<unsigned int, dim> elements;
		elements[0] = N[0];
		elements[1] = N[1];

		//std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
		std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, elements);

		typedef Grid::LeafGridView GV;

		GV gv = grid->leafGridView();

#else ALUGRID
		typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid;
		auto ll = Dune::FieldVector<Grid::ctype, dim>{{0, 0}};
		auto ur = Dune::FieldVector<Grid::ctype, dim>{{L[0], L[1]}};
		std::array<unsigned int, dim> elements;
		elements[0] = N[0];
		elements[1] = N[1];

		//std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createSimplexGrid(ll, ur, elements);
		std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(ll, ur, elements); // load balance the grid
		
		
		typedef Grid::LeafGridView GV;
		GV gv = grid->leafGridView();

  		// Transfer partitioning from ParMETIS to our grid
  		grid->loadBalance();

		
#endif
		proj_Hydrate_SimplexDG(gv, dt, t_END, t_OP);
	}
	catch (Dune::Exception &e)
	{
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
