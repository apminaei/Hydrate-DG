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
#include <exception>
#include <chrono>
#include <stdio.h>

#include "dune/Hydrate-DG/IncludesDUNE.hh"
#include "dune/Hydrate-DG/akerbp2D_pockmark/include_problem.hh"


//#define ALUGRID		already defined in IncludeDUNE.hh
//#define PARALLEL		already defined in IncludeDUNE.hh

int main(int argc, char **argv)
{
	try
	{
		// Maybe initialize MPI
		 Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
	    if(helper.size()==1){
	    std::cout << "This is a test of project hydrate DG. " << std::endl;
	    }
	    if(Dune::MPIHelper::isFake){
	      std::cout<< "This is a sequential program." << std::endl;
	    }
	    else {
	    	if(helper.size()==1){
	    		std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
	    	}
	    }

		if (argc != 2)
		{
			if (helper.rank() == 0)
			{
				std::cout << "usage: ./hydrateDG-Ex1 <input_file.ini> " << std::endl;
			}
			return 1;
		}

		int level = 1;

		char input[80];
	    sscanf(argv[1],"%39s", input);
	    std::string input_file = "/home/peiravim/dune/Hydrate-DG/dune/Hydrate-DG/akerbp2D_pockmark/inputs/";
	    input_file += input;
	    std::cout<< "input file: " << input_file << std::endl ;

	    Dune::ParameterTree ptree;
	    Dune::ParameterTreeParser ptreeparser;
	    ptreeparser.readINITree(input_file,ptree);
	    ptreeparser.readOptions(argc,argv,ptree);

		/**************************************************************************************************/
		// MESH
	    MeshParameters<Dune::ParameterTree> mesh(ptree);
	    const int dim = mesh.dimension;

		/*____________________________________________*/

		
		const double X = mesh.X_length;
		//const double Y = mesh.Y_length;
		const double Z = mesh.Z_length;
		/*____________________________________________*/

		
		Dune::FieldVector<double, dim> L(0.0); // L represents the right top node of the rectangular/cuboidal domain
		L[0] = X;
		if (dim == 2)
		{
			L[1] = Z;
		}
		// else if (dim == 3)
		// {
		// 	L[1] = Y;
		// 	L[2] = Z;
		// }
		std::array<int, dim> N(Dune::filledArray<dim, int>(1));
		N[0] = mesh.X_cells * level; //gridParams.gridScaleRatio;
		if (dim == 2)
			N[1] = mesh.Z_cells * level;
		// else if (dim == 3)
		// {
		// 	N[1] = parameters.Y_cells; // * level;
		// 	N[2] = parameters.Z_cells; // * level;
		// }
#ifdef YASP
		typedef Dune::YaspGrid<dim> Grid;
		std::bitset<dim> periodic(false);

		int overlap = 1;
		std::shared_ptr<Grid> grid = std::shared_ptr<Grid>(new Grid(L, N, periodic, overlap, Dune::MPIHelper::getCollectiveCommunication()));
		grid->refineOptions(false); // keep overlap in cells

		typedef Grid::LeafGridView GV;
		const GV &gv = grid->leafGridView();
		grid->loadBalance();
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
		grid->loadBalance();
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
	// auto bxlambda = [&](const auto& x){ // Dirichlet for x-component of velocity
	// 	if (x[0]>X-1.e-6 || x[1]>Z-1.e-6 || x[0]<1.e-6 || x[1]<1.e-6) return true;
	// 	return false;
    // };
    // auto bx = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bxlambda);
		Dune::VTKWriter<GV> vtkWriter(gv);
  		vtkWriter.write(std::string("gridviews"));
		driver(gv, ptree, helper);
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
