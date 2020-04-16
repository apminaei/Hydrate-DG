/*
 * proj_Hydrate_SimplexDG.hh
 *
 *  Created on: Sep 30, 2016
 *      Author: shubhangi
 */

#ifndef PROJ_HYDRATE_SIMPLEXDG_HH_
#define PROJ_HYDRATE_SIMPLEXDG_HH_

template <class GV>
void proj_Hydrate_SimplexDG(const GV &gv, // GridView
							double dt,	  // (sec) Time-step size
							double t_END, // (sec) End time
							double t_OP	  // (sec) Interval at which output must be plotted,
)
{

	//	INCLUDE EXTERNAL CLASSES
	IncludeClasses paramclass;
	//	PostProcess postprocess;

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	const int dim = GV::dimension;
	/*Non-dimensionalize time prams*/
	dt *= 1. / paramclass.characteristicValue.t_c;
	t_END *= 1. / paramclass.characteristicValue.t_c;
	t_OP *= 1. / paramclass.characteristicValue.t_c;
	Real dt_min = paramclass.timeStepControl.mindt / paramclass.characteristicValue.t_c;
	Real dt_max = paramclass.timeStepControl.maxdt / paramclass.characteristicValue.t_c;

	Real time = 0.0;
	Real dtstart = dt;
	Real time_op = time;

	int maxAllowableIterations = 6;
	int minAllowableIterations = 3;

	const int degree_S = 1;
	const int degree_P = 1;
	const int degree_T = 1;

	//	GFS
#ifdef PARALLEL
	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON0;
#else
	typedef Dune::PDELab::NoConstraints CON0;	// pure Neumann: no constraints
#endif									
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE0;	// default block size: 1
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_P, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_P; 
	FEM_P fem_P;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_S, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_S;// basis function
	FEM_S fem_S;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_T, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_T; 
	FEM_T fem_T;
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P, CON0, VBE0> GFS_P; // gfs
	GFS_P gfs_P(gv, fem_P);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_S, CON0, VBE0> GFS_S; // gfs
	GFS_S gfs_S(gv, fem_S);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_T, CON0, VBE0> GFS_T; // gfs
	GFS_T gfs_T(gv, fem_T);

	//	COMPOSITE GFS
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE; //  block size -> numOfFlowPVs

	// gfs for composite system Pg,Pc,Sw,Sh,T
	typedef Dune::PDELab::CompositeGridFunctionSpace<VBE,
													 Dune::PDELab::EntityBlockedOrderingTag,
													 GFS_P,
													 GFS_P,
													 GFS_S,
													 GFS_S,
													 GFS_T>
		GFS;
	GFS gfs(gfs_P, gfs_P, gfs_S, gfs_S, gfs_T);
	//PowerGridFunctionSpace< GFS_P,
	// 											  Indices::numOfFlowPVs,
	// 											  VBE,
	// 											  Dune::PDELab::EntityBlockedOrderingTag > GFS;
	// GFS gfs(gfs_P);
	typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
	CC cc;
	//cc.clear();
	gfs.update(); // initializing the gfs
	std::cout << "degrees of freedom: " << gfs.globalSize() << std::endl;

	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
	//	primary variables: Pg, Pc, Sw, Sh, T
	typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Pg>> SUBGFS_Pg; //
	SUBGFS_Pg subgfs_Pg(gfs);
	typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Pc>> SUBGFS_Pc; //
	SUBGFS_Pc subgfs_Pc(gfs);
	typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Sw>> SUBGFS_Sw; //
	SUBGFS_Sw subgfs_Sw(gfs);
	typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Sh>> SUBGFS_Sh; //
	SUBGFS_Sh subgfs_Sh(gfs);
	typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_T>> SUBGFS_T; //
	SUBGFS_T subgfs_T(gfs);

	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS, double>;
	U uold(gfs, 0.0);
	
	U unew(gfs, 0.0);


	//	MAKE FUNCTION FOR INITIAL VALUES
	typedef Pg_Initial<GV, Real> Pg_InitialType;
	Pg_InitialType Pg_initial(gv);

	typedef Pc_Initial<GV, Real> Pc_InitialType;
	Pc_InitialType Pc_initial(gv);
	typedef Sw_Initial<GV, Real> Sw_InitialType;
	Sw_InitialType Sw_initial(gv);
	typedef Sh_Initial<GV, Real> Sh_InitialType;
	Sh_InitialType Sh_initial(gv);
	typedef T_Initial<GV, Real> T_InitialType;
	T_InitialType T_initial(gv);
	typedef Dune::PDELab::CompositeGridFunction<Pg_InitialType,
												Pc_InitialType,
												Sw_InitialType,
												Sh_InitialType,
												T_InitialType>
		InitialType;
	InitialType initial(Pg_initial, Pc_initial, Sw_initial, Sh_initial, T_initial);

	Dune::PDELab::interpolate(initial, gfs, uold); // Initialize the solution at t=0 (uold) with the given initial values

	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	ConvectionDiffusionDGMethod::Type method_g = ConvectionDiffusionDGMethod::NIPG;
	ConvectionDiffusionDGMethod::Type method_w = ConvectionDiffusionDGMethod::NIPG;
	ConvectionDiffusionDGMethod::Type method_T = ConvectionDiffusionDGMethod::NIPG;
	double alpha_g = 10.;
	double alpha_w = 10.;
	double alpha_s = 10.;
	double alpha_T = 10000.;

	typedef FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T> LOP; // spatial part
	LOP lop(gv, &unew, gfs, &time, &dt, 6, method_g, method_w, method_T, alpha_g, alpha_w, alpha_s, alpha_T);

	typedef FLOW_TimeOperator TLOP; // temporal part
	TLOP tlop;

	typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
	MBE mbe(20);

	typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GOLOP;
	GOLOP goLOP(gfs, cc, gfs, cc, lop, mbe);

	// How well did we estimate the number of entries per matrix row?
	// => print Jacobian pattern statistics
	typename GOLOP::Traits::Jacobian jac(goLOP);
	//std::cout << jac.patternStatistics() << std::endl;
	//if (helper.rank() == 1)
	std::cout << " LOP DONE ! " << std::endl;

	typedef Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, Real, Real, Real, CC, CC> GOTLOP;
	GOTLOP goTLOP(gfs, cc, gfs, cc, tlop, mbe);

	//std::cout << " TLOP DONE ! " << std::endl;

	typedef Dune::PDELab::OneStepGridOperator<GOLOP, GOTLOP> IGO;
	IGO igo(goLOP, goTLOP);
	//if (helper.rank() == 1)
	std::cout << " IGO DONE ! " << std::endl;

	// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL

	//make vector consistent NEW IN PARALLEL
	Dune::PDELab::ISTL::ParallelHelper<GFS> helper(gfs);
	helper.maskForeignDOFs(uold);
	Dune::PDELab::AddDataHandle<GFS, U> adddh(gfs, uold);
	if (gfs.gridView().comm().size() > 1)
		gfs.gridView().communicate(adddh, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

	//typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS, CC> LS;// works
	//LS ls(gfs, cc, 1000, 2);

	//typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; //should be checked
	//LS ls(gfs, 100, 1, true, true);

	//typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	//LS ls(gfs, 1000, 1, true, true);

	typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS; //works
	LS ls(gfs,5000,1,true,true);

	//typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<GFS, CC> LS; //works
	//LS ls(gfs, cc);

	// typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; // should be checked
	// int verbose = 0;
	// if (gfs.gridView().comm().rank() == 0)
	// 	verbose = 1;
	// LS ls(gfs, 100, verbose);
	std::cout << " PARALLEL LS DONE ! " << std::endl;
#else
	//typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	//LS ls(1000, true);

	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;
	std::cout << " LS DONE ! " << std::endl;
#endif

	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	typedef Dune::PDELab::NewtonMethod<IGO, LS> PDESOLVER;
	PDESOLVER pdesolver(igo, ls);
	// 	select control parameters for non-linear PDE-solver
	typedef Dune::PDELab::LineSearchHackbuschReusken<PDESOLVER> lineSearchStrategy;
	lineSearchStrategy linesearchstrategy(pdesolver);
	pdesolver.setVerbosityLevel(2);
	pdesolver.setReduction(1e-6);
	pdesolver.setMinLinearReduction(1e-6);
	pdesolver.setAbsoluteLimit(1e-6);
	//if (helper.rank() == 1)
	std::cout << " PDESOLVER DONE ! " << std::endl;

	// SELECT TIME-STEPPER
	Dune::PDELab::ImplicitEulerParameter<Real> method1;
	Dune::PDELab::OneStepThetaParameter<Real> method2(0.0); //Crank-Nicholson -> 0.5, Implicit Euler -> 1.0, Explicit Euler -> 0.0
	Dune::PDELab::Alexander2Parameter<Real> method3;
	Dune::PDELab::Alexander3Parameter<Real> method4;
	Dune::PDELab::FractionalStepParameter<Real> method5;
	Dune::PDELab::HeunParameter<Real> method6;
	Dune::PDELab::Shu3Parameter<Real> method7;
	Dune::PDELab::RK4Parameter<Real> method8;

	Dune::PDELab::TimeSteppingParameterInterface<Real> *pmethod = &method1;

	Dune::PDELab::OneStepMethod<Real, IGO, PDESOLVER, U, U> osm(*pmethod, igo, pdesolver);
	osm.setVerbosityLevel(2);
	//if (helper.rank() == 1)
	std::cout << " OSM DONE ! " << std::endl;

	//	GRAPHICS FOR INITIAL GUESS
	// primary variables
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pg, U> DGF_Pg;
	DGF_Pg dgf_Pg(subgfs_Pg, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pc, U> DGF_Pc;
	DGF_Pc dgf_Pc(subgfs_Pc, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sw, U> DGF_Sw;
	DGF_Sw dgf_Sw(subgfs_Sw, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sh, U> DGF_Sh;
	DGF_Sh dgf_Sh(subgfs_Sh, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_T, U> DGF_T;
	DGF_T dgf_T(subgfs_T, uold);

	// secondary variables
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_Pw, U_PP> DGFPP_Pw;
	// DGFPP_Pw dgfpp_Pw(subgfspp_Pw, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_Sg, U_PP> DGFPP_Sg;
	// DGFPP_Sg dgfpp_Sg(subgfspp_Sg, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_Pc, U_PP> DGFPP_Pc;
	// DGFPP_Pc dgfpp_Pc(subgfspp_Pc, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_K, U_PP> DGFPP_K;
	// DGFPP_K dgfpp_K(subgfspp_K, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_krw, U_PP> DGFPP_krw;
	// DGFPP_krw dgfpp_krw(subgfspp_krw, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_krg, U_PP> DGFPP_krg;
	// DGFPP_krg dgfpp_krg(subgfspp_krg, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_zCH4, U_PP> DGFPP_zCH4;
	// DGFPP_zCH4 dgfpp_zCH4(subgfspp_zCH4, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_XCH4, U_PP> DGFPP_XCH4;
	// DGFPP_XCH4 dgfpp_XCH4(subgfspp_XCH4, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_XH2O, U_PP> DGFPP_XH2O;
	// DGFPP_XH2O dgfpp_XH2O(subgfspp_XH2O, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_YCH4, U_PP> DGFPP_YCH4;
	// DGFPP_YCH4 dgfpp_YCH4(subgfspp_YCH4, upp);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBGFSPP_YH2O, U_PP> DGFPP_YH2O;
	// DGFPP_YH2O dgfpp_YH2O(subgfspp_YH2O, upp);

	//	VTK
	auto pathName = paramclass.problemSpecs.getPathName();
	auto fileName = paramclass.problemSpecs.getFileName();
	const std::string str = "";
	Dune::PDELab::FilenameHelper fn(pathName + fileName);

	int subsampling = 0;
	typedef Dune::SubsamplingVTKWriter<GV> VTKWRITER;
	VTKWRITER vtkwriter(gv, subsampling, (int)0);
	typedef Dune::VTKSequenceWriter<GV> VTKSEQUENCEWRITER;
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), fileName, pathName, "");

	// add data field for all components of the space to the VTK writer
	// primary variables
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pg>>(dgf_Pg, "Pg"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pc>>(dgf_Pc, "Pc"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_Sh, "Sh"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sw>>(dgf_Sw, "Sw"));

	vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
	vtkSequenceWriter.clear();

	//	INITIALIZE
	unew = uold;
	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;

	bool exceptionCaught = false;

	int newton_iterations = 0;

	//	BEGIN TIME LOOP
	while (time < t_END - 1e-8)
	{
		std::cout << "_____________________________________________________" << std::endl;
		//std::cout<< " current opcount = " << opcount - 1 << std::endl;

		try
		{
			std::cout << "****************************" << std::endl;
			std::cout << "  CALLING osm.apply() !" << std::endl;
			std::cout << "****************************" << std::endl;
			osm.apply(time, dt, uold, unew);
			newton_iterations = osm.getPDESolver().result().iterations;
			std::cout << "****************************" << std::endl;
			std::cout << " osm.apply() DONE !" << std::endl;
			std::cout << "****************************" << std::endl;
			//slp_proj.apply(uproj);

			exceptionCaught = false;
		}
		catch (Dune::Exception &e)
		{
			exceptionCaught = true;
			if (dt > 1e-6)
			{
				std::cout << "Catched Error, Dune reported error: " << e << std::endl;

				unew = uold;

				dt *= 0.5;
				continue;
			}
			else
			{
				std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				exit(0);
			}
		}

		std::cout << "DONE" << std::endl;
		std::cout << "_____________________________________________________" << std::endl;

		// GRAPHICS FOR NEW OUTPUT
		// primary variables
		typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pg, U> DGF_Pg;
		DGF_Pg dgf_Pg(subgfs_Pg, unew);
		typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pc, U> DGF_Pc;
		DGF_Pc dgf_Pc(subgfs_Pc, unew);
		typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sw, U> DGF_Sw;
		DGF_Sw dgf_Sw(subgfs_Sw, unew);
		typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sh, U> DGF_Sh;
		DGF_Sh dgf_Sh(subgfs_Sh, unew);
		typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_T, U> DGF_T;
		DGF_T dgf_T(subgfs_T, unew);

		/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
		if ((time + dt > t_OP * opcount - dt_min) and (time + dt < t_OP * opcount + 1.e-6))
		{
			//vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			// primary variables
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pg>>(dgf_Pg, "Pg"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pc>>(dgf_Pc, "Pc"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_Sh, "Sh"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sw>>(dgf_Sw, "Sw"));

			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			vtkSequenceWriter.clear();
			std::cout << " ******************************************************************* " << std::endl;
			std::cout << " OUTPUT WRITTEN " << opcount << " ----processor: " << std::endl;
			std::cout << " ******************************************************************* " << std::endl;

			timecount = time;
			opcount = opcount + 1;
		}

		//		PREPARE FOR NEXT TIME INTEGRATION
		//		1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		uold = unew;
		//		2. ADVANCE TIME:
		time += dt;
		std::cout << " " << std::endl;
		std::cout << " time = " << time * paramclass.characteristicValue.t_c;

		if (paramclass.timeStepControl.adaptiveTimeStepControl)
		{
			if (newton_iterations > maxAllowableIterations)
			{
				dt = std::max(dt * 0.9, dt_min);
			}
			else if (newton_iterations <= minAllowableIterations)
			{
				dt = std::min(dt * 2, dt_max);
			}
		}
		else
		{
			dt = dtstart;
		}

		std::cout << " , time+dt = " << (time + dt) * paramclass.characteristicValue.t_c
				  << " , opTime = " << t_OP * opcount * paramclass.characteristicValue.t_c;

		if (time + dt > t_OP * opcount - 1.e-6)
		{
			dtLast = dt;
			dt = t_OP * opcount - time;

			std::cout << " , because timeNext > opNext , dt set to : " << dt * paramclass.characteristicValue.t_c << std::endl;
			dtFlag = 0;
		}
		dtFlag += 1;

		if (opcount > 1 and dtFlag == 2)
		{
			dt = std::max(dt, dtLast * 0.6);
		}
		std::cout << " , dt  : " << dt * paramclass.characteristicValue.t_c << std::endl;
		std::cout << " " << std::endl;

		std::cout << " READY FOR NEXT ITERATION. " << std::endl;
	}
};

#endif /* PROJ_HYDRATE_SIMPLEXDG_HH_ */
