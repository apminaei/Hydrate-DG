/*
 * driver.hh
 *
 */

#ifndef PROJ_HYDRATE_SIMPLEXDG_HH_
#define PROJ_HYDRATE_SIMPLEXDG_HH_

template <class GV, class PTree>
void driver(const GV &gv, // GridView
			const PTree& ptree, Dune::MPIHelper& helper)
{

	//	INCLUDE EXTERNAL CLASSES
	typedef Properties<GV,PTree> Properties;

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	const int dim = GV::dimension;
	Real time = 0.0;
	Real dt = 0.0;
	
	Properties property(gv,ptree);
	/* Non-dimensionalize time prams */
	dt = ptree.get("time.dt_initial",(double)1.)/ property.characteristicValue.t_c;
	Real t_END = ptree.get("time.time_end",(double)86400.) / property.characteristicValue.t_c;
	Real t_OP = ptree.get("output.time_interval",(double)1.) / property.characteristicValue.t_c;
	Real dt_min = ptree.get("adaptive_time_control.dt_min",(double)0.001) / property.characteristicValue.t_c;
	Real dt_max = ptree.get("adaptive_time_control.dt_max",(double)1.) / property.characteristicValue.t_c;
	bool adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	
	Real dtstart = dt;
	Real time_op = time;
	//std::cout << " time = " << time << " t_end = " << t_END <<  " dt_min = " << dt_min << " dt = " << dt << std::endl;
	//exit(0);
	int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)12);
	int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);

	const int degree_S = 1;
	const int degree_P = 1;
	const int degree_T = 1;
	const int degree_X = 1;
	const int degree_Y = 1;
	//	GFS
#ifdef PARALLEL
	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON0;
#else
	typedef Dune::PDELab::ConformingDirichletConstraints CON0;	// pure Neumann: no constraints
#endif									
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE0;	// default block size: 1
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_P, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_P;// 
	FEM_P fem_P;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_S, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_S;// basis function
	FEM_S fem_S;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_T, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_T; 
	FEM_T fem_T;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_X, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_X; 
	FEM_X fem_x;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Y, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_Y; 
	FEM_Y fem_y;
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P, CON0, VBE0> GFS_P; // gfs
	GFS_P gfs_P(gv, fem_P);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_S, CON0, VBE0> GFS_S; // gfs
	GFS_S gfs_S(gv, fem_S);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_T, CON0, VBE0> GFS_T; // gfs
	GFS_T gfs_T(gv, fem_T);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_X, CON0, VBE0> GFS_X; // gfs
	GFS_X gfs_x(gv, fem_x);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Y, CON0, VBE0> GFS_Y; // gfs
	GFS_Y gfs_y(gv, fem_y);

	//	COMPOSITE GFS
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE; //  block size -> numOfPVs

	// gfs for composite system Pw,Pc,Sg,Sh,T,XCH4,YH2O
	typedef Dune::PDELab::CompositeGridFunctionSpace<VBE, Dune::PDELab::EntityBlockedOrderingTag,
													 GFS_P, GFS_P, GFS_S, GFS_S, GFS_T, GFS_X,GFS_Y,GFS_X>
		GFS;
	GFS gfs(gfs_P, gfs_P, gfs_S, gfs_S, gfs_T, gfs_x, gfs_y,gfs_x);

	// typedef PowerGridFunctionSpace< GFS_P, Indices::numOfPVs, VBE, Dune::PDELab::LexicographicOrderingTag > GFS;
	// GFS gfs(gfs_P);
	
	
	// BCTypeParam0<GV> u0_bctype(gv);
    // BCTypeParam1<GV> u1_bctype(gv);
    // typedef Dune::PDELab::CompositeConstraintsParameters<BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>,BCTypeParam0<GV>> U_BCTypeParam;
	
	// U_BCTypeParam u_bctype( u0_bctype, u0_bctype,u0_bctype,u0_bctype,u0_bctype,u0_bctype,u0_bctype,u0_bctype );
    typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    CC cc;                       
	gfs.update(); // initializing the gfs
    //Dune::PDELab::constraints( u_bctype, gfs, cc);  // to artificial boundaries
    std::cout << "constrained dofs = " << cc.size() << " of " << gfs.globalSize() << std::endl;

	

	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
	using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pw>>;
    using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;
    SUBGFS_Pw    subgfs_Pw(gfs);
	using PathPc = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pc>>;
    using SUBGFS_Pc = Dune::PDELab::GridFunctionSubSpace<GFS,PathPc>;
    SUBGFS_Pc    subgfs_Pc(gfs);
	using PathSg = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sg>>;
    using SUBGFS_Sg = Dune::PDELab::GridFunctionSubSpace<GFS,PathSg>;
    SUBGFS_Sg    subgfs_Sg(gfs);
	using PathSh = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sh>>;
    using SUBGFS_Sh = Dune::PDELab::GridFunctionSubSpace<GFS,PathSh>;
    SUBGFS_Sh    subgfs_Sh(gfs);
	using PathT = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_T>>;
    using SUBGFS_T = Dune::PDELab::GridFunctionSubSpace<GFS,PathT>;
    SUBGFS_T    subgfs_T(gfs);
	using PathXCH4 = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_XCH4>>;
    using SUBGFS_XCH4 = Dune::PDELab::GridFunctionSubSpace<GFS,PathXCH4>;
    SUBGFS_XCH4    subgfs_XCH4(gfs);
	using PathYH2O = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_YH2O>>;
    using SUBGFS_YH2O = Dune::PDELab::GridFunctionSubSpace<GFS,PathYH2O>;
    SUBGFS_YH2O    subgfs_YH2O(gfs);
	using Pathc = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_C>>;
    using SUBGFS_XC = Dune::PDELab::GridFunctionSubSpace<GFS,Pathc>;
    SUBGFS_XC    subgfs_XC(gfs);

	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS, double>;
	U uold(gfs, 0.0);
	
	U unew(gfs, 0.0);


	//	MAKE FUNCTION FOR INITIAL VALUES   Which must be nondim 
	typedef Pw_Initial<GV,Properties,Real> Pw_InitialType;
	Pw_InitialType Pw_initial(gv,property); /* ndim */

	typedef Pc_Initial<GV,Properties,Real> Pc_InitialType;
	Pc_InitialType Pc_initial(gv,property); /* ndim */
	typedef Sg_Initial<GV,Properties,Real> Sg_InitialType;
	Sg_InitialType Sg_initial(gv,property);
	typedef Sh_Initial<GV,Properties,Real> Sh_InitialType;
	Sh_InitialType Sh_initial(gv,property);
	typedef T_Initial<GV,Properties,Real> T_InitialType;
	T_InitialType T_initial(gv,property); /* ndim */
	typedef XCH4_Initial<GV,Properties,Real> XCH4_InitialType;
	XCH4_InitialType XCH4_initial(gv,property);
	typedef YH2O_Initial<GV,Properties,Real> YH2O_InitialType;
	YH2O_InitialType YH2O_initial(gv,property);
	typedef XC_Initial<GV,Properties,Real> XC_InitialType;
	XC_InitialType XC_initial(gv,property);
	typedef Dune::PDELab::CompositeGridFunction<Pw_InitialType,
												Pc_InitialType,
												Sg_InitialType,
												Sh_InitialType,
												T_InitialType,
												XCH4_InitialType,
												YH2O_InitialType, XC_InitialType>
		InitialType;
	InitialType initial(Pw_initial, Pc_initial, Sg_initial, Sh_initial, T_initial, XCH4_initial, YH2O_initial, XC_initial);

	Dune::PDELab::interpolate(initial, gfs, uold); // Initialize the solution at t=0 (uold) with the given initial values

	//	BOUNDARY CONDITIONS
	
	// const Dune::FieldVector<double, dim> xtest(0.1);// xtest[0] = 0.15;xtest[1]=0.25;
	// std::cout << " bct[indices.PVId_Pw ] = " <<  bc.type( xtest, 9.e5, dt)[Indices::PVId_Pw] << std::endl;
	// std::cout << " bcvalue[indices.PVId_Pw ] = " <<  bc.value( xtest, 9.e5, dt)[Indices::PVId_Pw] << std::endl;
	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	ConvectionDiffusionDGMethod::Type method_g = ConvectionDiffusionDGMethod::IIPG;
	ConvectionDiffusionDGMethod::Type method_w = ConvectionDiffusionDGMethod::IIPG;
	ConvectionDiffusionDGMethod::Type method_T = ConvectionDiffusionDGMethod::IIPG;
	ConvectionDiffusionDGMethod::Type method_x = ConvectionDiffusionDGMethod::IIPG;
	ConvectionDiffusionDGMethod::Type method_y = ConvectionDiffusionDGMethod::IIPG;
	double alpha_g = 0.e1;
	double alpha_w = 0.e1;
	double alpha_s = 0.e1;
	double alpha_T = 0.e1;
	double alpha_x = 0.e1;
	double alpha_y = 0.e1;

	typedef LocalOperator<GV, Properties, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y> LOP; // spatial part
	//time *= (1./property.characteristicValue.t_c);
	LOP lop(gv, property, &unew, gfs, &time, &dt, 6, method_g, method_w, method_T, method_x, method_y, alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y);

	typedef TimeOperator<GV, Properties> TLOP; // temporal part
	TLOP tlop(gv, property);

	typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
	MBE mbe(100);

	typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GOLOP;
	GOLOP goLOP(gfs, cc, gfs, cc, lop, mbe);

	// How well did we estimate the number of entries per matrix row?
	// => print Jacobian pattern statistics
	typename GOLOP::Traits::Jacobian jac(goLOP);
	std::cout << " LOP DONE ! " << std::endl;

	typedef Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, Real, Real, Real, CC, CC> GOTLOP;
	GOTLOP goTLOP(gfs, cc, gfs, cc, tlop, mbe);

	typedef Dune::PDELab::OneStepGridOperator<GOLOP, GOTLOP> IGO;
	IGO igo(goLOP, goTLOP);
	std::cout << " IGO DONE ! " << std::endl;

	// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL

	//make vector consistent NEW IN PARALLEL
	Dune::PDELab::ISTL::ParallelHelper<GFS> phelper(gfs);
	phelper.maskForeignDOFs(uold);
	Dune::PDELab::AddDataHandle<GFS, U> adddh(gfs, uold);
	if (gfs.gridView().comm().size() > 1)
		gfs.gridView().communicate(adddh, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

	typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS, CC> LS;// works
	LS ls(gfs, cc, 100, 2);

	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	// LS ls(gfs, 100, 1, true, true);

	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS; //works
	// LS ls(gfs,5000,1,true,true);

	//typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<GFS, CC> LS; //works
	//LS ls(gfs, cc);

	// typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; // should be checked
	// int verbose = 0;
	// if (gfs.gridView().comm().rank() == 0)
	// 	verbose = 1;
	// LS ls(gfs, 100, verbose);
	std::cout << " PARALLEL LS DONE ! " << std::endl;
#else
	// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	// LS ls(1000, true);

	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;
	std::cout << " LS DONE ! " << std::endl;
#endif

    //    SELECT SOLVER FOR NON-LINEAR PROBLEM
    using PDESOLVER = Dune::PDELab::Newton< IGO, LS, U >;
    PDESOLVER pdesolver( igo, ls );
    //     select control parameters for non-linear PDE-solver
	pdesolver.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    //pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	pdesolver.setReassembleThreshold(0.0);
    pdesolver.setVerbosityLevel(2);
    pdesolver.setReduction(1e-5);
    pdesolver.setMinLinearReduction(1e-5);
	pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    pdesolver.setForceIteration(true);
	pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-5)); 

	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	// typedef Dune::PDELab::NewtonMethod<IGO, LS> PDESOLVER;
	// PDESOLVER pdesolver(igo, ls);
	
	// 	select control parameters for non-linear PDE-solver
	//typedef Dune::PDELab::LineSearchNone<PDESOLVER> lineSearchStrategy;
	
	//typedef Dune::PDELab::LineSearchHackbuschReusken<PDESOLVER> lineSearchStrategy;
	//lineSearchStrategy linesearchstrategy(pdesolver);
	
	//pdesolver.setParameters(ptree.sub("newton"));
	// pdesolver.setVerbosityLevel(2);
	// pdesolver.setReduction(1e-3);
	// pdesolver.setMinLinearReduction(1e-6);
	// pdesolver.setAbsoluteLimit(1e-6);
	
	std::cout << " PDESOLVER DONE ! " << std::endl;

	// SELECT TIME-STEPPER
	Dune::PDELab::ImplicitEulerParameter<Real> method1;
	Dune::PDELab::OneStepThetaParameter<Real> method2(0.5); //Crank-Nicholson -> 0.5, Implicit Euler -> 1.0, Explicit Euler -> 0.0
	Dune::PDELab::Alexander2Parameter<Real> method3;
	Dune::PDELab::Alexander3Parameter<Real> method4;
	Dune::PDELab::FractionalStepParameter<Real> method5;
	Dune::PDELab::HeunParameter<Real> method6;
	Dune::PDELab::Shu3Parameter<Real> method7;
	Dune::PDELab::RK4Parameter<Real> method8; // didnot work

	Dune::PDELab::TimeSteppingParameterInterface<Real> *pmethod = &method1;

	Dune::PDELab::OneStepMethod<Real, IGO, PDESOLVER, U, U> osm(*pmethod, igo, pdesolver);
	osm.setVerbosityLevel(2);
	std::cout << " OSM DONE ! " << std::endl;

	//	GRAPHICS FOR INITIAL GUESS
	// primary variables
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pw, U> DGF_Pw;
	DGF_Pw dgf_Pw(subgfs_Pw, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pc, U> DGF_Pc;
	DGF_Pc dgf_Pc(subgfs_Pc, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sg, U> DGF_Sg;
	DGF_Sg dgf_Sg(subgfs_Sg, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sh, U> DGF_Sh;
	DGF_Sh dgf_Sh(subgfs_Sh, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_T, U> DGF_T;
	DGF_T dgf_T(subgfs_T, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_XCH4, U> DGF_XCH4;
	DGF_XCH4 dgf_XCH4(subgfs_XCH4, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_YH2O, U> DGF_YH2O;
	DGF_YH2O dgf_YH2O(subgfs_YH2O, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_XC, U> DGF_XC;
	DGF_XC dgf_XC(subgfs_XC, uold);


	//	VTK
	auto pathName = "/home/amir/dune-master/Hydrate-DG/dune/Hydrate-DG/akerbp2D_pockmark/outputs/";
	auto fileName = ptree.get("output.file_name",(std::string)"test");
	const std::string str = "";
	//Dune::PDELab::FilenameHelper fn(pathName + fileName);

	int subsampling = 1;
	Dune::RefinementIntervals RefInt(subsampling);

	using VTKWRITER = Dune::SubsamplingVTKWriter<GV> ;
	VTKWRITER vtkwriter(gv, RefInt, false, Dune::VTK::Precision::float32);
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV> ;
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), fileName, pathName, "");

	// add data field for all components of the space to the VTK writer
	// primary variables
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw>>(dgf_Pw, "Pw"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pc>>(dgf_Pc, "Pc"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_Sh, "Sh"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sg>>(dgf_Sg, "Sg"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XCH4>>(dgf_XCH4, "XCH4"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_YH2O>>(dgf_YH2O, "YH2O"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XC>>(dgf_XC, "XC"));

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
	while ( time < t_END - 1e-8)
	{
		std::cout << "_____________________________________________________" << std::endl;
		
		try
		{
			std::cout << "****************************" << std::endl;
			std::cout << "  CALLING osm.apply() !" << std::endl;
			std::cout << "****************************" << std::endl;
			osm.apply(time, dt, uold, unew);
			
			newton_iterations = osm.getPDESolver().result().iterations;
			std::cout << "****************************" << std::endl;
		
			exceptionCaught = true;
		}
		catch(Dune::Exception &e)
		{
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
		DGF_Pw dgf_Pw(subgfs_Pw, unew);
		DGF_Pc dgf_Pc(subgfs_Pc, unew);
		DGF_Sg dgf_Sg(subgfs_Sg, unew);
		DGF_Sh dgf_Sh(subgfs_Sh, unew);
		DGF_T dgf_T(subgfs_T, unew);
		DGF_XCH4 dgf_XCH4(subgfs_XCH4, unew);
		DGF_YH2O dgf_YH2O(subgfs_YH2O, unew);
		DGF_XC dgf_XC(subgfs_XC, unew);


		/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
		if ((time + dt > t_OP * opcount - dt_min) and (time + dt < t_OP * opcount + 1.e-6))
		{
			//vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			// primary variables
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw>>(dgf_Pw, "Pw"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pc>>(dgf_Pc, "Pc"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_Sh, "Sh"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sg>>(dgf_Sg, "Sg"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XCH4>>(dgf_XCH4, "XCH4"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_YH2O>>(dgf_YH2O, "YH2O"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XC>>(dgf_XC, "XC"));

			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			vtkSequenceWriter.clear();
			std::cout << " ******************************************************************* " << std::endl;
			std::cout << " OUTPUT WRITTEN " << opcount << " ----processor: " << helper.rank() << std::endl;
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
		std::cout << " time = " << time ;

		if (adaptive_time_control)
		{
			if (newton_iterations > maxAllowableIterations)
			{
				dt = std::max(dt*0.6 , dt_min);
			}
			else if (newton_iterations <= minAllowableIterations)
			{
				dt = std::min(dt * 1.1, dt_max);
			}
		}
		else
		{
			dt = dtstart;
		}

		std::cout << " , time+dt = " << (time + dt) 
				  << " , opTime = " << t_OP * opcount ;

		if (time + dt > t_OP * opcount - 1.e-7)
		{
			dtLast = dt;
			dt = t_OP * opcount - time;

			std::cout << " , because timeNext > opNext , dt set to : " << dt  << std::endl;
			dtFlag = 0;
		}
		dtFlag += 1;

		if (opcount > 1 and dtFlag == 2)
		{
			dt = std::max(dt, dtLast * 0.5);
		}
		std::cout << " , dt  : " << dt << std::endl;
		std::cout << " " << std::endl;

		std::cout << " READY FOR NEXT ITERATION. " << std::endl;
	}
};

#endif /* PROJ_HYDRATE_SIMPLEXDG_HH_ */
