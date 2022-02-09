/*
 * driver.hh
 *
 */

#ifndef HYDRATE_DG_HH_
#define HYDRATE_DG_HH_


template <class GV, class PTree>
void driver(const GV &gv, // GridView
			const PTree& ptree, //ParameterTree 
			Dune::MPIHelper& helper //MPI-Helper 
		   )
{

	//	INCLUDE EXTERNAL CLASSES


	typedef Properties<GV, PTree> Properties;

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;
	const int dim = GV::dimension;
	Real time = 0.0;
	Real dt = 0.0;
	
	Properties property(gv,ptree);

	

	/* Non-dimensionalize time prams */
	/*in input file dt, dt_initprob, t_END, t_END_initprob, dt_min, dt_max, t_OP are specified in kilo-annum*/
	//INITIAL PROBLEM time and dt
	auto Xc_t = property.characteristicValue.t_c;
	auto t_day_sec = 12.*30.*24*60.*60.;
	// double dt_initprob  = ptree.get("initial_problem.dt_initial",(double)0.0001);
	// dt_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// dt_initprob *= 1./Xc_t; /*ndim*/
	// double t_END_initprob  = ptree.get("initial_problem.time_end",(double)1);
	// t_END_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// t_END_initprob 	 *= 1./Xc_t; /*ndim*/
	//MAIN PROBLEM time and dt
	dt  = ptree.get("time.dt_initial",(double)0.001);
	dt *= t_day_sec/Xc_t; /*convert to seconds then nondim*/
	double t_END  = ptree.get("time.time_end",(double)300);
	t_END *= t_day_sec/Xc_t; /*convert to seconds*/
	//t_END += t_END_initprob; //Total simulation time
	// output time interval
	double t_OP   = ptree.get("output.time_interval",(double)1000);
	t_OP *= t_day_sec/Xc_t; /*convert to seconds*/
	//adaptive time control
	bool adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6);
	dt_min *= t_day_sec/Xc_t; /*convert to seconds*/
	double dt_max = ptree.get("adaptive_time_control.dt_max",(double)10.);
	dt_max *= t_day_sec/Xc_t; /*convert to seconds*/


	
	// Real T_end = property.parameter.T_END();
	// double time_fraction = t_END / T_end;
	
	
	Real dtstart = dt;
	// Real time_op = time;
	Real clock_time_elapsed = 0.;

	int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)10);
	int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);
	int max_linear_iteration = ptree.get("newton.MaxLinearIteration",(int)10);
	const int degree_Sg = 1;
	const int degree_P = 1;
	const int degree_Sh = 1;
	const int degree_T = 1;
	const int degree_X = 1;
	const int degree_C = 1;
	const int degree_Y = 1;

	const int degree_PP = 1;

	// const int degree_pen = 1;

	//	GFS
#ifdef PARALLEL
	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON0;
#else
	typedef Dune::PDELab::ConformingDirichletConstraints CON0;	// pure Neumann: no constraints
#endif		
// #ifdef PARALLEL
// 	using CON0 = Dune::PDELab::P0ParallelConstraints;
// #else
// 	using CON0 = Dune::PDELab::NoConstraints;
// #endif							
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE0;	// default block size: 1
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_P,dim,Dune::GeometryType::simplex > FEM_P;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_Sg,dim,Dune::GeometryType::simplex > FEM_Sg;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_T,dim,Dune::GeometryType::simplex > FEM_T;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_X,dim,Dune::GeometryType::simplex > FEM_X;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_Y,dim,Dune::GeometryType::simplex > FEM_Y;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_Sh,dim,Dune::GeometryType::simplex > FEM_Sh;
	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_X,dim,Dune::GeometryType::simplex > FEM_C;

	// typedef OPBLocalFiniteElementMap<Coord,Real,degree_PP,dim,Dune::GeometryType::simplex > FEM_PP;

	auto const xc_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	auto const xch4_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	auto const yh2o_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	auto const pw_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;// ptree.get("basis.pw",);
	auto const sg_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	auto const sh_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	auto const t_basis_polynomial = Dune::PDELab::QkDGBasisPolynomial::lagrange;
	
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_P, dim, pw_basis_polynomial > FEM_P;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Sg, dim, sg_basis_polynomial > FEM_Sg;
	
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_T, dim, t_basis_polynomial > FEM_T; 
	
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_X, dim, xch4_basis_polynomial > FEM_X; 

	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Y, dim, yh2o_basis_polynomial > FEM_Y;

	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Sh, dim, sh_basis_polynomial > FEM_Sh; 

	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_C, dim, xc_basis_polynomial > FEM_C;  

	// typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_PP, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_PP;//l2orthonormal 

	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_P> FEM_P;
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_Sg> FEM_Sg;
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_Sh> FEM_Sh;
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_T> FEM_T; 
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_X> FEM_X; 
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_Y> FEM_Y;
	
	// typedef Dune::PDELab::PkLocalFiniteElementMap<GV, Coord, Real, degree_C> FEM_C; 

	typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, Real, degree_PP> FEM_PP;

	FEM_X fem_x;//(gv);
	FEM_T fem_T;//(gv);
	FEM_Sg fem_Sg;//(gv);
	FEM_P fem_P;//(gv);
	FEM_Y fem_y;//(gv);
	FEM_C fem_c;//(gv);
	FEM_Sh fem_Sh;//(gv);
	
	
	// if (fem_x.polynomial()==Dune::PDELab::QkDGBasisPolynomial::lagrange){
	// 	std::cout << "fem_x.polynomial()" << std::endl;
	// 	exit(0);
	// }
	FEM_PP fem_PP(gv);

	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_P, CON0, VBE0> GFS_P; // gfs
	GFS_P gfs_P(gv, fem_P);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Sg, CON0, VBE0> GFS_Sg; // gfs
	GFS_Sg gfs_Sg(gv, fem_Sg);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Sh, CON0, VBE0> GFS_Sh; // gfs
	GFS_Sh gfs_Sh(gv, fem_Sh);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_T, CON0, VBE0> GFS_T; // gfs
	GFS_T gfs_T(gv, fem_T);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_X, CON0, VBE0> GFS_X; // gfs
	GFS_X gfs_x(gv, fem_x);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Y, CON0, VBE0> GFS_Y; // gfs
	GFS_Y gfs_y(gv, fem_y);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_C, CON0, VBE0> GFS_C; // gfs
	GFS_C gfs_c(gv, fem_c);

	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_PP, CON0, VBE0> GFS_CC; // gfs
	GFS_CC gfs_cc(gv, fem_PP);
	// typedef Dune::PDELab::GridFunctionSpace<GV, FEM_pen, CON0, VBE0> GFS_pen; // gfs
	// GFS_pen gfs_pen(gv, fem_pen);
	//	COMPOSITE GFS
	using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;//  block size -> numOfPVs

	// gfs for composite system Pw,Sg,Sh,T,XCH4,YH2O,XC
	typedef Dune::PDELab::CompositeGridFunctionSpace<VBE, 
													 Dune::PDELab::EntityBlockedOrderingTag,
													 GFS_P, GFS_Sg, GFS_Sh, GFS_T, GFS_X,GFS_Y,GFS_C> 
													 GFS;
	// typedef Dune::PDELab::PowerGridFunctionSpace< GFS_P,
	// 												  Indices::numOfPVs,
	// 												  VBE0,
	// 												  Dune::PDELab::LexicographicOrderingTag > GFS;
	GFS gfs(gfs_P,  gfs_Sg, gfs_Sh, gfs_T, gfs_x, gfs_y,gfs_c);
	
    typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    CC cc;
	cc.clear();
	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS, double>;
	U uold(gfs, 0.0);
	U unew(gfs, 0.0);


	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
	using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pw>>;
    using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;
    SUBGFS_Pw    subgfs_Pw(gfs);
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

	//  EVALUATION FUNCTIONS FOR PRIMARY VARIABLES
	Dune::PDELab::Evaluation<SUBGFS_Pw	,U> evaluation_Pw(	 subgfs_Pw	,&unew );
	Dune::PDELab::Evaluation<SUBGFS_Sg	,U> evaluation_Sg(	 subgfs_Sg	,&unew );
	Dune::PDELab::Evaluation<SUBGFS_Sh	,U> evaluation_Sh(	 subgfs_Sh	,&unew );
	Dune::PDELab::Evaluation<SUBGFS_T,U> evaluation_T( subgfs_T,&unew );
	Dune::PDELab::Evaluation<SUBGFS_YH2O,U> evaluation_YH2O( subgfs_YH2O,&unew );
	Dune::PDELab::Evaluation<SUBGFS_XCH4,U> evaluation_XCH4( subgfs_XCH4,&unew );
	Dune::PDELab::Evaluation<SUBGFS_XC,U> evaluation_XC( subgfs_XC,&unew );

	Dune::PDELab::Evaluation<SUBGFS_Pw	,U> evaluation_Pw_old(	 subgfs_Pw	,&uold );
	Dune::PDELab::Evaluation<SUBGFS_Sg	,U> evaluation_Sg_old(	 subgfs_Sg	,&uold );
	Dune::PDELab::Evaluation<SUBGFS_Sh	,U> evaluation_Sh_old(	 subgfs_Sh	,&uold );
	Dune::PDELab::Evaluation<SUBGFS_T,U> evaluation_T_old( subgfs_T,&unew );
	Dune::PDELab::Evaluation<SUBGFS_YH2O,U> evaluation_YH2O_old( subgfs_YH2O,&uold );
	Dune::PDELab::Evaluation<SUBGFS_XCH4,U> evaluation_XCH4_old( subgfs_XCH4,&uold );
	Dune::PDELab::Evaluation<SUBGFS_XC,U> evaluation_XC_old( subgfs_XC,&uold );

	// using Pathp = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_pen>>;
    // using SUBGFS_pen = Dune::PDELab::GridFunctionSubSpace<GFS,Pathp>;
    // SUBGFS_pen    subgfs_pen(gfs);

	

	//	MAKE FUNCTION FOR INITIAL VALUES   Which must be nondim 
	typedef Pw_Initial<GV,Properties,Real> Pw_InitialType;
	Pw_InitialType Pw_initial(gv,property); /* ndim */
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
	// auto pen_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 1.;};
	// auto pen = Dune::PDELab::makeGridFunctionFromCallable(gv, pen_Initiallamda);
	typedef Dune::PDELab::CompositeGridFunction<Pw_InitialType,
												Sg_InitialType,
												Sh_InitialType,
												T_InitialType,
												XCH4_InitialType,
												YH2O_InitialType, XC_InitialType> InitialType;
	InitialType initial(Pw_initial, Sg_initial, Sh_initial, T_initial, XCH4_initial, YH2O_initial, XC_initial);
	Dune::PDELab::interpolate(initial, gfs, uold); 
    // Initialize the solution at t=0 (uold) with the given initial values
	
	// ProblemInitialConditions<GV,Properties> Icvalue;
	// Icvalue icvalue(gv, property);
	// auto Pw_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){
	// 	return (1.5e7 + 1030.21 * 9.81 * (0.-x[1])*property.characteristicValue.x_c)/property.characteristicValue.P_c;};
	// auto Pw = Dune::PDELab::makeGridFunctionFromCallable(gv, Pw_Initiallamda);
	// auto Sg_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 0.;};
	// auto Sg = Dune::PDELab::makeGridFunctionFromCallable(gv, Sg_Initiallamda);
	// auto Sh_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){double Sh = 0.0 ;
	// 	double GHSZ_width = -3.20 - (-4.);
		
	// 	if( x[1]<= -3.20 && x[1]>=-4.){
	// 		Sh = 1.2 * (x[1]-(-4.))/GHSZ_width * (x[1]-(-3.2))/(-GHSZ_width);//* (rand()%2) + 0.001;//
	// 	}
	// 	return Sh;
	// };
	// auto Sh = Dune::PDELab::makeGridFunctionFromCallable(gv, Sh_Initiallamda);
	// auto T_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){
	// 	return (4.65 + 273.15 - 0.035 * x[1]*property.characteristicValue.x_c)/property.characteristicValue.T_c;};
	// auto T = Dune::PDELab::makeGridFunctionFromCallable(gv, T_Initiallamda);
	// auto XCH4_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 0.;};
	// auto XCH4 = Dune::PDELab::makeGridFunctionFromCallable(gv, XCH4_Initiallamda);
	// // auto YH2O_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 6.e-5;};
	// // auto YH2O = Dune::PDELab::makeGridFunctionFromCallable(gv, YH2O_Initiallamda);
	// auto XC_Initiallamda = [&](const Dune::FieldVector<double,dim>& x){return 0.035 * (property.gas.MolarMass()/property.salt.MolarMass());};
	// auto XC = Dune::PDELab::makeGridFunctionFromCallable(gv, XC_Initiallamda);
	// using InitialTypelamda = Dune::PDELab::CompositeGridFunction<decltype(Pw), decltype(Sg), decltype(Sh),
	// 									decltype(T), decltype(XCH4),
	// 									YH2O_InitialType, decltype(XC)> ;
	// InitialTypelamda initialamda(Pw, Sg, Sh, T, XCH4, YH2O_initial, XC);
	// Dune::PDELab::interpolate(initialamda, gfs, uold); // Initialize the solution at t=0 (uold) with the given initial values
	

	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	double method_g = 0.;//ConvectionDiffusionDGMethod::SIPG;
	double method_w = 0.;//ConvectionDiffusionDGMethod::SIPG;
	double method_T = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double method_x = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double method_y = 0.;//ConvectionDiffusionDGMethod::NIPG;
	double alpha_g = ptree.get("penalty_coeff.Sg",(double)1.e3); // 1.e3; //
	double alpha_w = ptree.get("penalty_coeff.Pw",(double)1.e3); // 1.e1; //
	double alpha_s = ptree.get("penalty_coeff.Sh",(double)1.e3); // 1.e1; //
	double alpha_T = ptree.get("penalty_coeff.T",(double)1.e3); // 1.e1; //
	double alpha_x = ptree.get("penalty_coeff.XC",(double)1.e3); // 1.e1; //
	double alpha_y = ptree.get("penalty_coeff.YH2O",(double)1.e3); // 1.e1; //
	double intorder= ptree.get("quadrature.order",(int)6);

	typedef ProblemBoundaryConditions<GV,Properties> BoundaryConditions ;
	BoundaryConditions bc( gv,property ) ;

	typedef LocalOperator<GV, Properties, BoundaryConditions, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C> LOP; // spatial part
	LOP lop(gv, property, bc, gfs, &time, &dt, intorder, method_g, method_w, method_T, method_x, method_y, alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y);

	// typedef LocalOperator<GV, Properties, U, GFS,
	// 	 FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C,
	// 	Dune::PDELab::Evaluation<SUBGFS_Pw	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_Sg	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_Sh	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_T	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_XCH4 ,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_YH2O ,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_XC	,U>> LOP; // spatial part
	// LOP lop(gv, property, &unew, gfs,
	// 						&evaluation_Pw_old,
	// 						&evaluation_Sg_old,
	// 						&evaluation_Sh_old,
	// 						&evaluation_T_old,
	// 						&evaluation_XCH4_old,
	// 						&evaluation_YH2O_old,
	// 						&evaluation_XC_old,
	// 						 &time, &dt, intorder, method_g, method_w, method_T, method_x, method_y, alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y);

	
	typedef TimeOperator<GV, Properties> TLOP; // temporal part
	TLOP tlop(gv, property, intorder);
	// typedef TimeOperator<GV, Properties,
	// 	Dune::PDELab::Evaluation<SUBGFS_Pw	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_Sg	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_Sh	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_T	,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_XCH4 ,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_YH2O ,U>,
	// 	Dune::PDELab::Evaluation<SUBGFS_XC	,U>> TLOP; // temporal part
	// TLOP tlop(gv, property,
	// 						&evaluation_Pw_old,
	// 						&evaluation_Sg_old,
	// 						&evaluation_Sh_old,
	// 						&evaluation_T_old,
	// 						&evaluation_XCH4_old,
	// 						&evaluation_YH2O_old,
	// 						&evaluation_XC_old, intorder);

	typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
	MBE mbe(150);

	typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GOLOP;
	GOLOP goLOP(gfs, cc, gfs, cc, lop, mbe);

	// How well did we estimate the number of entries per matrix row?
	// => print Jacobian pattern statistics0
	typename GOLOP::Traits::Jacobian jac(goLOP);
	std::cout << " LOP DONE ! " << std::endl;

	typedef Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, Real, Real, Real, CC, CC> GOTLOP;
	GOTLOP goTLOP(gfs, cc, gfs, cc, tlop, mbe);

	typedef Dune::PDELab::OneStepGridOperator<GOLOP, GOTLOP> IGO;
	IGO igo(goLOP, goTLOP);
	std::cout << " IGO DONE ! " << std::endl;

	// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS, CC> LS; //doesn't work
	// LS ls(gfs, cc, 100, 2);

	// typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFS,CC> LS;
	// LS ls(gfs,cc,100,2,max_linear_iteration);

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SSORk<GFS,CC> LS;
	// LS ls(gfs, cc, max_linear_iteration, 10,1);

	typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	LS ls(gfs, max_linear_iteration, 1, false, true);
	/* 	NOTES:
		LINEAR SOLVER STATISTICS
		res.iterations = i;
		res.reduction = def/def0;
		res.conv_rate  = pow(res.reduction,1.0/i);
		res.elapsed = watch.elapsed();
		// final print
		if (_verbose>0)
		{
			std::cout << "=== rate=" << res.conv_rate
					<< ", T=" << res.elapsed
					<< ", TIT=" << res.elapsed/i
					<< ", IT=" << i << std::endl;
		}
	*/

	

	// using LS = Dune::PDELab::ISTLBackend_OVLP_AMG_4_DG< IGO, CC, GFS, CC,
	// Dune::PDELab::CG2DGProlongation,
	// Dune::SeqSSOR,
	// Dune::BiCGSTABSolver>;
	// LS ls(igo, cc, gfs, cc, max_linear_iteration,2,true,true);
	// // set parameters for AMG in CG-subspace
	// Dune::Amg::Parameters params = ls.parameters();
	// params.setCoarsenTarget(10000);
	// params.setMaxLevel(20);
	// params.setProlongationDampingFactor(0.5);
	// params.setNoPreSmoothSteps(6);
	// params.setNoPostSmoothSteps(6);
	// params.setGamma(1);
	// params.setAdditive(false);
	// ls.setParameters(params);
	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS;//doesn't work
	// LS ls(gfs,max_linear_iteration,2,true,true);

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<GFS, CC> LS; //doesn't work
	// LS ls(gfs, cc);

	// typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; // should be checked
	// int verbose = 0;
	// // if (gfs.gridView().comm().rank() == 0)
	// // 	verbose = 1;
	// LS ls(gfs, 100, verbose);

	auto param = ls.parameters();
	// param.setMaxLevel(3); // max number of coarsening levels
	param.setCoarsenTarget(100000); // max DoF at coarsest level
	ls.setParameters(param);

	std::cout << " PARALLEL LS DONE ! " << std::endl;

#else
	// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;//doesn't work
	// LS ls(max_linear_iteration, true);

	// using LS = Dune::PDELab::ISTLBackend_SEQ_AMG_4_DG<IGO,GFS,Dune::PDELab::CG2DGProlongation,Dune::SeqSSOR,Dune::BiCGSTABSolver>; //works
	// LS ls(igo,gfs,max_linear_iteration,2,true,true);
	// // // // set parameters for AMG in CG-subspace
	// Dune::Amg::Parameters params = ls.parameters();
	// params.setCoarsenTarget(10000);
	// params.setMaxLevel(20);
	// params.setProlongationDampingFactor(1.8);
	// params.setNoPreSmoothSteps(3);
	// params.setNoPostSmoothSteps(3);
	// params.setGamma(1);
	// params.setAdditive(false);
	// ls.setParameters(params);

	
	
	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;
	// std::cout << " LS DONE ! " << std::endl;
#endif

    //    SELECT SOLVER FOR NON-LINEAR PROBLEM
    using PDESOLVER = Dune::PDELab::Newton< IGO, LS, U >;
    PDESOLVER pdesolver( igo, ls );
    // select control parameters for non-linear PDE-solver
	pdesolver.setLineSearchStrategy(ptree.get("newton.LineSearchStrategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    //pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	pdesolver.setReassembleThreshold(0.0);
    pdesolver.setVerbosityLevel(2);
    pdesolver.setReduction(ptree.get("newton.Reduction",(double)1e-5));
    pdesolver.setMinLinearReduction(ptree.get("newton.MinLinearReduction",(double)1.e-9));
	pdesolver.setMaxIterations(ptree.get("newton.MaxIterations",(int)15));
    pdesolver.setForceIteration(ptree.get("newton.ForceIteration",(bool)true));
	pdesolver.setAbsoluteLimit(ptree.get("newton.AbsoluteLimit",(double)1.e-4)); 

	//TODO: CHECK NEW NEWTON PARAMS
	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	// typedef Dune::PDELab::NewtonMethod<IGO, LS> PDESOLVER;
	// PDESOLVER pdesolver(igo, ls);

	// 	select control parameters for non-linear PDE-solver
	// typedef Dune::PDELab::LineSearchNone<PDESOLVER> lineSearchStrategy;
	// typedef Dune::PDELab::LineSearchHackbuschReusken<PDESOLVER> lineSearchStrategy;
	// lineSearchStrategy linesearchstrategy(pdesolver);
	// pdesolver.setParameters(ptree.sub("newton"));
	//linesearchstrategy.setParameters(ptree.sub("newton.line_search"));
	// pdesolver.setVerbosityLevel(ptree.get("newton.VerbosityLevel",(int)3));
	// pdesolver.setReduction(ptree.get("newton.Reduction",(double)1e-5));
    // // pdesolver.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver.setAbsoluteLimit(ptree.get("newton.AbsoluteLimit",(double)1.e-4)); 
	// auto reduction = ptree.get("newton.Reduction",(double)1e-5);
	// auto absolutlimit = ptree.get("newton.AbsoluteLimit",(double)1e-5);
	
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
	unew = uold;
	

	

	// STEP 1:
		// L2-PROJECTION OF P0-PRESSURES ONTO Q1-SPACE
		const int degree = 1; //polynomial degree
		typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, double, degree> FEMPROJ;
		FEMPROJ femproj(gv);
		typedef Dune::PDELab::GridFunctionSpace<GV, FEMPROJ, CON0, VBE0> GFSPROJ0;
		GFSPROJ0 gfsproj0(gv, femproj);
#ifdef PARALLEL
		using VBEPROJ = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFSPROJ0,
													  Indices::numOfPVs,
													  VBEPROJ,
													  Dune::PDELab::EntityBlockedOrderingTag > GFSPROJ;
		GFSPROJ gfsproj(gfsproj0);
#else
		typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none > VBEPROJ;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFSPROJ0,
													  Indices::numOfPVs,
													  VBEPROJ,
													  Dune::PDELab::LexicographicOrderingTag > GFSPROJ;
		GFSPROJ gfsproj(gfsproj0);
#endif
		typedef typename GFSPROJ::template ConstraintsContainer<double>::Type CCPROJ;
		CCPROJ ccproj;
		ccproj.clear();

		using SUBPROJ_Pw = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pw>> >;
		SUBPROJ_Pw subproj_Pw(gfsproj);
		using SUBPROJ_Sg = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sg>> >;
		SUBPROJ_Sg subproj_Sg(gfsproj);
		using SUBPROJ_Sh = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sh>> >;
		SUBPROJ_Sh subproj_Sh(gfsproj);
		using SUBPROJ_T = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_T>> >;
		SUBPROJ_T subproj_T(gfsproj);
		using SUBPROJ_XC = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_C>> >;
		SUBPROJ_XC subproj_XC(gfsproj);
		using SUBPROJ_XCH4 = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_XCH4>> >;
		SUBPROJ_XCH4 subproj_XCH4(gfsproj);
		using SUBPROJ_YH2O = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_YH2O>> >;
		SUBPROJ_YH2O subproj_YH2O(gfsproj);

		using UPROJ = Dune::PDELab::Backend::Vector<GFSPROJ,double>;
		UPROJ u_proj(gfsproj); // u_proj = [Pg_q1, Pw_q1]^T

		Dune::PDELab::Evaluation<SUBPROJ_Sg	,UPROJ> evaluation_projSg(	 subproj_Sg	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_Pw	,UPROJ> evaluation_projPw(	 subproj_Pw	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_Sh	,UPROJ> evaluation_projSh(	 subproj_Sh	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_T	,UPROJ> evaluation_projT(	 subproj_T	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_XC	,UPROJ> evaluation_projXC(	 subproj_XC	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_XCH4	,UPROJ> evaluation_projXCH4(	 subproj_XCH4	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_YH2O	,UPROJ> evaluation_projYH2O(	 subproj_YH2O	,&u_proj );

		typedef LocalOperatorPROJ< GV,Properties,GFSPROJ,UPROJ,
		Dune::PDELab::Evaluation<SUBGFS_Pw	,U>,
		Dune::PDELab::Evaluation<SUBGFS_Sg	,U>,
		Dune::PDELab::Evaluation<SUBGFS_Sh	,U>,
		Dune::PDELab::Evaluation<SUBGFS_T	,U>,
		Dune::PDELab::Evaluation<SUBGFS_XCH4	,U>,
		Dune::PDELab::Evaluation<SUBGFS_YH2O	,U>,
		Dune::PDELab::Evaluation<SUBGFS_XC	,U>
		 > LOPPROJ;	// spatial part
		LOPPROJ lopproj(gv, property, gfsproj, &u_proj, &evaluation_Pw,
												 &evaluation_Sg,
												 &evaluation_Sh,
												 &evaluation_T,
												 &evaluation_XCH4,
												 &evaluation_YH2O,
												 &evaluation_XC);

		MBE mbeproj(10);
		typedef Dune::PDELab::GridOperator< GFSPROJ, GFSPROJ, LOPPROJ, MBE, double, double, double, CCPROJ, CCPROJ> GOPROJ;
		GOPROJ goproj(gfsproj, ccproj, gfsproj, ccproj, lopproj, mbeproj);

#ifdef PARALLEL
	#ifdef ALUGRID
			typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFSPROJ,CCPROJ> LSPROJ;
			LSPROJ lsproj(gfsproj,ccproj,100,1,200);
	#else 
			typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOPROJ> LSPROJ;
			LSPROJ lsproj(gfsproj,100,3,false,true);
	#endif
#else
		typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LSPROJ;
		LSPROJ lsproj;
#endif

		typedef Dune::PDELab::StationaryLinearProblemSolver<GOPROJ,LSPROJ,UPROJ> SLPPROJ;
		SLPPROJ slpproj(goproj,lsproj,u_proj,1e-10);
		// slpproj.apply(); // get initial values of u_proj


	

		
	//  POST-PROCESS
		typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed> VBE_PP;
		typedef Dune::PDELab::CompositeGridFunctionSpace<VBE, 
													 Dune::PDELab::EntityBlockedOrderingTag,
													 GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC,
													 GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC,
													 GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC,
													 GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC, GFS_CC,
													 GFS_CC, GFS_CC, GFS_CC, GFS_CC>  GFS_PP;
		// typedef Dune::PDELab::PowerGridFunctionSpace< GFS_CC,
		// 											  Indices::numOfSVs,
		// 											  VBE0,
		// 											  Dune::PDELab::LexicographicOrderingTag > GFS_PP;
		GFS_PP gfs_pp(	gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc,
						gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc,
						gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc,
						gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc, gfs_cc,
						 gfs_cc, gfs_cc, gfs_cc, gfs_cc);
		
		typedef typename GFS_PP::template ConstraintsContainer<Real>::Type CC_PP;
		CC_PP cc_pp;
		cc_pp.clear();

		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pg>	  > > SUBPP_Pg;
		SUBPP_Pg 	subpp_Pg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pw>	  > > SUBPP_Pw;
		SUBPP_Pw 	subpp_Pw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pc>	  > > SUBPP_Pc;
		SUBPP_Pc 	subpp_Pc(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sg>	  > > SUBPP_Sg;
		SUBPP_Sg 	subpp_Sg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sh>	  > > SUBPP_Sh;
		SUBPP_Sh 	subpp_Sh(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sw>	  > > SUBPP_Sw;
		SUBPP_Sw 	subpp_Sw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_T>  > > SUBPP_T;
		SUBPP_T 	subpp_T(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XCH4>  > > SUBPP_XCH4;
		SUBPP_XCH4 	subpp_XCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XC>  > > SUBPP_XC;
		SUBPP_XC 	subpp_XC(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XH2O>  > > SUBPP_XH2O;
		SUBPP_XH2O 	subpp_XH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YCH4 > > > SUBPP_YCH4;
		SUBPP_YCH4 	subpp_YCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YH2O  >> > SUBPP_YH2O;
		SUBPP_YH2O 	subpp_YH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_rhow  >> > SUBPP_rhow;
		SUBPP_rhow 	subpp_rhow(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_rhog  > >> SUBPP_rhog;
		SUBPP_rhog 	subpp_rhog(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_K	  > >> SUBPP_K;
		SUBPP_K 	subpp_K(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_krw	  > >> SUBPP_krw;
		SUBPP_krw 	subpp_krw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_krg	  > > >SUBPP_krg;
		SUBPP_krg 	subpp_krg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_muw>	  > > SUBPP_muw;
		SUBPP_muw 	subpp_muw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_mug>	  > > SUBPP_mug;
		SUBPP_mug 	subpp_mug(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_zCH4>  > > SUBPP_zCH4;
		SUBPP_zCH4 	subpp_zCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_por	>  > > SUBPP_por;
		SUBPP_por 	subpp_por(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_DH2O  >> > SUBPP_DH2O;
		SUBPP_DH2O 	subpp_DH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_DCH4  >> > SUBPP_DCH4;
		SUBPP_DCH4 	subpp_DCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pwsat > >> SUBPP_Pwsat;
		SUBPP_Pwsat 	subpp_Pwsat(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_HCH4  > >> SUBPP_HCH4;
		SUBPP_HCH4 	subpp_HCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_tau	  > >> SUBPP_tau;
		SUBPP_tau 	subpp_tau(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Peq	  > >> SUBPP_Peq;
		SUBPP_Peq 	subpp_Peq(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_HS > >> SUBPP_HS;
		SUBPP_HS 	subpp_HS(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Vwx > >> SUBPP_Vwx;
		SUBPP_Vwx 	subpp_Vwx(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Vwy > >> SUBPP_Vwy;
		SUBPP_Vwy 	subpp_Vwy(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Vgx > >> SUBPP_Vgx;
		SUBPP_Vgx 	subpp_Vgx(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Vgy > >> SUBPP_Vgy;
		SUBPP_Vgy 	subpp_Vgy(gfs_pp);

		
		using U_PP = Dune::PDELab::Backend::Vector<GFS_PP,double>;
		U_PP u_pp(gfs_pp, 0.0);
		

		PostProcess<GV, Properties, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C,
					Dune::PDELab::Evaluation<SUBGFS_Pw	,U> ,
					Dune::PDELab::Evaluation<SUBGFS_Sg	,U> ,
					Dune::PDELab::Evaluation<SUBGFS_Sh	,U> ,
					Dune::PDELab::Evaluation<SUBGFS_T,U> ,
					Dune::PDELab::Evaluation<SUBGFS_XCH4,U> ,
					Dune::PDELab::Evaluation<SUBGFS_YH2O,U> ,
					Dune::PDELab::Evaluation<SUBGFS_XC,U>,
					 GFS_PP, U_PP > postprocess( gv, property, fem_P, fem_Sg, fem_Sh, fem_T, fem_x, fem_y, fem_c, 
												 &evaluation_Pw,
												 &evaluation_Sg,
												 &evaluation_Sh,
												 &evaluation_T,
												 &evaluation_XCH4,
												 &evaluation_YH2O,
												 &evaluation_XC,
												 gfs_pp, &u_pp );
		postprocess.evaluate();

		std::cout << " Post process init DONE ! " << std::endl;
	//	GRAPHICS FOR INITIAL GUESS
	// primary variables
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_Pw, UPROJ> DGF_Pw;
	// DGF_Pw dgf_Pw(subproj_Pw, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_Sg, UPROJ> DGF_Sg;
	// DGF_Sg dgf_Sg(subproj_Sg, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_Sh, UPROJ> DGF_Sh;
	// DGF_Sh dgf_Sh(subproj_Sh, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_T, UPROJ> DGF_T;
	// DGF_T dgf_T(subproj_T, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_XCH4, UPROJ> DGF_XCH4;
	// DGF_XCH4 dgf_XCH4(subproj_XCH4, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_YH2O, UPROJ> DGF_YH2O;
	// DGF_YH2O dgf_YH2O(subproj_YH2O, u_proj);
	// typedef Dune::PDELab::DiscreteGridFunction<SUBPROJ_XC, UPROJ> DGF_XC;
	// DGF_XC dgf_XC(subproj_XC, u_proj);

	// secondary variables
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pg, U_PP > DGF_Pg;
		DGF_Pg dgf_Pg( subpp_Pg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pw, U_PP > DGF_Pw;
		DGF_Pw dgf_Pw( subpp_Pw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sg, U_PP > DGF_Sg;
		DGF_Sg dgf_Sg( subpp_Sg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sh, U_PP > DGF_Sh;
		DGF_Sh dgf_Sh( subpp_Sh, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_T, U_PP > DGF_T;
		DGF_T dgf_T( subpp_T, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XCH4, U_PP > DGF_XCH4;
		DGF_XCH4 dgf_XCH4( subpp_XCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XC, U_PP > DGF_XC;
		DGF_XC dgf_XC( subpp_XC, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YH2O, U_PP > DGF_YH2O;
		DGF_YH2O dgf_YH2O( subpp_YH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pc, U_PP > DGF_Pc;
		DGF_Pc dgf_Pc( subpp_Pc, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sw, U_PP > DGF_Sw;
		DGF_Sw dgf_Sw( subpp_Sw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XH2O, U_PP > DGF_XH2O;
		DGF_XH2O dgf_XH2O( subpp_XH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YCH4, U_PP > DGF_YCH4;
		DGF_YCH4 dgf_YCH4( subpp_YCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_rhow, U_PP > DGF_rhow;
		DGF_rhow dgf_rhow( subpp_rhow, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_rhog, U_PP > DGF_rhog;
		DGF_rhog dgf_rhog( subpp_rhog, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_K, U_PP > DGF_K;
		DGF_K dgf_K( subpp_K, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_krw, U_PP > DGF_krw;
		DGF_krw dgf_krw( subpp_krw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_krg, U_PP > DGF_krg;
		DGF_krg dgf_krg( subpp_krg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_muw, U_PP > DGF_muw;
		DGF_muw dgf_muw( subpp_muw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_mug, U_PP > DGF_mug;
		DGF_mug dgf_mug( subpp_mug, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_zCH4, U_PP > DGF_zCH4;
		DGF_zCH4 dgf_zCH4( subpp_zCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_por, U_PP > DGF_por;
		DGF_por dgf_por( subpp_por, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_DH2O, U_PP > DGF_DH2O;
		DGF_DH2O dgf_DH2O( subpp_DH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_DCH4, U_PP > DGF_DCH4;
		DGF_DCH4 dgf_DCH4( subpp_DCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pwsat, U_PP > DGF_Pwsat;
		DGF_Pwsat dgf_Pwsat( subpp_Pwsat, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_HCH4, U_PP > DGF_HCH4;
		DGF_HCH4 dgf_HCH4( subpp_HCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_tau, U_PP > DGF_tau;
		DGF_tau dgf_tau( subpp_tau, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Peq, U_PP > DGF_Peq;
		DGF_Peq dgf_Peq( subpp_Peq, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_HS, U_PP > DGF_HS;
		DGF_HS dgf_HS( subpp_HS, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Vwx, U_PP > DGF_Vwx;
		DGF_Vwx dgf_Vwx( subpp_Vwx, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Vwy, U_PP > DGF_Vwy;
		DGF_Vwy dgf_Vwy( subpp_Vwy, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Vgx, U_PP > DGF_Vgx;
		DGF_Vgx dgf_Vgx( subpp_Vgx, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Vgy, U_PP > DGF_Vgy;
		DGF_Vgy dgf_Vgy( subpp_Vgy, u_pp );
	//*****************************


	//*******************************	

	 
	//	VTK
	std::string fileName = ptree.get("output.file_name",(std::string)"test");
	std::string pathName = ptree.get("output.path_name",(std::string)"test");
	pathName += "outputs/";
	pathName += fileName ;
	std::string fileNameDefects = fileName;
	std::string pathNameDefects = pathName+"/"+fileName;
	std::time_t now = std::time(0);
	struct tm *tstruct;
	char buf [80];
	std::time(&now);
	tstruct = std::localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", tstruct );
	std::string timeStr(buf); 
	pathNameDefects.append("_"+timeStr);
    std::string jacPath = pathNameDefects;
		jacPath +="/";
		jacPath +=fileName;

	if(helper.rank()==0){
		
		std::filesystem::create_directory(pathNameDefects);
	}
	// pathName += "/" ;
	// const std::string str = "";
	//Dune::PDELab::FilenameHelper fn(pathName + fileName);

	int subsampling = 1;
	Dune::RefinementIntervals RefInt(subsampling);
	using VTKWRITER = Dune::VTKWriter<GV>;
	VTKWRITER vtkwriter(gv,Dune::VTK::nonconforming);
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
//	VTKSEQUENCEWRITER vtkSequenceWriter( gv,fileName,pathName,"",Dune::VTK::nonconforming);
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");
	

	// using VTKWRITER = Dune::SubsamplingVTKWriter<GV> ;
	// VTKWRITER vtkwriter(gv, RefInt, false, Dune::VTK::Precision::float32); // vtk nonconforming
	// using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV> ;
	// VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), fileName, pathName, "");

	// add data field for all components of the space to the VTK writer
	// primary variables
	// Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,uold,Dune::PDELab::vtk::DefaultFunctionNameGenerator("u"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw>>(dgf_Pw, "Pw"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_Sh, "Sh"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sg>>(dgf_Sg, "Sg"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_T, "T"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XCH4>>(dgf_XCH4, "XCH4"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_YH2O>>(dgf_YH2O, "YH2O"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XC>>(dgf_XC, "XC"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pg>>(dgf_Pg, "Pg"));
	vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pc>>(dgf_Pc, "Pc"));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sw    > >( dgf_Sw    , "Sw"   ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_XH2O  > >( dgf_XH2O  , "XH2O" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_YCH4  > >( dgf_YCH4  , "YCH4" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_rhow  > >( dgf_rhow  , "rhow" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_rhog  > >( dgf_rhog  , "rhog" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_K     > >( dgf_K     , "K"    ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_krw   > >( dgf_krw   , "krw"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_krg   > >( dgf_krg   , "krg"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_muw   > >( dgf_muw   , "muw"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_mug   > >( dgf_mug   , "mug"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_zCH4  > >( dgf_zCH4  , "z"    ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_por   > >( dgf_por   , "porosity"));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_DH2O  > >( dgf_DH2O  , "DH2O" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_DCH4  > >( dgf_DCH4  , "DCH4" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pwsat > >( dgf_Pwsat , "Pwsat"));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_HCH4  > >( dgf_HCH4  , "HCH4" ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_tau   > >( dgf_tau   , "tau"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Peq   > >( dgf_Peq   , "Peq"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_HS > >( dgf_HS   , "HS"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Vwx   > >( dgf_Vwx   , "Vwx"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Vwy   > >( dgf_Vwy   , "Vwy"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Vgx   > >( dgf_Vgx   , "Vgx"  ));
	vtkSequenceWriter.addVertexData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Vgy   > >( dgf_Vgy   , "Vgy"  ));



	vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
	// vtkSequenceWriter.clear();
	
	std::string dgmethod_g = std::__cxx11::to_string(method_g);
	std::string dgmethod_w = std::__cxx11::to_string(method_w);
	std::string dgmethod_T = std::__cxx11::to_string(method_T);
	std::string dgmethod_x = std::__cxx11::to_string(method_x);
	std::string dgmethod_y = std::__cxx11::to_string(method_y);
	double dissCoeff = property.parameter.HydrateDissociationRateConstant();
	double formCoeff = property.parameter.HydrateFormationRateConstant();
	if(helper.rank()==0){
		std::string parameters_file = pathNameDefects;
		parameters_file +="/";
		parameters_file +=fileName;
		parameters_file +="_parameters";
		parameters_file += ".txt";
		property.ReportParameters( 	parameters_file,
									dgmethod_g, dgmethod_w, dgmethod_T, dgmethod_x, dgmethod_y,
									alpha_g, alpha_w, alpha_s, alpha_T, alpha_x, alpha_y, dissCoeff, formCoeff,
									property.characteristicValue.permeability_c,
									property.characteristicValue.X_gravity,
									property.characteristicValue.X_source_mass,
									property.characteristicValue.X_source_heat,
									property.characteristicValue.dispersivity_c,
									property.characteristicValue.specificheat_c);
	}

	//	INITIALIZE
	unew = uold;
	int opcount = 1;
	// double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;

	bool exceptionCaught = false;

	int newton_iterations = 0;
	double newton_first_defect = 0.;
	double newton_defect = 0.;
	//Dune::BlockVector<Dune::FieldVector<Real, 1>> newton_defects;
	auto time_in_year = time * Xc_t / t_day_sec;
	auto dt_in_year = dt * Xc_t / t_day_sec;
	
	//	BEGIN TIME LOOP
	while ( time < (t_END - 1e-9/Xc_t))
	{
		if( exceptionCaught==false ){
				dt = std::max(dt,dt_min);
		}

		if(helper.rank()==0){
			std::cout<< "_____________________________________________________" <<std::endl;
			std::cout<< " current opcount = " << opcount - 1 << std::endl;
		}

		clock_t start = clock();
		try{
			if(helper.rank()==0){
			std::cout<<"****************************" << std::endl;
			std::cout<<"  CALLING osm.apply() !"	  << std::endl;
			std::cout<<"****************************" << std::endl;
			}
			time_in_year = time * Xc_t / t_day_sec;
			dt_in_year = dt * Xc_t / t_day_sec;
			osm.apply( time, dt, uold, unew );
			
			newton_iterations = osm.getPDESolver().result().iterations;
			newton_first_defect = osm.getPDESolver().result().first_defect;
			newton_defect = osm.getPDESolver().result().defect;
			
            // auto jacobian = osm.getPDESolver().getJacobian();
			// if(helper.rank()==0 &&  (opcount%100==0)){
			// Dune::writeMatrixToMatlab ( Dune::PDELab::Backend::native(jacobian), jacPath+"jacobian");
			// Dune::writeVectorToMatlab(Dune::PDELab::Backend::native(unew),jacPath+"solution");
			// }
			auto newton_defects = osm.getPDESolver().result().defects;
			auto u_norm_two = unew.two_norm();
			auto u_norm_one = unew.one_norm();
			auto u_norm_infinity = unew.infinity_norm();
			if(helper.rank()==0 &&  (newton_iterations>1)){//((time+dt )/(t_OP * opcount) > (1.-1.e-6)) and ((time+dt ) / (t_OP * opcount)< (1. + 1.e-6))
				std::string s_OP = std::__cxx11::to_string(time_in_year);
				std::string parameters_file = pathNameDefects;
				parameters_file +="/";
				parameters_file +=fileNameDefects;
				parameters_file +="_";
				parameters_file +=s_OP;
				parameters_file += ".txt";
				property.ReportNewton( parameters_file,
							time_in_year /*s*/,
							dt_in_year /*s*/,
							newton_iterations, newton_defects,
							u_norm_two, u_norm_one,  u_norm_infinity);
			}
			//std::cout << "  ================////////=============== " <<  newton_defects[newton_iterations] << std::endl;
			/*
				std::cout << "  Newton iteration " << std::setw(2)(Sets the field width to be used on output operations) << this->res_.iterations
                           << ".  New defect: "
                           << std::setw(12) << std::setprecision(4) << std::scientific (Sets the floatfield format flag for the str stream to scientific)
                           << this->res_.defect
                           << ".  Reduction (this): "
                           << std::setw(12) << std::setprecision(4) << std::scientific
                           << this->res_.defect/this->prev_defect_
                           << ".  Reduction (total): "
                           << std::setw(12) << std::setprecision(4) << std::scientific
                           << this->res_.reduction << std::endl;
			*/

			exceptionCaught = false;

		}catch ( Dune::Exception &e ) {
			exceptionCaught = true;
			if( dt > (2. * dt_min) ){

				if(helper.rank()==0){
					std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}

				unew = uold;
				
				dt = 0.5 * dt;
				
				newton_iterations = 0;
				dtLast = dt;
					continue;
			}
			else
			{
				if(helper.rank()==0){
					std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
				}
				exit(0);
			}
		}
		clock_t end = clock();
		double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
		clock_time_elapsed += clock_time_this_step;
		if(helper.rank()==0 ){//&& newton_defect >  absolutlimit && newton_defect > (newton_first_defect * reduction)
			std::string statistics_file = pathNameDefects;
			statistics_file +="/";
			statistics_file +=fileName;
			statistics_file +="_statistics";
			statistics_file += ".txt";
					
			property.ReportStatistics( 	statistics_file,
											time_in_year,
											dt_in_year,
											newton_iterations,
											newton_first_defect,
											newton_defect,
											clock_time_elapsed );
		}

		if(helper.rank()==0){
			std::cout<<"DONE"<<std::endl;
			std::cout<<"_____________________________________________________"<<std::endl;
		}

		/*********************************************************************************************
		 * OUTPUT
		 *********************************************************************************************/
		/* At each time step: **Statistics**
		 * t_new,
		 * dt,
		 * newton iterations per fixed point iteration,
		 * total newton iterations
		 */
		
		

		postprocess.evalProjectedSolution(); 
		// auto uoldtmp = uold;
		// uoldtmp =0.;
		// uoldtmp.axpy(1,  unew);
		// uoldtmp.axpy(0.,  uold);
		uold = unew;
		
		// GRAPHICS FOR NEW OUTPUT
		// primary variables
		// DGF_Pw dgf_Pw(subgfs_Pw, uold);
		// DGF_Sg dgf_Sg(subgfs_Sg, uold);
		// DGF_Sh dgf_Sh(subgfs_Sh, uold);
		// DGF_T dgf_T(subgfs_T, uold);
		// DGF_XCH4 dgf_XCH4(subgfs_XCH4, uold);
		// DGF_YH2O dgf_YH2O(subgfs_YH2O, uold);
		// DGF_XC dgf_XC(subgfs_XC, uold);
		// DGF_pen dgf_pen(subgfs_pen, unew);

		/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
		if (((time + dt)/(t_OP * opcount) > (1.-1.e-9)) and ((time + dt) / (t_OP * opcount)< (1. + 1.e-9)))
		{
			// primary variables
			// slpproj.apply();
			postprocess.evaluate();
			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << " ----processor: " << helper.rank() << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
			}
			// timecount = time;
			opcount = opcount + 1;
		}

		//		PREPARE FOR NEXT TIME INTEGRATION
		//		1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		
		//		2. ADVANCE TIME:
		time += dt;
		if(helper.rank()==0){
			std::cout<<" "<< std::endl;
			std::cout<< " time = " << time*Xc_t ;
			std::cout<< std::flush;
		}
		if (adaptive_time_control)
		{
			if (newton_iterations > maxAllowableIterations)
			{
				dt = std::max(dt*0.9 , dt_min);
			}
			else if (newton_iterations <= minAllowableIterations)
			{
				dt = std::min(dt * 1.2, dt_max);
			}
			if (dtFlag == -1)
			{
				dt = dtLast;//std::max(dt, dtLast);
			}
			dtFlag = 0;
		}
		else
		{
			dt = dtstart;
		}
		if(helper.rank()==0){
			std::cout << " , time+dt = " << (time + dt)*Xc_t
					  << " , opTime = "  << t_OP * opcount * Xc_t ;
			std::cout<< std::flush;
		}
		dtLast = dt;
		if ((time + dt) > (t_OP * opcount + 1.e-9) and time < (t_OP * opcount - 1.e-9))
		{
			
			dt = t_OP * opcount - time;

			if(helper.rank()==0){
				std::cout<< " , because timeNext > opNext , dt set to : " << dt*Xc_t << std::endl;
				std::cout<< std::flush;
			}
			dtFlag = -1;
		}

		if(helper.rank()==0){
			std::cout<< " , dt  : " << dt*Xc_t << std::endl;
			std::cout<<" "<< std::endl;
			std::cout << " READY FOR NEXT ITERATION. " << std::endl;
			std::cout<< std::flush;
		}
	}
};

#endif /* HYDRATE_DG_HH_ */