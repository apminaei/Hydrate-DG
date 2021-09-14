/*
 * driver.hh
 *
 */

#ifndef HYDRATE_DG_HH_
#define HYDRATE_DG_HH_
// struct ConvectionDiffusionDGMethod
// {
//   enum Type { NIPG, SIPG, IIPG };
// };
template <class GV, class PTree>
void driver_Sh(const GV &gv, // GridView
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
	auto Xc_t = property.characteristicValue.t_c;
	auto t_year_sec = 12.*30.*24*60.*60.;
	// double dt_initprob  = ptree.get("initial_problem.dt_initial",(double)0.0001);
	// dt_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// dt_initprob *= 1./Xc_t; /*ndim*/
	// double t_END_initprob  = ptree.get("initial_problem.time_end",(double)1);
	// t_END_initprob *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
	// t_END_initprob 	 *= 1./Xc_t; /*ndim*/
	//MAIN PROBLEM time and dt
	dt  = ptree.get("time.dt_initial",(double)0.001);
	dt *= t_year_sec; /*convert to seconds*/
	dt *= 1./Xc_t; /*ndim*/
	double t_END  = ptree.get("time.time_end",(double)300);
	t_END *= t_year_sec; /*convert to seconds*/
	t_END *= 1./Xc_t; /*ndim*/
	//t_END += t_END_initprob; //Total simulation time
	// output time interval
	double t_OP   = ptree.get("output.time_interval",(double)1000);
	t_OP *= t_year_sec; /*convert to seconds*/
	t_OP *= 1./Xc_t; /*ndim*/
	//adaptive time control
	bool adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
	double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6);
	dt_min *= t_year_sec; /*convert to seconds*/
	dt_min *= 1./Xc_t; /*ndim*/
	double dt_max = ptree.get("adaptive_time_control.dt_max",(double)10.);
	dt_max *= t_year_sec; /*convert to seconds*/
	dt_max *= 1./Xc_t; /*ndim*/

	
	Real dtstart = dt;
	Real time_op = time;
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

	//	GFS
#ifdef PARALLEL
	typedef Dune::PDELab::OverlappingConformingDirichletConstraints CON0;
#else
	typedef Dune::PDELab::ConformingDirichletConstraints CON0;	// pure Neumann: no constraints
#endif									
	typedef Dune::PDELab::ISTL::VectorBackend<> VBE0;	// default block size: 1
	//	COMPOSITE GFS
	using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;//  block size -> numOfPVs

	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_P, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_Pw;// FEM_P;
	FEM_Pw fem_pw;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Sh, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_Sh;
	FEM_Sh fem_sh;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Sg, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_Sg;
	FEM_Sg fem_sg;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_T, dim, Dune::PDELab::QkDGBasisPolynomial::legendre> FEM_T; 
	FEM_T fem_t;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_XCH4, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_XCH4; 
	FEM_XCH4 fem_xch4;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_XC, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_XC; 
	FEM_XC fem_xc;
	typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord, Real, degree_Y, dim, Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM_YH2O; 
	FEM_YH2O fem_yh2o;

	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Pw, CON0, VBE0> GFS_Pw; // gfs
	GFS_Pw gfs_pw(gv, fem_pw);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Sh, CON0, VBE0> GFS_Sh; // gfs
	GFS_Sh gfs_sh(gv, fem_sh);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_Sg, CON0, VBE0> GFS_Sg; // gfs
	GFS_Sg gfs_sg(gv, fem_sg);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_T, CON0, VBE0> GFS_T; // gfs
	GFS_T gfs_t(gv, fem_t);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_XCH4, CON0, VBE0> GFS_XCH4; // gfs
	GFS_XCH4 gfs_xch4(gv, fem_xch4);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_XC, CON0, VBE> GFS_XC; // gfs
	GFS_XC gfs_xc(gv, fem_xc);
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM_YH2O, CON0, VBE0> GFS_YH2O; // gfs
	GFS_YH2O gfs_yh2o(gv, fem_yh2o);

	using VBE1 = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;//  block size -> numOfPVs

	// gfs for composite system Pw,Sg,Sh,T,XCH4,YH2O,XC
	typedef Dune::PDELab::CompositeGridFunctionSpace<VBE1, 
													 Dune::PDELab::EntityBlockedOrderingTag,
													 GFS_Pw, GFS_Sg, GFS_Sh, GFS_T, GFS_XCH4,GFS_YH2O> 
													 GFS;
	GFS gfs(gfs_pw,  gfs_sg, gfs_sh, gfs_t, gfs_xch4, gfs_yh2o);
	
	// typedef typename GFS_Sh::template ConstraintsContainer<Real>::Type CC_Sh;
    // CC_Sh cc_sh;
    // typedef typename GFS_Pw::template ConstraintsContainer<Real>::Type CC_Pw;
    // CC_Pw cc_pw;
	// typedef typename GFS_Sg::template ConstraintsContainer<Real>::Type CC_Sg;
    // CC_Sg cc_sg;
	// typedef typename GFS_T::template ConstraintsContainer<Real>::Type CC_T;
    // CC_T cc_t;
    // typedef typename GFS_XCH4::template ConstraintsContainer<Real>::Type CC_XCH4;
    // CC_XCH4 cc_xch4;
	typedef typename GFS_XC::template ConstraintsContainer<Real>::Type CC_XC;
    CC_XC cc_xc;
	// typedef typename GFS_YH2O::template ConstraintsContainer<Real>::Type CC_YH2O;
    // CC_YH2O cc_yh2o;
	typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    CC cc;

	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
	using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Pw>>;
    using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;
    SUBGFS_Pw    subgfs_Pw(gfs);
	using PathSg = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Sg>>;
    using SUBGFS_Sg = Dune::PDELab::GridFunctionSubSpace<GFS,PathSg>;
    SUBGFS_Sg    subgfs_Sg(gfs);
	using PathSh = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Sh>>;
    using SUBGFS_Sh = Dune::PDELab::GridFunctionSubSpace<GFS,PathSh>;
    SUBGFS_Sh    subgfs_Sh(gfs);
	using PathT = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_T>>;
    using SUBGFS_T = Dune::PDELab::GridFunctionSubSpace<GFS,PathT>;
    SUBGFS_T    subgfs_T(gfs);
	using PathXCH4 = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XCH4>>;
    using SUBGFS_XCH4 = Dune::PDELab::GridFunctionSubSpace<GFS,PathXCH4>;
    SUBGFS_XCH4    subgfs_XCH4(gfs);
	using PathYH2O = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_YH2O>>;
    using SUBGFS_YH2O = Dune::PDELab::GridFunctionSubSpace<GFS,PathYH2O>;
    SUBGFS_YH2O    subgfs_YH2O(gfs);
	// using Pathc = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XC>>;
    // using SUBGFS_XC = Dune::PDELab::GridFunctionSubSpace<GFS,Pathc>;
    // SUBGFS_XC    subgfs_XC(gfs);

	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	// using U_Sh = Dune::PDELab::Backend::Vector<GFS_Sh, double>;
	// U_Sh uold_sh(gfs_sh, 0.0);
	// U_Sh unew_sh(gfs_sh, 0.0);
	// using U_Sg = Dune::PDELab::Backend::Vector<GFS_Sg, double>;
	// U_Sg uold_sg(gfs_sg, 0.0);
	// U_Sg unew_sg(gfs_sg, 0.0);
	// using U_Pw = Dune::PDELab::Backend::Vector<GFS_Pw, double>;
	// U_Pw uold_pw(gfs_pw, 0.0);
	// U_Pw unew_pw(gfs_pw, 0.0);
	// using U_T = Dune::PDELab::Backend::Vector<GFS_T, double>;
	// U_T uold_t(gfs_t, 0.0);
	// U_T unew_t(gfs_t, 0.0);
	// using U_XCH4 = Dune::PDELab::Backend::Vector<GFS_XCH4, double>;
	// U_XCH4 uold_xch4(gfs_xch4, 0.0);
	// U_XCH4 unew_xch4(gfs_xch4, 0.0);
	// using U_YH2O = Dune::PDELab::Backend::Vector<GFS_YH2O, double>;
	// U_YH2O uold_yh2o(gfs_yh2o, 0.0);
	// U_YH2O unew_yh2o(gfs_yh2o, 0.0);
	using U_XC = Dune::PDELab::Backend::Vector<GFS_XC, double>;
	U_XC uold_xc(gfs_xc, 0.0);
	U_XC unew_xc(gfs_xc, 0.0);
	using U = Dune::PDELab::Backend::Vector<GFS, double>;
	U uold(gfs, 0.0);
	U unew(gfs, 0.0);

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

	typedef Dune::PDELab::CompositeGridFunction<Pw_InitialType,
												Sg_InitialType,
												Sh_InitialType,
												T_InitialType,
												XCH4_InitialType,
												YH2O_InitialType> InitialType;
	InitialType initial(Pw_initial, Sg_initial, Sh_initial, T_initial, XCH4_initial, YH2O_initial, XC_initial);
	Dune::PDELab::interpolate(initial, gfs, uold); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(Pw_initial, gfs_pw, uold_pw); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(Sg_initial, gfs_sg, uold_sg); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(Sh_initial, gfs_sh, uold_sh); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(T_initial, gfs_t, uold_t); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(XCH4_initial, gfs_xch4, uold_xch4); // Initialize the solution at t=0 (uold) with the given initial values
	Dune::PDELab::interpolate(XC_initial, gfs_xc, uold_xc); // Initialize the solution at t=0 (uold) with the given initial values
	// Dune::PDELab::interpolate(YH2O_initial, gfs_yh2o, uold_yh2o); // Initialize the solution at t=0 (uold) with the given initial values
	
	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	double method_g = -1.;
	double method_w = -1.;
	double method_T = 1.;
	double method_x = 1.;
	double method_y = 1.;
	double alpha_g = 1.e0;
	double alpha_w = 1.e0;
	double alpha_s = 1.e0;
	double alpha_T = 1.e0;
	double alpha_x = 1.e0;
	double alpha_y = 1.e0;
	double intorder=4;

	typedef ProblemBoundaryConditions<GV,Properties> BoundaryConditions ;
	BoundaryConditions bc( gv,property ) ;



	// typedef LocalOperator_Sh<GV, Properties, BoundaryConditions, U_Sh, GFS_Sh, U, GFS, U_T, GFS_T, FEM_Sh> LOP_Sh; // spatial part
	// LOP_Sh lop_sh(gv, property, bc, &unew_sh, gfs_sh, uold, gfs, uold_t, gfs_t,
    //             &time, &dt, intorder);

	// typedef LocalOperator_Sg<GV, Properties, BoundaryConditions, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_Sg> LOP_Sg; // spatial part
	// LOP_Sg lop_sg(gv, property, bc, uold_sh, gfs_sh, &unew_sg, gfs_sg, uold_pw,  gfs_pw,
    //                 uold_t, gfs_t, uold_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, &time, &dt, intorder, method_g, alpha_g);

	// typedef LocalOperator_Pw<GV, Properties, BoundaryConditions, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_Pw> LOP_Pw; // spatial part
	// LOP_Pw lop_pw(gv, property, bc, uold_sh, gfs_sh, uold_sg, gfs_sg, &unew_pw,  gfs_pw,
    //                 uold_t, gfs_t, uold_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, &time, &dt, intorder, method_w, alpha_w);

	// typedef LocalOperator_T<GV, Properties, BoundaryConditions, U_Sh, GFS_Sh, U, GFS, U_T, GFS_T, FEM_T> LOP_T; // spatial part
	// LOP_T lop_t(gv, property, bc, uold_sh, gfs_sh, uold, gfs, &unew_t, gfs_t,
    //                 &time, &dt, intorder, method_T, alpha_T);

	typedef LocalOperator_XC<GV, Properties, BoundaryConditions, U, GFS, U_XC, GFS_XC, FEM_XC> LOP_XC; // spatial part
	LOP_XC lop_xc(gv, property, bc, uold, gfs, 
                    &unew_xc, gfs_xc, &time, &dt, intorder, method_x, alpha_x);

	// typedef LocalOperator_XCH4<GV, Properties, BoundaryConditions, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_XCH4> LOP_XCH4; // spatial part
	// LOP_XCH4 lop_xch4(gv, property, bc, uold_sh, gfs_sh, uold_sg, gfs_sg, uold_pw,  gfs_pw,
    //                 uold_t, gfs_t, &unew_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, &time, &dt, intorder, method_x, alpha_x);

	typedef LocalOperator_2comps<GV, Properties, BoundaryConditions, U, GFS, U_XC, GFS_XC, 
          					FEM_Pw, FEM_Sg, FEM_Sh, FEM_T,FEM_XCH4, FEM_YH2O> LOP; // spatial part
	LOP lop(gv, property, bc, &unew, gfs, uold_xc, gfs_xc,
			&time, &dt, intorder, method_g, method_w, method_x, method_y, alpha_g, alpha_w, alpha_x, alpha_y);

	// typedef TimeOperator_Sh<GV, Properties, U, GFS, U_T, GFS_T, FEM_Sh> TOP_Sh; // spatial part
	// TOP_Sh top_sh(gv, property, uold, gfs, uold_t, gfs_t, intorder);

	// typedef TimeOperator_Sg<GV, Properties, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_Sg> TOP_Sg; // spatial part
	// TOP_Sg top_sg(gv, property, uold_sh, gfs_sh, &unew_sg, gfs_sg, uold_pw,  gfs_pw,
    //                 uold_t, gfs_t, uold_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, intorder);

	// typedef TimeOperator_Pw<GV, Properties, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_Pw> TOP_Pw; // spatial part
	// TOP_Pw top_pw(gv, property, uold_sh, gfs_sh, uold_sg, gfs_sg, &unew_pw,  gfs_pw,
    //                 uold_t, gfs_t, uold_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, intorder);

	// typedef TimeOperator_T<GV, Properties, U_Sh, GFS_Sh, U, GFS, FEM_T> TOP_T; // spatial part
	// TOP_T top_t(gv, property, uold_sh, gfs_sh, uold, gfs, 
    //                  intorder);

	typedef TimeOperator_XC<GV, Properties, U, GFS, FEM_XC> TOP_XC; // spatial part
	TOP_XC top_xc(gv, property, uold, gfs, intorder);

	// typedef TimeOperator_XCH4<GV, Properties, U_Sh, GFS_Sh, U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T,
    //       					U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_XCH4> TOP_XCH4; // spatial part
	// TOP_XCH4 top_xch4(gv, property, uold_sh, gfs_sh, uold_sg, gfs_sg, uold_pw,  gfs_pw,
    //                 uold_t, gfs_t, &unew_xch4, gfs_xch4, uold_yh2o, gfs_yh2o,
    //                 uold_xc, gfs_xc, intorder);

	typedef TimeOperator_2comps<GV, Properties, U_XC, GFS_XC, FEM_YH2O> TOP; // spatial part
	TOP top(gv, property, uold_xc, gfs_xc,  intorder);


	typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
	MBE mbe(100);

	// typedef Dune :: PDELab :: GridOperator <
	// GFS , GFS , //  ansatz and test space 
	// LOP , 	//  local operator 
	// MBE , 	//  matrix backend 
	// RF , RF , RF , 	//  domain , range , jacobian field type 
	// CC , CC	//  constraints for ansatz and test space 
	// > GO ;

	// typedef Dune::PDELab::GridOperator<GFS_Sh, GFS_Sh, LOP_Sh, MBE, Real, Real, Real, CC_Sh, CC_Sh> GOLOP_Sh;
	// GOLOP_Sh goLOP_sh(gfs_sh, cc_sh, gfs_sh, cc_sh, lop_sh, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_Sg, GFS_Sg, LOP_Sg, MBE, Real, Real, Real, CC_Sg, CC_Sg> GOLOP_Sg;
	// GOLOP_Sg goLOP_sg(gfs_sg, cc_sg, gfs_sg, cc_sg, lop_sg, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_Pw, GFS_Pw, LOP_Pw, MBE, Real, Real, Real, CC_Pw, CC_Pw> GOLOP_Pw;
	// GOLOP_Pw goLOP_pw(gfs_pw, cc_pw, gfs_pw, cc_pw, lop_pw, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_T, GFS_T, LOP_T, MBE, Real, Real, Real, CC_T, CC_T> GOLOP_T;
	// GOLOP_T goLOP_t(gfs_t, cc_t, gfs_t, cc_t, lop_t, mbe);

	typedef Dune::PDELab::GridOperator<GFS_XC, GFS_XC, LOP_XC, MBE, Real, Real, Real, CC_XC, CC_XC> GOLOP_XC;
	GOLOP_XC goLOP_xc(gfs_xc, cc_xc, gfs_xc, cc_xc, lop_xc, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_XCH4, GFS_XCH4, LOP_XCH4, MBE, Real, Real, Real, CC_XCH4, CC_XCH4> GOLOP_XCH4;
	// GOLOP_XCH4 goLOP_xch4(gfs_xch4, cc_xch4, gfs_xch4, cc_xch4, lop_xch4, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_YH2O, GFS_YH2O, LOP_YH2O, MBE, Real, Real, Real, CC_YH2O, CC_YH2O> GOLOP_YH2O;
	// GOLOP_YH2O goLOP_yh2o(gfs_yh2o, cc_yh2o, gfs_yh2o, cc_yh2o, lop_yh2o, mbe);

	typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GOLOP;
	GOLOP goLOP(gfs, cc, gfs, cc, lop, mbe);
	// How well did we estimate the number of entries per matrix row?
	// => print Jacobian pattern statistics
	// typename GOLOP::Traits::Jacobian jac(goLOP);
	// std::cout << " LOP DONE ! " << std::endl;


	// typedef Dune::PDELab::GridOperator<GFS_Sh, GFS_Sh, TOP_Sh, MBE, Real, Real, Real, CC_Sh, CC_Sh> GOTOP_Sh;
	// GOTOP_Sh goTOP_sh(gfs_sh, cc_sh, gfs_sh, cc_sh, top_sh, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_Sg, GFS_Sg, TOP_Sg, MBE, Real, Real, Real, CC_Sg, CC_Sg> GOTOP_Sg;
	// GOTOP_Sg goTOP_sg(gfs_sg, cc_sg, gfs_sg, cc_sg, top_sg, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_Pw, GFS_Pw, TOP_Pw, MBE, Real, Real, Real, CC_Pw, CC_Pw> GOTOP_Pw;
	// GOTOP_Pw goTOP_pw(gfs_pw, cc_pw, gfs_pw, cc_pw, top_pw, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_T, GFS_T, TOP_T, MBE, Real, Real, Real, CC_T, CC_T> GOTOP_T;
	// GOTOP_T goTOP_t(gfs_t, cc_t, gfs_t, cc_t, top_t, mbe);

	typedef Dune::PDELab::GridOperator<GFS_XC, GFS_XC, TOP_XC, MBE, Real, Real, Real, CC_XC, CC_XC> GOTOP_XC;
	GOTOP_XC goTOP_xc(gfs_xc, cc_xc, gfs_xc, cc_xc, top_xc, mbe);

	// typedef Dune::PDELab::GridOperator<GFS_XCH4, GFS_XCH4, TOP_XCH4, MBE, Real, Real, Real, CC_XCH4, CC_XCH4> GOTOP_XCH4;
	// GOTOP_XCH4 goTOP_xch4(gfs_xch4, cc_xch4, gfs_xch4, cc_xch4, top_xch4, mbe);

	typedef Dune::PDELab::GridOperator<GFS, GFS, TOP, MBE, Real, Real, Real, CC, CC> GOTOP;
	GOTOP goTOP(gfs, cc, gfs, cc, top, mbe);


	// typedef Dune::PDELab::OneStepGridOperator<GOLOP_T, GOTOP_T> IGO_T;
	// IGO_T igo_t(goLOP_t, goTOP_t);

	typedef Dune::PDELab::OneStepGridOperator<GOLOP_XC, GOTOP_XC> IGO_XC;
	IGO_XC igo_xc(goLOP_xc, goTOP_xc);

	// typedef Dune::PDELab::OneStepGridOperator<GOLOP_Sh, GOTOP_Sh> IGO_Sh;
	// IGO_Sh igo_sh(goLOP_sh, goTOP_sh);

	// typedef Dune::PDELab::OneStepGridOperator<GOLOP_Sg, GOTOP_Sg> IGO_Sg;
	// IGO_Sg igo_sg(goLOP_sg, goTOP_sg);

	// typedef Dune::PDELab::OneStepGridOperator<GOLOP_Pw, GOTOP_Pw> IGO_Pw;
	// IGO_Pw igo_pw(goLOP_pw, goTOP_pw);

	// typedef Dune::PDELab::OneStepGridOperator<GOLOP_XCH4, GOTOP_XCH4> IGO_XCH4;
	// IGO_XCH4 igo_xch4(goLOP_xch4, goTOP_xch4);

	typedef Dune::PDELab::OneStepGridOperator<GOLOP, GOTOP> IGO;
	IGO igo(goLOP, goTOP);

	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	// LS ls(gfs, 100, 2, true, true);
	std::cout << " IGO DONE ! " << std::endl;
	
		// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_Sh, CC_Sh> LS_Sh;
	// LS_Sh ls_sh(gfs_sh, cc_sh, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO_Sh> LS_Sh;
	// LS_Sh ls_sh(gfs_sh, 100, 1, true, true);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_Sg, CC_Sg> LS_Sg;
	// LS_Sg ls_sg(gfs_sg, cc_sg, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_Pw, CC_Pw> LS_Pw;
	// LS_Pw ls_pw(gfs_pw, cc_pw, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_T, CC_T> LS_T;
	// LS_T ls_t(gfs_t, cc_t, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO_T> LS_T;
	// LS_T ls_t(gfs_t, 100, 1, true, true);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_XC, CC_XC> LS_XC;
	// LS_XC ls_xc(gfs_xc, cc_xc, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_XCH4, CC_XCH4> LS_XCH4;
	// LS_XCH4 ls_xch4(gfs_xch4, cc_xch4, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS_YH2O, CC_YH2O> LS_YH2O;
	// LS_YH2O ls_yh2o(gfs_yh2o, cc_yh2o, 100, 1);
	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_SuperLU<GFS, CC> LS;
	// LS ls(gfs, cc, 100, 1);
	ypedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO_XC> LS_XC; //works
	LS_XC ls_xc(gfs_xc, 100, 1, true, true);


	typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
	LS ls(gfs, 100, 1, true, true);

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

	// typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0< GFS, CC > LS;
	// LS ls(gfs, cc, 100, 1, 10);

	// typedef Dune::PDELab::ISTLBackend_BCGS_AMG_ILU0<IGO> LS; // does not work for more than 1 level coarsening from AMG
	// LS ls(gfs,500,2,true,true);

	// typedef Dune::PDELab::ISTLBackend_OVLP_BCGS_ILUn<GFS, CC> LS; //works
	// LS ls(gfs, cc);

	// typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<IGO> LS; // should be checked
	// int verbose = 0;
	// if (gfs.gridView().comm().rank() == 0)
	// 	verbose = 1;
	// LS ls(gfs, 100, verbose);

	auto param = ls.parameters();
	//param.setMaxLevel(3); // max number of coarsening levels
	param.setCoarsenTarget(100000); // max DoF at coarsest level
	ls.setParameters(param);
	auto param_t = ls_t.parameters();
	//param.setMaxLevel(3); // max number of coarsening levels
	// param_t.setCoarsenTarget(100000); // max DoF at coarsest level
	// ls_t.setParameters(param_t);
	// auto param_sh = ls_sh.parameters();
	// //param.setMaxLevel(3); // max number of coarsening levels
	// param_sh.setCoarsenTarget(100000); // max DoF at coarsest level
	// ls_sh.setParameters(param_sh);

	std::cout << " PARALLEL LS DONE ! " << std::endl;

#else
	// typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	// LS ls(1000, true);

	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS;
	LS ls;

	typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS_XC;
	LS_XC ls_xc;
	std::cout << " LS DONE ! " << std::endl;
#endif


	// Linear and PDESOLVER
<<<<<<< HEAD
	using PDESOLVER_Sh = Dune::PDELab::Newton< IGO_Sh, LS, U_Sh >;
    PDESOLVER_Sh pdesolver_sh( igo_sh, ls );
=======
	// using PDESOLVER_Sh = Dune::PDELab::Newton< IGO_Sh, LS, U_Sh >;
    // PDESOLVER_Sh pdesolver_sh( igo_sh, ls );
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
	// using PDESOLVER_Sg = Dune::PDELab::Newton< IGO_Sg, LS_Sg, U_Sg >;
    // PDESOLVER_Sg pdesolver_sg( igo_sg, ls_sg );
	// using PDESOLVER_Pw = Dune::PDELab::Newton< IGO_Pw, LS_Pw, U_Pw >;
    // PDESOLVER_Pw pdesolver_pw( igo_pw, ls_pw );
<<<<<<< HEAD
	using PDESOLVER_T = Dune::PDELab::Newton< IGO_T, LS, U_T >;
    PDESOLVER_T pdesolver_t( igo_t, ls );
	// using PDESOLVER_XC = Dune::PDELab::Newton< IGO_XC, LS_XC, U_XC >;
    // PDESOLVER_XC pdesolver_xc( igo_xc, ls_xc );
=======
	// using PDESOLVER_T = Dune::PDELab::Newton< IGO_T, LS, U_T >;
    // PDESOLVER_T pdesolver_t( igo_t, ls );
	using PDESOLVER_XC = Dune::PDELab::NewtonMethod< IGO_XC, LS_XC >;
    PDESOLVER_XC pdesolver_xc( igo_xc, ls_xc );
	pdesolver_xc.setParameters(ptree.sub("newton"));
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
	// using PDESOLVER_XCH4 = Dune::PDELab::Newton< IGO_XCH4, LS_XCH4, U_XCH4 >;
    // PDESOLVER_XCH4 pdesolver_xch4( igo_xch4, ls_xch4 );
	using PDESOLVER = Dune::PDELab::NewtonMethod< IGO, LS >;
    PDESOLVER pdesolver( igo, ls);
	pdesolver.setParameters(ptree.sub("newton"));
	
	// Dune::PDELab::OneStepMethod<Real, IGO_Sg, PDESOLVER, U, U> osm_sg(*pmethod, igo_sg, pdesolver);
	
	// Dune::PDELab::OneStepMethod<Real, IGO_Sh, PDESOLVER, U, U> osm_sh(*pmethod, igo_sh, pdesolver);

	// Dune::PDELab::OneStepMethod<Real, IGO_Pw, PDESOLVER, U, U> osm_pw(*pmethod, igo_pw, pdesolver);

	// Dune::PDELab::OneStepMethod<Real, IGO_T, PDESOLVER, U, U> osm_t(*pmethod, igo_t, pdesolver);

	// Dune::PDELab::OneStepMethod<Real, IGO_C, PDESOLVER, U, U> osm_c(*pmethod, igo_c, pdesolver);
	// std::cout << " OSM DONE ! " << std::endl;




    //    SELECT SOLVER FOR NON-LINEAR PROBLEM
    // using PDESOLVER = Dune::PDELab::Newton< IGO, LS, U >;
    // PDESOLVER pdesolver( igo, ls );
    // select control parameters for non-linear PDE-solver
	// pdesolver_sh.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_sh.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_sh.setReassembleThreshold(0.0);
    // pdesolver_sh.setVerbosityLevel(2);
    // pdesolver_sh.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_sh.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_sh.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_sh.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_sh.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4)); 
	
	// pdesolver_sg.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_sg.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_sg.setReassembleThreshold(0.0);
    // pdesolver_sg.setVerbosityLevel(2);
    // pdesolver_sg.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_sg.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_sg.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_sg.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_sg.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4)); 
	
	// pdesolver_pw.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_pw.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_pw.setReassembleThreshold(0.0);
    // pdesolver_pw.setVerbosityLevel(2);
    // pdesolver_pw.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_pw.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_pw.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_pw.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_pw.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));

	// pdesolver_t.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_t.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_t.setReassembleThreshold(0.0);
    // pdesolver_t.setVerbosityLevel(2);
    // pdesolver_t.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_t.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_t.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_t.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_t.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));

	// pdesolver_xc.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_xc.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_xc.setReassembleThreshold(0.0);
    // pdesolver_xc.setVerbosityLevel(2);
    // pdesolver_xc.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_xc.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_xc.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_xc.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_xc.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));

	// pdesolver_xch4.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver_xch4.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver_xch4.setReassembleThreshold(0.0);
    // pdesolver_xch4.setVerbosityLevel(2);
    // pdesolver_xch4.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver_xch4.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver_xch4.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver_xch4.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver_xch4.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));

	// pdesolver.setLineSearchStrategy(ptree.get("newton.line_search_strategy",(std::string)"noLineSearch"));//Strategy {  hackbuschReusken, hackbuschReuskenAcceptBest }
    // //pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::hackbuschReuskenAcceptBest);
	// pdesolver.setReassembleThreshold(0.0);
    // pdesolver.setVerbosityLevel(2);
    // pdesolver.setReduction(ptree.get("newton.reduction",(double)1e-5));
    // pdesolver.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
    // pdesolver.setForceIteration(ptree.get("newton.force_iterations",(bool)true));
	// pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));

	//TODO: CHECK NEW NEWTON PARAMS
	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	// typedef Dune::PDELab::NewtonMethod<IGO, LS> PDESOLVER;
	// PDESOLVER pdesolver(igo, ls);
	// //	select control parameters for non-linear PDE-solver
	// typedef Dune::PDELab::LineSearchNone<PDESOLVER> lineSearchStrategy;
	// //typedef Dune::PDELab::LineSearchHackbuschReusken<PDESOLVER> lineSearchStrategy;
	// lineSearchStrategy linesearchstrategy(pdesolver);
	// //pdesolver.setParameters(ptree.sub("newton"));
	// pdesolver.setVerbosityLevel(2);
	// pdesolver.setReduction(ptree.get("newton.reduction",(double)1e-5));
	// pdesolver.setMinLinearReduction(ptree.get("newton.min_linear_reduction",(double)1.e-9));
	// pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-4));
	
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
	Dune::PDELab::TimeSteppingParameterInterface<Real> *pmethod_sh = &method1;
	Dune::PDELab::TimeSteppingParameterInterface<Real> *pmethod_t = &method1;

	// Dune::PDELab::OneStepMethod<Real, IGO_Sh, PDESOLVER_Sh, U_Sh, U_Sh> osm_sh(*pmethod_sh, igo_sh, pdesolver_sh);
	// osm_sh.setVerbosityLevel(2);
	// Dune::PDELab::OneStepMethod<Real, IGO_Sg, PDESOLVER_Sg, U_Sg, U_Sg> osm_sg(*pmethod, igo_sg, pdesolver_sg);
	// osm_sg.setVerbosityLevel(2);
	// Dune::PDELab::OneStepMethod<Real, IGO_Pw, PDESOLVER_Pw, U_Pw, U_Pw> osm_pw(*pmethod, igo_pw, pdesolver_pw);
	// osm_pw.setVerbosityLevel(2);
	// Dune::PDELab::OneStepMethod<Real, IGO_T, PDESOLVER_T, U_T, U_T> osm_t(*pmethod_t, igo_t, pdesolver_t);
	// osm_t.setVerbosityLevel(2);
	Dune::PDELab::OneStepMethod<Real, IGO_XC, PDESOLVER_XC, U_XC, U_XC> osm_xc(*pmethod, igo_xc, pdesolver_xc);
	osm_xc.setVerbosityLevel(2);
	// Dune::PDELab::OneStepMethod<Real, IGO_XCH4, PDESOLVER_XCH4, U_XCH4, U_XCH4> osm_xch4(*pmethod, igo_xch4, pdesolver_xch4);
	// osm_xch4.setVerbosityLevel(2);
	Dune::PDELab::OneStepMethod<Real, IGO, PDESOLVER, U, U> osm(*pmethod, igo, pdesolver);
	osm.setVerbosityLevel(2);
	std::cout << " OSM DONE ! " << std::endl;

	//	GRAPHICS FOR INITIAL GUESS
	// primary variables  // Make a grid function out of it
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Pw, U> DGF_Pw;
	DGF_Pw dgf_pw(subgfs_Pw, uold);	
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sg, U> DGF_Sg;
	DGF_Sg dgf_sg(subgfs_Sg, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_Sh, U> DGF_Sh;
	DGF_Sh dgf_sh(subgfs_Sh, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_T, U> DGF_T;
	DGF_T dgf_t(subgfs_T, uold);	
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_XCH4, U> DGF_XCH4;
	DGF_XCH4 dgf_xch4(subgfs_XCH4, uold);
	typedef Dune::PDELab::DiscreteGridFunction<SUBGFS_YH2O, U> DGF_YH2O;
	DGF_YH2O dgf_yh2o(subgfs_YH2O, uold);
	typedef Dune::PDELab::DiscreteGridFunction<GFS_XC, U_XC> DGF_XC;
	DGF_XC dgf_xc(gfs_xc, uold_xc);
	//Dune::FieldVector<double, 1> y=0.;

	// for ( const auto & e : elements ( gv )) {
		
	// 	  auto geo = e.geometry();
	// 	  auto geo_ref = referenceElement(geo);
	// 	  for ( int i = 0; i < geo_ref.size(1); i++) {
	// 		//idxSet.subIndex (e , i , c ); // index of subentity
	// 		auto geo_ref_pos = geo.corner(i);
	// 		//x = geo.global(geo_ref_pos);
	// 		dgf_pw.evaluate(e, geo_ref_pos, y);
	// 		std::cout << "   --" << 
	// 		geo_ref_pos << "     " << y << "   ";
	// 	  }
	// 	  std::cout << " "<<std::endl;
	// }


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


	int subsampling = 1;
	Dune::RefinementIntervals RefInt(subsampling);

	using VTKWRITER = Dune::SubsamplingVTKWriter<GV> ;
	VTKWRITER vtkwriter(gv, RefInt, false, Dune::VTK::Precision::float32);
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV> ;
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter), fileName, pathName, "");

	// add data field for all components of the space to the VTK writer
	// primary variables
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw>>(dgf_pw, "Pw"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_sh, "Sh"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sg>>(dgf_sg, "Sg"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_t, "T"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XCH4>>(dgf_xch4, "XCH4"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_YH2O>>(dgf_yh2o, "YH2O"));
	vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XC>>(dgf_xc, "XC"));

	vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
	vtkSequenceWriter.clear();

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
	// unew_sh = uold_sh;
	// unew_sg = uold_sg;
	// unew_pw = uold_pw;
	// unew_t = uold_t;
	// unew_xch4 = uold_xch4;
	unew = uold;
	unew_xc = uold_xc;
	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;

	bool exceptionCaught = false;

	int newton_iterations = 0;
	double newton_first_defect = 0.;
	double newton_defect = 0.;
<<<<<<< HEAD
	
=======

>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
	//	BEGIN TIME LOOP
	while ( time < (t_END - 1e-3/Xc_t))
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
			// auto current_time = time;
			// auto current_dt = dt;
			// osm_sh.apply( time, dt, uold_sh, unew_sh );
			// newton_iterations = osm_sh.getPDESolver().result().iterations;
		
			// uold_sh = unew_sh;
			// // time = current_time;
			// // dt = current_dt; 
			
			// DGF_Sh dgf_sh(gfs_sh, unew_sh);
			
			
			// if(helper.rank()==0){
			// 	// std::cout << "current_time = " << current_time  << "   time = " << time<< std::endl;
			// 	// std::cout << "current_dt = " << current_dt  << "   dt = " << dt<< std::endl;
			// 	std::cout << "========== Sh DONE!" <<  " ======== " << std::endl;
			// }


			// osm_t.apply( time, dt, uold_t, unew_t );
			// newton_iterations = osm_t.getPDESolver().result().iterations;

			// uold_t = unew_t;
			// // time = current_time;
			// // dt = current_dt; 
			// DGF_T dgf_t(gfs_t, unew_t);
			
			// if(helper.rank()==0){
			// 	// std::cout << "current_time = " << current_time  << "   time = " << time<< std::endl;
			// 	// std::cout << "current_dt = " << current_dt  << "   dt = " << dt<< std::endl;
			// 	std::cout << "========== T DONE!" <<  " ======== " << std::endl;
			// }
			
			osm.apply( time, dt, uold, unew );
			newton_iterations = osm.getPDESolver().result().iterations;
			newton_first_defect = osm.getPDESolver().result().first_defect;
			newton_defect = osm.getPDESolver().result().defect;
<<<<<<< HEAD
		
			uold = unew;
			// time = current_time;
			// dt = current_dt; 
			// DGF_Pw dgf_pw(subgfs_Pw, unew);
			// DGF_Sg dgf_sg(subgfs_Sg, unew);
			// DGF_XCH4 dgf_xch4(subgfs_XCH4, unew);
			// DGF_YH2O dgf_yh2o(subgfs_YH2O, unew);
			// DGF_XC dgf_xc(subgfs_XC, unew);
=======

			newton_iterations = osm.getPDESolver().result().iterations;
			newton_first_defect = osm.getPDESolver().result().first_defect;
			newton_defect = osm.getPDESolver().result().defect;
            auto jacobian = osm.getPDESolver().getJacobian();
			if(helper.rank()==0 &&  (opcount%10==0)){
			Dune::writeMatrixToMatlab ( Dune::PDELab::Backend::native(jacobian), jacPath+"jacobian");
			Dune::writeVectorToMatlab(Dune::PDELab::Backend::native(unew),jacPath+"solution");
			}
			auto newton_defects = osm.getPDESolver().result().defects;
			auto u_norm_two = unew.two_norm();
			auto u_norm_one = unew.one_norm();
			auto u_norm_infinity = unew.infinity_norm();
			if(helper.rank()==0 &&  (newton_iterations>1)){//((time+dt )/(t_OP * opcount) > (1.-1.e-6)) and ((time+dt ) / (t_OP * opcount)< (1. + 1.e-6))
				std::string s_OP = std::__cxx11::to_string(time);
				std::string parameters_file = pathNameDefects;
				parameters_file +="/";
				parameters_file +=fileNameDefects;
				parameters_file +="_";
				parameters_file +=s_OP;
				parameters_file += ".txt";
				property.ReportNewton( parameters_file,
							time /*s*/,
							dt /*s*/,
							newton_iterations, newton_defects,
							u_norm_two, u_norm_one,  u_norm_infinity);
			}
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
			if(helper.rank()==0){
				// std::cout << "current_time = " << current_time  << "   time = " << time<< std::endl;
				// std::cout << "current_dt = " << current_dt  << "   dt = " << dt<< std::endl;
				std::cout << "========== 3 Comps DONE!" <<  " ======== " << std::endl;
			}
			osm_xc.apply( time, dt, uold_xc, unew_xc );
			// newton_iterations = osm_xc.getPDESolver().result().iterations;
			

			if(helper.rank()==0){
				// std::cout << "current_time = " << current_time  << "   time = " << time<< std::endl;
				// std::cout << "current_dt = " << current_dt  << "   dt = " << dt<< std::endl;
				std::cout << "========== XC DONE!" <<  " ======== " << std::endl;
			}	
			
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
			if( (dt) > (1e-3/Xc_t) ){

				if(helper.rank()==0){
					std::cout << "Catched Error, Dune reported error: " << e << std::endl;
				}

				// unew_sh = uold_sh;
				// unew_sg = uold_sg;
				// unew_pw = uold_pw;
				// unew_t = uold_t;
				// unew_xch4 = uold_xch4;
				unew = uold;
				unew_xc = uold_xc;

				newton_iterations = 0;

				dt *= 0.5;
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
		if(helper.rank()==0){
			std::string statistics_file = pathNameDefects;
			statistics_file +="/";
			statistics_file +=fileName;
			statistics_file +="_statistics";
			statistics_file += ".txt";
			property.ReportStatistics( 	statistics_file,
										time,
										dt,
										newton_iterations,
										newton_first_defect,
										newton_defect,
										clock_time_elapsed );
		}
		// GRAPHICS FOR NEW OUTPUT
		// primary variables
		DGF_Pw dgf_pw(subgfs_Pw, unew);
		DGF_Sg dgf_sg(subgfs_Sg, unew);
		DGF_Sh dgf_sh(subgfs_Sh, unew);
		DGF_T dgf_t(subgfs_T, unew);
		DGF_XCH4 dgf_xch4(subgfs_XCH4, unew);
		DGF_YH2O dgf_yh2o(subgfs_YH2O, unew);
		DGF_XC dgf_xc(subgfs_XC, unew_xc);

		/*********************************************************************************************
			 * OUTPUT 
			 *********************************************************************************************/
		if (((time + dt)/(t_OP * opcount) > (1.-1.e-3)) and ((time + dt) / (t_OP * opcount)< (1. + 1.e-3)))
		{
			// primary variables
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Pw>>(dgf_pw, "Pw"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sh>>(dgf_sh, "Sh"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_Sg>>(dgf_sg, "Sg"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_T>>(dgf_t, "T"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XCH4>>(dgf_xch4, "XCH4"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_YH2O>>(dgf_yh2o, "YH2O"));
			vtkSequenceWriter.addCellData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF_XC>>(dgf_xc, "XC"));

			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
			vtkSequenceWriter.clear();
			if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << " ----processor: " << helper.rank() << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
			}
			timecount = time;
			opcount = opcount + 1;
		}

		//		PREPARE FOR NEXT TIME INTEGRATION
		//		1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		// unew_sh = uold_sh;
		// unew_t = uold_t;
		unew = uold;
		unew_xc = uold_xc;
		//		2. ADVANCE TIME:
		time += dt;
		if(helper.rank()==0){
			std::cout<<" "<< std::endl;
			<< " time = " << time ;
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
<<<<<<< HEAD
			std::cout << " , time+dt = " 
			<< std::setw(12) << std::setprecision(8) << std::scientific
			<< (time + dt) * Xc_t
			<< std::setw(12) << std::setprecision(8) << std::scientific
			<< " , opTime = "  << (t_OP * opcount)* Xc_t  ;
			std::cout<< std::flush;
		}
		dtLast = dt;
		if ((time + dt) * Xc_t > (t_OP * opcount * Xc_t + 1.e-6) )
=======
			std::cout << " , time+dt = " << (time + dt)*Xc_t
					  << " , opTime = "  << t_OP * opcount * Xc_t ;
			std::cout<< std::flush;
		}
		dtLast = dt;
		if ((time + dt) > (t_OP * opcount + 1.e-6) and time < (t_OP * opcount - 1.e-6) )
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
		{
			
			dt = t_OP * opcount - time;

			if(helper.rank()==0){
<<<<<<< HEAD
				std::cout<< " , because timeNext > opNext , dt set to : " 
				<< std::setw(12) << std::setprecision(8) << std::scientific
				<< dt * Xc_t << std::endl;
=======
				std::cout<< " , because timeNext > opNext , dt set to : " << dt*Xc_t << std::endl;
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
				std::cout<< std::flush;
			}
			dtFlag = -1;
		}
 
		if(helper.rank()==0){
<<<<<<< HEAD
			std::cout<< " , dt  : "  
			<< std::setw(12) << std::setprecision(8) << std::scientific
			<< dt * Xc_t << std::endl;
=======
			std::cout<< " , dt  : " << dt*Xc_t << std::endl;
>>>>>>> 9ed57f5a8aa16f35fbe8f065875bdb4e2a910b38
			std::cout<<" "<< std::endl;
			std::cout << " READY FOR NEXT ITERATION. " << std::endl;
			std::cout<< std::flush;
		}
	}
};

#endif /* HYDRATE_DG_HH_ */
