/*
 * LocalOperator.hh
 *
 *
 */



using namespace Dune::PDELab;

// struct ConvectionDiffusionDGMethod
// {
//   enum Type { NIPG, SIPG, IIPG };
// };

template <class GV, typename Params, typename U, class GFS, class FEM_P, class FEM_Sg, class FEM_Sh, class FEM_T, class FEM_X, class FEM_Y, class FEM_C>
class LocalOperator :
                      public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::NumericalJacobianVolume<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::NumericalJacobianBoundary<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_Sg, FEM_Sh, FEM_T, FEM_X, FEM_Y, FEM_C>>,
                      public Dune::PDELab::FullSkeletonPattern, // matrix entries skeleton
                      public Dune::PDELab::FullVolumePattern,
                      public Dune::PDELab::LocalOperatorDefaultFlags,
                      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:


  const GV &gv;
  const Params&	  property;
  typedef ProblemBoundaryConditions<GV,Params> BC ;
  Dune::FieldVector<double,GV::dimension> gravity;
  U *unew;
  GFS gfs;
  double *time;
  double *dt;
  double alpha_g;
  double alpha_w;
  double alpha_s;
  double alpha_T;
  double alpha_x;
  double alpha_y;
  double method_g;
  double method_w;
  double method_T;
  double method_x;
  double method_y;
  double theta_g;
  double theta_w;
  double theta_T;
  double theta_x;
  double theta_y;
  constexpr static double eps = 1.0e-6;
  constexpr static double eps_ap	= 0.;
  constexpr static double pi = atan(1.) * 4;
  unsigned int intorder;
  double Xc_conv_m;
  double Xc_conv_h;
  double Xc_source_m;
  double Xc_source_h;
  double Xc_diff_m;
  double Xc_diff_h;
  double Xc_grav;
  double Xc_K;
  double Xc_mu;
  double Xc_rho;
  double Xc_kth;
  double Xc_v;
  double Xc_C;
  double Xc_P;
  double Xc_T;
  double Xc_t;
  double Xc_x;
  double T_ref;

public:
  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };


  typedef typename GV::IndexSet IndexSet;

  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U::template LocalView<LFSCache> VectorView;



  using LocalBasisType_Pw = typename FEM_P::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Pw = Dune::PDELab::LocalBasisCache<LocalBasisType_Pw>;
  using LocalBasisType_Sg = typename FEM_Sg::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sg = Dune::PDELab::LocalBasisCache<LocalBasisType_Sg>;
  using LocalBasisType_Sh = typename FEM_Sh::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sh = Dune::PDELab::LocalBasisCache<LocalBasisType_Sh>;
  using LocalBasisType_T = typename FEM_T::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_T = Dune::PDELab::LocalBasisCache<LocalBasisType_T>;
  using LocalBasisType_XCH4 = typename FEM_X::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_XCH4 = Dune::PDELab::LocalBasisCache<LocalBasisType_XCH4>;
  using LocalBasisType_YH2O = typename FEM_Y::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_YH2O = Dune::PDELab::LocalBasisCache<LocalBasisType_YH2O>;
  using LocalBasisType_XC = typename FEM_C::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_XC = Dune::PDELab::LocalBasisCache<LocalBasisType_XC>;
  // using LocalBasisType_pen = typename FEM_pen::Traits::FiniteElementType::Traits::LocalBasisType;
  // using Cache_pen = Dune::PDELab::LocalBasisCache<LocalBasisType_pen>;


  // In theory it is possible that one and the same local operator is
  // called first with a finite element of one type and later with a
  // finite element of another type.  Since finite elements of different
  // type will usually produce different results for the same local
  // coordinate they cannot share a cache.  Here we use a vector of caches
  // to allow for different orders of the shape functions, which should be
  // enough to support p-adaptivity.  (Another likely candidate would be
  // differing geometry types, i.e. hybrid meshes.)

  std::vector<Cache_Pw> cache_Pw;
  std::vector<Cache_Sg> cache_Sg;
  std::vector<Cache_Sh> cache_Sh;
  std::vector<Cache_T> cache_T;
  std::vector<Cache_XCH4> cache_XCH4;
  std::vector<Cache_YH2O> cache_YH2O;
  std::vector<Cache_XC> cache_XC;
  // std::vector<Cache_pen> cache_pen;

  // constructor stores parameters
  LocalOperator(const GV &gv_, const Params&	 property_,
                     U *unew_,
                     GFS gfs_,
                     double *time_,
                     double *dt_,
                     unsigned int intorder_ = 6,
                     double method_g_ = 0.,
                     double method_w_ = 0.,
                     double method_T_ = 0.,
                     double method_x_ = 0.,
                     double method_y_ = 0.,
                     double alpha_g_ = 1., double alpha_w_ = 1., double alpha_s_ = 1., double alpha_T_ = 1., double alpha_x_ = 1., double alpha_y_ = 1.)
      : gv(gv_), property( property_ ),
        unew(unew_),
        gfs(gfs_),
        time(time_),
        dt(dt_),
        intorder(intorder_),
        method_g(method_g_), method_w(method_w_), method_T(method_T_), method_x(method_x_), method_y(method_y_),
        alpha_g(alpha_g_), alpha_w(alpha_w_), alpha_s(alpha_s_),  alpha_T(alpha_T_), alpha_x(alpha_x_), alpha_y(alpha_y_),
        cache_Pw(20), cache_Sg(20), cache_Sh(20),  cache_T(20), cache_XCH4(20), cache_YH2O(20), cache_XC(20)
  {
    theta_g = method_g;
    // if (method_g == -1.)
    //   theta_g = -1.0;
    // if (method_g == 1.)
    //   theta_g = 1.0;

    theta_w = method_w;
    // if (method_w == -1.)
    //   theta_w = -1.0;
    // if (method_w == 1.)
    //   theta_w = 1.0;

    theta_T = method_T;
    // if (method_T == -1.)
    //   theta_T = -1.0;
    // if (method_T == 1.)
    //   theta_T = 1.0;

    theta_x = method_x;
    // if (method_x == ConvectionDiffusionDGMethod::SIPG)
    //   theta_x = -1.0;
    // if (method_w == ConvectionDiffusionDGMethod::NIPG)
    //   theta_x = 1.0;

    theta_y = method_y;
    // if (method_y == ConvectionDiffusionDGMethod::SIPG)
    //   theta_y = -1.0;
    // if (method_y == ConvectionDiffusionDGMethod::NIPG)
    //   theta_y = 1.0;

    Xc_conv_m = property.characteristicValue.X_convective_mass;
    Xc_conv_h = property.characteristicValue.X_convective_heat;
    Xc_source_m = property.characteristicValue.X_source_mass;
    Xc_source_h = property.characteristicValue.X_source_heat;
    Xc_diff_m = property.characteristicValue.X_diffusive_mass;
    Xc_diff_h = property.characteristicValue.X_diffusive_heat;
    Xc_grav = property.characteristicValue.X_gravity;

    Xc_K = property.characteristicValue.permeability_c;
    Xc_mu = property.characteristicValue.viscosity_c;
    Xc_rho = property.characteristicValue.density_c;
    Xc_kth = property.characteristicValue.thermalconductivity_c;
    Xc_C = property.characteristicValue.specificheat_c;
    Xc_v = property.characteristicValue.volumetricheat_c;
    Xc_P = property.characteristicValue.P_c;
    Xc_T = property.characteristicValue.T_c;
    Xc_t = property.characteristicValue.t_c;
    Xc_x = property.characteristicValue.x_c;
    T_ref = property.parameter.ReferenceTemperature()/Xc_T;
    gravity = property.parameter.g() / Xc_grav  ; /* ndim */

  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
  {
    // subspaces
    //water pressure
    const auto &lfsv_Pw = lfsv.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw = lfsu.template child<Indices::PVId_Pw>();

    //gas Saturation
    const auto &lfsv_Sg = lfsv.template child<Indices::PVId_Sg>();
    const auto &lfsu_Sg = lfsu.template child<Indices::PVId_Sg>();

    //Hydrate Saturation
    const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

    //Temperature
    const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

    //Methane (diss.) mole fraction
    const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    //H2O (vap.) mole fraction
    const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();

    //Salt (diss.) mole fraction
    const auto &lfsv_XC = lfsv.template child<Indices::PVId_C>();
    const auto &lfsu_XC = lfsu.template child<Indices::PVId_C>();



    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    typedef typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;

    // dimensions
    const int dim = EG::Entity::dimension;
    const int order_p = std::max(lfsu_Pw.finiteElement().localBasis().order(),
                               lfsv_Pw.finiteElement().localBasis().order());/* If different degrees are used for different functions ? */
    const int order_x = std::max(lfsu_XC.finiteElement().localBasis().order(),
                               lfsv_XC.finiteElement().localBasis().order());
    const int order_s = std::max(lfsu_Sg.finiteElement().localBasis().order(),
                               lfsv_Sg.finiteElement().localBasis().order());
    const int order_t = std::max(lfsu_T.finiteElement().localBasis().order(),
                               lfsv_T.finiteElement().localBasis().order());

    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pw(lfsu_Pw.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pw(lfsv_Pw.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sg(lfsu_Sg.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sg(lfsv_Sg.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sh(lfsu_Sh.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sh(lfsv_Sh.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T(lfsv_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4(lfsu_XCH4.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4(lfsv_XCH4.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O(lfsu_YH2O.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O(lfsv_YH2O.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC(lfsu_XC.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC(lfsv_XC.size());

    Dune::FieldVector<RF, dim> gradu_Pw(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh(0.0);
    Dune::FieldVector<RF, dim> gradu_T(0.0);
    Dune::FieldVector<RF, dim> Ktgradu_T(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O(0.0);
    Dune::FieldVector<RF, dim> gradu_XC(0.0);
    Dune::FieldVector<RF, dim> Kg(0.0);

    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_Pw = cache_Pw[order_p].evaluateFunction(ip.position(), lfsu_Pw.finiteElement().localBasis());
      auto &psi_Pw = cache_Pw[order_p].evaluateFunction(ip.position(), lfsv_Pw.finiteElement().localBasis());
      auto &phi_Sg = cache_Sg[order_s].evaluateFunction(ip.position(), lfsu_Sg.finiteElement().localBasis());
      auto &psi_Sg = cache_Sg[order_s].evaluateFunction(ip.position(), lfsv_Sg.finiteElement().localBasis());
      auto &phi_Sh = cache_Sh[order_s].evaluateFunction(ip.position(), lfsu_Sh.finiteElement().localBasis());
      auto &psi_Sh = cache_Sh[order_s].evaluateFunction(ip.position(), lfsv_Sh.finiteElement().localBasis());
      auto &phi_T  = cache_T[order_t].evaluateFunction(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &psi_T  = cache_T[order_t].evaluateFunction(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &phi_XCH4 = cache_XCH4[order_x].evaluateFunction(ip.position(), lfsu_XCH4.finiteElement().localBasis());
      auto &psi_XCH4 = cache_XCH4[order_x].evaluateFunction(ip.position(), lfsv_XCH4.finiteElement().localBasis());
      auto &phi_YH2O = cache_YH2O[order_x].evaluateFunction(ip.position(), lfsu_YH2O.finiteElement().localBasis());
      auto &psi_YH2O = cache_YH2O[order_x].evaluateFunction(ip.position(), lfsv_YH2O.finiteElement().localBasis());
      auto &phi_XC = cache_XC[order_x].evaluateFunction(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &psi_XC = cache_XC[order_x].evaluateFunction(ip.position(), lfsv_XC.finiteElement().localBasis());

      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);
      // evaluate Pw
      RF Pw = 0.0;
      for (size_type i = 0; i < lfsu_Pw.size(); i++)
        Pw += x(lfsu_Pw, i) * phi_Pw[i];

      // evaluate Sg
      RF Sg = 0.0;
      for (size_type i = 0; i < lfsu_Sg.size(); i++){
        Sg += x(lfsu_Sg, i) * phi_Sg[i];
      }

      // evaluate Sh
      RF Sh = 0.0;
      for (size_type i = 0; i < lfsu_Sh.size(); i++){
        Sh += x(lfsu_Sh, i) * phi_Sh[i];
      }

      // evaluate T
      RF T = 0.0;
      for (size_type i = 0; i < lfsu_T.size(); i++)
        T += x(lfsu_T, i) * phi_T[i];

      // evaluate XCH4
      RF XCH4 = 0.0;
      for (size_type i = 0; i < lfsu_XCH4.size(); i++)
        XCH4 += x(lfsu_XCH4, i) * phi_XCH4[i];
      // if (XCH4 < 0.)
      //   XCH4 = 0.;

      // evaluate YH2O
      RF YH2O = 0.0;
      for (size_type i = 0; i < lfsu_YH2O.size(); i++)
        YH2O += x(lfsu_YH2O, i) * phi_YH2O[i];

      // evaluate XC
      RF XC = 0.0;
      for (size_type i = 0; i < lfsu_XC.size(); i++)
        XC += x(lfsu_XC, i) * phi_XC[i];

      RF Sw = (1. - Sg - Sh);


      // evaluate Pg
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por = property.soil.SedimentPorosity(cell, ip_local) ;//* (1.-Sh)
      auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      RF Pg = Pw + Pc; /* ndim */
      RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh); /* ndim */

      auto Pw_dim = Pw * Xc_P;
      auto Pg_dim = Pg * Xc_P;
      auto T_dim = T * Xc_T;

      // evaluate gradient of basis functions
      auto &js_Pw = cache_Pw[order_p].evaluateJacobian(ip.position(), lfsu_Pw.finiteElement().localBasis());
      auto &js_v_Pw = cache_Pw[order_p].evaluateJacobian(ip.position(), lfsv_Pw.finiteElement().localBasis());
      auto &js_Sg = cache_Sg[order_s].evaluateJacobian(ip.position(), lfsu_Sg.finiteElement().localBasis());
      auto &js_v_Sg = cache_Sg[order_s].evaluateJacobian(ip.position(), lfsv_Sg.finiteElement().localBasis());
      auto &js_Sh = cache_Sh[order_s].evaluateJacobian(ip.position(), lfsu_Sh.finiteElement().localBasis());
      auto &js_v_Sh = cache_Sh[order_s].evaluateJacobian(ip.position(), lfsv_Sh.finiteElement().localBasis());
      auto &js_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &js_XCH4 = cache_XCH4[order_x].evaluateJacobian(ip.position(), lfsu_XCH4.finiteElement().localBasis());
      auto &js_v_XCH4 = cache_XCH4[order_x].evaluateJacobian(ip.position(), lfsv_XCH4.finiteElement().localBasis());
      auto &js_YH2O = cache_YH2O[order_x].evaluateJacobian(ip.position(), lfsu_YH2O.finiteElement().localBasis());
      auto &js_v_YH2O = cache_YH2O[order_x].evaluateJacobian(ip.position(), lfsv_YH2O.finiteElement().localBasis());
      auto &js_XC = cache_XC[order_x].evaluateJacobian(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &js_v_XC = cache_XC[order_x].evaluateJacobian(ip.position(), lfsv_XC.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (size_type i = 0; i < lfsu_Pw.size(); i++)
        jac.mv(js_Pw[i][0], gradphi_Pw[i]);
      for (size_type i = 0; i < lfsv_Pw.size(); i++)
        jac.mv(js_v_Pw[i][0], gradpsi_Pw[i]);


      for (size_type i = 0; i < lfsu_Sg.size(); i++)
        jac.mv(js_Sg[i][0], gradphi_Sg[i]);
      for (size_type i = 0; i < lfsv_Sg.size(); i++)
        jac.mv(js_v_Sg[i][0], gradpsi_Sg[i]);

      for (size_type i = 0; i < lfsu_Sh.size(); i++)
        jac.mv(js_Sh[i][0], gradphi_Sh[i]);
      for (size_type i = 0; i < lfsv_Sh.size(); i++)
        jac.mv(js_v_Sh[i][0], gradpsi_Sh[i]);

      for (size_type i = 0; i < lfsu_T.size(); i++)
        jac.mv(js_T[i][0], gradphi_T[i]);
      for (size_type i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_T[i][0], gradpsi_T[i]);

      for (size_type i = 0; i < lfsu_XCH4.size(); i++)
        jac.mv(js_XCH4[i][0], gradphi_XCH4[i]);
      for (size_type i = 0; i < lfsv_XCH4.size(); i++)
        jac.mv(js_v_XCH4[i][0], gradpsi_XCH4[i]);

      for (size_type i = 0; i < lfsu_YH2O.size(); i++)
        jac.mv(js_YH2O[i][0], gradphi_YH2O[i]);
      for (size_type i = 0; i < lfsv_YH2O.size(); i++)
        jac.mv(js_v_YH2O[i][0], gradpsi_YH2O[i]);


      for (size_type i = 0; i < lfsu_XC.size(); i++)
        jac.mv(js_XC[i][0], gradphi_XC[i]);
      for (size_type i = 0; i < lfsv_XC.size(); i++)
        jac.mv(js_v_XC[i][0], gradpsi_XC[i]);

      // compute gradient of Pw
      gradu_Pw = 0.0;
      for (size_type i = 0; i < lfsu_Pw.size(); i++)
        gradu_Pw.axpy(x(lfsu_Pw, i), gradphi_Pw[i]);

      // compute gradient of Sg
      gradu_Sg = 0.0;
      for (size_type i = 0; i < lfsu_Sg.size(); i++)
        gradu_Sg.axpy(x(lfsu_Sg, i), gradphi_Sg[i]);

      // compute gradient of Sh
      gradu_Sh = 0.0;
      for (size_type i = 0; i < lfsu_Sh.size(); i++)
        gradu_Sh.axpy(x(lfsu_Sh, i), gradphi_Sh[i]);

      // compute gradient of T
      gradu_T = 0.0;
      for (size_type i = 0; i < lfsu_T.size(); i++)
        gradu_T.axpy(x(lfsu_T, i), gradphi_T[i]);

      // compute gradient of XCH4
      gradu_XCH4 = 0.0;
      for (size_type i = 0; i < lfsu_XCH4.size(); i++)
        gradu_XCH4.axpy(x(lfsu_XCH4, i), gradphi_XCH4[i]);

      // compute gradient of YH2O
      gradu_YH2O = 0.0;
      for (size_type i = 0; i < lfsu_YH2O.size(); i++)
        gradu_YH2O.axpy(x(lfsu_YH2O, i), gradphi_YH2O[i]);

      // compute gradient of XCH4
      gradu_XC = 0.0;
      for (size_type i = 0; i < lfsu_XC.size(); i++)
        gradu_XC.axpy(x(lfsu_XC, i), gradphi_XC[i]);


      auto K = property.soil.SedimentPermeabilityTensor(cell, ip_local)
                    * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por ); /*ndim K from soil.hh*/
      K.mv(gravity, Kg);

      // compute K * gradient of Pw
      K.mv(gradu_Pw, Kgradu_Pw);

      // compute K * gradient of Sg
      K.mv(gradu_Sg, Kgradu_Sg);

      // compute K * gradient of Sh
      K.mv(gradu_Sh, Kgradu_Sh);

      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );

      auto tau = property.soil.Tortuosity(por);/*ndim tau from soil.hh*/
      auto DH2O_g = tau * por * property.mixture.DiffCoeffH2OInGas(T_dim, Pg_dim); /*ndim D from mixture.hh*/
      auto DCH4_w = tau * por * property.mixture.DiffCoeffCH4InLiquid(T_dim, Pw_dim); /*ndim D from mixture.hh*/
      auto DC_w = tau * por * property.salt.DiffCoeff(T_dim, Pw_dim); /*ndim D from salt.hh*/

      double S = XC * (property.salt.MolarMass()/property.gas.MolarMass());
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T_dim, Pg_dim);
      auto H_CH4_w = property.gas.SolubilityCoefficient(  T_dim/*K*/, S ); /*ndim */
      auto P_H_satu = property.water.SaturatedVaporPressure( T_dim /*K*/, S ); /*ndim */

      auto YCH4 =  property.mixture.YCH4(XCH4, T_dim, Pg_dim, XC, zCH4);
      auto XH2O =  property.mixture.XH2O(YH2O, T_dim, Pg_dim, XC);

      double eta = 1/BrooksCParams[1];
      auto Swe = property.hydraulicProperty.EffectiveSw(Sw,Sh,0., 0.);
      auto dPc_dSwe = property.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]); /*ndim */
      auto dSwe_dSw =  property.hydraulicProperty.dSwe_dSw(Sw,Sh,0., 0.);
      auto coeff_grad_Sw = dPc_dSwe * dSwe_dSw * property.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]);

      auto dPcSF1_dSh =  property.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh = property.hydraulicProperty.dSwe_dSh(Sw,Sh,0., 0.);
      auto coeff_grad_Sh = dPcSF1_dSh * std::pow( Swe , -eta ) * BrooksCParams[0] / Xc_P - Sg * coeff_grad_Sw *  dSwe_dSw  ;

      auto rho_g = property.gas.Density(T_dim, Pg_dim, zCH4); /*ndim density from CH4.hh; the input arguments are dimensional   */
      auto rho_w = property.water.Density(T_dim, Pw_dim, S); /*ndim density from H2O.hh; the input arguments are dimensional*/

      auto krW = property.hydraulicProperty.krw(cell, ip_local, Sw, Sh) / (property.water.DynamicViscosity(T_dim, Pw_dim, S) );
      auto krN = property.hydraulicProperty.krg(cell, ip_local, Sw, Sh) / (property.gas.DynamicViscosity(T_dim, Pg_dim));

      auto perm_dim = property.soil.SedimentPermeability(cell, ip_local ) * Xc_K;
      // compute source terms
			auto q_g  = property.kinetics.GasGenerationRate( T_dim, Pg_dim, Sh, Sw, XCH4, zCH4, S, por, perm_dim); /*dim*/
			auto q_w  = property.kinetics.WaterGenerationRate( q_g ); /*dim*/
			auto q_h  = property.kinetics.HydrateDissociationRate( q_g ); /*dim*/
			auto q_s = property.salt.Source(); /*dim------kg/m³s*/
			auto Q = property.kinetics.HeatOfDissociation( q_g, T_dim ); /*dim----[W/m³]*/

      auto Cp_g = property.gas.Cp(T_dim, Pg_dim, zCH4); /* ndim */
      auto Cp_w = property.water.Cp(T_dim, Pw_dim, S); /* ndim */
      auto kth_g = property.gas.ThermalConductivity(T_dim, Pg_dim); /* ndim */
      auto kth_w = property.water.ThermalConductivity(T_dim, Pw_dim, S); /* ndim */
      auto kth_h = property.hydrate.ThermalConductivity(T_dim, Peff * Xc_P); /* ndim */
      auto kth_s = property.soil.ThermalConductivity(); /* ndim */
      auto kth_eff = (1. - por) * kth_s + por * (Sg * kth_g + Sw * kth_w + Sh * kth_h); /* ndim */

      auto gradu_Pg = gradu_Pw  - coeff_grad_Sw * gradu_Sg + (coeff_grad_Sh ) * gradu_Sh;
      auto Kgradu_Pg = Kgradu_Pw - coeff_grad_Sw * Kgradu_Sg + (coeff_grad_Sh ) * Kgradu_Sh;
      
      auto convectiveflux_CH4_g = rho_g * (1. - YH2O) * krN * ( Xc_P/Xc_x * Kgradu_Pg - Xc_rho * rho_g * Kg); //(1. - YH2O)
      auto convectiveflux_CH4_w = rho_w * (XCH4) * krW * ( Xc_P/Xc_x * Kgradu_Pw - Xc_rho * rho_w * Kg);
      auto convectiveflux_H2O_g = rho_g * YH2O * krN * ( Xc_P/Xc_x * Kgradu_Pg - Xc_rho * rho_g * Kg);
      auto convectiveflux_H2O_w = rho_w * (1. - XC - XCH4) * krW * ( Xc_P/Xc_x * Kgradu_Pw - Xc_rho * rho_w * Kg);
      auto convectiveflux_SALT_w = rho_w * (XC) * krW * ( Xc_P/Xc_x * Kgradu_Pw - Xc_rho * rho_w * Kg);
      auto convectiveflux_Heat_w = rho_w * Cp_w * (T - T_ref) * krW * ( Xc_P/Xc_x * Kgradu_Pw - Xc_rho * rho_w * Kg);
      auto convectiveflux_Heat_g = rho_g * Cp_g * (T - T_ref) * krN *( Xc_P/Xc_x * Kgradu_Pg - Xc_rho * rho_g * Kg);

      auto j_H2O_g = rho_g * Sg * DH2O_g * gradu_YH2O;
      auto j_CH4_w = rho_w * Sw * DCH4_w * gradu_XCH4;
      auto j_SALT_w = rho_w * Sw * DC_w * gradu_XC;
      auto j_H2O_w = -j_CH4_w -j_SALT_w;
      auto j_CH4_g = -j_H2O_g;

      auto convectiveflux_CH4 = convectiveflux_CH4_g + convectiveflux_CH4_w;
      auto convectiveflux_H2O = convectiveflux_H2O_g + convectiveflux_H2O_w;
      auto convectiveflux_Heat = convectiveflux_Heat_g + convectiveflux_Heat_w;

      auto diffusiveflux_CH4 = j_CH4_g + j_CH4_w;
      auto diffusiveflux_H2O = j_H2O_g + j_H2O_w;
      auto diffusiveflux_SALT = j_SALT_w;
      auto diffusiveflux_Heat = kth_eff * gradu_T;

      
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      double tmp = 0.;
      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      for (size_type i = 0; i < lfsv_Sg.size(); i++)
      {
        r.accumulate(lfsv_Sg, i, (( Xc_conv_m * convectiveflux_CH4
                                  + Xc_diff_m * diffusiveflux_CH4 ) * gradpsi_Pw[i]
                                  - q_g *Xc_source_m * psi_Pw[i]) * factor);
      }

      for (size_type i = 0; i < lfsv_XC.size(); i++)
      {
        r.accumulate(lfsv_XC, i, (( Xc_conv_m * convectiveflux_SALT_w
                                  + Xc_diff_m * diffusiveflux_SALT ) * gradpsi_Pw[i]
                                  ) * factor);
      }

      for (size_type i = 0; i < lfsv_Pw.size(); i++)
      {
        r.accumulate(lfsv_Pw, i, (( Xc_conv_m * convectiveflux_H2O
                                  + Xc_diff_m * diffusiveflux_H2O ) * gradpsi_Pw[i]
                                  - (q_w ) * Xc_source_m * psi_Pw[i]) * factor);
      }

      for (size_type i = 0; i < lfsv_Sh.size(); i++)
      {
        r.accumulate(lfsv_Sh, i, (-q_h * Xc_source_m * psi_Sh[i]) * factor);
      }

      for (size_type i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, ((Xc_conv_h * convectiveflux_Heat +  Xc_diff_h * diffusiveflux_Heat ) * gradpsi_T[i]
                                  - Q * Xc_source_h * psi_T[i]) * factor);
      }

      //Integrals regarding the NCP1
			// RF max1 = std::max(0., (Sg -1. + YCH4 + YH2O));
			// for (size_type i=0; i<lfsv_YH2O.size(); i++){
			// 	r.accumulate(lfsv_YH2O,i,( (Sg - max1) * psi_YH2O[i]  *factor));
			// }

			// // Integrals regarding the NCP2
			// RF max2 = std::max(0., (Sw -1. + XC + XCH4 + XH2O ));
			// for (size_type i=0; i<lfsv_XCH4.size(); i++){
			// 	r.accumulate(lfsv_XCH4,i,((Sw - max2) * psi_XCH4[i]  *factor));
			// }

      
      // for (size_type i = 0; i < lfsv_XCH4.size(); i++)
      // {
      //   auto &phi_Pw0 = cache_Pw[order_p].evaluateFunction(geo.local(geo.corner(i)), lfsu_Pw.finiteElement().localBasis());
      //   auto &psi_Pw0 = cache_Pw[order_p].evaluateFunction(geo.local(geo.corner(i)), lfsv_Pw.finiteElement().localBasis());

      //   auto &phi_Sg0 = cache_Sg[order_s].evaluateFunction(geo.local(geo.corner(i)), lfsu_Sg.finiteElement().localBasis());
      //   auto &psi_Sg0 = cache_Sg[order_s].evaluateFunction(geo.local(geo.corner(i)), lfsv_Sg.finiteElement().localBasis());
      //   auto &phi_Sh0 = cache_Sh[order_s].evaluateFunction(geo.local(geo.corner(i)), lfsu_Sh.finiteElement().localBasis());
      //   auto &psi_Sh0 = cache_Sh[order_s].evaluateFunction(geo.local(geo.corner(i)), lfsv_Sh.finiteElement().localBasis());
      //   auto &phi_T0  = cache_T[order_t].evaluateFunction(geo.local(geo.corner(i)), lfsu_T.finiteElement().localBasis());
      //   auto &psi_T0  = cache_T[order_t].evaluateFunction(geo.local(geo.corner(i)), lfsv_T.finiteElement().localBasis());
      //   auto &phi_XCH40 = cache_XCH4[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsu_XCH4.finiteElement().localBasis());
      //   auto &psi_XCH40 = cache_XCH4[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsv_XCH4.finiteElement().localBasis());
      //   auto &phi_YH2O0 = cache_YH2O[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsu_YH2O.finiteElement().localBasis());
      //   auto &psi_YH2O0 = cache_YH2O[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsv_YH2O.finiteElement().localBasis());
      //   auto &phi_XC0 = cache_XC[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsu_XC.finiteElement().localBasis());
      //   auto &psi_XC0 = cache_XC[order_x].evaluateFunction(geo.local(geo.corner(i)), lfsv_XC.finiteElement().localBasis());
      //   // evaluate Pw0
      // RF Pw0 = 0.0;
      //   for (size_type i = 0; i < lfsu_Pw.size(); i++)
      //   Pw0 += x(lfsu_Pw, i) * phi_Pw0[i];

      // // evaluate Sg
      // RF Sg0 = 0.0;
      // for (size_type i = 0; i < lfsu_Sg.size(); i++){
      //   Sg0 += x(lfsu_Sg, i) * phi_Sg0[i];
      // }

      // // evaluate Sh
      // RF Sh0 = 0.0;
      // for (size_type i = 0; i < lfsu_Sh.size(); i++){
      //   Sh0 += x(lfsu_Sh, i) * phi_Sh0[i];
      // }

      // // evaluate T
      // RF T0 = 0.0;
      // for (size_type i = 0; i < lfsu_T.size(); i++)
      //   T0 += x(lfsu_T, i) * phi_T0[i];

      // // evaluate XCH4
      // RF XCH40 = 0.0;
      // for (size_type i = 0; i < lfsu_XCH4.size(); i++)
      //   XCH40 += x(lfsu_XCH4, i) * phi_XCH40[i];

      // // evaluate YH2O
      // RF YH2O0 = 0.0;
      // for (size_type i = 0; i < lfsu_YH2O.size(); i++)
      //   YH2O0 += x(lfsu_YH2O, i) * phi_YH2O0[i];

      // // evaluate XC
      // RF XC0 = 0.0;
      // for (size_type i = 0; i < lfsu_XC.size(); i++)
      //   XC0 += x(lfsu_XC, i) * phi_XC0[i];

      //   auto Sw0 = 1. - Sh0 - Sg0;
      //   auto Pc0 = property.hydraulicProperty.CapillaryPressure(cell, geo.local(geo.corner(i)), Sw0, Sh0, por) ; /* ndim */
      //   auto zCH40 = property.eos.EvaluateCompressibilityFactor(T0* Xc_T, (Pw0+Pc0)* Xc_P);
      //   auto YCH40 =  property.mixture.YCH4(XCH40, T0* Xc_T, (Pw0+Pc0)* Xc_P, XC0, zCH40);
      //   auto XH2O0 =  property.mixture.XH2O(YH2O0, T0* Xc_T, (Pw0+Pc0)* Xc_P, XC0);

      //   if(Sw0 - (1. - XC0 - XH2O0 - XCH40)> 0.)
      //      tmp = 1. - XC0 - XH2O0 - XCH40;
      //   else
      //      tmp= 1. - Sh0 - Sg0;
      //   r.accumulate(lfsv_XCH4 , i, +tmp * psi_XCH40[i] * geo.integrationElement(geo.local(geo.corner(i))));

      //   tmp = 0.;
      //   if(Sg0 -(1. - YH2O0 - YCH40)> 0.)
      //      tmp = 1. - YH2O0 - YCH40 ;
      //   else
      //      tmp= Sg0;
      //   r.accumulate(lfsv_YH2O , i, +tmp * psi_YH2O0[i] * geo.integrationElement(geo.local(geo.corner(i))));
      //   //exit(0);
      // }

      // NCP -> water phase
			tmp = 0.;
			//auto XH2O_alg = property.mixture.XH2O(YH2O,T_dim,Pg_dim,S);
			if( Sw - ( 1. - XCH4 - XH2O - XC ) > eps_ap ){//active set.std::max(0., std::min(1. ,  )
				tmp +=  1. - XCH4 - XH2O - XC ;//Active => phase is present => summation condition holds
        //	std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			}else {//std::max(0., std::min(1. , )if ( Sw - ( 1. - XCH4 - XH2O - XC )  < eps_ap )
				tmp += Sw; // inactive set. Inactive => phase is absent => Sw=0
			}
      for (size_type i = 0; i < lfsv_XCH4.size(); i++)
      {
			r.accumulate(lfsv_XCH4 , i, +tmp * psi_XCH4[i]  *factor);
      }

			// NCP -> gas phase
			tmp = 0.;
			//auto YCH4_alg = property.mixture.YCH4(XCH4,T_dim,Pg_dim,S,zCH4);
			if(  Sg - ( 1. - YCH4 - YH2O )  > eps_ap  ){ //active set.std::max(0., std::min(1. , )
				tmp +=   1. - YCH4 - YH2O;//Active => phase is present => summation condition holds
        //std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			}else {//std::max(0., std::min(1. ,  )if( Sg -( 1. - YCH4 - YH2O )  < eps_ap )
				tmp += Sg;// inactive set. Inactive => phase is absent => Sg=0
			}
      for (size_type i = 0; i < lfsv_YH2O.size(); i++)
      {
			r.accumulate(lfsv_YH2O , i, +tmp * psi_YH2O[i]  *factor);
      }




    } //End Quadrature Rule
  }  // End of alpha_volume



  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG &ig,
                      const LFSU &lfsu_s, const X &x_s, const LFSV &lfsv_s,
                      const LFSU &lfsu_n, const X &x_n, const LFSV &lfsv_n,
                      R &r_s, R &r_n) const
  {
    // subspaces
    //water pressure
    const auto &lfsv_Pw_s = lfsv_s.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_s = lfsu_s.template child<Indices::PVId_Pw>();
    const auto &lfsv_Pw_n = lfsv_n.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_n = lfsu_n.template child<Indices::PVId_Pw>();

    //gas Saturation
    const auto &lfsv_Sg_s = lfsv_s.template child<Indices::PVId_Sg>();
    const auto &lfsu_Sg_s = lfsu_s.template child<Indices::PVId_Sg>();
    const auto &lfsv_Sg_n = lfsv_n.template child<Indices::PVId_Sg>();
    const auto &lfsu_Sg_n = lfsu_n.template child<Indices::PVId_Sg>();

    //Hydrate Saturation
    const auto &lfsv_Sh_s = lfsv_s.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh_s = lfsu_s.template child<Indices::PVId_Sh>();
    const auto &lfsv_Sh_n = lfsv_n.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh_n = lfsu_n.template child<Indices::PVId_Sh>();

    //Temperature
    const auto &lfsv_T_s = lfsv_s.template child<Indices::PVId_T>();
    const auto &lfsu_T_s = lfsu_s.template child<Indices::PVId_T>();
    const auto &lfsv_T_n = lfsv_n.template child<Indices::PVId_T>();
    const auto &lfsu_T_n = lfsu_n.template child<Indices::PVId_T>();

    //Methane mole fraction
    const auto &lfsv_XCH4_s = lfsv_s.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4_s = lfsu_s.template child<Indices::PVId_XCH4>();
    const auto &lfsv_XCH4_n = lfsv_n.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4_n = lfsu_n.template child<Indices::PVId_XCH4>();

    //Water mole fraction
    const auto &lfsv_YH2O_s = lfsv_s.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O_s = lfsu_s.template child<Indices::PVId_YH2O>();
    const auto &lfsv_YH2O_n = lfsv_n.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O_n = lfsu_n.template child<Indices::PVId_YH2O>();

    //Salt mole fraction
    const auto &lfsv_XC_s = lfsv_s.template child<Indices::PVId_C>();
    const auto &lfsu_XC_s = lfsu_s.template child<Indices::PVId_C>();
    const auto &lfsv_XC_n = lfsv_n.template child<Indices::PVId_C>();
    const auto &lfsu_XC_n = lfsu_n.template child<Indices::PVId_C>();


    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;

    // dimensions
    const int dim= IG::Entity::dimension; // = GV::dimension
    const int dimension = GV::dimension;

    const int order_p = std::max(lfsu_Pw_s.finiteElement().localBasis().order(),
                                lfsv_Pw_s.finiteElement().localBasis().order());/* If different degrees are used for different functions ? */
    const int order_x = std::max(lfsu_XC_s.finiteElement().localBasis().order(),
                               lfsv_XC_s.finiteElement().localBasis().order());
    const int order_s = std::max(lfsu_Sg_s.finiteElement().localBasis().order(),
                               lfsv_Sg_s.finiteElement().localBasis().order());
    const int order_t = std::max(lfsu_T_s.finiteElement().localBasis().order(),
                               lfsv_T_s.finiteElement().localBasis().order());


    // References to inside and outside cells
    const auto &cell_inside = ig.inside();
    const auto &cell_outside = ig.outside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();
    auto geo_outside = cell_outside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();
    auto geo_in_outside = ig.geometryInOutside();
    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
	  auto ref_el_outside = referenceElement(geo_outside);
	  auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	  auto outside_cell_center_local 	= ref_el_outside.position(0,0);
    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = std::min(geo_inside.volume(), geo_outside.volume()) / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    omega_s = 0.5;
    omega_n = 0.5;
    harmonic_average = 1.0;

        //  if (true)
        //    {
        //      RF delta_s = (An_F_s*n_F);
        //      RF delta_n = (An_F_n*n_F);
        //      omega_s = delta_n/(delta_s+delta_n+1e-20);
        //      omega_n = delta_s/(delta_s+delta_n+1e-20);
        //      harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
        //    }
        //  else
        //    {
        //      omega_s = omega_n = 0.5;
        //      harmonic_average = 1.0;
        //    }

    // get polynomial degree
    auto order_i = lfsv_Pw_s.finiteElement().localBasis().order();
    auto order_o = lfsv_Pw_n.finiteElement().localBasis().order();
    auto degree = std::max(order_i, order_o);

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pw_s(lfsu_Pw_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pw_s(lfsv_Pw_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sg_s(lfsu_Sg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sg_s(lfsv_Sg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sh_s(lfsu_Sh_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sh_s(lfsv_Sh_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_s(lfsu_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_s(lfsv_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_s(lfsu_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_s(lfsv_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_s(lfsu_XC_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_s(lfsv_XC_s.size());

    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pw_n(lfsu_Pw_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pw_n(lfsv_Pw_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sg_n(lfsu_Sg_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sg_n(lfsv_Sg_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sh_n(lfsu_Sh_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sh_n(lfsv_Sh_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_n(lfsu_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_n(lfsv_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_n(lfsu_XCH4_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_n(lfsv_XCH4_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_n(lfsu_YH2O_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_n(lfsv_YH2O_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_n(lfsu_XC_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_n(lfsv_XC_n.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> gradu_T_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_n(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_n(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);

    Dune::FieldVector<RF, dim> Kg_s(0.0);
    Dune::FieldVector<RF, dim> Kg_n(0.0);

    
    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    auto intorder1 = intorder;
    if( dim == 1 )
      intorder1 = 1;
    for (const auto &ip : quadratureRule(geo, intorder1))
    {
      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto iplocal_n = geo_in_outside.global(ip.position());

      auto ip_global_s = geo_inside.global(iplocal_s);
      auto ip_global_n = geo_outside.global(iplocal_n);

      // evaluate basis functions
      auto &phi_Pw_s = cache_Pw[order_p].evaluateFunction(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &psi_Pw_s = cache_Pw[order_p].evaluateFunction(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &phi_Sg_s = cache_Sg[order_s].evaluateFunction(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &psi_Sg_s = cache_Sg[order_s].evaluateFunction(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &phi_Sh_s = cache_Sh[order_s].evaluateFunction(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &psi_Sh_s = cache_Sh[order_s].evaluateFunction(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &phi_XCH4_s = cache_XCH4[order_x].evaluateFunction(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &psi_XCH4_s = cache_XCH4[order_x].evaluateFunction(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &phi_YH2O_s = cache_YH2O[order_x].evaluateFunction(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &psi_YH2O_s = cache_YH2O[order_x].evaluateFunction(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &phi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &psi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      auto &phi_Pw_n = cache_Pw[order_p].evaluateFunction(iplocal_n, lfsu_Pw_n.finiteElement().localBasis());
      auto &psi_Pw_n = cache_Pw[order_p].evaluateFunction(iplocal_n, lfsv_Pw_n.finiteElement().localBasis());
      auto &phi_Sg_n = cache_Sg[order_s].evaluateFunction(iplocal_n, lfsu_Sg_n.finiteElement().localBasis());
      auto &psi_Sg_n = cache_Sg[order_s].evaluateFunction(iplocal_n, lfsv_Sg_n.finiteElement().localBasis());
      auto &phi_Sh_n = cache_Sh[order_s].evaluateFunction(iplocal_n, lfsu_Sh_n.finiteElement().localBasis());
      auto &psi_Sh_n = cache_Sh[order_s].evaluateFunction(iplocal_n, lfsv_Sh_n.finiteElement().localBasis());
      auto &phi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &psi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      auto &phi_XCH4_n = cache_XCH4[order_x].evaluateFunction(iplocal_n, lfsu_XCH4_n.finiteElement().localBasis());
      auto &psi_XCH4_n = cache_XCH4[order_x].evaluateFunction(iplocal_n, lfsv_XCH4_n.finiteElement().localBasis());
      auto &phi_YH2O_n = cache_YH2O[order_x].evaluateFunction(iplocal_n, lfsu_YH2O_n.finiteElement().localBasis());
      auto &psi_YH2O_n = cache_YH2O[order_x].evaluateFunction(iplocal_n, lfsv_YH2O_n.finiteElement().localBasis());
      auto &phi_XC_n = cache_XC[order_x].evaluateFunction(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &psi_XC_n = cache_XC[order_x].evaluateFunction(iplocal_n, lfsv_XC_n.finiteElement().localBasis());



      // evaluate Pw
      RF Pw_s = 0.0;
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        Pw_s += x_s(lfsu_Pw_s, i) * phi_Pw_s[i];
      RF Pw_n = 0.0;
      for (size_type i = 0; i < lfsu_Pw_n.size(); i++)
        Pw_n += x_n(lfsu_Pw_n, i) * phi_Pw_n[i];

      // evaluate Sg
      RF Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        Sg_s += x_s(lfsu_Sg_s, i) * phi_Sg_s[i];
      RF Sg_n = 0.0;
      for (size_type i = 0; i < lfsu_Sg_n.size(); i++)
        Sg_n += x_n(lfsu_Sg_n, i) * phi_Sg_n[i];

      // evaluate Sh
      RF Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        Sh_s += x_s(lfsu_Sh_s, i) * phi_Sh_s[i];
      RF Sh_n = 0.0;
      for (size_type i = 0; i < lfsu_Sh_n.size(); i++)
        Sh_n += x_n(lfsu_Sh_n, i) * phi_Sh_n[i];

      // evaluate T
      RF T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        T_s += x_s(lfsu_T_s, i) * phi_T_s[i];
      RF T_n = 0.0;
      for (size_type i = 0; i < lfsu_T_n.size(); i++)
        T_n += x_n(lfsu_T_n, i) * phi_T_n[i];

      // evaluate XCH4
      RF XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        XCH4_s += x_s(lfsu_XCH4_s, i) * phi_XCH4_s[i];
      RF XCH4_n = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_n.size(); i++)
        XCH4_n += x_n(lfsu_XCH4_n, i) * phi_XCH4_n[i];

      // evaluate YH2O
      RF YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        YH2O_s += x_s(lfsu_YH2O_s, i) * phi_YH2O_s[i];
      RF YH2O_n = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_n.size(); i++)
        YH2O_n += x_n(lfsu_YH2O_n, i) * phi_YH2O_n[i];

      // evaluate XC
      RF XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        XC_s += x_s(lfsu_XC_s, i) * phi_XC_s[i];
      RF XC_n = 0.0;
      for (size_type i = 0; i < lfsu_XC_n.size(); i++)
        XC_n += x_n(lfsu_XC_n, i) * phi_XC_n[i];

      RF Sw_s = (1. - Sg_s - Sh_s);
      RF Sw_n = (1. - Sg_n - Sh_n);


      // evaluate Pw
      auto BrooksCParams_s = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto BrooksCParams_n = property.hydraulicProperty.BrooksCoreyParameters(cell_outside, iplocal_n);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto por_n = property.soil.SedimentPorosity(cell_outside, iplocal_n);
      auto Pc_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */
      RF Pg_s = Pw_s + Pc_s;

      auto Pc_n = property.hydraulicProperty.CapillaryPressure(cell_outside, iplocal_n, Sw_n, Sh_n, por_n) ; /* ndim */
      RF Pg_n = Pw_n + Pc_n;/* ndim */

      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      auto Pw_s_dim = Pw_s * Xc_P;
      auto Pw_n_dim = Pw_n * Xc_P;
      auto Pg_s_dim = Pg_s * Xc_P;
      auto Pg_n_dim = Pg_n * Xc_P;
      auto T_s_dim = T_s * Xc_T;
      auto T_n_dim = T_n * Xc_T;



      // evaluate gradient of basis functions
      auto &js_Pw_s = cache_Pw[order_p].evaluateJacobian(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &js_v_Pw_s = cache_Pw[order_p].evaluateJacobian(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &js_Sg_s = cache_Sg[order_s].evaluateJacobian(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &js_v_Sg_s = cache_Sg[order_s].evaluateJacobian(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &js_Sh_s = cache_Sh[order_s].evaluateJacobian(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &js_v_Sh_s = cache_Sh[order_s].evaluateJacobian(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order_x].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order_x].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order_x].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order_x].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &js_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      auto &js_Pw_n = cache_Pw[order_p].evaluateJacobian(iplocal_n, lfsu_Pw_n.finiteElement().localBasis());
      auto &js_v_Pw_n = cache_Pw[order_p].evaluateJacobian(iplocal_n, lfsv_Pw_n.finiteElement().localBasis());
      auto &js_Sg_n = cache_Sg[order_s].evaluateJacobian(iplocal_n, lfsu_Sg_n.finiteElement().localBasis());
      auto &js_v_Sg_n = cache_Sg[order_s].evaluateJacobian(iplocal_n, lfsv_Sg_n.finiteElement().localBasis());
      auto &js_Sh_n = cache_Sh[order_s].evaluateJacobian(iplocal_n, lfsu_Sh_n.finiteElement().localBasis());
      auto &js_v_Sh_n = cache_Sh[order_s].evaluateJacobian(iplocal_n, lfsv_Sh_n.finiteElement().localBasis());
      auto &js_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &js_v_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      auto &js_XCH4_n = cache_XCH4[order_x].evaluateJacobian(iplocal_n, lfsu_XCH4_n.finiteElement().localBasis());
      auto &js_v_XCH4_n = cache_XCH4[order_x].evaluateJacobian(iplocal_n, lfsv_XCH4_n.finiteElement().localBasis());
      auto &js_YH2O_n = cache_YH2O[order_x].evaluateJacobian(iplocal_n, lfsu_YH2O_n.finiteElement().localBasis());
      auto &js_v_YH2O_n = cache_YH2O[order_x].evaluateJacobian(iplocal_n, lfsv_YH2O_n.finiteElement().localBasis());
      auto &js_XC_n = cache_XC[order_x].evaluateJacobian(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &js_v_XC_n = cache_XC[order_x].evaluateJacobian(iplocal_n, lfsv_XC_n.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        jac.mv(js_Pw_s[i][0], gradphi_Pw_s[i]);
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        jac.mv(js_v_Pw_s[i][0], gradpsi_Pw_s[i]);

      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        jac.mv(js_Sg_s[i][0], gradphi_Sg_s[i]);
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
        jac.mv(js_v_Sg_s[i][0], gradpsi_Sg_s[i]);

      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        jac.mv(js_Sh_s[i][0], gradphi_Sh_s[i]);
      for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
        jac.mv(js_v_Sh_s[i][0], gradpsi_Sh_s[i]);

      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        jac.mv(js_XCH4_s[i][0], gradphi_XCH4_s[i]);
      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
        jac.mv(js_v_XCH4_s[i][0], gradpsi_XCH4_s[i]);

      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        jac.mv(js_YH2O_s[i][0], gradphi_YH2O_s[i]);
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
        jac.mv(js_v_YH2O_s[i][0], gradpsi_YH2O_s[i]);

      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        jac.mv(js_XC_s[i][0], gradphi_XC_s[i]);
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        jac.mv(js_v_XC_s[i][0], gradpsi_XC_s[i]);

      jac = geo_outside.jacobianInverseTransposed(iplocal_n);
      for (size_type i = 0; i < lfsu_Pw_n.size(); i++)
        jac.mv(js_Pw_n[i][0], gradphi_Pw_n[i]);
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
        jac.mv(js_v_Pw_n[i][0], gradpsi_Pw_n[i]);

      for (size_type i = 0; i < lfsu_Sg_n.size(); i++)
        jac.mv(js_Sg_n[i][0], gradphi_Sg_n[i]);
      for (size_type i = 0; i < lfsv_Sg_n.size(); i++)
        jac.mv(js_v_Sg_n[i][0], gradpsi_Sg_n[i]);

      for (size_type i = 0; i < lfsu_Sh_n.size(); i++)
        jac.mv(js_Sh_n[i][0], gradphi_Sh_n[i]);
      for (size_type i = 0; i < lfsv_Sh_n.size(); i++)
        jac.mv(js_v_Sh_n[i][0], gradpsi_Sh_n[i]);

      for (size_type i = 0; i < lfsu_T_n.size(); i++)
        jac.mv(js_T_n[i][0], gradphi_T_n[i]);
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
        jac.mv(js_v_T_n[i][0], gradpsi_T_n[i]);

      for (size_type i = 0; i < lfsu_XCH4_n.size(); i++)
        jac.mv(js_XCH4_n[i][0], gradphi_XCH4_n[i]);
      for (size_type i = 0; i < lfsv_XCH4_n.size(); i++)
        jac.mv(js_v_XCH4_n[i][0], gradpsi_XCH4_n[i]);

      for (size_type i = 0; i < lfsu_YH2O_n.size(); i++)
        jac.mv(js_YH2O_n[i][0], gradphi_YH2O_n[i]);
      for (size_type i = 0; i < lfsv_YH2O_n.size(); i++)
        jac.mv(js_v_YH2O_n[i][0], gradpsi_YH2O_n[i]);

      for (size_type i = 0; i < lfsu_XC_n.size(); i++)
        jac.mv(js_XC_n[i][0], gradphi_XC_n[i]);
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
        jac.mv(js_v_XC_n[i][0], gradpsi_XC_n[i]);

      // compute gradient of Pw
      gradu_Pw_s = 0.0;
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        gradu_Pw_s.axpy(x_s(lfsu_Pw_s, i), gradphi_Pw_s[i]);
      gradu_Pw_n = 0.0;
      for (size_type i = 0; i < lfsu_Pw_n.size(); i++)
        gradu_Pw_n.axpy(x_n(lfsu_Pw_n, i), gradphi_Pw_n[i]);

      // compute gradient of Sg
      gradu_Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        gradu_Sg_s.axpy(x_s(lfsu_Sg_s, i), gradphi_Sg_s[i]);
      gradu_Sg_n = 0.0;
      for (size_type i = 0; i < lfsu_Sg_n.size(); i++)
        gradu_Sg_n.axpy(x_n(lfsu_Sg_n, i), gradphi_Sg_n[i]);

      // compute gradient of Sh
      gradu_Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        gradu_Sh_s.axpy(x_s(lfsu_Sh_s, i), gradphi_Sh_s[i]);
      gradu_Sh_n = 0.0;
      for (size_type i = 0; i < lfsu_Sh_n.size(); i++)
        gradu_Sh_n.axpy(x_n(lfsu_Sh_n, i), gradphi_Sh_n[i]);

      // compute gradient of T
      gradu_T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x_s(lfsu_T_s, i), gradphi_T_s[i]);
      gradu_T_n = 0.0;
      for (size_type i = 0; i < lfsu_T_n.size(); i++)
        gradu_T_n.axpy(x_n(lfsu_T_n, i), gradphi_T_n[i]);

      // compute gradient of XCH4
      gradu_XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        gradu_XCH4_s.axpy(x_s(lfsu_XCH4_s, i), gradphi_XCH4_s[i]);
      gradu_XCH4_n = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_n.size(); i++)
        gradu_XCH4_n.axpy(x_n(lfsu_XCH4_n, i), gradphi_XCH4_n[i]);

      // compute gradient of YH2O
      gradu_YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        gradu_YH2O_s.axpy(x_s(lfsu_YH2O_s, i), gradphi_YH2O_s[i]);
      gradu_YH2O_n = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_n.size(); i++)
        gradu_YH2O_n.axpy(x_n(lfsu_YH2O_n, i), gradphi_YH2O_n[i]);

      // compute gradient of XC
      gradu_XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        gradu_XC_s.axpy(x_s(lfsu_XC_s, i), gradphi_XC_s[i]);
      gradu_XC_n = 0.0;
      for (size_type i = 0; i < lfsu_XC_n.size(); i++)
        gradu_XC_n.axpy(x_n(lfsu_XC_n, i), gradphi_XC_n[i]);

      auto K_s = property.soil.SedimentPermeabilityTensor(cell_inside, iplocal_s)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s); /* ndim */
      auto K_n = property.soil.SedimentPermeabilityTensor(cell_outside, iplocal_n)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_outside,iplocal_n, Sh_n, por_n); /* ndim */
      // auto K_s = K_s1;//
      auto permeability_s = property.soil.SedimentPermeability(cell_inside, iplocal_s)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s); /* ndim */
      auto permeability_n = property.soil.SedimentPermeability(cell_outside, iplocal_n)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_outside,iplocal_n, Sh_n, por_n); /* ndim */
      // permeability_s = 2. * (permeability_s*permeability_n)/(permeability_s+permeability_n);//
      // // K_s[1][1] = 2. * (K_s[1][1]*K_n[1][1])/(K_s[1][1]+K_n[1][1]);//
      // permeability_n = permeability_s;//

      K_s.mv(gravity, Kg_s);
      K_n.mv(gravity, Kg_n);

      // compute K * gradient of Pw
      K_s.mv(gradu_Pw_s, Kgradu_Pw_s);
      K_n.mv(gradu_Pw_n, Kgradu_Pw_n);

      // compute K * gradient of Sg
      K_s.mv(gradu_Sg_s, Kgradu_Sg_s);
      K_n.mv(gradu_Sg_n, Kgradu_Sg_n);

      // compute K * gradient of Sh
      K_s.mv(gradu_Sh_s, Kgradu_Sh_s);
      K_n.mv(gradu_Sh_n, Kgradu_Sh_n);

      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F_local, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K_n.mv(n_F_local, Kn_F_n);

      double eta_s = 1/BrooksCParams_s[1];
      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s, Sh_s, 0., 0.);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams_s[0], BrooksCParams_s[1]);/* ndim */
      auto dSwe_dSw_s = property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0., 0.);
      auto coeff_grad_Sw_s = dPc_dSwe_s * dSwe_dSw_s * property.hydraulicProperty.PcSF1(Sh_s, BrooksCParams_s[1], BrooksCParams_s[4]);

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams_s[1], BrooksCParams_s[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0., 0.);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s * std::pow( Swe_s , -eta_s ) * BrooksCParams_s[0] / Xc_P - Sg_s * coeff_grad_Sw_s *  dSwe_dSw_s ;

      double eta_n = 1/BrooksCParams_n[1];
      auto Swe_n = property.hydraulicProperty.EffectiveSw(Sw_n,Sh_n, 0., 0.);
      auto dPc_dSwe_n =  property.hydraulicProperty.dPc_dSwe(Swe_n, BrooksCParams_n[0], BrooksCParams_n[1]);/* ndim */
      auto dSwe_dSw_n = property.hydraulicProperty.dSwe_dSw(Sw_n, Sh_n, 0., 0.);
      auto coeff_grad_Sw_n = dPc_dSwe_n * dSwe_dSw_n * property.hydraulicProperty.PcSF1(Sh_n, BrooksCParams_n[1], BrooksCParams_n[4]);

      auto dPcSF1_dSh_n =  property.hydraulicProperty.dPcSF1_dSh( Sh_n, BrooksCParams_n[1], BrooksCParams_n[4]);
      auto dSwe_dSh_n = property.hydraulicProperty.dSwe_dSh(Sw_n, Sh_n, 0., 0.);
      auto coeff_grad_Sh_n = dPcSF1_dSh_n * std::pow( Swe_n , -eta_n ) * BrooksCParams_n[0] / Xc_P - Sg_n * coeff_grad_Sw_n *  dSwe_dSw_n ;


      auto Kgradu_Pg_s = Kgradu_Pw_s - coeff_grad_Sw_s * Kgradu_Sg_s + (coeff_grad_Sh_s ) * Kgradu_Sh_s;
      auto gradu_Pg_s = gradu_Pw_s - coeff_grad_Sw_s * gradu_Sg_s + (coeff_grad_Sh_s ) * gradu_Sh_s;

      auto Kgradu_Pg_n = Kgradu_Pw_n - coeff_grad_Sw_n * Kgradu_Sg_n + (coeff_grad_Sh_n ) * Kgradu_Sh_n;
      auto gradu_Pg_n = gradu_Pw_n - coeff_grad_Sw_n * gradu_Sg_n + (coeff_grad_Sh_n ) * gradu_Sh_n;

      

      double S_s = XC_s * (property.salt.MolarMass()/property.gas.MolarMass());
      double S_n = XC_n * (property.salt.MolarMass()/property.gas.MolarMass());

      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s);
      auto mu_w_s = property.water.DynamicViscosity(T_s_dim, Pw_s_dim, S_s) ; /* ndim */
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s);
      auto mu_g_s = property.gas.DynamicViscosity(T_s_dim, Pg_s_dim) ; /* ndim */

      auto krW_n = property.hydraulicProperty.krw(cell_outside, iplocal_n, Sw_n, Sh_n);
      auto mu_w_n = property.water.DynamicViscosity(T_n_dim, Pw_n_dim, S_n) ; /* ndim */
      auto krN_n = property.hydraulicProperty.krg(cell_outside, iplocal_n, Sw_n, Sh_n);
      auto mu_g_n = property.gas.DynamicViscosity(T_n_dim, Pg_n_dim); /* ndim */

      // if (krW_s != 0. || krW_n != 0.){
      //   krW_s = 2. * (krW_s*krW_n)/(krW_s+krW_n);
      // }

      // if (krN_s != 0. || krN_n != 0.){
      //   krN_s = 2. * (krN_s*krN_n)/(krN_s+krN_n);
      // }


      // auto mu_w_inverse = 0.5 * (mu_w_s+mu_w_n)/(mu_w_s*mu_w_n);
      // auto mu_g_inverse = 0.5 * (mu_g_s+mu_g_n)/(mu_g_s*mu_g_n);

      krW_s = krW_s / mu_w_s;; /* ndim */
      krN_s = krN_s / mu_g_s; /* ndim */

      krW_n = krW_n / mu_w_n;
      krN_n = krN_n / mu_g_n;


      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s);
      auto tau_n = property.soil.Tortuosity(por_n);
      // tau_s = 0.5 * (tau_s+tau_n);
      // tau_n = tau_s;
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s_dim, Pg_s_dim); /* ndim */
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s_dim, Pw_s_dim); /* ndim */
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s_dim, Pw_s_dim); /* ndim */
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s_dim, Pg_s_dim);
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n_dim, Pg_n_dim);
      // zCH4_s = 0.5*(zCH4_s+zCH4_n);
      // zCH4_n = zCH4_s;
      auto H_CH4_w_s = property.gas.SolubilityCoefficient(  T_s_dim/*K*/, S_s ); /*ndim */
      auto P_H_satu_s = property.water.SaturatedVaporPressure( T_s_dim /*K*/, S_s ); /*ndim */



      auto Cp_g_s = property.gas.Cp(T_s_dim, Pg_s_dim, zCH4_s) ; /* ndim */
      auto Cp_w_s = property.water.Cp(T_s_dim, Pw_s_dim, S_s) ; /* ndim */
      auto kth_g_s = property.gas.ThermalConductivity(T_s_dim, Pg_s_dim) ; /* ndim */
      auto kth_w_s = property.water.ThermalConductivity(T_s_dim, Pw_s_dim, S_s) ; /* ndim */
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s_dim, Peff_s * Xc_P); /* ndim */
      auto kth_s_s = property.soil.ThermalConductivity() ; /* ndim */
      // auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * ( Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s); /* ndim */



      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n_dim, Pg_n_dim); /* ndim */
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n_dim, Pw_n_dim); /* ndim */
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n_dim, Pw_n_dim); /* ndim */
      auto H_CH4_w_n = property.gas.SolubilityCoefficient(  T_n_dim/*K*/, S_n ); /*ndim */
      auto P_H_satu_n = property.water.SaturatedVaporPressure( T_n_dim /*K*/, S_n ); /*ndim */

      auto rho_g_s = property.gas.Density(T_s_dim, Pg_s_dim, zCH4_s) ; /* ndim */
      auto rho_w_s = property.water.Density(T_s_dim, Pw_s_dim, S_s) ; /* ndim */

      auto rho_g_n = property.gas.Density(T_n_dim, Pg_n_dim, zCH4_n) ; /* ndim */
      auto rho_w_n = property.water.Density(T_n_dim, Pw_n_dim, S_n) ; /* ndim */
      // rho_g_s = 0.5 * (rho_g_s+rho_g_n);
      // rho_g_n = rho_g_s;
      // rho_w_s = 0.5 * (rho_w_s+rho_w_n);
      // rho_w_n = rho_w_s;

      auto YCH4_s = property.mixture.YCH4(XCH4_s, T_s_dim, Pg_s_dim, XC_s, zCH4_s);
      auto XH2O_s = property.mixture.XH2O(YH2O_s, T_s_dim, Pg_s_dim, XC_s);

      auto YCH4_n = property.mixture.YCH4(XCH4_n, T_n_dim, Pg_n_dim, XC_n, zCH4_n);
      auto XH2O_n = property.mixture.XH2O(YH2O_n, T_n_dim, Pg_n_dim, XC_n);

      auto Cp_g_n = property.gas.Cp(T_n_dim, Pg_n_dim, zCH4_n) ; /* ndim */
      auto Cp_w_n = property.water.Cp(T_n_dim, Pw_n_dim, S_n) ; /* ndim */
      auto kth_g_n = property.gas.ThermalConductivity(T_n_dim, Pg_n_dim) ; /* ndim */
      auto kth_w_n = property.water.ThermalConductivity(T_n_dim, Pw_n_dim, S_n) ; /* ndim */
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n_dim, Peff_n * Xc_P) ; /* ndim */
      auto kth_s_n = property.soil.ThermalConductivity() ; /* ndim */
      // auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * (Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);

      //  adding terms regarding components

      // DH2O_g_s = 0.5 * (DH2O_g_s + DH2O_g_n); /* ndim */
      // DH2O_g_n = DH2O_g_s;
      // DCH4_w_s = 0.5*(DCH4_w_s+DCH4_w_n); /* ndim */
      // DCH4_w_n = DCH4_w_s;
      // DC_w_s = 0.5*(DC_w_s+DC_w_n); /* ndim */
      // DC_w_n = DC_w_s;



      // Cp_g_s = 0.5*(Cp_g_s+Cp_g_n) ; /* ndim */
      // Cp_g_n = Cp_g_s;
      // Cp_w_s = 0.5*(Cp_w_s+Cp_w_n) ; /* ndim */
      // Cp_w_n = Cp_w_s;

      // kth_g_s = 0.5*(kth_g_s+kth_g_n); /* ndim */
      // kth_g_n = kth_g_s;
      // kth_w_s = 0.5*(kth_w_s+kth_w_n); /* ndim */
      // kth_w_n = kth_w_s;
      // kth_h_s = 0.5*(kth_h_s+kth_h_n); /* ndim */
      // kth_h_n = kth_h_s;
      // kth_s_s = 0.5*(kth_s_s+kth_s_n); /* ndim */
      // kth_s_n = kth_s_s;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * ( Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s); /* ndim */
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * ( Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);

			for(int i = 0;i<dim;i++){
        v_g[i] = -(omega_s * krN_s * ( Xc_P/Xc_x * Kgradu_Pg_s[i] - Xc_rho * rho_g_s * Kg_s[i])
                  +omega_n * krN_n * ( Xc_P/Xc_x * Kgradu_Pg_n[i] - Xc_rho *  rho_g_n * Kg_n[i]));

        v_w[i] = -(omega_s * krW_s * ( Xc_P/Xc_x * Kgradu_Pw_s[i] - Xc_rho *  rho_w_s * Kg_s[i])
                  +omega_n * krW_n * ( Xc_P/Xc_x * Kgradu_Pw_n[i] - Xc_rho *  rho_w_n * Kg_n[i]));
      }
      auto velocity_g_s = krN_s * ( Xc_P/Xc_x * Kgradu_Pg_s - Xc_rho * rho_g_s * Kg_s);
      auto velocity_w_s = krW_s * ( Xc_P/Xc_x * Kgradu_Pw_s - Xc_rho * rho_w_s * Kg_s);

      auto velocity_g_n = krN_n * ( Xc_P/Xc_x * Kgradu_Pg_n - Xc_rho * rho_g_n * Kg_n);
      auto velocity_w_n = krW_n * ( Xc_P/Xc_x * Kgradu_Pw_n - Xc_rho * rho_w_n * Kg_n);

      double normalflux_g = (v_g ) * n_F_local;
      double normalflux_w = (v_w) * n_F_local;
      double normalflux_x = (omega_s * gradu_XC_s + omega_n * gradu_XC_n) * n_F_local;
      double normalflux_T = (omega_s * gradu_T_s + omega_n * gradu_T_n) * n_F_local;
      double normaldiff_xch4 = (omega_s * gradu_XCH4_s + omega_n * gradu_XCH4_n)* n_F_local;
      double normaldiff_yh2o = (omega_s * gradu_YH2O_s + omega_n * gradu_YH2O_n)* n_F_local;


      // upwinding wrt gas-phase velocity
      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g>=0.0)
      {
        omegaup_g_s = 1.0;
        omegaup_g_n = 0.0;
      }
      else
      {
        omegaup_g_s = 0.0;
        omegaup_g_n = 1.0;
      }
      // upwinding wrt water-phase velocity
      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w>=0.0)
      {
        omegaup_w_s = 1.0;
        omegaup_w_n = 0.0;
      }
      else
      {
        omegaup_w_s = 0.0;
        omegaup_w_n = 1.0;
      }
      // if (normalflux_g>=0.0)
      // {
      //   omegaup_g_s = 0.5;
      //   omegaup_g_n = 0.5;
      // }
      // else
      // {
      //   omegaup_g_s = 0.5;
      //   omegaup_g_n = 0.5;
      // }
      // // upwinding wrt water-phase velocity
      // RF omegaup_w_s, omegaup_w_n;
      // if (normalflux_w>=0.0)
      // {
      //   omegaup_w_s = 0.5;
      //   omegaup_w_n = 0.5;
      // }
      // else
      // {
      //   omegaup_w_s = 0.5;
      //   omegaup_w_n = 0.5;
      // }

      RF omegaup_xch4_s, omegaup_xch4_n;
      if (normaldiff_xch4>=0.0)
      {
        omegaup_xch4_s = 0.5;
        omegaup_xch4_n = 0.5;
      }
      else
      {
        omegaup_xch4_s = 0.5;
        omegaup_xch4_n = 0.5;
      }
      // // upwinding wrt water-phase diff
      RF omegaup_yh2o_s, omegaup_yh2o_n;
      if (normaldiff_yh2o>=0.0)
      {
        omegaup_yh2o_s = 0.5;
        omegaup_yh2o_n = 0.5;
      }
      else
      {
        omegaup_yh2o_s = 0.5;
        omegaup_yh2o_n = 0.5;
      }

      RF omegaup_x_s, omegaup_x_n;
      if (normalflux_x>=0.0)
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      else
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>=0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }

      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);

      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      
      
      //   fluxes and diff. flux
      auto convectiveflux_CH4_g_s = omegaup_g_s * rho_g_s * (1. - YH2O_s) * velocity_g_s;//(1. - YH2O_s)
      auto convectiveflux_CH4_w_s = omegaup_w_s * rho_w_s * (XCH4_s) * velocity_w_s;
      auto convectiveflux_H2O_g_s = omegaup_g_s * rho_g_s * YH2O_s * velocity_g_s;
      auto convectiveflux_H2O_w_s = omegaup_w_s * rho_w_s * (1. - XC_s - XCH4_s) * velocity_w_s;
      auto convectiveflux_SALT_w_s = omegaup_w_s * rho_w_s * (XC_s) * velocity_w_s;
      auto convectiveflux_Heat_w_s = omegaup_w_s * rho_w_s * Cp_w_s * (T_s - T_ref) * velocity_w_s;
      auto convectiveflux_Heat_g_s = omegaup_g_s * rho_g_s * Cp_g_s * (T_s - T_ref) * velocity_g_s;

      auto j_H2O_g_s = rho_g_s * Sg_s * DH2O_g_s * gradu_YH2O_s * n_F_local;
      auto j_CH4_w_s = rho_w_s * Sw_s * DCH4_w_s * gradu_XCH4_s * n_F_local;
      auto j_SALT_w_s = rho_w_s * Sw_s * DC_w_s * (gradu_XC_s * n_F_local);//std::max(0., std::min(2.*,2.));
      auto j_H2O_w_s = - j_CH4_w_s - j_SALT_w_s;
      auto j_CH4_g_s = - j_H2O_g_s;


      auto convectiveflux_CH4_s =  convectiveflux_CH4_g_s + convectiveflux_CH4_w_s;
      auto convectiveflux_H2O_s = convectiveflux_H2O_g_s + convectiveflux_H2O_w_s;
      auto convectiveflux_Heat_s = convectiveflux_Heat_g_s + convectiveflux_Heat_w_s;

      auto diffusiveflux_CH4_s = j_CH4_g_s + j_CH4_w_s;
      auto diffusiveflux_H2O_s = j_H2O_g_s + j_H2O_w_s;
      auto diffusiveflux_SALT_s = j_SALT_w_s;
      auto diffusiveflux_Heat_s = omegaup_T_s * kth_eff_s * gradu_T_s;
      // *******************   //
      // auto conv_coeff_CH4_g = 0.5 * (  +  ) ;
      // auto conv_coeff_CH4_w = 0.5 * (  +  )  ;
      // auto conv_coeff_H2O_g = 0.5 * (  +  ) ;
      // auto conv_coeff_H2O_w = 0.5 * ( +  ) ;
      // auto conv_coeff_SALT_w = 0.5 * ( +  ) ;
      // auto conv_coeff_Heat_w = 0.5 * ( +  );
      // auto conv_coeff_Heat_g = 0.5 * ( +  );


      
      auto convectiveflux_CH4_g_n = omegaup_g_n * rho_g_n * (1. - YH2O_n) * velocity_g_n;//(1. - YH2O_n)
      auto convectiveflux_CH4_w_n = omegaup_w_n * rho_w_n * (XCH4_n) * velocity_w_n;
      auto convectiveflux_H2O_g_n = omegaup_g_n * rho_g_n * YH2O_n * velocity_g_n;
      auto convectiveflux_H2O_w_n = omegaup_w_n * rho_w_n * (1. - XC_n - XCH4_n) * velocity_w_n;
      auto convectiveflux_SALT_w_n = omegaup_w_n * rho_w_n * (XC_n) * velocity_w_n;
      auto convectiveflux_Heat_w_n = omegaup_w_n * rho_w_n * Cp_w_n * (T_n - T_ref) * velocity_w_n;
      auto convectiveflux_Heat_g_n = omegaup_g_n * rho_g_n * Cp_g_n * (T_n - T_ref) * velocity_g_n;


      auto j_H2O_g_n = rho_g_n * Sg_n * DH2O_g_n * gradu_YH2O_n * n_F_local;
      auto j_CH4_w_n = rho_w_n * Sw_n * DCH4_w_n * gradu_XCH4_n * n_F_local;
      auto j_SALT_w_n = rho_w_n * Sw_n * DC_w_n * (gradu_XC_n * n_F_local);//std::max(0., std::min(2.*,2.));

      auto j_H2O_w_n = - j_CH4_w_n - j_SALT_w_n;
      auto j_CH4_g_n = - j_H2O_g_n;


      auto convectiveflux_CH4_n = convectiveflux_CH4_g_n + convectiveflux_CH4_w_n;
      auto convectiveflux_H2O_n = convectiveflux_H2O_g_n + convectiveflux_H2O_w_n;
      auto convectiveflux_Heat_n = convectiveflux_Heat_g_n + convectiveflux_Heat_w_n;

      auto diffusiveflux_CH4_n = j_CH4_g_n + j_CH4_w_n;
      auto diffusiveflux_H2O_n = j_H2O_g_n + j_H2O_w_n;
      auto diffusiveflux_SALT_n = j_SALT_w_n;
      auto diffusiveflux_Heat_n = omegaup_T_n * kth_eff_n * gradu_T_n;


      auto convectiveflux_CH4 = -( convectiveflux_CH4_s + convectiveflux_CH4_n) * n_F_local;
      auto convectiveflux_H2O = -( convectiveflux_H2O_s + convectiveflux_H2O_n) * n_F_local;
      auto convectiveflux_SALT = -( convectiveflux_SALT_w_s + convectiveflux_SALT_w_n) * n_F_local;
      auto convectiveflux_Heat = -( convectiveflux_Heat_s + convectiveflux_Heat_n) * n_F_local;

      auto diffusiveflux_CH4 = -(omegaup_xch4_s * diffusiveflux_CH4_s +  omegaup_xch4_n * diffusiveflux_CH4_n );// - 0.5 * penalty_factor_x * (gradu_XCH4_n - gradu_XCH4_s)* n_F_local;
      auto diffusiveflux_H2O = -(omegaup_yh2o_s * diffusiveflux_H2O_s  + omegaup_yh2o_n * diffusiveflux_H2O_n) ;//- 0.5 * penalty_factor_x * (gradu_YH2O_n - gradu_YH2O_s)* n_F_local;
      auto diffusiveflux_SALT = -(omegaup_x_s * diffusiveflux_SALT_s + omegaup_x_n * diffusiveflux_SALT_n );//- 0.5 * penalty_factor_x  * (gradu_XC_n - gradu_XC_s)* n_F_local ;
      auto diffusiveflux_Heat = -(diffusiveflux_Heat_s + diffusiveflux_Heat_n) * n_F_local;

      

      /*ACCCUMULATE RESIDUALS*/
			double tmp=0.;
      // CH4-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4 ;

      double term_nipg_g = -theta_g * (Sg_s - Sg_n);
      double term_penalty_sg =  penalty_factor_g * (Sg_s - Sg_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sg_s, i, tmp * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sg_n, i, tmp * -psi_Pw_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sg_s, i,   term_nipg_g * krN_s * omegaup_g_s * Xc_conv_m * rho_g_s *Xc_P/Xc_x
                                    * (1. - YH2O_s) * (- coeff_grad_Sw_s) * Kn_F_s * gradpsi_Sg_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sg_n, i,  term_nipg_g * krN_n * omegaup_g_n * Xc_conv_m * rho_g_n * Xc_P/Xc_x
                                    * (1. - YH2O_n) * (- coeff_grad_Sw_n) * Kn_F_n * gradpsi_Sg_n[i] * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sg_s, i, term_penalty_sg * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sg_n, i, term_penalty_sg * -psi_Pw_n[i] * factor);
      }


      // SALT-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;

      double term_nipg_c_x = -theta_x * (XC_s - XC_n);

      double term_penalty_c =  penalty_factor_x * (XC_s-XC_n);//(diffusiveflux_SALT_n - diffusiveflux_SALT_s); //  0.5 * (rho_w_s * Sw_s * DC_w_s + rho_w_n * Sw_n * DC_w_n )*
      // // diffusion term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, tmp * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, tmp * -psi_Pw_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i,  omegaup_x_s *  term_nipg_c_x * Xc_diff_m * rho_w_s * Sw_s * DC_w_s * gradpsi_XC_s[i] * n_F_local * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i,  omegaup_x_n *  term_nipg_c_x * Xc_diff_m * rho_w_n * Sw_n * DC_w_n * gradpsi_XC_n[i] * n_F_local * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i,  term_penalty_c * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, term_penalty_c  * -psi_Pw_n[i] * factor);
      }



      // H2O-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O ;

      double term_nipg_w = -theta_w * (Pw_s - Pw_n);
      double term_penalty_w =  penalty_factor_w * (Pw_s - Pw_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i, tmp * psi_Pw_s[i] * factor);
      }

      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i, tmp * -psi_Pw_n[i] * factor);
      }
      // (non-)symmetric IP term

      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i,  term_nipg_w * (omegaup_w_s * krW_s * Xc_conv_m * rho_w_s * (1. - XC_s - XCH4_s) * Xc_P/Xc_x
                                                                   * n_F_local * gradpsi_Pw_n[i] ) * factor);//+ krN_s * omegaup_g_s * rho_g_s * YH2O_s  - Xc_rho * rho_w_s * gravity
      }
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i,  term_nipg_w * (omegaup_w_n * krW_n * Xc_conv_m * rho_w_n * (1. - XC_n - XCH4_n) * Xc_P/Xc_x
                                                                   * n_F_local * gradpsi_Pw_n[i] ) * factor );// + krN_n * omegaup_g_n * rho_g_n * YH2O_n - Xc_rho * rho_w_n * gravity
      }
      //standard IP term integral
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i, term_penalty_w * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i, term_penalty_w *-psi_Pw_n[i] * factor);
      }




      // double term_penalty_sh = penalty_factor_s * (Sh_s - Sh_n);
      // for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
      // {
      //   r_s.accumulate(lfsv_Sh_s, i, term_penalty_sh * psi_Sh_s[i] * factor);
      // }
      // for (size_type i = 0; i < lfsv_Sh_n.size(); i++)
      // {
      //   r_n.accumulate(lfsv_Sh_n, i, term_penalty_sh * -psi_Sh_n[i] * factor);
      // }


      // double term_penalty_XCH4 =   penalty_factor_g *  (XCH4_s-XCH4_n);
      // for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      // {
      //   r_s.accumulate(lfsv_XCH4_s, i, term_penalty_XCH4 * psi_XCH4_s[i] * factor);
      // }
      // for (size_type i = 0; i < lfsv_XCH4_n.size(); i++)
      // {
      //   r_n.accumulate(lfsv_XCH4_n, i,   term_penalty_XCH4 * -psi_XCH4_n[i] * factor);
      // }

      // double term_penalty_YH2O =  penalty_factor_y *  (YH2O_s-YH2O_n);
      // for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      // {
      //   r_s.accumulate(lfsv_YH2O_s, i,  term_penalty_YH2O * psi_YH2O_s[i] * factor);
      // }
      // for (size_type i = 0; i < lfsv_YH2O_n.size(); i++)
      // {
      //   r_n.accumulate(lfsv_YH2O_n, i,  term_penalty_YH2O * -psi_YH2O_n[i] * factor);
      // }

      // ENERGY balance
      tmp =  Xc_conv_h * convectiveflux_Heat +  Xc_diff_h * diffusiveflux_Heat;
      double term_nipg_T = -theta_T * (T_s - T_n);
      double term_penalty_T =  penalty_factor_T * (T_s - T_n) ;
      // diffusion term
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, tmp * psi_T_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, tmp * -psi_T_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, omegaup_T_s *  term_nipg_T * Xc_diff_h * kth_eff_s * n_F_local * gradpsi_T_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, omegaup_T_n *  term_nipg_T * Xc_diff_h * kth_eff_n * n_F_local * gradpsi_T_n[i] * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_penalty_T * kth_eff_s * psi_T_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_penalty_T * kth_eff_n * -psi_T_n[i] * factor);
      }

    } //End Quadrature Rule
  }//End of alpha_skeleton

  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG &ig,
                      const LFSU &lfsu,
                      const X &x,
                      const LFSV &lfsv,
                      R &r) const
  {
    // subspaces
    //Water pressure
    const auto &lfsv_Pw_s = lfsv.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_s = lfsu.template child<Indices::PVId_Pw>();

    //Gas Saturation
    const auto &lfsv_Sg_s = lfsv.template child<Indices::PVId_Sg>();
    const auto &lfsu_Sg_s = lfsu.template child<Indices::PVId_Sg>();

    //Hydrate Saturation
    const auto &lfsv_Sh_s = lfsv.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh_s = lfsu.template child<Indices::PVId_Sh>();

    //Temperature
    const auto &lfsv_T_s = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T_s = lfsu.template child<Indices::PVId_T>();

    //Hydrate mole fraction
    const auto &lfsv_XCH4_s = lfsv.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4_s = lfsu.template child<Indices::PVId_XCH4>();

    //Water mole fraction
    const auto &lfsv_YH2O_s = lfsv.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O_s = lfsu.template child<Indices::PVId_YH2O>();

    //Salt mole fraction
    const auto &lfsv_XC_s = lfsv.template child<Indices::PVId_C>();
    const auto &lfsu_XC_s = lfsu.template child<Indices::PVId_C>();

    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;

    // dimensions
    const int dimension = GV::dimension;
    const int dim = IG::Entity::dimension;
    const int order_p = std::max(lfsu_Pw_s.finiteElement().localBasis().order(),
                                lfsv_Pw_s.finiteElement().localBasis().order());/* If different degrees are used for different functions ? */
    const int order_x = std::max(lfsu_XC_s.finiteElement().localBasis().order(),
                               lfsv_XC_s.finiteElement().localBasis().order());
    const int order_s = std::max(lfsu_Sg_s.finiteElement().localBasis().order(),
                               lfsv_Sg_s.finiteElement().localBasis().order());
    const int order_t = std::max(lfsu_T_s.finiteElement().localBasis().order(),
                               lfsv_T_s.finiteElement().localBasis().order());

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();

    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
    auto inside_cell_center_local 	= ref_el_inside.position(0,0);
    auto inside_cell_center_global 	= geo_inside.center();

    // face geometry
    auto ref_el = referenceElement(geo);
    auto face_center_local = ref_el.position(0,0);
    auto face_center_global = geo.center();

    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = geo_inside.volume() / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    harmonic_average = 1.0;

    // get polynomial degree
    auto degree = lfsv_Pw_s.finiteElement().localBasis().order();

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pw_s(lfsu_Pw_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pw_s(lfsv_Pw_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sg_s(lfsu_Sg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sg_s(lfsv_Sg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Sh_s(lfsu_Sh_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Sh_s(lfsv_Sh_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_s(lfsu_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_s(lfsv_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_s(lfsu_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_s(lfsv_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_s(lfsu_XC_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_s(lfsv_XC_s.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> gradu_T_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_n(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_n(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);
    Dune::FieldVector<RF, dim> Kg(0.0);

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    auto intorder1 = intorder;
    if( dim == 1 )
      intorder1 = 1;
    // loop over quadrature points
    for (const auto &ip : quadratureRule(geo, intorder1))
    {
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto ip_global_s = geo_inside.global(iplocal_s);
      // std::cout << iplocal_s << "  " << ip_global_s << std::endl;
      // exit(0);
      BC bc( gv,property ) ;

      // evaluate boundary condition types for {Pw,Sg} or {Fw,Fg}
			auto bctype = bc.type(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
      auto veltype = bc.velType(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg}
			auto bcvalue = bc.value(ig, ip.position(), (*time)*Xc_t, (*dt) *Xc_t) ;
      auto velvalue = bc.velValue(ig, ip.position(), (*time)*Xc_t, (*dt) *Xc_t) ;

      // evaluate basis functions at local quadrature points
      auto &psi_Pw_s = cache_Pw[order_p].evaluateFunction(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &psi_Sg_s = cache_Sg[order_s].evaluateFunction(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &psi_Sh_s = cache_Sh[order_s].evaluateFunction(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &psi_XCH4_s = cache_XCH4[order_x].evaluateFunction(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &psi_YH2O_s = cache_YH2O[order_x].evaluateFunction(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &psi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      auto &phi_Pw_s = cache_Pw[order_p].evaluateFunction(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &phi_Sg_s = cache_Sg[order_s].evaluateFunction(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &phi_Sh_s = cache_Sh[order_s].evaluateFunction(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &phi_XCH4_s = cache_XCH4[order_x].evaluateFunction(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &phi_YH2O_s = cache_YH2O[order_x].evaluateFunction(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &phi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());

      // evaluate Pw
      RF Pw_s = 0.0;
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        Pw_s += x(lfsu_Pw_s, i) * phi_Pw_s[i];
      RF Pw_n = Pw_s;
      if (bctype[Indices::PVId_Pw] == Indices::BCId_dirichlet)
      {
        Pw_n = bcvalue[Indices::PVId_Pw] ;
      }

      // evaluate T
      RF T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        T_s += x(lfsu_T_s, i) * phi_T_s[i];

      RF T_n = T_s;
      if (bctype[Indices::PVId_T] == Indices::BCId_dirichlet)
      {
        T_n = bcvalue[Indices::PVId_T] ;
      }

      // evaluate Sh
      RF Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        Sh_s += x(lfsu_Sh_s, i) * phi_Sh_s[i];
      RF Sh_n = Sh_s;

      // evaluate Sg
      RF Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        Sg_s += x(lfsu_Sg_s, i) * phi_Sg_s[i];

      RF Sg_n = Sg_s ;
      if (bctype[Indices::PVId_Sg] == Indices::BCId_dirichlet)
      {
        Sg_n = bcvalue[Indices::PVId_Sg];
      }


      RF Sw_s = (1. - Sg_s - Sh_s);
      RF Sw_n = (1. - Sg_n - Sh_n);

      // evaluate XC
      RF XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        XC_s += x(lfsu_XC_s, i) * phi_XC_s[i];

      RF XC_n = XC_s ;
      if (bctype[Indices::PVId_C] == Indices::BCId_dirichlet)
      {
        XC_n = bcvalue[Indices::PVId_C];
      }

      // evaluate XCH4
      RF XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        XCH4_s += x(lfsu_XCH4_s, i) * phi_XCH4_s[i];
      RF XCH4_n = XCH4_s;

      if (bctype[Indices::PVId_XCH4] == Indices::BCId_dirichlet)
      {
        XCH4_n = bcvalue[Indices::PVId_XCH4];
      }

      // evaluate YH2O
      RF YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        YH2O_s += x(lfsu_YH2O_s, i) * phi_YH2O_s[i];

      RF YH2O_n = YH2O_s;
      if (bctype[Indices::PVId_YH2O] == Indices::BCId_dirichlet)
      {
        YH2O_n = bcvalue[Indices::PVId_YH2O];
      }

       // evaluate Pg
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto Pc_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */

      RF Pg_s = Pw_s + Pc_s;
      auto por_n = property.soil.SedimentPorosity(cell_inside, iplocal_s) ;
      auto Pc_n = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_n, Sh_n, por_n) ; /* ndim */

      RF Pg_n = Pw_n + Pc_n;
      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      auto Pw_s_dim = Pw_s * Xc_P;
      auto Pw_n_dim = Pw_n * Xc_P;
      auto Pg_s_dim = Pg_s * Xc_P;
      auto Pg_n_dim = Pg_n * Xc_P;
      auto T_s_dim = T_s * Xc_T;
      auto T_n_dim = T_n * Xc_T;

      double S_s = XC_s * (property.salt.MolarMass()/property.gas.MolarMass());
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s_dim, Pg_s_dim);
      auto H_CH4_w_s = property.gas.SolubilityCoefficient(  T_s_dim/*K*/, S_s ); /*ndim */
      auto P_H_satu_s = property.water.SaturatedVaporPressure( T_s_dim /*K*/, S_s ); /*ndim */
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n_dim, Pg_n_dim);
      double S_n = XC_n * (property.salt.MolarMass()/property.gas.MolarMass());
      auto H_CH4_w_n = property.gas.SolubilityCoefficient(  T_n_dim/*K*/, S_n ); /*ndim */
      auto P_H_satu_n = property.water.SaturatedVaporPressure( T_n_dim /*K*/, S_n ); /*ndim */
      auto YCH4_s = property.mixture.YCH4(XCH4_s, T_s_dim, Pg_s_dim, XC_s, zCH4_s);
      auto XH2O_s = property.mixture.XH2O(YH2O_s, T_s_dim, Pg_s_dim, XC_s);

      auto YCH4_n = property.mixture.YCH4(XCH4_n, T_n_dim, Pg_n_dim, XC_n, zCH4_n);
      auto XH2O_n = property.mixture.XH2O(YH2O_n, T_n_dim, Pg_n_dim, XC_n);


      auto K = property.soil.SedimentPermeabilityTensor(cell_inside,  iplocal_s)
      * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s);
      auto permeability = property.soil.SedimentPermeability(cell_inside,  iplocal_s)
      * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s);

      
      double eta_s = 1/BrooksCParams[1];
      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s, Sh_s, 0., 0.);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams[0], BrooksCParams[1]);/* ndim */
      auto dSwe_dSw_s = property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0., 0.);
      auto coeff_grad_Sw_s = dPc_dSwe_s * dSwe_dSw_s * property.hydraulicProperty.PcSF1(Sh_s, BrooksCParams[1], BrooksCParams[4]);

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0., 0.);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s * std::pow( Swe_s , -eta_s ) * BrooksCParams[0] / Xc_P - Sg_s * coeff_grad_Sw_s *  dSwe_dSw_s ;


      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s_dim, Pw_s_dim, S_s));
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.gas.DynamicViscosity(T_s_dim, Pg_s_dim) );

      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s);
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s_dim, Pg_s_dim);
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s_dim, Pw_s_dim);
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s_dim, Pw_s_dim);
      auto rho_g_s = property.gas.Density(T_s_dim, Pg_s_dim, zCH4_s) ;
      auto rho_w_s = property.water.Density(T_s_dim, Pw_s_dim, S_s);

      auto Cp_g_s = property.gas.Cp(T_s_dim, Pg_s_dim, zCH4_s);
      auto Cp_w_s = property.water.Cp(T_s_dim, Pw_s_dim, S_s);
      auto kth_g_s = property.gas.ThermalConductivity(T_s_dim, Pg_s_dim) ;
      auto kth_w_s = property.water.ThermalConductivity(T_s_dim, Pw_s_dim, S_s);
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s_dim, Peff_s * Xc_P);
      auto kth_s_s = property.soil.ThermalConductivity() ;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * ( Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s);
      auto h_g_s =  Cp_g_s * (T_s-T_ref) ;
      auto h_w_s =  Cp_w_s * (T_s-T_ref) ;


      double eta_n = 1/BrooksCParams[1];
      auto Swe_n = property.hydraulicProperty.EffectiveSw(Sw_n,Sh_n, 0., 0.);
      auto dPc_dSwe_n =  property.hydraulicProperty.dPc_dSwe(Swe_n, BrooksCParams[0], BrooksCParams[1]);/* ndim */
      auto dSwe_dSw_n = property.hydraulicProperty.dSwe_dSw(Sw_n, Sh_n, 0., 0.);
      auto coeff_grad_Sw_n = dPc_dSwe_n * dSwe_dSw_n * property.hydraulicProperty.PcSF1(Sh_n, BrooksCParams[1], BrooksCParams[4]);

      auto dPcSF1_dSh_n =  property.hydraulicProperty.dPcSF1_dSh( Sh_n, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_n = property.hydraulicProperty.dSwe_dSh(Sw_n, Sh_n, 0., 0.);
      auto coeff_grad_Sh_n = dPcSF1_dSh_n * std::pow( Swe_n , -eta_n ) * BrooksCParams[0] / Xc_P - Sg_n * coeff_grad_Sw_n *  dSwe_dSw_n ;


      auto krW_n = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n_dim, Pw_n_dim, S_n));
      auto krN_n = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.gas.DynamicViscosity(T_n_dim, Pg_n_dim) );

      auto tau_n = property.soil.Tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n_dim, Pg_n_dim);
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n_dim, Pw_n_dim);
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n_dim, Pw_n_dim);

      auto rho_g_n = property.gas.Density(T_n_dim, Pg_n_dim, zCH4_n) ;
      auto rho_w_n = property.water.Density(T_n_dim, Pw_n_dim, S_n);

      auto Cp_g_n = property.gas.Cp(T_n_dim, Pg_n_dim, zCH4_n);
      auto Cp_w_n = property.water.Cp(T_n_dim, Pw_n_dim, S_n);
      auto kth_g_n = property.gas.ThermalConductivity(T_n_dim, Pg_n_dim) ;
      auto kth_w_n = property.water.ThermalConductivity(T_n_dim, Pw_n_dim, S_n);
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n_dim, Peff_n * Xc_P);
      auto kth_s_n = property.soil.ThermalConductivity() ;
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * ( Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
      auto h_g_n =  Cp_g_n * (T_n-T_ref) ;
      auto h_w_n =  Cp_w_n * (T_n-T_ref) ;

      omega_s = 0.5;
      omega_n = 0.5;

			auto normalgravity = gravity * n_F_local;

      // evaluate gradient of basis functions
      auto &js_Pw_s = cache_Pw[order_p].evaluateJacobian(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &js_v_Pw_s = cache_Pw[order_p].evaluateJacobian(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &js_Sg_s = cache_Sg[order_s].evaluateJacobian(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &js_v_Sg_s = cache_Sg[order_s].evaluateJacobian(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &js_Sh_s = cache_Sh[order_s].evaluateJacobian(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &js_v_Sh_s = cache_Sh[order_s].evaluateJacobian(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order_x].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order_x].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order_x].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order_x].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &js_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        jac.mv(js_Pw_s[i][0], gradphi_Pw_s[i]);
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        jac.mv(js_v_Pw_s[i][0], gradpsi_Pw_s[i]);

      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        jac.mv(js_Sg_s[i][0], gradphi_Sg_s[i]);
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
        jac.mv(js_v_Sg_s[i][0], gradpsi_Sg_s[i]);

      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        jac.mv(js_Sh_s[i][0], gradphi_Sh_s[i]);
      for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
        jac.mv(js_v_Sh_s[i][0], gradpsi_Sh_s[i]);

      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        jac.mv(js_XCH4_s[i][0], gradphi_XCH4_s[i]);
      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
        jac.mv(js_v_XCH4_s[i][0], gradpsi_XCH4_s[i]);

      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        jac.mv(js_YH2O_s[i][0], gradphi_YH2O_s[i]);
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
        jac.mv(js_v_YH2O_s[i][0], gradpsi_YH2O_s[i]);

      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        jac.mv(js_XC_s[i][0], gradphi_XC_s[i]);
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        jac.mv(js_v_XC_s[i][0], gradpsi_XC_s[i]);

      // compute gradient of Pw
      gradu_Pw_s = 0.0;
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        gradu_Pw_s.axpy(x(lfsu_Pw_s, i), gradphi_Pw_s[i]);
      gradu_Pw_n = gradu_Pw_s;

      // compute gradient of Sg
      gradu_Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        gradu_Sg_s.axpy(x(lfsu_Sg_s, i), gradphi_Sg_s[i]);
      gradu_Sg_n = gradu_Sg_s;

      // compute gradient of Sh
      gradu_Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        gradu_Sh_s.axpy(x(lfsu_Sh_s, i), gradphi_Sh_s[i]);
      gradu_Sh_n = gradu_Sh_s;

      // compute gradient of T
      gradu_T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x(lfsu_T_s, i), gradphi_T_s[i]);
      gradu_T_n = gradu_T_s;

      // compute gradient of XCH4
      gradu_XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        gradu_XCH4_s.axpy(x(lfsu_XCH4_s, i), gradphi_XCH4_s[i]);
      gradu_XCH4_n = gradu_XCH4_s;

      // compute gradient of YH2O
      gradu_YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        gradu_YH2O_s.axpy(x(lfsu_YH2O_s, i), gradphi_YH2O_s[i]);
      gradu_YH2O_n = gradu_YH2O_s;

      // compute gradient of XC
      gradu_XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        gradu_XC_s.axpy(x(lfsu_XC_s, i), gradphi_XC_s[i]);
      gradu_XC_n = gradu_XC_s;

      K.mv(gravity, Kg);

      // compute K * gradient of Pw
      K.mv(gradu_Pw_s, Kgradu_Pw_s);
      K.mv(gradu_Pw_n, Kgradu_Pw_n);

      // compute K * gradient of Sg
      K.mv(gradu_Sg_s, Kgradu_Sg_s);
      K.mv(gradu_Sg_n, Kgradu_Sg_n);

      // compute K * gradient of Sh
      K.mv(gradu_Sh_s, Kgradu_Sh_s);
      K.mv(gradu_Sh_n, Kgradu_Sh_n);

      Dune::FieldVector<RF, dim> Kn_F_s;
      K.mv(n_F_local, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K.mv(n_F_local, Kn_F_n);

	    // evaluate normal flux of Pw i.e. grad_Pw.n
      RF grad_Pw_s = gradu_Pw_s * n_F_local;
      RF grad_Pw_n = grad_Pw_s;
      // if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
      // {
      //   grad_Pw_n = bcvalue[Indices::PVId_Pw];//(-1./(K*krW_n)) * velvalue[Indices::BCId_water] + rho_w_n * normalgravity;//
      // }
      // if (veltype[Indices::BCId_water] == Indices::BCId_neumann)
      // {
      //   grad_Pw_n = (-1./(K*krW_n)) * velvalue[Indices::BCId_water] + rho_w_n * normalgravity;//
      // }

      // evaluate normal flux of Sh
      RF grad_Sh_s = gradu_Sh_s * n_F_local;
      RF grad_Sh_n = grad_Sh_s;

      // evaluate normal flux of Sg
      RF grad_Sg_s = gradu_Sg_s * n_F_local;
      RF grad_Sg_n = grad_Sg_s;

      // if (veltype[Indices::BCId_gas] == Indices::BCId_neumann)
      // {
      // //   std::cout << coeff_grad_Sw_n << " " << dPc_dSwe_n << " " << dSwe_dSw_n << std::endl;
      //   grad_Sg_n = 0.0;
      //   if (krN_n > 0.){
      //     grad_Sg_n = ((1./(K*krN_n)) * velvalue[Indices::BCId_gas] + grad_Pw_n - rho_g_n * normalgravity
      //     + (coeff_grad_Sh_n - coeff_grad_Sw_n) * grad_Sh_n) / coeff_grad_Sw_n;// NOTE: put the correct coefficients K krg and Mug instead of 1.
      //   }
      // }

      // evaluate normal flux of T
      RF grad_T_s = gradu_T_s * n_F_local;
      RF grad_T_n = grad_T_s;
      if (veltype[Indices::BCId_heat] == Indices::BCId_neumann)
      {
        grad_T_n = velvalue[Indices::BCId_heat];
      }

      // evaluate normal flux of XCH4
      RF grad_XCH4_s = gradu_XCH4_s * n_F_local;
      RF grad_XCH4_n = grad_XCH4_s;
      if (bctype[Indices::PVId_XCH4] == Indices::BCId_neumann)
      {
        grad_XCH4_n = bcvalue[Indices::PVId_XCH4];
      }
      // evaluate normal flux of YH2O
      RF grad_YH2O_s = gradu_YH2O_s * n_F_local;
      RF grad_YH2O_n = grad_YH2O_s;
      if (bctype[Indices::PVId_YH2O] == Indices::BCId_neumann)
      {
        grad_XCH4_n = bcvalue[Indices::PVId_YH2O];
      }
      // evaluate normal flux of XC
      RF grad_XC_s = gradu_XC_s * n_F_local;
      RF grad_XC_n = grad_XC_s;
      if (veltype[Indices::BCId_salt] == Indices::BCId_neumann)
      {
        grad_XC_n = velvalue[Indices::BCId_salt];
      }

      auto gradu_Pg_s = gradu_Pw_s - coeff_grad_Sw_s * gradu_Sg_s + (coeff_grad_Sh_s ) * gradu_Sh_s;
      auto gradu_Pg_n = gradu_Pw_n - coeff_grad_Sw_n * gradu_Sg_n + (coeff_grad_Sh_n) * gradu_Sh_n;
      auto Kgradu_Pg_s = Kgradu_Pw_s - coeff_grad_Sw_s * Kgradu_Sg_s + (coeff_grad_Sh_s ) * Kgradu_Sh_s;
      auto Kgradu_Pg_n = Kgradu_Pw_n - coeff_grad_Sw_n * Kgradu_Sg_n + (coeff_grad_Sh_n) * Kgradu_Sh_n;
      auto grad_Pg_s = grad_Pw_s - coeff_grad_Sw_s * grad_Sg_s + (coeff_grad_Sh_s ) * grad_Sh_s;
      auto grad_Pg_n = grad_Pw_n - coeff_grad_Sw_n * grad_Sg_n + (coeff_grad_Sh_n) * grad_Sh_n;
      // if (veltype[Indices::BCId_gas] == Indices::BCId_neumann)
      // {
      //   grad_Pg_n = 0.0;
      //   if (krN_n > 0.){
      //   grad_Pg_n = (-1./(K*krN_n)) * velvalue[Indices::BCId_gas] + rho_g_n * normalgravity;//velvalue[Indices::BCId_gas];
      //   }
      // }



			double tmp = 0.;
      auto normalvelocity_g_s = krN_s * (Xc_P/Xc_x * Kgradu_Pg_s - Xc_rho * rho_g_s * Kg) * n_F_local;

      auto normalvelocity_w_s = krW_s * (Xc_P/Xc_x * Kgradu_Pw_s - Xc_rho * rho_w_s * Kg) * n_F_local;

      auto normalvelocity_g_n = krN_n * (Xc_P/Xc_x * Kgradu_Pg_n - Xc_rho * rho_g_n * Kg) * n_F_local;
      if (veltype[Indices::BCId_gas] == Indices::BCId_neumann){
        normalvelocity_g_n = velvalue[Indices::BCId_gas];
        //std::cout << "neuman gas velocity"<< std::endl;
      }

      auto normalvelocity_w_n = krW_n * (Xc_P/Xc_x * Kgradu_Pw_n - Xc_rho * rho_w_n * Kg) * n_F_local;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann){
        normalvelocity_w_n = velvalue[Indices::BCId_water];
        //std::cout << "neuman water velocity"<< std::endl;
      }

      double normalflux_g = -(omega_s * normalvelocity_g_s + omega_n * normalvelocity_g_n); //
      double normalflux_w = -(omega_s * normalvelocity_w_s + omega_n * normalvelocity_w_n); //
      double normalflux_x = ( omega_s * grad_XC_s +omega_n * grad_XC_n); //
      double normalflux_T = ( omega_s * grad_T_s +omega_n * grad_T_n); //
      double normaldiff_xch4 = (omega_s * grad_XCH4_s +omega_n * grad_XCH4_n); //
      double normaldiff_yh2o = (omega_s * grad_YH2O_s +omega_n * grad_YH2O_n); //
      // upwinding wrt gas-phase velocity
      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g>=0.0)
      {
        omegaup_g_s = 0.5;
        omegaup_g_n = 0.5;
      }
      else
      {
        omegaup_g_s = 0.5;
        omegaup_g_n = 0.5;
      }
      // upwinding wrt water-phase velocity
      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w>=0.0)
      {
        omegaup_w_s = 0.5;
        omegaup_w_n = 0.5;
      }
      else
      {
        omegaup_w_s = 0.5;
        omegaup_w_n = 0.5;
      }
      RF omegaup_xch4_s, omegaup_xch4_n;
      if (normaldiff_xch4>=0.0)
      {
        omegaup_xch4_s = 0.5;
        omegaup_xch4_n = 0.5;
      }
      else
      {
        omegaup_xch4_s = 0.5;
        omegaup_xch4_n = 0.5;
      }
      // // upwinding wrt water-phase diff
      RF omegaup_yh2o_s, omegaup_yh2o_n;
      if (normaldiff_yh2o>=0.0)
      {
        omegaup_yh2o_s = 0.5;
        omegaup_yh2o_n = 0.5;
      }
      else
      {
        omegaup_yh2o_s = 0.5;
        omegaup_yh2o_n = 0.5;
      }

      RF omegaup_x_s, omegaup_x_n;
      if (normalflux_x>0.0)
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      else
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>=0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }


      //   fluxes and diff. flux
      auto convectiveflux_CH4_g_s = omegaup_g_s * rho_g_s * (1. - YH2O_s) * normalvelocity_g_s;//(1. - YH2O_s)
      auto convectiveflux_CH4_w_s = omegaup_w_s * rho_w_s * (XCH4_s) * normalvelocity_w_s;
      auto convectiveflux_H2O_g_s = omegaup_g_s * rho_g_s * YH2O_s * normalvelocity_g_s;
      auto convectiveflux_H2O_w_s = omegaup_w_s * rho_w_s * (1. - XC_s - XCH4_s) * normalvelocity_w_s;
      auto convectiveflux_SALT_w_s = omegaup_w_s * rho_w_s * (XC_s) * normalvelocity_w_s;
      auto convectiveflux_Heat_w_s = omegaup_w_s * rho_w_s * Cp_w_s * (T_s - T_ref) * normalvelocity_w_s;
      auto convectiveflux_Heat_g_s = omegaup_g_s * rho_g_s * Cp_g_s * (T_s - T_ref) * normalvelocity_g_s;

      auto j_H2O_g_s = omegaup_yh2o_s * rho_g_s * Sg_s * DH2O_g_s * grad_YH2O_s;
      auto j_CH4_w_s = omegaup_xch4_s * rho_w_s * Sw_s * DCH4_w_s * grad_XCH4_s;
      auto j_SALT_w_s = omegaup_x_s * rho_w_s * Sw_s * DC_w_s * grad_XC_s;//std::max(0., std::min(2*,2.));
      auto j_H2O_w_s = - j_CH4_w_s - j_SALT_w_s;
      auto j_CH4_g_s = - j_H2O_g_s;

      auto diffusiveflux_CH4_s = j_CH4_g_s + j_CH4_w_s;
      auto diffusiveflux_H2O_s = j_H2O_g_s + j_H2O_w_s;
      auto diffusiveflux_SALT_s = j_SALT_w_s;
      auto diffusiveflux_Heat_s = omegaup_T_s * kth_eff_s * grad_T_s; // k_eff will be harmonic_average of k_eff_s and k_eff_n

      //   *******************   //
      auto convectiveflux_CH4_g_n = omegaup_g_n * rho_g_n * Sg_n * (1. - YH2O_n) * normalvelocity_g_n;//(1. - YH2O_n)
      auto convectiveflux_CH4_w_n = omegaup_w_n * rho_w_n * Sw_n * (XCH4_n) * normalvelocity_w_n;
      auto convectiveflux_H2O_g_n = omegaup_g_n * rho_g_n * Sg_n * YH2O_n * normalvelocity_g_n;
      auto convectiveflux_H2O_w_n = omegaup_w_n * rho_w_n * Sw_n * (1. - XC_n - XCH4_n) * normalvelocity_w_n;
      auto convectiveflux_SALT_w_n = omegaup_w_n * rho_w_n * Sw_n * (XC_n) * normalvelocity_w_n;
      auto convectiveflux_Heat_w_n = omegaup_w_n * rho_w_n * Sw_n * Cp_w_n * (T_n - T_ref) * normalvelocity_w_n;
      auto convectiveflux_Heat_g_n = omegaup_g_n * rho_g_n * Sg_n * Cp_g_n * (T_n - T_ref) * normalvelocity_g_n;

      auto j_H2O_g_n =  omegaup_yh2o_n * rho_g_n * Sg_n * DH2O_g_n * grad_YH2O_n;
      auto j_CH4_w_n = omegaup_xch4_n * rho_w_n * Sw_n * DCH4_w_n * grad_XCH4_n;
      auto j_SALT_w_n = omegaup_x_n * rho_w_n * Sw_n * DC_w_n * grad_XC_n;//std::max(0., std::min(2*,2.));
      auto j_H2O_w_n = - j_CH4_w_n - j_SALT_w_n;
      auto j_CH4_g_n = - j_H2O_g_n;


      auto diffusiveflux_CH4_n = j_CH4_g_n + j_CH4_w_n;
      auto diffusiveflux_H2O_n = j_H2O_g_n + j_H2O_w_n;
      auto diffusiveflux_SALT_n = j_SALT_w_n;
      auto diffusiveflux_Heat_n = omegaup_T_n * kth_eff_n * grad_T_n; // k_eff will be harmonic_average of k_eff_s and k_eff_n

      auto convectiveflux_CH4_g = convectiveflux_CH4_g_s +  convectiveflux_CH4_g_n;
      if (veltype[Indices::BCId_gas] == Indices::BCId_neumann){
        convectiveflux_CH4_g =  ( rho_g_n * (1. - YH2O_n)) * normalvelocity_g_n; //(1. - YH2O_n)  rho_g_s * (1. - YH2O_s) +
      }
      auto convectiveflux_CH4_w = convectiveflux_CH4_w_s + convectiveflux_CH4_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann || bctype[Indices::PVId_Pw] == Indices::BCId_neumann){
        convectiveflux_CH4_w = (  rho_w_n * (XCH4_n)) * normalvelocity_w_n; // rho_w_s * (XCH4_s) +
      }

      auto convectiveflux_H2O_g =  convectiveflux_H2O_g_s + convectiveflux_H2O_g_n;
      if (veltype[Indices::BCId_gas] == Indices::BCId_neumann){
        convectiveflux_H2O_g =  ( rho_g_n * (YH2O_n)) * normalvelocity_g_n; // rho_g_s * (YH2O_s) +
      }
      auto convectiveflux_H2O_w = convectiveflux_H2O_w_s + convectiveflux_H2O_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann || bctype[Indices::PVId_Pw] == Indices::BCId_neumann){
        convectiveflux_H2O_w =  ( rho_w_n * (1. - XC_n - XCH4_n)) * normalvelocity_w_n;// rho_w_s * (1. - XC_s - XCH4_s) +
      }

      auto convectiveflux_SALT_w =  convectiveflux_SALT_w_s +  convectiveflux_SALT_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann || bctype[Indices::PVId_Pw] == Indices::BCId_neumann){
        convectiveflux_SALT_w =  ( rho_w_n * (XC_n)) * normalvelocity_w_n; //  rho_w_s * ( XC_s ) +
      }

      auto convectiveflux_Heat_g = convectiveflux_Heat_g_s + convectiveflux_Heat_g_n;
      if (veltype[Indices::BCId_gas] == Indices::BCId_neumann){
        convectiveflux_Heat_g =  (  rho_g_n * Cp_g_n * (T_n - T_ref)) * normalvelocity_g_n; // rho_g_s * Cp_g_s * (T_s - T_ref) +
      }

      auto convectiveflux_Heat_w = convectiveflux_Heat_w_s + convectiveflux_Heat_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann || bctype[Indices::PVId_Pw] == Indices::BCId_neumann){
        convectiveflux_Heat_w =  ( rho_w_n * Cp_w_n * (T_n - T_ref)) * normalvelocity_w_n; //rho_w_s * Cp_w_s * (T_s - T_ref) +// omegaup_w_n *
      }

      auto convectiveflux_CH4 = -(convectiveflux_CH4_g + convectiveflux_CH4_w);
      auto diffusiveflux_CH4 =  -(diffusiveflux_CH4_s +  diffusiveflux_CH4_n) ;//- 0.5 * penalty_factor_x / Xc_diff_m* (grad_XCH4_n - grad_XCH4_s);

      auto convectiveflux_H2O = -(convectiveflux_H2O_g + convectiveflux_H2O_w);
      auto diffusiveflux_H2O =  -(diffusiveflux_H2O_s +  diffusiveflux_H2O_n) ;//- 0.5 * penalty_factor_x / Xc_diff_m* (grad_YH2O_n - grad_YH2O_s);

      auto convectiveflux_SALT = -convectiveflux_SALT_w;//(omegaup_w_s * convectiveflux_SALT_w_s + omegaup_w_n * convectiveflux_SALT_w_n);
      auto diffusiveflux_SALT = -(diffusiveflux_SALT_s +  diffusiveflux_SALT_n);// - 0.5 * penalty_factor_x / Xc_diff_m* (grad_XC_n - grad_XC_s);
      if (veltype[Indices::BCId_salt] == Indices::BCId_neumann){
        diffusiveflux_SALT =    -( rho_w_n * Sw_n * DC_w_n) * grad_XC_n; //rho_w_s * Sw_s * DC_w_s +
      }

      auto convectiveflux_Heat = -(convectiveflux_Heat_g + convectiveflux_Heat_w);
      auto diffusiveflux_Heat =  -(diffusiveflux_Heat_s + diffusiveflux_Heat_n);
      if (veltype[Indices::BCId_heat] == Indices::BCId_neumann){
        diffusiveflux_Heat =   -(kth_eff_n * grad_T_n);
      }

      //  ACCCUMULATE RESIDUALS  //
			tmp=0.;

      // CH4-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4 ;
      double term_nipg_g = -method_g * (Sg_s - Sg_n);
      double term_penalty_sg = penalty_factor_g * (Sg_s - Sg_n);
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r.accumulate(lfsv_Sg_s, i, tmp * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r.accumulate(lfsv_Sg_s, i,  term_nipg_g * permeability * (omegaup_g_n * krN_n * Xc_conv_m * rho_g_n * Xc_P/Xc_x
                                        * (1. - YH2O_n) * (- coeff_grad_Sw_n)) * n_F_local * gradpsi_Sg_s[i] * factor); //+ omegaup_g_n * krN_n * rho_g_n  * (1. - YH2O_n) * (- coeff_grad_Sw_n)
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r.accumulate(lfsv_Sg_s, i, term_penalty_sg *  psi_Pw_s[i] * factor);
      }

      // SALT-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;
      double term_nipg_c_x = -method_x * (XC_s  - XC_n  );
      double term_penalty_c =   penalty_factor_x *  (XC_s-XC_n);//(diffusiveflux_SALT_n - diffusiveflux_SALT_s); //  0.5 * (rho_w_s * Sw_s * DC_w_s + rho_w_n * Sw_n * DC_w_n )*
      // diffusion term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i, tmp * psi_Pw_s[i] * factor);
      }

      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i,  term_nipg_c_x * Xc_diff_m *
            (  rho_w_n * Sw_n * DC_w_n  )* gradpsi_XC_s[i] * n_F_local * factor);//+ omegaup_x_n *Sw_n * rho_w_n * DC_w_n
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i,  term_penalty_c  * psi_Pw_s[i] * factor);
      }


      // H2O-component-wise mass-balance
      tmp =  Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O ;
      double term_nipg_w = -method_w * (Pw_s - Pw_n);
      double term_penalty_w = penalty_factor_w * (Pw_s - Pw_n) ;
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r.accumulate(lfsv_Pw_s, i, tmp * psi_Pw_s[i] * factor);
      }

      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r.accumulate(lfsv_Pw_s, i, term_nipg_w * permeability * (omegaup_w_n * krW_n * Xc_conv_m * rho_w_n * (1. - XC_n - XCH4_n) * Xc_P/Xc_x   )
                                                                   * n_F_local * (gradpsi_Pw_s[i]) * factor);//+ omegaup_w_n * krW_n * rho_w_n  * (1. - XC_n - XCH4_n) + omegaup_g_s * krN_s * rho_g_s *  YH2O_s

      }
      //standard IP term integral
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r.accumulate(lfsv_Pw_s, i, term_penalty_w *  psi_Pw_s[i] * factor);
      }

      // double term_penalty_sh = penalty_factor_s * (Sh_s - Sh_n);
      // // standard IP term integral
      // for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
      // {
      //   r.accumulate(lfsv_Sh_s, i, term_penalty_sh * psi_Sh_s[i] * factor);
      // }

      // double term_penalty_XCH4 =   Xc_diff_m * penalty_factor_x  * (XCH4_s - XCH4_n);
      // // standard IP term integral
      // for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      // {
      //   r.accumulate(lfsv_XCH4_s, i, term_penalty_XCH4 * psi_XCH4_s[i] * factor);
      // }

      // double term_penalty_YH2O =  Xc_diff_m * penalty_factor_y * (YH2O_s - YH2O_n);
      // // standard IP term integral
      // for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      // {
      //   r.accumulate(lfsv_YH2O_s, i,  term_penalty_YH2O * psi_YH2O_s[i] * factor);
      // }

      // ENERGY balance
      tmp =  Xc_conv_h * convectiveflux_Heat +  Xc_diff_h * diffusiveflux_Heat;
      double term_nipg_T = -method_T * (T_s - T_n);
      double term_penalty_T =  penalty_factor_T * (T_s - T_n) ;

      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i, tmp * psi_T_s[i] * factor);
      }

      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i,  Xc_diff_h * kth_eff_s * term_nipg_T * n_F_local * gradpsi_T_s[i] * factor); // in the run testAveragingXC-T there is no upwinding for sym terms
      }

      // standard IP term integral
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i, term_penalty_T * kth_eff_s * psi_T_s[i] * factor);
      }

      //Integrals regarding the NCP1
			// RF max1 = 0.5 * std::max(0., (Sg_s -1. + YCH4_s + YH2O_s)) + 0.5 * std::max(0., (Sg_n -1. + YCH4_n + YH2O_n));
			// for (size_type i=0; i<lfsv_YH2O_s.size(); i++){
			// 	r.accumulate(lfsv_YH2O_s,i,( ((0.5*Sg_s+0.5*Sg_n) - max1) * psi_YH2O_s[i]  *factor));
			// }
      // Integrals regarding the NCP1
			// RF max1 = std::max(0., (Sg_s -1. + YCH4_s + YH2O_s)) - std::max(0., (Sg_n -1. + YCH4_n + YH2O_n));
			// for (size_type i=0; i<lfsv_YH2O_s.size(); i++){
			// 	r.accumulate(lfsv_YH2O_s,i,( ((Sg_s-Sg_n) - max1) * psi_YH2O_s[i]  *factor));
			// }

			// // // Integrals regarding the NCP2
			// // RF max2 = 0.5 * std::max(0., (Sw_s -1. + XC_s + XCH4_s + XH2O_s )) + 0.5 * std::max(0., (Sw_n -1. + XC_n + XCH4_n + XH2O_n ));
			// // for (size_type i=0; i<lfsv_XCH4_s.size(); i++){
			// // 	r.accumulate(lfsv_XCH4_s,i,(((0.5*Sw_s+0.5*Sw_n) - max2) * psi_XCH4_s[i]  *factor));
			// // }
      // Integrals regarding the NCP2
			// RF max2 = std::max(0., (Sw_s -1. + XC_s + XCH4_s + XH2O_s )) - std::max(0., (Sw_n -1. + XC_n + XCH4_n + XH2O_n ));
			// for (size_type i=0; i<lfsv_XCH4_s.size(); i++){
			// 	r.accumulate(lfsv_XCH4_s,i,(((Sw_s-Sw_n) - max2) * psi_XCH4_s[i]  *factor));
			// }

      tmp = 0;

      if(  Sw_s -  (1. - XCH4_s - XH2O_s - XC_s )> eps_ap){//active set. std::max(0., std::min(1. ,  )
				tmp += 1. - XCH4_s - XH2O_s - XC_s;//Active => phase is present => summation condition holds
        //	std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			}else if( Sw_s - ( 1. - XCH4_s - XH2O_s - XC_s ) < eps_ap ){//std::max(0., std::min(1. , ) )
				tmp += Sw_s; // inactive set. Inactive => phase is absent => Sw=0
			}
      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      {
			r.accumulate(lfsv_XCH4_s , i, + tmp * psi_XCH4_s[i]  *factor);
      }

      tmp = 0;

      if( Sw_n - ( 1. - XCH4_n - XH2O_n - XC_n )  > eps_ap  ){//active set.std::max(0., std::min(1. , )
				tmp += 1. - XCH4_n - XH2O_n - XC_n;//Active => phase is present => summation condition holds
        //	std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			}else if( Sw_n - ( 1. - XCH4_n - XH2O_n - XC_n )  < eps_ap ){//std::max(0., std::min(1. ,  )
				tmp += Sw_n; // inactive set. Inactive => phase is absent => Sw=0
			}

      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      {
			r.accumulate(lfsv_XCH4_s , i, -tmp * psi_XCH4_s[i]  *factor);
      }

			// NCP -> gas phase
      tmp = 0.;
			if(  Sg_s -( 1. - YCH4_s - YH2O_s ) > eps_ap ){ //active set.std::max(0., std::min(1. ,  ))
				tmp +=  1. - YCH4_s - YH2O_s ;//Active => phase is present => summation condition holds
        //std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			}else if( Sg_s - ( 1. - YCH4_s - YH2O_s ) < eps_ap ){//std::max(0., std::min(1. , ) )
				tmp += Sg_s;// inactive set. Inactive => phase is absent => Sg=0
			}
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      {
			r.accumulate(lfsv_YH2O_s , i, tmp * psi_YH2O_s[i]  *factor);
      }

      tmp = 0.;
			if( Sg_n - ( 1. - YCH4_n - YH2O_n )> eps_ap ){ //active set.std::max(0., std::min(1. ,  ))
				tmp +=  1. - YCH4_n - YH2O_n ;//Active => phase is present => summation condition holds
        //std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			}else if( Sg_n - ( 1. - YCH4_n - YH2O_n ) < eps_ap ){//std::max(0., std::min(1. , ) )
				tmp += Sg_n;// inactive set. Inactive => phase is absent => Sg=0
			}
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      {
			r.accumulate(lfsv_YH2O_s , i, -tmp * psi_YH2O_s[i]  *factor);
      }



    } // end of quadrature rule
  } // end of alpha_boundary



};
