/*
 * FLOW_LocalOperator.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef FLOW_LOCALOPERATOR_HH_
#define FLOW_LOCALOPERATOR_HH_

using namespace Dune::PDELab;

struct ConvectionDiffusionDGMethod
{
  enum Type
  {
    NIPG,
    SIPG,
    IIPG
  };
};

template <class GV, typename U, class GFS, class FEM_P, class FEM_S, class FEM_T, class FEM_X, class FEM_Y>
class FLOW_LocalOperator : public Dune::PDELab::NumericalJacobianApplyVolume<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianVolume<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianApplySkeleton<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianSkeleton<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianApplyBoundary<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianBoundary<FLOW_LocalOperator<GV, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::FullSkeletonPattern, // matrix entries skeleton
                           public Dune::PDELab::FullVolumePattern,
                           public Dune::PDELab::LocalOperatorDefaultFlags,
                           public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
  IncludeClasses paramclass;
  const GV &gv;
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
  ConvectionDiffusionDGMethod::Type method_g;
  ConvectionDiffusionDGMethod::Type method_w;
  ConvectionDiffusionDGMethod::Type method_T;
  ConvectionDiffusionDGMethod::Type method_x;
  ConvectionDiffusionDGMethod::Type method_y;
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
  double Xc_C;
  double Xc_P;
  double Xc_T;
  double Xc_X;
  double Xc_Y;

public:
  // pattern assembly flags
  enum
  {
    doPatternVolume = true
  };
  enum
  {
    doPatternSkeleton = true
  };

  // residual assembly flags
  enum
  {
    doAlphaVolume = true
  };
  enum
  {
    doAlphaSkeleton = true
  }; // assemble skeleton term
  enum
  {
    doAlphaBoundary = true
  };

  typedef typename GV::IndexSet IndexSet;

  typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U::template LocalView<LFSCache> VectorView;

  typedef typename LFS::template Child<Indices::PVId_Pg>::Type LFS_Pg;
  typedef typename LFS::template Child<Indices::PVId_Pc>::Type LFS_Pc;
  typedef typename LFS::template Child<Indices::PVId_Sw>::Type LFS_Sw;
  typedef typename LFS::template Child<Indices::PVId_Sh>::Type LFS_Sh;
  typedef typename LFS::template Child<Indices::PVId_T>::Type LFS_T;
  typedef typename LFS::template Child<Indices::PVId_XCH4>::Type LFS_X;
  typedef typename LFS::template Child<Indices::PVId_YH2O>::Type LFS_Y;

  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Pg>> GFS_Pg; //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Pc>> GFS_Pc; //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Sw>> GFS_Sw; //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_Sh>> GFS_Sh; //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_T>> GFS_T;   //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_XCH4>> GFS_X; //
  typedef typename Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<Indices::PVId_YH2O>> GFS_Y;   //

  using LocalBasisType_Pg = typename FEM_P::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Pg = Dune::PDELab::LocalBasisCache<LocalBasisType_Pg>;
  using LocalBasisType_Pc = typename FEM_P::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Pc = Dune::PDELab::LocalBasisCache<LocalBasisType_Pc>;
  using LocalBasisType_Sw = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sw = Dune::PDELab::LocalBasisCache<LocalBasisType_Sw>;
  using LocalBasisType_Sh = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sh = Dune::PDELab::LocalBasisCache<LocalBasisType_Sh>;
  using LocalBasisType_T = typename FEM_T::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_T = Dune::PDELab::LocalBasisCache<LocalBasisType_T>;
  using LocalBasisType_XCH4 = typename FEM_X::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_XCH4 = Dune::PDELab::LocalBasisCache<LocalBasisType_XCH4>;
  using LocalBasisType_YH2O = typename FEM_Y::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_YH2O = Dune::PDELab::LocalBasisCache<LocalBasisType_YH2O>;

  //    using LocalBasisType = typename FiniteElementMap::Traits::FiniteElementType::Traits::LocalBasisType;

  // In theory it is possible that one and the same local operator is
  // called first with a finite element of one type and later with a
  // finite element of another type.  Since finite elements of different
  // type will usually produce different results for the same local
  // coordinate they cannot share a cache.  Here we use a vector of caches
  // to allow for different orders of the shape functions, which should be
  // enough to support p-adaptivity.  (Another likely candidate would be
  // differing geometry types, i.e. hybrid meshes.)

  std::vector<Cache_Pg> cache_Pg;
  std::vector<Cache_Pc> cache_Pc;
  std::vector<Cache_Sw> cache_Sw;
  std::vector<Cache_Sh> cache_Sh;
  std::vector<Cache_T> cache_T;
  std::vector<Cache_XCH4> cache_XCH4;
  std::vector<Cache_YH2O> cache_YH2O;

  // constructor stores parameters
  FLOW_LocalOperator(const GV &gv_,
                     U *unew_,
                     GFS gfs_,
                     double *time_,
                     double *dt_,
                     unsigned int intorder_ = 6,
                     ConvectionDiffusionDGMethod::Type method_g_ = ConvectionDiffusionDGMethod::NIPG,
                     ConvectionDiffusionDGMethod::Type method_w_ = ConvectionDiffusionDGMethod::NIPG,
                     ConvectionDiffusionDGMethod::Type method_T_ = ConvectionDiffusionDGMethod::NIPG,
                     ConvectionDiffusionDGMethod::Type method_x_ = ConvectionDiffusionDGMethod::NIPG,
                     ConvectionDiffusionDGMethod::Type method_y_ = ConvectionDiffusionDGMethod::NIPG,
                     double alpha_g_ = 1., double alpha_w_ = 1., double alpha_s_ = 1., double alpha_T_ = 1., double alpha_x_ = 1., double alpha_y_ = 1.)
      : gv(gv_),
        unew(unew_),
        gfs(gfs_),
        time(time_),
        dt(dt_),
        intorder(intorder_),
        method_g(method_g_), method_w(method_w_), method_T(method_T_), method_x(method_x_), method_y(method_y_),
        alpha_g(alpha_g_), alpha_w(alpha_w_), alpha_s(alpha_s_), alpha_T(alpha_T_), alpha_x(alpha_x_), alpha_y(alpha_y_),
        cache_Pg(20), cache_Pc(20), cache_Sw(20), cache_Sh(20), cache_T(20), cache_XCH4(20), cache_YH2O(20)
  {
    theta_g = 1.0;
    if (method_g == ConvectionDiffusionDGMethod::SIPG)
      theta_g = -1.0;
    if (method_g == ConvectionDiffusionDGMethod::IIPG)
      theta_g = 0.0;

    theta_w = 1.0;
    if (method_w == ConvectionDiffusionDGMethod::SIPG)
      theta_w = -1.0;
    if (method_w == ConvectionDiffusionDGMethod::IIPG)
      theta_w = 0.0;

    theta_T = 1.0;
    if (method_T == ConvectionDiffusionDGMethod::SIPG)
      theta_T = -1.0;
    if (method_T == ConvectionDiffusionDGMethod::IIPG)
      theta_T = 0.0;

    theta_x = 1.0;
    if (method_x == ConvectionDiffusionDGMethod::SIPG)
      theta_x = -1.0;
    if (method_w == ConvectionDiffusionDGMethod::IIPG)
      theta_x = 0.0;

    theta_y = 1.0;
    if (method_y == ConvectionDiffusionDGMethod::SIPG)
      theta_y = -1.0;
    if (method_y == ConvectionDiffusionDGMethod::IIPG)
      theta_y = 0.0;

    Xc_conv_m = paramclass.characteristicValue.X_convective_mass;
    Xc_conv_h = paramclass.characteristicValue.X_convective_heat;
    Xc_source_m = paramclass.characteristicValue.X_source_mass;
    Xc_source_h = paramclass.characteristicValue.X_source_heat;
    Xc_diff_m = paramclass.characteristicValue.X_diffusive_mass;
    Xc_diff_h = paramclass.characteristicValue.X_diffusive_heat;
    Xc_grav = paramclass.characteristicValue.X_gravity;

    Xc_K = paramclass.characteristicValue.permeability_c;
    Xc_mu = paramclass.characteristicValue.viscosity_c;
    Xc_rho = paramclass.characteristicValue.density_c;
    Xc_kth = paramclass.characteristicValue.thermalconductivity_c;
    Xc_C = paramclass.characteristicValue.specificheat_c;
    Xc_P = paramclass.characteristicValue.P_c;
    Xc_T = paramclass.characteristicValue.T_c;
    Xc_X = paramclass.characteristicValue.x_c;
    Xc_Y = paramclass.characteristicValue.x_c;
  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
  {

    // subspaces
    //Gas pressure
    const auto &lfsv_Pg = lfsv.template child<Indices::PVId_Pg>();
    const auto &lfsu_Pg = lfsu.template child<Indices::PVId_Pg>();

    //Capillary Pressure
    const auto &lfsv_Pc = lfsv.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc = lfsu.template child<Indices::PVId_Pc>();

    //Water Saturation
    const auto &lfsv_Sw = lfsv.template child<Indices::PVId_Sw>();
    const auto &lfsu_Sw = lfsu.template child<Indices::PVId_Sw>();

    //Hydrate Saturation
    const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

    //Temperature
    const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

    //Hydrate Saturation
    const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    //Temperature
    const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();
    

    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    typedef typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::JacobianType JacobianType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::SizeType;

    // dimensions
    const int dim = EG::Entity::dimension;
    const int order = std::max(lfsu_Pg.finiteElement().localBasis().order(),
                               lfsv_Pg.finiteElement().localBasis().order());

    // Get cell
    const auto &cell = eg.entity();

    // Get geometry
    auto geo = eg.geometry();

    // evaluate diffusion tensor at cell center, assume it is constant over elements
    auto ref_el = referenceElement(geo);
    auto localcenter = ref_el.position(0, 0);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pg(lfsu_Pg.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pg(lfsv_Pg.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc(lfsu_Pc.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc(lfsv_Pc.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T(lfsv_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4(lfsu_XCH4.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4(lfsv_XCH4.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O(lfsu_YH2O.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O(lfsv_YH2O.size());

    Dune::FieldVector<RF, dim> gradu_Pg(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pg(0.0);
    Dune::FieldVector<RF, dim> gradu_Pc(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc(0.0);
    Dune::FieldVector<RF, dim> gradu_T(0.0);
    Dune::FieldVector<RF, dim> Ktgradu_T(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O(0.0);

    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_Pg = cache_Pg[order].evaluateFunction(ip.position(), lfsu_Pg.finiteElement().localBasis());
      auto &psi_Pg = cache_Pg[order].evaluateFunction(ip.position(), lfsv_Pg.finiteElement().localBasis());
      auto &phi_Sw = cache_Sw[order].evaluateFunction(ip.position(), lfsu_Sw.finiteElement().localBasis());
      auto &psi_Sw = cache_Sw[order].evaluateFunction(ip.position(), lfsv_Sw.finiteElement().localBasis());
      auto &phi_Sh = cache_Sh[order].evaluateFunction(ip.position(), lfsu_Sh.finiteElement().localBasis());
      auto &psi_Sh = cache_Sh[order].evaluateFunction(ip.position(), lfsv_Sh.finiteElement().localBasis());
      auto &phi_Pc = cache_Pc[order].evaluateFunction(ip.position(), lfsu_Pc.finiteElement().localBasis());
      auto &psi_Pc = cache_Pc[order].evaluateFunction(ip.position(), lfsv_Pc.finiteElement().localBasis());
      auto &phi_T = cache_T[order].evaluateFunction(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &psi_T = cache_T[order].evaluateFunction(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &phi_XCH4 = cache_XCH4[order].evaluateFunction(ip.position(), lfsu_XCH4.finiteElement().localBasis());
      auto &psi_XCH4 = cache_XCH4[order].evaluateFunction(ip.position(), lfsv_XCH4.finiteElement().localBasis());
      auto &phi_YH2O = cache_YH2O[order].evaluateFunction(ip.position(), lfsu_YH2O.finiteElement().localBasis());
      auto &psi_YH2O = cache_YH2O[order].evaluateFunction(ip.position(), lfsv_YH2O.finiteElement().localBasis());

      auto ip_global = geo.global(ip.position());

      // evaluate Pg
      RF Pg = 0.0;
      for (size_type i = 0; i < lfsu_Pg.size(); i++)
        Pg += x(lfsu_Pg, i) * phi_Pg[i];

      // evaluate Sw
      RF Sw = 0.0;
      for (size_type i = 0; i < lfsu_Sw.size(); i++)
        Sw += x(lfsu_Sw, i) * phi_Sw[i];

      // evaluate Sh
      RF Sh = 0.0;
      for (size_type i = 0; i < lfsu_Sh.size(); i++)
        Sh += x(lfsu_Sh, i) * phi_Sh[i];

      // evaluate Pc
      RF Pc = 0.0;
      for (size_type i = 0; i < lfsu_Pc.size(); i++)
        Pc += x(lfsu_Pc, i) * phi_Pc[i];

      // evaluate T
      RF T = 0.0;
      for (size_type i = 0; i < lfsu_T.size(); i++)
        T += x(lfsu_T, i) * phi_T[i];
      
      // evaluate XCH4
      RF XCH4 = 0.0;
      for (size_type i = 0; i < lfsu_XCH4.size(); i++)
        XCH4 += x(lfsu_XCH4, i) * phi_XCH4[i];

      // evaluate YH2O
      RF YH2O = 0.0;
      for (size_type i = 0; i < lfsu_YH2O.size(); i++)
        YH2O += x(lfsu_YH2O, i) * phi_YH2O[i];

      // evaluate Pw
      RF Pw = Pg - Pc;
      RF Peff = (Pg * (1. - Sw - Sh) + Pw * Sw) / (1. - Sh);

      // evaluate gradient of basis functions
      auto &js_Pg = cache_Pg[order].evaluateJacobian(ip.position(), lfsu_Pg.finiteElement().localBasis());
      auto &js_v_Pg = cache_Pg[order].evaluateJacobian(ip.position(), lfsv_Pg.finiteElement().localBasis());
      auto &js_Pc = cache_Pc[order].evaluateJacobian(ip.position(), lfsu_Pc.finiteElement().localBasis());
      auto &js_v_Pc = cache_Pc[order].evaluateJacobian(ip.position(), lfsv_Pc.finiteElement().localBasis());
      auto &js_T = cache_T[order].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &js_XCH4 = cache_XCH4[order].evaluateJacobian(ip.position(), lfsu_XCH4.finiteElement().localBasis());
      auto &js_v_XCH4 = cache_XCH4[order].evaluateJacobian(ip.position(), lfsv_XCH4.finiteElement().localBasis());
      auto &js_YH2O = cache_YH2O[order].evaluateJacobian(ip.position(), lfsu_YH2O.finiteElement().localBasis());
      auto &js_v_YH2O = cache_YH2O[order].evaluateJacobian(ip.position(), lfsv_YH2O.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (size_type i = 0; i < lfsu_Pg.size(); i++)
        jac.mv(js_Pg[i][0], gradphi_Pg[i]);
      for (size_type i = 0; i < lfsv_Pg.size(); i++)
        jac.mv(js_v_Pg[i][0], gradpsi_Pg[i]);

      for (size_type i = 0; i < lfsu_Pc.size(); i++)
        jac.mv(js_Pc[i][0], gradphi_Pc[i]);
      for (size_type i = 0; i < lfsv_Pc.size(); i++)
        jac.mv(js_v_Pc[i][0], gradpsi_Pc[i]);

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
      for (size_type i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_YH2O[i][0], gradpsi_YH2O[i]);

      // compute gradient of Pg
      gradu_Pg = 0.0;
      for (size_type i = 0; i < lfsu_Pg.size(); i++)
        gradu_Pg.axpy(x(lfsu_Pg, i), gradphi_Pg[i]);

      // compute gradient of Pc
      gradu_Pc = 0.0;
      for (size_type i = 0; i < lfsu_Pc.size(); i++)
        gradu_Pc.axpy(x(lfsu_Pc, i), gradphi_Pc[i]);

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

      auto K = paramclass.problemSpecs.SedimentPermeabilityTensor(ip_global);
      K *= 1. / Xc_K; /*ndim K*/
      K *= Xc_conv_m;

      // compute K * gradient of Pg
      K.mv(gradu_Pg, Kgradu_Pg);

      // compute K * gradient of Pc
      K.mv(gradu_Pc, Kgradu_Pc);

      auto por = paramclass.problemSpecs.SedimentPorosity(ip_global);
      auto rho_g = paramclass.methane.density(T * Xc_T, Pg * Xc_P, 1.) / Xc_rho;
      auto rho_w = paramclass.water.density(T * Xc_T, Pw * Xc_P) / Xc_rho;
      auto rho_h = paramclass.hydrate.density() / Xc_rho;

      //  adding terms regarding components
      auto Sg = 1. - Sw - Sh;
      auto tau = paramclass.soil.tortuosity(por);
      auto DH2O_g = tau * por * paramclass.mixture.binaryDiffCoeffInGas(T * Xc_T, Pg * Xc_P);
      auto DCH4_w = tau * por * paramclass.mixture.binaryDiffCoeffInLiquid(T * Xc_T, Pw * Xc_P);
      //auto XCH4 = paramclass.mixture.mole_x_CH4(T * Xc_T, Pg * Xc_P);
      auto YCH4 = paramclass.mixture.mole_y_CH4(T * Xc_T, Pg * Xc_P);
      auto XH2O = paramclass.mixture.mole_x_H2O(T * Xc_T, Pg * Xc_P);
      //auto YH2O = paramclass.mixture.mole_y_H2O(T * Xc_T, Pg * Xc_P);
      auto H_M_w = paramclass.methane.henrysConstant(T * Xc_T);
      auto P_H_sat = paramclass.water.saturatedVaporPressure(T * Xc_T);
      auto zCH4 = paramclass.eos.evaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      //std::cout << "----" << zCH4 << std::endl;
      //  end of terms regarding components

      //  std::cout<< "densities: (g,w,h) " << rho_g <<", "<< rho_w<<", "<<rho_h<<std::endl;
      auto krW = paramclass.hydraulicProperty.krW(Sw, Sh) / (paramclass.water.dynamicViscosity(T * Xc_T, Pw * Xc_P) / Xc_mu);
      auto krN = paramclass.hydraulicProperty.krNW(Sw, Sh) / (paramclass.methane.dynamicViscosity(T * Xc_T, Pg * Xc_P) / Xc_mu);
      //          std::cout<<"lambdas: (w,n) " <<krW<<", "<<krN<<std::endl;
      auto suctionPressure = paramclass.hydraulicProperty.suctionPressure(Sw, Sh) / Xc_P;
      auto PcSF1 = paramclass.hydraulicProperty.PcSF1(Sh);
      //          std::cout<<"suction pressure: " << suctionPressure*PcSF1<<std::endl;
      auto q_g = Xc_source_m * paramclass.reactionKinetics.gasGenerationRate(paramclass.problemSpecs.T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      auto q_w = Xc_source_m * paramclass.reactionKinetics.waterGenerationRate(paramclass.problemSpecs.T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      auto q_h = Xc_source_m * paramclass.reactionKinetics.hydrateDissociationRate(paramclass.problemSpecs.T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      auto Q = Xc_source_h * paramclass.reactionKinetics.heatOfDissociation(paramclass.problemSpecs.T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      //          std::cout<< "source terms: (g,w,h) " << g_g <<", "<<g_w<<", "<<g_h<<std::endl;
      auto Cp_g = paramclass.methane.Cp(T * Xc_T, Pg * Xc_P, 1.) / Xc_C;
      auto Cp_w = paramclass.water.Cp(T * Xc_T, Pw * Xc_P) / Xc_C;
      auto kth_g = paramclass.methane.thermalConductivity(T * Xc_T, Pg * Xc_P) / Xc_kth;
      auto kth_w = paramclass.water.thermalConductivity(T * Xc_T, Pw * Xc_P) / Xc_kth;
      auto kth_h = paramclass.hydrate.thermalConductivity(T * Xc_T, Peff * Xc_P) / Xc_kth;
      auto kth_s = paramclass.soil.thermalConductivity(T * Xc_T, Peff * Xc_P) / Xc_kth;
      auto kth_eff = (1. - por) * kth_s + por * ((1. - Sw - Sh) * kth_g + Sw * kth_w + Sh * kth_h);
      kth_eff *= Xc_diff_h;

      //Methane
      auto Concoeff_m_g = rho_g * krN * YCH4;
      auto Concoeff_m_w = rho_w * krW * XCH4;

      auto Diffcoeff_m_g_p = DH2O_g * rho_g * Sg * -YCH4  / Pg;
      auto Diffcoeff_m_g_x = DH2O_g * rho_g * Sg * H_M_w / ( zCH4 * Pg);
      auto Diffcoeff_m_w = DCH4_w * rho_w * Sw ;

      
      
      //Water
      auto Concoeff_w_g = rho_g * krN * YH2O;
      auto Concoeff_w_w = rho_w * krW * XH2O;

      auto Diffcoeff_w_w_p = DCH4_w * rho_w * Sw * YH2O / P_H_sat ;
      auto Diffcoeff_w_w_x = DCH4_w * rho_w * Sw * Pg /  P_H_sat;
      auto Diffcoeff_w_g = DH2O_g * rho_g * Sg ;
      //auto Diffcoeff_w = -Diffcoeff_m;

      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      // rho_g * krN * (Kgradu_Pg * gradpsi_Pg[i])
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      for (size_type i = 0; i < lfsv_Pg.size(); i++)
      {
        r.accumulate(lfsv_Pg, i, ((Concoeff_m_g * (Kgradu_Pg * gradpsi_Pg[i]) + Concoeff_m_w * (Kgradu_Pg * gradpsi_Pg[i] 
                                                  - Kgradu_Pc * gradpsi_Pg[i])) - Diffcoeff_m_g_p * (gradu_Pg * gradpsi_Pg[i]) 
                                                  - Diffcoeff_m_g_x * (gradu_XCH4 * gradpsi_Pg[i]) 
                                                  - Diffcoeff_m_w * (gradu_XCH4 * gradpsi_Pg[i])- q_g * psi_Pg[i]) * factor);
      }
      // rho_w * krW * (Kgradu_Pg * gradpsi_Pc[i] - Kgradu_Pc * gradpsi_Pc[i])
      for (size_type i = 0; i < lfsv_Pc.size(); i++)
      {
        r.accumulate(lfsv_Pc, i, ((Concoeff_w_g * (Kgradu_Pg * gradpsi_Pc[i]) + Concoeff_w_w * (Kgradu_Pg * gradpsi_Pc[i] 
                                                  - Kgradu_Pc * gradpsi_Pc[i])) - Diffcoeff_w_w_p * (gradu_Pg * gradpsi_Pc[i])
                                                  - Diffcoeff_w_w_x * (gradu_YH2O * gradpsi_Pc[i])
                                                  - Diffcoeff_w_g * (gradu_YH2O * gradpsi_Pc[i]) 
                                                  - q_w * psi_Pc[i]) * factor);
      }
      for (size_type i = 0; i < lfsv_Sh.size(); i++)
      {
        r.accumulate(lfsv_Sh, i, (-q_h * psi_Sh[i]) * factor);
      }
      //Integrals regarding the NCP
			RF max1 = std::max(0., (Sg -1. + YCH4 + YH2O));
			for (size_type i=0; i<lfsv_XCH4.size(); i++){
				r.accumulate(lfsv_XCH4,i,( (Sg - max1) * psi_XCH4[i]  *factor));
			}

			// Integrals regarding the NCP
			RF max2 = std::max(0., (-Sg + XCH4 + XH2O ));
			for (size_type i=0; i<lfsv_YH2O.size(); i++){
				r.accumulate(lfsv_YH2O,i,((Sw - max2) * psi_YH2O[i]  *factor));
			}

      // // NCP -> water phase
			// double tmp = 0.;
			// //auto XH2O_alg = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,S);
			// if( ( Sw - ( 1. - XCH4 - XH2O - Xc_X ) ) > eps_ap ){//active set.
			// 	tmp += 1. - XCH4 - XH2O - Xc_X;//Active => phase is present => summation condition holds
      //   //	std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			// }else{
			// 	tmp += Sw; // inactive set. Inactive => phase is absent => Sw=0
			// }
      // for (size_type i = 0; i < lfsv_XCH4.size(); i++)
      // {
			// r.accumulate(lfsv_XCH4 , i, +tmp * psi_XCH4[i]  *factor);
      // }
			// // NCP -> gas phase
			// tmp = 0.;
			// //auto YCH4_alg = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,S,zCH4);
			// if( ( Sg - ( 1. - YCH4 - YH2O ) ) > eps_ap ){ //active set.			
			// 	tmp +=  1. - YCH4 - YH2O ;//Active => phase is present => summation condition holds
      //   //std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			// }else{
			// 	tmp += Sg;// inactive set. Inactive => phase is absent => Sg=0
			// }
      // for (size_type i = 0; i < lfsv_YH2O.size(); i++)
      // {
			// r.accumulate(lfsv_YH2O , i, +tmp * psi_YH2O[i]  *factor);
      // }

      for (size_type i = 0; i < lfsv_Sw.size(); i++)
      {
        r.accumulate(lfsv_Sw, i, (Pc - suctionPressure * PcSF1) * psi_Sw[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, (rho_g * krN * (Kgradu_Pg * gradpsi_T[i]) * Cp_g * T + rho_w * krW * (Kgradu_Pg * gradpsi_T[i] 
                                              - Kgradu_Pc * gradpsi_T[i]) * Cp_w * T + kth_eff * (gradu_T * gradpsi_T[i]) 
                                              - Q * psi_T[i]) * factor);
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
    //Gas pressure
    const auto &lfsv_Pg_s = lfsv_s.template child<Indices::PVId_Pg>();
    const auto &lfsu_Pg_s = lfsu_s.template child<Indices::PVId_Pg>();
    const auto &lfsv_Pg_n = lfsv_n.template child<Indices::PVId_Pg>();
    const auto &lfsu_Pg_n = lfsu_n.template child<Indices::PVId_Pg>();

    //Capillary Pressure
    const auto &lfsv_Pc_s = lfsv_s.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_s = lfsu_s.template child<Indices::PVId_Pc>();
    const auto &lfsv_Pc_n = lfsv_n.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_n = lfsu_n.template child<Indices::PVId_Pc>();

    //Water Saturation
    const auto &lfsv_Sw_s = lfsv_s.template child<Indices::PVId_Sw>();
    const auto &lfsu_Sw_s = lfsu_s.template child<Indices::PVId_Sw>();
    const auto &lfsv_Sw_n = lfsv_n.template child<Indices::PVId_Sw>();
    const auto &lfsu_Sw_n = lfsu_n.template child<Indices::PVId_Sw>();

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

    //Hydrate mole fraction
    const auto &lfsv_XCH4_s = lfsv_s.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4_s = lfsu_s.template child<Indices::PVId_XCH4>();
    const auto &lfsv_XCH4_n = lfsv_n.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4_n = lfsu_n.template child<Indices::PVId_XCH4>();
    //Water mole fraction
    const auto &lfsv_YH2O_s = lfsv_s.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O_s = lfsu_s.template child<Indices::PVId_YH2O>();
    const auto &lfsv_YH2O_n = lfsv_n.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O_n = lfsu_n.template child<Indices::PVId_YH2O>();

    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::SizeType;

    // dimensions
    const int dim = IG::Entity::dimension;
    //std::cout << dim << std::endl;
    const int order = std::max(
        std::max(lfsu_Pg_s.finiteElement().localBasis().order(),
                 lfsu_Pg_n.finiteElement().localBasis().order()),
        std::max(lfsv_Pg_s.finiteElement().localBasis().order(),
                 lfsv_Pg_n.finiteElement().localBasis().order()));

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();
    const auto &cell_outside = ig.outside();

    // Get geometries
    auto geo = ig.geometry();
    auto geo_inside = cell_inside.geometry();
    auto geo_outside = cell_outside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();
    auto geo_in_outside = ig.geometryInOutside();

    // evaluate permeability tensors
    auto ref_el_inside = referenceElement(geo_inside);
    auto ref_el_outside = referenceElement(geo_outside);
    auto local_inside = ref_el_inside.position(0, 0);
    auto local_outside = ref_el_outside.position(0, 0);

    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = std::min(geo_inside.volume(), geo_outside.volume()) / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    omega_s = omega_n = 0.5;
    harmonic_average = 1.0;

    //      if (weights==ConvectionDiffusionDGWeights::weightsOn)
    //        {
    //          RF delta_s = (An_F_s*n_F);
    //          RF delta_n = (An_F_n*n_F);
    //          omega_s = delta_n/(delta_s+delta_n+1e-20);
    //          omega_n = delta_s/(delta_s+delta_n+1e-20);
    //          harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
    //        }
    //      else
    //        {
    //          omega_s = omega_n = 0.5;
    //          harmonic_average = 1.0;
    //        }

    // get polynomial degree
    auto order_s = lfsv_Pg_s.finiteElement().localBasis().order();
    auto order_n = lfsv_Pg_n.finiteElement().localBasis().order();
    auto degree = std::max(order_s, order_n);

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pg_s(lfsu_Pg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pg_s(lfsv_Pg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_s(lfsu_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_s(lfsv_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_s(lfsu_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_s(lfsv_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_s(lfsu_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_s(lfsv_YH2O_s.size());

    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pg_n(lfsu_Pg_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pg_n(lfsv_Pg_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_n(lfsu_Pc_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_n(lfsv_Pc_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_n(lfsu_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_n(lfsv_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_n(lfsu_XCH4_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_n(lfsv_XCH4_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_n(lfsu_YH2O_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_n(lfsv_YH2O_n.size());

    Dune::FieldVector<RF, dim> gradu_Pg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Pg_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pg_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Pc_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_n(0.0);
    Dune::FieldVector<RF, dim> gradu_T_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_n(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_n(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd+quadrature_factor*order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto iplocal_n = geo_in_outside.global(ip.position());

      // evaluate basis functions

      auto &phi_Pg_s = cache_Pg[order].evaluateFunction(iplocal_s, lfsu_Pg_s.finiteElement().localBasis());
      auto &psi_Pg_s = cache_Pg[order].evaluateFunction(iplocal_s, lfsv_Pg_s.finiteElement().localBasis());
      auto &phi_Sw_s = cache_Sw[order].evaluateFunction(iplocal_s, lfsu_Sw_s.finiteElement().localBasis());
      auto &psi_Sw_s = cache_Sw[order].evaluateFunction(iplocal_s, lfsv_Sw_s.finiteElement().localBasis());
      auto &phi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &psi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &phi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &psi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &phi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &psi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &phi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &psi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());

      auto &phi_Pg_n = cache_Pg[order].evaluateFunction(iplocal_n, lfsu_Pg_n.finiteElement().localBasis());
      auto &psi_Pg_n = cache_Pg[order].evaluateFunction(iplocal_n, lfsv_Pg_n.finiteElement().localBasis());
      auto &phi_Sw_n = cache_Sw[order].evaluateFunction(iplocal_n, lfsu_Sw_n.finiteElement().localBasis());
      auto &psi_Sw_n = cache_Sw[order].evaluateFunction(iplocal_n, lfsv_Sw_n.finiteElement().localBasis());
      auto &phi_Sh_n = cache_Sh[order].evaluateFunction(iplocal_n, lfsu_Sh_n.finiteElement().localBasis());
      auto &psi_Sh_n = cache_Sh[order].evaluateFunction(iplocal_n, lfsv_Sh_n.finiteElement().localBasis());
      auto &phi_Pc_n = cache_Pc[order].evaluateFunction(iplocal_n, lfsu_Pc_n.finiteElement().localBasis());
      auto &psi_Pc_n = cache_Pc[order].evaluateFunction(iplocal_n, lfsv_Pc_n.finiteElement().localBasis());
      auto &phi_T_n = cache_T[order].evaluateFunction(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &psi_T_n = cache_T[order].evaluateFunction(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      auto &phi_XCH4_n = cache_XCH4[order].evaluateFunction(iplocal_n, lfsu_XCH4_n.finiteElement().localBasis());
      auto &psi_XCH4_n = cache_XCH4[order].evaluateFunction(iplocal_n, lfsv_XCH4_n.finiteElement().localBasis());
      auto &phi_YH2O_n = cache_YH2O[order].evaluateFunction(iplocal_n, lfsu_YH2O_n.finiteElement().localBasis());
      auto &psi_YH2O_n = cache_YH2O[order].evaluateFunction(iplocal_n, lfsv_YH2O_n.finiteElement().localBasis());

      auto ip_global_s = geo_inside.global(iplocal_s);
      auto ip_global_n = geo_outside.global(iplocal_n);
      //          std::cout<<"  ip_global_s = "<<ip_global_s<<"  iplocal_s = "<<iplocal_s
      //        		  <<"  ip_global_n = "<<ip_global_n<<"  iplocal_n = "<<iplocal_n<<std::endl;

      // evaluate Pg
      RF Pg_s = 0.0;
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        Pg_s += x_s(lfsu_Pg_s, i) * phi_Pg_s[i];
      RF Pg_n = 0.0;
      for (size_type i = 0; i < lfsu_Pg_n.size(); i++)
        Pg_n += x_n(lfsu_Pg_n, i) * phi_Pg_n[i];

      // evaluate Sw
      RF Sw_s = 0.0;
      for (size_type i = 0; i < lfsu_Sw_s.size(); i++)
        Sw_s += x_s(lfsu_Sw_s, i) * phi_Sw_s[i];
      RF Sw_n = 0.0;
      for (size_type i = 0; i < lfsu_Sw_n.size(); i++)
        Sw_n += x_n(lfsu_Sw_n, i) * phi_Sw_n[i];

      // evaluate Sh
      RF Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        Sh_s += x_s(lfsu_Sh_s, i) * phi_Sh_s[i];
      RF Sh_n = 0.0;
      for (size_type i = 0; i < lfsu_Sh_n.size(); i++)
        Sh_n += x_n(lfsu_Sh_n, i) * phi_Sh_n[i];

      // evaluate Pc
      RF Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        Pc_s += x_s(lfsu_Pc_s, i) * phi_Pc_s[i];
      RF Pc_n = 0.0;
      for (size_type i = 0; i < lfsu_Pc_n.size(); i++)
        Pc_n += x_n(lfsu_Pc_n, i) * phi_Pc_n[i];

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

      // evaluate Pw
      RF Pw_s = Pg_s - Pc_s;
      RF Pw_n = Pg_n - Pc_n;
      RF Peff_s = (Pg_s * (1. - Sw_s - Sh_s) + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * (1. - Sw_n - Sh_n) + Pw_n * Sw_n) / (1. - Sh_n);

      // evaluate gradient of basis functions
      auto &js_Pg_s = cache_Pg[order].evaluateJacobian(iplocal_s, lfsu_Pg_s.finiteElement().localBasis());
      auto &js_v_Pg_s = cache_Pg[order].evaluateJacobian(iplocal_s, lfsv_Pg_s.finiteElement().localBasis());
      auto &js_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &js_v_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());

      auto &js_Pg_n = cache_Pg[order].evaluateJacobian(iplocal_n, lfsu_Pg_n.finiteElement().localBasis());
      auto &js_v_Pg_n = cache_Pg[order].evaluateJacobian(iplocal_n, lfsv_Pg_n.finiteElement().localBasis());
      auto &js_Pc_n = cache_Pc[order].evaluateJacobian(iplocal_n, lfsu_Pc_n.finiteElement().localBasis());
      auto &js_v_Pc_n = cache_Pc[order].evaluateJacobian(iplocal_n, lfsv_Pc_n.finiteElement().localBasis());
      auto &js_T_n = cache_T[order].evaluateJacobian(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &js_v_T_n = cache_T[order].evaluateJacobian(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      auto &js_XCH4_n = cache_XCH4[order].evaluateJacobian(iplocal_n, lfsu_XCH4_n.finiteElement().localBasis());
      auto &js_v_XCH4_n = cache_XCH4[order].evaluateJacobian(iplocal_n, lfsv_XCH4_n.finiteElement().localBasis());
      auto &js_YH2O_n = cache_YH2O[order].evaluateJacobian(iplocal_n, lfsu_YH2O_n.finiteElement().localBasis());
      auto &js_v_YH2O_n = cache_YH2O[order].evaluateJacobian(iplocal_n, lfsv_YH2O_n.finiteElement().localBasis());

      // transform gradients of shape functions to real element

      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        jac.mv(js_Pg_s[i][0], gradphi_Pg_s[i]);
      for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        jac.mv(js_v_Pg_s[i][0], gradpsi_Pg_s[i]);
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        jac.mv(js_Pc_s[i][0], gradphi_Pc_s[i]);
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        jac.mv(js_v_Pc_s[i][0], gradpsi_Pc_s[i]);
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

      jac = geo_outside.jacobianInverseTransposed(iplocal_n);
      for (size_type i = 0; i < lfsu_Pg_n.size(); i++)
        jac.mv(js_Pg_n[i][0], gradphi_Pg_n[i]);
      for (size_type i = 0; i < lfsv_Pg_n.size(); i++)
        jac.mv(js_v_Pg_n[i][0], gradpsi_Pg_n[i]);
      for (size_type i = 0; i < lfsu_Pc_n.size(); i++)
        jac.mv(js_Pc_n[i][0], gradphi_Pc_n[i]);
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
        jac.mv(js_v_Pc_n[i][0], gradpsi_Pc_n[i]);
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

      // compute gradient of Pg
      gradu_Pg_s = 0.0;
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        gradu_Pg_s.axpy(x_s(lfsu_Pg_s, i), gradphi_Pg_s[i]);
      gradu_Pg_n = 0.0;
      for (size_type i = 0; i < lfsu_Pg_n.size(); i++)
        gradu_Pg_n.axpy(x_s(lfsu_Pg_n, i), gradphi_Pg_n[i]);

      // compute gradient of Pc
      gradu_Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        gradu_Pc_s.axpy(x_s(lfsu_Pc_s, i), gradphi_Pc_s[i]);
      gradu_Pc_n = 0.0;
      for (size_type i = 0; i < lfsu_Pc_n.size(); i++)
        gradu_Pc_n.axpy(x_n(lfsu_Pc_n, i), gradphi_Pc_n[i]);

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

      auto K_s = paramclass.problemSpecs.SedimentPermeabilityTensor(ip_global_s);
      K_s *= 1. / Xc_K;
      K_s *= Xc_conv_m;
      auto K_n = paramclass.problemSpecs.SedimentPermeabilityTensor(ip_global_n);
      K_n *= 1. / Xc_K;
      K_n *= Xc_conv_m;

      // compute K * gradient of Pg
      K_s.mv(gradu_Pg_s, Kgradu_Pg_s);
      K_n.mv(gradu_Pg_n, Kgradu_Pg_n);

      // compute K * gradient of Pc
      K_s.mv(gradu_Pc_s, Kgradu_Pc_s);
      K_n.mv(gradu_Pc_n, Kgradu_Pc_n);

      auto n_F = ig.centerUnitOuterNormal();
      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K_n.mv(n_F, Kn_F_n);

      for (int i = 0; i < paramclass.problemSpecs.dimension; i++)
      {
        v_g[i] = (omega_s * Kgradu_Pg_s[i] + omega_n * Kgradu_Pg_n[i]);
        v_w[i] = (omega_s * (Kgradu_Pg_s[i] - Kgradu_Pc_s[i]) + omega_n * (Kgradu_Pg_n[i] - Kgradu_Pc_n[i]));
      }

      double normalflux_g = -1. * (v_g * n_F_local);
      double normalflux_w = -1. * (v_w * n_F_local);

      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g >= 0.0)
      {
        omegaup_g_s = 1.0;
        omegaup_g_n = 0.0;
      }
      else
      {
        omegaup_g_s = 0.0;
        omegaup_g_n = 1.0;
      }

      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w >= 0.0)
      {
        omegaup_w_s = 1.0;
        omegaup_w_n = 0.0;
      }
      else
      {
        omegaup_w_s = 0.0;
        omegaup_w_n = 1.0;
      }

      auto por_s = paramclass.problemSpecs.SedimentPorosity(ip_global_s);
      auto rho_g_s = paramclass.methane.density(T_s * Xc_T, Pg_s * Xc_P, 1.) / Xc_rho;
      auto rho_w_s = paramclass.water.density(T_s * Xc_T, Pw_s * Xc_P) / Xc_rho;
      auto krW_s = paramclass.hydraulicProperty.krW(Sw_s, Sh_s) / (paramclass.water.dynamicViscosity(T_s * Xc_T, Pw_s * Xc_P) / Xc_mu);
      auto krN_s = paramclass.hydraulicProperty.krNW(Sw_s, Sh_s) / (paramclass.methane.dynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) / Xc_mu);
      //          std::cout<< "densities: (g,w)_s " << rho_g_s <<", "<< rho_w_s<<std::endl;
      //          std::cout<<"lambdas: (w,n)_s " <<krW_s<<", "<<krN_s<<std::endl;

      //  adding terms regarding components
      auto Sg_s = 1. - Sw_s - Sh_s;
      auto tau_s = paramclass.soil.tortuosity(por_s);
      auto DH2O_g_s = tau_s * por_s * paramclass.mixture.binaryDiffCoeffInGas(T_s * Xc_T, Pg_s * Xc_P);
      auto DCH4_w_s = tau_s * por_s * paramclass.mixture.binaryDiffCoeffInLiquid(T_s * Xc_T, Pw_s * Xc_P);
      //auto XCH4_s = paramclass.mixture.mole_x_CH4(T_s * Xc_T, Pw_s * Xc_P);
      auto YCH4_s = paramclass.mixture.mole_y_CH4(T_s * Xc_T, Pw_s * Xc_P);
      auto XH2O_s = paramclass.mixture.mole_x_H2O(T_s * Xc_T, Pw_s * Xc_P);
      //auto YH2O_s = paramclass.mixture.mole_y_H2O(T_s * Xc_T, Pw_s * Xc_P);
      auto H_M_w_s = paramclass.methane.henrysConstant(T_s * Xc_T);
      auto P_H_sat_s = paramclass.water.saturatedVaporPressure(T_s * Xc_T);
      auto zCH4_s = paramclass.eos.evaluateCompressibilityFactor(T_s * Xc_T, Pw_s * Xc_P);
      //std::cout << "----" << zCH4 << std::endl;
      //  end of terms regarding components

      auto Cp_g_s = paramclass.methane.Cp(T_s * Xc_T, Pg_s * Xc_P, 1.) / Xc_C;
      auto Cp_w_s = paramclass.water.Cp(T_s * Xc_T, Pw_s * Xc_P) / Xc_C;
      auto kth_g_s = paramclass.methane.thermalConductivity(T_s * Xc_T, Pg_s * Xc_P) / Xc_kth;
      auto kth_w_s = paramclass.water.thermalConductivity(T_s * Xc_T, Pw_s * Xc_P) / Xc_kth;
      auto kth_h_s = paramclass.hydrate.thermalConductivity(T_s * Xc_T, Peff_s * Xc_P) / Xc_kth;
      auto kth_s_s = paramclass.soil.thermalConductivity(T_s * Xc_T, Peff_s * Xc_P) / Xc_kth;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * ((1. - Sw_s - Sh_s) * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s);
      kth_eff_s *= Xc_diff_h;

      auto por_n = paramclass.problemSpecs.SedimentPorosity(ip_global_n);
      auto rho_g_n = paramclass.methane.density(T_n * Xc_T, Pg_n * Xc_P, 1.) / Xc_rho;
      auto rho_w_n = paramclass.water.density(T_n * Xc_T, Pw_n * Xc_P) / Xc_rho;
      auto krW_n = paramclass.hydraulicProperty.krW(Sw_n, Sh_n) / (paramclass.water.dynamicViscosity(T_n * Xc_T, Pw_n * Xc_P) / Xc_mu);
      auto krN_n = paramclass.hydraulicProperty.krNW(Sw_n, Sh_n) / (paramclass.methane.dynamicViscosity(T_n * Xc_T, Pg_n * Xc_P) / Xc_mu);
      //          std::cout<< "densities: (g,w)_n " << rho_g_n <<", "<< rho_w_n<<std::endl;
      //          std::cout<<"lambdas: (w,n)_n " <<krW_n<<", "<<krN_n<<std::endl;
      //  adding terms regarding components
      auto Sg_n = 1. - Sw_n - Sh_n;
      auto tau_n = paramclass.soil.tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * paramclass.mixture.binaryDiffCoeffInGas(T_n * Xc_T, Pg_n * Xc_P);
      auto DCH4_w_n = tau_n * por_n * paramclass.mixture.binaryDiffCoeffInLiquid(T_n * Xc_T, Pg_n * Xc_P);
      //auto XCH4_n = paramclass.mixture.mole_x_CH4(T_n * Xc_T, Pg_n * Xc_P);
      auto YCH4_n = paramclass.mixture.mole_y_CH4(T_n * Xc_T, Pg_n * Xc_P);
      auto XH2O_n = paramclass.mixture.mole_x_H2O(T_n * Xc_T, Pg_n * Xc_P);
      //auto YH2O_n = paramclass.mixture.mole_y_H2O(T_n * Xc_T, Pg_n * Xc_P);
      auto H_M_w_n = paramclass.methane.henrysConstant(T_n * Xc_T);
      auto P_H_sat_n = paramclass.water.saturatedVaporPressure(T_n * Xc_T);
      auto zCH4_n = paramclass.eos.evaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);
      //std::cout << "----" << zCH4 << std::endl;
      //  end of terms regarding components

      auto Cp_g_n = paramclass.methane.Cp(T_n * Xc_T, Pg_n * Xc_P, 1.) / Xc_C;
      auto Cp_w_n = paramclass.water.Cp(T_n * Xc_T, Pw_n * Xc_P) / Xc_C;
      auto kth_g_n = paramclass.methane.thermalConductivity(T_n * Xc_T, Pg_n * Xc_P) / Xc_kth;
      auto kth_w_n = paramclass.water.thermalConductivity(T_n * Xc_T, Pw_n * Xc_P) / Xc_kth;
      auto kth_h_n = paramclass.hydrate.thermalConductivity(T_n * Xc_T, Peff_n * Xc_P) / Xc_kth;
      auto kth_s_n = paramclass.soil.thermalConductivity(T_n * Xc_T, Peff_n * Xc_P) / Xc_kth;
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * ((1. - Sw_n - Sh_n) * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);
      kth_eff_n *= Xc_diff_h;

      auto krN = omegaup_g_s * krN_s + omegaup_g_n * krN_n;
      auto krW = omegaup_w_s * krW_s + omegaup_w_n * krW_n;
      auto h_g = omegaup_g_s * Cp_g_s * T_s + omegaup_g_n * Cp_g_n * T_n;
      auto h_w = omegaup_w_s * Cp_w_s * T_s + omegaup_w_n * Cp_w_n * T_n;
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);

      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      //Methane
      auto DH2O_g = omega_s * DH2O_g_s + omega_n * DH2O_g_n; // = DCH4_g
      auto Diffcoeff_m_g_p_s = rho_g_s * Sg_s * -YCH4_s / Pg_s;
      auto Diffcoeff_m_g_x_s = rho_g_s * Sg_s * H_M_w_s / ( zCH4_s * Pg_s);
      auto Diffcoeff_m_g_p_n = rho_g_n * Sg_n * -YCH4_n  / Pg_n;
      auto Diffcoeff_m_g_x_n = rho_g_n * Sg_n * H_M_w_n / ( zCH4_n * Pg_n);
      
      auto DCH4_w = omega_s * DCH4_w_s + omega_n * DCH4_w_n; // = DH2O_w
      auto Diffcoeff_m_w_s = rho_w_s * Sw_s ;
      auto Diffcoeff_m_w_n = rho_w_n * Sw_n ;
      //auto Diffcoeff_m_g_s = rho_g_s * Sg_s * (-YCH4_s + H_M_w_s / (H_M_w_s - zCH4_s * P_H_sat_s)) / Pg_s;
      //auto Diffcoeff_m_g_n = rho_g_n * Sg_n * (-YCH4_n + H_M_w_n / (H_M_w_n - zCH4_n * P_H_sat_n)) / Pg_n;

      
      //auto Diffcoeff_m_w_s = rho_w_s * Sw_s * (zCH4_s / (H_M_w_s - zCH4_s * P_H_sat_s));
      //auto Diffcoeff_m_w_n = rho_w_n * Sw_n * (zCH4_n / (H_M_w_n - zCH4_n * P_H_sat_n));

      double term_conv_m_g = -krN * (omega_s * rho_g_s * YCH4_s * (Kgradu_Pg_s * n_F_local) 
                                    + omega_n * rho_g_n * YCH4_n * (Kgradu_Pg_n * n_F_local));
      double term_conv_m_w = -krW * (omega_s * rho_w_s * XCH4_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local)) 
                                    + omega_n * rho_w_n * XCH4_n * ((Kgradu_Pg_n * n_F_local) - (Kgradu_Pc_n * n_F_local)));

      double term_diffusion_m = term_conv_m_g + term_conv_m_w + DH2O_g * (omega_s * (Diffcoeff_m_g_p_s * (gradu_Pg_s * n_F_local) 
                                              + Diffcoeff_m_g_x_s * (gradu_XCH4_s * n_F_local))
                                              + omega_n * (Diffcoeff_m_g_p_n * (gradu_Pg_n * n_F_local) 
                                              + Diffcoeff_m_g_x_n * (gradu_XCH4_n * n_F_local)))
                                              + DCH4_w * (omega_s * Diffcoeff_m_w_s * (gradu_XCH4_s * n_F_local)
                                              + omega_n * Diffcoeff_m_w_n * (gradu_XCH4_n * n_F_local)) ;

      double term_nipg_g = theta_g * (Pg_s - Pg_n);
      double term_nipg_w = theta_w * ((Pg_s - Pc_s) - (Pg_n - Pc_n));
      double term_nipg_m_x = theta_x * (XCH4_s - XCH4_n);

      double term_penalty_g = penalty_factor_g * (Pg_s - Pg_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pg_s, i, term_diffusion_m * psi_Pg_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pg_n, i, term_diffusion_m * -psi_Pg_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pg_s, i, (term_nipg_g * krN * omega_s * rho_g_s * YCH4_s * (Kn_F_s * gradpsi_Pg_s[i]) 
                                      + term_nipg_w * krW * omega_s * rho_w_s * XCH4_s  * (Kn_F_s * gradpsi_Pg_s[i])
                                      + omega_s * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_s 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_s
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_s) * (n_F_local* gradpsi_Pg_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_Pg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pg_n, i, (term_nipg_g * krN * omega_n * rho_g_n * YCH4_n * (Kn_F_n * gradpsi_Pg_n[i]) 
                                      + term_nipg_w * krW * omega_n * rho_w_n * XCH4_n  * (Kn_F_n * gradpsi_Pg_n[i])
                                      + omega_n * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_n 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_n
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_n) * (n_F_local* gradpsi_Pg_n[i])) * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pg_s, i, term_penalty_g * psi_Pg_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pg_n, i, term_penalty_g * -psi_Pg_n[i] * factor);
      }
      //Water
      
      //auto DH2O_g = omega_s * DH2O_g_s + omega_n * DH2O_g_n; // = DCH4_g
      auto Diffcoeff_w_w_p_s = rho_w_s * Sw_s * YH2O_s / P_H_sat_s;
      auto Diffcoeff_w_w_y_s = rho_w_s * Sw_s * Pg_s / P_H_sat_s;
      auto Diffcoeff_w_w_p_n = rho_w_n * Sw_n * YH2O_n  / P_H_sat_n;
      auto Diffcoeff_w_w_y_n = rho_w_n * Sw_n * Pg_n / P_H_sat_n;
      
      //auto DCH4_w = omega_s * DCH4_w_s + omega_n * DCH4_w_n; // = DH2O_w
      auto Diffcoeff_w_g_s = rho_g_s * Sg_s ;
      auto Diffcoeff_w_g_n = rho_g_n * Sg_n ;

      double term_conv_w_g = -krN * ( omega_s * rho_g_s * YH2O_s * (Kgradu_Pg_s * n_F_local) 
                                    + omega_n * rho_g_n * YH2O_n * (Kgradu_Pg_n * n_F_local));

      double term_conv_w_w = -krW * ( omega_s * rho_w_s * XH2O_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local)) 
                                    + omega_n * rho_w_n * XH2O_n * ((Kgradu_Pg_n * n_F_local) - (Kgradu_Pc_n * n_F_local)));

      double term_diffusion_w = term_conv_w_g + term_conv_w_w + DH2O_g * (omega_s * Diffcoeff_w_g_s * (gradu_YH2O_s * n_F_local) 
                              + omega_n * Diffcoeff_w_g_n * (gradu_YH2O_n * n_F_local)) 
                              + DCH4_w * (omega_s * (Diffcoeff_w_w_p_s * (gradu_Pg_s * n_F_local)
                              + Diffcoeff_w_w_y_s * (gradu_YH2O_s * n_F_local))
                              + omega_n * (Diffcoeff_w_w_p_n * (gradu_Pg_n * n_F_local)
                              + Diffcoeff_w_w_y_n * (gradu_YH2O_n * n_F_local)));


      
      double term_nipg_w_y = theta_y * (YH2O_s - YH2O_n);
      double term_penalty_w = penalty_factor_w * ((Pg_s - Pc_s) - (Pg_n - Pc_n));
      // diffusion term
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pc_s, i, term_diffusion_w * psi_Pc_s[i] * factor);
      }

      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pc_n, i, term_diffusion_w * -psi_Pc_n[i] * factor);
      }
      // (non-)symmetric IP term
                                      
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pc_s, i, (term_nipg_g * krN * omega_s * rho_g_s * YH2O_s * (Kn_F_s * gradpsi_Pc_s[i]) 
                                      + term_nipg_w * krW * omega_s * rho_w_s * XH2O_s  * (Kn_F_s * gradpsi_Pc_s[i])
                                      + omega_s * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_s 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_s
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_s) * (n_F_local* gradpsi_Pc_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pc_n, i, (term_nipg_g * krN * omega_n * rho_g_n * YH2O_n * (Kn_F_n * gradpsi_Pc_n[i]) 
                                      + term_nipg_w * krW * omega_n * rho_w_n * XH2O_n  * (Kn_F_n * gradpsi_Pc_n[i])
                                      + omega_n * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_n 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_n
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_n) * (n_F_local* gradpsi_Pc_n[i])) * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pc_s, i, term_penalty_w * psi_Pc_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pc_n, i, term_penalty_w * -psi_Pc_n[i] * factor);
      }

      double term_penalty_sw = penalty_factor_s * (Sw_s - Sw_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sw_s, i, term_penalty_sw * psi_Sw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sw_n, i, term_penalty_sw * -psi_Sw_n[i] * factor);
      }

      double term_penalty_sh = penalty_factor_s * (Sh_s - Sh_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sh_s, i, term_penalty_sh * psi_Sh_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sh_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sh_n, i, term_penalty_sh * -psi_Sh_n[i] * factor);
      }

      double term_penalty_XCH4 = penalty_factor_x * (XCH4_s - XCH4_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      {
        r_s.accumulate(lfsv_XCH4_s, i, term_penalty_XCH4 * psi_XCH4_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_XCH4_n.size(); i++)
      {
        r_n.accumulate(lfsv_XCH4_n, i, term_penalty_XCH4 * -psi_XCH4_n[i] * factor);
      }

      double term_penalty_YH2O = penalty_factor_y * (YH2O_s - YH2O_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      {
        r_s.accumulate(lfsv_YH2O_s, i, term_penalty_YH2O * psi_YH2O_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_YH2O_n.size(); i++)
      {
        r_n.accumulate(lfsv_YH2O_n, i, term_penalty_YH2O * -psi_YH2O_n[i] * factor);
      }

      double term_diffusion_T_1 = -h_g * krN * (omega_s * rho_g_s * (Kgradu_Pg_s * n_F_local) + omega_n * rho_g_n * (Kgradu_Pg_n * n_F_local));
      double term_diffusion_T_2 = -h_w * krW * (omega_s * rho_w_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local)) + omega_n * rho_w_n * ((Kgradu_Pg_n * n_F_local) - (Kgradu_Pc_n * n_F_local)));
      double term_diffusion_T_3 = -kth_eff * (omega_s * (gradu_T_s * n_F_local) + omega_n * (gradu_T_n * n_F_local));
      double term_diffusion_T = term_diffusion_T_1 + term_diffusion_T_2 + term_diffusion_T_3;
      double term_nipg_T = theta_T * (T_s - T_n);
      double term_penalty_T = penalty_factor_T * (T_s - T_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_diffusion_T * psi_T_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_diffusion_T * -psi_T_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_nipg_T * (kth_eff * omega_s * (n_F * gradpsi_T_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_nipg_T * (kth_eff * omega_n * (n_F * gradpsi_T_n[i])) * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_penalty_T * -psi_T_n[i] * factor);
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
    //Gas pressure
    const auto &lfsv_Pg_s = lfsv.template child<Indices::PVId_Pg>();
    const auto &lfsu_Pg_s = lfsu.template child<Indices::PVId_Pg>();

    //Capillary Pressure
    const auto &lfsv_Pc_s = lfsv.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_s = lfsu.template child<Indices::PVId_Pc>();

    //Water Saturation
    const auto &lfsv_Sw_s = lfsv.template child<Indices::PVId_Sw>();
    const auto &lfsu_Sw_s = lfsu.template child<Indices::PVId_Sw>();

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

    // define types
    using RF = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
        Traits::LocalBasisType::Traits::RangeFieldType;
    using size_type = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::SizeType;

    // dimensions
    const int dim = IG::Entity::dimension;
    const int order = std::max(lfsu_Pg_s.finiteElement().localBasis().order(),
                               lfsv_Pg_s.finiteElement().localBasis().order());

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();

    // Get geometries
    auto geo = ig.geometry();
    auto geo_inside = cell_inside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();

    // evaluate permeability tensors
    auto ref_el_inside = referenceElement(geo_inside);
    auto local_inside = ref_el_inside.position(0, 0);

    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = geo_inside.volume() / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF harmonic_average(0.0);
    omega_s = 1.;
    harmonic_average = 1.0;

    //      if (weights==ConvectionDiffusionDGWeights::weightsOn)
    //        {
    //          RF delta_s = (An_F_s*n_F);
    //          RF delta_n = (An_F_n*n_F);
    //          omega_s = delta_n/(delta_s+delta_n+1e-20);
    //          omega_n = delta_s/(delta_s+delta_n+1e-20);
    //          harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
    //        }
    //      else
    //        {
    //          omega_s = omega_n = 0.5;
    //          harmonic_average = 1.0;
    //        }

    // get polynomial degree
    auto order_s = lfsv_Pg_s.finiteElement().localBasis().order();
    auto degree = lfsv_Pg_s.finiteElement().localBasis().order();

    // penalty factor
    auto penalty_factor_g = (alpha_g / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_w = (alpha_w / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_s = (alpha_s / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_T = (alpha_T / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);
    auto penalty_factor_y = (alpha_y / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pg_s(lfsu_Pg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pg_s(lfsv_Pg_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_s(lfsu_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_s(lfsv_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XCH4_s(lfsu_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XCH4_s(lfsv_XCH4_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_YH2O_s(lfsu_YH2O_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_YH2O_s(lfsv_YH2O_s.size());

    Dune::FieldVector<RF, dim> gradu_Pg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd+quadrature_factor*order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto ip_global_s = geo_inside.global(iplocal_s);

      // evaluate basis functions

      auto bctype = paramclass.problemSpecs.ProblemBCTypes(ip_global_s);
      auto bcvalue = paramclass.problemSpecs.ProblemBCValues(ip_global_s, *(time));

      auto &psi_Pg_s = cache_Pg[order].evaluateFunction(iplocal_s, lfsv_Pg_s.finiteElement().localBasis());
      auto &psi_Sw_s = cache_Sw[order].evaluateFunction(iplocal_s, lfsv_Sw_s.finiteElement().localBasis());
      auto &psi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &psi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &psi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &psi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());

      if (bctype[Indices::PVId_Pg] == Indices::BCId_neumann)
      {

        double flux_g = bcvalue[Indices::PVId_Pg];
        for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        {
          r.accumulate(lfsv_Pg_s, i, flux_g * psi_Pg_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_Pc] == Indices::BCId_neumann)
      {

        double flux_w = bcvalue[Indices::PVId_Pc];
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, flux_w * psi_Pc_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_T] == Indices::BCId_neumann)
      {

        double flux_T = bcvalue[Indices::PVId_T];
        for (size_type i = 0; i < lfsv_T_s.size(); i++)
        {
          r.accumulate(lfsv_T_s, i, flux_T * psi_T_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_T] == Indices::BCId_neumann and
          bctype[Indices::PVId_Pc] == Indices::BCId_neumann and
          bctype[Indices::PVId_Pg] == Indices::BCId_neumann and
          bctype[Indices::PVId_XCH4] == Indices::BCId_neumann and
          bctype[Indices::PVId_YH2O] == Indices::BCId_neumann)
      {
        continue;
      }

      auto &phi_Pg_s = cache_Pg[order].evaluateFunction(iplocal_s, lfsu_Pg_s.finiteElement().localBasis());
      auto &phi_Sw_s = cache_Sw[order].evaluateFunction(iplocal_s, lfsu_Sw_s.finiteElement().localBasis());
      auto &phi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &phi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &phi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &phi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      

      // evaluate Pg
      RF Pg_s = 0.0;
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        Pg_s += x(lfsu_Pg_s, i) * phi_Pg_s[i];
      RF Pg_n = 0.0;
      if (bctype[Indices::PVId_Pg] == Indices::BCId_neumann)
      {
        Pg_n = Pg_s;
      }
      else if (bctype[Indices::PVId_Pg] == Indices::BCId_dirichlet)
      {
        Pg_n = bcvalue[Indices::PVId_Pg] / Xc_P;
      }

      // evaluate Sh
      RF Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        Sh_s += x(lfsu_Sh_s, i) * phi_Sh_s[i];
      RF Sh_n = Sh_s;

      // evaluate Sw
      RF Sw_s = 0.0;
      for (size_type i = 0; i < lfsu_Sw_s.size(); i++)
        Sw_s += x(lfsu_Sw_s, i) * phi_Sw_s[i];
      RF Sw_n = 0.0;
      if (bctype[Indices::PVId_Sw] == Indices::BCId_neumann)
      {
        Sw_n = Sw_s;
      }
      else if (bctype[Indices::PVId_Sw] == Indices::BCId_dirichlet)
      {
        Sw_n=bcvalue[Indices::PVId_Sw];
        //Sw_n = Sw_s;
      }

      // evaluate Pc
      RF Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        Pc_s += x(lfsu_Pc_s, i) * phi_Pc_s[i];
      RF Pc_n = 0.;
      if (bctype[Indices::PVId_Pc] == Indices::BCId_neumann)
      {
        Pc_n = Pc_s;
      }
      else if (bctype[Indices::PVId_Pc] == Indices::BCId_dirichlet)
      {
        //auto PcSF1_n = paramclass.hydraulicProperty.PcSF1(Sh_n);
        //Pc_n=paramclass.hydraulicProperty.suctionPressure(Sw_n,Sh_n)*PcSF1_n;
        Pc_n = bcvalue[Indices::PVId_Pc] * paramclass.hydraulicProperty.PcSF1(Sh_n) / Xc_P;
      }

      // evaluate T
      RF T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        T_s += x(lfsu_T_s, i) * phi_T_s[i];
      RF T_n = 0.0;
      if (bctype[Indices::PVId_T] == Indices::BCId_neumann)
      {
        T_n = T_s;
      }
      else if (bctype[Indices::PVId_T] == Indices::BCId_dirichlet)
      {
        T_n = bcvalue[Indices::PVId_T] / Xc_T;
      }

      // evaluate XCH4
      RF XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        XCH4_s += x(lfsu_XCH4_s, i) * phi_XCH4_s[i];
      RF XCH4_n = 0.0;
      if (bctype[Indices::PVId_XCH4] == Indices::BCId_neumann)
      {
        XCH4_n = XCH4_s;
      }
      else if (bctype[Indices::PVId_XCH4] == Indices::BCId_dirichlet)
      {
        XCH4_n = bcvalue[Indices::PVId_XCH4] / Xc_X;
      }

      // evaluate YH2O
      RF YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        YH2O_s += x(lfsu_YH2O_s, i) * phi_YH2O_s[i];
      RF YH2O_n = 0.0;
      if (bctype[Indices::PVId_YH2O] == Indices::BCId_neumann)
      {
        YH2O_n = YH2O_s;
      }
      else if (bctype[Indices::PVId_YH2O] == Indices::BCId_dirichlet)
      {
        YH2O_n = bcvalue[Indices::PVId_YH2O] / Xc_Y;
      }


      // evaluate Pw
      RF Pw_s = Pg_s - Pc_s;
      RF Pw_n = Pg_n - Pc_n;
      RF Peff_s = (Pg_s * (1. - Sw_s - Sh_s) + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * (1. - Sw_n - Sh_n) + Pw_n * Sw_n) / (1. - Sh_n);

      // evaluate gradient of basis functions
      auto &js_Pg_s = cache_Pg[order].evaluateJacobian(iplocal_s, lfsu_Pg_s.finiteElement().localBasis());
      auto &js_v_Pg_s = cache_Pg[order].evaluateJacobian(iplocal_s, lfsv_Pg_s.finiteElement().localBasis());
      auto &js_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &js_v_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      

      // transform gradients of shape functions to real element

      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        jac.mv(js_Pg_s[i][0], gradphi_Pg_s[i]);
      for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        jac.mv(js_v_Pg_s[i][0], gradpsi_Pg_s[i]);
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        jac.mv(js_Pc_s[i][0], gradphi_Pc_s[i]);
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        jac.mv(js_v_Pc_s[i][0], gradpsi_Pc_s[i]);
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

      // compute gradient of Pg
      gradu_Pg_s = 0.0;
      for (size_type i = 0; i < lfsu_Pg_s.size(); i++)
        gradu_Pg_s.axpy(x(lfsu_Pg_s, i), gradphi_Pg_s[i]);

      // compute gradient of Pc
      gradu_Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        gradu_Pc_s.axpy(x(lfsu_Pc_s, i), gradphi_Pc_s[i]);

      // compute gradient of T
      gradu_T_s = 0.0;
      for (size_type i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x(lfsu_T_s, i), gradphi_T_s[i]);

      // compute gradient of XCH4
      gradu_XCH4_s = 0.0;
      for (size_type i = 0; i < lfsu_XCH4_s.size(); i++)
        gradu_XCH4_s.axpy(x(lfsu_XCH4_s, i), gradphi_XCH4_s[i]);

      // compute gradient of YH2O
      gradu_YH2O_s = 0.0;
      for (size_type i = 0; i < lfsu_YH2O_s.size(); i++)
        gradu_YH2O_s.axpy(x(lfsu_YH2O_s, i), gradphi_YH2O_s[i]);

      auto K_s = paramclass.problemSpecs.SedimentPermeabilityTensor(ip_global_s);
      K_s *= 1. / Xc_K;
      K_s *= Xc_conv_m;
      auto n_F = ig.centerUnitOuterNormal();
      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F, Kn_F_s);

      // compute K * gradient of Pg
      K_s.mv(gradu_Pg_s, Kgradu_Pg_s);

      // compute K * gradient of Pc
      K_s.mv(gradu_Pc_s, Kgradu_Pc_s);

      for (int i = 0; i < paramclass.problemSpecs.dimension; i++)
      {
        v_g[i] = omega_s * Kgradu_Pg_s[i];
        v_w[i] = omega_s * (Kgradu_Pg_s[i] - Kgradu_Pc_s[i]);
      }

      double normalflux_g = -1. * (v_g * n_F_local);
      double normalflux_w = -1. * (v_w * n_F_local);

      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g >= 0.0)
      {
        omegaup_g_s = 1.0;
        omegaup_g_n = 0.0;
      }
      else
      {
        omegaup_g_s = 0.0;
        omegaup_g_n = 1.0;
      }

      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w >= 0.0)
      {
        omegaup_w_s = 1.0;
        omegaup_w_n = 0.0;
      }
      else
      {
        omegaup_w_s = 0.0;
        omegaup_w_n = 1.0;
      }

      auto rho_g_s = paramclass.methane.density(T_s * Xc_T, Pg_s * Xc_P, 1.) / Xc_rho;
      auto rho_w_s = paramclass.water.density(T_s * Xc_T, Pw_s * Xc_P) / Xc_rho;
      auto por_s = paramclass.problemSpecs.SedimentPorosity(iplocal_s);
      auto krW_s = paramclass.hydraulicProperty.krW(Sw_s, Sh_s) / (paramclass.water.dynamicViscosity(T_s * Xc_T, Pw_s * Xc_P) / Xc_mu);
      auto krN_s = paramclass.hydraulicProperty.krNW(Sw_s, Sh_s) / (paramclass.methane.dynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) / Xc_mu);
      //          std::cout<< "densities: (g,w)_s " << rho_g_s <<", "<< rho_w_s<<std::endl;
      //          std::cout<<"lambdas: (w,n)_s " <<krW_s<<", "<<krN_s<<std::endl;
      //  adding terms regarding components
      auto Sg_s = 1. - Sw_s - Sh_s;
      auto tau_s = paramclass.soil.tortuosity(por_s);
      auto DH2O_g_s = tau_s * por_s * paramclass.mixture.binaryDiffCoeffInGas(T_s * Xc_T, Pg_s * Xc_P);
      auto DCH4_w_s = tau_s * por_s * paramclass.mixture.binaryDiffCoeffInLiquid(T_s * Xc_T, Pw_s * Xc_P);
      //auto XCH4_s = paramclass.mixture.mole_x_CH4(T_s * Xc_T, Pw_s * Xc_P);
      auto YCH4_s = paramclass.mixture.mole_y_CH4(T_s * Xc_T, Pw_s * Xc_P);
      auto XH2O_s = paramclass.mixture.mole_x_H2O(T_s * Xc_T, Pw_s * Xc_P);
      //auto YH2O_s = paramclass.mixture.mole_y_H2O(T_s * Xc_T, Pw_s * Xc_P);
      auto H_M_w_s = paramclass.methane.henrysConstant(T_s * Xc_T);
      auto P_H_sat_s = paramclass.water.saturatedVaporPressure(T_s * Xc_T);
      auto zCH4_s = paramclass.eos.evaluateCompressibilityFactor(T_s * Xc_T, Pw_s * Xc_P);
      //std::cout << "----" << zCH4 << std::endl;
      //  end of terms regarding components

      auto Cp_g_s = paramclass.methane.Cp(T_s * Xc_T, Pg_s * Xc_P, 1.) / Xc_C;
      auto Cp_w_s = paramclass.water.Cp(T_s * Xc_T, Pw_s * Xc_P) / Xc_C;
      auto kth_g_s = paramclass.methane.thermalConductivity(T_s * Xc_T, Pg_s * Xc_P) / Xc_kth;
      auto kth_w_s = paramclass.water.thermalConductivity(T_s * Xc_T, Pw_s * Xc_P) / Xc_kth;
      auto kth_h_s = paramclass.hydrate.thermalConductivity(T_s * Xc_T, Peff_s * Xc_P) / Xc_kth;
      auto kth_s_s = paramclass.soil.thermalConductivity(T_s * Xc_T, Peff_s * Xc_P) / Xc_kth;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * ((1. - Sw_s - Sh_s) * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s);

      auto rho_g_n = paramclass.methane.density(T_n * Xc_T, Pg_n * Xc_P, 1.) / Xc_rho;
      auto rho_w_n = paramclass.water.density(T_n * Xc_T, Pw_n * Xc_P) / Xc_rho;
      auto por_n = paramclass.problemSpecs.SedimentPorosity(iplocal_s);
      auto krW_n = paramclass.hydraulicProperty.krW(Sw_n, Sh_n) / (paramclass.water.dynamicViscosity(T_n * Xc_T, Pw_n * Xc_P) / Xc_mu);
      auto krN_n = paramclass.hydraulicProperty.krNW(Sw_n, Sh_n) / (paramclass.methane.dynamicViscosity(T_n * Xc_T, Pg_n * Xc_P) / Xc_mu);
      //          std::cout<< "densities: (g,w)_n " << rho_g_n <<", "<< rho_w_n<<std::endl;
      //          std::cout<<"lambdas: (w,n)_n " <<krW_n<<", "<<krN_n<<std::endl;
      //  adding terms regarding components
      auto Sg_n = 1. - Sw_n - Sh_n;
      auto tau_n = paramclass.soil.tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * paramclass.mixture.binaryDiffCoeffInGas(T_n * Xc_T, Pg_n * Xc_P);
      auto DCH4_w_n = tau_n * por_n * paramclass.mixture.binaryDiffCoeffInLiquid(T_n * Xc_T, Pg_n * Xc_P);
      //auto XCH4_n = paramclass.mixture.mole_x_CH4(T_n * Xc_T, Pg_n * Xc_P);
      auto YCH4_n = paramclass.mixture.mole_y_CH4(T_n * Xc_T, Pg_n * Xc_P);
      auto XH2O_n = paramclass.mixture.mole_x_H2O(T_n * Xc_T, Pg_n * Xc_P);
      //auto YH2O_n = paramclass.mixture.mole_y_H2O(T_n * Xc_T, Pg_n * Xc_P);
      auto H_M_w_n = paramclass.methane.henrysConstant(T_n * Xc_T);
      auto P_H_sat_n = paramclass.water.saturatedVaporPressure(T_n * Xc_T);
      auto zCH4_n = paramclass.eos.evaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);
      //std::cout << "----" << zCH4 << std::endl;
      //  end of terms regarding components
      auto Cp_g_n = paramclass.methane.Cp(T_n * Xc_T, Pg_n * Xc_P, 1.) / Xc_C;
      auto Cp_w_n = paramclass.water.Cp(T_n * Xc_T, Pw_n * Xc_P) / Xc_C;
      auto kth_g_n = paramclass.methane.thermalConductivity(T_n * Xc_T, Pg_n * Xc_P) / Xc_kth;
      auto kth_w_n = paramclass.water.thermalConductivity(T_n * Xc_T, Pw_n * Xc_P) / Xc_kth;
      auto kth_h_n = paramclass.hydrate.thermalConductivity(T_n * Xc_T, Peff_n * Xc_P) / Xc_kth;
      auto kth_s_n = paramclass.soil.thermalConductivity(T_n * Xc_T, Peff_n * Xc_P) / Xc_kth;
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * ((1. - Sw_n - Sh_n) * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);

      auto krN = omegaup_g_s * krN_s + omegaup_g_n * krN_n;
      auto krW = omegaup_w_s * krW_s + omegaup_w_n * krW_n;
      auto h_g = omegaup_g_s * Cp_g_s * T_s + omegaup_g_n * Cp_g_n * T_n;
      auto h_w = omegaup_w_s * Cp_w_s * T_s + omegaup_w_n * Cp_w_n * T_n;
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
      kth_eff *= Xc_diff_h;

      auto DH2O_g = omega_s * DH2O_g_s + omega_s * DH2O_g_n; // = DCH4_g
      auto DCH4_w = omega_s * DCH4_w_s + omega_s * DCH4_w_n; // = DH2O_w


      if (bctype[Indices::PVId_Pg] == Indices::BCId_dirichlet)
      {
        //Methane
        //auto DH2O_g = omega_s * DH2O_g_s + omega_n * DH2O_g_n; // = DCH4_g
        auto Diffcoeff_m_g_p_s = rho_g_s * Sg_s * -YCH4_s / Pg_s;
        auto Diffcoeff_m_g_x_s = rho_g_s * Sg_s * H_M_w_s / ( zCH4_s * Pg_s);
      
        //auto DCH4_w = omega_s * DCH4_w_s + omega_n * DCH4_w_n; // = DH2O_w
        auto Diffcoeff_m_w_s = rho_w_s * Sw_s ;

        double term_conv_m_g = -krN * omega_s * rho_g_s * YCH4_s * (Kgradu_Pg_s * n_F_local); 
        double term_conv_m_w = -krW * omega_s * rho_w_s * XCH4_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local)); 

        double term_diffusion_m = term_conv_m_g + term_conv_m_w + DH2O_g * (omega_s * (Diffcoeff_m_g_p_s * (gradu_Pg_s * n_F_local) 
                                              + Diffcoeff_m_g_x_s * (gradu_XCH4_s * n_F_local)))
                                              + DCH4_w * omega_s * Diffcoeff_m_w_s * (gradu_XCH4_s * n_F_local);
                                             
        double term_nipg_g = theta_g * (Pg_s - Pg_n);
        double term_nipg_w = theta_w * ((Pg_s - Pc_s) - (Pg_n - Pc_n));
        double term_nipg_m_x = theta_x * (XCH4_s - XCH4_n);

        double term_penalty_g = penalty_factor_g * (Pg_s - Pg_n);
        
        // diffusion term
        for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        {
          r.accumulate(lfsv_Pg_s, i, term_diffusion_m * psi_Pg_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        {
          r.accumulate(lfsv_Pg_s, i, (term_nipg_g * krN * omega_s * rho_g_s * YCH4_s * (Kn_F_s * gradpsi_Pg_s[i]) 
                                      + term_nipg_w * krW * omega_s * rho_w_s * XCH4_s  * (Kn_F_s * gradpsi_Pg_s[i])
                                      + omega_s * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_s 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_s
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_s) * (n_F_local* gradpsi_Pg_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Pg_s.size(); i++)
        {
          r.accumulate(lfsv_Pg_s, i, term_penalty_g * psi_Pg_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_Pc] == Indices::BCId_dirichlet)
      {
        //auto DH2O_g = omega_s * DH2O_g_s + omega_n * DH2O_g_n; // = DCH4_g
        auto Diffcoeff_w_w_p_s = rho_w_s * Sw_s * YH2O_s / P_H_sat_s;
        auto Diffcoeff_w_w_y_s = rho_w_s * Sw_s * Pg_s / P_H_sat_s;
      
        //auto DCH4_w = omega_s * DCH4_w_s + omega_n * DCH4_w_n; // = DH2O_w
        auto Diffcoeff_w_g_s = rho_g_s * Sg_s ;
     

        double term_conv_w_g = -krN *  omega_s * rho_g_s * YH2O_s * (Kgradu_Pg_s * n_F_local);

        double term_conv_w_w = -krW * omega_s * rho_w_s * XH2O_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local));

        double term_diffusion_w = term_conv_w_g + term_conv_w_w + DH2O_g * omega_s * Diffcoeff_w_g_s * (gradu_YH2O_s * n_F_local) 
                                + DCH4_w * (omega_s * (Diffcoeff_w_w_p_s * (gradu_Pg_s * n_F_local)
                                + Diffcoeff_w_w_y_s * (gradu_YH2O_s * n_F_local)));


        double term_nipg_g = theta_g * (Pg_s - Pg_n);
        double term_nipg_w = theta_w * ((Pg_s - Pc_s) - (Pg_n - Pc_n));
        
        double term_nipg_w_y = theta_y * (YH2O_s - YH2O_n);
        double term_penalty_w = penalty_factor_w * ((Pg_s - Pc_s) - (Pg_n - Pc_n));
        
        // diffusion term
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, term_diffusion_w * psi_Pc_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, (term_nipg_g * krN * omega_s * rho_g_s * YH2O_s * (Kn_F_s * gradpsi_Pc_s[i]) 
                                      + term_nipg_w * krW * omega_s * rho_w_s * XH2O_s  * (Kn_F_s * gradpsi_Pc_s[i])
                                      + omega_s * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_s 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_s
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_s)  * (n_F_local* gradpsi_Pc_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, term_penalty_w * psi_Pc_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_Sw] == Indices::BCId_dirichlet)
      {

        double term_penalty_sw = penalty_factor_s * (Sw_s - Sw_n);
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Sw_s.size(); i++)
        {
          r.accumulate(lfsv_Sw_s, i, term_penalty_sw * psi_Sw_s[i] * factor);
        }
      }

      double term_penalty_sh = penalty_factor_s * (Sh_s - Sh_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sh_s.size(); i++)
      {
        r.accumulate(lfsv_Sh_s, i, term_penalty_sh * psi_Sh_s[i] * factor);
      }

      double term_penalty_XCH4 = penalty_factor_x * (XCH4_s - XCH4_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_XCH4_s.size(); i++)
      {
        r.accumulate(lfsv_XCH4_s, i, term_penalty_XCH4 * psi_XCH4_s[i] * factor);
      }
      
      double term_penalty_YH2O = penalty_factor_y * (YH2O_s - YH2O_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_YH2O_s.size(); i++)
      {
        r.accumulate(lfsv_YH2O_s, i, term_penalty_YH2O * psi_YH2O_s[i] * factor);
      }
      

      if (bctype[Indices::PVId_T] == Indices::BCId_dirichlet)
      {

        //        	  double term_diffusion_T_1 = - h_g * krN * ( omega_s*(Kgradu_Pg_s*n_F_local) );
        double term_diffusion_T_1 = -h_g * krN * (omega_s * rho_g_s * (Kgradu_Pg_s * n_F_local));
        //              double term_diffusion_T_2 = - h_w * krW * ( omega_s*( (Kgradu_Pg_s*n_F_local) - (Kgradu_Pc_s*n_F_local) ) );
        double term_diffusion_T_2 = -h_w * krW * (omega_s * rho_w_s * ((Kgradu_Pg_s * n_F_local) - (Kgradu_Pc_s * n_F_local)));
        if (bctype[Indices::PVId_Pg] == Indices::BCId_neumann)
        {
          term_diffusion_T_1 = -bcvalue[Indices::PVId_Pg] * h_g;
        }
        if (bctype[Indices::PVId_Pc] == Indices::BCId_neumann)
        {
          term_diffusion_T_2 = -bcvalue[Indices::PVId_Pc] * h_w;
        }
        double term_diffusion_T_3 = -kth_eff * (omega_s * (gradu_T_s * n_F_local));
        double term_diffusion_T = term_diffusion_T_1 + term_diffusion_T_2 + term_diffusion_T_3;
        double term_nipg_T = theta_T * (T_s - T_n);
        double term_penalty_T = penalty_factor_T * (T_s - T_n);

        // diffusion term
        for (size_type i = 0; i < lfsv_T_s.size(); i++)
        {
          r.accumulate(lfsv_T_s, i, term_diffusion_T * psi_T_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_T_s.size(); i++)
        {
          r.accumulate(lfsv_T_s, i, term_nipg_T * (kth_eff * omega_s * (n_F * gradpsi_T_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_T_s.size(); i++)
        {
          r.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
        }
      }

    } //End Quadrature Rule
  }

  //  template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
  //  void jacobian_boundary (const IG& ig,
  //                          const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
  //                          M& mat_ss) const
  //  {
  //      // define types
  //      using RF = typename LFSV::Traits::FiniteElementType::
  //        Traits::LocalBasisType::Traits::RangeFieldType;
  //      using size_type = typename LFSV::Traits::SizeType;
  //
  //      // dimensions
  //      const int dim = IG::dimension;
  //      const int order = std::max(
  //          lfsu_s.finiteElement().localBasis().order(),
  //          lfsv_s.finiteElement().localBasis().order()
  //          );
  //
  //      // References to the inside cell
  //      const auto& cell_inside = ig.inside();
  //
  //      // Get geometries
  //      auto geo = ig.geometry();
  //      auto geo_inside = cell_inside.geometry();
  //
  //      // Get geometry of intersection in local coordinates of cell_inside
  //      auto geo_in_inside = ig.geometryInInside();
  //
  //      // evaluate permeability tensors
  //      auto ref_el_inside = referenceElement(geo_inside);
  //      auto local_inside = ref_el_inside.position(0,0);
  //      auto A_s = param.A(cell_inside,local_inside);
  //
  //      // face diameter
  //      auto h_F = geo_inside.volume()/geo.volume(); // Houston!
  //
  //      // compute weights
  //      auto n_F = ig.centerUnitOuterNormal();
  //      Dune::FieldVector<RF,dim> An_F_s;
  //      A_s.mv(n_F,An_F_s);
  //      RF harmonic_average;
  //      if (weights==ConvectionDiffusionDGWeights::weightsOn)
  //        harmonic_average = An_F_s*n_F;
  //      else
  //        harmonic_average = 1.0;
  //
  //      // get polynomial degree
  //      auto order_s = lfsu_s.finiteElement().localBasis().order();
  //      auto degree = order_s;
  //
  //      // penalty factor
  //      auto penalty_factor = (alpha/h_F) * harmonic_average * degree*(degree+dim-1);
  //
  //      // Initialize vectors outside for loop
  //      std::vector<Dune::FieldVector<RF,dim> > tgradphi_s(lfsu_s.size());
  //
  //      // Transformation matrix
  //      typename IG::Entity::Geometry::JacobianInverseTransposed jac;
  //
  //      // loop over quadrature points
  //      auto intorder = intorderadd+quadrature_factor*order;
  //      for (const auto& ip : quadratureRule(geo,intorder))
  //        {
  //          auto bctype = param.bctype(ig.intersection(),ip.position());
  //
  //          if (bctype == ConvectionDiffusionBoundaryConditions::None ||
  //              bctype == ConvectionDiffusionBoundaryConditions::Neumann)
  //            continue;
  //
  //          // position of quadrature point in local coordinates of elements
  //          auto iplocal_s = geo_in_inside.global(ip.position());
  //
  //          // local normal
  //          auto n_F_local = ig.unitOuterNormal(ip.position());
  //
  //          // evaluate basis functions
  //          auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,lfsu_s.finiteElement().localBasis());
  //
  //          // integration factor
  //          auto factor = ip.weight() * geo.integrationElement(ip.position());
  //
  //          // evaluate velocity field and upwinding, assume H(div) velocity field => choose any side
  //          auto b = param.b(cell_inside,iplocal_s);
  //          auto normalflux = b*n_F_local;
  //
  //          if (bctype == ConvectionDiffusionBoundaryConditions::Outflow)
  //            {
  //              if (normalflux<-1e-30)
  //                DUNE_THROW(Dune::Exception,
  //                  "Outflow boundary condition on inflow! [b("
  //                  << geo.global(ip.position()) << ") = "
  //                  << b << ")" << n_F_local << " " << normalflux);
  //
  //              // convection term
  //              for (size_type j=0; j<lfsu_s.size(); j++)
  //                for (size_type i=0; i<lfsu_s.size(); i++)
  //                  mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * normalflux * factor * phi_s[i]);
  //
  //              continue;
  //            }
  //
  //          // evaluate gradient of basis functions
  //          auto& gradphi_s = cache[order_s].evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
  //
  //          // transform gradients of shape functions to real element
  //          jac = geo_inside.jacobianInverseTransposed(iplocal_s);
  //          for (size_type i=0; i<lfsu_s.size(); i++) jac.mv(gradphi_s[i][0],tgradphi_s[i]);
  //
  //          // upwind
  //          RF omegaup_s, omegaup_n;
  //          if (normalflux>=0.0)
  //            {
  //              omegaup_s = 1.0;
  //              omegaup_n = 0.0;
  //            }
  //          else
  //            {
  //              omegaup_s = 0.0;
  //              omegaup_n = 1.0;
  //            }
  //
  //          // convection term
  //          for (size_type j=0; j<lfsu_s.size(); j++)
  //            for (size_type i=0; i<lfsu_s.size(); i++)
  //              mat_ss.accumulate(lfsu_s,i,lfsu_s,j,omegaup_s * phi_s[j] * normalflux * factor * phi_s[i]);
  //
  //          // diffusion term
  //          for (size_type j=0; j<lfsu_s.size(); j++)
  //            for (size_type i=0; i<lfsu_s.size(); i++)
  //              mat_ss.accumulate(lfsu_s,i,lfsu_s,j,-(An_F_s*tgradphi_s[j]) * factor * phi_s[i]);
  //
  //          // (non-)symmetric IP term
  //          for (size_type j=0; j<lfsu_s.size(); j++)
  //            for (size_type i=0; i<lfsu_s.size(); i++)
  //              mat_ss.accumulate(lfsu_s,i,lfsu_s,j,phi_s[j] * factor * theta * (An_F_s*tgradphi_s[i]));
  //
  //          // standard IP term
  //          for (size_type j=0; j<lfsu_s.size(); j++)
  //            for (size_type i=0; i<lfsu_s.size(); i++)
  //              mat_ss.accumulate(lfsu_s,i,lfsu_s,j,penalty_factor * phi_s[j] * phi_s[i] * factor);
  //        }
  //    }
};

#endif /* FLOW_LOCALOPERATOR_HH_ */
