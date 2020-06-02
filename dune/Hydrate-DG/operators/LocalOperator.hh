/*
 * LocalOperator.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef FLOW_LOCALOPERATOR_HH_
#define FLOW_LOCALOPERATOR_HH_

using namespace Dune::PDELab;

struct ConvectionDiffusionDGMethod
{
  enum Type { NIPG, SIPG, IIPG };
};

template <typename GV, typename Params, typename U, class GFS, class FEM_P, class FEM_S, class FEM_T, class FEM_X, class FEM_Y>
class LocalOperator : public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianVolume<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::NumericalJacobianBoundary<LocalOperator<GV, Params, U, GFS, FEM_P, FEM_S, FEM_T, FEM_X, FEM_Y>>,
                           public Dune::PDELab::FullSkeletonPattern, // matrix entries skeleton
                           public Dune::PDELab::FullVolumePattern,
                           public Dune::PDELab::LocalOperatorDefaultFlags,
                           public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
  
	
  const GV &gv;
  const Params&	  property;
  typedef ProblemBoundaryConditions<GV,Params> BC ;
	

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
  double Xc_t;

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
  using LocalBasisType_Pc = typename FEM_P::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Pc = Dune::PDELab::LocalBasisCache<LocalBasisType_Pc>;
  using LocalBasisType_Sg = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sg = Dune::PDELab::LocalBasisCache<LocalBasisType_Sg>;
  using LocalBasisType_Sh = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_Sh = Dune::PDELab::LocalBasisCache<LocalBasisType_Sh>;
  using LocalBasisType_T = typename FEM_T::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_T = Dune::PDELab::LocalBasisCache<LocalBasisType_T>;
  using LocalBasisType_XCH4 = typename FEM_X::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_XCH4 = Dune::PDELab::LocalBasisCache<LocalBasisType_XCH4>;
  using LocalBasisType_YH2O = typename FEM_Y::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_YH2O = Dune::PDELab::LocalBasisCache<LocalBasisType_YH2O>;
  using LocalBasisType_XC = typename FEM_X::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_XC = Dune::PDELab::LocalBasisCache<LocalBasisType_XC>;
  

  // In theory it is possible that one and the same local operator is
  // called first with a finite element of one type and later with a
  // finite element of another type.  Since finite elements of different
  // type will usually produce different results for the same local
  // coordinate they cannot share a cache.  Here we use a vector of caches
  // to allow for different orders of the shape functions, which should be
  // enough to support p-adaptivity.  (Another likely candidate would be
  // differing geometry types, i.e. hybrid meshes.)

  std::vector<Cache_Pw> cache_Pw;
  std::vector<Cache_Pc> cache_Pc;
  std::vector<Cache_Sg> cache_Sg;
  std::vector<Cache_Sh> cache_Sh;
  std::vector<Cache_T> cache_T;
  std::vector<Cache_XCH4> cache_XCH4;
  std::vector<Cache_YH2O> cache_YH2O;
  std::vector<Cache_XC> cache_XC;

  // constructor stores parameters
  LocalOperator(const GV &gv_, const Params&	 property_,
                     U *unew_,
                     GFS gfs_,
                     double *time_,
                     double *dt_,
                     unsigned int intorder_ = 6,
                     ConvectionDiffusionDGMethod::Type method_g_ = ConvectionDiffusionDGMethod::IIPG,
                     ConvectionDiffusionDGMethod::Type method_w_ = ConvectionDiffusionDGMethod::IIPG,
                     ConvectionDiffusionDGMethod::Type method_T_ = ConvectionDiffusionDGMethod::IIPG,
                     ConvectionDiffusionDGMethod::Type method_x_ = ConvectionDiffusionDGMethod::IIPG,
                     ConvectionDiffusionDGMethod::Type method_y_ = ConvectionDiffusionDGMethod::IIPG,
                     double alpha_g_ = 1., double alpha_w_ = 1., double alpha_s_ = 1., double alpha_T_ = 1., double alpha_x_ = 1., double alpha_y_ = 1.)
      : gv(gv_), property( property_ ),
        unew(unew_),
        gfs(gfs_),
        time(time_),
        dt(dt_),
        intorder(intorder_),
        method_g(method_g_), method_w(method_w_), method_T(method_T_), method_x(method_x_), method_y(method_y_),
        alpha_g(alpha_g_), alpha_w(alpha_w_), alpha_s(alpha_s_),  alpha_T(alpha_T_), alpha_x(alpha_x_), alpha_y(alpha_y_),
        cache_Pw(20), cache_Pc(20),  cache_T(20), cache_Sg(20), cache_Sh(20), cache_XCH4(20), cache_YH2O(20), cache_XC(20)
  {
    theta_g = 0.0;
    if (method_g == ConvectionDiffusionDGMethod::SIPG)
      theta_g = -1.0;
    if (method_g == ConvectionDiffusionDGMethod::NIPG)
      theta_g = 1.0;

    theta_w = 0.0;
    if (method_w == ConvectionDiffusionDGMethod::SIPG)
      theta_w = -1.0;
    if (method_w == ConvectionDiffusionDGMethod::NIPG)
      theta_w = 1.0;

    theta_T = 0.0;
    if (method_T == ConvectionDiffusionDGMethod::SIPG)
      theta_T = -1.0;
    if (method_T == ConvectionDiffusionDGMethod::NIPG)
      theta_T = 1.0;

    theta_x = 0.0;
    if (method_x == ConvectionDiffusionDGMethod::SIPG)
      theta_x = -1.0;
    if (method_w == ConvectionDiffusionDGMethod::NIPG)
      theta_x = 1.0;

    theta_y = 0.0;
    if (method_y == ConvectionDiffusionDGMethod::SIPG)
      theta_y = -1.0;
    if (method_y == ConvectionDiffusionDGMethod::NIPG)
      theta_y = 1.0;

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
    Xc_P = property.characteristicValue.P_c;
    Xc_T = property.characteristicValue.T_c;
    Xc_X = property.characteristicValue.x_c;
    Xc_Y = property.characteristicValue.x_c;
    Xc_t = property.characteristicValue.t_c;
  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
  {

    // subspaces
    //Gas pressure
    const auto &lfsv_Pw = lfsv.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw = lfsu.template child<Indices::PVId_Pw>();

    //Capillary Pressure
    const auto &lfsv_Pc = lfsv.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc = lfsu.template child<Indices::PVId_Pc>();

    //Water Saturation
    const auto &lfsv_Sg = lfsv.template child<Indices::PVId_Sg>();
    const auto &lfsu_Sg = lfsu.template child<Indices::PVId_Sg>();

    //Hydrate Saturation
    const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
    const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

    // //Temperature
    const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
    const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

    //Methane mole fraction
    const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    //H2O mole fraction
    const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();

    //Salt mole fraction
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
    const int order = std::max(lfsu_Pw.finiteElement().localBasis().order(),
                               lfsv_Pw.finiteElement().localBasis().order());

    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pw(lfsu_Pw.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pw(lfsv_Pw.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc(lfsu_Pc.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc(lfsv_Pc.size());
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
    Dune::FieldVector<RF, dim> gradu_Pc(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh(0.0);
    Dune::FieldVector<RF, dim> gradu_T(0.0);
    Dune::FieldVector<RF, dim> Ktgradu_T(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O(0.0);
    Dune::FieldVector<RF, dim> gradu_XC(0.0);

    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_Pw = cache_Pw[order].evaluateFunction(ip.position(), lfsu_Pw.finiteElement().localBasis());
      auto &psi_Pw = cache_Pw[order].evaluateFunction(ip.position(), lfsv_Pw.finiteElement().localBasis());
      auto &phi_Sg = cache_Sg[order].evaluateFunction(ip.position(), lfsu_Sg.finiteElement().localBasis());
      auto &psi_Sg = cache_Sg[order].evaluateFunction(ip.position(), lfsv_Sg.finiteElement().localBasis());
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
      auto &phi_XC = cache_XC[order].evaluateFunction(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &psi_XC = cache_XC[order].evaluateFunction(ip.position(), lfsv_XC.finiteElement().localBasis());

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
        //std::cout<< " Sg = " << phi_Sg[i] << std::endl;
       
      }
      //exit(0);
      // evaluate Sh
      RF Sh = 0.0;
      for (size_type i = 0; i < lfsu_Sh.size(); i++){
        Sh += x(lfsu_Sh, i) * phi_Sh[i];
      }
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

      // evaluate XC
      RF XC = 0.0;
      for (size_type i = 0; i < lfsu_XC.size(); i++)
        XC += x(lfsu_XC, i) * phi_XC[i];

      RF Sw = 1. - Sg - Sh;
      // evaluate Pw
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por = property.soil.SedimentPorosity(cell, ip_local);
      auto suctionPressure = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      auto PcSF1 = property.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]);
      
      //auto Pc = suctionPressure * PcSF1;
      RF Pg = Pw + Pc;
      RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);

      // evaluate gradient of basis functions
      auto &js_Pw = cache_Pw[order].evaluateJacobian(ip.position(), lfsu_Pw.finiteElement().localBasis());
      auto &js_v_Pw = cache_Pw[order].evaluateJacobian(ip.position(), lfsv_Pw.finiteElement().localBasis());
      auto &js_Pc = cache_Pc[order].evaluateJacobian(ip.position(), lfsu_Pc.finiteElement().localBasis());
      auto &js_v_Pc = cache_Pc[order].evaluateJacobian(ip.position(), lfsv_Pc.finiteElement().localBasis());
      auto &js_Sg = cache_Sg[order].evaluateJacobian(ip.position(), lfsu_Sg.finiteElement().localBasis());
      auto &js_v_Sg = cache_Sg[order].evaluateJacobian(ip.position(), lfsv_Sg.finiteElement().localBasis());
      auto &js_Sh = cache_Sh[order].evaluateJacobian(ip.position(), lfsu_Sh.finiteElement().localBasis());
      auto &js_v_Sh = cache_Sh[order].evaluateJacobian(ip.position(), lfsv_Sh.finiteElement().localBasis());
      auto &js_T = cache_T[order].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());
      auto &js_XCH4 = cache_XCH4[order].evaluateJacobian(ip.position(), lfsu_XCH4.finiteElement().localBasis());
      auto &js_v_XCH4 = cache_XCH4[order].evaluateJacobian(ip.position(), lfsv_XCH4.finiteElement().localBasis());
      auto &js_YH2O = cache_YH2O[order].evaluateJacobian(ip.position(), lfsu_YH2O.finiteElement().localBasis());
      auto &js_v_YH2O = cache_YH2O[order].evaluateJacobian(ip.position(), lfsv_YH2O.finiteElement().localBasis());
      auto &js_XC = cache_XC[order].evaluateJacobian(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &js_v_XC = cache_XC[order].evaluateJacobian(ip.position(), lfsv_XC.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (size_type i = 0; i < lfsu_Pw.size(); i++)
        jac.mv(js_Pw[i][0], gradphi_Pw[i]);
      for (size_type i = 0; i < lfsv_Pw.size(); i++)
        jac.mv(js_v_Pw[i][0], gradpsi_Pw[i]);

      for (size_type i = 0; i < lfsu_Pc.size(); i++)
        jac.mv(js_Pc[i][0], gradphi_Pc[i]);
      for (size_type i = 0; i < lfsv_Pc.size(); i++)
        jac.mv(js_v_Pc[i][0], gradpsi_Pc[i]);

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

      // compute gradient of Pc
      gradu_Pc = 0.0;
      for (size_type i = 0; i < lfsu_Pc.size(); i++)
        gradu_Pc.axpy(x(lfsu_Pc, i), gradphi_Pc[i]);

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


      auto K = property.soil.SedimentPermeabilityTensor(cell, ip_local); /*ndim K from soil.hh*/
      //K *= 1. / Xc_K; /*ndim K*/
      //K *= Xc_conv_m;

      // compute K * gradient of Pw
      K.mv(gradu_Pw, Kgradu_Pw);

      // compute K * gradient of Pc
      K.mv(gradu_Pc, Kgradu_Pc);

      // compute K * gradient of Sg
      K.mv(gradu_Sg, Kgradu_Sg);

      // compute K * gradient of Sh
      K.mv(gradu_Sh, Kgradu_Sh);
      
      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );
      //auto rho_h = property.hydrate.Density() / Xc_rho;

      //  adding terms regarding components
      
      auto tau = property.soil.Tortuosity(por);/*ndim tau from soil.hh*/
      auto DH2O_g = tau * por * property.mixture.DiffCoeffH2OInGas(T * Xc_T, Pg * Xc_P); /*ndim D from mixture.hh*/
      auto DCH4_w = tau * por * property.mixture.DiffCoeffCH4InLiquid(T * Xc_T, Pw * Xc_P); /*ndim D from mixture.hh*/
      auto DC_w = tau * por * property.salt.DiffCoeff(T * Xc_T, Pw * Xc_P); /*ndim D from salt.hh*/
      
      double S = XC * (property.salt.MolarMass()/property.water.MolarMass());
      auto H_M_w = property.gas.SolubilityCoefficient(T * Xc_T, S ); /*ndim H_M_w from CH4.hh*/
      auto P_H_sat = property.water.SaturatedVaporPressure(T * Xc_T, S ); /*ndim P_H_sat from H2O.hh*/
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      //  end of terms regarding components
      auto VLequil = property.mixture.EquilibriumMoleFractions( T * Xc_T, Pg * Xc_P, XC, zCH4);
      auto YCH4 = VLequil[Indices::compId_YCH4];//property.mixture.YCH4(XCH4, T * Xc_T, Pg * Xc_P, XC, zCH4);
      auto XH2O = VLequil[Indices::compId_XH2O];//property.mixture.XH2O(YH2O, T * Xc_T, Pg * Xc_P, XC);
      
      auto Swe = property.hydraulicProperty.EffectiveSw(Sw,Sh,0.0,0.0);
      auto dPc_dSwe = property.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]);
      auto dSwe_dSg = - property.hydraulicProperty.dSwe_dSw(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sg = dPc_dSwe * dSwe_dSg ;

      auto dPcSF1_dSh =  property.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh = property.hydraulicProperty.dSwe_dSh(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sh = dPcSF1_dSh + dPc_dSwe * dSwe_dSh ;

      auto rho_g = property.gas.Density(T * Xc_T, Pg * Xc_P, zCH4); /*ndim density from CH4.hh; the input arguments are dimensional   */
      auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S); /*ndim density from H2O.hh; the input arguments are dimensional*/
      
      auto krW = property.hydraulicProperty.krw(cell, ip_local, Sw, Sh) / (property.water.DynamicViscosity(T * Xc_T, Pw * Xc_P, S) );
      auto krN = property.hydraulicProperty.krg(cell, ip_local, Sw, Sh) / (property.gas.DynamicViscosity(T * Xc_T, Pg * Xc_P));
      
      // compute source terms
			auto q_g  = property.kinetics.GasGenerationRate( T*Xc_T,
														    Pg*Xc_P,
														    Sh,
														    Sw,
														    XCH4,
														    zCH4,
														    S,
														    por,
														    permeability*Xc_K); /*[kg/m³s]*/
			auto q_w  = property.kinetics.WaterGenerationRate( q_g ); /*[kg/m³s]*/
			auto q_h  = property.kinetics.HydrateDissociationRate( q_g ); /*[kg/m³s]*/
			auto q_s = property.salt.Source(); /*kg/m³s*/
			auto Q = property.kinetics.HeatOfDissociation( q_g, T*Xc_T ); /*[W/m³]*/
      // auto q_g = Xc_source_m * property.reactionKinetics.gasGenerationRate( T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global); // property.propertyeter.ReferenceTemperature()
      // auto q_w = Xc_source_m * property.reactionKinetics.waterGenerationRate(T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      // auto q_h = Xc_source_m * property.reactionKinetics.hydrateDissociationRate(T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      // auto Q = Xc_source_h * property.reactionKinetics.heatOfDissociation(T * Xc_T, Pg * Xc_P, Sh, Sw, por, ip_global);
      
      auto Cp_g = property.gas.Cp(T * Xc_T, Pg * Xc_P, zCH4); /* ndim */
      auto Cp_w = property.water.Cp(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_g = property.gas.ThermalConductivity(T * Xc_T, Pg * Xc_P); /* ndim */
      auto kth_w = property.water.ThermalConductivity(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_h = property.hydrate.ThermalConductivity(T * Xc_T, Peff * Xc_P); /* ndim */
      auto kth_s = property.soil.ThermalConductivity(); /* ndim */
      auto kth_eff = (1. - por) * kth_s + por * (Sg * kth_g + Sw * kth_w + Sh * kth_h); /* ndim */
      //kth_eff *= Xc_diff_h;


      //Methane
      auto Concoeff_m_g = rho_g * krN * YCH4;
      auto Concoeff_m_w = rho_w * krW * XCH4;
      auto Concoeff_c_w = rho_w * krW * XC;

      auto Diffcoeff_m_g_p = DH2O_g * rho_g * Sg * -YCH4  / Pg;
      auto Diffcoeff_m_g_x = DH2O_g * rho_g * Sg * H_M_w / ( zCH4 * Pg);
      auto Diffcoeff_m_w = DCH4_w * rho_w * Sw ;
      auto Diffcoeff_c_w = DC_w * rho_w * Sw ;

      
      
      //Water
      auto Concoeff_w_g = rho_g * krN * YH2O;
      auto Concoeff_w_w = rho_w * krW * XH2O;

      auto Diffcoeff_w_w_p = DCH4_w * rho_w * Sw * YH2O / P_H_sat ;
      auto Diffcoeff_w_w_x = DCH4_w * rho_w * Sw * Pg /  P_H_sat;
      auto Diffcoeff_w_g = DH2O_g * rho_g * Sg ;

      auto gradu_Pg = gradu_Pw + coeff_grad_Sg * gradu_Sg + coeff_grad_Sh * gradu_Sh;
      
      auto Kgradu_Pg = Kgradu_Pw + coeff_grad_Sg * Kgradu_Sg + coeff_grad_Sh * Kgradu_Sh;

      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      for (size_type i = 0; i < lfsv_Pw.size(); i++)
      {
        r.accumulate(lfsv_Pw, i, ((Xc_conv_m*Concoeff_m_g * (Kgradu_Pg) * gradpsi_Pw[i] + Xc_conv_m*Concoeff_m_w * (Kgradu_Pw * gradpsi_Pw[i] )
                                                  - Xc_diff_m*Diffcoeff_m_g_p * (gradu_Pg) * gradpsi_Pw[i]) 
                                                  - Xc_diff_m*Diffcoeff_m_g_x * (gradu_XCH4 * gradpsi_Pw[i]) 
                                                  - Xc_diff_m*Diffcoeff_m_w * (gradu_XCH4 * gradpsi_Pw[i])- Xc_source_m*q_g * psi_Pw[i]) * factor);
      }
     
      for (size_type i = 0; i < lfsv_XC.size(); i++)
      {
        r.accumulate(lfsv_XC, i, (Xc_conv_m*Concoeff_c_w * (Kgradu_Pw  * gradpsi_XC[i])
                                                  - Xc_diff_m*Diffcoeff_c_w * (gradu_XC * gradpsi_XC[i])) * factor);
      }

      for (size_type i = 0; i < lfsv_Pc.size(); i++)
      {
        r.accumulate(lfsv_Pc, i, (Xc_conv_m*Concoeff_w_g * (Kgradu_Pg) * gradpsi_Pc[i] + Xc_conv_m*Concoeff_w_w * (Kgradu_Pw * gradpsi_Pc[i]) 
                                                  - Xc_diff_m*Diffcoeff_w_w_p * (gradu_Pg) * gradpsi_Pc[i]
                                                  - Xc_diff_m*Diffcoeff_w_w_x * (gradu_YH2O * gradpsi_Pc[i])
                                                  - Xc_diff_m*Diffcoeff_w_g * (gradu_YH2O * gradpsi_Pc[i]) 
                                                  - Xc_source_m*q_w * psi_Pc[i]) * factor);
      }
      for (size_type i = 0; i < lfsv_Sh.size(); i++)
      {
        r.accumulate(lfsv_Sh, i, (-Xc_source_m*q_h * psi_Sh[i]) * factor);
      }
      
      //Integrals regarding the NCP
			RF max1 = std::max(0., (Sg -1. + YCH4 + YH2O));
			for (size_type i=0; i<lfsv_YH2O.size(); i++){
				r.accumulate(lfsv_YH2O,i,( (Sg - max1) * psi_YH2O[i]  *factor));
			}

			// Integrals regarding the NCP
			RF max2 = std::max(0., (Sw -1. + XC + XCH4 + XH2O ));
			for (size_type i=0; i<lfsv_XCH4.size(); i++){
				r.accumulate(lfsv_XCH4,i,((Sw - max2) * psi_XCH4[i]  *factor));
			}
      
      
      //NCP -> water phase
			// double tmp = 0.;
			// //auto XH2O_alg = property.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,S);
			// if( ( Sw - ( 1. - XCH4 - XH2O - XC ) ) > 1.e-8 ){//active set.
			// 	tmp += 1. - XCH4 - XH2O - XC;//Active => phase is present => summation condition holds
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
			// //auto YCH4_alg = property.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,S,zCH4);
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

      for (size_type i = 0; i < lfsv_Sg.size(); i++)
      {
        r.accumulate(lfsv_Sg, i, (Pc - suctionPressure * PcSF1) * psi_Sg[i] * factor);
      }
      for (size_type i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, (Xc_conv_h*rho_g * krN * ((Kgradu_Pg) * gradpsi_T[i]) * Cp_g * T 
                                + Xc_conv_h*rho_w * krW * (Kgradu_Pw * gradpsi_T[i]) * Cp_w * T + Xc_diff_h*kth_eff * (gradu_T * gradpsi_T[i]) 
                                              - Xc_source_h*Q * psi_T[i]) * factor);
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
    const auto &lfsv_Pw_s = lfsv_s.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_s = lfsu_s.template child<Indices::PVId_Pw>();
    const auto &lfsv_Pw_n = lfsv_n.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_n = lfsu_n.template child<Indices::PVId_Pw>();

    //Capillary Pressure
    const auto &lfsv_Pc_s = lfsv_s.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_s = lfsu_s.template child<Indices::PVId_Pc>();
    const auto &lfsv_Pc_n = lfsv_n.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_n = lfsu_n.template child<Indices::PVId_Pc>();

    //Water Saturation
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
    const int dim= IG::Entity::dimension;
    const int dimension = GV::dimension;
    //std::cout << "grid dimension " << dimension << "  element dim " << dim << std::endl;
    const int order = std::max(
        std::max(lfsu_Pw_s.finiteElement().localBasis().order(),
                 lfsu_Pw_n.finiteElement().localBasis().order()),
        std::max(lfsv_Pw_s.finiteElement().localBasis().order(),
                 lfsv_Pw_n.finiteElement().localBasis().order()));

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
    auto order_s = lfsv_Pw_s.finiteElement().localBasis().order();
    auto order_n = lfsv_Pw_n.finiteElement().localBasis().order();
    auto degree = std::max(order_s, order_n);

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
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_s(lfsu_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_s(lfsv_Pc_s.size());
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
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_n(lfsu_Pc_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_n(lfsv_Pc_n.size());
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
    Dune::FieldVector<RF, dim> gradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_s(0.0);
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
    Dune::FieldVector<RF, dim> gradu_Pc_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_n(0.0);
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

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    // auto intorder = intorderadd+quadrature_factor*order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto iplocal_n = geo_in_outside.global(ip.position());

      auto ip_global_s = geo_inside.global(iplocal_s);
      auto ip_global_n = geo_outside.global(iplocal_n);

      //std::cout << "ip_global_s " << ip_global_s << "  local " << iplocal_s <<  "  ip position "<< ip.position() <<  std::endl;
      // evaluate basis functions

      auto &phi_Pw_s = cache_Pw[order].evaluateFunction(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &psi_Pw_s = cache_Pw[order].evaluateFunction(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &phi_Sg_s = cache_Sg[order].evaluateFunction(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &psi_Sg_s = cache_Sg[order].evaluateFunction(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
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
      auto &phi_XC_s = cache_XC[order].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &psi_XC_s = cache_XC[order].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      auto &phi_Pw_n = cache_Pw[order].evaluateFunction(iplocal_n, lfsu_Pw_n.finiteElement().localBasis());
      auto &psi_Pw_n = cache_Pw[order].evaluateFunction(iplocal_n, lfsv_Pw_n.finiteElement().localBasis());
      auto &phi_Sg_n = cache_Sg[order].evaluateFunction(iplocal_n, lfsu_Sg_n.finiteElement().localBasis());
      auto &psi_Sg_n = cache_Sg[order].evaluateFunction(iplocal_n, lfsv_Sg_n.finiteElement().localBasis());
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
      auto &phi_XC_n = cache_XC[order].evaluateFunction(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &psi_XC_n = cache_XC[order].evaluateFunction(iplocal_n, lfsv_XC_n.finiteElement().localBasis());

      
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
      RF T_n = 277.15;
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

      RF Sw_s = 1. - Sg_s - Sh_s;
      RF Sw_n = 1. - Sg_n - Sh_n;
      // evaluate Pw
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto suctionPressure_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */
      auto PcSF1_s = property.hydraulicProperty.PcSF1(Sh_s, BrooksCParams[1], BrooksCParams[4]);
      
      //auto Pc_s = suctionPressure_s * PcSF1_s;
      RF Pg_s = Pw_s + Pc_s;
     
      auto por_n = property.soil.SedimentPorosity(cell_outside, iplocal_n);
      auto suctionPressure_n = property.hydraulicProperty.CapillaryPressure(cell_outside, iplocal_n, Sw_n, Sh_n, por_n) ; /* ndim */
      auto PcSF1_n = property.hydraulicProperty.PcSF1(Sh_n, BrooksCParams[1], BrooksCParams[4]);
      
      //auto Pc_n = suctionPressure_n * PcSF1_n;
      RF Pg_n = Pw_n + Pc_n;

      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      // evaluate gradient of basis functions
      auto &js_Pw_s = cache_Pw[order].evaluateJacobian(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &js_v_Pw_s = cache_Pw[order].evaluateJacobian(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &js_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &js_v_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &js_Sg_s = cache_Sg[order].evaluateJacobian(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &js_v_Sg_s = cache_Sg[order].evaluateJacobian(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &js_Sh_s = cache_Sh[order].evaluateJacobian(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &js_v_Sh_s = cache_Sh[order].evaluateJacobian(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &js_XC_s = cache_XC[order].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      auto &js_Pw_n = cache_Pw[order].evaluateJacobian(iplocal_n, lfsu_Pw_n.finiteElement().localBasis());
      auto &js_v_Pw_n = cache_Pw[order].evaluateJacobian(iplocal_n, lfsv_Pw_n.finiteElement().localBasis());
      auto &js_Pc_n = cache_Pc[order].evaluateJacobian(iplocal_n, lfsu_Pc_n.finiteElement().localBasis());
      auto &js_v_Pc_n = cache_Pc[order].evaluateJacobian(iplocal_n, lfsv_Pc_n.finiteElement().localBasis());
      auto &js_Sg_n = cache_Sg[order].evaluateJacobian(iplocal_n, lfsu_Sg_n.finiteElement().localBasis());
      auto &js_v_Sg_n = cache_Sg[order].evaluateJacobian(iplocal_n, lfsv_Sg_n.finiteElement().localBasis());
      auto &js_Sh_n = cache_Sh[order].evaluateJacobian(iplocal_n, lfsu_Sh_n.finiteElement().localBasis());
      auto &js_v_Sh_n = cache_Sh[order].evaluateJacobian(iplocal_n, lfsv_Sh_n.finiteElement().localBasis());
      auto &js_T_n = cache_T[order].evaluateJacobian(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &js_v_T_n = cache_T[order].evaluateJacobian(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      auto &js_XCH4_n = cache_XCH4[order].evaluateJacobian(iplocal_n, lfsu_XCH4_n.finiteElement().localBasis());
      auto &js_v_XCH4_n = cache_XCH4[order].evaluateJacobian(iplocal_n, lfsv_XCH4_n.finiteElement().localBasis());
      auto &js_YH2O_n = cache_YH2O[order].evaluateJacobian(iplocal_n, lfsu_YH2O_n.finiteElement().localBasis());
      auto &js_v_YH2O_n = cache_YH2O[order].evaluateJacobian(iplocal_n, lfsv_YH2O_n.finiteElement().localBasis());
      auto &js_XC_n = cache_XC[order].evaluateJacobian(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &js_v_XC_n = cache_XC[order].evaluateJacobian(iplocal_n, lfsv_XC_n.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        jac.mv(js_Pw_s[i][0], gradphi_Pw_s[i]);
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        jac.mv(js_v_Pw_s[i][0], gradpsi_Pw_s[i]);
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        jac.mv(js_Pc_s[i][0], gradphi_Pc_s[i]);
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        jac.mv(js_v_Pc_s[i][0], gradpsi_Pc_s[i]);
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
        jac.mv(js_XC_s[i][0], gradphi_XCH4_s[i]);
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        jac.mv(js_v_XC_s[i][0], gradpsi_XC_s[i]);

      jac = geo_outside.jacobianInverseTransposed(iplocal_n);
      for (size_type i = 0; i < lfsu_Pw_n.size(); i++)
        jac.mv(js_Pw_n[i][0], gradphi_Pw_n[i]);
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
        jac.mv(js_v_Pw_n[i][0], gradpsi_Pw_n[i]);
      for (size_type i = 0; i < lfsu_Pc_n.size(); i++)
        jac.mv(js_Pc_n[i][0], gradphi_Pc_n[i]);
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
        jac.mv(js_v_Pc_n[i][0], gradpsi_Pc_n[i]);
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
        gradu_Pw_n.axpy(x_s(lfsu_Pw_n, i), gradphi_Pw_n[i]);

      // compute gradient of Pc
      gradu_Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        gradu_Pc_s.axpy(x_s(lfsu_Pw_s, i), gradphi_Pc_s[i]);
      gradu_Pc_n = 0.0;
      for (size_type i = 0; i < lfsu_Pc_n.size(); i++)
        gradu_Pc_n.axpy(x_s(lfsu_Pc_n, i), gradphi_Pc_n[i]);

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
      gradu_Sg_n = 0.0;
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

      auto K_s = property.soil.SedimentPermeabilityTensor(cell_inside, iplocal_s); /* ndim */
      // K_s *= 1. / Xc_K;
      // K_s *= Xc_conv_m;
      auto K_n = property.soil.SedimentPermeabilityTensor(cell_outside, iplocal_n); /* ndim */
      // K_n *= 1. / Xc_K;
      // K_n *= Xc_conv_m;
      
      
      // compute K * gradient of Pw
      K_s.mv(gradu_Pw_s, Kgradu_Pw_s);
      K_n.mv(gradu_Pw_n, Kgradu_Pw_n);

      // compute K * gradient of Pc
      K_s.mv(gradu_Pc_s, Kgradu_Pc_s);
      K_n.mv(gradu_Pc_n, Kgradu_Pc_n);

      // compute K * gradient of Sg
      K_s.mv(gradu_Sg_s, Kgradu_Sg_s);
      K_n.mv(gradu_Sg_n, Kgradu_Sg_n);
      
      // compute K * gradient of Sh
      K_s.mv(gradu_Sh_s, Kgradu_Sh_s);
      K_n.mv(gradu_Sh_n, Kgradu_Sh_n);

      auto n_F = ig.centerUnitOuterNormal();
      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K_n.mv(n_F, Kn_F_n);

      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s,Sh_s,0.0,0.0);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams[0], BrooksCParams[1]);
      auto dSwe_dSg_s = -property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sg_s = dPc_dSwe_s * dSwe_dSg_s ;

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s + dPc_dSwe_s * dSwe_dSh_s ;

      auto Swe_n = property.hydraulicProperty.EffectiveSw(Sw_n,Sh_n,0.0,0.0);
      auto dPc_dSwe_n =  property.hydraulicProperty.dPc_dSwe(Swe_n, BrooksCParams[0], BrooksCParams[1]);
      auto dSwe_dSg_n = -property.hydraulicProperty.dSwe_dSw(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sw_n = dPc_dSwe_n * dSwe_dSg_n ;

      auto dPcSF1_dSh_n =  property.hydraulicProperty.dPcSF1_dSh( Sh_n, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_n = property.hydraulicProperty.dSwe_dSh(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sh_n = dPcSF1_dSh_n + dPc_dSwe_n * dSwe_dSh_n ;


      auto Kgradu_Pg_s = Kgradu_Pw_s + coeff_grad_Sg_s * Kgradu_Sg_s + coeff_grad_Sh_s * Kgradu_Sh_s;
      auto gradu_Pg_s = gradu_Pw_s + coeff_grad_Sg_s * gradu_Sg_s + coeff_grad_Sh_s * gradu_Sh_s;

      auto Kgradu_Pg_n = Kgradu_Pw_n + coeff_grad_Sw_n * Kgradu_Sg_n + coeff_grad_Sh_n * Kgradu_Sh_n;
      auto gradu_Pg_n = gradu_Pw_n + coeff_grad_Sw_n * gradu_Sg_n + coeff_grad_Sh_n * gradu_Sh_n;
      


      for (int i = 0; i < dim; i++)
      {
        v_g[i] = (omega_s * (gradu_Pg_s[i]) + omega_n * (gradu_Pg_n[i]));
        v_w[i] = (omega_s * (gradu_Pw_s[i] ) + omega_n * (gradu_Pw_n[i] ));
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

      
      double S_s = XC_s * (property.salt.MolarMass()/property.water.MolarMass());
      //auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s) ); /* ndim */
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.gas.DynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) ); /* ndim */
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s); 
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s * Xc_T, Pg_s * Xc_P); /* ndim */
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s * Xc_T, Pw_s * Xc_P); /* ndim */
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P); /* ndim */

      auto H_M_w_s = property.gas.SolubilityCoefficient(T_s * Xc_T, S_s ); /* ndim */
      auto P_H_sat_s = property.water.SaturatedVaporPressure(T_s * Xc_T, S_s ); /* ndim */
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pw_s * Xc_P); 
      //  end of terms regarding components
      auto VLequil_s = property.mixture.EquilibriumMoleFractions( T_s * Xc_T, Pg_s * Xc_P, XC_s, zCH4_s);
      auto YCH4_s = property.mixture.YCH4(XCH4_s, T_s * Xc_T, Pg_s * Xc_P, XC_s, zCH4_s);
      auto XH2O_s = property.mixture.XH2O(YH2O_s, T_s * Xc_T, Pg_s * Xc_P, XC_s);
      
      auto rho_g_s = property.gas.Density(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ; /* ndim */
      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
      
      auto Cp_g_s = property.gas.Cp(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ; /* ndim */
      auto Cp_w_s = property.water.Cp(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
      auto kth_g_s = property.gas.ThermalConductivity(T_s * Xc_T, Pg_s * Xc_P) ; /* ndim */
      auto kth_w_s = property.water.ThermalConductivity(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s * Xc_T, Peff_s * Xc_P); /* ndim */
      auto kth_s_s = property.soil.ThermalConductivity() ; /* ndim */
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * (Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s); /* ndim */
      //kth_eff_s *= Xc_diff_h;


      double S_n = XC_n * (property.salt.MolarMass()/property.water.MolarMass());
      //auto por_n = property.soil.SedimentPorosity(cell_outside, iplocal_n);
      auto krW_n = property.hydraulicProperty.krw(cell_outside, iplocal_n, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_n) ); /* ndim */
      auto krN_n = property.hydraulicProperty.krg(cell_outside, iplocal_n, Sw_n, Sh_n) / (property.gas.DynamicViscosity(T_n * Xc_T, Pg_n * Xc_P) ); /* ndim */
      
      //  adding terms regarding components
      
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n * Xc_T, Pg_n * Xc_P); /* ndim */
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n * Xc_T, Pw_n * Xc_P); /* ndim */
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P); /* ndim */

      auto H_M_w_n = property.gas.SolubilityCoefficient(T_n * Xc_T, S_n ); /* ndim */
      auto P_H_sat_n = property.water.SaturatedVaporPressure(T_n * Xc_T, S_n ); /* ndim */
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n * Xc_T, Pw_n * Xc_P); 
      //  end of terms regarding components
      auto YCH4_n = property.mixture.YCH4(XCH4_n, T_n * Xc_T, Pg_n * Xc_P, XC_n, zCH4_n);
      auto XH2O_n = property.mixture.XH2O(YH2O_n, T_n * Xc_T, Pg_n * Xc_P, XC_n);
  
      auto rho_g_n = property.gas.Density(T_n * Xc_T, Pg_n * Xc_P, zCH4_n) ; /* ndim */
      auto rho_w_n = property.water.Density(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      
      auto Cp_g_n = property.gas.Cp(T_n * Xc_T, Pg_n * Xc_P, zCH4_n) ; /* ndim */
      auto Cp_w_n = property.water.Cp(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      auto kth_g_n = property.gas.ThermalConductivity(T_n * Xc_T, Pg_n * Xc_P) ; /* ndim */
      auto kth_w_n = property.water.ThermalConductivity(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n * Xc_T, Peff_n * Xc_P) ; /* ndim */
      auto kth_s_n = property.soil.ThermalConductivity() ; /* ndim */
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * (Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);
      //kth_eff_n *= Xc_diff_h;

      auto krN = omega_s * krN_s + omega_n * krN_n;
      auto krW = omega_s * krW_s + omega_n * krW_n;
      auto h_g = omega_s * Cp_g_s * T_s + omega_n * Cp_g_n * T_n;
      auto h_w = omega_s * Cp_w_s * T_s + omega_n * Cp_w_n * T_n;
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
       //std::cout << "H_M_w_n = " << H_M_w_n << " P_H_sat_n = " << P_H_sat_n << " Cp_g_n = " << Cp_g_n << " rho_g_n = " << rho_g_n <<std::endl;
    //     " source_m = " << Xc_source_m << " source_h = " << Xc_source_h << std::endl;
        // exit(0);
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
      
      double term_conv_m_g = -krN * (omega_s * rho_g_s * YCH4_s * ((Kgradu_Pg_s) * n_F_local) 
                                    + omega_n * rho_g_n * YCH4_n * ((Kgradu_Pg_n) * n_F_local));
      double term_conv_m_w = -krW * (omega_s * rho_w_s * XCH4_s * (Kgradu_Pw_s * n_F_local) 
                                    + omega_n * rho_w_n * XCH4_n * (Kgradu_Pw_n * n_F_local));

      double term_diffusion_m = Xc_conv_m*term_conv_m_g + Xc_conv_m*term_conv_m_w + Xc_diff_m*DH2O_g * (omega_s * (Diffcoeff_m_g_p_s * ((gradu_Pg_s) * n_F_local) 
                                              + Diffcoeff_m_g_x_s * (gradu_XCH4_s * n_F_local))
                                              + omega_n * (Diffcoeff_m_g_p_n * ((gradu_Pg_n) * n_F_local) 
                                              + Diffcoeff_m_g_x_n * (gradu_XCH4_n * n_F_local)))
                                              + Xc_diff_m * DCH4_w * (omega_s * Diffcoeff_m_w_s * (gradu_XCH4_s * n_F_local)
                                              + omega_n * Diffcoeff_m_w_n * (gradu_XCH4_n * n_F_local)) ;

      double term_nipg_g = theta_g * ((Pg_s) - (Pg_n));
      double term_nipg_w = theta_w * (Pw_s - Pw_n);
      double term_nipg_m_x = theta_x * (XCH4_s - XCH4_n);

      double term_penalty_w = penalty_factor_w * (Pw_s - Pw_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i, term_diffusion_m * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i, term_diffusion_m * -psi_Pw_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i, (Xc_conv_m*term_nipg_g * krN * omega_s * rho_g_s * YCH4_s * (Kn_F_s * gradpsi_Pw_s[i]) 
                                      + Xc_conv_m*term_nipg_w * krW * omega_s * rho_w_s * XCH4_s  * (Kn_F_s * gradpsi_Pw_s[i])
                                      + Xc_diff_m*omega_s * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_s 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_s
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_s) * (n_F_local* gradpsi_Pw_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i, (Xc_conv_m*term_nipg_g * krN * omega_n * rho_g_n * YCH4_n * (Kn_F_n * gradpsi_Pw_n[i]) 
                                      + Xc_conv_m*term_nipg_w * krW * omega_n * rho_w_n * XCH4_n  * (Kn_F_n * gradpsi_Pw_n[i])
                                      + Xc_diff_m*omega_n * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_n 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_n
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_n) * (n_F_local* gradpsi_Pw_n[i])) * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pw_s, i, term_penalty_w * psi_Pw_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pw_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pw_n, i, term_penalty_w * -psi_Pw_n[i] * factor);
      }
      
      // Salt
      
      
      auto DC_w = omega_s * DC_w_s + omega_n * DC_w_n; 
      auto Diffcoeff_c_w_s = rho_w_s * Sw_s ;
      auto Diffcoeff_c_w_n = rho_w_n * Sw_n ;
      
      
      double term_conv_c_w = -krW * (omega_s * rho_w_s * XC_s * Kgradu_Pw_s * n_F_local  
                                    + omega_n * rho_w_n * XC_n * Kgradu_Pw_n * n_F_local);

      double term_diffusion_c = Xc_conv_m*term_conv_c_w + Xc_diff_m* DC_w * (omega_s * Diffcoeff_c_w_s * (gradu_XC_s * n_F_local)
                                              + omega_n * Diffcoeff_c_w_n * (gradu_XC_n * n_F_local)) ;

      
      double term_nipg_c_x = theta_x * (XC_s - XC_n);

      double term_penalty_c = penalty_factor_x * (XC_s - XC_n);
      // diffusion term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, term_diffusion_c * psi_XC_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, term_diffusion_c * -psi_XC_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, ( Xc_conv_m*term_nipg_w * krW * omega_s * rho_w_s * XC_s  * (Kn_F_s * gradpsi_XC_s[i])
                                      + Xc_diff_m*omega_s * ( term_nipg_c_x * DC_w * Diffcoeff_c_w_s) * (n_F_local* gradpsi_XC_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, (Xc_conv_m*term_nipg_w * krW * omega_n * rho_w_n * XC_n  * (Kn_F_n * gradpsi_XC_n[i])
                                      + Xc_diff_m*omega_n * ( term_nipg_c_x * DC_w * Diffcoeff_c_w_n) * (n_F_local* gradpsi_XC_n[i])) * factor);
      }
      // standard IP term integral
      for (size_type i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, term_penalty_c * psi_XC_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, term_penalty_c * -psi_XC_n[i] * factor);
      }
     
      
      //Water
      //auto DH2O_g = DCH4_g
      auto Diffcoeff_w_w_p_s = rho_w_s * Sw_s * YH2O_s / P_H_sat_s;
      auto Diffcoeff_w_w_y_s = rho_w_s * Sw_s * Pg_s / P_H_sat_s;
      auto Diffcoeff_w_w_p_n = rho_w_n * Sw_n * YH2O_n  / P_H_sat_n;
      auto Diffcoeff_w_w_y_n = rho_w_n * Sw_n * Pg_n / P_H_sat_n;
      
      //auto DCH4_w = DH2O_w
      auto Diffcoeff_w_g_s = rho_g_s * Sg_s ;
      auto Diffcoeff_w_g_n = rho_g_n * Sg_n ;

      double term_conv_w_g = -krN * ( omega_s * rho_g_s * YH2O_s * (Kgradu_Pg_s) * n_F_local 
                                    + omega_n * rho_g_n * YH2O_n * (Kgradu_Pg_n) * n_F_local);

      double term_conv_w_w = -krW * ( omega_s * rho_w_s * XH2O_s * Kgradu_Pw_s * n_F_local 
                                    + omega_n * rho_w_n * XH2O_n * Kgradu_Pw_n * n_F_local );

      double term_diffusion_w = Xc_conv_m*term_conv_w_g + Xc_conv_m*term_conv_w_w + Xc_diff_m*DH2O_g * (omega_s * Diffcoeff_w_g_s * (gradu_YH2O_s * n_F_local) 
                              + omega_n * Diffcoeff_w_g_n * (gradu_YH2O_n * n_F_local)) 
                              + Xc_diff_m*DCH4_w * (omega_s * (Diffcoeff_w_w_p_s * ((gradu_Pg_s ) * n_F_local)
                              + Diffcoeff_w_w_y_s * (gradu_YH2O_s * n_F_local))
                              + omega_n * (Diffcoeff_w_w_p_n * ((gradu_Pg_n ) * n_F_local)
                              + Diffcoeff_w_w_y_n * (gradu_YH2O_n * n_F_local)));


      
      double term_nipg_w_y = theta_y * (YH2O_s - YH2O_n);
      double term_penalty_p = penalty_factor_g * ((Pg_s) - (Pg_n));
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
        r_s.accumulate(lfsv_Pc_s, i, (Xc_conv_m*term_nipg_g * krN * omega_s * rho_g_s * YH2O_s * (Kn_F_s * gradpsi_Pc_s[i]) 
                                      + Xc_conv_m*term_nipg_w * krW * omega_s * rho_w_s * XH2O_s  * (Kn_F_s * gradpsi_Pc_s[i])
                                      + Xc_diff_m*omega_s * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_s 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_s
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_s) * (n_F_local* gradpsi_Pc_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pc_n, i, (Xc_conv_m*term_nipg_g * krN * omega_n * rho_g_n * YH2O_n * (Kn_F_n * gradpsi_Pc_n[i]) 
                                      + Xc_conv_m*term_nipg_w * krW * omega_n * rho_w_n * XH2O_n  * (Kn_F_n * gradpsi_Pc_n[i])
                                      + Xc_diff_m*omega_n * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_n 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_n
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_n) * (n_F_local* gradpsi_Pc_n[i])) * factor);
      }
      //standard IP term integral
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
      {
        r_s.accumulate(lfsv_Pc_s, i, term_penalty_p * psi_Pc_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Pc_n.size(); i++)
      {
        r_n.accumulate(lfsv_Pc_n, i, term_penalty_p * -psi_Pc_n[i] * factor);
      }

      
      double term_penalty_sg = penalty_factor_s * (Sg_s - Sg_n);
      // standard IP term integral
      for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
      {
        r_s.accumulate(lfsv_Sg_s, i, term_penalty_sg * psi_Sg_s[i] * factor);
      }
      for (size_type i = 0; i < lfsv_Sg_n.size(); i++)
      {
        r_n.accumulate(lfsv_Sg_n, i, term_penalty_sg * -psi_Sg_n[i] * factor);
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

      double term_diffusion_T_1 = - h_g * krN * (omega_s * rho_g_s * ((Kgradu_Pg_s) * n_F_local) + omega_n * rho_g_n * ((Kgradu_Pg_n) * n_F_local));
      double term_diffusion_T_2 = -h_w * krW * (omega_s * rho_w_s * (Kgradu_Pw_s ) * n_F_local + omega_n * rho_w_n * (Kgradu_Pw_n) * n_F_local);
      double term_diffusion_T_3 = -kth_eff * (omega_s * (gradu_T_s * n_F_local) + omega_n * (gradu_T_n * n_F_local));
      double term_diffusion_T = Xc_conv_h*(term_diffusion_T_1 + term_diffusion_T_2) + Xc_diff_h*term_diffusion_T_3;
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
        r_s.accumulate(lfsv_T_s, i, Xc_diff_h*term_nipg_T * (kth_eff * omega_s * (n_F * gradpsi_T_s[i])) * factor);
      }
      for (size_type i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, Xc_diff_h*term_nipg_T * (kth_eff * omega_n * (n_F * gradpsi_T_n[i])) * factor);
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
    // std::cout << "conv_m = " << Xc_conv_m << " diff_m = " << Xc_diff_m << " conv_h = " << Xc_conv_h << " diff_h = " << Xc_diff_h <<
    //     " source_m = " << Xc_source_m << " source_h = " << Xc_source_h << std::endl;
    //     exit(0);
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
    const auto &lfsv_Pw_s = lfsv.template child<Indices::PVId_Pw>();
    const auto &lfsu_Pw_s = lfsu.template child<Indices::PVId_Pw>();

    //Capillary Pressure
    const auto &lfsv_Pc_s = lfsv.template child<Indices::PVId_Pc>();
    const auto &lfsu_Pc_s = lfsu.template child<Indices::PVId_Pc>();

    //Water Saturation
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
    const int order = std::max(lfsu_Pw_s.finiteElement().localBasis().order(),
                               lfsv_Pw_s.finiteElement().localBasis().order());

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();
    
    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();

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
    std::vector<Dune::FieldVector<RF, dim>> gradphi_Pc_s(lfsu_Pc_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_Pc_s(lfsv_Pc_s.size());
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
    Dune::FieldVector<RF, dim> gradu_Pc_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pc_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // auto intorder = intorderadd+quadrature_factor*order;
    // loop over quadrature points
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto ip_global_s = geo_inside.global(iplocal_s);

      BC bc( gv,property ) ;
      
      // evaluate boundary condition types for {Pw,Sg} or {Fw,Fg} 
			auto bctype = bc.type(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t ) ;

			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg} 
			auto bcvalue = bc.value(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t ) ;
      //auto bctype = property.problemSpecs.ProblemBCTypes();
      //auto bcvalue = property.problemSpecs.ProblemBCValues(ip_global_s, *(time));
      
      // evaluate basis functions
      auto &psi_Pw_s = cache_Pw[order].evaluateFunction(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &psi_Sg_s = cache_Sg[order].evaluateFunction(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &psi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &psi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &psi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &psi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &psi_XC_s = cache_XC[order].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
      {

        double flux_w = bcvalue[Indices::PVId_Pw];
        for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        {
          r.accumulate(lfsv_Pw_s, i, flux_w * psi_Pw_s[i] * factor);
        }
      }

      if (bctype[Indices::PVId_Pc] == Indices::BCId_neumann)
      {

        double flux_g = bcvalue[Indices::PVId_Pc];
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, flux_g * psi_Pc_s[i] * factor);
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
          bctype[Indices::PVId_Sg] == Indices::BCId_neumann and
          bctype[Indices::PVId_Pw] == Indices::BCId_neumann and
          bctype[Indices::PVId_XCH4] == Indices::BCId_neumann and
          bctype[Indices::PVId_C] == Indices::BCId_neumann and
          bctype[Indices::PVId_YH2O] == Indices::BCId_neumann)
      {
        continue;
      }

      auto &phi_Pw_s = cache_Pw[order].evaluateFunction(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &phi_Sg_s = cache_Sg[order].evaluateFunction(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &phi_Sh_s = cache_Sh[order].evaluateFunction(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &phi_Pc_s = cache_Pc[order].evaluateFunction(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &phi_XCH4_s = cache_XCH4[order].evaluateFunction(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &phi_YH2O_s = cache_YH2O[order].evaluateFunction(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &phi_XC_s = cache_XC[order].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      

      // evaluate Pw
      RF Pw_s = 0.0;
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        Pw_s += x(lfsu_Pw_s, i) * phi_Pw_s[i];
      RF Pw_n = 0.0;
      if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
      {
        Pw_n = Pw_s;
      }
      else if (bctype[Indices::PVId_Pw] == Indices::BCId_dirichlet)
      {
        Pw_n = bcvalue[Indices::PVId_Pw] / Xc_P;
      }

      // evaluate Sh
      RF Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        Sh_s += x(lfsu_Sh_s, i) * phi_Sh_s[i];
      RF Sh_n = 0.;
      if (bctype[Indices::PVId_Sh] == Indices::BCId_neumann)
      {
        Sh_n = Sh_s;
      }
      else if (bctype[Indices::PVId_Sh] == Indices::BCId_dirichlet)
      {
        Sh_n = bcvalue[Indices::PVId_Sh] ;
      }

      // evaluate Sg
      RF Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        Sg_s += x(lfsu_Sg_s, i) * phi_Sg_s[i];
      RF Sg_n = 0.0;
      if (bctype[Indices::PVId_Sg] == Indices::BCId_neumann)
      {
        Sg_n = Sg_s;
      }
      else if (bctype[Indices::PVId_Sg] == Indices::BCId_dirichlet)
      {
        Sg_n=bcvalue[Indices::PVId_Sg];
        //Sg_n = Sg_s;
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
        //auto PcSF1_n = propertyclass.hydraulicProperty.PcSF1(Sh_n);
        //Pc_n=propertyclass.hydraulicProperty.suctionPressure(Sw_n,Sh_n)*PcSF1_n;
        Pc_n = bcvalue[Indices::PVId_Pc];// * property.hydraulicProperty.PcSF1(Sh_n) / Xc_P;
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
        XCH4_n = bcvalue[Indices::PVId_XCH4] ;
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
        YH2O_n = bcvalue[Indices::PVId_YH2O] ;
      }

      // evaluate XC
      RF XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        XC_s += x(lfsu_XC_s, i) * phi_XC_s[i];
      RF XC_n = 0.0;
      if (bctype[Indices::PVId_C] == Indices::BCId_neumann)
      {
        XC_n = XC_s;
      }
      else if (bctype[Indices::PVId_C] == Indices::BCId_dirichlet)
      {
        XC_n = bcvalue[Indices::PVId_C] ;
      }

      RF Sw_s = 1. - Sg_s - Sh_s;
      RF Sw_n = 1. - Sg_n - Sh_n;
      // evaluate Pw
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto suctionPressure_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */
      auto PcSF1_s = property.hydraulicProperty.PcSF1(Sh_s, BrooksCParams[1], BrooksCParams[4]);
      
      //auto Pc_s = suctionPressure_s * PcSF1_s;
      RF Pg_s = Pw_s + Pc_s;

      auto por_n = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto suctionPressure_n = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_n, Sh_n, por_n) ; /* ndim */
      auto PcSF1_n = property.hydraulicProperty.PcSF1(Sh_n, BrooksCParams[1], BrooksCParams[4]);
      
      //auto Pc_n = suctionPressure_n * PcSF1_n;
      RF Pg_n = Pw_n + Pc_n;
      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      // evaluate gradient of basis functions
      auto &js_Pw_s = cache_Pw[order].evaluateJacobian(iplocal_s, lfsu_Pw_s.finiteElement().localBasis());
      auto &js_v_Pw_s = cache_Pw[order].evaluateJacobian(iplocal_s, lfsv_Pw_s.finiteElement().localBasis());
      auto &js_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsu_Pc_s.finiteElement().localBasis());
      auto &js_v_Pc_s = cache_Pc[order].evaluateJacobian(iplocal_s, lfsv_Pc_s.finiteElement().localBasis());
      auto &js_Sg_s = cache_Sg[order].evaluateJacobian(iplocal_s, lfsu_Sg_s.finiteElement().localBasis());
      auto &js_v_Sg_s = cache_Sg[order].evaluateJacobian(iplocal_s, lfsv_Sg_s.finiteElement().localBasis());
      auto &js_Sh_s = cache_Sh[order].evaluateJacobian(iplocal_s, lfsu_Sh_s.finiteElement().localBasis());
      auto &js_v_Sh_s = cache_Sh[order].evaluateJacobian(iplocal_s, lfsv_Sh_s.finiteElement().localBasis());
      auto &js_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsu_XCH4_s.finiteElement().localBasis());
      auto &js_v_XCH4_s = cache_XCH4[order].evaluateJacobian(iplocal_s, lfsv_XCH4_s.finiteElement().localBasis());
      auto &js_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsu_YH2O_s.finiteElement().localBasis());
      auto &js_v_YH2O_s = cache_YH2O[order].evaluateJacobian(iplocal_s, lfsv_YH2O_s.finiteElement().localBasis());
      auto &js_XC_s = cache_XC[order].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);
      for (size_type i = 0; i < lfsu_Pw_s.size(); i++)
        jac.mv(js_Pw_s[i][0], gradphi_Pw_s[i]);
      for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        jac.mv(js_v_Pw_s[i][0], gradpsi_Pw_s[i]);
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        jac.mv(js_Pc_s[i][0], gradphi_Pc_s[i]);
      for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        jac.mv(js_v_Pc_s[i][0], gradpsi_Pc_s[i]);
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

      // compute gradient of Pc
      gradu_Pc_s = 0.0;
      for (size_type i = 0; i < lfsu_Pc_s.size(); i++)
        gradu_Pc_s.axpy(x(lfsu_Pc_s, i), gradphi_Pc_s[i]);

      // compute gradient of Sg
      gradu_Sg_s = 0.0;
      for (size_type i = 0; i < lfsu_Sg_s.size(); i++)
        gradu_Sg_s.axpy(x(lfsu_Sg_s, i), gradphi_Sg_s[i]);
      

      // compute gradient of Sh
      gradu_Sh_s = 0.0;
      for (size_type i = 0; i < lfsu_Sh_s.size(); i++)
        gradu_Sh_s.axpy(x(lfsu_Sh_s, i), gradphi_Sh_s[i]);
     
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

      // compute gradient of XC
      gradu_XC_s = 0.0;
      for (size_type i = 0; i < lfsu_XC_s.size(); i++)
        gradu_XC_s.axpy(x(lfsu_XC_s, i), gradphi_XC_s[i]);


      auto K_s = property.soil.SedimentPermeabilityTensor(cell_inside,  iplocal_s);
      // K_s *= 1. / Xc_K;
      // K_s *= Xc_conv_m;
      auto n_F = ig.centerUnitOuterNormal(); // n_F = n_F_local


      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F, Kn_F_s);

      // compute K * gradient of Pw
      K_s.mv(gradu_Pw_s, Kgradu_Pw_s);

      // compute K * gradient of Pc
      K_s.mv(gradu_Pc_s, Kgradu_Pc_s);

      // compute K * gradient of Sg
      K_s.mv(gradu_Sg_s, Kgradu_Sg_s);

      // compute K * gradient of Sh
      K_s.mv(gradu_Sh_s, Kgradu_Sh_s);


      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s,Sh_s,0.0,0.0);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams[0], BrooksCParams[1]);
      auto dSwe_dSg_s = -property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sg_s = dPc_dSwe_s * dSwe_dSg_s ;

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s + dPc_dSwe_s * dSwe_dSh_s ;

      

      auto Kgradu_Pg_s = Kgradu_Pw_s + coeff_grad_Sg_s * Kgradu_Sg_s + coeff_grad_Sh_s * Kgradu_Sh_s;
      auto gradu_Pg_s = gradu_Pw_s + coeff_grad_Sg_s * gradu_Sg_s + coeff_grad_Sh_s * gradu_Sh_s;

      for (int i = 0; i < dim; i++)
      {
        v_w[i] = omega_s * gradu_Pw_s[i];
        v_g[i] = omega_s * (gradu_Pg_s[i]);
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

      //auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      double S_s = XC_s * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s));
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.gas.DynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) );
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s);
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s * Xc_T, Pg_s * Xc_P);
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s * Xc_T, Pw_s * Xc_P);
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P);
      
      
      auto H_M_w_s = property.gas.SolubilityCoefficient(T_s * Xc_T, S_s );
      auto P_H_sat_s = property.water.SaturatedVaporPressure(T_s * Xc_T, S_s );
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pw_s * Xc_P);
      //  end of terms regarding components
      auto YCH4_s = property.mixture.YCH4(XCH4_s, T_s * Xc_T, Pg_s * Xc_P, XC_s, zCH4_s);
      auto XH2O_s = property.mixture.XH2O(YH2O_s, T_s * Xc_T, Pg_s * Xc_P, XC_s);
      
      auto rho_g_s = property.gas.Density(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ;
      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s);
      
      auto Cp_g_s = property.gas.Cp(T_s * Xc_T, Pg_s * Xc_P, zCH4_s);
      auto Cp_w_s = property.water.Cp(T_s * Xc_T, Pw_s * Xc_P, S_s);
      auto kth_g_s = property.gas.ThermalConductivity(T_s * Xc_T, Pg_s * Xc_P) ;
      auto kth_w_s = property.water.ThermalConductivity(T_s * Xc_T, Pw_s * Xc_P, S_s);
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s * Xc_T, Peff_s * Xc_P);
      auto kth_s_s = property.soil.ThermalConductivity() ;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * (Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s);

      //auto por_n = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto krW_n = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_s));
      auto krN_n = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.gas.DynamicViscosity(T_n * Xc_T, Pg_n * Xc_P));
      
      //  adding terms regarding components
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n * Xc_T, Pg_n * Xc_P);
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n * Xc_T, Pw_n * Xc_P);
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P);
      //  end of terms regarding components
      auto Cp_g_n = property.gas.Cp(T_n * Xc_T, Pg_n * Xc_P, zCH4_s) ;
      auto Cp_w_n = property.water.Cp(T_n * Xc_T, Pw_n * Xc_P, S_s);
      auto kth_g_n = property.gas.ThermalConductivity(T_n * Xc_T, Pg_n * Xc_P) ;
      auto kth_w_n = property.water.ThermalConductivity(T_n * Xc_T, Pw_n * Xc_P, S_s);
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n * Xc_T, Peff_n * Xc_P);
      auto kth_s_n = property.soil.ThermalConductivity() ;
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * (Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);

      auto krN = omega_s * krN_s ;//+ omega_n * krN_n;
      auto krW = omega_s * krW_s ;//+ omega_n * krW_n;
      auto h_g = omega_s * Cp_g_s * T_s ;//+ omega_n * Cp_g_n * T_n;
      auto h_w = omega_s * Cp_w_s * T_s ;//+ omega_n * Cp_w_n * T_n;
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
      //kth_eff *= Xc_diff_h;

      auto DH2O_g = omega_s * DH2O_g_s ;//+ omega_n * DH2O_g_n; // = DCH4_g
      auto DCH4_w = omega_s * DCH4_w_s ;//+ omega_n * DCH4_w_n; // = DH2O_w
      auto DC_w = omega_s * DC_w_s ;//+ omega_n * DC_w_n; // = DH2O_w

      if (bctype[Indices::PVId_Pw] == Indices::BCId_dirichlet)
      {
        //Methane
        auto Diffcoeff_m_g_p_s = rho_g_s * Sg_s * -YCH4_s / Pg_s;
        auto Diffcoeff_m_g_x_s = rho_g_s * Sg_s * H_M_w_s / ( zCH4_s * Pg_s);
      
        auto Diffcoeff_m_w_s = rho_w_s * Sw_s ;

        double term_conv_m_g = -krN * omega_s * rho_g_s * YCH4_s * (Kgradu_Pg_s) * n_F_local; 
        double term_conv_m_w = -krW * omega_s * rho_w_s * XCH4_s * (Kgradu_Pw_s * n_F_local) ; 

        double term_diffusion_m = Xc_conv_m*term_conv_m_g + Xc_conv_m*term_conv_m_w + Xc_diff_m*DH2O_g * (omega_s * (Diffcoeff_m_g_p_s * ((gradu_Pg_s) * n_F_local) 
                                              + Diffcoeff_m_g_x_s * (gradu_XCH4_s * n_F_local)))
                                              + Xc_diff_m*DCH4_w * omega_s * Diffcoeff_m_w_s * (gradu_XCH4_s * n_F_local);
                                             
        double term_nipg_g = theta_g * (Pg_s - Pg_n);
        double term_nipg_w = theta_w * (Pw_s - Pw_n);
        double term_nipg_m_x = theta_x * (XCH4_s - XCH4_n);

        double term_penalty_w = penalty_factor_w * (Pw_s - Pw_n);
        
        // diffusion term
        for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        {
          r.accumulate(lfsv_Pw_s, i, term_diffusion_m * psi_Pw_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        {
          r.accumulate(lfsv_Pw_s, i, (Xc_conv_m*term_nipg_g * krN * omega_s * rho_g_s * YCH4_s * (Kn_F_s * gradpsi_Pw_s[i]) 
                                      + Xc_conv_m*term_nipg_w * krW * omega_s * rho_w_s * XCH4_s  * (Kn_F_s * gradpsi_Pw_s[i])
                                      + Xc_diff_m*omega_s * (term_nipg_g * DH2O_g *  Diffcoeff_m_g_p_s 
                                      + term_nipg_m_x * DH2O_g * Diffcoeff_m_g_x_s
                                      + term_nipg_m_x * DCH4_w * Diffcoeff_m_w_s) * (n_F_local* gradpsi_Pw_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Pw_s.size(); i++)
        {
          r.accumulate(lfsv_Pw_s, i, term_penalty_w * psi_Pw_s[i] * factor);
        }
      }


      if (bctype[Indices::PVId_C] == Indices::BCId_dirichlet)
      {
        // Salt
        auto Diffcoeff_c_w_s = rho_w_s * Sw_s ;      
        
        double term_conv_c_w = -krW * omega_s * rho_w_s * XC_s * Kgradu_Pw_s  * n_F_local;

        double term_diffusion_c = Xc_conv_m*term_conv_c_w + Xc_diff_m*DC_w * omega_s * Diffcoeff_c_w_s * (gradu_XC_s * n_F_local) ;

        
        double term_nipg_c_x = theta_x * (XC_s - XC_n);
        double term_nipg_w = theta_w * (Pw_s - Pw_n);
        double term_penalty_c = penalty_factor_x * (XC_s - XC_n);
        // diffusion term
        for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        {
          r.accumulate(lfsv_XC_s, i, term_diffusion_c * psi_XC_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        {
          r.accumulate(lfsv_XC_s, i, ( Xc_conv_m*term_nipg_w * krW * omega_s * rho_w_s * XC_s  * (Kn_F_s * gradpsi_XC_s[i])
                                        +Xc_diff_m* omega_s * ( term_nipg_c_x * DC_w * Diffcoeff_c_w_s) * (n_F_local* gradpsi_XC_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_XC_s.size(); i++)
        {
          r.accumulate(lfsv_XC_s, i, term_penalty_c * psi_XC_s[i] * factor);
        }
       
      }
     

      if (bctype[Indices::PVId_Pc] == Indices::BCId_dirichlet)
      {
        
        auto Diffcoeff_w_w_p_s = rho_w_s * Sw_s * YH2O_s / P_H_sat_s;
        auto Diffcoeff_w_w_y_s = rho_w_s * Sw_s * Pg_s / P_H_sat_s;
      
        auto Diffcoeff_w_g_s = rho_g_s * Sg_s ;
     
        double term_conv_w_g = -krN *  omega_s * rho_g_s * YH2O_s * ((Kgradu_Pg_s) * n_F_local);

        double term_conv_w_w = -krW * omega_s * rho_w_s * XH2O_s * Kgradu_Pw_s * n_F_local;

        double term_diffusion_w = Xc_conv_m*term_conv_w_g + Xc_conv_m*term_conv_w_w + Xc_diff_m*(DH2O_g * omega_s * Diffcoeff_w_g_s * (gradu_YH2O_s * n_F_local) 
                                + DCH4_w * (omega_s * (Diffcoeff_w_w_p_s * ((gradu_Pg_s) * n_F_local)
                                + Diffcoeff_w_w_y_s * (gradu_YH2O_s * n_F_local))));


        double term_nipg_g = theta_g * (Pg_s - Pg_n);
        double term_nipg_w = theta_w * (Pw_s - Pw_n);
        
        double term_nipg_w_y = theta_y * (YH2O_s - YH2O_n);
        double term_penalty_g = penalty_factor_g * (Pg_s - Pg_n);
        
        // diffusion term
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, term_diffusion_w * psi_Pc_s[i] * factor);
        }
        // (non-)symmetric IP term
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, (Xc_conv_m*term_nipg_g * krN * omega_s * rho_g_s * YH2O_s * (Kn_F_s * gradpsi_Pc_s[i]) 
                                      +Xc_conv_m* term_nipg_w * krW * omega_s * rho_w_s * XH2O_s  * (Kn_F_s * gradpsi_Pc_s[i])
                                      + Xc_diff_m*omega_s * (term_nipg_g * DCH4_w *  Diffcoeff_w_w_p_s 
                                      + term_nipg_w_y * DCH4_w * Diffcoeff_w_w_y_s
                                      + term_nipg_w_y * DH2O_g * Diffcoeff_w_g_s)  * (n_F_local* gradpsi_Pc_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Pc_s.size(); i++)
        {
          r.accumulate(lfsv_Pc_s, i, term_penalty_g * psi_Pc_s[i] * factor);
        }
         
      }

      if (bctype[Indices::PVId_Sg] == Indices::BCId_dirichlet)
      {

        double term_penalty_sg = penalty_factor_s * (Sg_s - Sg_n);
        // standard IP term integral
        for (size_type i = 0; i < lfsv_Sg_s.size(); i++)
        {
          r.accumulate(lfsv_Sg_s, i, term_penalty_sg * psi_Sg_s[i] * factor);
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

        
        double term_diffusion_T_1 = -Xc_conv_h*h_g * krN * omega_s * rho_g_s * (Kgradu_Pg_s) * n_F_local;
        double term_diffusion_T_2 = -Xc_conv_h*h_w * krW * omega_s * rho_w_s * (Kgradu_Pw_s * n_F_local);
        if (bctype[Indices::PVId_Pc] == Indices::BCId_neumann)
        {
          term_diffusion_T_1 = -bcvalue[Indices::PVId_Pc] * h_g;//(bcvalue[Indices::PVId_Pc]+bcvalue[Indices::PVId_Pw])
        }
        if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
        {
          term_diffusion_T_2 = -bcvalue[Indices::PVId_Pw] * h_w;
        }
        double term_diffusion_T_3 = -Xc_diff_h*kth_eff * (omega_s * (gradu_T_s * n_F_local));
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
          r.accumulate(lfsv_T_s, i, Xc_diff_h*term_nipg_T * (kth_eff * omega_s * (n_F * gradpsi_T_s[i])) * factor);
        }
        // standard IP term integral
        for (size_type i = 0; i < lfsv_T_s.size(); i++)
        {
          r.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
        }
      }

    } //End Quadrature Rule
  } //End of alpha_boundary

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
  //      auto A_s = property.A(cell_inside,local_inside);
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
  //          auto bctype = property.bctype(ig.intersection(),ip.position());
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
  //          auto b = property.b(cell_inside,iplocal_s);
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
