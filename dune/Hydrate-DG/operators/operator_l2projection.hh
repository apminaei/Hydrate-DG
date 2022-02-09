/*
 * operator_l2projetcion.hh
 *
 *  Created on: Nov 14, 2019
 *      Author: demo
 */

template <class GV, class PARAMS, class GFS, class U,
			class Evaluation_Pw,
		  class Evaluation_Sg,
		  class Evaluation_Sh,
		  class Evaluation_T,
		  class Evaluation_XCH4,
		  class Evaluation_YH2O,
		  class Evaluation_XC>
class LocalOperatorPROJ :  
	public Dune::PDELab::NumericalJacobianApplyVolume	<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::NumericalJacobianVolume		<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::NumericalJacobianApplySkeleton	<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::NumericalJacobianSkeleton		<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary	<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::NumericalJacobianBoundary		<LocalOperatorPROJ<GV,PARAMS,GFS,U, Evaluation_Pw, Evaluation_Sg, Evaluation_Sh, Evaluation_T, Evaluation_XCH4, Evaluation_YH2O, Evaluation_XC> >,
	public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
		    const GV& gv;
		    const PARAMS& param;
			GFS gfs;
			Evaluation_Pw    *evaluation_Pw;
			Evaluation_Sg    *evaluation_Sg;
			Evaluation_Sh    *evaluation_Sh;
			Evaluation_T  	*evaluation_T;
			Evaluation_XCH4  *evaluation_XCH4;
			Evaluation_YH2O  *evaluation_YH2O;
			Evaluation_XC  *evaluation_XC;
			U *ufv;
			unsigned int intorder ;
			double epsilon;
			double Xc_K ;
			double Xc_mu ;
			double Xc_rho;
			double Xc_kth;
			double Xc_C;
			double Xc_D;
			double Xc_P;
			double Xc_T;
			double Xc_t;
			double Xc_x;

public:
			// pattern assembly flags
			enum { doPatternVolume 	 = true  };
			enum { doPatternSkeleton = false };

			// residual assembly flags
			enum { doAlphaVolume  	= true 	};
			enum { doLambdaVolume  	= true 	};
			enum { doAlphaSkeleton  = false };
			enum { doAlphaBoundary  = false };
			enum { doLambdaBoundary = false };

			typedef typename GV::IndexSet IndexSet;

			typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
		    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
			typedef typename U::template LocalView<LFSCache> VectorView;

			// constructor stores parameters
			LocalOperatorPROJ(	const GV& gv_,
								const PARAMS& param_,
								GFS 	gfs_,
								U		*ufv_, 
								Evaluation_Pw 	*evaluation_Pw_,
								Evaluation_Sg 	*evaluation_Sg_,
								Evaluation_Sh 	*evaluation_Sh_,
								Evaluation_T *evaluation_T_,
								Evaluation_XCH4 *evaluation_XCH4_,
								Evaluation_YH2O *evaluation_YH2O_,
								Evaluation_XC *evaluation_XC_,
								unsigned int 	intorder_ = 2,
								double 			epsilon_ = 1.e-6)
			: gv(gv_),
			  param(param_),
			  gfs(gfs_),
			  ufv(ufv_),
			  evaluation_Pw(evaluation_Pw_),
			evaluation_Sg(evaluation_Sg_),
			evaluation_Sh(evaluation_Sh_),
			evaluation_T(evaluation_T_),
			evaluation_XCH4(evaluation_XCH4_),
			evaluation_YH2O(evaluation_YH2O_),
			evaluation_XC(evaluation_XC_),
			  intorder( intorder_ ),
			  epsilon( epsilon_ )
			{
				  Xc_K 		= param.characteristicValue.permeability_c;
				  Xc_mu 	= param.characteristicValue.viscosity_c;
				  Xc_rho 	= param.characteristicValue.density_c;
				  Xc_kth 	= param.characteristicValue.thermalconductivity_c;
				  Xc_C 		= param.characteristicValue.specificheat_c;
				  Xc_D		= param.characteristicValue.dispersivity_c;
				  Xc_P 		= param.characteristicValue.P_c;
				  Xc_T 		= param.characteristicValue.T_c;
				  Xc_t 		= param.characteristicValue.t_c;
				  Xc_x 		= param.characteristicValue.x_c;
			 }

			// volume integral depending on test and ansatz functions
			template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
			void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
			{
				// define types
				using RF = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
				using RangeType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
				using JacobianType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

				// dimensions
				const int dim = EG::Entity::dimension;

				// Get cell
				const auto& cell = eg.entity();

				// Get geometry
				auto geo = eg.geometry();

				// evaluate diffusion tensor at cell center, assume it is constant over elements
				auto ref_el = referenceElement(geo);
				auto localcenter = ref_el.position(0,0);

				// Transformation matrix
				typename EG::Geometry::JacobianInverseTransposed jac;

				// loop over quadrature points
				for (const auto& ip : quadratureRule(geo,intorder))
				{
					// evaluate basis functions
					std::vector< std::vector<RangeType> > phi( Indices::numOfPVs );
			        for( int i=0; i<phi.size(); i++ ){
			        	 phi[i] = std::vector<RangeType> ( lfsu.child(i).size() );
			        	 lfsu.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),phi[i]);
			        }

					std::vector< std::vector<RangeType> > psi( Indices::numOfPVs );
			        for( int i=0; i<psi.size(); i++ ){
			        	 psi[i] = std::vector<RangeType> ( lfsv.child(i).size() );
			        	 lfsv.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),psi[i]);
			        }

					auto ip_global = geo.global(ip.position());

					// evaluate u (u->phase pressures)
			        std::vector<double> u( Indices::numOfPVs ,0.);
			        for(int i=0; i<u.size(); i++ ){
			        	for (int j=0; j<lfsu.child(i).size(); j++)
			        		u[i] += x(lfsu.child(i),j) * phi[i][j];
			        }

					// evaluate gradient of basis functions
			        // std::vector< std::vector<JacobianType> > jsu( Indices::numOfPVs );
			        // for( int i=0; i<jsu.size(); i++ ){
			        // 	 jsu[i] = std::vector<JacobianType> ( lfsu.child(i).size() );
			        // 	 lfsu.child(i).finiteElement().localBasis().evaluateJacobian(ip.position(),jsu[i]);
			        // }

			        // std::vector< std::vector<JacobianType> > jsv( Indices::numOfPVs );
			        // for( int i=0; i<jsv.size(); i++ ){
			        // 	 jsv[i] = std::vector<JacobianType> ( lfsv.child(i).size() );
			        // 	 lfsv.child(i).finiteElement().localBasis().evaluateJacobian(ip.position(),jsv[i]);
			        // }

					// // transform gradients of shape functions to real element
					// jac = geo.jacobianInverseTransposed(ip.position());

					// // evaluade gradients of shape fncs.
			        // std::vector< std::vector < Dune::FieldVector<RF,dim > > > gradphi( Indices::numOfPVs );
			        // for( int i=0; i<gradphi.size(); i++ ){
			        // 	gradphi[i] = std::vector<Dune::FieldVector<RF,dim>> ( lfsu.child(i).size() );
			        // 	for (int j=0; j<lfsu.child(i).size(); j++){
			        // 	          jac.mv(jsu[i][j][0],gradphi[i][j]);
			        // 	}
			        // }

			        // std::vector< std::vector < Dune::FieldVector<RF,dim > > > gradpsi( Indices::numOfPVs );
			        // for( int i=0; i<gradpsi.size(); i++ ){
			        // 	gradpsi[i] = std::vector<Dune::FieldVector<RF,dim>> ( lfsv.child(i).size() );
			        // 	for (int j=0; j<lfsv.child(i).size(); j++){
			        // 	          jac.mv(jsv[i][j][0],gradpsi[i][j]);
			        // 	}
			        // }

					// // compute gradient of u (phase pressures)
			        // std::vector< Dune::FieldVector<RF,dim > > gradu( Indices::numOfPVs );
					// for( int i=0; i<gradu.size(); i++ ){
					// 	gradu[i] =  Dune::FieldVector<RF,dim>(0.0);
					// 	for (int j=0; j<lfsu.child(i).size(); j++){
					// 		gradu[i].axpy(x(lfsu.child(i),j),gradphi[i][j]);
					// 	}
					// }

					// integrate
					// FV --> FEM
					//|| u_FV - u_FE || --> Min
					//
					RF factor = ip.weight() * geo.integrationElement(ip.position());
					for (int i=0; i<Indices::numOfPVs; i++){
						double tmp = 0.;
						for (int j=0; j<lfsv.child(i).size(); j++){
							tmp =  u[i] * psi[i][j] ;//+ gradu[i] * gradpsi[i][j];
							r.accumulate(lfsv.child(i),j,tmp*factor );
						}
					}

					

				}//End Quadrature Rule

			}


			// volume integral depending only on test functions
			template<typename EG, typename LFSV, typename R>
			void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
			{
				// define types
				using RF = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
				using RangeType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
				using JacobianType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

				// dimensions
				const int dim = EG::Entity::dimension;

				// Get cell
				const auto& cell = eg.entity();

				// Get geometry
				auto geo = eg.geometry();

				// evaluate diffusion tensor at cell center, assume it is constant over elements
				// auto ref_el = referenceElement(geo);
				// auto localcenter = ref_el.position(0,0);
		        // auto globalcenter = cell.geometry().global(localcenter);

				
				
				


				// loop over quadrature points
				for (const auto& ip : quadratureRule(geo,intorder))
				{
					// evaluate shape functions
					std::vector<std::vector<RangeType>> psi( Indices::numOfPVs );
			        for( int i=0; i<psi.size(); i++ ){
			        	 psi[i] = std::vector<RangeType> ( lfsv.child(i).size() );
			        	 lfsv.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),psi[i]);
			        }

					auto ip_global = geo.global(ip.position());

					

					// evaluate gradient of basis functions
			        // std::vector< std::vector<JacobianType> > jsv( Indices::numOfPVs );
			        // for( int i=0; i<jsv.size(); i++ ){
			        // 	 jsv[i] = std::vector<JacobianType> ( lfsv.child(i).size() );
			        // 	 lfsv.child(i).finiteElement().localBasis().evaluateJacobian(ip.position(),jsv[i]);
			        // }

					// // transform gradients of shape functions to real element
					// auto jac = geo.jacobianInverseTransposed(ip.position());

					// // evaluade gradients of shape fncs.
			        // std::vector< std::vector < Dune::FieldVector<RF,dim > > > gradpsi( Indices::numOfPVs );
			        // for( int i=0; i<gradpsi.size(); i++ ){
			        // 	gradpsi[i] = std::vector<Dune::FieldVector<RF,dim>> ( lfsv.child(i).size() );
			        // 	for (int j=0; j<lfsv.child(i).size(); j++){
			        // 	          jac.mv(jsv[i][j][0],gradpsi[i][j]);
			        // 	}
			        // }

					RF Pw=0.;
					evaluation_Pw->evalFunction(cell, ip.position(), &Pw);
					RF Sg=0.;
					evaluation_Sg->evalFunction(cell,ip.position(),&Sg);
					RF Sh=0.;
					evaluation_Sh->evalFunction(cell,ip.position(),&Sh);
					RF T=0.;
					evaluation_T->evalFunction(cell,ip.position(),&T);
					RF YH2O=0.;
					evaluation_YH2O->evalFunction(cell,ip.position(),&YH2O);
					RF XCH4=0.;
					evaluation_XCH4->evalFunction(cell,ip.position(),&XCH4);
					RF XC=0.;
					evaluation_XC->evalFunction(cell,ip.position(),&XC);
					
					// Dune::FieldVector<RF, dim> grad_Pw(0.0);
					// evaluation_Pw->evalGradient(cell, ip.position(), &grad_Pw);
					// Dune::FieldVector<RF, dim> grad_Sg(0.0);
					// evaluation_Sg->evalGradient(cell,ip.position(),&grad_Sg);
					// Dune::FieldVector<RF, dim> grad_Sh(0.0);
					// evaluation_Sh->evalGradient(cell,ip.position(),&grad_Sh);
					// Dune::FieldVector<RF, dim> grad_T(0.0);
					// evaluation_T->evalGradient(cell,ip.position(),&grad_T);
					// Dune::FieldVector<RF, dim> grad_YH2O(0.0);
					// evaluation_YH2O->evalGradient(cell,ip.position(),&grad_YH2O);
					// Dune::FieldVector<RF, dim> grad_XCH4(0.0);
					// evaluation_XCH4->evalGradient(cell,ip.position(),&grad_XCH4);
					// Dune::FieldVector<RF, dim> grad_XC(0.0);
					// evaluation_XC->evalGradient(cell,ip.position(),&grad_XC);

					// integrate
					RF factor = ip.weight() * geo.integrationElement(ip.position());
					RF tmp = 0;
					for (int j=0; j<lfsv.child(Indices::PVId_Pw).size(); j++){
						//tmp = grad_Pw * gradpsi[Indices::PVId_Pw][j];
						r.accumulate(lfsv.child(Indices::PVId_Pw),j,( -Pw * psi[Indices::PVId_Pw][j]  )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_Sg).size(); j++){
						//tmp = grad_Sg * gradpsi[Indices::PVId_Sg][j];
						r.accumulate(lfsv.child(Indices::PVId_Sg),j,( -Sg * psi[Indices::PVId_Sg][j]    )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_Sh).size(); j++){
						//tmp = grad_Sh * gradpsi[Indices::PVId_Sh][j];
						r.accumulate(lfsv.child(Indices::PVId_Sh),j,( -Sh * psi[Indices::PVId_Sh][j]   )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_T).size(); j++){
						//tmp = grad_T * gradpsi[Indices::PVId_T][j];
						r.accumulate(lfsv.child(Indices::PVId_T),j,( -T * psi[Indices::PVId_T][j]   )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_C).size(); j++){
						//tmp = grad_XC * gradpsi[Indices::PVId_C][j];
						r.accumulate(lfsv.child(Indices::PVId_C),j,( -XC * psi[Indices::PVId_C][j]   )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_XCH4).size(); j++){
						//tmp = grad_XCH4 * gradpsi[Indices::PVId_XCH4][j];
						r.accumulate(lfsv.child(Indices::PVId_XCH4),j,( -XCH4 * psi[Indices::PVId_XCH4][j]   )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::PVId_YH2O).size(); j++){
						//tmp = grad_YH2O * gradpsi[Indices::PVId_YH2O][j];
						r.accumulate(lfsv.child(Indices::PVId_YH2O),j,( -YH2O * psi[Indices::PVId_YH2O][j]   )*factor);
					}

				}//END: quadrature rule

			}//END: lambda_volume

	  };

