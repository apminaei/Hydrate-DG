/*********************************************************
 * EVALUATE OUTPUT VARIABLES
 *********************************************************/

template< class GV,
		  class Params, class FEM_P,  class FEM_Sg,  class FEM_Sh,  class FEM_T,  class FEM_X,  class FEM_Y,  class FEM_C,
		  class Evaluation_Pw,
		  class Evaluation_Sg,
		  class Evaluation_Sh,
		  class Evaluation_T,
		  class Evaluation_XCH4,
		  class Evaluation_YH2O,
		  class Evaluation_XC,
		  class GFS_PP, typename U_pp >
class PostProcess{
private:
	const GV&  gv;
	const Params& param; 
	const FEM_P &fem_P;
	const FEM_Sg &fem_Sg;
	const FEM_Sh &fem_Sh;
	const FEM_T &fem_T;
	const FEM_X &fem_x;
	const FEM_Y &fem_y;
	const FEM_C &fem_c;
	Evaluation_Pw    *evaluation_Pw;
	Evaluation_Sg    *evaluation_Sg;
	Evaluation_Sh    *evaluation_Sh;
	Evaluation_T  	*evaluation_T;
	Evaluation_XCH4  *evaluation_XCH4;
	Evaluation_YH2O  *evaluation_YH2O;
	Evaluation_XC  *evaluation_XC;
	GFS_PP gfs_pp;
	U_pp *u_pp;
	double eps=1e-6;
	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_D;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;

	double Xc_gravity;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PostProcess(	const GV& gv_,
					const Params& param_, const FEM_P& fem_P_, const FEM_Sg& fem_Sg_, const FEM_Sh& fem_Sh_,
					const FEM_T& fem_T_, const FEM_X& fem_x_, const FEM_Y& fem_y_, const FEM_C& fem_c_, 
					Evaluation_Pw 	*evaluation_Pw_,
					Evaluation_Sg 	*evaluation_Sg_,
					Evaluation_Sh 	*evaluation_Sh_,
					Evaluation_T *evaluation_T_,
					Evaluation_XCH4 *evaluation_XCH4_,
					Evaluation_YH2O *evaluation_YH2O_,
					Evaluation_XC *evaluation_XC_,
					GFS_PP	gfs_pp_,
					U_pp	*u_pp_ )
	: gv(gv_),
	  param(param_), fem_P(fem_P_), fem_Sg(fem_Sg_), fem_Sh(fem_Sh_), fem_T(fem_T_), fem_x(fem_x_), fem_y(fem_y_), fem_c(fem_c_),
	  evaluation_Pw(evaluation_Pw_),
	  evaluation_Sg(evaluation_Sg_),
	  evaluation_Sh(evaluation_Sh_),
	  evaluation_T(evaluation_T_),
	  evaluation_XCH4(evaluation_XCH4_),
	  evaluation_YH2O(evaluation_YH2O_),
	  evaluation_XC(evaluation_XC_),
	  gfs_pp(gfs_pp_),
	  u_pp(u_pp_)
	{
		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_D		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T 		= param.characteristicValue.T_c;
		  Xc_t 		= param.characteristicValue.t_c;
		  Xc_x 		= param.characteristicValue.x_c;
		  Xc_gravity= param.characteristicValue.X_gravity;
	  }

	virtual ~PostProcess()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS_PP > LFS_PP;
		LFS_PP lfs_pp(gfs_pp);
		
		typedef Dune::PDELab::LFSIndexCache<LFS_PP> LFSCache_PP;
		LFSCache_PP lfs_cache_pp(lfs_pp);

		typedef typename LFS_PP::template Child<0>::Type::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();
		auto gravity = param.parameter.g()  ; /* ndim */
		auto dof= gv.size(param.mesh.getdim());
		std::vector<double> repeat(dof,0.);
		Dune::PDELab::Backend::native(*u_pp) = 0.;
		
		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
			
			for (int i=0; i <  geo.corners(); ++i){
				auto cornervalue = geo.corner(i);
				
				auto ip_local = geo.local(cornervalue);
				auto globalVertexIndex = indexSet.subIndex(cell , i , dim );
				repeat[globalVertexIndex]+=1.;
				
				Dune::FieldVector<RF, dim> grad_Pw(0.0);
				Dune::FieldVector<RF, dim> Kgrad_Pw(0.0);
				Dune::FieldVector<RF, dim> grad_Sg(0.0);
				Dune::FieldVector<RF, dim> Kgrad_Sg(0.0);
				Dune::FieldVector<RF, dim> grad_Sh(0.0);
				Dune::FieldVector<RF, dim> Kgrad_Sh(0.0);
				Dune::FieldVector<RF, dim> Kg(0.0);
				
				RF Pw=0.;
				evaluation_Pw->evalFunction(cell, ip_local, &Pw);
				evaluation_Pw->evalGradient(cell, ip_local, &grad_Pw);
				RF Sg=0.;
				evaluation_Sg->evalFunction(cell,ip_local,&Sg);
				evaluation_Sg->evalGradient(cell, ip_local, &grad_Sg);
				RF Sh=0.;
				evaluation_Sh->evalFunction(cell,ip_local,&Sh);
				evaluation_Sh->evalGradient(cell, ip_local, &grad_Sh);
				RF T=0.;
				evaluation_T->evalFunction(cell,ip_local,&T);
				RF YH2O=0.;
				evaluation_YH2O->evalFunction(cell,ip_local,&YH2O);
				RF XCH4=0.;
				evaluation_XCH4->evalFunction(cell,ip_local,&XCH4);
				RF XC=0.;
				evaluation_XC->evalFunction(cell,ip_local,&XC);
				
				
				
				auto BrooksCParams = param.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
				RF porosity = param.soil.SedimentPorosity( cell,ip_local );
				RF permeability = param.soil.SedimentPermeability( cell,ip_local )
					* param.hydraulicProperty.PermeabilityScalingFactor( cell,ip_local, Sh, porosity );

				auto K = param.soil.SedimentPermeabilityTensor(cell, ip_local)
						* param.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, porosity ); /*ndim K from soil.hh*/
				
				K.mv(gravity, Kg);

				// compute K * gradient of Pw
				K.mv(grad_Pw, Kgrad_Pw);

				// compute K * gradient of Sg
				K.mv(grad_Sg, Kgrad_Sg);

				// compute K * gradient of Sh
				K.mv(grad_Sh, Kgrad_Sh);

				RF S = XC * (param.salt.MolarMass()/param.gas.MolarMass());
				// RF Xc = param.parameter.ReferenceSaltConcentration();
				RF T_ref = param.parameter.ReferenceTemperature()/Xc_T; /*ndim*/
				RF Sw = 1.- Sg - Sh;
				RF Pc = param.hydraulicProperty.CapillaryPressure( cell,ip_local, Sw, Sh, porosity );
				RF Pg = Pw + Pc;
				
				RF zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );
				
				RF XH2O = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,XC);
				RF YCH4 = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,XC,zCH4);
				
				RF rhog = param.gas.Density( T*Xc_T, Pg*Xc_P, zCH4 );
				RF rhow = param.water.Density( T*Xc_T, Pw*Xc_P, S );
				
				RF krg = param.hydraulicProperty.krg( cell, ip_local, Sw, Sh );
				RF krw = param.hydraulicProperty.krw( cell, ip_local, Sw, Sh );
				
				RF mug = param.gas.DynamicViscosity( T*Xc_T, Pg*Xc_P );
				RF muw = param.water.DynamicViscosity( T*Xc_T, Pw*Xc_P, S );
				
				RF tau = param.soil.Tortuosity( porosity );
				RF DH2O_g	= param.mixture.DiffCoeffH2OInGas( T*Xc_T,Pg*Xc_P );
				RF DCH4_w	= param.mixture.DiffCoeffCH4InLiquid( T*Xc_T,Pw*Xc_P );

				RF Pwsat = param.water.SaturatedVaporPressure( T*Xc_T,S );
				RF Hch4 = param.gas.SolubilityCoefficient(T*Xc_T,S);

				RF Coeff = 16.32;//4.4824;//14.543;//
				RF Peq = 1.e3 * exp( 38.592 - (8533.8/ (T * Xc_T) ) + Coeff*S );
				
				RF HS =  ((Peq/ Xc_P-Pg)/abs(Peq/ Xc_P-Pg) +1.) /2.;
				// if((Peq-Pg*Xc_P)<= 1.e-2 && (Peq-Pg*Xc_P)>= -1.e-2){HS=1.;}//

				double eta = 1/BrooksCParams[1];
				auto Swe = param.hydraulicProperty.EffectiveSw(Sw,Sh,0., 0.);
				auto dPc_dSwe = param.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]); /*ndim */
				auto dSwe_dSw =  param.hydraulicProperty.dSwe_dSw(Sw,Sh,0., 0.);
				auto coeff_grad_Sw = dPc_dSwe * dSwe_dSw * param.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]);

				auto dPcSF1_dSh =  param.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
				auto dSwe_dSh = param.hydraulicProperty.dSwe_dSh(Sw,Sh,0., 0.);
				auto coeff_grad_Sh = dPcSF1_dSh * std::pow( Swe , -eta ) * BrooksCParams[0] / Xc_P - Sg * coeff_grad_Sw *  dSwe_dSw  ;

							//   + param.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]) * dPc_dSwe * (dSwe_dSh);//
				auto Kgrad_Pg = Kgrad_Pw - coeff_grad_Sw * Kgrad_Sg + coeff_grad_Sh * Kgrad_Sh;

				auto Vwx = - krw/(muw*Xc_mu)*Xc_K*(Kgrad_Pw[0]*Xc_P/Xc_x-rhow*Xc_rho*Kg[0]);
				auto Vwy = - krw/(muw*Xc_mu)*Xc_K*(Kgrad_Pw[1]*Xc_P/Xc_x-rhow*Xc_rho*Kg[1]);

				auto Vgx = - krg/(mug*Xc_mu)*Xc_K*(Kgrad_Pg[0]*Xc_P/Xc_x-rhog*Xc_rho*Kg[0]);
				auto Vgy = - krg/(mug*Xc_mu)*Xc_K*(Kgrad_Pg[1]*Xc_P/Xc_x-rhog*Xc_rho*Kg[1]);
				

				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Pg]	  += Pg*Xc_P;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Pw] 	  += Pw*Xc_P ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Pc]	  += Pc*Xc_P;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Sg]	  += Sg;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Sh]	  += Sh ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Sw]	  += Sw;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_T] 	  += T * Xc_T;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_XCH4]  += XCH4;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_XC]  += XC ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_XH2O]  += XH2O ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_YCH4] += YCH4;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_YH2O]  += YH2O;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_rhow]  += rhow*Xc_rho;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_rhog]  += rhog*Xc_rho;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_K]	  += permeability*Xc_K;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_krw]   += krw;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_krg]   += krg ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_muw]   += muw*Xc_mu;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_mug]   += mug*Xc_mu ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_zCH4]  += zCH4;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_por]  += porosity ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_DH2O]  += DH2O_g*Xc_D;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_DCH4]  += DCH4_w*Xc_D;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Pwsat] += Pwsat*Xc_P;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_HCH4]  += Hch4*Xc_P;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_tau]   += tau;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Peq]   += Peq;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_HS]   += HS ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Vwx]   += Vwx;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Vwy]   += Vwy ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Vgx]   += Vgx ;
				Dune::PDELab::Backend::native(*u_pp)[globalVertexIndex][Indices::SVId_Vgy]   += Vgy ;
			
			}
			
		}//END:iterate over each volume

		for(int i=0; i<gv.size(param.mesh.getdim()); ++i){
			auto corners_number = 1./repeat[i];
			Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Pg]	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Pw] 	  *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Pc]	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Sg]	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Sh]	  *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Sw]	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_T] 	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_XCH4]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_XC]  *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_XH2O]  *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_YCH4] *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_YH2O]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_rhow]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_rhog]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_K]	  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_krw]   *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_krg]   *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_muw]   *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_mug]   *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_zCH4]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_por]  *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_DH2O]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_DCH4]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Pwsat] *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_HCH4]  *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_tau]   *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Peq]   *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_HS]   *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Vwx]   *= corners_number;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Vwy]   *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Vgx]   *= corners_number ;
	        Dune::PDELab::Backend::native(*u_pp)[i][Indices::SVId_Vgy]   *= corners_number ;
		}
		
	}
	void evalProjectedSolution() {

		
		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
	        // evaluation_Pw->updateSolution(cell);
	        if(fem_Sg.polynomial()==Dune::PDELab::QkDGBasisPolynomial::lagrange){
	        	evaluation_Sg->updateSolution(cell);
			}else{
				evaluation_Sg->updateLegendreSolution(cell);
			}
	        if(fem_Sh.polynomial()==Dune::PDELab::QkDGBasisPolynomial::lagrange){
	        	evaluation_Sh->updateSolution(cell);
				// std::cout << "lagrange   " << std::endl; 
			}else{
				evaluation_Sh->updateLegendreSolution(cell);
				// std::cout << "legendre   " << std::endl;
			}
	        // evaluation_T->updateSolution(cell);
	        
	        // evaluation_YH2O->updateSolution(cell);
	        if(fem_x.polynomial()==Dune::PDELab::QkDGBasisPolynomial::lagrange){
	        	evaluation_XCH4->updateSolution(cell);
			}else{
				evaluation_XCH4->updateLegendreSolution(cell);
			}
	       

		}//END:iterate over each volume
	}
};
