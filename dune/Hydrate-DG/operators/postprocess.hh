/*********************************************************
 * EVALUATE OUTPUT VARIABLES
 *********************************************************/

template< class GV,
		  class Params,
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
	Evaluation_Pw    *evaluation_Pw;
	Evaluation_Sg    *evaluation_Sg;
	Evaluation_Sh    *evaluation_Sh;
	Evaluation_T  	*evaluation_T;
	Evaluation_XCH4  *evaluation_XCH4;
	Evaluation_YH2O  *evaluation_YH2O;
	Evaluation_XC  *evaluation_XC;
	GFS_PP gfs_pp;
	U_pp *u_pp;

	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_D;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;

	double Xc_gravity;

//    typedef typename GV::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator LeafIterator;
	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PostProcess(	const GV& gv_,
					const Params& param_,
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
	  param(param_),
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
		typedef typename LFS_PP::template Child<Indices::SVId_Pg>::Type LFS_PP_Pg;
		const LFS_PP_Pg& lfs_pp_Pg = lfs_pp.template child<Indices::SVId_Pg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pw>::Type LFS_PP_Pw;
		const LFS_PP_Pw& lfs_pp_Pw = lfs_pp.template child<Indices::SVId_Pw>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pc>::Type LFS_PP_Pc;
		const LFS_PP_Pc& lfs_pp_Pc = lfs_pp.template child<Indices::SVId_Pc>();

		typedef typename LFS_PP::template Child<Indices::SVId_Sg>::Type LFS_PP_Sg;
		const LFS_PP_Sg& lfs_pp_Sg = lfs_pp.template child<Indices::SVId_Sg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Sh>::Type LFS_PP_Sh;
		const LFS_PP_Sh& lfs_pp_Sh = lfs_pp.template child<Indices::SVId_Sh>();
		typedef typename LFS_PP::template Child<Indices::SVId_Sw>::Type LFS_PP_Sw;
		const LFS_PP_Sw& lfs_pp_Sw = lfs_pp.template child<Indices::SVId_Sw>();

		typedef typename LFS_PP::template Child<Indices::SVId_T>::Type LFS_PP_T;
		const LFS_PP_T& lfs_pp_T = lfs_pp.template child<Indices::SVId_T>();

		typedef typename LFS_PP::template Child<Indices::SVId_XC>::Type LFS_PP_XC;
		const LFS_PP_XC& lfs_pp_XC = lfs_pp.template child<Indices::SVId_XC>();
		typedef typename LFS_PP::template Child<Indices::SVId_XCH4>::Type LFS_PP_XCH4;
		const LFS_PP_XCH4& lfs_pp_XCH4 = lfs_pp.template child<Indices::SVId_XCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_XH2O>::Type LFS_PP_XH2O;
		const LFS_PP_XH2O& lfs_pp_XH2O = lfs_pp.template child<Indices::SVId_XH2O>();
		typedef typename LFS_PP::template Child<Indices::SVId_YCH4>::Type LFS_PP_YCH4;
		const LFS_PP_YCH4& lfs_pp_YCH4 = lfs_pp.template child<Indices::SVId_YCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_YH2O>::Type LFS_PP_YH2O;
		const LFS_PP_YH2O& lfs_pp_YH2O = lfs_pp.template child<Indices::SVId_YH2O>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_rhow>::Type LFS_PP_rhow;
		const LFS_PP_rhow& lfs_pp_rhow = lfs_pp.template child<Indices::SVId_rhow>();
		typedef typename LFS_PP::template Child<Indices::SVId_rhog>::Type LFS_PP_rhog;
		const LFS_PP_rhog& lfs_pp_rhog = lfs_pp.template child<Indices::SVId_rhog>();

		typedef typename LFS_PP::template Child<Indices::SVId_K>::Type LFS_PP_K;
		const LFS_PP_K& lfs_pp_K = lfs_pp.template child<Indices::SVId_K>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_krw>::Type LFS_PP_krw;
		const LFS_PP_krw& lfs_pp_krw = lfs_pp.template child<Indices::SVId_krw>();
		typedef typename LFS_PP::template Child<Indices::SVId_krg>::Type LFS_PP_krg;
		const LFS_PP_krg& lfs_pp_krg = lfs_pp.template child<Indices::SVId_krg>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_muw>::Type LFS_PP_muw;
		const LFS_PP_muw& lfs_pp_muw = lfs_pp.template child<Indices::SVId_muw>();
		typedef typename LFS_PP::template Child<Indices::SVId_mug>::Type LFS_PP_mug;
		const LFS_PP_mug& lfs_pp_mug = lfs_pp.template child<Indices::SVId_mug>();

		typedef typename LFS_PP::template Child<Indices::SVId_zCH4>::Type LFS_PP_zCH4;
		const LFS_PP_zCH4& lfs_pp_zCH4 = lfs_pp.template child<Indices::SVId_zCH4>();

		typedef typename LFS_PP::template Child<Indices::SVId_por>::Type LFS_PP_por;
		const LFS_PP_por& lfs_pp_por = lfs_pp.template child<Indices::SVId_por>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_DH2O>::Type LFS_PP_DH2O;
		const LFS_PP_DH2O& lfs_pp_DH2O = lfs_pp.template child<Indices::SVId_DH2O>();
		typedef typename LFS_PP::template Child<Indices::SVId_DCH4>::Type LFS_PP_DCH4;
		const LFS_PP_DCH4& lfs_pp_DCH4 = lfs_pp.template child<Indices::SVId_DCH4>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_Pwsat>::Type LFS_PP_Pwsat;
		const LFS_PP_Pwsat& lfs_pp_Pwsat = lfs_pp.template child<Indices::SVId_Pwsat>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_HCH4>::Type LFS_PP_HCH4;
		const LFS_PP_HCH4& lfs_pp_HCH4 = lfs_pp.template child<Indices::SVId_HCH4>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_tau>::Type LFS_PP_tau;
		const LFS_PP_tau& lfs_pp_tau = lfs_pp.template child<Indices::SVId_tau>();

		typedef typename LFS_PP::template Child<Indices::SVId_Peq>::Type LFS_PP_Peq;
		const LFS_PP_Peq& lfs_pp_Peq = lfs_pp.template child<Indices::SVId_Peq>();

		typedef typename LFS_PP::template Child<Indices::SVId_HS>::Type LFS_PP_HS;
		const LFS_PP_HS& lfs_pp_HS = lfs_pp.template child<Indices::SVId_HS>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_Vwx>::Type LFS_PP_Vwx;
		const LFS_PP_Vwx& lfs_pp_Vwx = lfs_pp.template child<Indices::SVId_Vwx>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_Vwy>::Type LFS_PP_Vwy;
		const LFS_PP_Vwy& lfs_pp_Vwy = lfs_pp.template child<Indices::SVId_Vwy>();

		typedef typename LFS_PP::template Child<Indices::SVId_Vgx>::Type LFS_PP_Vgx;
		const LFS_PP_Vgx& lfs_pp_Vgx = lfs_pp.template child<Indices::SVId_Vgx>();
		
		typedef typename LFS_PP::template Child<Indices::SVId_Vgy>::Type LFS_PP_Vgy;
		const LFS_PP_Vgy& lfs_pp_Vgy = lfs_pp.template child<Indices::SVId_Vgy>();
		

		typedef Dune::PDELab::LFSIndexCache<LFS_PP> LFSCache_PP;
		LFSCache_PP lfs_cache_pp(lfs_pp);
		typedef typename U_pp::template LocalView<LFSCache_PP> VectorView_PP;
		VectorView_PP u_pp_view( (*u_pp) );

		typedef typename LFS_PP::template Child<0>::Type::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();
		auto gravity = param.parameter.g()  ; /* ndim */

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

			lfs_pp.bind(*self);
			lfs_cache_pp.update();
			u_pp_view.bind(lfs_cache_pp);
	        std::vector<double> ul_pp(lfs_pp.size());
			
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0) ;//barycentric position
	        auto cell_volume = geo.volume();
			auto ip_global = geo.global(cell_center_local);
      		auto ip_local = geo.local(ip_global);
			// auto cornervalue = ip_local;
			// for (int i=0; i <  geo.corners(); i++){
			// 	cornervalue = geo.corner(i);
			// 	ip_local = geo.local(cornervalue);
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
			// std::cout << "lfs_pp_Pg.localIndex(i)  ----------------   "<< std::endl;
			// std::cout << i << "  " << Pw << std::endl;
			// std::cout << "cornervalue  ----------------   "<< std::endl;
	        
	        for(int i = 0. ; i < lfs_pp.size() ; i++){
	        	ul_pp[i] = 0.;
	        }
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
			
			RF HS = (15.e0 * (Peq/ Xc_P-Pg)/abs(Peq/ Xc_P-Pg)+15.e0 ) ;
			// if((Peq-Pg*Xc_P)<= 1.e-2 && (Peq-Pg*Xc_P)>= -1.e-2){HS=1.;}//

			double eta = 1/BrooksCParams[1];
			auto Swe = param.hydraulicProperty.EffectiveSw(Sw,Sh,0., 0.);
			auto dPc_dSwe = param.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]); /*ndim */
			auto dSwe_dSw =  param.hydraulicProperty.dSwe_dSw(Sw,Sh,0., 0.);
			auto coeff_grad_Sw = dPc_dSwe * dSwe_dSw * param.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4]);

			auto dPcSF1_dSh =  param.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
			auto dSwe_dSh = param.hydraulicProperty.dSwe_dSh(Sw,Sh,0., 0.);
			auto coeff_grad_Sh = dPcSF1_dSh * std::pow( Swe , -eta ) * BrooksCParams[0] / Xc_P - Sg * coeff_grad_Sw *  dSwe_dSw  ;

      		auto Kgrad_Pg = Kgrad_Pw - coeff_grad_Sw * Kgrad_Sg + (coeff_grad_Sh ) * Kgrad_Sh;

			auto Vwx = - krw/(muw*Xc_mu)*Xc_K*(Kgrad_Pw[0]*Xc_P/Xc_x-rhow*Xc_rho*Kg[0]);
			auto Vwy = - krw/(muw*Xc_mu)*Xc_K*(Kgrad_Pw[1]*Xc_P/Xc_x-rhow*Xc_rho*Kg[1]);

			auto Vgx = - krg/(mug*Xc_mu)*Xc_K*(Kgrad_Pg[0]*Xc_P/Xc_x-rhog*Xc_rho*Kg[0]);
			auto Vgy = - krg/(mug*Xc_mu)*Xc_K*(Kgrad_Pg[1]*Xc_P/Xc_x-rhog*Xc_rho*Kg[1]);
			

			ul_pp[lfs_pp_Pg.localIndex(0)]	  = Pg*Xc_P ;
	        ul_pp[lfs_pp_Pw.localIndex(0)] 	  = Pw*Xc_P ;
	        ul_pp[lfs_pp_Pc.localIndex(0)]	  = Pc*Xc_P ;
	        ul_pp[lfs_pp_Sg.localIndex(0)]	  = Sg ;
	        ul_pp[lfs_pp_Sh.localIndex(0)]	  = Sh ;
	        ul_pp[lfs_pp_Sw.localIndex(0)] 	  = Sw ;
	        ul_pp[lfs_pp_T.localIndex(0)] 	  = T * Xc_T;
	        ul_pp[lfs_pp_XCH4.localIndex(0)]  = XCH4 ;
	        ul_pp[lfs_pp_XC.localIndex(0)]  = XC ;
	        ul_pp[lfs_pp_XH2O.localIndex(0)]  = XH2O ;
	        ul_pp[lfs_pp_YCH4.localIndex(0)]  = YCH4 ;
	        ul_pp[lfs_pp_YH2O.localIndex(0)]  = YH2O ;
	        ul_pp[lfs_pp_rhow.localIndex(0)]  = rhow*Xc_rho ;
	        ul_pp[lfs_pp_rhog.localIndex(0)]  = rhog*Xc_rho ;
	        ul_pp[lfs_pp_K.localIndex(0)] 	  = permeability*Xc_K ;
	        ul_pp[lfs_pp_krw.localIndex(0)]   = krw ;
	        ul_pp[lfs_pp_krg.localIndex(0)]   = krg ;
	        ul_pp[lfs_pp_muw.localIndex(0)]   = muw*Xc_mu ;
	        ul_pp[lfs_pp_mug.localIndex(0)]   = mug*Xc_mu ;
	        ul_pp[lfs_pp_zCH4.localIndex(0)]  = zCH4 ;
	        ul_pp[lfs_pp_por.localIndex(0)]   = porosity ;
	        ul_pp[lfs_pp_DH2O.localIndex(0)]  = DH2O_g*Xc_D ;
	        ul_pp[lfs_pp_DCH4.localIndex(0)]  = DCH4_w*Xc_D ;
	        ul_pp[lfs_pp_Pwsat.localIndex(0)] = Pwsat*Xc_P ;
	        ul_pp[lfs_pp_HCH4.localIndex(0)]  = Hch4*Xc_P ;
	        ul_pp[lfs_pp_tau.localIndex(0)]   = tau ;
	        ul_pp[lfs_pp_Peq.localIndex(0)]   = Peq ;
	        ul_pp[lfs_pp_HS.localIndex(0)]   = HS ;
	        ul_pp[lfs_pp_Vwx.localIndex(0)]   = Vwx ;
	        ul_pp[lfs_pp_Vwy.localIndex(0)]   = Vwy ;
	        ul_pp[lfs_pp_Vgx.localIndex(0)]   = Vgx ;
	        ul_pp[lfs_pp_Vgy.localIndex(0)]   = Vgy ;

			
			u_pp_view.write( ul_pp );
			u_pp_view.commit();
			u_pp_view.unbind();
			// }
			// exit(0);
			

		}//END:iterate over each volume
		
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
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;

			
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0) ;
	        auto cell_volume = geo.volume();
			auto ip_global = geo.global(cell_center_local);
      		auto ip_local = geo.local(ip_global);

	        // evaluation_Pw->updateSolution(cell);
	        
	        evaluation_Sg->updateSolution(cell);
	        
	        evaluation_Sh->updateSolution(cell);
	        
	        // evaluation_T->updateSolution(cell);
	        
	        // evaluation_YH2O->updateSolution(cell);
	        
	        evaluation_XCH4->updateSolution(cell);
	        
	        // evaluation_XC->updateSolution(cell);
			
			

		}//END:iterate over each volume
	}
};
