template<typename GV, typename Properties>
class ProblemInitialConditions
{
private:
	  const GV& gv;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;


public:

		//! construct from grid view
		ProblemInitialConditions (const GV& gv_, const Properties& property_)
		: gv( gv_ ),
			property(property_)
		{}

		/* Initial Conditions */
		std::vector< double >
		evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
			  	const Dune::FieldVector<double,dim>& xlocal) const {

			auto xglobal = element.geometry().global(xlocal);
			
			std::vector< double > icvalue(Indices::numOfPVs,0.);

			/******************************************************************************/
			// SATURATIONS
			double Sg = property.parameter.InitialSg(xglobal);
			double Sh = property.parameter.InitialSh(xglobal);
			double Sw = 1.-Sg-Sh;
			
			/******************************************************************************/
			// PRESSURES
			//double porosity = property.soil.SedimentPorosity(xglobal);
			double Pc = property.hydraulicProperty.suctionPressure(Sw, Sh)* property.hydraulicProperty.PcSF1(Sh)
						* property.characteristicValue.P_c; /*Pa*/
			double Pw = property.parameter.InitialPw(xglobal);  /*Pa*/
			double Pg = Pw+Pc;//  /*Pa*/
			
			double T = property.parameter.InitialT(xglobal)+273.15; /*K*/
			/******************************************************************************/
			// MOLE FRACTIONS
			double XC = property.parameter.InitialXC(xglobal);
			auto z_CH4 = property.eos.evaluateCompressibilityFactor( T,Pg );
		
			/******************************************************************************/
			
			double XCH4 = property.mixture.mole_x_CH4(T, Pg ,z_CH4, XC);
			double YH2O = property.mixture.mole_y_H2O(T, Pg ,z_CH4, XC);
			

			icvalue[Indices::PVId_Pw] = Pw ; //P_w + P_c ;
			icvalue[Indices::PVId_Sg] = Sg ;
			icvalue[Indices::PVId_Sh] = Sh ;
			icvalue[Indices::PVId_Pc] = Pc ;
			icvalue[Indices::PVId_T ] = T  ;
			icvalue[Indices::PVId_XCH4] = XCH4 ;
			icvalue[Indices::PVId_YH2O ] = YH2O  ;
			icvalue[Indices::PVId_C] = XC ;

			
		  /******************************************************************************/

		  return icvalue; /*ndim*/
	  	}

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};
