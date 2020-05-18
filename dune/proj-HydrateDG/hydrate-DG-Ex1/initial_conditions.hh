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
			double Sg = 0.;//property.parameter.InitialSg(xglobal);
			double Sh = 0.3;//property.parameter.InitialSh(xglobal);
			//double Sw = 1.-Sg-Sh;
			
			/******************************************************************************/
			// PRESSURES
			//double porosity = property.soil.SedimentPorosity(xglobal);
			//double Pc = property.hydraulicProperty.CapillaryPressure(element,xlocal,Sw,porosity)
			//			* property.characteristicValue.P_c; /*Pa*/
			double Pw = 2.e6;//property.parameter.InitialPw(xglobal);  /*Pa*/
			//double Pg = property.parameter.InitialPg(xglobal);  /*Pa*/
			
			/******************************************************************************/
			// MOLE FRACTIONS
			// auto S = property.parameter.ReferenceSalinity();
			// auto T = property.parameter.ReferenceTemperature();
			// auto zCH4 = property.eos.EvaluateCompressibilityFactor( T,Pg );
			// auto Xf = property.mixture.EquilibriumMoleFractions( T, Pg, S, zCH4);

			// auto xch4 = Xf[Indices::compId_XCH4];
			// auto yh2o = Xf[Indices::compId_YH2O];

			/******************************************************************************/
			//double Pw = 2.*1.e6; /*Pa*/
			
			double T = 4.+273.15;//property.parameter.InitialT(xglobal);; /*K*/
			double XCH4 = 0.;//property.parameter.InitialXCH4(xglobal);;
			double YH2O = 0.0005;//property.parameter.InitialYH2O(xglobal);;
			double XC = 5.5e-3;//property.parameter.InitialXC(xglobal);;
			
			double Pc = 8.48e4;//property.hydraulicProperty.suctionPressure(Sw,Sh) * property.hydraulicProperty.PcSF1(Sh);
			

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
