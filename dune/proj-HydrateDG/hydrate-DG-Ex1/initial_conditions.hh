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
			  	const Dune::FieldVector<double,dim>& xglobal) const {

			//auto xglobal = element.geometry().global(xlocal);
			
			std::vector< double > icvalue(Indices::numOfPVs,0.);

			/******************************************************************************/
			// SATURATIONS
			double Sg = property.parameter.InitialSg(xglobal);
			double Sh = property.parameter.InitialSh(xglobal);
			double Sw = 1.-Sg-Sh;

			/******************************************************************************/
			// PRESSURES
			//double porosity = property.soil.SedimentPorosity(xglobal);
			//double Pc = property.hydraulicProperty.CapillaryPressure(element,xlocal,Sw,porosity)
			//			* property.characteristicValue.P_c; /*Pa*/
			double Pw = property.parameter.InitialPw(xglobal);  /*Pa*/
			//double Pg = Pw + Pc; /*Pa*/
			
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
			
			double T = property.parameter.InitialT(xglobal);; /*K*/
			double XCH4 = property.parameter.InitialXCH4(xglobal);;
			double YH2O = property.parameter.InitialYH2O(xglobal);;
			double XC = property.parameter.InitialXC(xglobal);;
			HydraulicProperties hydraulicProperty;
			double Pc = hydraulicProperty.suctionPressure(Sw,Sh) * hydraulicProperty.PcSF1(Sh);
			double Pg = Pw + Pc; /*Pa*/

			icvalue[Indices::PVId_Pw] = Pg ; //P_w + P_c ;
			icvalue[Indices::PVId_Sg] = Sw ;
			icvalue[Indices::PVId_Sh] = Sh ;
			//icvalue[Indices::PVId_Pc] = Pc ;
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
