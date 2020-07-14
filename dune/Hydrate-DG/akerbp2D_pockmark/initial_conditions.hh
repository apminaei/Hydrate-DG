/* All Values are dimensional and transfer to nondim in Initial.hh */
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

			auto prop_L = property.parameter.layer_properties();
			auto Pentry = prop_L[0][2]; // Characteristic value for Pressure
			std::vector< double > icvalue(Indices::numOfPVs,0.);

			/******************************************************************************/
			// SATURATIONS
			double Sg = property.parameter.InitialSg(xglobal);
			double Sh = property.parameter.InitialSh(xglobal);
			double Sw = 1.-Sg-Sh;
			
			/******************************************************************************/
			// PRESSURES
			auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(element, xlocal);
			double por = property.soil.SedimentPorosity(element, xlocal);
			double Pc = property.hydraulicProperty.CapillaryPressure(element, xlocal, Sw, Sh, por)
						* property.hydraulicProperty.PcSF1(Sh, BrooksCParams[1], BrooksCParams[4])
						* property.characteristicValue.P_c; /*Pa*/
			double Pw = property.parameter.InitialPw(xglobal);  /*Pa*/
			double Pg = Pw+Pc;//  /*Pa*/
			//std::cout << " Pg = " << Pg << std::endl;
			double T = property.parameter.InitialT(xglobal) + 273.15; /*K*/
			/******************************************************************************/
			// MOLE FRACTIONS
			double XC = property.parameter.InitialXC(xglobal);
			auto z_CH4 = property.eos.EvaluateCompressibilityFactor( T,Pg );
		
			/******************************************************************************/
			auto VLEqui = property.mixture.EquilibriumMoleFractions( T/*K*/, Pg/*Pa*/, Sg, Sw, XC,  z_CH4 );
			double XCH4 = property.parameter.InitialXCH4(xglobal); //0.0
			double YH2O = VLEqui[Indices::compId_YH2O];// change 0.0

			icvalue[Indices::PVId_Pw] = Pw ; //P_w + P_c ;
			icvalue[Indices::PVId_Sg] = Sg ;
			icvalue[Indices::PVId_Sh] = Sh ;
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
