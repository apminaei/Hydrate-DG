/*
 * parameters.hh
 *
 * 
 */



template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	MeshParameters<PTree> mesh;
//	const double pi = boost::math::constants::pi<double>();
	const double pi = 3.14159265358979323846;
	const static int dim = MeshParameters<PTree>::dimension;
	CharacteristicValues X_c;
	double *time;
	double *dt;

	double T_t0;
	double Sg_t0;
	double Pw_t0;
	double Pg_t0;
	double Sh_t0;
	double XCH4_t0;
	double YH2O_t0;
	double XC_t0;
	double Sg_x0;
	double Pw_x0;
	double Sgin_x0;

	double time_initial_problem;
	double kd_initial_problem;
	double kf_initial_problem;
	double kd;
	double kf;
	int numMaterials;
	int numProps;
	std::vector<std::vector<double> > prop;

	double ref_salinity;
	double ref_saltconcentration;
	double ref_temperature;

	bool gravity_flag;
	double g_magnitude;

public:

  //! constructor
  Parameters (const PTree& ptree_)
  :ptree(ptree_),
   mesh(ptree_)
  {
		Sg_t0 = ptree.get("initial.Sg",(double)0.0);
		Pw_t0 = ptree.get("initial.Pw",(double)2.e6);
		Pg_t0 = ptree.get("initial.Pg",(double)2.0848e6);
		T_t0 = ptree.get("initial.T",(double)4.) ; // in Celcius
		Sh_t0 = ptree.get("initial.Sh",(double)0.3);
		YH2O_t0 = ptree.get("initial.YH2O",(double)0.0005);
		XCH4_t0 = ptree.get("initial.XCH4",(double)0.);
		XC_t0 = ptree.get("initial.XC",(double)5.5e-3);

		Pw_x0 = ptree.get("boundary.Pw_at_left",(double)2.e6);
		Sgin_x0 = ptree.get("boundary.Sg_at_inlet",(double)0.0);
		
		numMaterials = ptree.get("sediment.number_of_materials",(int)1);

		numProps = 8;
		prop = std::vector<std::vector<double> > (numMaterials,std::vector<double>(numProps, 0.));
		for(int n_mat=0; n_mat<numMaterials; n_mat++ ){
			std::string name = "sediment.material"+std::to_string(n_mat);
			prop[n_mat][0] = ptree.get(name+".por",	(double)0.3);
			prop[n_mat][1] = ptree.get(name+".K",	(double)1.e-12);
			prop[n_mat][2] = ptree.get(name+".pentry",(double)5.e4);
			prop[n_mat][3] = ptree.get(name+".lambda",(double)1.2);
			prop[n_mat][4] = ptree.get(name+".swr",	(double)0.);
			prop[n_mat][5] = ptree.get(name+".sgr",	(double)0.);
			prop[n_mat][6] = ptree.get(name+".m",	(double)3.);
			prop[n_mat][7] = ptree.get(name+".beta",(double)1.);
		}

		//reference state
		ref_salinity = ptree.get("reference_state.salinity",(double)0.);
		ref_saltconcentration = ref_salinity * (18.0/58.4); /*MolarMass_H2O/MolarMass_salt*/
		ref_temperature = (273.15+ptree.get("reference_state.temperature",(double)0.))/X_c.T_c;


		kd = ptree.get("hydrate_phase_change.dissociation_rate",(double)1.e-12);/*mol/m².Pa.s*/
		kf = ptree.get("hydrate_phase_change.formation_rate",(double)1.e-12);/*mol/m².Pa.s*/
		time_initial_problem = ptree.get("initial_problem.time_end",(double)5.); /*years*/
		time_initial_problem *= (1000.*364.*24.*60.*60.); /*convert to seconds*/
		kf_initial_problem = ptree.get("initial_problem.hydrate_formation_rate",(double)1.e-12);/*mol/m².Pa.s*/
		kd_initial_problem = ptree.get("initial_problem.hydrate_dissociation_rate",(double)1.e-12);/*mol/m².Pa.s*/

		//gravity
		gravity_flag = ptree.get("gravity.flag",(bool)true);
		g_magnitude = ptree.get("gravity.magnitude",(double)9.81);
  }

	/**********************************************************************
	 * INPUTS
	 **********
	 * z_domain : height of the computational domain [m]
	 * z_cells	: no. of cells along Z-axis
	 * 
	 *
	 *
	 *
	 **********************************************************************/

	//2. Initial Values

	double InitialSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sg_t0;
		return Sg;
	}

	double InitialPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_t0;
		return Pw; /* Pa */
	}

	double InitialPg(Dune::FieldVector< double,dim > xglobal) const {
		double Pg = Pg_t0;
		return Pg;/* Pa */
	}

	double InitialT(Dune::FieldVector< double,dim > xglobal) const {
		double T = T_t0;
		return T; /* Celcius */
	}
	double InitialSh(Dune::FieldVector< double,dim > xglobal) const {
		double Sh = Sh_t0;
		return Sh;
	}

	double InitialXCH4(Dune::FieldVector< double,dim > xglobal) const {
		double XCH4 = XCH4_t0;
		return XCH4;
	}
	double InitialYH2O(Dune::FieldVector< double,dim > xglobal) const {
		double YH2O = YH2O_t0;
		return YH2O;
	}

	double InitialXC(Dune::FieldVector< double,dim > xglobal) const {
		double XC = XC_t0;
		return XC;
	}

	//3. Boundary values
	double InletSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sgin_x0;
		return Sg;
	}

	double LeftPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_x0;
		return Pw; /* Pa */
	}


	//4. Material properties
	std::vector< std::vector<double> > layer_properties() const {
		return prop;
	}

	/**********************************************************************/
	/* REFERENCE STATE */
	double ReferenceSalinity() const {
		return ref_salinity; /*kg/kg*/
	}
	double ReferenceSaltConcentration() const {
		return ref_saltconcentration;
	}
	double ReferenceTemperature() const {
		return ref_temperature; /*K*/
	}

    /**********************************************************************/
	Dune::FieldVector<double,dim>
	SedimentVelocity ( double time, double dt ) const {

		Dune::FieldVector<double,dim> vs( 0. );
		vs[dim-1] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}

	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = g_magnitude;
		gravity[dim-1] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}

	// //HYDRATE REACTION KINETIC CONSTANTS
	// double HydrateDissociationRateConstant() const {
	// 	double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
	// 	if( time_new<time_initial_problem )
	// 		return kd_initial_problem;
	// 	else return kd;
	// 	/*mol/m².Pa.s*/
	// }
	// double HydrateFormationRateConstant() const {
	// 	double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
	// 	if( time_new<time_initial_problem ){
	// 		//std::cout<< "rho = " << kf_initial_problem << " T = " << time  << std::endl;
	// 		return kf_initial_problem;
	// 	}
	// 	else{ 
	// 		//std::cout<< "rho = " << kf << " T = " << time  << std::endl;
	// 		return kf;
	// 	}

		
			//exit(0);
		/*mol/m².Pa.s*/
	//}


	/**********************************************************************/

};


