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
	double *time;
	double *dt;

	MeshParameters<PTree> mesh;
	const double pi = 3.14159265358979323846;
	const static int dim = MeshParameters<PTree>::dimension;
	CharacteristicValues X_c;

	double Zmax;
	double Xmax;
	double FGP_thickness;
	double FGP_width;
	double FGP_offset;
	double BSR_depth;
	double fracture_diameter;
	double fracture_depth;
	std::vector<double> z;
	std::vector<std::vector<double> > prop;
	std::vector<std::vector<double> > prop_initial_problem;
	std::vector<double> prop_fracture;
	bool gravity_flag;
	double gravity_magnitude;
	double burial_rate;
	double kd;
	double kf;
	double time_initial_problem;
	double kd_initial_problem;
	double kf_initial_problem;
	double salinity_initial_problem;
	double FGP_Sg;
	double PSF_Pw;
	double PSF_T;
	double thermal_gradient;
	double BSR_z;
	
public:

  //! constructor
  Parameters (const PTree& ptree_, double *time_/*ndim*/, double *dt_/*ndim*/)
  :ptree(ptree_),
   time(time_),
   dt(dt_),
   mesh(ptree_)
  {
	/* layer numbering: example with numLayers=3
	Zmax___________
			|
			|
			|	layer0
			|
		____|______ z0
			|
			|
			|	layer1
			|
		____|______ z1
			|
			|
			|	layer2
			|
	z=0	____|______ z2
	*/

	Zmax=ptree.get("grid.LZ",(double)1000.); /*m*/
	Xmax=ptree.get("grid.LX",(double)1000.); /*m*/

	FGP_thickness=ptree.get("free_gas_pocket.thickness",(double)100.); /*m*/
	FGP_width=ptree.get("free_gas_pocket.width",(double)500.); /*m*/
	FGP_offset=ptree.get("free_gas_pocket.depth_below_BSR",(double)0.); /*m*/
	BSR_depth=ptree.get("paleo_conditions.BSR_depth",(double)360.); /*m*/

	fracture_diameter=ptree.get("fracture.diameter",(double)10.); /*m*/
	fracture_depth=ptree.get("fracture.depth",(double)200.); /*m*/

	int numLayers = ptree.get("sediment.number_of_layers",(int)1);
	/*z -> z-location of the bottom of a layer*/
	z = std::vector<double> (numLayers,0.);
	for(int n_layer=0; n_layer<numLayers; n_layer++ ){
		std::string name = "sediment.layer"+std::to_string(n_layer);
		z[n_layer] = ptree.get(name+".z",(double)0.); /*m*/
	}

	prop = std::vector<std::vector<double> > (numLayers,std::vector<double>(8, 0.));
	for(int n_layer=0; n_layer<numLayers; n_layer++ ){
		std::string name = "sediment.layer"+std::to_string(n_layer);
		prop[n_layer][0] = ptree.get(name+".porosity",	(double)0.5);
		prop[n_layer][1] = ptree.get(name+".abs_permeability",(double)1.e-15);/*m^2*/
		prop[n_layer][2] = ptree.get(name+".entry_pressure",(double)1.e5); /*Pa*/
		prop[n_layer][3] = ptree.get(name+".lambda",(double)1.2);
		prop[n_layer][4] = ptree.get(name+".swr",	(double)0.);
		prop[n_layer][5] = ptree.get(name+".sgr",	(double)0.);
		prop[n_layer][6] = ptree.get(name+".m",		(double)3.);
		prop[n_layer][7] = ptree.get(name+".beta",	(double)1.);
	}
	prop_initial_problem = std::vector<std::vector<double> > (numLayers,std::vector<double>(8, 0.));
	for(int n_layer=0; n_layer<numLayers; n_layer++ ){
		prop_initial_problem[n_layer][0] = ptree.get("initial_problem.porosity",(double)0.5);
		prop_initial_problem[n_layer][1] = ptree.get("initial_problem.abs_permeability",(double)1.e-15); /*m^2*/
		prop_initial_problem[n_layer][2] = ptree.get("initial_problem.entry_pressure",(double)1.e5); /*Pa*/
		prop_initial_problem[n_layer][3] = ptree.get("initial_problem.lambda",(double)1.2);
		prop_initial_problem[n_layer][4] = ptree.get("initial_problem.swr",	  (double)0.);
		prop_initial_problem[n_layer][5] = ptree.get("initial_problem.sgr",	  (double)0.);
		prop_initial_problem[n_layer][6] = ptree.get("initial_problem.m",	  (double)3.);
		prop_initial_problem[n_layer][7] = ptree.get("initial_problem.beta",  (double)1.);
	}
	prop_fracture = std::vector<double> (8,0.);
	prop_fracture[0] = ptree.get("fracture.porosity",(double)0.75);
	prop_fracture[1] = ptree.get("fracture.abs_permeability",(double)1.e-10); /*m^2*/
	prop_fracture[2] = ptree.get("fracture.entry_pressure",(double)0.); /*Pa*/
	prop_fracture[3] = ptree.get("fracture.lambda",(double)1.2);
	prop_fracture[4] = ptree.get("fracture.swr",(double)0.);
	prop_fracture[5] = ptree.get("fracture.sgr",(double)0.);
	prop_fracture[6] = ptree.get("fracture.m",(double)3.);
	prop_fracture[7] = ptree.get("fracture.beta",(double)1.);

	gravity_flag = ptree.get("gravity.flag",(bool)true);
	gravity_magnitude = ptree.get("gravity.magnitude",(double)9.81);

	burial_rate = ptree.get("sedimentation.burial_rate",(double)0.01); /*cm/year*/
	burial_rate *= 0.01/(364.*24.*60.*60.); /*convert to m/s*/

	kd = ptree.get("hydrate_phase_change.dissociation_rate",(double)1.e-18);/*mol/m².Pa.s*/
	kf = ptree.get("hydrate_phase_change.formation_rate",(double)1.e-18);/*mol/m².Pa.s*/

	time_initial_problem = ptree.get("initial_problem.time_end",(double)5.); /*years*/
	time_initial_problem *= (1000.*364.*24.*60.*60.); /*convert to seconds*/

	kf_initial_problem = ptree.get("initial_problem.hydrate_formation_rate",(double)1.e-18);/*mol/m².Pa.s*/
	kd_initial_problem = ptree.get("initial_problem.hydrate_dissociation_rate",(double)0.);/*mol/m².Pa.s*/

	salinity_initial_problem = ptree.get("initial_problem.salinity",(double)0.025);

	FGP_Sg = ptree.get("free_gas_pocket.gas_saturation",(double)0.);

	PSF_Pw = ptree.get("paleo_conditions.sea_floor_pressure",(double)10); /*MPa*/
	PSF_Pw *= 1.e6; /*convert to Pa*/

	PSF_T  = ptree.get("paleo_conditions.sea_floor_temperature",(double)4.0); /*degC*/
	PSF_T += 273.15; /*convert to K*/

	thermal_gradient = ptree.get("paleo_conditions.regional_temperature_gradient",(double)35.0); /*degC/km*/
	thermal_gradient *= 1./1000.; /*degC/m */

	BSR_z  = ptree.get("paleo_conditions.BSR_depth",(double)360.); /*m*/

  }
  
  	// FRACTURE-ZONE
	bool isFracture( Dune::FieldVector< double, dim > globalPos /*ndim*/ ) const{
		if( 	globalPos[0]< Xmax/X_c.x_c/2.+fracture_diameter/X_c.x_c/2.
			and globalPos[0]> Xmax/X_c.x_c/2.-fracture_diameter/X_c.x_c/2.
			and globalPos[1]> Zmax/X_c.x_c/2.+fracture_depth/X_c.x_c ){
			return true;
		}else return false;
	}
			//and globalPos[1]< Zmax/X_c.x_c/2.+FGP_thickness/X_c.x_c/2.){ //+BSR_depth/X_c.x_c+FGP_offset/X_c.x_c

	// INITIAL FREE-GAS-ZONE
	bool isInitialFreeGasPocket( Dune::FieldVector< double, dim > globalPos /*ndim*/ )const{
		if(  globalPos[1]> -0.5-FGP_thickness/X_c.x_c//+BSR_depth/X_c.x_c+FGP_offset/X_c.x_c
			and globalPos[1]< -0.5+FGP_thickness/X_c.x_c){ //+BSR_depth/X_c.x_c+FGP_offset/X_c.x_c
			return true;
		}else return false;
	}

	//SEDIMENT LAYERS AND HYDRAULIC PROPERTIES
	std::vector<double> layer_ztop() const {
		return z;
	}
	std::vector< std::vector<double> > layer_properties() const {
		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		if( time_new<time_initial_problem )
			return prop_initial_problem;
		else return prop;
		/*values returned WITH dimensions*/
	}
	std::vector<double> fracture_layer_properties() const {
		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		if( time_new<time_initial_problem )
			return prop_initial_problem[0];
		else return prop_fracture;
		/*values returned WITH dimensions*/
	}
	
	//HYDRATE REACTION KINETIC CONSTANTS
	double HydrateDissociationRateConstant() const {
		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		if( time_new<time_initial_problem )
			return kd_initial_problem;
		else return kd;
		/*mol/m².Pa.s*/
	}
	double HydrateFormationRateConstant() const {
		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		if( time_new<time_initial_problem )
			return kf_initial_problem;
		else return kf;
		/*mol/m².Pa.s*/
	}
	
	//INITIAL PALEO CONDITIONS
	double PSF_Pressure() const {
		return PSF_Pw; /*Pa*/
	}
	double PSF_Temperature() const {
		return PSF_T; /*K*/
	}
	double RegionalThermalGradient() const {
		return thermal_gradient; /*degC/m*/
	}
	double BSR_Depth() const {
		return BSR_z; /*m*/
	}
	double PaleoSalinity() const {
		return salinity_initial_problem; /**/
	}
	
	//FREE GAS POCKET
	double FGP_GasSaturation() const {
		return FGP_Sg; /*-*/
	}
	
    /* SEDIMENT BURIAL RATE */
	Dune::FieldVector<double,dim>
	SedimentationVelocity ( double x /*m*/) const {

		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		double vs0 = burial_rate; /*m/s*/
		
		Dune::FieldVector<double,dim> vs( 0. );
		vs[ dim - 1 ] = -vs0;
		vs[0] = 0.;

		return vs; /*m/s*/
	}
	
	double SedimentationDepth( double x /*m*/ ) const {

		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/

		double sedimentation_depth = 0.;

        double vs /*m/s*/ = (-1.) * SedimentationVelocity(x)[dim-1];
        
        if(time_new>time_initial_problem){
        	sedimentation_depth /*m*/ = vs * (time_new-time_initial_problem);
        }

		return sedimentation_depth; /*m*/
	}

	// COMPACTION
	double CompactionFunction( Dune::FieldVector<double,dim> globalpos /*ndim,ndim*/ ) const {

		double z = globalpos[dim-1]*X_c.x_c; /*m*/
		double sedimentation_depth /*m*/ = SedimentationDepth(globalpos[0]*X_c.x_c/*m*/);
		double beta = 1./3000.;
		double compaction_factor = std::exp( -beta*(Zmax-z + sedimentation_depth) );

		return compaction_factor;
	}
	
	// SEDIMENT DEFORMTION/FLOW VELOCITY
	Dune::FieldVector<double,dim>
	SedimentVelocity () const {
		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		Dune::FieldVector<double,dim> vs( 0. );
		vs[ dim - 1 ] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}
	
	
	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = gravity_magnitude;
		gravity[ dim - 1 ] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}


	/* REFERENCE STATE VALUES */
	double ReferenceSalinity() const {
		return 0.; /*kg/kg*/
	}
	double ReferenceTemperature() const {
		return 273.15; /*K*/
	}
	double ReferencePressure() const {
		return 1.01e5; /*Pa*/
	}
#ifdef STATEINDEPENDENTPROPERTIES
	// (P,T,sal) reference state for CASE1
	double RefP() const {
		return 20.*ReferencePressure(); /*Pa*/
	}

	double ReferenceTemperature() const {
		return ReferenceTemperature() + 8.0; /*K*/
	}

	double RefSal() const {
		return 0.02; /*kg/kg*/;
	}
#endif

};
