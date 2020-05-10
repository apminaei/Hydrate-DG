template<typename PTree>
class MeshParameters{
private:
	const PTree& ptree;
	constexpr static double eps = 1.0e-6;
	CharacteristicValues Xc;

	double Zmax;
	double Xmax;
	double dZ;
	double dX;
	double Xminlens;
	double Zminlens;
	double Xmaxlens;
	double Zmaxlens;
	double Zmininlet;
	double Zmaxinlet;

public:
	
	//! constructor
	MeshParameters (const PTree& ptree_)
	:ptree(ptree_)
	{
		Zmax = ptree.get("grid.yasp.LZ",(double)1.)/Xc.x_c; //m
		dZ = ptree.get("grid.yasp.NZ",(int)10);
		Xmax = ptree.get("grid.yasp.LX",(double)1.)/Xc.x_c; //m
		dX = ptree.get("grid.yasp.NX",(int)10)	;
		Xminlens = ptree.get("grid.yasp.lens.Xmin",(double)1.)/Xc.x_c; //m
		Zminlens = ptree.get("grid.yasp.lens.Zmin",(double)1.)/Xc.x_c; //m
		Xmaxlens = ptree.get("grid.yasp.lens.Xmax",(double)1.)/Xc.x_c; //m
		Zmaxlens = ptree.get("grid.yasp.lens.Zmax",(double)1.)/Xc.x_c; //m
		Zmininlet = ptree.get("grid.yasp.inlet.Zmin",(double)1.)/Xc.x_c; //m
		Zmaxinlet = ptree.get("grid.yasp.inlet.Zmax",(double)1.)/Xc.x_c; //m
	}

	/*
	 * 2D -> X and Z
	 */
	const static int dimension = 2;
	const double Z_length = ptree.get("grid.yasp.LZ",(double)1.)/Xc.x_c; //m 
	const int Z_cells = ptree.get("grid.yasp.NZ",(int)10);
	const double X_length = ptree.get("grid.yasp.LX",(double)1.)/Xc.x_c; //m
	const int X_cells = ptree.get("grid.yasp.NX",(int)10);

	bool isRight( Dune::FieldVector<double,dimension> globalPos ) const {
		if( globalPos[0] > Xmax - eps ){
			return true;
		}
		else return false;
	}

	bool isLeft( Dune::FieldVector<double,dimension> globalPos ) const {
		if( globalPos[0] < 0. + eps ){
			return true;
		}
		else return false;
	}

	bool isInlet( Dune::FieldVector<double,dimension> globalPos ) const {
		if( isLeft(globalPos) and (globalPos[1] < Zmaxinlet and globalPos[1] > Zmininlet) ){
			return true;
		}
		else return false;
	}

	bool isLens( Dune::FieldVector<double,dimension> globalPos ) const {
		if(    (globalPos[1] < Zmaxlens and globalPos[1] > Zminlens)
		   and (globalPos[0] < Xmaxlens and globalPos[0] > Xminlens) ){
			return true;
		}
		else return false;
	}
	
};
