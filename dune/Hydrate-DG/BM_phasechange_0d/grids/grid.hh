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

public:
	
	//! constructor
	MeshParameters (const PTree& ptree_)
	:ptree(ptree_)
	{
		Zmax = ptree.get("grid.yasp.LZ",(double)1.)/Xc.x_c; //m
		dZ = ptree.get("grid.yasp.NZ",(int)10);
		Xmax = ptree.get("grid.yasp.LX",(double)1.)/Xc.x_c; //m
		dX = ptree.get("grid.yasp.NX",(int)10)	;
	}

	/*
	 * 2D -> X and Z
	 */

	const static int dimension 			= 2		;
	constexpr static double origin 		= 0.	;

	const double Z_length = ptree.get("grid.yasp.LZ",(double)1.)/Xc.x_c; //m 
	const int Z_cells = ptree.get("grid.yasp.NZ",(int)10);
	const double X_length = ptree.get("grid.yasp.LX",(double)1.)/Xc.x_c; //m
	const int X_cells = ptree.get("grid.yasp.NX",(int)10);
	const double Y_length = 1.;	//only if dim=3
	const int Y_cells = 1.; //only for dim==3

	int getdim( ) const{
		return dimension;
	}
	
};
