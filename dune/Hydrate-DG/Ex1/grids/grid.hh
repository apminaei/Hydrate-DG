template <typename PTree>
class MeshParameters
{
private:
	const PTree &ptree;

	constexpr static double eps = 1.0e-6;
	CharacteristicValues Xc;

	double Zmin;
	double Xmax;
	double dZ;
	double dX;

public:
	//! constructor
	MeshParameters(const PTree &ptree_)
		: ptree(ptree_)
	{
		Zmin = ptree.get("grid.LZ", (double)1.) / Xc.x_c; // ndim
		dZ = ptree.get("grid.NZ", (int)10);
		Xmax = ptree.get("grid.LX", (double)1.) / Xc.x_c; // ndim
		dX = ptree.get("grid.NX", (int)10);
	}

	/*
	 * 2D -> X and Z
	 *					X_min,Z_max = 0,0 	----------------------------------------- X_max,Z_max = X_length,0
											|				|	|					|
											|				|	|					|
											|				-----					|
											|---------------------------------------|
											|---------------------------------------|
	 * 										|										|
	 * 										|										|
	 * 			X_min,Z_min = 0,Z_length	----------------------------------------- X_max,Z_min = X_length,-Z_length
	 *
	 *
	 *
	 */

	const static int dimension = 1; // ptree.get("grid.worlddim",(int)2);
	constexpr static double origin = 0.;

	const double Z_length = ptree.get("grid.LZ", (double)1.) / Xc.x_c; // m
	const int Z_cells = ptree.get("grid.NZ", (int)10);
	const double X_length = ptree.get("grid.LX", (double)1.) / Xc.x_c; // m
	const int X_cells = ptree.get("grid.NX", (int)10);
	const double Y_length = 1.; // only if dim=3
	const int Y_cells = 1.;		// only for dim==3

	const double Z_GHSZ_bottom = ptree.get("grid.ghsz.Zbottom", (double)1.) / Xc.x_c; // ndim
	const double Z_GHSZ_top = ptree.get("grid.ghsz.Ztop", (double)1.) / Xc.x_c;		  // ndim
	const double X_GHSZ_left = ptree.get("grid.ghsz.Xleft", (double)1.) / Xc.x_c;	  // ndim
	const double X_GHSZ_right = ptree.get("grid.ghsz.Xright", (double)1.) / Xc.x_c;	  // ndim

	const double Zmin_lenz = ptree.get("grid.lenz.Zmin", (double)1.) / Xc.x_c; // m
	const double Zmax_lenz = ptree.get("grid.lenz.Zmax", (double)1.) / Xc.x_c; // m
	const double Xmin_lenz = ptree.get("grid.lenz.Xmin", (double)1.) / Xc.x_c; // m
	const double Xmax_lenz = ptree.get("grid.lenz.Xmax", (double)1.) / Xc.x_c; // m
	const double X0_lenz = ptree.get("grid.lenz.X0", (double)1.) / Xc.x_c;	   // m
	const double X1_lenz = ptree.get("grid.lenz.X1", (double)1.) / Xc.x_c;	   // m

	double volumeFactor(double r /*radial component of cell center*/) const
	{

		/* Volume = 2 PI r_center dr dz
		 * */
		double factor = 1.; //
		return factor;		/*m^3*/
	}

	double areaFactor(double r_c /*radial component of cell center*/,
					  double r_f /*radial component of face center*/,
					  Dune::FieldVector<double, dimension> normal) const
	{

		/* Surface area = 2 PI ( r_center dr i_hat + r_face dz j_hat )
		 * */
		double factor = 1.;
		return abs(factor); /*m^2*/
	}

	bool isLeftBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{

		if (dimension == 2)
		{
			if (globalPos[0] == origin)
			{
				return true;
			}
			else
				return false;
		}
		else if (dimension == 3)
		{
			if (globalPos[0] < origin + eps)
			{
				return true;
			}
			else
				return false;
		}
	}

	bool isRightBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{

		if (dimension == 2)
		{
			if (globalPos[0] == X_length)
			{
				return true;
			}
			else
				return false;
		}
		else if (dimension == 3)
		{
			if (globalPos[0] > X_length - eps)
			{
				return true;
			}
			else
				return false;
		}
	}

	bool isBottomBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{
		if (dimension == 1)
		{
			if (globalPos[0] == Z_length)
			{
				return true;
			}
			else
				return false;
		}
		else if (dimension == 2)
		{
			if (globalPos[1] == Z_length && globalPos[0] != origin && globalPos[0] != X_length)
			{ //+eps
				return true;
			}
			else
				return false;
		}
		else if (dimension == 3)
		{
			if (globalPos[2] < Z_length + eps)
			{
				return true;
			}
			else
				return false;
		}
	}

	bool isTopBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{
		if (dimension == 1)
		{
			if (globalPos[0] == origin)
			{ //-eps
				return true;
			}
			else
				return false;
		}
		else if (dimension == 2)
		{
			if (globalPos[1] == origin && globalPos[0] != origin && globalPos[0] != X_length)
			{ //-eps
				return true;
			}
			else
				return false;
		}
		else if (dimension == 3)
		{
			if (globalPos[2] > origin - eps)
			{
				return true;
			}
			else
				return false;
		}
	}

	bool isFrontBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{
		if (dimension == 3)
		{
			if (globalPos[1] < origin + eps)
			{
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}

	bool isBackBoundary(Dune::FieldVector<double, dimension> globalPos) const
	{
		if (dimension == 3)
		{
			if (globalPos[1] > Y_length - eps)
			{
				return true;
			}
			else
				return false;
		}
		else
			return false;
	}

	//##################################################################################
	bool isGHSZ(Dune::FieldVector<double, dimension> globalPos) const
	{

		if (dimension == 1)
		{
			if ((Z_GHSZ_bottom) <= globalPos[0] && globalPos[0] <= (Z_GHSZ_top))
				return true;
			else
				return false;
		}
		else if (dimension == 2)
		{
			if (Z_GHSZ_bottom <= globalPos[1] && globalPos[1] <= Z_GHSZ_top)
				return true;
			else
				return false;
		}
		else if (dimension == 3)
		{
			if ((Z_GHSZ_bottom - eps) < globalPos[2] && globalPos[2] < (Z_GHSZ_top + eps))
				return true;
			else
				return false;
		}
	}
	//
	//##################################################################################
	bool isLenz1(Dune::FieldVector<double, dimension> globalPos) const
	{
		if (
			(150 / 150 * (globalPos[0] - 30 / Xc.x_c) + globalPos[1]) >= (-200 / Xc.x_c) && (150 / 150 * (globalPos[0] - 40 / Xc.x_c) + globalPos[1]) <= (-190 / Xc.x_c) && (-15 / 15 * (globalPos[0] - 30 / Xc.x_c) + globalPos[1]) <= (-200 / Xc.x_c) && (-15 / 15 * (globalPos[0] - 170 / Xc.x_c) + globalPos[1]) >= (-340 / Xc.x_c))
		{
			return true;
		}
		else
			return false;
	}

	//##################################################################################
	bool isWell(Dune::FieldVector<double, dimension> globalPos) const
	{
		if ((globalPos[0] == origin) && (globalPos[1] > Z_GHSZ_bottom - eps))
		{
			return true;
		}
		else
			return false;
	}

	int getdim() const
	{
		return dimension;
	}
};
