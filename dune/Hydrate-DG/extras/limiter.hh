#ifndef DUNE_PDELAB_LIMITER_HH
#define DUNE_PDELAB_LIMITER_HH

#include <limits>

namespace Dune
{
  namespace PDELab
  {

    template <typename GFS, class X>
    class Limiter
    {
    private:
      // extract useful types

      typedef typename GFS::Traits::GridViewType GV;
      enum
      {
        dim = GV::dimension
      };
      typedef typename GV::IndexSet IndexSet;
      typedef typename GV::Grid::ctype DF;
      typedef typename GV::Traits::template Codim<0>::Entity Element;
      typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;
      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename X::template LocalView<LFSCache> VectorView;
      typedef typename LFS::Traits::SizeType size_type;
      typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
      typedef typename LBTraits::RangeFieldType RF;
      typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;

      X *xglobal;

    private:
      /** \brief A class that is to passed to interpolate
       */
      class ConstantInterpolationArgument
      {
      public:
        ConstantInterpolationArgument(RF constant_) : constant(constant_) {}

        template <typename DT, typename RT>
        inline void evaluate(const DT &x, RT &y) const
        {
          y = constant;
          return;
        }

      private:
        RF constant;
      };

      /** \brief A class that is to passed to interpolate
       */
      template <typename Element, typename RF>
      class LinearInterpolationArgument
      {
        enum
        {
          dim = Element::Geometry::coorddimension
        };
        typedef typename Element::Geometry::ctype DF;

      public:
        LinearInterpolationArgument(Element e_, RF coeff_[dim + 1]) : e(e_)
        {
          for (int i = 0; i <= dim; i++)
            coeff[i] = coeff_[i];
          center = e.geometry().center();
        }

        inline RF operator()(const Dune::FieldVector<DF, dim> &x) const
        {
          Dune::FieldVector<DF, dim> global = e.geometry().global(x);
          auto y = coeff[0];
          for (int i = 0; i < dim; i++)
            y += coeff[i + 1] * (global[i] - center[i]);
          return y;
        }

      private:
        Element e;
        RF coeff[dim + 1];
        Dune::FieldVector<DF, dim> center;
      };

      template <typename RF>
      RF minmod(RF a, RF b)
      {
        if (a >= 0.0 && b >= 0.0)
          return theta * std::min(a, b);
        if (a <= 0.0 && b <= 0.0)
          return theta * std::max(a, b);
        return 0.0;
      }

      template <typename RF>
      RF minmod(RF a, RF b, RF c)
      {
        return minmod(a, minmod(b, c));
      }

      const GFS &gfs;
      LFS lfs;
      LFS lfsnb;
      LFSCache lfs_cache;
      LFSCache lfsnb_cache;
      RF theta;

    public:
      Limiter(const GFS &gfs_, X *xglobal_, RF theta_ = 1.0) : gfs(gfs_),
                                                               xglobal(xglobal_),
                                                               lfs(gfs),
                                                               lfsnb(gfs),
                                                               lfs_cache(lfs),
                                                               lfsnb_cache(lfsnb),
                                                               theta(theta_) {}

      //! factor for the timestep in order to fulfill the CFL condition
      inline RF timestepFactor() const
      {
        return 1. / 2.; /// static_cast<RF>(dim);
      }

      void prestage(X &x)
      {
      }

      void poststage(X &x)
      {

        typename X::template LocalView<LFSCache> x_view(x);
        // loop over all interior cells
        for (const auto &cell : elements(gfs.gridView(), Dune::Partitions::interior))
        {
          // extraction of linear part
          RF c[dim + 1];
          for (int i = 0; i <= dim; i++)
            c[i] = 0.0;
          RF b[dim + 1];
          for (int i = 0; i <= dim; i++)
            b[i] = 0.0;

          lfs.bind(cell);
          lfs_cache.update();
          std::vector<RF> xl(lfs.size());
          x_view.bind(lfs_cache);
          x_view.read(xl);
          x_view.unbind();

          const int order = lfs.finiteElement().localBasis().order(); // order of dg
          Dune::GeometryType gt = cell.geometry().type();
          Dune::FieldVector<DF, dim> localcenter = Dune::ReferenceElements<DF, dim>::general(gt).position(0, 0);
          Dune::FieldVector<DF, dim> center = cell.geometry().center();
          const Dune::QuadratureRule<DF, dim> &rule = Dune::QuadratureRules<DF, dim>::rule(gt, 2 * order + 2);
          for (const auto &qip : rule)
          {
            // global position of quadrature point
            Dune::FieldVector<DF, dim> global = cell.geometry().global(qip.position());

            // evaluate basis functions
            std::vector<RangeType> phi(lfs.size());
            lfs.finiteElement().localBasis().evaluateFunction(qip.position(), phi);

            // evaluate u
            RF u = 0.0;
            for (size_type i = 0; i < lfs.size(); i++)
              u += xl[i] * phi[i];

            // integrate
            RF factor = qip.weight() * cell.geometry().integrationElement(qip.position());
            b[0] += u * factor;
            for (int i = 1; i <= dim; i++)
            {
              b[i] += u * (global[i - 1] - center[i - 1]) * factor;
              c[i] += (global[i - 1] - center[i - 1]) * (global[i - 1] - center[i - 1]) * factor;
            }
          }
          c[0] = b[0] / cell.geometry().volume();
          for (int i = 1; i <= dim; i++)
            c[i] = b[i] / c[i];

          // loop over cell neighbors and compute cell averages
          RF avgup[dim];
          RF avgdn[dim];
          for (int i = 0; i < dim; i++)
          {
            avgup[i] = 0;
            avgdn[i] = 0;
          }
          for (const auto &intersection : intersections(gfs.gridView(), cell))
          {
            if (intersection.neighbor())
            {
              // compute cell averages in neighbors
              lfsnb.bind(intersection.outside());
              lfsnb_cache.update();
              std::vector<RF> xlnb(lfsnb.size());
              x_view.bind(lfsnb_cache);
              x_view.read(xlnb);
              x_view.unbind();

              const Dune::QuadratureRule<DF, dim> &rulenb = Dune::QuadratureRules<DF, dim>::rule(gt, 2 * order + 2);
              RF sum = 0.0;
              for (typename Dune::QuadratureRule<DF, dim>::const_iterator qit = rulenb.begin(); qit != rulenb.end(); ++qit)
              {
                // evaluate basis functions
                std::vector<RangeType> phi(lfsnb.size());
                lfsnb.finiteElement().localBasis().evaluateFunction(qit->position(), phi);

                // evaluate u
                RF u = 0.0;
                for (size_type i = 0; i < lfsnb.size(); i++)
                  u += xlnb[i] * phi[i];

                // integrate
                RF factor = qit->weight() * cell.geometry().integrationElement(qit->position());
                sum += u * factor;
              }
              sum /= intersection.outside().geometry().volume();

              // find coordinate direction
              const Dune::FieldVector<DF, dim> n_F = intersection.centerUnitOuterNormal();
              for (int i = 0; i < dim; i++)
              {
                if (n_F[i] > 0.5)
                  avgup[i] = sum;
                if (n_F[i] < -0.5)
                  avgdn[i] = sum;
              }
            }
            if (intersection.boundary())
            {
              const Dune::FieldVector<DF, dim> n_F = intersection.centerUnitOuterNormal();
              for (int i = 0; i < dim; i++)
              {
                if (n_F[i] > 0.5)
                  avgup[i] = c[0];
                if (n_F[i] < -0.5)
                  avgdn[i] = c[0];
              }
            }
          }

          // do the limiting per direction
          Dune::FieldMatrix<DF, dim, dim> jac;
          jac = cell.geometry().jacobianTransposed(localcenter);
          bool limityes = false;
          for (int i = 0; i < dim; i++)
          {
            RF cnew = minmod(c[i + 1], (avgup[i] - c[0]) / jac[i][i], (c[0] - avgdn[i]) / jac[i][i]);
            if (std::abs(cnew - c[i + 1]) > 1.e-12)
            {
              limityes = true;
              c[i + 1] = cnew;
            }
          }

          // set the function on this cell
          if (limityes)
          {
            LinearInterpolationArgument<Element, RF> f(cell, c);
            std::vector<RF> coeff(lfs.size());
            lfs.finiteElement().localInterpolation().interpolate(f, coeff);
            lfs.bind(cell);
            lfs_cache.update();
            x_view.bind(lfs_cache);
            x_view.write(coeff);
            x_view.unbind();
          }
        }

        // communicate limited function
        // Dune::PDELab::CopyDataHandle<GFS,X> dh(gfs,x);
        // if (gfs.gridView().comm().size()>1)
        //   gfs.gridView().communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
      }
    };

    //! second order flux reconstruction with minmod slope limiter in cartesian grids
    template <typename GFS>
    class LimiterFV
    {
    private:
      typedef typename GFS::Traits::GridViewType GV;
      enum
      {
        dim = GV::dimension
      };
      typedef typename GV::IndexSet IndexSet;
      typedef typename GV::Grid::ctype DF;
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;
      typedef typename GV::Traits::template Codim<0>::template Partition<Dune::Interior_Partition>::Iterator ElementIterator;
      typedef typename GV::IntersectionIterator IntersectionIterator;
      typedef typename IntersectionIterator::Intersection Intersection;
      typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
      typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
      typedef typename LFS::Traits::SizeType size_type;
      typedef typename GFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
      typedef typename LBTraits::RangeFieldType RangeFieldType;
      using X = Dune::PDELab::Backend::Vector<GFS, RangeFieldType>;
      typedef typename LBTraits::RangeType RangeType;
      typedef typename LBTraits::DomainType DomainType;
      typedef Dune::PDELab::MonomLocalFiniteElementMap<DF, RangeFieldType, dim, 1> rFEM;
      typedef Dune::PDELab::GridFunctionSpace<GV, rFEM> rGFS;
      using rVEC = Dune::PDELab::Backend::Vector<rGFS, RangeFieldType>;
      typedef Dune::PDELab::LocalFunctionSpace<rGFS> rLFS;
      typedef Dune::PDELab::LFSIndexCache<rLFS> rLFSCache;
      typedef Dune::PDELab::DiscreteGridFunction<rGFS, rVEC> rDGF;

      static GeometryType getCube()
      {
        return GeometryType(GeometryTypes::cube(dim));
      }

    public:
      LimiterFV(const GFS &gfs_, RangeFieldType theta_ = 1.0, const std::string model_ = "minmod") : gfs(gfs_),
                                                                                                     lfs(gfs),
                                                                                                     lfsnb(gfs),
                                                                                                     lfs_cache(lfs),
                                                                                                     lfsnb_cache(lfsnb),
                                                                                                     rfem(getCube()),
                                                                                                     rgfs(gfs.gridView(), rfem),
                                                                                                     rlfs(rgfs),
                                                                                                     rlfs_cache(rlfs),
                                                                                                     rvec(rgfs),
                                                                                                     rvec_view(rvec),
                                                                                                     rdgf(rgfs, rvec),
                                                                                                     theta(theta_),
                                                                                                     model(model_)
      {

        std::cout << "second order reconstruction with theta " << theta << std::endl;
      }

      //! set time of current time step
      void setTime(RangeFieldType t)
      {
        time2 = t;
      }

      //! update the construction values for a given timestep
      inline void updateReconstruction(RangeFieldType time, int stage)
      {
      }

      /** \brief A class that is to passed to interpolate
       */
      template <typename ElementIterator, typename RF>
      class LinearInterpolationArgument
      {
        enum
        {
          dim = ElementIterator::Entity::Geometry::coorddimension
        };
        typedef typename ElementIterator::Entity::Geometry::ctype DF;

      public:
        LinearInterpolationArgument(ElementIterator it_, RF coeff_[dim + 1]) : it(it_)
        {
          for (int i = 0; i <= dim; i++)
            coeff[i] = coeff_[i];
          center = it->geometry().center();
        }

        template <typename DT, typename RT>
        inline void evaluate(const DT &x, RT &y) const
        {
          Dune::FieldVector<DF, dim> global = it->geometry().global(x);
          y = coeff[0];
          for (int i = 0; i < dim; i++)
            y += coeff[i + 1] * (global[i] - center[i]);
          return;
        }

      private:
        ElementIterator it;
        RF coeff[dim + 1];
        Dune::FieldVector<DF, dim> center;
      };

      //! update the construction values for a given timestep
      void prestage(X &x)
      {

        // store slopes
        Dune::FieldVector<Dune::FieldVector<RangeFieldType, 1>, 2 * dim + 1> values(0.0);
        Dune::FieldVector<RangeFieldType, dim> slopes(0.0);

        typename X::template LocalView<LFSCache> x_view(x);

        // for each cell
        ElementIterator it = gfs.gridView().template begin<0, Dune::Interior_Partition>();
        ElementIterator endit = gfs.gridView().template end<0, Dune::Interior_Partition>();
        for (; it != endit; ++it)
        {
          lfs.bind(*it);
          lfs_cache.update();
          std::vector<RangeFieldType> xl(lfs.size());
          x_view.bind(lfs_cache);
          x_view.read(xl);
          x_view.unbind();

          /*
            as we have structured grid, the direction and the axis are computed as follows:

            axis = face/2
            dir = face%2

            where
            face \in [0, 2*dim[
            axis \in [0, dim[
            dir \in {0, 1}.

            the cell value is stored in values[0]
            the neighboring values are stored as values[1+axis*2+dir]
          */

          // compute cell value
          values = std::numeric_limits<RangeFieldType>::quiet_NaN();
          // evaluate basis functions

          Dune::GeometryType gt = it->geometry().type();
          Dune::FieldVector<DF, dim> localcenter = Dune::ReferenceElements<DF, dim>::general(gt).position(0, 0);

          std::vector<RangeType> phi(lfs.size());
          lfs.finiteElement().localBasis().evaluateFunction(it->geometry().center(), phi);
          // evaluate u
          RangeFieldType u = 0.0;
          for (size_type i = 0; i < lfs.size(); i++)
            u += xl[i] * phi[i];
          values[0] = u;
          // compute neighbor values
          typedef typename GV::IntersectionIterator IntersectionIterator;
          IntersectionIterator endit = gfs.gridView().iend(*it);
          IntersectionIterator iit = gfs.gridView().ibegin(*it);
          for (; iit != endit; ++iit)
          {

            int face = iit->indexInInside();
            int axis = face / 2;
            int dir = face % 2;

            if (iit->neighbor())
            {
              lfsnb.bind(*(iit->outside()));
              lfsnb_cache.update();
              std::vector<RangeFieldType> xlnb(lfsnb.size());
              x_view.bind(lfsnb_cache);
              x_view.read(xlnb);
              x_view.unbind();

              std::vector<RangeType> phi(lfsnb.size());
              lfsnb.finiteElement().localBasis().evaluateFunction(it->geometry().center(), phi);

              // evaluate u
              RangeFieldType u = 0.0;
              for (size_type i = 0; i < lfsnb.size(); i++)
                u += xlnb[i] * phi[i];
              values[1 + axis * 2 + dir] = u;
            }
            else
              values[1 + axis * 2 + dir] = values[0];
          }

          // compute slopes for each neighbor
          // and apply slope limiter
          Dune::FieldMatrix<DF, dim, dim> jac;
          jac = it->geometry().jacobianTransposed(localcenter);

          for (int a = 0; a < dim; a++)
          {
            if (model == "minmod")
            {
              slopes[a] = minmod(theta * (values[0] - values[1 + 2 * a]) / jac[a][a], theta * (values[1 + 2 * a + 1] - values[0]) / jac[a][a]);
              // std::cout << "slope " << slopes[a] << " " << values[1+2*a+1] << " " << values[1+2*a] << std::endl;
            }
            else if (model == "minmodchanged")
            {
              slopes[a] = 1. / jac[a][a] * minmod(theta * (values[0] - values[1 + 2 * a]), (values[1 + 2 * a + 1] - values[1 + 2 * a]) / 2.0, theta * (values[1 + 2 * a + 1] - values[0]));
              // std::cout << "slope minch" << slopes[a] << " " << jac[a][a] << std::endl;
            }
            else if (model == "superbee")
            {
              RangeFieldType s1 = minmod(2 * (values[0] - values[1 + 2 * a]), values[1 + 2 * a + 1] - values[0]);
              RangeFieldType s2 = minmod((values[0] - values[1 + 2 * a]), 2 * (values[1 + 2 * a + 1] - values[0]));
              slopes[a] = 1. / jac[a][a] * maxmod(s1, s2);
              //  std::cout << "slope sb" << slopes[a] << " " << values[1+2*a+1] << " " << values[1+2*a] << std::endl;
            }
            else
              DUNE_THROW(Dune::Exception, "limiter model " << model << " is not known.");
          }

          // project onto p1-dg
          rlfs.bind(*it);
          rlfs_cache.update();
          assert(rlfs.size() == dim + 1);
          RangeFieldType xh[dim + 1];
          xh[0] = values[0];
          for (int a = 0; a < dim; a++)
          {
            xh[a + 1] = slopes[a];
            // here I am not sure! should be ok because of interpolation argument!
            // xh[0] -= 0.5 * slopes[a];
          }

          LinearInterpolationArgument<ElementIterator, RangeFieldType> f(it, xh);
          std::vector<RangeFieldType> coeff(rlfs.size());
          rlfs.finiteElement().localInterpolation().interpolate(f, coeff);
          rvec_view.bind(rlfs_cache);
          rvec_view.write(coeff);
          rvec_view.unbind();
        }

        // communicate limited function
        Dune::PDELab::CopyDataHandle<rGFS, rVEC> dh(rgfs, rvec);
        if (rgfs.gridView().comm().size() > 1)
          rgfs.gridView().communicate(dh, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
      }

      //! update the construction values for a given timestep
      void poststage(X &x)
      {
      }

      /*! \brief evaulate the reconstructed flux with in a certain element

        \param e entity to evaluate in
        \param x local position in the entity
        \param c center value in this cell
        \param y result
      */
      inline void evaluate(const ElementType &e,
                           const DomainType &x,
                           const RangeFieldType &c,
                           RangeFieldType &y) const
      {
        FieldVector<RangeFieldType, 1> vn;
        rdgf.evaluate(e, x, vn);
        y = vn[0];
      }

      inline void evaluate(const ElementType &e,
                           const DomainType &x,
                           RangeFieldType &y) const
      {
        // default is the center value
        FieldVector<RangeFieldType, 1> vn;
        rdgf.evaluate(e, x, vn);
        y = vn[0];
      }

      //! factor for the timestep in order to fulfill the CFL condition
      inline RangeFieldType timestepFactor() const
      {
        // need revision
        if (theta > 0)
          return 1.0 / 2.0; // / Traits::dimDomain;
        else
          return 1.0;
      }

    private:
      // minmod slope limiter
      RangeFieldType minmod(RangeFieldType a, RangeFieldType b) const
      {
        if (a * b < 0.0)
          return 0.0;
        if (a > 0.0)
          return std::min(a, b);
        if (a < 0.0)
          return std::max(a, b);
        return 0;
      }

      // minmod slope limiter
      RangeFieldType maxmod(RangeFieldType a, RangeFieldType b) const
      {
        if (std::abs(a) > std::abs(b))
          return a;
        else
          return b;
      }

      // minmod slope limiter
      RangeFieldType minmod(RangeFieldType a, RangeFieldType b, RangeFieldType c) const
      {
        return minmod(a, minmod(b, c));
      }

      // function which is to be reconstructed
      const GFS &gfs;
      LFS lfs;
      LFS lfsnb;
      LFSCache lfs_cache;
      LFSCache lfsnb_cache;
      // functionspace related variables
      rFEM rfem;
      rGFS rgfs;
      rLFS rlfs;
      rLFSCache rlfs_cache;
      // data vector and function
      rVEC rvec;
      typename rVEC::template LocalView<rLFSCache> rvec_view;
      rDGF rdgf;

      // just needed for assertions
      RangeFieldType time1;
      RangeFieldType time2;

      // limiter parameters
      RangeFieldType theta;
      std::string model;
    };

  } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_LIMITER_HH
