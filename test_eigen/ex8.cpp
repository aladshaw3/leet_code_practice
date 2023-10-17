/**
  This example shows how to use the LM solvers from the
  'unsupported' section of Eigen.

  The LevenbergMarquardt method is usually for non-square
  optimization problems (such as curve-fitting), but can
  also be used for square problems.

  Like with the NonLinearOptimization methods, user must create
  an instance of a 'Functor' struct. However, this time the user has 2 existing structs:
  'DenseFunctor' and 'SparseFunctor' to start from. This object also allows the user
  to create and utilize Sparse Matrices coded by hand in order to solve very large
  problems. This method appears to be much more polished than the other in ex7.
*/

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/LevenbergMarquardt>
#include <time.h>

struct lmdif_drx_functor : Eigen::DenseFunctor<double>
{
    //Added new data for dynamically tracking model size
    unsigned int vars, funs;

    lmdif_drx_functor(void) : Eigen::DenseFunctor<double>() {}
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        double dx = 1.0/(double)funs;
        double k = 0.5;
        double D = 0.001;
        // D*d^2x/dz^2 - k*x^2
        for (unsigned int i=1; i<funs-1; i++)
        {
            fvec[i] = D*(x[i-1] - 2.*x[i] + x[i+1])/dx/dx - k*x[i]*x[i];
        }
        // BCs
        fvec[0] = x[0] - 1.;             //Dirchlet
        fvec[funs-1] = x[funs-1] - x[funs-2];  //1st Order Neumann
        return 0;
    }

    //Added new functions for setting size
    void setVars(int v) {vars = v;}
    void setFuns(int f) {funs = f;}

    // REQUIRED: Override 'inputs()' and 'values()' to get correct outputs
    int inputs() const { return vars; }
    int values() const { return funs; }
};


int testLmdif_drx(int n)
{
  int info;
  double fnorm;
  Eigen::VectorXd x(n);

  /* the following starting values provide a rough fit. */
  x.setConstant(n, 1.);

  // do the computation
  lmdif_drx_functor functor;
  functor.setVars(n);
  functor.setFuns(n);
  //std::cout << functor.values() << std::endl;
  //std::cout << functor.inputs() << std::endl;

  Eigen::NumericalDiff<lmdif_drx_functor > numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lmdif_drx_functor > > lm(numDiff);
  lm.setMaxfev(n*100);
  lm.setXtol(1e-6);
  lm.setFtol(1e-6);
  lm.setGtol(1e-6);
  lm.setEpsilon(1e-6);
  lm.setFactor(100);
  info = lm.minimize(x);
  EIGEN_UNUSED_VARIABLE(info)

  // check return values
  std::cout << "status=" << info << std::endl << std::endl;

  // check norm
  fnorm = lm.fnorm();
  std::cout << "fnorm=" << fnorm << std::endl << std::endl;

  // check x
  if (x.size() < 20)
  {
    std::cout << "x=" << std::endl;
    std::cout << x << std::endl << std::endl;
  }
  else
  {
    std::cout << "iter=" << lm.iterations() << std::endl;
    std::cout << "feval=" << lm.nfev() << std::endl;
    std::cout << "jeval=" << lm.njev() << std::endl;
  }

  return info;

}



template<int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct lmdif_drx_functor_sparse : Eigen::SparseFunctor<double, int>
{
    //Added new data for dynamically tracking model size
    unsigned int vars, funs;

    lmdif_drx_functor_sparse(void) : Eigen::SparseFunctor<double, int>(NX, NY)
      { vars=NX; funs=NY;}
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        double dx = 1.0/(double)funs;
        double k = 0.5;
        double D = 0.001;
        // D*d^2x/dz^2 - k*x^2
        for (unsigned int i=1; i<funs-1; i++)
        {
            fvec[i] = D*(x[i-1] - 2.*x[i] + x[i+1])/dx/dx - k*x[i]*x[i];
        }
        // BCs
        fvec[0] = x[0] - 1.;             //Dirchlet
        fvec[funs-1] = x[funs-1] - x[funs-2];  //1st Order Neumann
        return 0;
    }

    // Helper function to initialize the SparseMatrix (call ONLY once)
    int df_init(const Eigen::VectorXd &x, Eigen::SparseMatrix<double, Eigen::ColMajor, int> &fjac) const
    {
        fjac.conservativeResize(funs,funs);
        //Create space for 3 non-zeros per row/col
        fjac.reserve( Eigen::VectorXi::Constant(funs,3) );

        double dx = 1.0/(double)funs;
        double k = 0.5;
        double D = 0.001;
        //Loop over i functions
        for (unsigned int i = 1; i < funs-1; i++)
        {
            fjac.insert(i,i-1) = D*(1.)/dx/dx;
            fjac.insert(i,i) = D*(- 2.)/dx/dx - 2.*k*x[i];
            fjac.insert(i,i+1) = D*(1.)/dx/dx;
        }
        fjac.insert(0,0) = 1.;
        fjac.insert(funs-1,funs-1) = 1.;
        fjac.insert(funs-1,funs-2) = -1.;
        fjac.makeCompressed();
        return 0;
    }

    // This function will assume the sparcity pattern of the Jacobian has
    //  been previously established
    int df(const Eigen::VectorXd &x, Eigen::SparseMatrix<double, Eigen::ColMajor, int> &fjac) const
    {
        double dx = 1.0/(double)funs;
        double k = 0.5;
        double D = 0.001;
        //Loop over i functions
        /*
            NOTE: May be more efficient to access from valuePtr(),
            innerIndexPtr(), and outerIndexPtr().

            see https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html
        */
        for (unsigned int i = 1; i < funs-1; i++)
        {
            fjac.coeffRef(i,i-1) = D*(1.)/dx/dx;
            fjac.coeffRef(i,i) = D*(- 2.)/dx/dx - 2.*k*x[i];
            fjac.coeffRef(i,i+1) = D*(1.)/dx/dx;
        }
        fjac.coeffRef(0,0) = 1.;
        fjac.coeffRef(funs-1,funs-1) = 1.;
        fjac.coeffRef(funs-1,funs-2) = -1.;

        fjac.makeCompressed();
        return 0;
    }


    //Added new functions for setting size
    void setVars(int v) {vars = v;}
    void setFuns(int f) {funs = f;}

    // REQUIRED: Override 'inputs()' and 'values()' to get correct outputs
    int inputs() const { return vars; }
    int values() const { return funs; }

};


int testLmdif_drx_sparse(int n)
{
  int info;
  double fnorm;
  Eigen::VectorXd x(n), scale(n);

  info = 0;

  /* the following starting values provide a rough fit. */
  x.setConstant(n, 1.);

  // Use this to set a scaling factor for each equation
  //    Here, each equation is scaled to 1.
  scale.setConstant(n, 1.);

  // do the computation
  lmdif_drx_functor_sparse<Eigen::Dynamic,Eigen::Dynamic> functor;
  functor.setVars(n);
  functor.setFuns(n);

  Eigen::LevenbergMarquardt<lmdif_drx_functor_sparse<Eigen::Dynamic,Eigen::Dynamic> > lm(functor);

  //Jacobian reference accessible at lm.jacobian()
  //    Use to initialize matrix size and memory
  info = functor.df_init(x,lm.jacobian());

  //Quick test of df function
  /*
  if (n < 20)
  {
    std::cout << lm.jacobian() << std::endl;
    x.setConstant(n, 0.);
    info = functor.df(x,lm.jacobian());
    std::cout << lm.jacobian() << std::endl;
  }
  */

  // For a Sparse Functor, you CANNOT use NumericalDiff
  lm.setMaxfev(n*100);  // Maximum calls to f()
  lm.setXtol(1e-8);     // Cutoff tolerance for change in x
  lm.setFtol(1e-6);     // Cutoff tolerance for norm in f
  lm.setGtol(0.);       // Cutoff tolerance for change in df()
  lm.setEpsilon(0.);    // Error percision
  lm.setFactor(100);    // Step bound for diagonal shift

  // Use these options for user defined scaling
  lm.setExternalScaling(true);
  lm.diag() = scale;

  // Call the solve
  info = lm.minimize(x);
  EIGEN_UNUSED_VARIABLE(info)

  // check return values
  std::cout << "status=" << info << std::endl << std::endl;

  // check norm
  fnorm = lm.fnorm();
  std::cout << "fnorm=" << fnorm << std::endl;
  std::cout << "gnorm=" << lm.gnorm() << std::endl;
  std::cout << "lm_param=" << lm.lm_param() << std::endl << std::endl;

  // check x
  if (x.size() < 20)
  {
    std::cout << "x=" << std::endl;
    std::cout << x << std::endl << std::endl;
  }
  else
  {
    std::cout << "iter=" << lm.iterations() << std::endl;
    std::cout << "feval=" << lm.nfev() << std::endl;
    std::cout << "jeval=" << lm.njev() << std::endl;
  }

  return info;

}



int main()
{
  int success = 0;
  double time;

  // output status is good for
  //  [1-4, 6-8]
  for (int i=10; i<=50; i+=5)
  {
    time = clock();
    int dn = (i)*10;
    success = testLmdif_drx(dn);
    time = clock() - time;
    std::cout << dn << " equ - Dense Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
    if (success != 1 && success != 2 && success != 3 &&
        success != 4 && success != 6 && success != 7 && success != 8)
    {
      std::cout << "ERROR in Numerical Jacobian Solve for 'LM' method\n\n";
      return -1;
    }
    std::cout << "\n\n";

    time = clock();
    int sn = (i)*100;
    success = testLmdif_drx_sparse(sn);
    time = clock() - time;
    std::cout << sn << " equ - Sparse Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;
    if (success != 1 && success != 2 && success != 3 &&
        success != 4 && success != 6 && success != 7 && success != 8)
    {
      std::cout << "ERROR in Sparse Jacobian Solve for 'LM' method\n\n";
      return -1;
    }
  }

  success = 0;
  std::cout << "Finished!\n\n";
  return success;
}
