/**
  This example shows how to use the Non-linear solvers from the
  'unsupported' section of Eigen.

  NOTE: The NonLinearOptimization Lib REQUIRES the use of Dense Matrices
*/

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <unsupported/Eigen/NonLinearOptimization>
#include <time.h>

// Functor (template of a function to be used by Non-Linear Solvers)
// --------------------------------------------------------------------------------
template<typename Scalar_, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct Functor
{
  typedef Scalar_ Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

  // REQUIRED override
  virtual int operator()(const InputType &x, ValueType& fvec) const;

  // Only needed if you want to formulate own Jacobian
  int df(const InputType &x, JacobianType& fjac) const;
};
//================================================================================



// ============== Functor with hand-coded Jacobian =========================
struct hybrj_functor : Functor<double>
{
    hybrj_functor(void) : Functor<double>(9,9) {}

    // Evaluation Function (F(x) = ...)
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        double temp, temp1, temp2;
        const Eigen::VectorXd::Index n = x.size();
        assert(fvec.size()==n);
        for (Eigen::VectorXd::Index k = 0; k < n; k++)
        {
            temp = (3. - 2.*x[k])*x[k];
            temp1 = 0.;
            if (k) temp1 = x[k-1];
            temp2 = 0.;
            if (k != n-1) temp2 = x[k+1];
            fvec[k] = temp - temp1 - 2.*temp2 + 1.;
        }
        return 0;
    }

    // Jacobian Function (dF/dx = ...)
    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
    {
        const Eigen::VectorXd::Index n = x.size();
        assert(fjac.rows()==n);
        assert(fjac.cols()==n);
        for (Eigen::VectorXd::Index k = 0; k < n; k++)
        {
            for (Eigen::VectorXd::Index j = 0; j < n; j++)
                fjac(k,j) = 0.;
            fjac(k,k) = 3.- 4.*x[k];
            if (k) fjac(k,k-1) = -1.;
            if (k != n-1) fjac(k,k+1) = -2.;
        }
        return 0;
    }
};
// =============================================================================

// ================== Function to test the above solve ===========================
int testHybrj()
{
  const int n=9;
  int info;
  Eigen::VectorXd x(n);

  /* the following starting values provide a rough fit. */
  x.setConstant(n, -1.);

  // do the computation
  hybrj_functor functor;
  Eigen::HybridNonLinearSolver<hybrj_functor> solver(functor);
  info = solver.solve(x);
  EIGEN_UNUSED_VARIABLE(info)

  // check return value
  // VERIFY_IS_EQUAL(info, 1);

  // check norm
  //VERIFY_IS_APPROX(solver.fvec.blueNorm(), 1.192636e-08);
  std::cout << "fnorm="<< solver.fvec.blueNorm() << std::endl;


// check x
  Eigen::VectorXd x_ref(n);
  x_ref <<
     -0.5706545,    -0.6816283,    -0.7017325,
     -0.7042129,     -0.701369,    -0.6918656,
     -0.665792,    -0.5960342,    -0.4164121;
  //VERIFY_IS_APPROX(x, x_ref);
  std::cout << "x=" << std::endl;
  std::cout << x << std::endl << std::endl;
  std::cout << "x_ref=" << std::endl;
  std::cout << x_ref << std::endl;
  std::cout << "status=";
  std::cout << info << std::endl;

  return info;
}
//==============================================================================


// ============== Functor WITHOUT hand-coded Jacobian =========================
struct hybrd_functor : Functor<double>
{
    hybrd_functor(void) : Functor<double>(9,9) {}
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        double temp, temp1, temp2;
        const Eigen::VectorXd::Index n = x.size();

        assert(fvec.size()==n);
        for (Eigen::VectorXd::Index k=0; k < n; k++)
        {
            temp = (3. - 2.*x[k])*x[k];
            temp1 = 0.;
            if (k) temp1 = x[k-1];
            temp2 = 0.;
            if (k != n-1) temp2 = x[k+1];
            fvec[k] = temp - temp1 - 2.*temp2 + 1.;
        }
        return 0;
    }
};
//=============================================================================

// ================== Function to test the above solve ===========================
int testHybrd()
{
  const int n=9;
  int info;
  Eigen::VectorXd x;

  /* the following starting values provide a rough fit. */
  x.setConstant(n, -1.);

  // do the computation
  hybrd_functor functor;
  Eigen::HybridNonLinearSolver<hybrd_functor> solver(functor);
  solver.parameters.nb_of_subdiagonals = 1;
  solver.parameters.nb_of_superdiagonals = 1;
  solver.diag.setConstant(n, 1.);
  solver.useExternalScaling = true;
  info = solver.solveNumericalDiff(x); //Call this instead of '.solve(x)'
  EIGEN_UNUSED_VARIABLE(info)

  // check return value
  // VERIFY_IS_EQUAL(info, 1);
  //VERIFY(solver.nfev <= 14*LM_EVAL_COUNT_TOL);

  // check norm
  //VERIFY_IS_APPROX(solver.fvec.blueNorm(), 1.192636e-08);
  std::cout << "fnorm="<< solver.fvec.blueNorm() << std::endl;

  // check x
  Eigen::VectorXd x_ref(n);
  x_ref <<
      -0.5706545,    -0.6816283,    -0.7017325,
      -0.7042129,     -0.701369,    -0.6918656,
      -0.665792,    -0.5960342,    -0.4164121;
  //VERIFY_IS_APPROX(x, x_ref);
  std::cout << "x=" << std::endl;
  std::cout << x << std::endl << std::endl;
  std::cout << "x_ref=" << std::endl;
  std::cout << x_ref << std::endl;
  std::cout << "status=";
  std::cout << info << std::endl;

  return info;
}
//=============================================================================



int main()
{
  int success;

  success = testHybrj();
  if (success != 1)
  {
    std::cout << "ERROR in Jacobian Solve for 'doleg' method\n\n";
    return -1;
  }

  success = testHybrd();
  if (success != 1)
  {
    std::cout << "ERROR in Numerical Jacobian Solve for 'doleg' method\n\n";
    return -1;
  }

  std::cout << "Finished!\n\n";
  return 0;
}
