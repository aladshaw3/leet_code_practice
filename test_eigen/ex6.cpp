/**
  This example shows how to use the Matrix-Free linear solvers,
  including the 'Unofficially' supported methods like GMRES.
*/

/**
  Iterative solvers such as ConjugateGradient and BiCGSTAB can be used in a matrix
  free context. To this end, user must provide a wrapper class inheriting EigenBase<>
  and implementing the following methods:

    - Index rows() and Index cols(): returns number of rows and columns respectively
    - operator* with your type and an Eigen dense column vector
        (its actual implementation goes in a specialization of the
        internal::generic_product_impl class)

  Eigen::internal::traits<> must also be specialized for the wrapper type.

  Here is a complete example wrapping an Eigen::SparseMatrix:
*/

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <time.h>

class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
  namespace internal {
    // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
    //      NOTE: This is specifically for a Matrix like object of type 'double'
    //            That typdef is specific to this wrapper
    template<>
    struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
    {};
  }
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement>
{
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum
  {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  // Required method to override
  Index rows() const { return mp_mat->rows(); }
  Index cols() const { return mp_mat->cols(); }

  // Required operator to override
  //    NOTE: This is not the actual implementation section
  template<typename Rhs>
  Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const
  {
    return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // Custom API:
  MatrixReplacement() : mp_mat(0) {}

  void attachMyMatrix(const SparseMatrix<double> &mat) {
    mp_mat = &mat;
  }
  const SparseMatrix<double> my_matrix() const { return *mp_mat; }

private:
  const SparseMatrix<double> *mp_mat;
};


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
  namespace internal {

    // This is where the actual implementation of the Matrix Multiply action is done
    template<typename Rhs>
    struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
    : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
    {
      typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;

      template<typename Dest>
      static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
      {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
        assert(alpha==Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);

        // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
        // but let's do something fancier (and less efficient):
        // ------------- Manual multiply-------------
        //      NOTE: This is very, very slow
        /**
        for(Index i=0; i<lhs.cols(); ++i)
        {
          dst += rhs(i) * lhs.my_matrix().col(i);
        }
        */
        // ------------- Quick multiply (because it is just a wrapper to SparseMatrix) -------------
        dst.noalias() += lhs.my_matrix() * rhs;
      }
    };

  }
}


// Helper function to fill in the b vector for BCs of 3D Laplacian
void fillNatural3DLaplacian_BCs(unsigned int m, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& b)
{
  b.resize(m*m*m,1);
  b(0,0) = -1;
  b(m*m-1,0) = -1;
  b(m*m*m-1,0) = -1;
}

// Helper function to fill in a 3D sparse Laplacian
void fillNatural3DLaplacian(unsigned int m, Eigen::SparseMatrix<double>& S)
{
  // Setup the matrix size
  unsigned int N = m*m*m;
  S.conservativeResize(N,N);

  //Create space for 7 non-zeros per row/col
  S.reserve( Eigen::VectorXi::Constant(N,7) );

  //Create piviots for marking when things have changed in a row/col
  unsigned int r = m;
  unsigned int r2 = m*m;
  unsigned int r3 = r2;
  unsigned int r4 = 0;

  //Loop along each row
  for (unsigned int i=0; i<N; i++)
  {
    //If statements for tridiagonal portion
    S.insert(i, i) = 6;
    if (i == 0)
    {
      S.insert(i, i+1) = -1;
    }
    else if (i == N-1)
    {
      S.insert(i, i-1) = -1;
    }
    else if (i == r-1)
    {
      S.insert(i, i-1) = -1;
    }
    else if (i == r)
    {
      S.insert(i, i+1) = -1;
      r = r + m;
    }
    else
    {
      S.insert(i, i-1) = -1;
      S.insert(i, i+1) = -1;
    }

    //If statements for 2nd diagonal bands
    if (i > m-1)
    {
      if (i <= r3-1)
      {
        S.insert(i, i-m) = -1;
      }
      else if (i > r3-1)
      {
        r4 = r4+1;
        if (r4 == m-1)
        {
          r3 = r2;
          r4 = 0;
        }
      }
    }
    if (i <= N-m-1 && i <= r2-m-1)
    {
      S.insert(i, i+m) = -1;
    }
    if (i == r2-1)
    {
      r2 = r2+(m*m);
    }

    //If statements for 3rd diagonal bands
    if (i > (m*m)-1)
    {
      S.insert(i, i-(m*m)) = -1;
    }
    if (i <= N-(m*m)-1)
    {
      S.insert(i, i+(m*m)) = -1;
    }
  }

  // Lastly, compress the representation
  S.makeCompressed();
}




int main()
{
  double time;
  int n = 50;
  Eigen::SparseMatrix<double> S;
  fillNatural3DLaplacian(n, S);

  MatrixReplacement A;
  A.attachMyMatrix(S);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b, x;
  x.resize(n*n*n,1);
  fillNatural3DLaplacian_BCs(n, b);

  std::cout << "Mat Size: " << n*n*n << " x " << n*n*n << std::endl << std::endl;

  // Solve Ax = b using various iterative solver with matrix-free version:
  {
    time = clock();
    Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.setTolerance(1e-6);
    cg.compute(A);
    x = cg.solve(b);
    time = clock() - time;
    std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    std::cout << "MatFree CG Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( cg.error() > 1e-6)
    {
      std::cout << "Failed line 246" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> cg;
    cg.setTolerance(1e-6);
    cg.compute(S);
    x = cg.solve(b);
    time = clock() - time;
    std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
    std::cout << "Std CG Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( cg.error() > 1e-6)
    {
      std::cout << "Failed line 264" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
    bicg.setTolerance(1e-6);
    bicg.compute(A);
    x = bicg.solve(b);
    time = clock() - time;
    std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
    std::cout << "MatFree BiCGSTAB Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( bicg.error() > 1e-6)
    {
      std::cout << "Failed line 282" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> bicg;
    bicg.setTolerance(1e-6);
    bicg.compute(S);
    x = bicg.solve(b);
    time = clock() - time;
    std::cout << "BiCGSTAB: #iterations: " << bicg.iterations() << ", estimated error: " << bicg.error() << std::endl;
    std::cout << "Std BiCGSTAB Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( bicg.error() > 1e-6)
    {
      std::cout << "Failed line 300" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::GMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.setTolerance(1e-6);
    gmres.setMaxIterations(100);
    gmres.set_restart(50);
    gmres.compute(A);
    x = gmres.solve(b);
    time = clock() - time;
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "MatFree GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( gmres.error() > 1e-6)
    {
      std::cout << "Failed line 320" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> gmres;
    gmres.setTolerance(1e-6);
    gmres.setMaxIterations(100);
    gmres.set_restart(50);
    gmres.compute(S);
    x = gmres.solve(b);
    time = clock() - time;
    std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "Std GMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( gmres.error() > 1e-6)
    {
      std::cout << "Failed line 340" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::DGMRES<MatrixReplacement, Eigen::IdentityPreconditioner> gmres;
    gmres.setTolerance(1e-6);
    gmres.setMaxIterations(1000);
    gmres.set_restart(50);
    gmres.compute(A);
    x = gmres.solve(b);
    time = clock() - time;
    std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "MatFree DGMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( gmres.error() > 1e-6)
    {
      std::cout << "Failed line 360" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> gmres;
    gmres.setTolerance(1e-6);
    gmres.setMaxIterations(1000);
    gmres.set_restart(50);
    gmres.compute(S);
    x = gmres.solve(b);
    time = clock() - time;
    std::cout << "DGMRES:   #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    std::cout << "Std DGMRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( gmres.error() > 1e-6)
    {
      std::cout << "Failed line 380" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::MINRES<MatrixReplacement, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.setTolerance(1e-6);
    minres.compute(A);
    x = minres.solve(b);
    time = clock() - time;
    std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
    std::cout << "MatFree MINRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( minres.error() > 1e-6)
    {
      std::cout << "Failed line 398" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;

  {
    time = clock();
    Eigen::MINRES<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> minres;
    minres.setTolerance(1e-6);
    minres.compute(S);
    x = minres.solve(b);
    time = clock() - time;
    std::cout << "MINRES:   #iterations: " << minres.iterations() << ", estimated error: " << minres.error() << std::endl;
    std::cout << "Std MINRES Solve (s):\t" << (time / CLOCKS_PER_SEC) << std::endl;

    if ( minres.error() > 1e-6)
    {
      std::cout << "Failed line 416" << std::endl;
      return -1;
    }
  }
  std::cout << std::endl;
}
