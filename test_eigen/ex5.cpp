/**
  This example shows how to use the Basic Linear Algebra
  solvers available in Eigen for both Dense and Sparse
  Matrix objects.

  The Basic Linear Algebra solvers cover:
    (a) Direct solvers (LU, QR, etc)
    (b) Eigen Value solvers
    (c) Iterative solvers (simple)

  Advanced Solvers will include matrix-free methods.
*/

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include<Eigen/IterativeLinearSolvers>
#include <complex>
#include <vector>

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
  // ------------------- Test 1: Dense Matrix Solvers ----------------------
  //    How to solve Ax = b using LU and QR decomposition

  // Declare your matrices (any way you want)
  Eigen::Matrix<float, 3, 3> A;
  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> b;
  b.resize(3,1);
  Eigen::Matrix<float, 3, 1> x1, x2, x3, x4;

  // Fill in the matrices
  A << 1,2,3,  4,5,6,  7,8,10;
  b << 3, 3, 4;
  std::cout << "Here is the matrix A:\n" << A << std::endl;
  std::cout << "Here is the vector b:\n" << b << std::endl;

  // Call a solver
  //    (a) Column Pivioting QR - Slower, but more stable
  x1 = A.colPivHouseholderQr().solve(b);
  std::cout << "The solution is:\n" << x1 << std::endl;
  //    (b) Full Pivioting LU - Much slower, but more stable
  x2 = A.fullPivLu().solve(b);
  std::cout << "The solution is:\n" << x2 << std::endl;
  //    (c) Standard QR - faster, but less stable
  x3 = A.householderQr().solve(b);
  std::cout << "The solution is:\n" << x3 << std::endl;
  //    (b) Partial Pivioting LU - Much faster, but less stable
  x4 = A.partialPivLu().solve(b);
  std::cout << "The solution is:\n" << x4 << std::endl;

  // Validate solution
  if ( ( x1 - x2 ).norm() > 1e-3 )
  {
    std::cout << "Failed line 48" << std::endl;
    return -1;
  }
  if ( ( x2 - x3 ).norm() > 1e-3 )
  {
    std::cout << "Failed line 54" << std::endl;
    return -1;
  }
  if ( ( x3 - x4 ).norm() > 1e-3 )
  {
    std::cout << "Failed line 59" << std::endl;
    return -1;
  }

  // Several other solvers, but these are the most generic
  // Full list here: https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html



  // ------------------- Test 2: Dense Matrix Eigenvalues ----------------------
  //    How to find the Eigenvalues/Eigenvectors of A
  //        NOTE: The Eigenvalues && Eigenvectors on return are of type std::complex< T >
  Eigen::Matrix< std::complex<float> , 3, 1> lam;
  Eigen::Matrix< std::complex<float> , 3, 3> lam_mat;

  // Declare an EigenSolver object of same type as the matrix object
  // we are solving the eigenvalues of.
  //
  //    Other Options include:
  //          SelfAdjointEigenSolver - For matrices that are Self-adjoint/Symmetric Only
  //                (Fast, but conditional on symmetry)
  //          EigenSolver - For matrices that are Square and Real
  //                (Average, but real matrices only)
  //          GeneralizedSelfAdjointEigenSolver - For matrices that are Square
  //                (Fast, but only for A*v = lam*B*v)
  //          ComplexEigenSolver - For matrices that are Square
  //                (Slow, but most generic)

  // Either of these work for this problem...
  Eigen::EigenSolver< Eigen::Matrix<float, 3, 3> > e_solver;
  //Eigen::ComplexEigenSolver< Eigen::Matrix<float, 3, 3> > e_solver;

  // Call the 'compute' method passing the matrix to solve 'A'
  //  and boolean determining whether or not to also compute
  //  the eigenvectors. Select 'true' to get the eigenvectors.
  e_solver.compute(A, true);

  // After the solve, you can retrive the eigenvalues and
  // eigenvectors with `.eigenvalues()` and `.eigenvectors()`
  lam = e_solver.eigenvalues();
  std::cout << lam << std::endl << std::endl;
  lam_mat = e_solver.eigenvectors();
  std::cout << lam_mat << std::endl << std::endl;

  if ( fabs(lam.norm() - 16.7332) > 1e-3)
  {
    std::cout << "Failed line 91" << std::endl;
    return -1;
  }

  // ------------------- Test 3: Sparse Matrix Solvers ----------------------
  //    How to solve As*xs = bs using Sparse versions of LU and QR decomposition
  Eigen::SparseMatrix<double> As;
  int m = 4; // 3D size var
  fillNatural3DLaplacian(m, As);
  std::cout << As << std::endl;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bs, xs;
  xs.resize(m*m*m,1);
  fillNatural3DLaplacian_BCs(m, bs);
  //std::cout << bs << std::endl;

  //    (a) Solve with LU
  //        Object MUST contain the type of matrix we are solving, as well
  //        as the OrderMethod to use. The OrderMethods are as follows:
  //            COLAMDOrdering<int>  - Column Major [Default]
  //            AMDOrdering<int>     - Requires Symmetry
  //            NaturalOrdering<int> - Ordering of the matrix
  Eigen::SparseLU< Eigen::SparseMatrix<double> > lu_solver;
  // Analyze the pattern - DO ONLY ONCE if sparsity pattern doesn't change
  //      for this step the numerical values of A are not used
  lu_solver.analyzePattern(As);

  // This is the actual 'solve' step, the 'solve' grabs the solution (if available)
  lu_solver.factorize(As);
  if (lu_solver.info() != Eigen::ComputationInfo::Success)
  {
    std::cout << "LU decomposition Failed" << std::endl;
    return -1;
  }
  xs = lu_solver.solve(bs);
  if (lu_solver.info() != Eigen::ComputationInfo::Success)
  {
    std::cout << "LU solve Failed" << std::endl;
    return -1;
  }

  double norm = (As*xs - bs).norm();
  std::cout << norm << std::endl;
  if ( norm > 1e-6)
  {
    std::cout << "Failed line 251" << std::endl;
    return -1;
  }


  //    (b) Solve with QR
  //        Object MUST contain the type of matrix we are solving, as well
  //        as the OrderMethod to use. The OrderMethods are as follows:
  //            COLAMDOrdering<int>  - Column Major [Default]
  //            AMDOrdering<int>     - Requires Symmetry
  //            NaturalOrdering<int> - Ordering of the matrix
  Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > qr_solver;
  // Analyze the pattern - DO ONLY ONCE if sparsity pattern doesn't change
  //      for this step the numerical values of A are not used
  qr_solver.analyzePattern(As);

  // This is the actual 'solve' step, the 'solve' grabs the solution (if available)
  qr_solver.factorize(As);
  if (qr_solver.info() != Eigen::ComputationInfo::Success)
  {
    std::cout << "QR decomposition Failed" << std::endl;
    return -1;
  }
  xs = qr_solver.solve(bs);
  if (qr_solver.info() != Eigen::ComputationInfo::Success)
  {
    std::cout << "QR solve Failed" << std::endl;
    return -1;
  }

  norm = (As*xs - bs).norm();
  std::cout << norm << std::endl;
  if ( norm > 1e-6)
  {
    std::cout << "Failed line 285" << std::endl;
    return -1;
  }


  // ------------------- Test 4: Sparse Matrix Solvers ----------------------
  //    How to solve As*xs = bs using PCG, BiCGSTAB, and LSCG

  //      (a) Using PCG
  //            ConjugateGradient< Matrix type, TrigPart, Precon>
  //                TrigPart = Lower [default], Lower|Upper, Upper
  //                Precon = DiagonalPreconditioner [default]
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;

  // As before, you call a 'compute()' function before calling solve
  cg.compute(As);
  xs = cg.solve(bs);
  std::cout << "#iterations:     " << cg.iterations() << std::endl;
  std::cout << "estimated error: " << cg.error()      << std::endl;
  if ( cg.error() > 1e-6)
  {
    std::cout << "Failed line 307" << std::endl;
    return -1;
  }


  //      (b) Using BiCGSTAB
  //            BiCGSTAB< Matrix type, Precon>
  //                Precon = DiagonalPreconditioner [default]
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> bicg;

  // As before, you call a 'compute()' function before calling solve
  bicg.compute(As);
  xs = bicg.solve(bs);
  std::cout << "#iterations:     " << bicg.iterations() << std::endl;
  std::cout << "estimated error: " << bicg.error()      << std::endl;
  if ( bicg.error() > 1e-6)
  {
    std::cout << "Failed line 324" << std::endl;
    return -1;
  }


  //      (c) Using LSCG
  //            LeastSquaresConjugateGradient< Matrix type, Precon>
  //                Precon = DiagonalPreconditioner [default]
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> lscg;

  // As before, you call a 'compute()' function before calling solve
  lscg.compute(As);
  xs = lscg.solve(bs);
  std::cout << "#iterations:     " << lscg.iterations() << std::endl;
  std::cout << "estimated error: " << lscg.error()      << std::endl;
  if ( lscg.error() > 1e-6)
  {
    std::cout << "Failed line 341" << std::endl;
    return -1;
  }

}
