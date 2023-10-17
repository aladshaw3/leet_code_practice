/**
  This example shows how to use the Full Templated
  Matrix object. This is for more advanced control
  over how the matrix is formed and how it is defined
  at compile time.
*/

#include <iostream>
#include <complex>
#include <Eigen/Dense>

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::ColMajor;

int main()
{
  /**
    Below is the Syntax for the full templated Matrix
    -------------------------------------------------

    Matrix<typename Scalar,
         int RowsAtCompileTime,
         int ColsAtCompileTime,
         int Options = 0,
         int MaxRowsAtCompileTime = RowsAtCompileTime,
         int MaxColsAtCompileTime = ColsAtCompileTime>

    The 'Dynamic' Keyword can be used to replace `RowsAtCompileTime`
      and `ColsAtCompileTime`, so that this info can by dynamically
      allocated

    The `Options` is an int flag to determine how the data is stored
      in memory. Choices are: (i) `ColMajor` [default] and (ii) `RowMajor`.

  */
  // NOTE: Dense Matrix MAX size is 100x100
  Matrix<float, Dynamic, Dynamic, ColMajor, 100, 100> m;
  m.resize(10,10);

  if (m.rows() != 10)
  {
    std::cout << "Failed line 40" << std::endl;
    return -1;
  }

  if (m.cols() != 10)
  {
    std::cout << "Failed line 46" << std::endl;
    return -1;
  }

  // NOTE: ONLY the first 3 items are REQUIRED
  //    (i.e., provide the type, rows, and cols)
  Matrix<int, 3, 4> m2;

  if (m2.rows() != 3)
  {
    std::cout << "Failed line 56" << std::endl;
    return -1;
  }

  if (m2.cols() != 4)
  {
    std::cout << "Failed line 62" << std::endl;
    return -1;
  }

  //Can create matrices of complex values
  Matrix< std::complex<double>, 5, 5> m3;
  if (m3.rows() != 5)
  {
    std::cout << "Failed line 70" << std::endl;
    return -1;
  }

  if (m3.cols() != 5)
  {
    std::cout << "Failed line 77" << std::endl;
    return -1;
  }

  // Use the `real` function to get the real portion of
  //  the complex number
  if ( real( m3(0,0) ) != 0)
  {
    std::cout << "Failed line 85" << std::endl;
    return -1;
  }

  // Use the `imag` function to get the imaginary portion of
  //  the complex number
  if ( imag( m3(0,0) ) != 0)
  {
    std::cout << "Failed line 93" << std::endl;
    return -1;
  }

  //Fill in m3
  for (unsigned int i=0; i<m3.rows(); i++)
  {
    m3(i,i) = std::complex<double>(10.0*(double)i,1.0);
  }


  // Vectors can be created as matrices
  Matrix< std::complex<double>, 5, 1> v3, res;

  //Fill in v3
  for (unsigned int i=0; i<m3.rows(); i++)
  {
    v3(i,0) = std::complex<double>(1,0);
  }

  // Matrix multiplication
  res = m3 * v3;

  std::cout << res << std::endl;
  if ( real( res(0,0) ) != 0)
  {
    std::cout << "Failed line 119" << std::endl;
    return -1;
  }

  // Dot product
  std::complex<double> d;
  d = res.dot(v3);

  std::cout << d << std::endl;
  if ( real( d ) != 100)
  {
    std::cout << "Failed line 130" << std::endl;
    return -1;
  }

  std::complex<double> norm;
  norm = res.norm();

  std::cout << norm << std::endl;
  if ( fabs(real( norm ) - 54.8179) > 1e-3 )
  {
    std::cout << "Failed line 140" << std::endl;
    return -1;
  }

  if ( imag( norm ) != 0)
  {
    std::cout << "Failed line 146" << std::endl;
    return -1;
  }

  // Many, many other matrix operations available
  // (https://eigen.tuxfamily.org/dox/group__QuickRefPage.html)

  std::cout << "Finished" << std::endl;
  return 0;
}
