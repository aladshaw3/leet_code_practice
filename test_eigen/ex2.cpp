/**
  This example shows how to specify the sizes of
  matrices and vectors at run time (non-Dynamically),
  at compile time (defined by object names), and
  at run time (fully Dynamically).
*/
#include <iostream>
#include <Eigen/Dense>

// Importing the namespaces
using Eigen::MatrixXd; ///< *Xd = Set at runtime, of type double
using Eigen::VectorXd; ///< *Xd = Set at runtime, of type double

using Eigen::MatrixXf; ///< *Xf = Set at runtime, of type float
using Eigen::Matrix3f; ///< *3f = Set at compile time, of type float

using Eigen::Matrix3d; ///< *3d = Set at compile time, of type double
using Eigen::Vector3d; ///< *3d = Set at compile time, of type double

int main()
{
  //Create a 3x3 matrix of random numbers between -1 and 1
  MatrixXd m = MatrixXd::Random(3,3);

  //Inline creation of another matrix for adding to existing values
  //    Makes a 3x3 matrix with 1.2 in all positions
  m = (m + MatrixXd::Constant(3,3,1.2)) * 50;

  if (m.rows() != 3)
  {
    std::cout << "Failed line 25" << std::endl;
    return -1;
  }

  if (m.cols() != 3)
  {
    std::cout << "Failed line 31" << std::endl;
    return -1;
  }

  std::cout << "m =" << std::endl << m << std::endl;

  //Create a 3x1 vector
  VectorXd v(3);

  if (v.size() != 3)
  {
    std::cout << "Failed line 42" << std::endl;
    return -1;
  }

  //Specify values of vector using the `<<` operator
  v << 1, 2, 3;

  std::cout << "m * v =" << std::endl << m * v << std::endl;

  //Comma initializer
  //    Create a 5x5 matrix
  //    Fill in the first square section 3x3 with the `<<` operator
  //    Fill in the non-square sides with 0s
  //    Fill in remaining diagonal with 1s
  int rows=5, cols=5;
  MatrixXf m2(rows,cols);
  m2 << (Matrix3f() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished(),
     MatrixXf::Zero(3,cols-3),
     MatrixXf::Zero(rows-3,3),
     MatrixXf::Identity(rows-3,cols-3);
  std::cout << m2;
  std::cout << std::endl;
  /**
      1 2 3 0 0
      4 5 6 0 0
      7 8 9 0 0
      0 0 0 1 0
      0 0 0 0 1
  */
  if (m2(1,0) != 4)
  {
    std::cout << "Failed line 77" << std::endl;
    return -1;
  }


  //Create a matrix of blank size
  MatrixXf m3;

  //Dynamically set the size
  m3.resize(4,5);

  //Check
  if (m3.rows() != 4)
  {
    std::cout << "Failed line 91" << std::endl;
    return -1;
  }

  if (m3.cols() != 5)
  {
    std::cout << "Failed line 97" << std::endl;
    return -1;
  }


  // Compile time sizing
  Matrix3d m4 = Matrix3d::Random();
  m4 = (m4 + Matrix3d::Constant(1.2)) * 50;
  std::cout << "m4 =" << std::endl << m4 << std::endl;
  Vector3d v4(1,2,3);

  std::cout << "m4 * v4 =" << std::endl << m4 * v4 << std::endl;

  if (m4.rows() != 3)
  {
    std::cout << "Failed line 112" << std::endl;
    return -1;
  }

  if (m4.cols() != 3)
  {
    std::cout << "Failed line 118" << std::endl;
    return -1;
  }

  if (v4.size() != 3)
  {
    std::cout << "Failed line 124" << std::endl;
    return -1;
  }

  return 0;
}
