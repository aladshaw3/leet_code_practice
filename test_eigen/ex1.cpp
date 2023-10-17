/**
  This example shows how to create a basic matrix.
*/

#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);

  // Test each value to make sure it is accurate
  if (m(0,0) != 3)
  {
    std::cout << "Failed line 15" << std::endl;
    return -1;
  }

  if (m(1,0) != 2.5)
  {
    std::cout << "Failed line 21" << std::endl;
    return -1;
  }

  if (m(0,1) != -1)
  {
    std::cout << "Failed line 27" << std::endl;
    return -1;
  }

  if (m(1,1) != 1.5)
  {
    std::cout << "Failed line 33" << std::endl;
    return -1;
  }

  return 0;
}
