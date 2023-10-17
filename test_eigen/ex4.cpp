/**
  This example shows how to use the SparseMatrix object.
  For most practical problems, we will likely want to
  use SparseMatrix for the sake of memory efficiency.
  However, the Dense Matrix class is more operationally
  efficient (if the matrices are relatively small).

  Vectors will generally be dense and can be used in
  conjunction with SparseMatrix objects.
*/


#include <Eigen/Sparse>
#include <vector>
#include <iostream>

// Defining T as Eigen::Triplet<double> for convienience
typedef Eigen::Triplet<double> T;

// Inserts coefficients into b
void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs,
                       Eigen::VectorXd& b, const Eigen::VectorXd& boundary)
{
  int n = int(boundary.size());
  int id1 = i+j*n;

  if(i==-1 || i==n)
    b(id) -= w * boundary(j); // constrained coefficient
  else if(j==-1 || j==n)
    b(id) -= w * boundary(i); // constrained coefficient
  else
    coeffs.push_back(T(id,id1,w)); // unknown coefficient
}

// Builds the b vector
void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n)
{
  b.setZero();
  Eigen::ArrayXd boundary = Eigen::ArrayXd::LinSpaced(n, 0,3.14159).sin().pow(2);
  for(int j=0; j<n; ++j)
  {
    for(int i=0; i<n; ++i)
    {
      int id = i+j*n;
      insertCoefficient(id, i-1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i+1,j, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j-1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j+1, -1, coefficients, b, boundary);
      insertCoefficient(id, i,j,    4, coefficients, b, boundary);
    }
  }
}

int main()
{

  // ------------ Test 1: Setting up and Solving Large Sparse System ----------------------------------
  int n = 300;  // size of the image
  int m = n*n;  // number of unknowns (=number of pixels)

  // Assembly: Collect RHS of Ax = b
  std::vector<T> coefficients;            // list of non-zeros coefficients
  Eigen::VectorXd b(m);                   // the right hand side-vector resulting from the constraints
  buildProblem(coefficients, b, n);

  // Formulate the Sparse version of A
  //  NOTE: This matrix is far too large to create a dense form of
  Eigen::SparseMatrix<double> A(m,m);

  //This forms a Tridiagonal Matrix
  A.setFromTriplets(coefficients.begin(), coefficients.end());

  // Solving: Using Cholesky factorization on a SparseMatrix
  Eigen::SimplicialCholesky< Eigen::SparseMatrix<double> > chol(A);  // performs a Cholesky factorization of A
  Eigen::VectorXd x = chol.solve(b);         // use the factorization to solve for the given right hand side


  // ------------------- Test 2: Basic Matrix Creation ----------------------------------

  // SparseMatrix templated class has 3 args (last 2 are optional)
  /**
        SparseMatrix< type, Options, StorageIndex >

            type: float, double, complex, etc.
            Options:  ColMajor [default], RowMajor
            StorageIndex: int [default], short, etc.
  */

  // Declare a sparse matrix
  Eigen::SparseMatrix<double> B;

  // Size is intially zero
  std::cout << B.rows() << std::endl;
  std::cout << B.cols() << std::endl;
  std::cout << B << std::endl;
  if (B.rows() != 0 && B.cols() != 0)
  {
    std::cout << "Failed line 96" << std::endl;
    return -1;
  }

  // Specify size (without allocation)
  B.conservativeResize(10,10);
  std::cout << B.rows() << std::endl;
  std::cout << B.cols() << std::endl;
  std::cout << B << std::endl;

  if (B.rows() != 10 && B.cols() != 10)
  {
    std::cout << "Failed line 108" << std::endl;
    return -1;
  }

  // There are 2 Ways to insert data into a SparseMatrix

  //  (1) =============== Using 'insert' ============================
  //        This method is generally less efficient, but easy to use.
  B.insert(0,0) = 1;
  B.insert(1,1) = 1;

  // Use the .`coeff` function to grab value from matrix
  //    (if value does not exist, it returns 0)
  if (B.coeff(0,0) != 1 )
  {
    std::cout << "Failed line 121" << std::endl;
    return -1;
  }

  //        This can be improved by pre-allocating memory space
  //        for a known number of non-zeros.
  unsigned int cols = 10;
  unsigned int nonzeros = 3;
  //        Here, we are reserving 3 non-zeros per colomn
  B.reserve( Eigen::VectorXi::Constant(cols,nonzeros) );

  //        After reserving space, iterate through and insert
  //        remaining elements. (NOTE: `insert` is only valid
  //        if the element does NOT already exist).
  for (unsigned int i=2; i<cols; i++)
  {
    B.insert(i,i) = 1;
  }

  //      Lastly, you can re-compress the matrix after formation
  //      (During a typical insertion, more than what is needed
  //      usually allocated to make the insertion faster).
  B.makeCompressed();

  std::cout << B.rows() << std::endl;
  std::cout << B.cols() << std::endl;
  std::cout << B << std::endl;

  //  (2) =============== Using 'setFromTriplets' ============================
  //        This method is usually more complicated, but more efficient
  Eigen::SparseMatrix<double> C;

  // Specify size (without allocation)
  C.conservativeResize(10,10);

  // Recall from earlier, we defined `T` as `Eigen::Triplet<double>`
  //  We can use a list/vector of triplets to specify values.
  //  A Triplet is constructed with 3 args
  //        T( row_index, col_index, value )
  //  We create a list of these, then use `setFromTriplets` to fill
  //    the SparseMatrix...
  std::vector<T> list;
  list.resize(C.rows());
  for (unsigned int i=0; i<C.rows(); i++)
  {
    list[i] = T(i, i, 2); //Fill diagonal with 2s
  }
  C.setFromTriplets(list.begin(), list.end());
  std::cout << C << std::endl;
  if (C.coeff(0,0) != 2 )
  {
    std::cout << "Failed line 174" << std::endl;
    return -1;
  }


  // ------------------- Test 3: Iterating Through SparseMatrix ----------------------------------
  for (unsigned int k=0; k<C.outerSize(); ++k)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(C,k); it; ++it)
    {
      std::cout << "C(" << it.row() << "," << it.col() << ") = " << it.value() << std::endl;

      // NOTE: Inner Index...
      std::cout << "Inner Index: " << it.index() << std::endl << std::endl;
    }
  }


  // Other Notes:
  /**
      Addition of a SparseMatrix by a DenseMatrix results in a DenseMatrix
      Addition of a SparseMatrix by a SparseMatrix results in a SparseMatrix
      Multiplication of a SparseMatrix by a DenseMatrix results in a DenseMatrix
      Multiplication of a SparseMatrix by a SparseMatrix results in a SparseMatrix
      The 'prune' function removes values in SparseMatrix that are 'close to zero'
  */



  std::cout << "Finished" << std::endl;
  return 0;
}
