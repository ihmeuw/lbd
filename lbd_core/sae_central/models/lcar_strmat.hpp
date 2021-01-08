// Function for preparing a LCAR RE structure matrix given a neighborhood graph and rho parameter
template<class Type> 
SparseMatrix<Type> lcar_strmat(SparseMatrix<Type> graph, Type rho) {
  SparseMatrix<Type> K = rho * graph; 
  for (size_t i = 0; i < K.rows(); i++)
    K.coeffRef(i,i) += (1 - rho);
  return K; 
}
