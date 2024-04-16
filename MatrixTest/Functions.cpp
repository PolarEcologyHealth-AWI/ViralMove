// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
#include <cassert> 
#include <RcppArmadillo.h>

// Namespace for Parameters
namespace sdp {

template <class T>

class checkedVector : public std::vector<T>
{
public:
  
  checkedVector()
  { }
  checkedVector(size_t n, const T& value = T())
    : std::vector<T>(n,value)
    { }
  checkedVector(T* i, T* j)
    : std::vector<T>(i,j)
    { }
  
  T& operator [] (ptrdiff_t index)
  {
    assert (index >= 0 && index < static_cast <ptrdiff_t> (size()));
    return std::vector<T>::operator[](index);
  }
};

// ---
// 2d
// ---

template <class T>
class Matrix
  : public sdp::checkedVector< sdp::checkedVector<T> >
{
protected:
  size_t rows,
  columns;
public:
  Matrix(size_t x = 0, size_t y = 0)
    : sdp::checkedVector< sdp::checkedVector<T> > (x,
      sdp::checkedVector<T>(y)), rows(x), columns(y)
      {}
  
  size_t Rows() const {return rows;}
  size_t Columns() const {return columns;}
  
  void init(const T& Value) {
    for (size_t i=0; i< rows; ++i)
      for (size_t j=0; j < columns; ++j)
        sdp::Matrix<T>::operator[](i)(j) = Value;
  }
  
  void resize(size_t x, size_t y, T t=T()) {
    sdp::checkedVector< sdp::checkedVector<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix<T>::operator[](i).resize(y, t);
    rows = x ; columns = y;
  }
  
  Matrix<T>& I()
  {
    for(size_t i = 0; i < rows; ++i)
      for(size_t j = 0; j < columns; ++j)
        sdp::Matrix<T>::operator[](i)(j) = (i == j) ? T(1) : T(0);
    return *this;
  }
  
};

// ---
// 3d
// ---

template <class T>
class Matrix3D : public sdp::checkedVector<Matrix<T> >
{
protected:
  size_t rows,
  columns,
  zDim;
public:
  Matrix3D(size_t x = 0, size_t y = 0, size_t z = 0)
    : sdp::checkedVector<Matrix<T> >(x, Matrix<T>(y, z)),
      rows(x), columns(y), zDim(z)
      {}
  
  size_t Rows()   const {return rows;}
  size_t Columns() const {return columns;}
  size_t zDIM()   const {return zDim;}
  
  void init(const T& value)
  {
    for (size_t i = 0; i < rows; ++i)
      sdp::Matrix3D<T>::operator[](i).init(value);
  }
  
  void resize (size_t x, size_t y, size_t z, T t=T()) {
    sdp::checkedVector< Matrix<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix3D<T>::operator[](i).resize(y,z,t);
    rows = x ; columns = y; zDim = z;
  }
};


// ---
// 4d
// ---

template <class T>
class Matrix4D : public sdp::checkedVector<Matrix3D<T> >
{
protected:
  size_t rows,
  columns,
  zDim,
  z2Dim;
public:
  Matrix4D(size_t x = 0, size_t y = 0, size_t z = 0, size_t z2 = 0)
    : sdp:: checkedVector<Matrix3D<T> >(x, Matrix3D<T>(y, z, z2)),
      rows(x), columns(y), zDim(z), z2Dim(z2)
      {}
  
  size_t Rows()    const {return rows;}
  size_t Columns() const {return columns;}
  size_t zDIM()    const {return zDim;}
  size_t z2DIM()   const {return z2Dim;}
  
  void init(const T& value)
  {
    for (size_t i = 0; i < rows; ++i)
      sdp::Matrix4D<T>::operator[](i).init(value);
  }
  
  void resize (size_t x, size_t y, size_t z, size_t z2, T t=T())
  {
    sdp::checkedVector< Matrix3D<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix4D<T>::operator[](i).resize(y,z,z2,t);
    rows = x ; columns = y; zDim = z; z2Dim = z2;
  }
};

}



// [[Rcpp::export]]
sdp::Matrix4D<float> timesTwo(double time, double sites, double maxX, double infect) {
  sdp::Matrix4D<float>Matrix;
  Matrix.resize(time, sites, maxX, infect, 0.0);
  return Matrix;
}

