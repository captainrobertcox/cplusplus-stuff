#ifndef MATRIX_H
#define MATRIX_X

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring> // memcpy, memcmp
#include "vector3d.h"

using namespace std;

#define MAT(M, r, c) (M)[(r)*4 + (c)]
#define m11 MAT(Elements, 0, 0)
#define m12 MAT(Elements, 0, 1)
#define m13 MAT(Elements, 0, 2)
#define m21 MAT(Elements, 1, 0)
#define m22 MAT(Elements, 1, 1)
#define m23 MAT(Elements, 1, 2)
#define m31 MAT(Elements, 2, 0)
#define m32 MAT(Elements, 2, 1)
#define m33 MAT(Elements, 2, 2)

#define a11 MAT(A.Elements, 0, 0)
#define a12 MAT(A.Elements, 0, 1)
#define a13 MAT(A.Elements, 0, 2)
#define a21 MAT(A.Elements, 1, 0)
#define a22 MAT(A.Elements, 1, 1)
#define a23 MAT(A.Elements, 1, 2)
#define a31 MAT(A.Elements, 2, 0)
#define a32 MAT(A.Elements, 2, 1)
#define a33 MAT(A.Elements, 2, 2)

#define b11 MAT(B.Elements, 0, 0)
#define b12 MAT(B.Elements, 0, 1)
#define b13 MAT(B.Elements, 0, 2)
#define b21 MAT(B.Elements, 1, 0)
#define b22 MAT(B.Elements, 1, 1)
#define b23 MAT(B.Elements, 1, 2)
#define b31 MAT(B.Elements, 2, 0)
#define b32 MAT(B.Elements, 2, 1)
#define b33 MAT(B.Elements, 2, 2)

class Matrix3{

  friend ostream& operator<<(ostream& os, const Matrix3& A);
  friend istream& operator>>(istream& is, Matrix3& A);

  public:

    Matrix3();
    Matrix3(const Matrix3& A); // copy constructor
    Matrix3(double d11, double d12, double d13, double d21, double d22, double d23, double d31, double d32, double d33);
    //const Matrix3& operator=(const Matrix3& A); // overload assignment operator
    bool operator==(const Matrix3& A); // comparison
    void operator+=(const Matrix3& A);
    void operator-=(const Matrix3& A);
    void operator*=(const double d);
    void operator*=(const Matrix3& A);
    Matrix3 operator*(const double d) const;
    Matrix3 operator*(const Matrix3& A) const;
    Vec3d operator*(const Vec3d& V) const;
    void operator/=(const double d);
    double operator[](int element);
    void Set(double d11, double d12, double d13, double d21, double d22, double d23, double d31, double d32, double d33);
    void Reset(); // set Elements to zero
    void Identity(); // make this the identity matrix
    void RotateX(const double d); // make this a rotation matrix, radians around x axis
    void RotateY(const double d); // make this a rotation matrix, radians around y axis
    void RotateZ(const double d); // make this a rotation matrix, radians around z axis
    Matrix3& Transpose(); // transpose this matrix
    Matrix3& Transpose(const Matrix3& A); // return transpose of matrix A

  private:

    double Elements[16];

};

ostream& operator<<(ostream& os, const Matrix3& A){
  os << endl;
  os << setw(10) << a11 << " " << setw(10) << a12 << " " << setw(10) << a13 << endl;
  os << setw(10) << a21 << " " << setw(10) << a22 << " " << setw(10) << a23 << endl;
  os << setw(10) << a31 << " " << setw(10) << a32 << " " << setw(10) << a33 << endl;
  return os;
}

istream& operator>>(istream& is, Matrix3& A){
  is >> a11 >> a12 >> a13 >> a21 >> a22 >> a23 >> a31 >> a32 >> a33;
  return is;
}

Matrix3::Matrix3(){
  Reset();
}

Matrix3::Matrix3(const Matrix3& A){
  memcpy(Elements, A.Elements, sizeof(Elements));
}

Matrix3::Matrix3(double d11, double d12, double d13, double d21, double d22, double d23, double d31, double d32, double d33){
  m11 = d11; m12 = d12, m13 = d13;
  m21 = d21; m22 = d22; m23 = d23;
  m31 = d31; m32 = d32; m33 = d33;
}

void Matrix3::Reset(){
  memset(Elements, 0, sizeof(Elements));
}

void Matrix3::Identity(){
  Reset();
  m11 = m22 = m33 = 1.0;
}

bool Matrix3::operator==(const Matrix3& A){
  return memcmp(Elements, A.Elements, sizeof(Elements));
}

void Matrix3::operator+=(const Matrix3& A){
  m11 += a11; m12 += a12; m13 += a13;
  m21 += a21; m22 += a22; m23 += a23;
  m31 += a31; m32 += a32; m33 += a33;
}

void Matrix3::operator-=(const Matrix3& A){
  m11 -= a11; m12 -= a12; m13 -= a13;
  m21 -= a21; m22 -= a22; m23 -= a23;
  m31 -= a31; m32 -= a32; m33 -= a33;
}

void Matrix3::operator*=(const double d){
  m11 *= d; m12 *= d; m13 *= d;
  m21 *= d; m22 *= d; m23 *= d;
  m31 *= d; m32 *= d; m33 *= d;
}

void Matrix3::operator/=(const double d){
  m11 /= d; m12 /= d; m13 /= d;
  m21 /= d; m22 /= d; m23 /= d;
  m31 /= d; m32 /= d; m33 /= d;
}

Matrix3 Matrix3::operator*(const double d) const{
  Matrix3 A;
  a11 = m11*d; a12 = m12*d; a13 = m13*d;
  a21 = m21*d; a22 = m22*d; a23 = m23*d;
  a31 = m31*d; a32 = m32*d; a33 = m33*d;
  return A;
}

Matrix3 Matrix3::operator*(const Matrix3& A) const{
  Matrix3 B;
  b11=m11*a11 + m12*a21 + m13*a31; b12=m11*a12 + m12*a22 + m13*a32; b13=m11*a13 + m12*a23 + m13*a33;
  b21=m21*a11 + m22*a21 + m23*a31; b22=m21*a12 + m22*a22 + m23*a32; b23=m21*a13 + m22*a23 + m23*a33;
  b31=m31*a11 + m32*a21 + m33*a31; b32=m31*a12 + m32*a22 + m33*a32; b33=m31*a13 + m32*a23 + m33*a33;
  return B;
}

void Matrix3::operator*=(const Matrix3& A){
  Matrix3 B;
  memcpy(B.Elements, Elements, sizeof(Elements));
  m11 = b11*a11 + b12*a21 + b13*a31; m12 = b11*a12 + b12*a22 + b13*a32; m13 = b11*a13 + b12*a23 + b13*a33;
  m21 = b21*a11 + b22*a21 + b23*a31; m22 = b21*a12 + b22*a22 + b23*a32; m23 = b21*a13 + b22*a23 + b23*a33;
  m31 = b31*a11 + b32*a21 + b33*a31; m32 = b31*a12 + b32*a22 + b33*a32; m33 = b31*a13 + b32*a23 + b33*a33;
}

Vec3d Matrix3::operator*(const Vec3d& V) const{
  Vec3d A;
  A.x = m11*V.x + m12*V.y + m13*V.z;
  A.y = m21*V.x + m22*V.y + m32*V.z;
  A.z = m31*V.x + m32*V.y + m33*V.z;
  return A;
}

double Matrix3::operator[](int element){
  return Elements[element];
}

void Matrix3::Set(double d11, double d12, double d13, double d21, double d22, double d23, double d31, double d32, double d33){
  m11 = d11; m12 = d12; m13 = d13;
  m21 = d21; m22 = d22; m23 = d23;
  m31 = d31; m32 = d32; m33 = d33;
}

Matrix3& Matrix3::Transpose(){
  Matrix3 A;
  memcpy(A.Elements, Elements, sizeof(Elements));
  m11=a11; m12=a21; m13=a31;
  m21=a12; m22=a22; m23=a32;
  m31=a13; m32=a23; m33=a33;
  return *this;
}

Matrix3& Matrix3::Transpose(const Matrix3& A){
  m11=a11; m12=a21; m13=a31;
  m21=a12; m22=a22; m23=a32;
  m31=a13; m32=a23; m33=a33;
  return *this;
}

void Matrix3::RotateX(const double d){
  // d is in radians, if d is in degrees, convert it to radians by multiplying degrees by PI/180
  m11 = 1.0;
  m22 = cos(d);
  m33 = cos(d);
  m23 = sin(d);
  m32 = -sin(d);
  m12 = m13 = m21 = m31 = 0.0;
}

void Matrix3::RotateY(const double d){
  // d is in radians. If d is in degrees, convert it to radians by multiplying by PI/180
  m11 = cos(d);
  m33 = cos(d);
  m13 = -sin(d);
  m31 = sin(d);
  m22 = 1.0;
  m12 = m21 = m23 = m32 = 0.0;
}

void Matrix3::RotateZ(const double d){
  // d is in radians. If d is in degrees, convert to radians by multiplying by PI/180
  m11 = cos(d);
  m22 = cos(d);
  m12 = sin(d);
  m21 = -sin(d);
  m33 = 1.0;
  m13 = m23 = m31 = m32 = 0.0;
}

#endif

