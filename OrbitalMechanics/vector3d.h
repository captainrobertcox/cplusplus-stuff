#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
#include <cmath>

using namespace std;

class Vec3d {

  friend ostream& operator<<(ostream&, const Vec3d&);
  friend istream& operator>>(istream&, Vec3d&);

  public:

    Vec3d(double xi=0.0, double yi=0.0, double zi=0.0, double wi=1.0){ // Constructor
      x=xi; y=yi; z=zi; w=wi;
    }

    Vec3d(const Vec3d &V); // Copy Constructor

    void set(double xi, double yi, double zi, double wi=1.0);
    double length() const;
    void normalize();
    Vec3d& operator-(); // negate
    //void operator=(const Vec3d &V); // assignment overload (not necessary without pointers)
    void operator+=(const Vec3d &V);
    Vec3d operator+(const Vec3d &V) const;
    void operator+=(double d); // add scalar
    Vec3d operator+(double d) const;
    void operator-=(const Vec3d &V);
    Vec3d operator-(const Vec3d &V) const;
    void operator-=(double d); // subtract scalar
    Vec3d operator-(double d) const;
    void operator*=(double d); // scale up
    Vec3d operator*(double d) const;
    void operator/=(double d); // scale down
    Vec3d operator/(double d) const;
    bool operator==(const Vec3d &V) const; // compare

    double operator*(const Vec3d &V) const; // dot product
    Vec3d proj(const Vec3d &V) const; // projection
    void cross(const Vec3d &V1, const Vec3d &V2); // cross product
    Vec3d operator/(const Vec3d &V) const; // also cross product


    double x, y, z, w;
};

ostream& operator<<(ostream& os, const Vec3d& V){
  os << V.x << " " << V.y << " " << V.z;
  return os;
}

istream& operator>>(istream& is, const Vec3d& V){
  is >> V.x >> V.y >> V.z;
  return is;
}

Vec3d::Vec3d(const Vec3d &V){ // copy constructor
  x=V.x; y=V.y; z=V.z; w=V.w;
}

void Vec3d::set(double xi, double yi, double zi, double wi){
  x=xi; y=yi; z=zi; w=wi;
}

double Vec3d::length() const{
  return sqrt(x*x + y*y + z*z);
}

void Vec3d::normalize(){
  double d = sqrt(x*x + y*y + z*z);
  if (d != 0) { x/=d; y/=d; z/=d; }
}


Vec3d& Vec3d::operator-(){
  x=-x; y=-y; z=-z;
  return *this;
}

void Vec3d::operator+=(const Vec3d &V){
  x+=V.x; y+=V.y; z+=V.z;
}

Vec3d Vec3d::operator+(const Vec3d &V) const{
  return Vec3d(x+V.x, y+V.y, z+V.z);
}

void Vec3d::operator+=(double d){
  x+=d; y+=d; z+=d;
}

Vec3d Vec3d::operator+(double d) const{
  return Vec3d(x+d, y+d, z+d);
}

void Vec3d::operator-=(const Vec3d &V){
  x-=V.x; y-=V.y; z-=V.z;
}

Vec3d Vec3d::operator-(const Vec3d &V) const{
  return Vec3d(x-V.x, y-V.y, z-V.z);
}

void Vec3d::operator-=(double d){
  x-=d; y-=d; z-=d;
}

Vec3d Vec3d::operator-(double d) const{
  return Vec3d(x-d, y-d, z-d);
}

void Vec3d::operator*=(double d){
  x*=d; y*=d; z*=d;
}

Vec3d Vec3d::operator*(double d) const{
  return Vec3d(x*d, y*d, z*d);
}

void Vec3d::operator/=(double d){
  x/=d; y/=d; z/=d;
}

Vec3d Vec3d::operator/(double d) const{
  return Vec3d(x/d, y/d, z/d);
}

bool Vec3d::operator==(const Vec3d &V) const{
  return ((abs(x-V.x)<1.0e-20) && (abs(y-V.y)<1.0e-20) && (abs(z-V.z)<1.0e-20));
}

double Vec3d::operator*(const Vec3d &V) const { // dot product, perpendicular=0, <90deg=positive, dot of norm vecs is cosine
  return (x*V.x + y*V.y + z*V.z);
}

Vec3d Vec3d::proj(const Vec3d &V) const{
  double d = (V* *this) / (V*V);
  return Vec3d(V.x*d, V.y*d, V.z*d);
}

void Vec3d::cross(const Vec3d &V1, const Vec3d &V2){
  x = V1.y * V2.z - V1.z * V2.y;
  y = V1.z * V2.x - V1.x * V2.z;
  z = V1.x * V2.y - V1.y * V2.x;
}

Vec3d Vec3d::operator/(const Vec3d &V) const{
  return Vec3d(y*V.z - z * V.y, z*V.x - x*V.z, x*V.y - y*V.x);
}

#endif

