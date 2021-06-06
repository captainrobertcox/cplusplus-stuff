#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "matrix3.h"

//#define MAT(m, r, c) (m)[(r)*4+c]
#define d11 MAT(D, 0, 0)
#define d12 MAT(D, 0, 1)
#define d13 MAT(D, 0, 2)
#define d21 MAT(D, 1, 0)
#define d22 MAT(D, 1, 1)
#define d23 MAT(D, 1, 2)
#define d31 MAT(D, 2, 0)
#define d32 MAT(D, 2, 1)
#define d33 MAT(D, 2, 2)

using namespace std;

const double PI = 3.14159265;
const double deg = PI/180.0;
const double g0 = .00981; // km/s^2
const double Gconst = 6.67408e-20; // km^3/(kg*s^2)
const double MassEarth = 5.974e24;
const double RadiusEarth = 6378.13655;
const double EarthEquatorialRotationVel = .46510105; // km/s
const double KSC_Rotation_Vel = 0.408456; //0.40839; // km/s or 1470.2 km/h
const double KSC_Latitude = 28.5729; // deg N
const double KSC_longitude = 80.6490; // deg W
const double FLAT = 1.0/298.256421867;
// Mu = G * MassOfBody, ie. for Earth: 6.674e-20 * 5.974e24 = 398600 km^3/s^2
const double MU = 398600.4418; // Earth's Sphere of influece radius is 925000 km
const double MU_SUN = 1.327124e11;
const double MU_MOON = 4903.0; // Moon's SOI radius is 66200 km
const double MU_MARS = 42828.0; // Mars' SOI radius is 577000 km

struct orbitalElements{
  double h; // angular momentum (km^2/s)
  double e; // eccentricity
  double RA; // right ascension of the ascending node (deg)
  double incl; // inclination (deg)
  double w; // argument of perihelion (deg)
  double TA; // true anomaly (deg)
  double a; // semimajor axis (km)
} oe;

struct planetElements{
  orbitalElements oe;
  double w_hat; // longitude of perihelion (= RA + w) (deg)
  double L; // mean longitude (= w_hat + M)  (deg)
  double M; // mean anomaly (deg)
  double E; // eccentric anomaly (deg)
} pe;


string months[]={"January","February","March", "April","May","June","July","August","September","October","November","December"};
string planets[]={"Mercury","Venus","Earth","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"};

void month_planet_names(string& month, string& planet, int month_id, int planet_id);
double zero_to_360(double x){if(x<0.0) while (x<0.0) x+=360.0; else if (x>=360.0) while (x>=360.0) x-=360.0; return x;}
double kepler_E(double e, double M);
double kepler_H(double e, double M);
double stumpS(double z);
double stumpC(double z);
double S(double z){return stumpS(z);}
double C(double z){return stumpC(z);}
double kepler_U(double dt, double ro, double vro, double a, double mu=MU);
void f_and_g(double &f, double &g, double x, double t, double ro, double a, double mu=MU);
void fdot_and_gdot(double &fdot, double &gdot, double x, double r, double ro, double a, double mu=MU);
void RV_from_R0V0(Vec3d &R, Vec3d &V, Vec3d R0, Vec3d V0, double t, double mu=MU);
void displayOE(orbitalElements oe, double mu=MU);
double period(double a, double mu=MU){return 2*PI/sqrt(mu)*pow(a,1.5);}
orbitalElements oe_from_sv(Vec3d R, Vec3d V, double mu=MU);
void sv_from_oe(Vec3d& R, Vec3d& V, orbitalElements oe, double mu=MU);
bool gibbs(Vec3d& V2, Vec3d R1, Vec3d R2, Vec3d R3, double mu=MU);
double y(double z, double r1, double r2, double A);
double F(double z, double t, double r1, double r2, double A, double mu=MU);
double dFdz(double z, double r1, double r2, double A);
double periapseRadius(double h, double e, double mu=MU){return h*h/mu/(1.0 + e);}
void lambert(Vec3d& V1, Vec3d& V2, Vec3d R1, Vec3d R2, double t, string orbitDirection, double mu=MU);
double julian(int year, int month, int day, int hour=0, int min=0, double sec=0.0);
double lst(int year, int month, int day, int hour, int min, double sec, int lon_deg, int lon_min, double lon_sec);
void RV_from_observe(Vec3d& R, Vec3d& V, double rho, double rhodot, double A, double Adot, double a, double adot, double theta, double phi, double H, double f=FLAT, double Re=RadiusEarth, double wE=7.292115e-5, double mu=MU);
void gauss(Vec3d& R, Vec3d& V, Vec3d& R_old, Vec3d& V_old, double dec1, double RA1, double dec2, double RA2, double dec3, double RA3, double Altitude, double Latitude, double LST1, double LST2, double LST3, double t1, double t2, double t3, double Re=RadiusEarth, double flat=FLAT, double mu=MU);
void planet_elements(planetElements& pe, Vec3d& R, Vec3d& V, double& jd, int planet_id, int year, int month, int day, int hour, int minute, double second, double mu=1.327124e11);
void planet_data(double J200_oe[], double rates[], int planet_id);
void interplanetary(Vec3d& Rp1, Vec3d& Vp1, double& jd1, Vec3d& Rp2, Vec3d& Vp2, double& jd2, Vec3d& V1, Vec3d& V2, int depart[], int arrive[], double mu=MU_SUN);
double rocket(double e1, double isp1, double e2, double isp2, double e3, double isp3, double payload, double dV);
double payload(double mE1, double mp1, double isp1, double mE2, double mp2, double isp2, double Vbo);
double deltaV(double mE1, double mp1, double isp1, double mE2, double mp2, double isp2, double mE3, double mp3, double isp3, double mPL);
double escapeV(double altitude = 0.0, double radiusBody = RadiusEarth, double massBody = MassEarth){ return sqrt(2*Gconst*massBody/(radiusBody+altitude));}
double orbitalV(double altitude = 0.0, double radiusBody = RadiusEarth, double massBody = MassEarth){ return sqrt(Gconst*massBody/(radiusBody + altitude));}
double acheiveOrbitDeltaV(double orbitAltitude = 0.0, double latitude = 0.0, double radiusBody = RadiusEarth, double bodyRotationVel = EarthEquatorialRotationVel, double g = g0){return sqrt(g*(radiusBody + orbitAltitude)) - (latitude!=0?(bodyRotationVel*cos(latitude*deg)):0);}

int main(){
/*
  double E, F, x;
  E = kepler_E(0.37255, 3.6029);
  cout << "Eccentric anomaly (radians): " << E << endl;

  F = kepler_H(2.7696, 40.69);
  cout << "Hyperbolic eccentric anomaly (dimensionless): " << F << endl;

  x = kepler_U(3600.0, 10000.0, 3.0752, 1.0 / -19655.0);
  cout << "Universal anomaly (km^0.5): " << x << endl;
*/
  Vec3d R0(7000.0, -12124.0, 0.0);
  Vec3d V0(2.6679, 4.6210, 0.0);

  Vec3d R,V;

  RV_from_R0V0(R, V, R0, V0, 3600.0);
  cout << "Position Vector (km): " << R << ", Velocity Vector (km/s): " << V << endl;
/*
  R.set(-6045.0, -3490.0, 2500.0);
  V.set(-3.457, 6.618, 2.533);
  oe = oe_from_sv(R,V);
  displayOE(oe);
*/
  oe.h = 80000.0;
  oe.e = 1.4;
  oe.RA = 40.0*deg;
  oe.incl = 30.0*deg;
  oe.w = 60.0*deg;
  oe.TA = 30.0*deg;
  sv_from_oe(R, V, oe);
  cout << "State Vector R (km): " << R << endl;
  cout << "State Vector V (km/s): " << V << endl;
  oe = oe_from_sv(R,V);
  displayOE(oe);
/*
  Vec3d R1(-294.32, 4265.1, 5986.7);
  Vec3d R2(-1365.4, 3637.6, 6346.8);
  Vec3d R3(-2940.3, 2473.7, 6555.8);
  Vec3d V2;
  bool err;
  err = gibbs(V2, R1, R2, R3);
  cout << "Velocity Vector, V2 (km/s): " << V2 << endl;
  if (err) cout << "The vectors are not coplanar." << endl;
  else {
    oe = oe_from_sv(R2, V2);
    displayOE(oe);
  }

  R1.set(5000.0, 10000.0, 2100.0);
  R2.set(-14600.0, 2500.0, 7000.0);
  double dt = 3600.0;
  Vec3d V1;
  lambert(V1, V2, R1, R2, dt, "pro");
  oe = oe_from_sv(R1, V1);
  double TA1 = oe.TA;
  oe = oe_from_sv(R2, V2);
  double TA2 = oe.TA;
  cout << "Velocity vector V1 (km): " << V1 << endl;
  cout << "Velocity vector V2 (km): " << V2 << endl;
  cout << "Angular momentum (km^2/s): " << oe.h << endl;
  cout << "Eccentricity: " << oe.e << endl;
  cout << "Inclination (deg): " << oe.incl/deg << endl;
  cout << "RA of ascending node (deg): " << oe.RA/deg << endl;
  cout << "Argument of perigee (deg): " << oe.w/deg << endl;
  cout << "Initial True anomaly (deg): " << TA1/deg << endl;
  cout << "Final True anomaly (deg): " << TA2/deg << endl;
  cout << "Semimajor axis (km): " << oe.a << endl;
  cout << "Periapse radius (km): " << oe.h*oe.h/MU/(1.0 + oe.e) << endl;
  if (oe.e < 1.0){ //ellipse
    double T = 2*PI/sqrt(MU)*pow(oe.a,1.5);
    cout << "Period (s): " << T << endl;
    int t = static_cast<int>(T);
    cout << "Or " << t/24/3600 << " days, " << (t%86400)/3600 << " hours, " << (t%3600)/60 << " minutes, " << (t%60) << " seconds." << endl;
  }

  cout << fixed << showpoint << "Julian day for 5/12/2004 14:45:30.0: " << julian(2004,5,12,14,45,30.0) << endl;

  double LST = lst(2004,3,3,4,30,0.0,139,47,0.0);
  cout << "LST for 3/3/2004 4:30:00, longitude 139 deg, 47' 0\": " << LST << " deg, or " << LST/15.0 << " hr." << endl;

  RV_from_observe(R0, V0, 2551.0, 0.0, 90.0, 0.1130, 30.0, 0.05651, 300.0, 60.0, 0.0);
  oe = oe_from_sv(R0, V0);
  cout << "State Radius Vector (km): " << R0 << endl;
  cout << "State Velocity Vector (km/s): " << V0 << endl;
  displayOE(oe);

  // given lst's
  gauss(R1, V1, R0, V0, -8.78334, 43.5365, -12.0739, 54.4196, -15.1054, 64.3178, 1.0, 40, 44.5065, 45.0, 45.4992, 0.0, 118.104, 237.577);
  cout << "Without iterative improvement:" << endl << "Rold: " << R0 << ", Vold: " << V0 << endl;
  oe = oe_from_sv(R0, V0);
  displayOE(oe);
  cout << "With iterative improvement:" << endl << "R: " << R1 << ", V: " << V1 << endl;
  oe = oe_from_sv(R1, V1);
  displayOE(oe);

  double jd;
  planet_elements(pe, R0, V0, jd, 3, 2003, 8, 27, 12, 0, 0.0);
  cout << "Planet: " << planets[3-1] << endl;
  cout << "Julian date: " << jd << endl;
  cout << "Orbital Elements" << endl;
  displayOE(pe.oe, MU_SUN);
  cout << "Longitude of perihelion: " << pe.w_hat << endl;
  cout << "Mean longitude: " << pe.L << endl;
  cout << "Mean anomaly: " << pe.M << endl;
  cout << "Eccentric anomaly: " << pe.E << endl;
  cout << "State vectors: R: " << R0 << ", V: " << V0 << endl;
  cout << "Magnitude of R: " << R0.length() << ", Magnitude of V: " << V0.length() << endl;

  int depart[7]={3, 1996, 11, 7, 0, 0, 0}; // depart earth 11/7/96 at midnight
  int arrive[7]={4, 1997, 9, 12, 0, 0, 0}; // arrive mars 9/12/97 at midnight
  double jd1, jd2;
  Vec3d Rp1, Vp1, Rp2, Vp2;
  interplanetary(Rp1, Vp1, jd1, Rp2, Vp2, jd2, V1, V2, depart, arrive);
  double tof = jd2 - jd1;
  Vec3d Vinf1, Vinf2;
  Vinf1 = V1 - Vp1;
  Vinf2 = V2 - Vp2;
  oe = oe_from_sv(Rp1, V1, MU_SUN);
  cout << "Depart " << planets[depart[0]-1] << " on Julian day: " << jd1 << endl;
  cout << "Planet position: " << Rp1 << ", Magnitude: " << Rp1.length() << endl;
  cout << "Planet velocity: " << Vp1 << ", Magnitude: " << Vp1.length() << endl;
  cout << "Spacecraft velocity: " << V1 << ", Magnitude: " << V1.length() << ", V-inf at departure: " << Vinf1 << ", Magnitude: " << Vinf1.length() << endl;
  cout << "Time of flight (days): " << tof << endl;
  cout << "Arrive " << planets[arrive[0]-1] << " on Julian day: " << jd2 << endl;
  cout << "Planet position: " << Rp2 << ", Magnitude: " << Rp2.length() << endl;
  cout << "Planet velocity: " << Vp2 << ", Magnitude: " << Vp2.length() << endl;
  cout << "Spacecraft velocity: " << V2 << ", Magnitude: " << V2.length() << ", V-inf at arrival: " << Vinf2 << ", Magnitude: " << Vinf2.length() << endl;
  oe = oe_from_sv(Rp1, V1, MU_SUN);
  cout << "Orbital elements of trajectory at departure:" << endl;
  displayOE(oe, MU_SUN);
  oe = oe_from_sv(Rp2, V2, MU_SUN);
  cout << "True anomaly at arrival (deg): " << oe.TA/deg << endl;

  double m0;
  m0 = rocket(0.10, 400.0, 0.15, 350.0, 0.20, 300.0, 5000, 10);
  cout << "Total rocket mass: " << m0 << endl;
  m0 = rocket(0.2, 300.0, 0.3, 235.0, 0.0, 0.0, 10.0, 6.2);
  cout << "Total rocket mass: " << m0 << endl;
  //m0 = rocket(0.142857, 290.0, 0.047619, 450.0, 0.0, 0.0, 114000, 9.686);
  //cout << "Total rocket mass (2,424,000): " << m0 << endl;
  double mPL = payload(150000.0, 900000.0, 290.0, 60000.0, 1200000.0, 450.0, 8.094 + 2.0 - KSC_Rotation_Vel);
  cout << "Payload mass: " << mPL << endl;

  double Vbo = deltaV(150000.0, 900000.0, 290.0, 60000.0, 1200000.0, 450.0, 0.0, 0.0, 0.0, 114000.0);
  cout << "Delta-V at burnout: " << Vbo << endl;

  cout << "Escape velocity from Earth's surface: " << escapeV() << endl;
  cout << "Orbital velocity at 300 km: " << orbitalV(300.0) << endl;
  cout << "Delta-V to acheive 300 km orbit (add 20% for drag, gravity, and buffer): " << acheiveOrbitDeltaV(300) << endl;
  cout << "Delta-V to 300 km orbit, heading East from KSC (latitude 28.5729): " << acheiveOrbitDeltaV(300.0, KSC_Latitude) << endl;
*/
  return 0;
}

void month_planet_names(string& month, string& planet, int month_id, int planet_id){
  if (1 <= month_id && month_id <= 12) month = months[month_id-1];
  if (1 <= planet_id && planet_id <= 9) planet = planets[planet_id-1];
}

void displayOE(orbitalElements oe, double mu){
  cout << "Angular momentum, h (km^2/s): " << oe.h << endl;
  cout << "Eccentricity, e: " << oe.e << endl;
  cout << "Inclination, i (deg): " << oe.incl/deg << endl; // divide by deg to convert radians to deg
  cout << "Right Ascension of ascending node (Capital Omega) (deg): " << oe.RA/deg << endl;
  cout << "Argument of perigee, w (deg): " << oe.w/deg << endl;
  cout << "True anomaly (Capital theta) (deg): " << oe.TA/deg << endl;
  cout << "Semimajor axis, a (km): " << oe.a << endl;
  cout << "Periapse radius (km): " << periapseRadius(oe.h, oe.e, mu) << endl;
  if (oe.e < 1.0){ // ellipse
    double T = period(oe.a, mu);
    cout << "Period (s): " << T << endl;
    int t = static_cast<int>(T);
    cout << "Or " << t/24/3600 << " days, " << (t%86400)/3600 << " hours, " << (t%3600)/60 << " minutes, " << (t%60) << " seconds." << endl;
  }
}

double kepler_E (double e, double M){
  /* This function uses Newton's method to solve Kepler's equation E - e*sin(E) = M for the eccentric anomaly, given the eccentricity, e, and the mean anomaly, M.

  E - eccentric anomaly (radians)
  e - eccentricity
  M - mean anomaly (radians)
  */

  double error = 1.0e-8;
  double ratio = 1.0;
  double E;

  if (M < PI) E = M + e/2.0;
  else E = M - e/2.0;

  while (abs(ratio) > error){
    ratio = (E - e*sin(E) - M) / (1.0 - e*cos(E));
    E = E - ratio;
  }

  return E;
}

double kepler_H(double e, double M){
  /* This function uses Newton's Method to solve Kepler's equation for the hyperbola e*sinh(F) - F = M for the hyperbolic eccentric anomaly, given the eccentricity e, and the hyperbolic mean anomaly M.
  F - hyperbolic eccentric anomaly (radians)
  e - eccentricity
  M - hyperbolic mean anomaly (radians)
  */

  double error = 1.0e-8;
  double ratio = 1.0;
  double F;

  F = M;
  while (abs(ratio) > error){
    ratio = (e*sinh(F) - F - M) / (e*cosh(F) - 1.0);
    F = F - ratio;
  }

  return F;
}

double stumpS(double z){
  /* This function calculates the stump function S(z) */

  double S;

  if (z > 0) S = (sqrt(z) - sin(sqrt(z))) / pow(sqrt(z),3);
  else if (z < 0)  S = (sinh(sqrt(-z)) - sqrt(-z)) / pow(sqrt(-z), 3);
  else S = 1.0 / 6.0;

  return S;
}

double stumpC(double z){
  /* This function calculates the stump function C(z) */

  double C;

  if (z > 0)  C = (1.0 - cos(sqrt(z))) / z;
  else if (z < 0) C = (cosh(sqrt(-z)) - 1.0) / -z;
  else C = 1.0 / 2.0;

  return C;
}

double kepler_U(double dt, double ro, double vro, double a, double mu){
  /* This function uses Newton's Method to solve the universal Kepler equation for the universal anomaly.
  mu - gravitational parameter (km^3/s^2)
  x - universal anomaly (km^0.5)
  dt - time since x = 0 (s)
  ro - radial position at x = 0 (km)
  vro - radial velocity at x = 0 (km/s)
  a - reciprocal of the semi-major axis (1/km)
  z - auxiliary variable (z=a*x^2)
  C - value of stump function C(z)
  S - value of stump function S(z)
  n - number of iterations for convergence
  nMax - maximum allowable number of iterations
  */

  double error = 1.0e-8;
  int nMax = 1000;
  int n = 0;
  double ratio = 1.0;
  double x = sqrt(mu) * abs(a) * dt;
  double C, S, F, dFdx;

  while (abs(ratio) > error && n <= nMax){
    n++;
    C = stumpC(a*x*x);
    S = stumpS(a*x*x);
    F = ro*vro / sqrt(mu)*x*x*C + (1.0 - a*ro)*pow(x,3)*S + ro*x - sqrt(mu)*dt;
    dFdx = ro*vro / sqrt(mu)*x*(1.0 - a*x*x*S) + (1.0 - a*ro)*x*x*C + ro;
    ratio = F/dFdx;
    x = x - ratio;
  }

  if (n > nMax){
    cout << "Number of iterations: " << n << endl;
    cout << "F/dFdx: " << ratio << endl;
  }

  return x;
}

void f_and_g(double &f, double &g, double x, double t, double ro, double a, double mu){
  /* This function calculates the Lagrange f and g coefficients.
  mu - the gravitational parameter (km^3/s^2)
  a - reciprocal of the semi-major axis (1/km)
  ro - the radial position at time t (km)
  t - the time elapsed since t (s)
  x - the universal anomaly after time t (km^0.5)
  f - the Lagrange f coefficient (dimensionless)
  g - the Lagrange g coeffiecint (s)
  */

  double z = a*x*x;

  f = 1 - x*x/ro*stumpC(z);
  g = t - 1/sqrt(mu)*pow(x,3)*stumpS(z);

  return;
}

void fdot_and_gdot(double &fdot, double &gdot, double x, double r, double ro, double a, double mu){
  /* This function calculates the time derivatives of the Lagrange f and g coefficients.
  mu - the gravitational parameter (km^3/s^2)
  a - the reciprocal of the semimajor axis (1/km)
  ro - the radial position at time t (km)
  t - the time elapsed since initial state vector (s)
  r - the radial position after time t (km)
  x - the universal anomaly after time t (km^0.5)
  fdot - the time derivative of the Lagrange f coefficient (1/s)
  gdot - the time derivative of the Lagrange g coefficient (dimensionless)
  */

  double z = a*x*x;

  fdot = sqrt(mu)/r/ro*(z*stumpS(z) - 1)*x;
  gdot = 1 - x*x/r*stumpC(z);

  return;
}

void RV_from_R0V0(Vec3d &R, Vec3d &V, Vec3d R0, Vec3d V0, double t, double mu){
  /* This function computes the state vector (R, V) from the initial state vector (R0, V0) after elapsed time t.
  mu - gravitational parameter (km^3/s^2)
  R0 - initial position vector (km)
  V0 - initial velocity vector (km/s)
  t - elapsed time (s)
  R - final position vector (km)
  V - final velocity vector (km/s)
  */

  // Magnitudes of R0 and V0
  double r0 = R0.length();
  double v0 = V0.length();
  cout << "r0: " << r0 << ", v0: " << v0 << endl;

  // Initial radial velocity
  double vr0 = (R0*V0)/r0;
  cout << "vr0: " << vr0 << endl;

  // reciprocal of the semimajor axis from the energy eqn
  double alpha = 2/r0 - v0*v0/mu;
  cout << "alpha: " << alpha << endl;

  // get universal anomaly
  double x = kepler_U(t, r0, vr0, alpha);
  cout << "x: " << x << endl;

  // get f and g
  double f, g;
  f_and_g(f, g, x, t, r0, alpha);
  cout << "f: " << f << ", g: " << g << endl;

  // compute final position vector
  R = R0*f + V0*g;

  // magnitude of R
  double r = R.length();

  // get derivatives of f and g
  double fdot, gdot;
  fdot_and_gdot(fdot, gdot, x, r, r0, alpha);

  // compute final velocity
  V = R0*fdot + V0*gdot;

  return;
}

orbitalElements oe_from_sv(Vec3d R, Vec3d V, double mu){
  /* This function computes the classical orbital elements from the state vector (R, V).
  mu - gravitational parameter (km^3/s^2)
  R - position vector in the geocentric equatorial frame (km)
  V - velocity vector in the geocentric equatorial frame (km/s)
  r, v - the magnitudes of R and V
  vr - the raidal velocity component (km/s)
  H - the angular momentum vector (km^2/s)
  h - the magnitude of H (km^2/s)
  incl - inclination of the orbit (radians)
  N - the node line vector (km^2/s)
  n - the magnitude of N
  cp - the cross product of N and R
  RA - the Right Ascension of the ascending node (radians)
  E - the eccentricity vector
  e - the eccentricity (magnitude of E)
  eps - a small number below which the eccentricity is considered to be zero
  w - argument of perigee (radians)
  TA - the true anomaly (radians)
  a - semimajor axis (km)
  oe - the structure of orbital elements, including h, e, RA, incl, w, TA, and a
  */

  double eps = 1.0e-10;
  double r = R.length();
  double v = V.length();
  double vr = (R*V)/r;
  Vec3d H = R/V;
  //cout << "H: " << H << endl;
  double h = H.length();
  double incl = acos(H.z/h);
  Vec3d N = Vec3d(0.0, 0.0, 1.0)/H;
  //cout << "N: " << N << endl;
  double n = N.length();
  double RA;
  if (n!=0){
    RA = acos(N.x/n);
    if (N.y<0) RA = 2*PI - RA;
  }
  else RA = 0.0;
  Vec3d E = (R*(v*v - mu/r) - V*r*vr)/mu;
  //cout << "E: " << E << endl;
  double e = E.length();
  double w;
  if (n!=0){
    if (e>eps){
      w = acos((N*E)/n/e);
      if (E.z<0) w = 2*PI - w;
    }
    else w=0.0;
  }
  else w=0.0;
  double TA;
  Vec3d cp;
  if (e>eps){
    TA = acos((E*R)/e/r);
    if (vr<0) TA = 2*PI - TA;
  }
  else{
    cp = N/R;
    TA = acos((N*R)/n/r);
    if (cp.z < 0) TA = 2*PI - TA;
  }
  double a = h*h/mu/(1.0 - e*e);
  orbitalElements oe;
  oe.h = h;
  oe.e = e;
  oe.RA = RA;
  oe.incl = incl;
  oe.w = w;
  oe.TA = TA;
  oe.a = a;

  return oe;
}

void sv_from_oe(Vec3d& R, Vec3d& V, orbitalElements oe, double mu){
  /* This function computes the state vector (R, V) from the classical orbital elements.
  mu - gravitational parameter (km^3/s^2)
  oe - orbital elements, including h, e, RA, incl, w, and TA
  h - angular momentum (km^2/s)
  e - eccentricity
  RA - right ascension of the ascending node (radians)
  incl - inclination of the orbit (radians)
  w - argument of perigee (radians)
  TA - true anomaly (radians)
  R3_w - rotation matrix about the z-axis through the angle w
  R1_i - rotation matrix about the x-axis through the angle i
  R3_RA - rotation matrix about the z-axis through the angle RA
  Q_px - matrix of the transformation from perifocal to geocentric equatorial frame
  Rp - position vector in the perifocal frame (km)
  Vp - velocity vector in the perifocal frame (km/s)
  R - position vector in the geocentric equatorial frame (km)
  V - velocity vector in the geocentric equatorial frame (km/s)
  */

  double h = oe.h;
  double e = oe.e;
  double RA = oe.RA;
  double incl = oe.incl;
  double w = oe.w;
  double TA = oe.TA;
  Vec3d Rp((Vec3d(0.0, 1.0, 0.0) * sin(TA) + Vec3d(1.0, 0.0, 0.0) * cos(TA)) * h*h/mu / (1.0 + e*cos(TA)));
  Vec3d Vp((Vec3d(1.0, 0.0, 0.0) * -(sin(TA)) + Vec3d(0.0, 1.0, 0.0) * (e + cos(TA))) * mu/h);
  //cout << "Rp: " << Rp << endl;
  //cout << "Vp: " << Vp << endl;
  Matrix3 R3_RA, R1_i, R3_w, Q_px;
  R3_RA.RotateZ(RA);
  //cout << "R3_RA: " << R3_RA << endl;
  R1_i.RotateX(incl);
  //cout << "R1_i: " << R1_i << endl;
  R3_w.RotateZ(w);
  //cout << "R3_w: " << R3_w << endl;
  Q_px = R3_RA.Transpose() * R1_i.Transpose() * R3_w.Transpose();
  //cout << "Q_px: " << Q_px << endl;
  R = Q_px * Rp;
  V = Q_px * Vp;
}

bool gibbs(Vec3d& V2, Vec3d R1, Vec3d R2, Vec3d R3, double mu){

  /* This function uses the Gibbs method of orbit determination to compute the velocity corresponding to the second of three supplied position vectors.
  mu - gravitational parameter (km^3/s^2)
  R1, R2, R3 - three coplanar geocentric position vectors (km)
  r1, r2, r3 - the magnitudes of R1, R2, and R3 (km)
  C12, C23, C31 - three independent cross products of R1, R2, and R3
  N, D, S - vectors formed during the Gibbs procedure
  tol - tolerance for determining if R1, R2, and R3 are coplanar
  err - (0, or false) if R1, R2, and R3 are coplanar, otherwise (1, or true)
  V2 - the velocity corresponding to R2 (km/s)
  */

  double tol = 1.0e-4;
  bool err = false;
  double r1 = R1.length();
  double r2 = R2.length();
  double r3 = R3.length();
  Vec3d C12 = R1/R2;
  Vec3d C23 = R2/R3;
  Vec3d C31 = R3/R1;
  if (abs((R1*C23)/r1/(C23.length())) > tol) err = true;
  Vec3d N = C23*r1 + C31*r2 + C12*r3;
  Vec3d D = C12 + C23 + C31;
  Vec3d S = R1*(r2-r3) + R2*(r3-r1) + R3*(r1-r2);
  V2 = ((D/R2)/r2 + S) * sqrt(mu/N.length()/D.length());
  return err;
}

double y(double z, double r1, double r2, double A){
  // Equation 5.38
  return (r1 + r2 + A*(z*S(z) - 1.0)/sqrt(C(z)));
}

double F(double z, double t, double r1, double r2, double A, double mu){
  // Equation 5.40
  return (pow(y(z, r1, r2, A)/C(z),1.5)*S(z) + A*sqrt(y(z, r1, r2, A)) - sqrt(mu)*t);
}

double dFdz(double z, double r1, double r2, double A){
  // Equation 5.43
  if (z==0) return (sqrt(2.0)/40.0*pow(y(0.0, r1, r2, A),1.5) + A/8.0*(sqrt(y(0.0, r1, r2, A)) + A*sqrt(1.0/2.0/y(0, r1, r2, A))));
  else return (pow((y(z,r1,r2,A)/C(z)), 1.5)*(1.0/2.0/z*(C(z) - 3.0*S(z)/2.0/C(z)) + 3.0*pow(S(z), 2)/4.0/C(z)) + A/8.0*(3.0*S(z)/C(z)*sqrt(y(z, r1, r2, A)) + A*sqrt(C(z)/y(z,r1,r2,A))));
}

void lambert(Vec3d& V1, Vec3d& V2, Vec3d R1, Vec3d R2, double t, string orb, double mu){
  /* This function solves Lambert's problem
  mu - gravitational parameter (km^3/s^2)
  R1, R2 - initial and final position vectors (km)
  r1, r2 - magnitudes of R1 and R2
  t - time of flight from R1 to R2 (s)
  V1, V2 - initial and final velocity vectors (km/s)
  C12 - cross product of R1 into R2
  theta - angle between R1 and R2
  orb - 'prograde' or 'retrograde' (string)
  A - a constant given by Equation 5.35
  z - alpha*x^2, where alpha is the reciprocal of the semimajor axis and x is the universal anomaly
  y(z) - function of z given by Equation 5.38
  F(z,t) - function of z and t from Equation 5.40
  dFdz(z) - derivative of F(z,t) from Equation 5.43
  ratio - F/dFdz
  tol - tolerance on precision of convergence
  nmax - maximum number of iterations of Newton's method
  f, g - Lagrange coefficients
  gdot - time derivative of g
  C(z), S(z) - stumpff functions
  */

  double r1 = R1.length();
  double r2 = R2.length();
  Vec3d C12 = R1/R2;
  double theta = acos((R1*R2)/r1/r2);
  if (orb[0]=='r' || orb[0]=='R'){ // retrograde
    if (C12.z >= 0.0) theta = 2*PI - theta;
  }
  else { // assume prograde
    if (C12.z <= 0.0) theta = 2*PI - theta;
  }
  double A = sin(theta)*sqrt(r1*r2/(1.0 - cos(theta)));
  double z = -100.0; // determine approximately where F(z,t) changes sign, use that z to start Eqn 5.45
  while (y(z, r1, r2, A) < 0.0) z+=0.1; // C++ can't take sqrt of negative number. F(z,t) takes sqrt of y(z).
  while (F(z,t,r1,r2,A,mu) < 0.0) z += 0.1;
  double tol = 1.0e-8;
  int nmax = 5000;
  double ratio = 1.0;
  int n = 0;
  while (abs(ratio) > tol && n <= nmax){
    n++;
    ratio = F(z,t,r1,r2,A,mu)/dFdz(z,r1,r2,A);
    z -= ratio;
  }
  if (n >= nmax) cout << "Number of iterations exceeded " << nmax << "." << endl;
  double f = 1.0 - y(z, r1, r2, A)/r1;
  double g = A*sqrt(y(z, r1, r2, A)/mu);
  double gdot = 1.0 - y(z, r1, r2, A)/r2;
  V1 = (R2 - (R1*f))/g;
  V2 = ((R2*gdot) - R1)/g;
}

double julian(int year, int month, int day, int hour, int min, double sec){
  /* This function returns the julian day for any year between 1900 and 2100.
  j0 - julian day at 0 hr UCT
  year - range: 1901 - 2099
  month - range: 1 - 12
  day - range: 1 - 31
  hour - range: 0 - 23
  min - range: 0 - 59
  sec - range: 0 - 59.9999
  */

  int j = 367*year - (7*(year + ((month + 9)/12))/4) + (275*month/9) + day;
  double j0 = j + 1721013.5;
  double ut = hour + min/60.0 + sec/3600.0;
  double jd = j0 + ut/24.0;
  return jd;
}

double lst(int year, int month, int day, int hour, int min, double sec, int lon_deg, int lon_min, double lon_sec){
  /* This function determines the local siderial time.
  year - range: 1901 - 2099
  month - range: 1 - 12
  day - range: 1 - 31
  hour - range: 0 - 23
  min - range: 0 - 59
  sec - range: 0 - 59.9999
  lon_deg - range: -360 - 360
  lon_min - range: 0 - 59
  lon_sec - range: 0 - 59.9999
  */

  double j0 = julian(year, month, day);
  double j = (j0 - 2451545)/36525;
  double g00 = 100.4606184 + 36000.77004*j + 0.000387933*j*j - 2.583e-8*pow(j,3);
  if (g00 >= 360) g00 = g00 - (static_cast<int>(g00)/360)*360;
  else if (g00 < 360) g00 = g00 - (static_cast<int>(g00)/360 - 1)*360;
  double ut = hour + min/60.0 + sec/3600.0;
  double elon = lon_deg + lon_min/60.0 + lon_sec/3600.0;
  double gst = g00 + 360.98564724*ut/24.0;
  double lst = gst + elon;
  lst -= 360*static_cast<int>(lst/360);
  return lst;
}

void RV_from_observe(Vec3d& R, Vec3d& V, double rho, double rhodot, double A, double Adot, double a, double adot, double theta, double phi, double H, double f, double Re, double wE, double mu){
  /* This function calulates the geocentric equatorial position and velocity vectors of an object from radar observations of range, azimuth, elevation angle and their rates.
  Re - equatorial radius of the earth (km)
  f - earth's flattening factor
  wE - angular velocity of the earth (rad/s)
  omega - earth's angular velocity vector in the geocentric equatorial frame (rad/s)
  theta - local siderial time of tracking site (degrees)
  phi - geodetic latitude of the tracking site (degrees)
  H - elevation of the tracking site (km)
  Rs - geocentric equatorial position vector of the tracking site (km)
  Rsdot- inertial velocity of the tracking site (km/s)
  rho - slant range of the object (km)
  rhodot - range rate (km/s)
  A - azimuth of the object relative to the observation (tracking) site (degrees)
  Adot - time rate of change of the Azimuth (degrees/s)
  a - elevation angle of the object relative to the observation (tracking) site (degrees)
  adot - time rate of change of the elevation angle (deg/s)
  dec - topocentric equatorial declination of the object (rad)
  decdot - declination rate (rad/s)
  h - hour angle of the object (rad)
  RA - topocentric equatorial right ascension of the object (rad)
  RAdot - right ascension rate (rad/s)
  Rho - unit vector from site to object
  Rhodot - time rate of change of Rho (1/s)
  R - geocentric equatorial position vector of the object (km)
  V - geocentric equatorial velocity vector of the object (km)
  */

  Vec3d omega(0.0, 0.0, wE);
  // convert to radians
  A *= deg;
  Adot *= deg;
  a *= deg;
  adot *= deg;
  theta *= deg;
  phi *= deg;
  double x = (Re/sqrt(1.0 - (2.0*f - f*f)*sin(phi)*sin(phi)) + H)*cos(phi)*cos(theta);
  Vec3d Rs((Re/sqrt(1.0 - (2.0*f - f*f)*sin(phi)*sin(phi)) + H)*cos(phi)*cos(theta),
              (Re/sqrt(1.0 - (2.0*f - f*f)*sin(phi)*sin(phi)) + H)*cos(phi)*sin(theta),
              (Re*(1.0 - f)*(1.0 - f)/sqrt(1.0 - (2.0*f - f*f)*sin(phi)*sin(phi)) + H)*sin(phi));
  Vec3d Rsdot = omega/Rs;
  double dec = asin(cos(phi)*cos(A)*cos(a) + sin(phi)*sin(a));
  double h = acos((cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a))/cos(dec));
  if (A > 0.0 && A < PI) h = 2*PI - h;
  double RA = theta - h;
  Vec3d Rho(cos(RA)*cos(dec), sin(RA)*cos(dec), sin(dec));
  R = Rs + Rho*rho;
  double decdot = (-Adot*cos(phi)*sin(A)*cos(a) + adot*(sin(phi)*cos(a) - cos(phi)*cos(A)*sin(a)))/cos(dec);
  double RAdot = wE + (Adot*cos(A)*cos(a) - adot*sin(A)*sin(a) + decdot*sin(A)*cos(a)*tan(dec))/(cos(phi)*sin(a) - sin(phi)*cos(A)*cos(a));
  Vec3d Rhodot(-RAdot*sin(RA)*cos(dec) - decdot*cos(RA)*sin(dec),
                RAdot*cos(RA)*cos(dec) - decdot*sin(RA)*sin(dec), decdot*cos(dec));
  V = Rsdot + Rho*rhodot + Rhodot*rho;
}

void gauss(Vec3d& R, Vec3d& V, Vec3d& R_old, Vec3d& V_old, double dec1, double RA1, double dec2, double RA2, double dec3, double RA3, double altitude, double latitude, double LST1, double LST2, double LST3, double t1, double t2, double t3, double Re, double flat, double mu){
  /* This function uses the Gauss method with iterative improvement to calculate the state vectors of an orbiting body from angles-only observations at three closely spaced time.
  mu - gravitation parameter (km^3/s^2)
  dec1, dec2, dec3 - the declinations of the 3 observations (deg)
  RA1, RA2, RA3 - the Right Ascensions of the 3 observations (deg)
  altitude - the altitude of the observing site (km)
  latitude - the latitude of the observing site (deg)
  LST1, LST2, LST3 - the local siderial times of the observations at t1, t2, and t3 (deg)
  t1, t2, t3 - times of the observations (s)
  tau1, tau2, tau3 - time intervals between observations (s)
  R1, R2, R3 - the observation site position vectors at t1, t2, and t3 (km)
  Rho1, Rho2, Rho3 - the direction cosine vectors of the satellite at t1, t2, and t3
  p1, p2, p3 - cross products among the three direction cosine vectors
  D0 - scalar triple product of Rho1, Rho2, and Rho3
  D - Matrix of the nine scalar triple products of R1, R2, and R3 with p1, p2, and p3
  E - dot product of R2 and Rho2
  A, B - constants in the expression relating slant range to geocentric radius
  a, b, c - coefficients of the 8th order polynomial in the estimated geocentric radius x
  x - positive root of the 8th order polynomial
  rho1, rho2, rho3 - the slant ranges at t1, t2, t3 (km)
  r1, r2, r3 - the position vectors at t1, t2, t3 (km)
  R_old, V_old - the estimated state vectors at the end of algorithm 5.5 (km, km/s)
  rho1_old, rho2_old, rho3_old - the values of the slant ranges at t1, t2, t3 at the beginning of iterative improvement, algorithm 5.6 (km)
  diff1, diff2, diff3 - magnitudes of the differences between old and new slant ranges at the end of each iteration
  tol - the error tolerance determining convergence
  n - number of passes throngh the iterative improvement loop
  nmax - limit on the number of iterations
  ro, vo - magnitude of the position and velocity vectors (km, km/s)
  vro - radial velocity component (km)
  a - reciprocal of the semimajor axis (1/km)
  v2 - computer velocity at time t2 (km/s)
  R, V - the state vectors at the end of iterative improvement algorithm 5.6 (km, km/s)
  */

  // convert to radians
  double phi = latitude*deg;
  dec1 *= deg; dec2 *= deg; dec3 *= deg;
  RA1 *= deg; RA2 *= deg; RA3 *= deg;
  LST1 *= deg; LST2 *= deg; LST3 *= deg;

  double fac1 = Re/sqrt(1.0 - (2.0*flat - flat*flat)*sin(phi)*sin(phi));
  double fac2 = (Re*(1.0 - flat)*(1.0 - flat)/sqrt(1.0 - (2.0*flat - flat*flat)*sin(phi)*sin(phi)) + altitude)*sin(phi);
  Vec3d R1((fac1 + altitude)*cos(phi)*cos(LST1), (fac1 + altitude)*cos(phi)*sin(LST1), fac2);
  Vec3d R2((fac1 + altitude)*cos(phi)*cos(LST2), (fac1 + altitude)*cos(phi)*sin(LST2), fac2);
  Vec3d R3((fac1 + altitude)*cos(phi)*cos(LST3), (fac1 + altitude)*cos(phi)*sin(LST3), fac2);
  Vec3d Rho1(cos(dec1)*cos(RA1), cos(dec1)*sin(RA1), sin(dec1));
  Vec3d Rho2(cos(dec2)*cos(RA2), cos(dec2)*sin(RA2), sin(dec2));
  Vec3d Rho3(cos(dec3)*cos(RA3), cos(dec3)*sin(RA3), sin(dec3));
  double tau1 = t1 - t2;
  double tau3 = t3 - t2;
  double tau = tau3 - tau1;
  Vec3d p1 = Rho2/Rho3;
  Vec3d p2 = Rho1/Rho3;
  Vec3d p3 = Rho1/Rho2;
  double D0 = Rho1*p1;
  Matrix3 D(R1*p1, R1*p2, R1*p3, R2*p1, R2*p2, R2*p3, R3*p1, R3*p2, R3*p3);
  double E = R2*Rho2;
  double A = (-d12*tau3/tau + d22 + d32*tau1/tau)/D0;
  double B = (d12*(tau3*tau3 - tau*tau)*tau3/tau + d32*(tau*tau - tau1*tau1)*tau1/tau)/6/D0;
  double a = -(A*A + 2.0*A*E + R2.length()*R2.length());
  double b = -2.0*mu*B*(A + E);
  double c = -(mu*B*mu*B);
  double tol = 1.0e-8;
  int nmax = 5000;
  double ratio = 1.0;
  int n = 0;
  // Newton's method of finding roots - just start at 1
  double x = 10000;
  while (abs(ratio)>tol && n<=nmax){
    n++;
    ratio = (pow(x,8) + a*pow(x,6) + b*pow(x,3) + c)/(8*pow(x,7) + 6*a*pow(x,5) + 3*b*pow(x,2));
    x -= ratio;
  }
  double f1 = 1.0 - mu*tau1*tau1/pow(x,3)/2.0;
  double f3 = 1.0 - mu*tau3*tau3/pow(x,3)/2.0;
  double g1 = tau1 - mu*pow((tau1/x),3)/6.0;
  double g3 = tau3 - mu*pow((tau3/x),3)/6.0;
  double rho2 = A + mu*B/pow(x,3);
  double rho1 = ((6.0*(d31*tau1/tau3 + d21*tau/tau3)*pow(x,3) + mu*d31*(tau*tau - tau1*tau1)*tau1/tau3)/(6.0*pow(x,3) + mu*(tau*tau - tau3*tau3)) - d11)/D0;
  double rho3 = ((6.0*(d13*tau3/tau1 - d23*tau/tau1)*pow(x,3) + mu*d13*(tau*tau - tau3*tau3)*tau3/tau1)/(6.0*pow(x,3) + mu*(tau*tau - tau3*tau3)) - d33)/D0;
  Vec3d r1 = R1 + Rho1*rho1;
  Vec3d r2 = R2 + Rho2*rho2;
  Vec3d r3 = R3 + Rho3*rho3;
  Vec3d v2 = (r1*-f3 + r3*f1)/(f1*g3 - f3*g1);
  R_old = r2;
  V_old = v2;
  // End of algorithm 5.5, start algorithm 5.6
  double rho1_old = rho1;
  double rho2_old = rho2;
  double rho3_old = rho3;
  double diff1 = 1.0;
  double diff2 = 1.0;
  double diff3 = 1.0;
  n = 0;
  nmax = 1000;
  double ro, vo, vro, x1, x3, ff1, ff3, gg1, gg3, c1, c3;
  while (diff1>tol && diff2>tol && diff3>tol && n <= nmax){
    n++;
    ro = r2.length();
    vo = v2.length();
    vro = (v2*r2)/ro;
    a = 2.0/ro - vo*vo/mu;
    x1 = kepler_U(tau1, ro, vro, a);
    x3 = kepler_U(tau3, ro, vro, a);
    f_and_g(ff1, gg1, x1, tau1, ro, a);
    f_and_g(ff3, gg3, x3, tau3, ro, a);
    f1 = (f1 + ff1)/2.0;
    f3 = (f3 + ff3)/2.0;
    g1 = (g1 + gg1)/2.0;
    g3 = (g3 + gg3)/2.0;
    c1 = g3/(f1*g3 - f3*g1);
    c3 = -g1/(f1*g3 - f3*g1);
    rho1 = (-(d11) + d21/c1 - d31*c3/c1)/D0;
    rho2 = (d12*(-c1) + d22 - d32*c3)/D0;
    rho3 = (d13*(-c1)/c3 + d23/c3 - d33)/D0;
    r1 = R1 + Rho1*rho1;
    r2 = R2 + Rho2*rho2;
    r3 = R3 + Rho3*rho3;
    v2 = (r1*(-f3) + r3*f1)/(f1*g3 - f3*g1);
    diff1 = abs(rho1 - rho1_old);
    diff2 = abs(rho2 - rho2_old);
    diff3 = abs(rho3 - rho3_old);
    rho1_old = rho1;
    rho2_old = rho2;
    rho3_old = rho3;
  } // close while loop
  if (n>nmax) cout << "Number of iterations exceeds " << nmax << "." << endl;
  else cout << "Number of iterations: " << n << endl;
  R = r2;
  V = v2;
}

void planet_data(double J200_oe[], double rates[], int planet_id){
  /* gets planet's J200 orbital element data and centennial rates from table 8.1
  planet_id - 1 to 9
  J200_elements - 9 by 6 matrix of J200 orbital elements for the nine planets. The columns are:
    a - semimajor axis (AU)
    e - eccentricity
    i - inclination (deg)
    RA - right ascension of the ascending node (deg)
    w_hat - longitude of perihelion (deg)
    L - mean longitude (deg)
  cent_rates - 9 by 6 matrix of the time rates of change of the J200_elements per Julian century (Cy). The columns are:
    a_dot - (AU/Cy)
    e_dot - (1/Cy)
    i_dot - (arcseconds/Cy)
    RA_dot - (arcseconds/Cy)
    w_hat_dot - (arcseconds/Cy)
    L_dot - (arcseconds/Cy)
  J200_oe - row of J200_elements corresponding to the planet_id, with AU converted to km
  rates - row of cent_rates corresponding to the planet_id, with AU converted to km, and arcseconds converted to deg
  au - astronomical units (km)
  */

  double J200_elements[9][6] = {
    { 0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084 },
    { 0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973 },
    { 1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435 },
    { 1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332 },
    { 5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438 },
    { 9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432 },
    { 19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218 },
    { 30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003 },
    { 39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881} };
  double cent_rates[9][6] = {
    { 0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29 },
    { 0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06 },
    { -0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63 },
    { -0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78 },
    { 0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35 },
    { -0.00301530, -0.00036762, 6.11, -1591.05, -1948.89, 4401052.95 },
    { 0.00152025, -0.00019150, -2.09, -1681.4, 1312.56, 1542547.79 },
    { -0.00125196, 0.00002514, -3.64, -151.25, -844.43, 786449.21 },
    { -0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90 } };
  for(int i=0; i<6; i++){
    J200_oe[i] = J200_elements[planet_id-1][i];
    rates[i] = cent_rates[planet_id-1][i];
  }
  const double AU = 149597871.0;
  J200_oe[0] *= AU;
  rates[0] *= AU;
  for(int i=2; i<6; i++) rates[i]/=3600.0; // convert arcseconds to deg
}

void planet_elements(planetElements& pe, Vec3d& R, Vec3d& V, double& jd, int planet_id, int year, int month, int day, int hour, int minute, double second, double mu){
  /* This function calculates the orbital elements and state vectors of a planet at a date and time
  mu - gravitational parameter of the sun, 1.327124e11 (km^3/s^2)
  pe - planetary elements, including:
    h - angular momentum (km^2/s)
    e - eccentricity
    RA - right ascension of the ascending node (deg)
    incl - inclination (deg)
    w - argument of perihelion (deg)
    TA - true anomaly (deg)
    a - semimajor axis (km)
    w_hat - longitude of perihelion (deg)
    L - mean longitude (= w_hat + M) (deg)
    M - mean anomaly (deg)
    E - eccentric anomaly
  planet_id - 1 to 9, for Mercury to Pluto
  year - range: 1901 - 2099
  month - range: 1 - 12
  day - range: 1 - 31
  hour - range: 0 - 23
  minute - range: 0 - 59
  second - range: 0.0 - 59.9999
  j0 - Julian day of the date at 0 hr UCT
  ut - universal time in fractions of a day
  jd - Julian day of the date and time
  J200_oe - array of 6 orbital elements for a planet
  rates - array of 6 orbital element time rates of change for a planet
  t0 - Julian centuries between J200 and jd
  elements - orbital elements at jd
  R - heliocentric position vector
  V - heliocentric velocity vector
  */

  jd = julian(year, month, day, hour, minute, second);
  double J200_oe[6];
  double rates[6];
  planet_data(J200_oe, rates, planet_id);
  double t0 = (jd - 2451545.0)/36525.0;
  double elements[6];
  for(int i=0; i<6; i++) elements[i] = J200_oe[i] + rates[i]*t0;
  double a = elements[0];
  double e = elements[1];
  double h = sqrt(mu*a*(1.0 - e*e));
  double incl = elements[2];//zero_to_360(elements[2]);
  double RA = zero_to_360(elements[3]);
  double w_hat = zero_to_360(elements[4]);
  double L = zero_to_360(elements[5]);
  double w = zero_to_360(w_hat - RA);
  double M = zero_to_360(L- w_hat);
  double E = kepler_E(e, M*deg); // convert M to radians
  double TA = zero_to_360(2.0*atan(sqrt((1.0 + e)/(1.0 - e))*tan(E/2.0))/deg);
  pe.oe.h = h;
  pe.oe.e = e;
  pe.oe.RA = RA*deg;//convert to radians
  pe.oe.incl = incl*deg;
  pe.oe.w = w*deg;
  pe.oe.TA = TA*deg;
  pe.oe.a = a;
  pe.w_hat = w_hat;
  pe.L = L;
  pe. M = M;
  pe. E = E/deg;
  orbitalElements oe;
  oe.h = h;
  oe.e = e;
  oe.RA = RA*deg;
  oe.incl = incl*deg;
  oe.w = w*deg;
  oe.TA = TA*deg;
  sv_from_oe(R, V, oe, mu);
}

void interplanetary(Vec3d& Rp1, Vec3d& Vp1, double& jd1, Vec3d& Rp2, Vec3d& Vp2, double& jd2, Vec3d& V1, Vec3d& V2, int depart[], int arrive[], double mu){
  /* This function determines a spacecraft trajectory from the sphere of influence of planet 1 to that of planet 2.
  mu - gravitational parameter of the sun (km^3/s^2)
  Rp1, Vp1 - state vectors of planet 1 at departure (km, km/s)
  Rp2, Vp2 - state vectors of planet 2 at arrival (km, km/s)
  jd1, jd2 - Julian day numbers at departure and arrival
  tof - time of flight (s)
  depart - array of planet_id, year, month, day, hour, minute, second at departure
  arrive - array of planet_id, year, month, day, hour, minute, second at arrival
    planet_id - range: 1 - 9
    year - range: 1901 - 2099
    month - range: 1 - 12
    day - range: 1 - 31
    hour - range: 0 - 23
    minute - range: 0 - 59
    second - range: 0 - 59
  Rp1, V1 - heliocentric state vectors of spacecraft at departure (km, km/s)
  Rp2, V2 - heliocentric state vectors of spacecraft at arrival (km, km/s)
  V1, V2 - the trajectory of the spacecraft
  */

  planet_elements(pe, Rp1, Vp1, jd1, depart[0], depart[1], depart[2], depart[3], depart[4], depart[5], static_cast<double>(depart[6]));
  planet_elements(pe, Rp2, Vp2, jd2, arrive[0], arrive[1], arrive[2], arrive[3], arrive[4], arrive[5], static_cast<double>(arrive[6]));
  double tof = (jd2 - jd1)*24.0*3600.0;
  lambert(V1, V2, Rp1, Rp2, tof, "prograde", MU_SUN);
}

double rocket(double e1, double isp1, double e2, double isp2, double e3, double isp3, double payload, double dV){
  /* This function determines the optimal mass for a multi-stage launch vehicle.
  e1, e2, e3 - structural mass ratios for each stage
  isp1, isp2, isp3 - specific impulse for each stage
  payload - mass of payload
  dV - final velocity of booster
  */

  double const g = .00981; // km/s^2
  double c1 = isp1*g;
  double c2 = isp2*g;
  double c3 = isp3*g;
  // find roots
  int nmax = 5000;
  int n = 0;
  int tol = 1.0e-8;
  double ratio = 1.0;
  double x;
  if (e3>0.0) x = 1.0/(fmin(fmin(c1, c2), c3)) + 0.1;
  else if (e3==0.0 && e2>0.0) x = 1.0/(fmin(c1, c2)) + 0.1;
  else if (e3==0.0 && e2==0.0 & e1>0.0) x = 1.0/c1 + 0.1;
  else return 0;
  cout << "Starting x: " << x << endl;

  while (abs(ratio)>tol && n<=nmax){
    n++;
    if(e3>0.0) ratio = (c1*log(c1*x - 1.0) + c2*log(c2*x - 1.0) + c3*log(c3*x - 1.0) - log(x)*(c1+c2+c3) - (c1*log(c1*e1) + c2*log(c2*e2) + c3*log(c3*e3)) - dV)/(c1*c1/(c1*x - 1.0) + c2*c2/(c2*x - 1.0) + c3*c3/(c3*x - 1.0) - (c1 + c2 + c3)/x);
    else if (e3==0.0 && e2>0.0) ratio = (c1*log(c1*x - 1.0) + c2*log(c2*x - 1.0) - log(x)*(c1+c2) - (c1*log(c1*e1) + c2*log(c2*e2)) - dV)/(c1*c1/(c1*x - 1.0) + c2*c2/(c2*x - 1.0) - (c1 + c2)/x);
    else if (e3==0.0 && e2==0.0 && e1>0.0) ratio = (c1*log(c1*x - 1.0) - log(x)*(c1) - (c1*log(c1*e1)) - dV)/(c1*c1/(c1*x - 1.0) - (c1)/x);
    x -= ratio;
  }
  cout << "Number of iterations: " << n << ", x: " << x << endl;
  cout << "Ratio: " << abs(ratio) << endl;

  // mass ratio of each stage (mE + mp + mPL) / (mE + mPL)
  double n1 = (c1*x - 1.0)/(c1*e1*x);
  double n2 = 0.0;
  if (e2>0.0) n2 = (c2*x - 1.0)/(c2*e2*x);
  double n3 = 0.0;
  if (e3>0.0) n3 = (c3*x - 1.0)/(c3*e3*x);
  cout << "Optimum mass ratios: Stage 1: " << n1 << ", Stage 2: " << n2 << ", Stage 3: " << n3 << endl;

  // step masses of each stage mE + mp
  double m3 = 0.0;
  if (e3>0.0) m3 = (n3 - 1.0)/(1.0 - n3*e3)*payload;
  double m2 = 0.0;
  if (e2>0.0) m2 = (n2 - 1.0)/(1.0 - n2*e2)*(m3 + payload);
  double m1 = (n1 - 1.0)/(1.0 - n1*e1)*(m2 + m3 + payload);
  cout << "Step masses: Stage 1: " << m1 << ", Stage 2: " << m2 << ", Stage 3: " << m3 << endl;

  // empty masses of each stage, mE
  double mE1 = e1*m1;
  double mE2 = 0.0;
  if (e2>0.0) mE2 = e2*m2;
  double mE3 = 0.0;
  if (e3>0.0) mE3 = e3*m3;
  cout << "Empty masses: Stage 1: " << mE1 << ", Stage 2: " << mE2 << ", Stage 3: " << mE3 << endl;

  // propellant masses of each stage, mp
  double mp1 = m1 - mE1;
  double mp2 = m2 - mE2;
  double mp3 = m3 - mE3;
  cout << "Propellant masses: Stage 1: " << mp1 << ", Stage 2: " << mp2 << ", Stage 3: " << mp3 << endl;

  // payload ratios, pr
  double pr1 = (m2 + m3 + payload)/m1;
  double pr2 = 0.0;
  if (e2>0.0) pr2 = (m3 + payload)/m2;
  double pr3 = 0.0;
  if (e3>0.0) pr3 = payload/m3;
  cout << "Payload ratios: Stage 1: " << pr1 << ", Stage 2: " << pr2 << ", Stage 3: " << pr3 << endl;

  // total vehicle mass
  double m0 = m1 + m2 + m3 + payload;
  cout << "Total vehicle mass: " << m0 << endl;
  // overall payload fraction
  double Pf = payload/m0;
  cout << "For the " << payload << " kg payload, Overall payload fraction: " << Pf << endl;

  // make sure its a minimum
  bool min=false;
  if (e3>0.0){
    if (x*c1*(e1*n1 - 1.0)*(e1*n1 - 1.0) + 2.0*e1*n1 - 1.0 > 0.0 &&
      x*c2*(e2*n2 - 1.0)*(e2*n2 - 1.0) + 2.0*e2*n2 - 1.0 > 0.0 &&
      x*c3*(e3*n3 - 1.0)*(e3*n3 - 1.0) + 2.0*e3*n3 - 1.0 > 0.0) min = true;
  }
  else if (e2>0.0){
    if (x*c1*(e1*n1 - 1.0)*(e1*n1 - 1.0) + 2.0*e1*n1 - 1.0 > 0.0 &&
      x*c2*(e2*n2 - 1.0)*(e2*n2 - 1.0) + 2.0*e2*n2 - 1.0 > 0.0) min = true;
  }
  else if (x*c1*(e1*n1 - 1.0)*(e1*n1 - 1.0) + 2.0*e1*n1 - 1.0 > 0.0) min = true;

  if (min) return m0;
  else return 0;
}

double payload(double mE1, double mp1, double isp1, double mE2, double mp2, double isp2, double Vbo){
  /* Determine mass of payload for a burn-out velocity, given structural mass, propellant mass, and specific impulse for each stage.
  mE1, mE2 - structural (empty) mass of stage 1 and 2
  mp1, mp2 - propellant mass for stage 1 and 2
  isp1, isp2 - specific impulse of stage 1 and 2
  Vb0 - final velocity acheived when propellant is entirely burned
  */

  //const double g0 = 0.00981;
  double c1 = isp1*g0;
  double c2 = isp2*g0;
  double m1 = mE1 + mp1;
  double m2 = mE2 + mp2;
  double tol = 1.0e-8;
  int nmax = 5000;
  int n = 0;
  double ratio = 1.0;
  double mPL = 0.0;
  while (abs(ratio)>tol && n<nmax){
    n++;
    ratio = (c1*(log(m1 + m2 + mPL) - log(mE1 + m2 + mPL)) + c2*(log(m2 + mPL) - log(mE2 + mPL)) - Vbo)/(c1/(m1 + m2 + mPL) - c1/(mE1 + m2 + mPL) + c2/(m2 + mPL) - c2/(mE2 + mPL));
    mPL -= ratio;
  }
  cout << "Number of iterations: " << n << ", ratio: " << ratio << ", mPL: " << mPL << endl;
/*
  // now, a different way
  n = 0;
  ratio = 1.0;
  double x = .406; //1.0/fmin(c1, c2) + 0.1;
  double e1 = mE1/m1;
  double e2 = mE2/m2;
  cout << "c1: " << c1 << ", c2: " << c2 << ", e1: " << e1 << ", e2: " << e2 << endl;
  double num, dem;
  while (abs(ratio)>tol && n<nmax){
    n++;
    num = (c1*(log(c1*x - 1.0) - log(c1*e1*x)) + c2*(log(c2*x - 1.0) - log(c2*e2*x)) - Vbo);
    //num = c1*log((c1*x - 1.0)/(c1*e1*x)) + c2*log((c2*x - 1.0)/(c2*e2*x)) - Vbo;
    dem = (c1*c1/(c1*x - 1.0) - c1/x + c2*c2/(c2*x - 1.0) - c2/x);
    //dem = -c1/(x - c1*x*x) - c2/(x - c2*x*x);
    //cout << "num: " << num << ", dem: " << dem << endl;
    ratio = num/dem;
    x -= ratio;
    //cout << "ratio: " << ratio << ", x: " << x << endl;
  }
  cout << "Number of iterations: " << n << ", ratio: " << ratio << ", x: " << x << endl;
  double n1 = (c1*x - 1.0)/(c1*e1*x);
  double n2 = (c2*x - 1.0)/(c2*e2*x);
  double mPL1 = (n1*(mE1+m2) - m1 - m2)/(1.0 - n1);
  double mPL2 = (n2*mE2-m2)/(1.0 - n2);
  cout << "n1 (must be > 1): " << n1 << ", n2 (must be > 1): " << n2 << endl;
  if (n1>1.0) cout << "mPL1: " << mPL1 << endl;
  if (n2>1.0) cout << "mPL2: " << mPL2 << endl;
  cout << "test1 (must be > 0 for min): " << x*c1*(e1*n1 - 1.0)*(e1*n1 - 1.0) + 2.0*e1*n1 - 1.0 << endl;
  cout << "test2 (must be > 0 for min): " << x*c2*(e2*n2 - 1.0)*(e2*n2 - 1.0) + 2.0*e2*n2 - 1.0 << endl;

  // one more time
  n = 0;
  ratio = 1.0;
  x = 1.0/fmin(c1, c2) + 0.1;
  while (abs(ratio)>tol && n<nmax){
    n++;
    num = (c1*log(c1*x - 1.0) + c2*log(c2*x - 1.0) - log(x)*(c1+c2) - c1*log(c1*e1) - c2*log(c2*e2) - Vbo);
    dem = (c1*c1/(c1*x - 1.0) + c2*c2/(c2*x - 1.0) - (c1 + c2)/x);
    ratio = num/dem;
    x -= ratio;
    //cout << "num: " << num << ", dem: " << dem << endl;
    //cout << "ratio: " << ratio << ", x: " << x << endl;
  }
  cout << "Number of iterations: " << n << ", ratio: " << ratio << ", x: " << x << endl;
  n1 = (c1*x - 1.0)/(c1*e1*x);
  n2 = (c2*x - 1.0)/(c2*e2*x);
  mPL1 = (n1*(mE1+m2) - m1 - m2)/(1.0 - n1);
  mPL2 = (n2*mE2-m2)/(1.0 - n2);
  cout << "n1 (must be > 1): " << n1 << ", n2 (must be > 1): " << n2 << endl;
  if (n1>1.0) cout << "mPL1: " << mPL1 << endl;
  if (n2>1.0) cout << "mPL2: " << mPL2 << endl;
  cout << "test1 (must be > 0 for min): " << x*c1*(e1*n1 - 1.0)*(e1*n1 - 1.0) + 2.0*e1*n1 - 1.0 << endl;
  cout << "test2 (must be > 0 for min): " << x*c2*(e2*n2 - 1.0)*(e2*n2 - 1.0) + 2.0*e2*n2 - 1.0 << endl;
*/
  return mPL;
}

double deltaV(double mE1, double mp1, double isp1, double mE2, double mp2, double isp2, double mE3, double mp3, double isp3, double mPL){
  /* This function determines the delta-V of a rocket given its empty stage masses, propellant masses, and specific impulses of each stage.
  mE1, mE2, mE3 - empty structural masses of each stage
  mp1, mp2, mp3 - the propellant masses of each stage
  isp1, isp2, isp3 - the specific impulses of each stage
  */

  double m1 = mE1 + mp1;
  double m2 = mE2 + mp2;
  double m3 = mE3 + mp3;
  double c1 = isp1*g0;
  double c2 = isp2*g0;
  double c3 = isp3*g0;
  double Vbo = 0.0;
  if (mPL==0.0) cout << "Must have a payload mass." << endl;
  else Vbo = c1*(log(m1 + m2 + m3 + mPL) - log(mE1 + m2 + m3 + mPL)) + c2*(log(m2 + m3 + mPL) - log(mE2 + m3 + mPL)) + c3*(log(m3 + mPL) - log(mE3 + mPL));
  return Vbo;
}

