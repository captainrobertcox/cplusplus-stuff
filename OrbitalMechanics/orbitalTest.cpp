#include <iostream>
#include <iomanip>
#include <string>
#include "orbital.h"

using namespace std;

int main(){
//Orbital o;


  double E, F, x;
  orbitalElements oe;
  planetElements pe;

  E = kepler_E(0.37255, 3.6029);
  cout << endl << "eccentricity 0.37255, mean anomaly 3.6029 = eccentric anomaly 3.47942" << endl;
  cout << "Eccentric anomaly (radians): " << E << endl;

  F = kepler_H(2.7696, 40.69);
  cout << endl << "eccentricity 2.7696, hyperbolic mean anomaly 40.69 = hyperbolic eccentric anomaly 3.46309" << endl;
  cout << "Hyperbolic eccentric anomaly (dimensionless): " << F << endl;

  x = kepler_U(3600.0, 10000.0, 3.0752, 1.0 / -19655.0);
  cout << endl << "elapsed time 3600, radius 10000, radial vel 3.0752, semimajor -19655 = Universal anomaly 128.511" << endl;
  cout << "Universal anomaly (km^0.5): " << x << endl;

  Vec3d R0(7000.0, -12124.0, 0.0);
  Vec3d V0(2.6679, 4.6210, 0.0);

  Vec3d R,V;

  RV_from_R0V0(R, V, R0, V0, 3600.0);
  cout << endl << "init pos (7000.0, -12124.0, 0.0), init vel (2.6679, 4.6210, 0.0), after 3600 s" << endl;
  cout << "r = (-3297.77, 7413.4, 0), v = (-8.2976, -0.964045, 0)" << endl;
  cout << "Position Vector (km): " << R << ", Velocity Vector (km/s): " << V << endl;
  cout << endl << "init pos (7000.0, -12124.0, 0.0), init vel (2.6679, 4.6210, 0.0), after 0 s" << endl;
  RV_from_R0V0(R, V, R0, V0, 0.0);
  cout << "Position Vector (km): " << R << ", Velocity Vector (km/s): " << V << endl;

  R.set(-6045.0, -3490.0, 2500.0);
  V.set(-3.457, 6.618, 2.533);
  oe = oe_from_sv(R,V);
  cout << endl << "r(-6045.0, -3490.0, 2500.0), v(-3.457, 6.618, 2.533)" << endl;
  cout << "h = 58311, e = .171212, RA = 255.279, incl 153.249, arg prg 20.0683, TA 28.4456, semi 8788.1, T 8198.96" << endl;
  displayOE(oe);

  oe.h = 80000.0;
  oe.e = 1.4;
  oe.RA = 40.0*deg;
  oe.incl = 30.0*deg;
  oe.w = 60.0*deg;
  oe.TA = 30.0*deg;
  sv_from_oe(R, V, oe);
  cout << endl << "h 80000, e 1.4, RA 40, incl 30, w 60, TA 30 = r(-4039.9, 4814.56, 3628.62), v(-10.386, -4.77192, 1.74388)" << endl;
  cout << "State Vector R (km): " << R << endl;
  cout << "State Vector V (km/s): " << V << endl;
  oe = oe_from_sv(R,V);
  displayOE(oe);
  //for (int i = 0; i < 100000; i++){
    //sv_from_oe(R, V, oe);
    //oe = oe_from_sv(R, V);
  //}
  //cout << endl << "After 100000 cycles:" << endl;
  //displayOE(oe);

  Vec3d R1(-294.32, 4265.1, 5986.7);
  Vec3d R2(-1365.4, 3637.6, 6346.8);
  Vec3d R3(-2940.3, 2473.7, 6555.8);
  Vec3d V2;
  bool err;
  err = gibbs(V2, R1, R2, R3);
  cout << endl << "gibbs: r1(-294.32, 4265.1, 5986.7), r2(-1365.4, 3637.6, 6346.8), r3(-2940.3, 2473.7, 6555.8) = v2(-6.2176, -4.01237, 1.59915)" << endl;
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
  cout << endl << "lambert: r1(5000.0, 10000.0, 2100.0), r2(-14600.0, 2500.0, 7000.0), dt=3600 s, prograde" << endl;
  cout << "v1(-5.99249 1.92536 3.24564), v2(-3.31246 -4.19662 -0.385288)" << endl;
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

  cout << fixed << showpoint << endl << "Julian day for 5/12/2004 14:45:30.0 (2453138.115): " << julian(2004,5,12,14,45,30.0) << noshowpoint << endl;
  cout.unsetf(ios_base::fixed);

  double LST = lst(2004,3,3,4,30,0.0,139,47,0.0);
  cout << endl << "LST for 3/3/2004 4:30:00, longitude 139 deg, 47' 0\" (8.57688): " << LST << " deg, or (0.571792) " << LST/15.0 << " hr." << endl;

  RV_from_observe(R0, V0, 2551.0, 0.0, 90.0, 0.1130, 30.0, 0.05651, 300.0, 60.0, 0.0);
  oe = oe_from_sv(R0, V0);
  cout << endl << "observed:" << endl;
  cout << "range 2551km, range rate 0km/s, azimuth 90deg, az rate .113deg/s, elevation 30deg, elev rt .05651deg/s, lst 300deg, latitude 60, altitude 0 km:" << endl;
  cout << "r(3830.68, -2216.47, 6605.09), v(1.50357, -4.56099, -0.291536)" << endl;
  cout << "State Radius Vector (km): " << R0 << endl;
  cout << "State Velocity Vector (km/s): " << V0 << endl;
  displayOE(oe);

  // given lst's
  cout << endl << "gauss:" << endl;
  cout << "observation  time      Right Ascension   Declination   lst" << endl;
  cout << "1            0         43.5365           -8.7833       44.5065" << endl;
  cout << "2            118.1     54.4196           -12.0739      45.0000" << endl;
  cout << "3            237.6     64.3178           -15.1054      45.4992" << endl;
  cout << "latitude 40 deg, altitude 1 km" << endl;
  gauss(R1, V1, R0, V0, -8.78334, 43.5365, -12.0739, 54.4196, -15.1054, 64.3178, 1.0, 40, 44.5065, 45.0, 45.4992, 0.0, 118.104, 237.577);
  cout << "with iterative improvement: r(5662.04, 6537.95, 3269.05), v(-3.88542, 5.12141, -2.2434)" << endl;
  cout << "Without iterative improvement:" << endl << "Rold: " << R0 << ", Vold: " << V0 << endl;
  oe = oe_from_sv(R0, V0);
  displayOE(oe);
  cout << "With iterative improvement:" << endl << "R: " << R1 << ", V: " << V1 << endl;
  oe = oe_from_sv(R1, V1);
  displayOE(oe);

  double jd;
  planet_elements(pe, R0, V0, jd, 3, 2003, 8, 27, 12, 0, 0.0);
  cout << endl << "For Earth, 8-27-2003 at noon: r(1.35589e+08 -6.68029e+07 286.909) v(12.6804 26.61 -0.000212731)" << endl;
  cout << "and rMagnitude = 1.51152e+08, vMagnitude = 29.4769" << endl;
  cout << "Planet: " << planets[3-1] << endl;
  cout << "Julian date: " << showpoint << fixed << jd << noshowpoint << endl;
  cout.unsetf(ios_base::fixed);
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
  cout << endl << "Depart Earth 11-7-1996 at midnight, arrive Mars 9-12-1997 at midnight (309 days)" << endl;
  cout << "Departure spacecraft velocity magnitude:  32.7427, Arrival spacecraft vel mag: 22.1637" << endl;
  cout << "Arrival: planet vel(25.0386 -0.220288 -0.620623), spacecraft vel(22.1581 -0.19666 -0.457847)" << endl;
  cout << "Depart " << planets[depart[0]-1] << " on Julian day: " << showpoint << fixed << jd1 << noshowpoint << endl;
  cout.unsetf(ios_base::fixed);
  cout << "Planet position: " << Rp1 << ", Magnitude: " << Rp1.length() << endl;
  cout << "Planet velocity: " << Vp1 << ", Magnitude: " << Vp1.length() << endl;
  cout << "Spacecraft velocity: " << V1 << ", Magnitude: " << V1.length() << ", V-inf at departure: " << Vinf1 << ", Magnitude: " << Vinf1.length() << endl;
  cout << "Time of flight (days): " << tof << endl;
  cout << "Arrive " << planets[arrive[0]-1] << " on Julian day: " << showpoint << fixed << jd2 << noshowpoint << endl;
  cout.unsetf(ios_base::fixed);
  cout << "Planet position: " << Rp2 << ", Magnitude: " << Rp2.length() << endl;
  cout << "Planet velocity: " << Vp2 << ", Magnitude: " << Vp2.length() << endl;
  cout << "Spacecraft velocity: " << V2 << ", Magnitude: " << V2.length() << ", V-inf at arrival: " << Vinf2 << ", Magnitude: " << Vinf2.length() << endl;
  oe = oe_from_sv(Rp1, V1, MU_SUN);
  cout << "Orbital elements of trajectory at departure:" << endl;
  displayOE(oe, MU_SUN);
  oe = oe_from_sv(Rp2, V2, MU_SUN);
  cout << "True anomaly at arrival (deg): " << oe.TA/deg << endl << endl;

  double m0;
  m0 = rocket(0.10, 400.0, 0.15, 350.0, 0.20, 300.0, 5000, 10);
  cout << endl << "rocket stage 1 mass ratio .1, isp 400, stage 2 mass ratio .15, isp 350, stage 3 mass ratio .2, isp 300, payload 5000, dv 10: " << endl;
  cout << "Total rocket mass (191,200): " << m0 << endl;
  m0 = rocket(0.2, 300.0, 0.3, 235.0, 0.0, 0.0, 10.0, 6.2);
  cout << endl << "rocket stage 1 mass ratio .2, isp 300, stage 2 mass ratio .3, isp 235, payload 10kg, dV 6.2km/s:" << endl;
  cout << "Total rocket mass (1125): " << m0 << endl;
  //m0 = rocket(0.142857, 290.0, 0.047619, 450.0, 0.0, 0.0, 114000, 9.686);
  //cout << "Total rocket mass (2,424,000): " << m0 << endl;
  double mPL = payload(150000.0, 900000.0, 290.0, 60000.0, 1200000.0, 415.0, 8.094 + 2.0 - KSC_Rotation_Vel);
  cout << endl << "stage 1: 2 solids, each 525000 total mass, 450000 propellant mass, isp 290. stage 2: 2 liquids, each empty mass 30000, propellant mass 600000, isp 450. 300 km orbit, 2km/s drag loss, from KSC:" << endl;
  cout << "Payload mass (114,000): " << mPL << endl;

  double Vbo = deltaV(150000.0, 900000.0, 290.0, 60000.0, 1200000.0, 450.0, 0.0, 0.0, 0.0, 114000.0);
  cout << endl << "stage 1: empty mass 150000, propellant mass 900000, isp 290. stage 2: empty mass 60000, propellant mass 1200000, isp 450. payload 114000." << endl;
  cout << "Delta-V at burnout: " << Vbo << endl;

  cout << endl << "Escape velocity from Earth's surface: " << escapeV() << endl;
  cout << "Orbital velocity at 300 km: " << orbitalV(300.0) << endl;
  cout << "Delta-V to acheive 300 km orbit: " << acheiveOrbitDeltaV(300.0) << endl;
  cout << "dV to LEO (300 km): " << dV_LEO(300.0) << endl;
  cout << "Delta-V to acheive 300 km orbit (+ 2km/s for drag, gravity, and buffer): " << acheiveOrbitDeltaV(300.0) + 2.0 << endl;
  cout << "Delta-V to acheive 300 km orbit (+ 20% for drag, gravity, and buffer): " << acheiveOrbitDeltaV(300.0) * 1.2 << endl;
  cout << "Delta-V to 300 km orbit, heading East from KSC (latitude 28.5729): " << acheiveOrbitDeltaV(300.0, KSC_Latitude) << endl;
  cout << "dV to LEO from KSC, 300 km orbit, heading East: " << dV_LEO_from_KSC(300.0) << endl;
  cout << "Delta-V to 300 km orbit, heading East from KSC (latitude 28.5729) (+2km/s loss...): " << acheiveOrbitDeltaV(300.0, KSC_Latitude) + 2.0 << endl;
  cout << "Delta-V to 300 km orbit, heading East from Equator (latitude 0.001) (+2km/s loss...): " << acheiveOrbitDeltaV(300.0, 0.001) + 2.0 << endl;
  cout << "Delta-V to 200 km orbit, heading 45 azimuth from KSC (+2km/s loss): " << acheiveOrbitDeltaV(200.0, KSC_Latitude, 45.0) + 2.0 << endl;
  cout << endl << "R: (0.0, -cos(28.5729)*Re, sin(28.5729)*Re), V: (0.408456, 0.0, 0.0)" << endl;
  R.set(0.0, -cos(28.5729*3.1415927/180.0)*6378, sin(28.5729*3.1415927/180.0)*6378);
  V.set(0.408456, 0.0, 0.0);
  oe = oe_from_sv(R, V);
  displayOE(oe);

  return 0;
}
