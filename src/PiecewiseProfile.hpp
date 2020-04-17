/// \file PiecewiseProfile.hpp
/// \author khermans 
///
/// \brief
///
/// file to create piecewise profiles and read in data from tracer files

#ifndef PIECEWISE_PROFILE_HPP_
#define PIECEWISE_PROFILE_HPP_

#include <math.h>
#include <memory> 
#include <algorithm> 
#include <iostream> 
#include <string> 
#include <fstream> 
#include <sstream> 
#include <iomanip>

#include "Utilities/FunctionVsTime.hpp" 
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class PiecewiseProfile : public FunctionVsTime<double> { 
public: 
  PiecewiseProfile(const std::vector<double>& time, 
      const std::vector<double>& vals) : mTimes(time), mVals(vals){}
  
  double operator()(const double time) const {
    return interp(time, mTimes, mVals); 
  }
   
  double ReverseInterp(const double val) const { 
    return interp(val, mVals, mTimes); 
  } 
   
  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const {
    return std::unique_ptr<FunctionVsTime<double>>(new PiecewiseProfile(mTimes, mVals));
  }

protected: 
  std::vector<double> mTimes, mVals; 
  
  double interp(const double xin, const std::vector<double>& x, 
      const std::vector<double>& y) const {
    auto lbound = x.begin(); 
    if (x[0] > x[x.back()]) {
      lbound = std::lower_bound(x.begin(), x.end(), xin, 
          [](double a, double b){return a>b;}); 
    } else { 
      lbound = std::lower_bound(x.begin(), x.end(), xin, 
          [](double a, double b){return a<b;}); 
    } 
    int idx = lbound - x.begin();
    if (lbound == x.begin()) return y[0];  
    if (lbound == x.end()) return y.back();  
    double h = (xin - x[idx-1])/(x[idx] - x[idx-1]); 
    return y[idx]*h + y[idx-1]*(1.0-h); 
  }
}; 

// read in density from file
PiecewiseProfile ReadDensityProfileFromFile(std::string fname) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, rhos;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
  
    // Read the data 
    while (std::getline(f, line)) {
        double tc, T9c, rhoc;
        std::stringstream ss(line); 
        ss >> tc >> T9c >> rhoc;
        times.push_back(tc); 
        rhos.push_back(rhoc); 
    //    std::cout << "rhoc: " << rhoc << "\n"; 
    }
    return PiecewiseProfile(times, rhos); 
} 

// read in temperatures from file in T9
PiecewiseProfile ReadTempProfileFromFile(std::string fname) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, temps;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    while (std::getline(f, line)) {
        double tc, T9c, rhoc;
        std::stringstream ss(line); 
        ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc;
        times.push_back(tc); 
        temps.push_back(T9c/1000000000.0);
    //    std::cout << "tc: " << tc << "\tT9c: " << T9c << "\n"; 
    }
    return PiecewiseProfile(times, temps); 
} 

// read ye profile from file
PiecewiseProfile ReadYeProfileFromFile(std::string fname) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, yes;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    while (std::getline(f, line)) {
        double tc, T9c, rhoc, rc, yec;
        std::stringstream ss(line); 
        ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> rc >> yec;
        times.push_back(tc); 
        yes.push_back(yec);
    //    std::cout << "tc: " << tc << "\tT9c: " << T9c << "\n"; 
    }
    return PiecewiseProfile(times, yes); 
} 

// find the start time for files that go to NSE
double FindStartTFromFile(std::string fname) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    int i = 0;
    int start_NSE = 0;
    int start_t_index;
    double max_T9c = 0;
    while (std::getline(f, line)) {
        double tc, T9c, rhoc, rc, yec;
        std::stringstream ss(line); 
        ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> rc >> yec;
        times.push_back(tc); 
    //    std::cout << T9c << std::endl;
        if (T9c > 10000000000.0 )
            start_NSE = i;
        if (start_NSE>0 & T9c < 10000000000.0 )
        {
            start_t_index = i-1;
            break;
        }
        else if (T9c > max_T9c)
        {
            max_T9c = T9c; 
            start_t_index = i;
        }
        i++;
  //    std::cout << "tc: " << tc << "\tT9c: " << T9c << "\n"; 
    }
    return start_t_index; 
} 

// read in density 
std::vector<double> ReadDensityFromFile(std::string fname) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, rho;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    while (std::getline(f, line)) {
        double tc, T9c, rhoc, radc, yec, ENUEc, ENUBc, LNUEc, LNUBc;
        std::stringstream ss(line); 
        ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> radc >> yec >> ENUEc >> ENUBc >> LNUEc >> LNUBc;
        times.push_back(tc); 
        rho.push_back(rhoc);
    }
    return rho; 
} 

// read in radii 
std::vector<double> ReadRadiiFromFile(std::string fname, int t_start_index) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, radii;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    int line_no = 1;
    while (std::getline(f, line)) {
        if ( line_no > t_start_index ){
            double tc, T9c, rhoc, radc, yec, ENUEc, ENUBc, LNUEc, LNUBc;
            std::stringstream ss(line); 
            ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> radc >> yec >> ENUEc >> ENUBc >> LNUEc >> LNUBc;
            times.push_back(tc); 
            radii.push_back(radc);
        } else {
            line_no++;
        }
    }
    return radii; 
} 

// read in times 
std::vector<double> ReadTimesFromFile(std::string fname, int t_start_index) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, temp_ENU;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    int line_no = 1;
    while (std::getline(f, line)) {
        if ( line_no > t_start_index ){
            double tc, T9c, rhoc, radc, yec, ENUEc, ENUBc, LNUEc, LNUBc;
            std::stringstream ss(line); 
            ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> radc >> yec >> ENUEc >> ENUBc >> LNUEc >> LNUBc;
            times.push_back(tc); 
        } else {
            line_no++;
        }
    }
    return times; 
} 

// read in neutrino average energies
std::vector<std::vector<double>> ReadENUFromFile(std::string fname, int t_start_index) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, temp_ENU;  
    std::vector<std::vector<double>> ENU;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    double ConstToGK = 1.0 / 3.15137 / Constants::BoltzmannConstantInMeVPerGK ;
    int line_no = 1;
    while (std::getline(f, line)) {
        if ( line_no > t_start_index ){
            double tc, T9c, rhoc, radc, yec, ENUEc, ENUBc, LNUEc, LNUBc;
            std::stringstream ss(line); 
            ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> radc >> yec >> ENUEc >> ENUBc >> LNUEc >> LNUBc;
            times.push_back(tc); 
            double temp_ENUE =  ENUEc * ConstToGK ;
            double temp_ENUB = ENUBc * ConstToGK ;
            double mult = 1.0;
            std::vector<double> temp_ENU = { temp_ENUE * mult , temp_ENUB * mult };
            ENU.push_back(temp_ENU);
 //           std::cout << "ENUEc: " << temp_ENUE << "\tENUBc: " << temp_ENUB << "\n"; 
        } else {
            line_no++;
        }
    }
    return ENU; 
} 

// read in neutrino luminosities 
std::vector<std::vector<double>> ReadLNUFromFile(std::string fname, int t_start_index) { 
    std::ifstream f(fname); 
    std::string line;
    std::vector<double> times, temp_LNU;  
    std::vector<std::vector<double>> LNU;  
    // Read the first two header lines 
    std::getline(f, line);
    std::getline(f, line);
    // Read the data 
    int line_no = 1;
    while (std::getline(f, line)) {
        if ( line_no > t_start_index ){
            double tc, T9c, rhoc, radc, yec, ENUEc, ENUBc, LNUEc, LNUBc;
            std::stringstream ss(line); 
            ss >> std::setprecision(std::numeric_limits<double>::digits10) >> tc >> T9c >> rhoc >> radc >> yec >> ENUEc >> ENUBc >> LNUEc >> LNUBc;
            times.push_back(tc); 
            double area = 4 * Constants::Pi * radc * radc ;
            double mag = 1e51;
            std::vector<double> temp_LNU = { LNUEc * area * mag , LNUBc * area * mag };
//            std::cout << "LNUEc: " << LNUEc * area << "\tLNUBc: " << LNUBc * area << "\n"; 
            LNU.push_back(temp_LNU);
        } else {
            line_no++;
        }
    }
    return LNU; 
} 

#endif // PIECEWISE_PROFILE_HPP_
