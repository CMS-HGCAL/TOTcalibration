#ifndef ADCFITTER
#define ADCFITTER

#include <vector>
#include <cmath>
#include <iostream>

struct adcFitterResult{
adcFitterResult() : coeff(0.), pedestal(0.), ecoeff(0.), epedestal(0.), status(-1), transitionPoint(4500) {;}
  double coeff;
  double pedestal;
  double chi2;
  double edm;
  double ecoeff;
  double epedestal;
  int status;
  int ncalls;
  int transitionPoint;
};


class adcFitter{
 public:
  adcFitter(int printLevel=1, int maxIterations=100, double toleranceAtTP=0.05 );
  ~adcFitter(){;}
  void run(std::vector<double> &dac, std::vector<double> &adc, std::vector<double> &adcError, adcFitterResult &fit, int maxDAC=1000);
 private:
  int m_printLevel;
  int m_maxIterations;
  double m_toleranceAtTP;
};

#endif
