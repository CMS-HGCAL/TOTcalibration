#ifndef TOTFITTER
#define TOTFITTER

#include <vector>
#include <cmath>
#include <iostream>

struct totFitterResult{
totFitterResult() : thr(0.), coeff(0.), pedestal(0.), c(0.), power(0.),
    ethr(0.), ecoeff(0.), epedestal(0.), ec(0.), epower(0.), status(-1), transitionPoint(4500) {;}

  double thr;
  double coeff;
  double pedestal;
  double c;
  double power;
  double chi2;
  double edm;
  double ethr;
  double ecoeff;
  double epedestal;
  double ec;
  double epower;
  double echi2;
  int status;
  int ncalls;
  int transitionPoint;
  
  double totShape(float nmip)
  {
    return nmip<thr ? 0 : coeff*nmip+pedestal-c/(std::pow(nmip,power)-thr);
  }
  float findNumberOfMipsFromToT(int tot)
  {
    float nmin(0),nmax(2000);
    int index;
    float diff(2000);
    float centralValue(0);
    float totThr=totShape(thr+thr/10);
    if(tot<totThr) return -1;
    while(fabs(diff)>1e-2&&index<10000){
      centralValue=nmin+fabs(nmax-nmin)/2;
      float totApprox=totShape(centralValue);
      diff=tot-totApprox;
      if( diff>0 )
	nmin=centralValue;
      else
	nmax=centralValue;
      index++;
    }
    return centralValue;
  }
};


class totFitter{
 public:
  totFitter( int printLevel=1, int maxIterations=100, int totNoise=5, double toleranceAtTP=0.05 );
  ~totFitter(){;}
  void run(std::vector<double> &dac, std::vector<double> &tot, totFitterResult &fit);
 private:
  int m_printLevel;
  int m_maxIterations;
  double m_toleranceAtTP;
};

#endif
