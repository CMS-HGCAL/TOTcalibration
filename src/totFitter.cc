#include <totFitter.h>

#include <limits>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

int _totNoise=5.;
int _thr=4500;
static std::vector<int> _tot;
static std::vector<int> _dac;
// double totShape_fcn(int x, double thr, double coeff, double pedestal, double c, double p)
// {
//   return x<thr ? 0 : coeff*x+pedestal-c/(std::pow(x,p)-thr);
// }
double totShape_fcn(int x, double coeff, double pedestal, double c, double p)
{
  return x<_thr ? 0 : coeff*x+pedestal-c/(std::pow(x,p)-_thr);
}
double totShape_chi2(const double *x)
{
  double sum = 0.0;
  for(size_t i=0; i<_tot.size(); i++){
    if( _tot[i]>10 && _dac[i]<800 ){
      double zero = _tot[i]-totShape_fcn( _dac[i],x[0],x[1],x[2],x[3]);
      sum += zero * zero / _totNoise / _totNoise;
    }
  }
  return sum;
}

totFitter::totFitter( int printLevel, int maxIterations, int noise, double toleranceAtTP )
{							     
  m_printLevel=printLevel;
  m_maxIterations=maxIterations;
  m_toleranceAtTP=toleranceAtTP;
  _totNoise=noise;
}

void totFitter::run(std::vector<double> &dac, std::vector<double> &tot, totFitterResult &fit)
{
  _thr=4500;
  _tot.clear();
  _dac.clear();
  for( uint16_t i=0; i<dac.size(); i++ ){
    _dac.push_back(dac[i]);
    _tot.push_back(tot[i]);
    if(tot[i]>4&&_thr==4500){
      bool totStart=true;
      for( uint16_t j=i+1; j<i+5; j++ )
	if(tot[j]<5)
	  totStart=false;
      if( totStart )
	_thr=_dac[i];
    }
  }
  ROOT::Math::Minimizer* m = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  m->SetMaxFunctionCalls(m_maxIterations);
  m->SetMaxIterations(m_maxIterations);
  m->SetTolerance(0.1);
  m->SetPrintLevel(m_printLevel);
  ROOT::Math::Functor f(&totShape_chi2, 4);

  m->SetFunction(f);

  m->Clear(); // just a precaution


  m->SetVariable(0, "coeff", 1.3, 0.01);
  m->SetVariableLimits(0,0,2);
  m->SetVariable(1, "pedestal", 150, 1);
  m->SetVariableLimits(1,0,1000);
  m->SetVariable(2, "c", 7e5, 100);
  m->SetVariableLimits(2,0,1e10);
  m->SetVariable(3, "power", 2, 0.01);
  m->SetVariableLimits(3,0.,4);

  m->Minimize();
  
  const double *xm = m->X();
  const double *errors = m->Errors();
  fit.thr=_thr;
  fit.ethr=1;
  fit.coeff=xm[0];
  fit.pedestal=xm[1];
  fit.c=xm[2];
  fit.power=xm[3];
  fit.ecoeff=errors[0];
  fit.epedestal=errors[1];
  fit.ec=errors[2];
  fit.epower=errors[3];
  fit.chi2=m->MinValue();
  fit.edm=m->Edm();
  fit.status=m->Status();
  fit.ncalls=m->NCalls();

  bool findlimit=false;
  for( size_t ievt=0; ievt<_tot.size(); ievt++ ){
    double delta=fit.coeff*_dac[ievt]+fit.pedestal-_tot[ievt];
    if( delta/_tot[ievt]<m_toleranceAtTP ){
      int count=0;
      for( size_t jevt=ievt+1; jevt<ievt+10; jevt++ ){
	delta=fit.coeff*_dac[jevt]+fit.pedestal-_tot[jevt];
	if( delta/_tot[jevt]<m_toleranceAtTP ) count++;
      }
      if(count<3) continue;
      findlimit=true;
      fit.transitionPoint=dac[ievt];
      break;
    }
  }
  if(!findlimit)
    fit.transitionPoint=4500;
}
