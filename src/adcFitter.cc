#include <adcFitter.h>

#include <limits>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

static int _maxDAC=1000;
static std::vector<int> _adc;
static std::vector<int> _adcError;
static std::vector<int> _dac;
double adcShape_fcn(int x, double coeff, double pedestal )
{
  return coeff*x+pedestal;
}
double adcShape_chi2(const double *x)
{
  double sum = 0.0;
  for(size_t i=0; i<1000; i++){
    if( _dac[i]<_maxDAC && _adc[i]>0 &&_adc[i]<4000 ){
      double zero = _adc[i]-adcShape_fcn( _dac[i],x[0],x[1]);
      sum += zero * zero / _adcError[i] / _adcError[i];
    }
  }
  return sum;
}

adcFitter::adcFitter(int printLevel, int maxIterations, double toleranceAtTP)
{							     
  m_printLevel=printLevel;
  m_maxIterations=maxIterations;
  m_toleranceAtTP=toleranceAtTP;
}

void adcFitter::run(std::vector<double> &dac, std::vector<double> &adc, std::vector<double> &adcError, adcFitterResult &fit, int maxDAC)
{
  if( dac.size()!=adc.size() ){
    std::cout << "ERROR : we should have the same vector size in void adcFitter::run(std::vector<double> &dac, std::vector<double> &adc, adcFitterResult &fit) -> return without fitting" << std::endl;
    return;
  }
  _maxDAC=maxDAC;
  _dac.clear();
  _adc.clear();
  _adcError.clear();
  for( uint16_t i=0; i<dac.size(); i++ ){
    _dac.push_back(dac[i]);
    _adc.push_back(adc[i]);
    adcError[i]>1?_adcError.push_back(adcError[i]):_adcError.push_back(1);
  }

  ROOT::Math::Minimizer* m = ROOT::Math::Factory::CreateMinimizer("Minuit2", "MIGRAD");
  m->SetMaxFunctionCalls(m_maxIterations);
  m->SetMaxIterations(m_maxIterations);
  m->SetTolerance(0.01);
  m->SetPrintLevel(m_printLevel);
  ROOT::Math::Functor f(&adcShape_chi2, 2);

  m->SetFunction(f);

  m->Clear(); // just a precaution


  m->SetVariable(0, "coeff", 1, 0.1);
  m->SetVariableLimits(0,0,1000);
  m->SetVariable(1, "pedestal", 0, 0);
  m->SetVariableLimits(1,0,0);

  m->Minimize();
  
  const double *xm = m->X();
  const double *errors = m->Errors();
  fit.coeff=xm[0];
  fit.pedestal=xm[1];
  fit.ecoeff=errors[0];
  fit.epedestal=errors[1];
  fit.chi2=m->MinValue();
  fit.edm=m->Edm();
  fit.status=m->Status();
  fit.ncalls=m->NCalls();

  bool findlimit=false;
  for( int ievt=0; ievt<_dac.size(); ievt++ ){
    if(ievt<50) continue;
    double delta=fit.coeff*_dac[ievt]+fit.pedestal-_adc[ievt];
    if( delta/_adc[ievt]>m_toleranceAtTP ){
      int count=0;
      for( int jevt=ievt+1; jevt<ievt+10; jevt++ ){
	delta=fit.coeff*_dac[jevt]+fit.pedestal-_adc[jevt];
	if( delta/_adc[jevt]>m_toleranceAtTP ) count++;
      }
      if(count<3) continue;
      delta=fit.coeff*_dac[ievt]+fit.pedestal-_adc[ievt];
      findlimit=true;
      fit.transitionPoint=dac[ievt];
      break;
    }
  }
  if(!findlimit)
    fit.transitionPoint=4500;
}
