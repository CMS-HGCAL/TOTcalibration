#include <iostream>
#include <unistd.h>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TLatex.h>
#include <TApplication.h>

#include "tdrstyle.h"
#include "totFitter.h"

double myfunc(double* x, double* par)
{
  return x[0]<par[0] ? 0 : par[1]*x[0]+par[2]-par[3]/(std::pow(x[0],par[4])-par[0]);

}
int main(int argc,char** argv)
{
  int argc1=0;
  char* argv1=(char*)"";
  TApplication* app = new TApplication("toto",&argc1,&argv1);

  TH1D* h=new TH1D("","",1000,0,4500);
  h->GetXaxis()->SetTitle("Arbitrary unit of injected charge [DAC]");
  h->GetYaxis()->SetTitle("TOT slow");
  h->GetYaxis()->SetTitleOffset(1.3);
  h->SetLineWidth(2);
  //  h->SetMaximum(400);
  TF1 *f=new TF1("f",myfunc,0,4500,5);
  f->SetParameter(0,362);
  f->SetParLimits(0,350,450);
  f->SetParameter(1,0.0337);
  f->SetParameter(2,209);
  f->SetParameter(3,1.66e6);
  f->SetParameter(4,1.5);
  h->FillRandom("f",5000000);

  std::vector<double> dac,tot;
  TGraphErrors *gr=new TGraphErrors();
  gr->SetMarkerStyle(7);
  for( int bin=0; bin<h->GetNbinsX(); bin++ ){
    dac.push_back( h->GetBinCenter(bin) );
    tot.push_back( h->GetBinContent(bin) );
    gr->SetPoint(bin,h->GetBinCenter(bin),h->GetBinContent(bin));
    gr->SetPointError(bin,0.1,5);
  }
  totFitter fitter(1,10000,100);
  totFitterResult fitresult;
  fitter.run(dac, tot, fitresult);
  
  setTDRStyle();
  TCanvas *cc=new TCanvas();
  cc->SetWindowSize(600,600);
  h->Draw("axis");
  gr->Fit(f);
  gr->Draw("psame");
  f->Draw("same");
  cc->Update();
  cc->WaitPrimitive();

  return 0;
}
