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
#include "adcFitter.h"

struct channelData{
  std::vector<double> dacInj; //DAC unit (aribitrary)
  std::vector<double> chargeInj; //MIP unit (aribitrary)
  std::vector<double> tot;
  std::vector<double> low;
  std::vector<double> high;
  std::vector<double> errorlow;
  std::vector<double> errorhigh;
  int channelID;
  int skirocID;
  float adcToMIP;
};

void loadADCToMIP(std::string mipFileName, int moduleID, std::map<int,channelData> &chanMap)
{
  TFile* file=new TFile(mipFileName.c_str(),"READ");
  if( file->IsOpen() )
    file->Print();
  else
    std::cout << "can not open file " << mipFileName << std::endl;
  TTree* tree=(TTree*)file->Get("mipCalib");
  float MIP;
  int module,chip;
  tree->SetBranchAddress("MIP",&MIP);
  tree->SetBranchAddress("module",&module);
  tree->SetBranchAddress("chip",&chip);
  for( unsigned int ientry(0); ientry<tree->GetEntries(); ++ientry ){
    tree->GetEntry(ientry);
    if( module!=moduleID ) continue;
    for( std::map<int,channelData>::iterator it=chanMap.begin(); it!=chanMap.end(); ++it ){
      if( it->first/100!=chip ) continue;
      it->second.adcToMIP=MIP;
      std::cout << "module chip channel adcToMIP" << moduleID << " " << it->first/100 << " " <<   it->first%100 << " " << it->second.adcToMIP << std::endl;
    }
  }
  delete tree;
  file->Close();
  delete file;
}

double myfunc(double* x, double* par)
{
  return x[0]<par[0] ? 0 : par[1]*x[0]+par[2]-par[3]/(std::pow(x[0],par[4])-par[0]);

}
int main(int argc,char** argv)
{
  int argc1=0;
  char* argv1=(char*)"";
  TApplication* app = new TApplication("toto",&argc1,&argv1);
  std::string m_fileName,m_mipFileName;
  int m_module=0;
  std::vector<int> m_skirocs;
  std::vector<int> m_channels;
  try { 
    /** Define and parse the program options 
     */ 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("fileName", po::value<std::string>(&m_fileName)->default_value("toto.root"), "name of input root file")
      ("mipFileName", po::value<std::string>(&m_mipFileName)->default_value("~/public/data/2018/mipChipCalibration.root"), "name of root file containing ADC to MIP conversion")
      ("module", po::value<int>(&m_module)->default_value(63), "module number")
      ("skirocs", po::value<std::vector<int> >()->multitoken()->zero_tokens()->composing(), "skiroc Ids to calibrate")
      ("channels", po::value<std::vector<int> >()->multitoken()->zero_tokens()->composing(), "channel Ids to calibrate");
    
    po::variables_map vm; 
    try { 
      po::store(po::parse_command_line(argc, argv, desc),  vm); 
      if ( vm.count("help")  ) { 
        std::cout << desc << std::endl; 
        return 0; 
      } 
      po::notify(vm);
    }
    catch(po::error& e) { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << desc << std::endl; 
      return 1; 
    }
    if( vm.count("fileName") ) std::cout << "fileName = " << m_fileName << std::endl;
    if( vm.count("mipFileName") ) std::cout << "mipFileName = " << m_mipFileName << std::endl;
    if( vm.count("module") ) std::cout << "module = " << m_module << std::endl;
    
    if( !vm.count("skirocs") ){
      std::cerr << "skirocs not set -> execute for the 4 skirocs" << std::endl;
      m_skirocs.push_back(0);
      m_skirocs.push_back(1);
      m_skirocs.push_back(2);
      m_skirocs.push_back(3);
    }
    else
      m_skirocs=vm["skirocs"].as<std::vector<int> >();
    
    if( !vm.count("channels") ){
      std::cerr << "channels not set -> execute for channel 0" << std::endl;
      m_channels.push_back(0);
    }
    else
      m_channels=vm["channels"].as<std::vector<int> >();
  }
  catch(std::exception& e) { 
    std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << ", application will now exit" << std::endl; 
    return 2; 
  } 

  TFile *file=new TFile(m_fileName.c_str(),"READ");
  if( file->IsOpen() )
    file->Print();
  else
    std::cout << "can not open file " << m_fileName << std::endl;
  TDirectory *dir=(TDirectory*)file->Get("pulseshapeplotter");
  TTree *tree = (TTree*)dir->Get("tree");
  if (!tree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 0;
  }
  int eventID;
  std::vector<int> *skirocID;
  std::vector<int> *channelID;
  std::vector<float> *HighGainADC;
  std::vector<float> *HighGainTmax;
  std::vector<float> *HighGainErrorADC;
  std::vector<float> *HighGainErrorTmax;
  std::vector<float> *HighGainChi2;
  std::vector<float> *LowGainADC;
  std::vector<float> *LowGainTmax;
  std::vector<float> *LowGainErrorADC;
  std::vector<float> *LowGainErrorTmax;
  std::vector<float> *LowGainChi2;
  std::vector<int> *TotSlow;
  std::vector<int> *ToaRise;
  std::vector<int> *ToaFall;
  skirocID = 0;
  channelID = 0;
  HighGainADC = 0;
  HighGainTmax = 0;
  HighGainChi2 = 0;
  HighGainErrorADC = 0;
  HighGainErrorTmax = 0;
  // HighGainStatus = 0;
  // HighGainNCalls = 0;
  LowGainADC = 0;
  LowGainTmax = 0;
  LowGainChi2 = 0;
  LowGainErrorADC = 0;
  LowGainErrorTmax = 0;
  // LowGainStatus = 0;
  // LowGainNCalls = 0;
  TotSlow = 0;
  ToaRise = 0;
  ToaFall = 0;

  tree->SetBranchAddress("eventID", &eventID);
  tree->SetBranchAddress("skirocID", &skirocID);
  tree->SetBranchAddress("channelID", &channelID);
  tree->SetBranchAddress("HighGainADC",&HighGainADC);
  tree->SetBranchAddress("HighGainTmax",&HighGainTmax);
  tree->SetBranchAddress("HighGainErrorADC",&HighGainErrorADC);
  tree->SetBranchAddress("HighGainErrorTmax",&HighGainErrorTmax);
  tree->SetBranchAddress("HighGainChi2",&HighGainChi2);
  tree->SetBranchAddress("LowGainADC",&LowGainADC);
  tree->SetBranchAddress("LowGainTmax",&LowGainTmax);
  tree->SetBranchAddress("LowGainErrorADC",&LowGainErrorADC);
  tree->SetBranchAddress("LowGainErrorTmax",&LowGainErrorTmax);
  tree->SetBranchAddress("LowGainChi2",&LowGainChi2);
  tree->SetBranchAddress("TotSlow",&TotSlow);
  tree->SetBranchAddress("ToaRise",&ToaRise);
  tree->SetBranchAddress("ToaFall",&ToaFall);
  std::map<int,channelData> cmap;
  for( unsigned int ientry(0); ientry<tree->GetEntries(); ++ientry ){
    tree->GetEntry(ientry);
    if( eventID==0 ) continue;
    for(size_t ihit=0;ihit<HighGainADC->size(); ihit++){
      if( std::find(m_channels.begin(),m_channels.end(), channelID->at(ihit))==m_channels.end() ) continue;
      if( std::find(m_skirocs.begin(),m_skirocs.end(), skirocID->at(ihit))==m_skirocs.end() ) continue;
      int key=100*skirocID->at(ihit)+channelID->at(ihit);
      if( cmap.find(key)==cmap.end() ){
	channelData chData;
	chData.channelID=channelID->at(ihit);
	chData.skirocID=skirocID->at(ihit);
	cmap.insert(std::pair<int,channelData>(key,chData));
      }
      cmap[key].dacInj.push_back(int(eventID*4.095));
      cmap[key].tot.push_back(TotSlow->at(ihit));
      cmap[key].high.push_back(HighGainADC->at(ihit));
      cmap[key].low.push_back(LowGainADC->at(ihit));
      cmap[key].errorhigh.push_back(HighGainErrorADC->at(ihit));
      cmap[key].errorlow.push_back(LowGainErrorADC->at(ihit));
    }
  }
  file->Close();
  delete file;

  loadADCToMIP(m_mipFileName, m_module, cmap);

  std::map<int,TGraphErrors*> grtotmap;
  std::map<int, int> indexmap;
  for( unsigned int iski=0; iski<m_skirocs.size(); iski++ ){
    for( unsigned int ichan=0; ichan<m_channels.size(); ichan++ ){
      TGraphErrors* grtot=new TGraphErrors();
      grtot->SetMarkerStyle(20);
      grtot->SetMarkerColor(kBlack);
      grtot->SetLineColor(kBlack);
      grtotmap.insert( std::pair<int,TGraphErrors*>(m_skirocs[iski]*100+m_channels[ichan],grtot) );
      indexmap.insert( std::pair<int,int>(iski*100+ichan,int(0)) );
    }
  }

  setTDRStyle();
  TCanvas *cc=new TCanvas();
  cc->SetWindowSize(600,600);

  TLatex *tex=new TLatex();
  tex->SetTextSize(0.035);
  tex->SetLineColor(kBlack);

  TF1* f=new TF1("f",myfunc,0,4500,5);
  f->SetLineColor(kBlue);

  adcFitter fhigh(-1,1000,0.05);
  adcFitterResult rhigh;
  totFitter fitter(1,1000,2,0.03);
  totFitterResult fitresult;
  for( std::map<int,channelData>::iterator it=cmap.begin(); it!=cmap.end(); ++it ){
    int key=it->first;
    channelData chD=it->second;
    fhigh.run(chD.dacInj, chD.high, chD.errorhigh, rhigh, 100);
    float dacToMIP=rhigh.coeff/chD.adcToMIP;
    int index=0;
    for( unsigned int ievt=0; ievt<chD.dacInj.size(); ievt++ ){
      chD.chargeInj.push_back( dacToMIP*chD.dacInj.at(ievt) );
      index++;
      grtotmap[key]->SetPoint(index,chD.chargeInj.at(ievt),chD.tot.at(ievt));
      grtotmap[key]->SetPointError(index,0,1);
    }    
    fitter.run(chD.chargeInj, chD.tot, fitresult);
    f->SetParameters(fitresult.thr,fitresult.coeff,fitresult.pedestal,fitresult.c,fitresult.power);
    f->Draw("same");
    
    // for(unsigned int ievt=0; ievt<chD.tot.size(); ievt++){
    //   std::cout << "event = " << ievt << "\t tot = " << chD.tot.at(ievt) << "\t nmip = " << fitresult.findNumberOfMipsFromToT( chD.tot.at(ievt) ) << std::endl;
    // }
  
    TH1D* h=new TH1D("","",100,0,5000*dacToMIP);
    h->GetXaxis()->SetTitle("Injected charge [MIP]");
    h->GetYaxis()->SetTitle("ADC counts");
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->SetMaximum(1800);
  
    h->DrawCopy("axis");
    grtotmap[key]->Draw("psame");      
    f->Draw("same");

    tex->DrawLatex(200,1600,"TOT = a.x + b - #frac{c}{x^{#alpha}-t} if x>t");
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    ss << std::setprecision(3) ;  
    ss.str("");
    ss << "a = " << fitresult.coeff << " +- " << fitresult.ecoeff;
    tex->DrawLatex(200,1500,ss.str().c_str());

    ss.str("");
    ss << "b = " << fitresult.pedestal << " +- " << fitresult.epedestal;
    tex->DrawLatex(200,1400,ss.str().c_str());

    ss.str("");
    ss << "c = " << fitresult.c << " +- " << fitresult.ec;
    tex->DrawLatex(200,1300,ss.str().c_str());

    ss.str("");
    ss << "t = " << fitresult.thr << " +- " << fitresult.ethr;
    tex->DrawLatex(200,1200,ss.str().c_str());

    ss.str("");
    ss << "#alpha = " << fitresult.power << " +- " << fitresult.epower;
    tex->DrawLatex(200,1100,ss.str().c_str());
    cc->Update();
    cc->WaitPrimitive();
    ss.str("");
    ss << "plots/module" << m_module << "/TOT_Vs_Inj_Module" << m_module << "_Chip" << chD.skirocID << "_Channel" << chD.channelID << ".pdf";
    cc->SaveAs(ss.str().c_str());

    delete grtotmap[key];
    delete h;
  }
  delete cc;
  return 0;
}
