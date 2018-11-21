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
  float adcToMIP=45;
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
      ("channels", po::value<std::vector<int> >()->multitoken()->zero_tokens()->composing(), "channel Ids to calibrate") ;
    
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

  std::cout << "cmap.size() " << cmap.size() << std::endl;
  loadADCToMIP(m_mipFileName, m_module, cmap);

  std::map<int,TGraphErrors*> grhgmap;
  std::map<int,TGraphErrors*> grlgmap;
  std::map<int,TGraphErrors*> grtotmap;
  std::map<int, int> indexmap;
  
  for( unsigned int iski=0; iski<m_skirocs.size(); iski++ ){
    for( unsigned int ichan=0; ichan<m_channels.size(); ichan++ ){
      TGraphErrors *grhg=new TGraphErrors();
      grhg->SetMarkerStyle(20);
      grhg->SetMarkerColor(kBlack);
      grhg->SetLineColor(kBlack);
      grhgmap.insert( std::pair<int,TGraphErrors*>(m_skirocs[iski]*100+m_channels[ichan],grhg) );
      TGraphErrors* grlg=new TGraphErrors();
      grlg->SetMarkerStyle(20);
      grlg->SetMarkerColor(kBlue);
      grlg->SetLineColor(kBlue);
      grlgmap.insert( std::pair<int,TGraphErrors*>(m_skirocs[iski]*100+m_channels[ichan],grlg) );
      TGraphErrors* grtot=new TGraphErrors();
      grtot->SetMarkerStyle(20);
      grtot->SetMarkerColor(kRed);
      grtot->SetLineColor(kRed);
      grtotmap.insert( std::pair<int,TGraphErrors*>(m_skirocs[iski]*100+m_channels[ichan],grtot) );
      indexmap.insert( std::pair<int,int>(iski*100+ichan,int(0)) );
    }
  }

  adcFitter flow(1,1000,0.05);
  adcFitterResult rlow;
  adcFitter fhigh(1,1000,0.05);
  adcFitterResult rhigh;
  setTDRStyle();
  TCanvas *cc=new TCanvas();
  cc->SetWindowSize(600,600);

  TLatex *tex=new TLatex();
  tex->SetTextSize(0.035);
  tex->SetLineColor(kBlack);

  TF1* fhg=new TF1("fhg","pol1",0,300);
  fhg->SetLineColor(kBlack);
  TF1* flg=new TF1("flg","pol1",0,2000);
  flg->SetLineColor(kBlue);

  std::stringstream ss(std::stringstream::in | std::stringstream::out);
  ss << std::setprecision(3) ;

  for( std::map<int,channelData>::iterator it=cmap.begin(); it!=cmap.end(); ++it ){
    int key=it->first;
    channelData chD=it->second;
    fhigh.run(chD.dacInj, chD.high, chD.errorhigh, rhigh, 100);
    float dacToMIP=rhigh.coeff/chD.adcToMIP;
    int index=0;
    for( unsigned int ievt=0; ievt<chD.dacInj.size(); ievt++ ){
      chD.chargeInj.push_back( dacToMIP*chD.dacInj.at(ievt) );
      index++;
      grhgmap[key]->SetPoint(index,chD.chargeInj.at(ievt),chD.high.at(ievt));
      grhgmap[key]->SetPointError(index,0,5);
      grlgmap[key]->SetPoint(index,chD.chargeInj.at(ievt),chD.low.at(ievt));
      grlgmap[key]->SetPointError(index,0,2);
      grtotmap[key]->SetPoint(index,chD.chargeInj.at(ievt),chD.tot.at(ievt));
      grtotmap[key]->SetPointError(index,0,1);
    }    
    fhigh.run(chD.chargeInj, chD.high, chD.errorhigh, rhigh, 100*dacToMIP);
    flow.run(chD.chargeInj, chD.low, chD.errorlow, rlow, 400*dacToMIP);
    fhg->SetParameters(rhigh.pedestal,rhigh.coeff);
    flg->SetParameters(rlow .pedestal,rlow .coeff);

    std::cout << "Coeff lg hg dacToMIP:" << rlow.coeff << "\t" << rhigh.coeff << "\t" << dacToMIP << std::endl;
  
    TH1D* h=new TH1D("","",100,0,1000);
    h->GetXaxis()->SetTitle("Injected charge [MIP]");
    h->GetYaxis()->SetTitle("ADC counts");
    h->GetYaxis()->SetTitleOffset(1.3);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->SetMaximum(3500);
  
    h->DrawCopy("axis");
    grlgmap[key]->Draw("psame");
    flg->Draw("same");
    grhgmap[key]->Draw("psame");
    fhg->Draw("same");
    grtotmap[key]->Draw("psame");      

    TLegend* leg=new TLegend(0.6,0.12,0.9,0.3);
    leg->SetFillStyle(0);
    leg->AddEntry(grhgmap[key],"High gain","ep");
    leg->AddEntry(grlgmap[key],"Low gain","ep");
    leg->AddEntry(grtotmap[key],"TOT","ep");
    leg->Draw();
    cc->Update();
    cc->WaitPrimitive();
    //sleep(2);
    ss.str("");
    ss << "plots/module" << m_module << "/hglgtot_Vs_Inj_Module" << m_module << "_Chip" << chD.skirocID << "_Channel" << chD.channelID << ".pdf";
    cc->SaveAs(ss.str().c_str());
    delete leg;
    delete grhgmap[key];
    delete grlgmap[key];
    delete grtotmap[key];
    delete h;
  }
  delete fhg;
  delete flg;
  delete tex;
  delete cc;
  return 0;
}
