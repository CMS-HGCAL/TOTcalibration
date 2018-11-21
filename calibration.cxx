#include <iostream>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <sys/stat.h>
#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>

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

inline bool exists_file (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc,char** argv)
{
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

  std::ostringstream os(std::ostringstream::ate);
  os.str("");
  os << "./calibrations/Calibration_Module" << m_module << ".root";
  TFile* outfile;
  TTree* outtree;
  int _module,_chip,_channel;
  double _hgcoeff,_lgcoeff;
  int _hgsat,_lgsat;
  double _totthr,_totcoeff,_totped,_totnorm,_totpower;
  double _hgchi2,_hgedm,_lgchi2,_lgedm,_totchi2,_totedm;
  int _hgncalls,_lgncalls,_totncalls;
  if( !exists_file ( os.str() ) ){
    outfile=new TFile(os.str().c_str(),"RECREATE");
    outtree=new TTree("calib","hgcal calibrations tree");
    outtree->Branch("module",&_module);
    outtree->Branch("chip",&_chip);
    outtree->Branch("channel",&_channel);
    outtree->Branch("hgcoeff",&_hgcoeff);
    outtree->Branch("hgsat",&_hgsat);
    outtree->Branch("hgchi2",&_hgchi2);
    outtree->Branch("hgedm",&_hgedm);
    outtree->Branch("hgncalls",&_hgncalls);
    outtree->Branch("lgcoeff",&_lgcoeff);
    outtree->Branch("lgsat",&_lgsat);
    outtree->Branch("lgchi2",&_lgchi2);
    outtree->Branch("lgedm",&_lgedm);
    outtree->Branch("lgncalls",&_lgncalls);
    outtree->Branch("totthr",&_totthr);
    outtree->Branch("totcoeff",&_totcoeff);
    outtree->Branch("totped",&_totped);
    outtree->Branch("totnorm",&_totnorm);
    outtree->Branch("totpower",&_totpower);
    outtree->Branch("totchi2",&_totchi2);
    outtree->Branch("totedm",&_totedm);
    outtree->Branch("totncalls",&_totncalls);
  }
  else{
    outfile=new TFile(os.str().c_str(),"UPDATE");
    outtree=(TTree*)outfile->Get("calib");
    outtree->SetBranchAddress("module",&_module);
    outtree->SetBranchAddress("chip",&_chip);
    outtree->SetBranchAddress("channel",&_channel);
    outtree->SetBranchAddress("hgcoeff",&_hgcoeff);
    outtree->SetBranchAddress("hgsat",&_hgsat);
    outtree->SetBranchAddress("hgchi2",&_hgchi2);
    outtree->SetBranchAddress("hgedm",&_hgedm);
    outtree->SetBranchAddress("hgncalls",&_hgncalls);
    outtree->SetBranchAddress("lgcoeff",&_lgcoeff);
    outtree->SetBranchAddress("lgsat",&_lgsat);
    outtree->SetBranchAddress("lgchi2",&_lgchi2);
    outtree->SetBranchAddress("lgedm",&_lgedm);
    outtree->SetBranchAddress("lgncalls",&_lgncalls);
    outtree->SetBranchAddress("totthr",&_totthr);
    outtree->SetBranchAddress("totcoeff",&_totcoeff);
    outtree->SetBranchAddress("totped",&_totped);
    outtree->SetBranchAddress("totnorm",&_totnorm);
    outtree->SetBranchAddress("totpower",&_totpower);
    outtree->SetBranchAddress("totchi2",&_totchi2);
    outtree->SetBranchAddress("totedm",&_totedm);
    outtree->SetBranchAddress("totncalls",&_totncalls);
  }

  totFitter ftot(0,1000,2,0.03);
  totFitterResult rtot;
  adcFitter flow(0,1000,0.05);
  adcFitterResult rlow;
  adcFitter fhigh(0,1000,0.05);
  adcFitterResult rhigh;
  for( std::map<int,channelData>::iterator it=cmap.begin(); it!=cmap.end(); ++it ){
    channelData chD=it->second;
    fhigh.run(chD.dacInj, chD.high, chD.errorhigh, rhigh, 100);
    float dacToMIP=rhigh.coeff/chD.adcToMIP;
    for( unsigned int ievt=0; ievt<chD.dacInj.size(); ievt++ )
      chD.chargeInj.push_back( dacToMIP*chD.dacInj.at(ievt) );

    fhigh.run(chD.chargeInj, chD.high, chD.errorhigh, rhigh, 100*dacToMIP);
    flow.run(chD.chargeInj, chD.low, chD.errorlow, rlow, 400*dacToMIP);
    ftot.run(chD.chargeInj, chD.tot, rtot);
    
    std::cout << "Coeff lg hg dacToMIP:" << rlow.coeff << "\t" << rhigh.coeff << "\t" << dacToMIP << std::endl;

    _module=m_module;
    _chip=chD.skirocID;
    _channel=chD.channelID;
    _hgcoeff=rhigh.coeff;
    _hgsat=rhigh.transitionPoint;
    _lgcoeff=rlow.coeff;
    _lgsat=rlow.transitionPoint;
    _totthr=rtot.thr;
    _totcoeff=rtot.coeff;
    _totped=rtot.pedestal;
    _totnorm=rtot.c;
    _totpower=rtot.power;

    _hgchi2=rhigh.chi2;
    _hgedm=rhigh.edm;
    _hgncalls=rhigh.ncalls;
    _lgchi2=rlow.chi2;
    _lgedm=rlow.edm;
    _lgncalls=rlow.ncalls;
    _totchi2=rtot.chi2;
    _totedm=rtot.edm;
    _totncalls=rtot.ncalls;

    std::cout << "ToT chi2 = " << _totchi2 << std::endl;

    outtree->Fill();
  }
  outtree->Write();
  outfile->Close();
  return 0;
}
