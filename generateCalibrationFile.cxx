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

struct channelCalib{
  int module;
  int chip;
  int channel;
  float hgcoeff,lgcoeff;
  float hgsat,lgsat;
  float totcoeff,totped,totthr,totnorm,totpower;
  float hgedm,lgedm,totedm;
  float hgchi2,lgchi2,totchi2;

  bool operator==(const channelCalib& c) const
  {
    return (channel == c.channel && module==c.module && chip==c.chip);
  } 

  void print()
  {
    std::cout << module << " " << chip << " " << channel << " " << hgcoeff << " " << lgcoeff << " " << totcoeff << " " << totchi2*totedm << std::endl;
  }  
};

class SortChannelCalib
{
  public : 
  SortChannelCalib(){;}
  ~SortChannelCalib(){;}
  static bool sortByTOT(channelCalib a, channelCalib b){ return a.totcoeff < b.totcoeff; }
  static bool sortByHG(channelCalib a, channelCalib b){ return a.hgcoeff < b.hgcoeff; }
  static bool sortByLG(channelCalib a, channelCalib b){ return a.lgcoeff < b.lgcoeff; }
};

channelCalib getMedian(std::vector<channelCalib> vec,int gain/*0 high gain,1 low gain, 2 tot*/,float &sigma)
{
  if(gain==0) std::sort( vec.begin(), vec.end(), SortChannelCalib::sortByHG );
  else if(gain==1) std::sort( vec.begin(), vec.end(), SortChannelCalib::sortByLG );
  else if(gain==2) std::sort( vec.begin(), vec.end(), SortChannelCalib::sortByTOT );
  else{ std::cout<<"Wrong options for getMedian -> exit"<<gain<<std::endl; exit(1);}
  int sigma_1_Index = int( 0.16*(vec.size()-1) );
  int sigma_3_Index = int( 0.84*(vec.size()-1) );
  if(gain==0) sigma=0.5*( vec[sigma_3_Index].hgcoeff-vec[sigma_1_Index].hgcoeff );
  else if(gain==1) sigma=0.5*( vec[sigma_3_Index].lgcoeff-vec[sigma_1_Index].lgcoeff );
  else if(gain==2) sigma=0.5*( vec[sigma_3_Index].totcoeff-vec[sigma_1_Index].totcoeff );
  return vec[0.5*(vec.size()-1)];
}

inline bool exists_file (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc,char** argv)
{
  std::string m_fileName;
  try { 
    /** Define and parse the program options 
     */ 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("fileName", po::value<std::string>(&m_fileName)->default_value("calibrations/Calibration.root"), "name of input root file");
    
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
  TTree *tree = (TTree*)file->Get("calib");
  if (!tree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 0;
  }

  int _module,_chip,_channel;
  double _hgcoeff,_lgcoeff;
  int _hgsat,_lgsat;
  double _totthr,_totcoeff,_totped,_totnorm,_totpower;
  double _hgchi2,_hgedm,_lgchi2,_lgedm,_totchi2,_totedm;
  int _hgncalls,_lgncalls,_totncalls;
  tree->SetBranchAddress("module",&_module);
  tree->SetBranchAddress("chip",&_chip);
  tree->SetBranchAddress("channel",&_channel);
  tree->SetBranchAddress("hgcoeff",&_hgcoeff);
  tree->SetBranchAddress("hgsat",&_hgsat);
  tree->SetBranchAddress("hgchi2",&_hgchi2);
  tree->SetBranchAddress("hgedm",&_hgedm);
  tree->SetBranchAddress("hgncalls",&_hgncalls);
  tree->SetBranchAddress("lgcoeff",&_lgcoeff);
  tree->SetBranchAddress("lgsat",&_lgsat);
  tree->SetBranchAddress("lgchi2",&_lgchi2);
  tree->SetBranchAddress("lgedm",&_lgedm);
  tree->SetBranchAddress("lgncalls",&_lgncalls);
  tree->SetBranchAddress("totthr",&_totthr);
  tree->SetBranchAddress("totcoeff",&_totcoeff);
  tree->SetBranchAddress("totped",&_totped);
  tree->SetBranchAddress("totnorm",&_totnorm);
  tree->SetBranchAddress("totpower",&_totpower);
  tree->SetBranchAddress("totchi2",&_totchi2);
  tree->SetBranchAddress("totedm",&_totedm);
  tree->SetBranchAddress("totncalls",&_totncalls);
  
  std::map< int,std::vector<channelCalib> > calibmap;
  
  for( unsigned int ientry(0); ientry<tree->GetEntries(); ++ientry ){
    tree->GetEntry(ientry);
    std::map< int,std::vector<channelCalib> >::iterator it=calibmap.find(_module*10+_chip);
    if( it==calibmap.end() ){
      channelCalib c;
      c.module=_module;
      c.chip=_chip;
      c.channel=_channel;
      c.hgcoeff=_hgcoeff;
      c.hgsat=_hgsat;
      c.hgchi2=_hgchi2;
      c.hgedm=_hgedm;
      c.lgcoeff=_lgcoeff;
      c.lgsat=_lgsat;
      c.lgchi2=_lgchi2;
      c.lgedm=_lgedm;
      c.totthr=_totthr;
      c.totcoeff=_totcoeff;
      c.totped=_totped;
      c.totnorm=_totnorm;
      c.totpower=_totpower;
      c.totchi2=_totchi2;
      c.totedm=_totedm;
      std::vector<channelCalib> vec; vec.push_back(c);
      calibmap.insert( std::pair< int,std::vector<channelCalib> >(_module*10+_chip,vec) );
    }
    else{
      channelCalib c;
      c.module=_module;
      c.chip=_chip;
      c.channel=_channel;
      std::vector<channelCalib>::iterator jt = std::find(it->second.begin(), it->second.end(), c);
      if( jt==it->second.end() ){
	c.hgcoeff=_hgcoeff;
	c.hgsat=_hgsat;
	c.hgchi2=_hgchi2;
	c.hgedm=_hgedm;
	c.lgcoeff=_lgcoeff;
	c.lgsat=_lgsat;
	c.lgchi2=_lgchi2;
	c.lgedm=_lgedm;
	c.totthr=_totthr;
	c.totcoeff=_totcoeff;
	c.totped=_totped;
	c.totnorm=_totnorm;
	c.totpower=_totpower;
	c.totchi2=_totchi2;
	c.totedm=_totedm;
	it->second.push_back(c);
      }
      else{
	(*jt).hgcoeff=_hgcoeff;
	(*jt).hgsat=_hgsat;
	(*jt).hgchi2=_hgchi2;
	(*jt).hgedm=_hgedm;
	(*jt).lgcoeff=_lgcoeff;
	(*jt).lgsat=_lgsat;
	(*jt).lgchi2=_lgchi2;
	(*jt).lgedm=_lgedm;
	(*jt).totthr=_totthr;
	(*jt).totcoeff=_totcoeff;
	(*jt).totped=_totped;
	(*jt).totnorm=_totnorm;
	(*jt).totpower=_totpower;
	(*jt).totchi2=_totchi2;
	(*jt).totedm=_totedm;
      }
    }
  }  
  file->Close();
  delete file;
  
  std::ostringstream os(std::ostringstream::ate);
  os.str("");
  os << "./calibrations/ADCCalibration.root";
  TFile* outfile=new TFile(os.str().c_str(),"RECREATE");
  TTree* outtree=new TTree("calib","hgcal calibrations tree");
  outtree->Branch("module",&_module);
  outtree->Branch("chip",&_chip);
  outtree->Branch("channel",&_channel);
  outtree->Branch("hgcoeff",&_hgcoeff);
  outtree->Branch("hgsat",&_hgsat);
  outtree->Branch("lgcoeff",&_lgcoeff);
  outtree->Branch("lgsat",&_lgsat);
  outtree->Branch("totthr",&_totthr);
  outtree->Branch("totcoeff",&_totcoeff);
  outtree->Branch("totped",&_totped);
  outtree->Branch("totnorm",&_totnorm);
  outtree->Branch("totpower",&_totpower);

  for(std::map< int,std::vector<channelCalib> >::iterator it=calibmap.begin(); it!=calibmap.end(); ++it){
    std::vector<channelCalib> vec=it->second;
    float sigmahg(0),sigmalg(0),sigmatot(0);
    channelCalib cmedhg=getMedian(vec,0,sigmahg);
    channelCalib cmedlg=getMedian(vec,1,sigmalg);
    channelCalib cmedtot=getMedian(vec,2,sigmatot);
    for(std::vector<channelCalib>::iterator jt=vec.begin(); jt!=vec.end(); ++jt){
      (*jt).print();
      if( (*jt).hgedm*(*jt).hgchi2<1e7 && (*jt).hgcoeff>20 && fabs( (*jt).hgcoeff-cmedhg.hgcoeff )<3*sigmahg ){
	_hgcoeff=(*jt).hgcoeff;
	_hgsat=(*jt).hgsat;
      }
      else{
	_hgcoeff=cmedhg.hgcoeff;
	_hgsat=cmedhg.hgsat;
      }
      if( (*jt).lgedm*(*jt).lgchi2<1e7  && (*jt).lgcoeff>2 && fabs( (*jt).lgcoeff-cmedlg.lgcoeff )<3*sigmalg ){
	_lgcoeff=(*jt).lgcoeff;
	_lgsat=(*jt).lgsat;
      }
      else{
	_lgcoeff=cmedlg.lgcoeff;
	_lgsat=cmedlg.lgsat;
      }
      if( (*jt).totedm*(*jt).totchi2<1e7 && (*jt).totcoeff>0.2 && (*jt).totthr<200 && fabs( (*jt).totcoeff-cmedtot.totcoeff )<3*sigmatot ){
	_totcoeff=(*jt).totcoeff;
	_totped=(*jt).totped;
	_totthr=(*jt).totthr;
	_totpower=(*jt).totpower;	
	_totnorm=(*jt).totnorm;	
      }
      else{
	_totcoeff=cmedtot.totcoeff;
	_totped=cmedtot.totped;
	_totthr=cmedtot.totthr;
	_totpower=cmedtot.totpower;	
	_totnorm=cmedtot.totnorm;	
      }
      _module=(*jt).module;
      _chip=(*jt).chip;
      _channel=(*jt).channel;
      outtree->Fill();
    }
    std::cout << "CHPI MEDIAN : " << _module << " " << _chip << "\t" 
	      << cmedhg.hgcoeff << "+-" << sigmahg << "\t"
	      << cmedlg.lgcoeff << "+-" << sigmalg << "\t"
      	      << cmedtot.totcoeff << "+-" << sigmatot << std::endl;
    std::cout << std::endl;

  }
  outfile->Write();
  outfile->Close();
  return 0;
}
