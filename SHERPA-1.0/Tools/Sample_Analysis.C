#include "Sample_Analysis.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"


#include "Final_Selector.H"
#include "Primitive_Detector.H"
#include "Jet_Observables.H"
#include "One_Particle_Observables.H"
#include "Two_Particle_Observables.H"
#include "Four_Particle_Observables.H"

//#include <ctype.h>

using namespace SHERPA;
using namespace ANALYSIS;
using namespace ATOOLS;

// Output treatment: 1 master (ME+PS ... -> einlesbar) -> n specials (ME, PS, ...)

Observable_Data::Observable_Data(std::string _type) : type(_type) {}

void Observable_Data::Output() {
  Flavour flav1,flav2;
  if (Specify()==1) {
    flav1 = Flavour(kf::code(abs(ints[0])));
    if (ints[0]<0) flav1=flav1.Bar();
    msg.Out()<<"Obs : "<<type<<" for "<<flav1<<" : "
	     <<"Range : "<<numbers[0]<<" ... "<<numbers[1]<<" in "<<ints[1]<<" bins,"
	     <<" extra : "<<keywords[0]<<std::endl;
  }
  if (Specify()==2) {
    flav1 = Flavour(kf::code(abs(ints[0])));
    if (ints[0]<0) flav1=flav1.Bar();
    flav2 = Flavour(kf::code(abs(ints[1])));
    if (ints[1]<0) flav2=flav2.Bar();
    msg.Out()<<"Obs : "<<type<<" for "<<flav1<<" "<<flav2<<" : "
	     <<"Range : "<<numbers[0]<<" ... "<<numbers[1]<<" in "<<ints[2]<<" bins,"
	     <<" extra : "<<keywords[0]<<std::endl;
  }
  if (Specify()==10) {
    flav1 = Flavour(kf::code(abs(ints[0])));
    if (ints[0]<0) flav1=flav1.Bar();
    msg.Out()<<"Obs : "<<type<<"("<<ints[2]<<", min :"<<ints[3]<<" jets, max : "<<ints[4]<<"jets), "
	     <<"Range : "<<numbers[0]<<" ... "<<numbers[1]<<" in "<<ints[1]<<" bins,"
	     <<" extra : "<<keywords[0]<<std::endl;
  }
}

int Observable_Data::Specify() {
  if (type==std::string("ET") || type==std::string("PT") ||
      type==std::string("Eta") || type==std::string("E") ||
      type==std::string("EVis"))                                 return 1;
  if (type==std::string("Mass") || type==std::string("PT2") ||
      type==std::string("Eta2"))                                 return 2;
  if (type==std::string("JetPT") || type==std::string("JetEta") || 
      type==std::string("JetE") || type==std::string("JetDR")) return 10;
  return -1;
}



Sample_Analysis::Sample_Analysis(std::ifstream * readin, std::string _phase, const std::string & prefix) :
  m_phase(_phase), m_outputpath(std::string("./")+_phase), m_prefix(prefix), p_analysis(NULL)
{  
  m_type = std::string("Perturbative");
  std::cout<<"Initialize new Sample_Analysis for "<<_phase<<std::endl;
  std::string phasemode;
  int  mode  = ANALYSIS::fill_all|ANALYSIS::splitt_jetseeds;
  bool split = false;
  while (_phase.length()>0) {
    if (_phase[0]==' ' || _phase[0]=='+') _phase = _phase.substr(1);
    else { 
      phasemode = _phase.substr(0,_phase.find(std::string("+")));
      if (phasemode==std::string("ME")) {
	mode  = mode|ANALYSIS::do_me;
	if (split) mode = mode|ANALYSIS::splitt_phase;
	else split = true;
	_phase = _phase.substr(2);
      }
      if (phasemode==std::string("Showers")) {
	mode = mode|ANALYSIS::do_shower;
	if (split) mode = mode|ANALYSIS::splitt_phase;
	else split = true;
	_phase = _phase.substr(7);
      }
      if (phasemode==std::string("Hadrons"))  {
	mode = mode|ANALYSIS::do_hadron;
	if (split) mode = mode|ANALYSIS::splitt_phase;
	else split = true;
	_phase = _phase.substr(7);
	m_type = std::string("Hadronization");
      }
    }
  }
  p_analysis            = new Primitive_Analysis(m_phase,mode);
  Final_Selector * fsel = new Final_Selector("FinalState","Analysed",(rpa.gen.Beam1()==Flavour(kf::e)));
  ReadInFinalSelectors(readin,fsel);
  p_analysis->AddObservable(fsel);
  ReadInObservables(readin);
  SetUpObservables();

  std::cout<<"Initialized new Sample_Analysis for "<<m_phase<<","<<mode<<std::endl;
  fsel->Output();
}

Sample_Analysis::~Sample_Analysis() 
{
  if (p_analysis) { delete p_analysis; p_analysis = NULL; }
  for (int i=m_obsdata.size()-1;i>=0;i--) {
    if (m_obsdata[i]) { delete m_obsdata[i]; m_obsdata.pop_back(); }
  }
}

void Sample_Analysis::SetOutputPath(const std::string & path)
{
  m_outputpath = path+std::string("/")+m_phase;
}

void Sample_Analysis::DoAnalysis(ATOOLS::Blob_List * const blist, double weight) 
{
  p_analysis->DoAnalysis(blist,weight);
}

void Sample_Analysis::Finish() 
{
  int  mode_dir = 448;
  mkdir(m_outputpath.c_str(),mode_dir); 
  p_analysis->FinishAnalysis(m_prefix+m_outputpath);
}

/*--------------------------------------------------------------------------------------

Read in routines.

---------------------------------------------------------------------------------------*/

void Sample_Analysis::ReadInFinalSelectors(std::ifstream * readin,Final_Selector * fsel)
{  
  std::string buffer, arg;
  unsigned int pos;
  int    kfc;
  double number;
  std::vector<int>    kfcs; 
  std::vector<double> numbers; 
  Flavour flav,flav2;
  Final_Selector_Data fd;

  for (;;) {
    if (readin->eof()) return;
    getline(*readin,buffer);
    buffer += std::string(" ");
    if (buffer[0] != '%' && buffer[0] != '!' && buffer[0] != '#' && buffer.length()>0) {
      if (buffer.find("END_ANALYSIS_PHASE")!=std::string::npos) return;
      if (buffer.find("OBSERVABLES")!=std::string::npos) return;
      if (buffer.find("OUTPUT =")!=std::string::npos) {
	buffer=buffer.substr(buffer.find("OUTPUT =")+8);
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else {
	    pos = buffer.find(std::string(" "));
	    if (pos!=std::string::npos) {
	      m_outputpath = buffer.substr(0,pos);
	    }
	    break;
	  }
	}
      }
      else {
	kfcs.clear();
	numbers.clear();
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else {
	    pos = buffer.find(std::string(" "));
	    if (pos!=std::string::npos) {
	      arg = buffer.substr(0,pos);
	      if (arg.find(".")==std::string::npos) {
		kfc = std::atoi(arg.c_str());
		kfcs.push_back(kfc);
	      }
	      else {
		number = std::atof(arg.c_str());
		numbers.push_back(number);
	      }
	      buffer = buffer.substr(pos);
	    }
	    else break;
	  }
	}
	if (kfcs.size()==1 && numbers.size()>=3) {
	  flav       = Flavour(kf::code(abs(kfcs[0])));
	  if (kfcs[0]<0) flav = flav.Bar();
	  fd.pt_min  =   0.;
	  fd.eta_min = -20.;
	  fd.eta_max = +20.;
	  if (numbers.size()>=0) fd.pt_min  = numbers[0];
	  if (numbers.size()>=1) fd.eta_min = numbers[1];
	  if (numbers.size()>=2) fd.eta_max = numbers[2];
	  if (numbers.size()>=3 && kfcs[0]==93) fd.r_min = numbers[3];
	  if (numbers.size()>=4 && kfcs[0]==93) fd.bf = (bool)numbers[4];
	  fsel->AddSelector(flav,fd);
	}
	else if (kfcs.size()==2 && numbers.size()>=1) {
	  flav       = Flavour(kf::code(abs(kfcs[0])));
	  if (kfcs[0]<0) flav  = flav.Bar();
	  flav2      = Flavour(kf::code(abs(kfcs[1])));
	  if (kfcs[0]<0) flav2 = flav2.Bar();
	  fd.r_min  =   0.;
	  if (numbers.size()>=0) fd.r_min  = numbers[0];
	  fsel->AddSelector(flav,flav2,fd);
	}
	else if (kfcs.size()==3  && numbers.size()==0) {
	  flav       = Flavour(kf::code(abs(kfcs[0])));
	  if (kfcs[0]<0) flav = flav.Bar();
	  fsel->AddSelector(flav,kfcs[1],kfcs[2]);
	}
      }
    }
  }
}

void Sample_Analysis::ReadInObservables(std::ifstream * readin)
{
  std::string buffer, arg;
  unsigned int pos;
  int    kfc;
  double number;
  Observable_Data * obs;

  for (;;) {
    if (readin->eof()) return;
    getline(*readin,buffer);
    buffer += std::string(" ");
    if (buffer[0] != '%' && buffer[0] != '!' && buffer[0] != '#' && buffer.length()>0) {
      if (buffer.find("END_ANALYSIS_PHASE")!=std::string::npos) return;
      if (buffer.find("OUTPUT =")!=std::string::npos) {
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else {
	    pos = buffer.find(std::string(" "));
	    if (pos!=std::string::npos) m_outputpath = buffer.substr(0,pos);
	    else break;
	  }
	}
      }
      else {
	obs = new Observable_Data();
	while(buffer.length()>0) {
	  if (buffer[0]==' ') buffer = buffer.substr(1);
	  else {
	    pos = buffer.find(std::string(" "));
	    if (pos!=std::string::npos) {
	      arg = buffer.substr(0,pos);
	      if (isalpha(arg[0])) {
		if (obs->type==std::string("")) obs->type = arg;
		else obs->keywords.push_back(arg);
	      }
	      else {
		if (arg.find(".")==std::string::npos) {
		  kfc = std::atoi(arg.c_str());
		  obs->ints.push_back(kfc);
		}
		else {
		  number = std::atof(arg.c_str());
		  obs->numbers.push_back(number);
		}
	      }
	      buffer = buffer.substr(pos);
	    }
	    else break;
	  }
	}
	m_obsdata.push_back(obs);
      }
    }
  }
}

void Sample_Analysis::SetUpObservables()
{
  Primitive_Observable_Base * obs;
  Observable_Data * od;
  int               odn,linlog;
  Flavour           flav,flav2;
  std::string type, listname;
  Final_Selector *  fsel;
  bool              found;

  for (size_t i=0;i<m_obsdata.size();++i) {
    msg.Tracking()<<"Try to initialize another observable from read in :"<<std::endl;
    od   = m_obsdata[i];
    od->Output();
    odn  = od->Specify();
    type = od->type;

    switch (odn) {
    case 1: // One-Particle Observables
      if (!(od->ints.size()==2 && od->numbers.size()==2 && od->keywords.size()==1)) {
	msg.Error()<<"Potential Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		   <<"   One particle observable with "
		   <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		   <<". Continue and hope for the best."<<std::endl;
      }
      flav   = Flavour(kf::code(abs(od->ints[0])));
      if (od->ints[0]<0) flav = flav.Bar();
      linlog = 0;
      if (od->keywords[0]==std::string("Log")) linlog = 10;
      
      if (type==std::string("ET")) { 
	obs = new One_Particle_ET(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1]); 
	break; 
      }
      if (type==std::string("PT")) { 
	obs = new One_Particle_PT(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1]); 
	break; 
      }
      if (type==std::string("Eta")) { 
	obs = new One_Particle_Eta(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1]); 
	break; 
      }
      if (type==std::string("E")) { 
	obs = new One_Particle_E(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1]); 
	break; 
      }
      if (type==std::string("EVis")) { 
	obs = new One_Particle_EVis(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1]); 
	break; 
      }
    case 2:
      if (!(od->ints.size()==3 && od->numbers.size()==2 && od->keywords.size()==1)) {
	msg.Error()<<"Potential Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		   <<"   One particle observable with "
		   <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		   <<". Continue and hope for the best."<<std::endl;
      }
      flav   = Flavour(kf::code(abs(od->ints[0])));
      if (od->ints[0]<0) flav = flav.Bar();
      flav2  = Flavour(kf::code(abs(od->ints[1])));
      if (od->ints[1]<0) flav2 = flav2.Bar();
      linlog = 0;
      if (od->keywords[0]==std::string("Log")) linlog = 10;
      
      if (type==std::string("Mass")) { 
	obs = new Two_Particle_Mass(flav,flav2,linlog,od->numbers[0],od->numbers[1],od->ints[2]); 
	break; 
      }
      if (type==std::string("PT2")) { 
	obs = new Two_Particle_PT(flav,flav2,linlog,od->numbers[0],od->numbers[1],od->ints[2]); 
	break; 
      }
      if (type==std::string("Eta2"))  { 
	obs = new Two_Particle_Eta(flav,flav2,linlog,od->numbers[0],od->numbers[1],od->ints[2]); 
	break; 
      }
    case 10:
      if (!(od->ints.size()==5 && od->numbers.size()==2 && od->keywords.size()==1)) {
	msg.Error()<<"Potential Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		   <<"   One particle observable with "
		   <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		   <<". Continue and hope for the best."<<std::endl;
      }
      flav   = Flavour(kf::code(abs(od->ints[0])));
      if (od->ints[0]<0) flav = flav.Bar();
      if (flav==Flavour(kf::jet)) listname = std::string("Analysed");
      linlog = 0;
      if (od->keywords[0]==std::string("Log")) linlog = 10;
      //type,  xmin,xmax,nbins,mode,minn,maxn,listname 
      //linlog,n0,  n1,  i1,   i2,  i3,  i4

      found = false;
      if (m_subsamples.size()>0) {
	for (std::map<Flavour,std::string>::iterator flit=m_subsamples.begin();
	     flit!=m_subsamples.end();flit++) {
	  if (flit->first==flav && flit->second==listname) { 
	    msg.Tracking()<<"List "<<listname<<" already to be projected out."<<std::endl;
	    found = true; break; 
	  }
	}
      }
      if (!found) {
	msg.Tracking()<<"List "<<listname<<" added to be projected out."<<std::endl;
	fsel = new Final_Selector("Analysed",listname,(rpa.gen.Beam1()==Flavour(kf::e)));
	fsel->AddKeepFlavour(flav);
	p_analysis->AddObservable(fsel);
	m_subsamples.insert(std::make_pair(flav,listname));
      }
      if (type==std::string("JetPT"))  { 
	obs = new Jet_PT_Distribution(linlog,od->numbers[0],od->numbers[1],
				      od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname);
	break; 
      }
      if (type==std::string("JetEta"))  { 
	obs = new Jet_Eta_Distribution(linlog,od->numbers[0],od->numbers[1],
				       od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	break; 
      }
      if (type==std::string("JetE"))  { 
	obs = new Jet_E_Distribution(linlog,od->numbers[0],od->numbers[1],
				       od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	break; 
      }
      if (type==std::string("JetDR"))  { 
	obs = new Jet_Differential_Rates(linlog,od->numbers[0],od->numbers[1],
				       od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	break; 
      }
    default:
      msg.Error()<<"Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		 <<"   "<<odn<<"-Particle Observables not yet realized for "<<type<<"."<<std::endl;
      abort();
    }
    p_analysis->AddObservable(obs);
  }
  SetUpSubSamples();
} 

void Sample_Analysis::SetUpSubSamples()
{
  // This is the place to add specific observables ....
} 



