#include "Analysis_Handler.H"
#include "MyStrStream.H"
#include "Run_Parameter.H"
#include "Message.H"

#define ANALYSE__ISR_Statistics
#define ANALYSE__Shower_Statistics

#include "Final_Selector.H"
#include "Primitive_Calorimeter.H"
#include "Primitive_Observable.H"
#include "Jet_Cone_Distribution.H"
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
  if (Specify()==20) {
    flav1 = Flavour(kf::code(abs(ints[0])));
    if (ints[0]<0) flav1=flav1.Bar();
    msg.Out()<<"Obs : "<<type<<" for "<<flav1<<" : "
	     <<"Range : "<<numbers[0]<<" ... "<<numbers[1]<<" in "<<ints[1]<<" bins,"
	     <<" extra : "<<keywords[0]<<std::endl;
  }
}

int Observable_Data::Specify() {
  if (type==std::string("ET") || type==std::string("PT") ||
      type==std::string("Eta") || type==std::string("E") ||
      type==std::string("EVis"))                                 return 1;
  if (type==std::string("Mass") || type==std::string("PT2") ||
      type==std::string("Eta2") || type==std::string("SPT2") ||
      type==std::string("Angles"))                               return 2;
  if (type==std::string("JetPT") || type==std::string("JetEta") || 
      type==std::string("JetE") || type==std::string("DiffJet") || 
      type==std::string("JetDR") || type==std::string("JetDEta") || 
      type==std::string("JetDPhi"))                              return 10;
  if (type==std::string("JetCone") || type==std::string("JetConeShape") ||
      type==std::string("JetConeDep"))                           return 20; 
  return -1;
}

Analysis_Handler::Analysis_Handler(std::ifstream * readin, std::string _phase,
				   const std::string & prefix) :
  m_phase(_phase), m_outputpath(std::string("./")+_phase), m_prefix(prefix), 
  p_detector(NULL), p_analysis(NULL)
{  
  msg.Info()<<"Initialize new Analysis_Handler for "<<_phase<<std::endl;
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
      if (phasemode==std::string("MI")) {
	mode  = mode|ANALYSIS::do_mi;
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
      }
    }
  }
  p_analysis            = new Primitive_Analysis(m_phase,mode);
  int jet_mode          = 0;
  if (rpa.gen.Beam1()==Flavour(kf::e)) jet_mode = 1;
  ReadInDetector(readin);
  if (p_detector) p_detector->Print();

  Final_Selector * fsel; 
  if (p_detector!=NULL) fsel = new Final_Selector(std::string("FinalState"),
						  std::string("Analysed"));
                   else fsel = new Final_Selector(std::string("FinalState"),
						  std::string("Analysed"),jet_mode);
  ReadInFinalSelectors(readin,fsel);
  p_analysis->AddObservable(fsel);
  ReadInObservables(readin);
  SetUpObservables();
  //  SetUpSubSamples();
  
  ATOOLS::Exception_Handler::AddTerminatorObject(this);

  if (msg.LevelIsInfo()) {
    msg.Out()<<"Initialized new Analysis_Handler for "<<m_phase<<","<<mode<<std::endl;
    fsel->Output();
  }
}

Analysis_Handler::~Analysis_Handler() 
{
  if (p_detector) { delete p_detector; p_detector = NULL; }
  if (p_analysis) { delete p_analysis; p_analysis = NULL; }
  for (int i=m_obsdata.size()-1;i>=0;i--) {
    if (m_obsdata[i]) { delete m_obsdata[i]; m_obsdata.pop_back(); }
  }
}

void Analysis_Handler::SetOutputPath(const std::string & path)
{
  m_outputpath = path+std::string("/")+m_phase;
}

void Analysis_Handler::DoAnalysis(ATOOLS::Blob_List * const blist, double weight) 
{
  if (p_detector) p_detector->Fill(blist);
  p_analysis->DoAnalysis(blist,weight);
  if (p_detector) p_detector->Reset();
}

void Analysis_Handler::Finish() 
{
  int  mode_dir = 448;
  mkdir(m_outputpath.c_str(),mode_dir); 
  p_analysis->FinishAnalysis(m_prefix+m_outputpath);
}

/*--------------------------------------------------------------------------------------

Read in routines.

---------------------------------------------------------------------------------------*/

void Analysis_Handler::ReadInDetector(std::ifstream * readin)
{
  std::string         buffer, arg;
  size_t              pos;
  bool                detector_on=false;
  int                 pixel;
  double              spec;
  std::vector<int>    pixels; 
  std::vector<double> specs; 
  for (;;) {
    if (readin->eof()) return;
    getline(*readin,buffer);
    buffer += std::string(" ");
    if (buffer[0] != '%' && buffer[0] != '!' && buffer[0] != '#' && buffer.length()>0) {
      if (buffer.find("END_ANALYSIS_PHASE")!=std::string::npos) return;
      if (buffer.find("OBSERVABLES")!=std::string::npos)        return;
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
	return;
      }
      if (buffer.find("DETECTOR")!=std::string::npos) {
	detector_on = true;
	continue;
      }
      if (!detector_on) return;
      pixels.clear();
      specs.clear();
      while(buffer.length()>0) {
	if (buffer[0]==' ') buffer = buffer.substr(1);
	else {
	  pos = buffer.find(std::string(" "));
	  if (pos!=std::string::npos) {
	    arg = buffer.substr(0,pos);
	    if (arg.find(".")==std::string::npos) {
	      pixel = std::atoi(arg.c_str());
	      pixels.push_back(pixel);
	    }
	    else {
	      spec = std::atof(arg.c_str());
	      specs.push_back(spec);
	    }
	    buffer = buffer.substr(pos);
	  }
	  else break;
	}
      }
      if (pixels.size()>0 && specs.size()>0) {
	if (p_detector==NULL) p_detector = new Primitive_Detector();
	switch(pixels[0]) {
	case 1: 
	  // Primitive_Calorimeter
	  if (specs.size()<2||pixels.size()<3) {
	    msg.Error()<<"Error in Analysis_Handler::ReadInDetector."<<std::endl
		       <<"   Cannot specify calorimenter. Too little input :"
		       <<specs.size()<<"/"<<pixels.size()<<"."<<std::endl
		       <<"   Abort the run."<<std::endl;
	    abort();
	  }
	  p_detector->Add(new Primitive_Calorimeter(specs[0],specs[1],pixels[1],pixels[2]));
	  break;
	default:continue;
	}
      }
    }
  }
}

void Analysis_Handler::ReadInFinalSelectors(std::ifstream * readin,
					    Final_Selector *& fsel)
{  
  std::string buffer, arg;
  size_t              pos;
  int                 kfc;
  double              number;
  std::vector<int>    kfcs; 
  std::vector<double> numbers; 
  Flavour             flav,flav2;
  Final_Selector_Data fd;

  bool mode = (p_detector==NULL);
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
	  if (!mode && (flav==Flavour(kf::jet) || flav==Flavour(kf::bjet))) {
	    Primitive_Calorimeter * Hcal = 
	      dynamic_cast<Primitive_Calorimeter * >
	      (p_detector->GetElement(std::string("Hadronic Calorimeter")));
	    Calorimeter_Cone * Hcal_jetfinder = 
	      new Calorimeter_Cone(numbers[0],numbers[3],Hcal);
	    Hcal_jetfinder->SetEtaRangeForJets(numbers[1],numbers[2],int(numbers[4]>0.));
	    fd.eta_min = numbers[1];
	    fd.eta_max = numbers[2];
	    fd.pt_min  = numbers[0];
	    fd.r_min   = numbers[3];
	    fsel->AddSelector(flav,fd,Hcal_jetfinder);
	  }
	  else {
	    fd.pt_min  =   0.;
	    fd.eta_min = -20.;
	    fd.eta_max = +20.;
	    if (numbers.size()>0) fd.pt_min  = numbers[0];
	    if (numbers.size()>1) fd.eta_min = numbers[1];
	    if (numbers.size()>2) fd.eta_max = numbers[2];
	    if (numbers.size()>3 && kfcs[0]==93) fd.r_min = numbers[3];
	    if (numbers.size()>4 && kfcs[0]==93) fd.bf = (bool)numbers[4];
	    fsel->AddSelector(flav,fd);
	  }
	}
	else if (kfcs.size()==2 && numbers.size()>=1) {
	  flav       = Flavour(kf::code(abs(kfcs[0])));
	  if (kfcs[0]<0) flav  = flav.Bar();
	  flav2      = Flavour(kf::code(abs(kfcs[1])));
	  if (kfcs[0]<0) flav2 = flav2.Bar();
	  fd.r_min  =   0.;
	  if (numbers.size()>0) fd.r_min  = numbers[0];
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

void Analysis_Handler::ReadInObservables(std::ifstream * readin)
{
  std::string buffer, arg;
  size_t pos;
  int    kfc;
  double number;
  Observable_Data * obs=NULL;

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

void Analysis_Handler::SetUpObservables()
{
  Primitive_Observable_Base * obs = NULL;
  Observable_Data           * od  = NULL;
  int                         odn,linlog;
  Flavour                     flav,flav2;
  std::string                 type;
  
  std::string                 listname = std::string("Analysed");

  for (size_t i=0;i<m_obsdata.size();++i) {
    msg.Tracking()<<"Try to initialize another observable from read in :"<<std::endl;
    od   = m_obsdata[i];
    od->Output();
    odn  = od->Specify();
    type = od->type;

    switch (odn) {
    case 1: // One-Particle Observables
      if (!(od->ints.size()==2 && od->numbers.size()==2 && od->keywords.size()==1)) {
	msg.Error()<<"Potential Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
		   <<"   One particle observable with "
		   <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		   <<". Continue and hope for the best."<<std::endl;
      }
      flav   = Flavour(kf::code(abs(od->ints[0])));
      if (od->ints[0]<0) flav = flav.Bar();
      linlog = 0;
      if (od->keywords[0]==std::string("Log")) linlog = 10;
      
      if (type==std::string("ET")) { 
	obs = new One_Particle_ET(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1],listname); 
	break; 
      }
      if (type==std::string("PT")) { 
	obs = new One_Particle_PT(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1],listname); 
	break; 
      }
      if (type==std::string("Eta")) { 
	obs = new One_Particle_Eta(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1],listname); 
	break; 
      }
      if (type==std::string("E")) { 
	obs = new One_Particle_E(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1],listname); 
	break; 
      }
      if (type==std::string("EVis")) { 
	obs = new One_Particle_EVis(flav,linlog,od->numbers[0],od->numbers[1],od->ints[1],listname); 
	break; 
      }
    case 2:
      if (!(od->ints.size()==3 && od->numbers.size()==2 && od->keywords.size()==1)) {
	msg.Error()<<"Potential Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
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
	obs = new Two_Particle_Mass(flav,flav2,linlog,od->numbers[0],od->numbers[1],
				    od->ints[2],listname); 
	break; 
      }
      if (type==std::string("PT2")) { 
	obs = new Two_Particle_PT(flav,flav2,linlog,od->numbers[0],od->numbers[1],
				  od->ints[2],listname); 
	break; 
      }
      if (type==std::string("SPT2")) { 
	obs = new Two_Particle_Scalar_PT(flav,flav2,linlog,od->numbers[0],od->numbers[1],
					 od->ints[2],listname); 
	break; 
      }
      if (type==std::string("Angles")) { 
#ifdef ROOT_SUPPORT
	obs = new Two_Particle_Angles(flav,flav2,linlog,od->numbers[0],od->numbers[1],od->ints[2]); 
#endif
	break; 
      }
      if (type==std::string("Eta2"))  { 
	obs = new Two_Particle_Eta(flav,flav2,linlog,od->numbers[0],od->numbers[1],
				   od->ints[2],listname); 
	break; 
      }
    case 10:
      {
	if (!(od->ints.size()==5 && od->numbers.size()==2 && od->keywords.size()==1)) {
	  msg.Error()<<"Potential Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		     <<"   One particle observable with "
		     <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		     <<". Continue and hope for the best."<<std::endl;
	}
	
	flav   = Flavour(kf::code(abs(od->ints[0])));
	if (od->ints[0]<0) flav = flav.Bar();
	
	linlog = 0;
	if (od->keywords[0]==std::string("Log")) linlog = 10;

	if (type!="DiffJet") {
	  // jet observables have a "particle list" containing only "jets"
	  if (flav==Flavour(kf::jet)) listname = std::string("AnalysedJets");
	  else listname = std::string("Analysed")+flav.Name();
	  //type,  xmin,xmax,nbins,mode,minn,maxn,listname 
	  //linlog,n0,  n1,  i1,   i2,  i3,  i4

	  bool	found = false;
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
	    Final_Selector * fsel = 
	      new Final_Selector("Analysed",listname,(rpa.gen.Beam1()==Flavour(kf::e)));
	    fsel->AddKeepFlavour(flav);
	    p_analysis->AddObservable(fsel);
	    m_subsamples.insert(std::make_pair(flav,listname));
	  }
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
	if (type==std::string("DiffJet"))  { 
	  obs = new Jet_Differential_Rates(linlog,od->numbers[0],od->numbers[1],
					   od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	  break; 
	}
	if (type==std::string("JetDR"))  {
	  obs = new Jet_DeltaR_Distribution(linlog,od->numbers[0],od->numbers[1],
					   od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	  break; 
	}
	if (type==std::string("JetDEta"))  {
	  obs = new Jet_DeltaEta_Distribution(linlog,od->numbers[0],od->numbers[1],
					   od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	  break; 
	}
	if (type==std::string("JetDPhi"))  {
	  obs = new Jet_DeltaPhi_Distribution(linlog,od->numbers[0],od->numbers[1],
					   od->ints[1],od->ints[2],od->ints[3],od->ints[4],listname); 
	  break; 
	}
      }
    case 20:
      {
 	if (!(od->ints.size()==2 && od->numbers.size()==3 && od->keywords.size()==1) ||
	    !(od->ints.size()==2 && od->numbers.size()==5 && od->keywords.size()==1) ||
	    !(od->ints.size()==4 && od->numbers.size()==6 && od->keywords.size()==1) ||
	    !(od->ints.size()==4 && od->numbers.size()==3 && od->keywords.size()==1) ||
	    !(od->ints.size()==4 && od->numbers.size()==5 && od->keywords.size()==1) ) {
	  msg.Error()<<"Potential Error in Sample_Analysis::SetUpSubObservables()"<<std::endl
		     <<"   Detector observable with "
		     <<od->ints.size()<<" "<<od->numbers.size()<<" "<<od->keywords.size()
		     <<". Continue and hope for the best."<<std::endl;
	}
	
	flav   = Flavour(kf::code(abs(od->ints[0])));
	if (od->ints[0]<0) flav = flav.Bar();
	
	linlog = 0;
	if (od->keywords[0]==std::string("Log")) linlog = 10;
	if (type==std::string("JetConeShape"))  {
	  if (flav!=Flavour(kf::jet)) {
	    msg.Error()<<"Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
		       <<"   Jet cone number not defined for non-jets."<<std::endl
		       <<"   Don't initialize an observable."<<std::endl;
	    continue;
	  }
	  Primitive_Calorimeter * Hcal = 
	    dynamic_cast<Primitive_Calorimeter * >
	    (p_detector->GetElement(std::string("Hadronic Calorimeter")));
	  Calorimeter_Cone * Hcal_jetfinder = 
	    new Calorimeter_Cone(od->numbers[0],od->numbers[3],Hcal);
	  Hcal_jetfinder->SetEtaRangeForJets(od->numbers[1],od->numbers[2],1);
	  obs = new Jet_Cone_Shape(linlog,od->numbers[4],od->numbers[5],
				   od->ints[1],od->ints[2],od->ints[3],Hcal_jetfinder);
	  break; 
	}
	if (type==std::string("JetConeDep"))  {
	  if (flav!=Flavour(kf::jet)) {
	    msg.Error()<<"Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
		       <<"   Jet cone dependence not defined for non-jets."<<std::endl
		       <<"   Don't initialize an observable."<<std::endl;
	    continue;
	  }
	  Primitive_Calorimeter * Hcal = 
	    dynamic_cast<Primitive_Calorimeter * >
	    (p_detector->GetElement(std::string("Hadronic Calorimeter")));
	  if (od->numbers.size()==3) {
	    obs = new Jet_Cone_Dependence(linlog,od->numbers[0],
					  od->numbers[1],od->numbers[2],
					  od->ints[1],od->ints[2],od->ints[3],Hcal);
	  }
	  else if (od->numbers.size()==5) {
	    obs = new Jet_Cone_Dependence(linlog,od->numbers[0],
					  od->numbers[1],od->numbers[2],
					  od->numbers[3],od->numbers[4],
					  od->ints[1],od->ints[2],od->ints[3],Hcal);
	  }
	  break;
	}
	if (type==std::string("JetCone"))  {
	  if (flav!=Flavour(kf::jet)) {
	    msg.Error()<<"Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
		       <<"   Jet cone number not defined for non-jets."<<std::endl
		       <<"   Don't initialize an observable."<<std::endl;
	    continue;
	  }
	  Primitive_Calorimeter * Hcal = 
	    dynamic_cast<Primitive_Calorimeter * >
	    (p_detector->GetElement(std::string("Hadronic Calorimeter")));
	  if (od->numbers.size()==3) {
	    obs = new Jet_Cone_Distribution(linlog,od->numbers[0],
					    od->numbers[1],od->numbers[2],
					    od->ints[1],Hcal);
	  }
	  else if (od->numbers.size()==5) {
	    obs = new Jet_Cone_Distribution(linlog,od->numbers[0],
					    od->numbers[1],od->numbers[2],
					    od->numbers[3],od->numbers[4],
					    od->ints[1],Hcal);
	  }
	  break;
	}
      }
    default:
      msg.Error()<<"Error in Analysis_Handler::SetUpSubObservables()"<<std::endl
		 <<"   "<<odn<<"-Particle Observables not yet realized for "<<type<<"."<<std::endl;
      abort();
    }
    p_analysis->AddObservable(obs);
  }
  //SetUpSubSamples();
} 

void Analysis_Handler::SetUpSubSamples()
{
  // This is the place to add specific observables ....

//   Final_Selector * fsel =  new Final_Selector("FinalState","Analysed");
//   // pure teilchen eigenschaften
//   Final_Selector_Data fd;
//   fd.rmin   =1.;
//   fd.pt_min =15.;
//   fd.eta_min=-2.;
//   fd.eta_max=+2.;
//   fsel->AddSelector(Flavour(kf::jet),fd);
//   Final_Selector_Data fcorr;
//   fcorr.rmin = 0.7;
//   fsel->AddSelector(Flavour(kf::jet),Flavour(kf::lepton),fcorr);
//   p_analysis->AddObservable(fsel);


  // --- Andreas Analysis ---
  //
  // special analysis   W + jets
  //


  // all jets
  Final_Selector * fsel = new Final_Selector("Analysed","Jets",-1);
  fsel->AddKeepFlavour(Flavour(kf::jet));
  p_analysis->AddObservable(fsel);

  // all lepton
  fsel = new Final_Selector("Analysed","Leptons",-1);
  fsel->AddKeepFlavour(Flavour(kf::lepton));
  p_analysis->AddObservable(fsel);


  // leptons in excl. 1 Jet events
  fsel = new Final_Selector("Analysed","Leptons(1j)",-1);
  fsel->AddKeepFlavour(Flavour(kf::lepton));
  fsel->AddSelector(Flavour(kf::jet),1,1);
  p_analysis->AddObservable(fsel);

  // leptons in excl. 2 Jet events
  fsel = new Final_Selector("Analysed","Leptons(2j)",-1);
  fsel->AddKeepFlavour(Flavour(kf::lepton));
  fsel->AddSelector(Flavour(kf::jet),2,2);
  p_analysis->AddObservable(fsel);

  // leptons in incl. 3 Jet events
  fsel = new Final_Selector("Analysed","Leptons(3j+)",-1);
  fsel->AddKeepFlavour(Flavour(kf::lepton));
  fsel->AddSelector(Flavour(kf::jet),3,-1);
  p_analysis->AddObservable(fsel);


  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"FinalState"));
  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"Jets"));
  p_analysis->AddObservable(new One_Particle_PT(Flavour(kf::e),
						00,0.,200.,100,"FinalState"));
  p_analysis->AddObservable(new One_Particle_Eta(Flavour(kf::e),
						 00,-5.,5.,50,"FinalState"));
  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"Analysed"));
  p_analysis->AddObservable(new One_Particle_PT(Flavour(kf::e),
						00,0.,200.,100,"Analysed"));
  p_analysis->AddObservable(new One_Particle_Eta(Flavour(kf::e),
						 00,-5.,5.,50,"Analysed"));
  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"Leptons(1j)"));
  p_analysis->AddObservable(new One_Particle_PT(Flavour(kf::e),
						00,0.,200.,100,"Leptons(1j)"));
  p_analysis->AddObservable(new One_Particle_Eta(Flavour(kf::e),
						 00,-5.,5.,50,"Leptons(1j)"));
  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"Leptons(2j)"));
  p_analysis->AddObservable(new One_Particle_PT(Flavour(kf::e),
						00,0.,200.,100,"Leptons(2j)"));
  p_analysis->AddObservable(new One_Particle_Eta(Flavour(kf::e),
						 00,-5.,5.,50,"Leptons(2j)"));
  p_analysis->AddObservable(new Multiplicity(00,-0.5,50.5,51,"Leptons(3j+)"));
  p_analysis->AddObservable(new One_Particle_PT(Flavour(kf::e),
						00,0.,200.,100,"Leptons(3j+)"));
  p_analysis->AddObservable(new One_Particle_Eta(Flavour(kf::e),
						 00,-5.,5.,50,"Leptons(3j+)"));

  p_analysis->AddObservable(new Two_Particle_PT(Flavour(kf::e),Flavour(kf::nue).Bar(),
						00,0.,200.,100,"FinalState"));
//   p_analysis->AddObservable(new Two_Particle_PT(Flavour(kf::e),Flavour(kf::nue).Bar(),
// 						00,0.,100.,200,"FinalState"));
  p_analysis->AddObservable(new Two_Particle_Eta(Flavour(kf::e),Flavour(kf::nue).Bar(),
						00,-5.,5.,50,"FinalState"));

  /*
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,84,1,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,84,1,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,84,1,3,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,200.,84,2,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,200.,84,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,200.,84,2,3,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.,2.,10,2,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.,2.,10,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.,2.,10,2,3,5,"Jets"));
  */
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,1,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,1,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,1,3,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,2,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_PT_Distribution(00,0.,210.,42,2,3,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,1,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,1,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,1,3,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,2,1,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_Eta_Distribution(00,-2.5,2.5,20,2,3,5,"Jets"));


  p_analysis->AddObservable(new Jet_Differential_Rates(10,2.e-2,200.,100,1,1,5,"Analysed")); 

  p_analysis->AddObservable(new Jet_DeltaR_Distribution(00,0.,3.7,37,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_DeltaPhi_Distribution(00,0.,3.7,37,2,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_DeltaR_Distribution(00,0.,3.7,37,1,2,5,"Jets"));
  p_analysis->AddObservable(new Jet_DeltaPhi_Distribution(00,0.,3.7,37,1,2,5,"Jets"));
  

} 

void Analysis_Handler::PrepareTerminate()
{
  Finish();
}


