#include "Output_LHEF.H"
#include "CXXFLAGS.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"

#include <iomanip>
#include <stdio.h>
#include <cassert>

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

void Output_LHEF::Header()
{
 
  std::string path(rpa.gen.Variable("SHERPA_DAT_PATH")+"/");
  std::string file(rpa.gen.Variable("RUN_DATA_FILE"));
  
  size_t sep = file.find("|");
  if (sep!=std::string::npos) {
    file.erase(sep);
  }
    
  m_outstream<<"<LesHouchesEvents>"<<std::endl;
  m_outstream<<"<header>"<<std::endl;
  m_outstream<<"<!-- "<<std::endl; 
  m_outstream<<"# created by SHERPA "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION
             <<endl;
  
  Data_Reader dr(" ",";","!","=");
  dr.SetInputPath(path);
  dr.SetInputFile(file);
  
  if (dr.OpenInFile()) {
    m_outstream<<"# Run data extracted from : "<<file<<std::endl; 
    m_outstream<<"--> "<<std::endl; 
    m_outstream<<"<SHRunCard> "<<std::endl; 
    std::vector<std::vector<std::string> > helpsvv;
    if (dr.MatrixFromFile(helpsvv,"")) {
      for (size_t i(0);i<helpsvv.size();++i) {
	for (size_t j(0);j<helpsvv[i].size();++j) {
	  m_outstream<<helpsvv[i][j]<<" ";
	}
	m_outstream<<std::endl;
      }
    }
  }
  m_outstream<<"</SHRunCard> "<<std::endl; 
  m_outstream<<"</header>"<<std::endl;
  
  m_outstream<<"<init>"<<std::endl;
  //run info to be dumped here
  Flavour Beam1 = rpa.gen.Beam1();
  Flavour Beam2 = rpa.gen.Beam2();
  int IDBMUP1 = Beam1.HepEvt();
  int IDBMUP2 = Beam2.HepEvt();
  double EBMUP1 = rpa.gen.PBeam(0)[0];
  double EBMUP2 = rpa.gen.PBeam(1)[0];
  
  int IDWTUP(ToType<int>(rpa.gen.Variable("EVENT_GENERATION_MODE"))==0?1:3);
  int NPRUP = 1;
  int PDFGUP1 = 0;
  int PDFGUP2 = 0;
  int PDFSUP1 =10041;
  int PDFSUP2 =10041;

  m_outstream<<std::setprecision(10);
  m_outstream<<std::setw(6)<<IDBMUP1<<" "
	     <<std::setw(6)<<IDBMUP2<<" "
	     <<std::setw(11)<<EBMUP1<<" "
	     <<std::setw(11)<<EBMUP2<<" "
	     <<std::setw(3)<<PDFGUP1<<" "
	     <<std::setw(3)<<PDFGUP2<<" "
	     <<std::setw(6)<<PDFSUP1<<" "
	     <<std::setw(6)<<PDFSUP2<<" "
	     <<std::setw(4)<<IDWTUP<<" "
	     <<std::setw(4)<<NPRUP<<std::endl;

  //process information
  m_outstream<<std::setw(18)<<m_xs<<" "
	     <<std::setw(18)<<m_xserr<<" "
	     <<std::setw(18)<<m_max<<" "
	     <<std::setw(4)<<NPRUP<<std::endl; 
  m_outstream<<"</init>"<<std::endl;
}

void Output_LHEF::Output(Blob_List* blobs, const double weight) 
{
  m_outstream<<"<event>"<<std::endl;
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit){
    if ((*blit)->Type()==ATOOLS::btp::Signal_Process) {
      //LHE event information
      int NUP = (*blit)->NInP()+(*blit)->NOutP();
      int IDPRUP = 1;
      double XWGTUP = weight;
      double SCALUP = sqrt((*(*blit))["Factorisation_Scale"]->Get<double>());
      double AQEDUP = -1.;
      double AQCDUP = -1.;
      m_outstream<<std::setprecision(10);
      m_outstream<<std::setiosflags(std::ios::scientific);
      m_outstream<<std::setw(4)<<NUP<<" "
		 <<std::setw(4)<<IDPRUP<<" "
		 <<std::setw(18)<<XWGTUP<<" "
		 <<std::setw(18)<<SCALUP<<" "
		 <<std::setw(18)<<AQEDUP<<" "
		 <<std::setw(18)<<AQCDUP<<std::endl;
      for (int i=0;i<(*blit)->NInP();i++)
	m_outstream<<std::setw(8)<<(*blit)->InParticle(i)->Flav().HepEvt()<<" -1  0  0 "
		   <<std::setw(4)<<(*blit)->InParticle(i)->GetFlow(1)<<" "
		   <<std::setw(4)<<(*blit)->InParticle(i)->GetFlow(2)<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->Momentum()[1]<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->Momentum()[2]<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->Momentum()[3]<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->Momentum()[0]<<" "
		   <<std::setw(18)<<(*blit)->InParticle(i)->FinalMass()<<" "
		   <<" 0  9"<<std::endl;
      for (int i=0;i<(*blit)->NOutP();i++)
	m_outstream<<std::setw(8)<<(*blit)->OutParticle(i)->Flav().HepEvt()<<"  1  1  2 "
		   <<std::setw(4)<<(*blit)->OutParticle(i)->GetFlow(1)<<" "
		   <<std::setw(4)<<(*blit)->OutParticle(i)->GetFlow(2)<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->Momentum()[1]<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->Momentum()[2]<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->Momentum()[3]<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->Momentum()[0]<<" "
		   <<std::setw(18)<<(*blit)->OutParticle(i)->FinalMass()<<" "
		   <<" 0  9"<<std::endl;
      m_outstream<<std::resetiosflags(std::ios::fixed);
    }
  }
  m_outstream<<"</event>"<<std::endl;
}

void Output_LHEF::Footer(std::string number)
{
  string newfile=m_basename+"."+number+m_ext;
  string footer = std::string("</LesHouchesEvents>");
  m_outstream<<footer<<std::endl;
}

void Output_LHEF::SetXS(const double& xs, const double& xserr) {
  //std::cout<<" set xs to "<<xs<<"  +/- "<<xserr<<std::endl; 
  m_xs = xs;
  m_xserr = xserr;
  m_max =1.;
}

