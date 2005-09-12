#include "ISR_Statistics.H"

using namespace ANALYSIS;

DECLARE_GETTER(ISR_Statistics_Getter,"ISRStats",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base * 
ISR_Statistics_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new ISR_Statistics(listname);
}

void ISR_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"
#include "ISR_Info.H"
#include "MyStrStream.H"
#include "Shell_Tools.H"

#include <fstream>

using namespace ATOOLS;

ISR_Statistics::ISR_Statistics(const std::string & listname, int type)
{
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-200.,400.,200));
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-10000.,1000.,200));
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-10.,10.,200));
  m_observables.push_back(new ATOOLS::Histogram(0,0.,20000.,200));
  m_observables.push_back(new ATOOLS::Histogram(10,0.01,100.,200));
  m_observables.push_back(new ATOOLS::Histogram(0,0.,1.,200));
  m_observables.push_back(new ATOOLS::Histogram(10,0.01,100.,200));
  m_name  = "ISR_Statistics_";
  m_type  = type;
  m_listname    = listname;
  m_splitt_flag = true;
}

void ISR_Statistics::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
  Particle_List * pl_cut = p_ana->GetParticleList(m_listname);
  if (pl_cut && pl_cut->size()==0) {
    for (size_t i=0;i<m_observables.size();++i) m_observables[i]->Insert(1.,0.,ncount);
    return;
  }
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type()==btp::Signal_Process) {
      Blob_Data_Base * info=(*(*bit))["ISR_Info_cms"];
      if (info) {
	std::vector<double> isrinfo=info->Get<std::vector<double> >();
	m_observables[0]->Insert(isrinfo[PDF::iic::E_1],weight,ncount);
	m_observables[2]->Insert(isrinfo[PDF::iic::t_1],weight,ncount);
	m_observables[4]->Insert(isrinfo[PDF::iic::Em_1],weight,ncount);
	m_observables[6]->Insert(isrinfo[PDF::iic::mu_1],weight,ncount);
	m_observables[7]->Insert(isrinfo[PDF::iic::mu_1]/
				 isrinfo[PDF::iic::mu_2],weight,ncount);
      }
      info=(*(*bit))["ISR_Info_lab"];
      if (info) {
	std::vector<double> isrinfo=info->Get<std::vector<double> >();
	m_observables[1]->Insert(isrinfo[PDF::iic::E_1],weight,ncount);
	m_observables[3]->Insert(isrinfo[PDF::iic::t_1],weight,ncount);
	m_observables[5]->Insert(isrinfo[PDF::iic::Em_1],weight,ncount);
	m_observables[8]->Insert(isrinfo[PDF::iic::z_1],weight,ncount);
	m_observables[9]->Insert(isrinfo[PDF::iic::z_1]/
				 isrinfo[PDF::iic::z_2],weight,ncount);
      }
    }
  }
}

Primitive_Observable_Base * ISR_Statistics::Copy() const 
{
  return new ISR_Statistics(m_listname,m_type);
}

void ISR_Statistics::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_observables.size();++i) {
    m_observables[i]->Finalize();
    if (scale!=1.) m_observables[i]->Scale(scale);
    m_observables[i]->Output();
  }
}

void ISR_Statistics::Output(const std::string & pname) {
  int  mode_dir = 448;
  ATOOLS::MakeDir((pname).c_str(),mode_dir); 
  for (size_t i=0; i<m_observables.size();++i) {
    std::string fname;
    MyStrStream s;
    s<<i;
    s<<".dat"; 
    s>>fname;
    m_observables[i]->Output((pname+std::string("/")+m_name+fname).c_str());
  }
}

Primitive_Observable_Base & ISR_Statistics::operator+=(const Primitive_Observable_Base & ob)
{
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    msg.Error()<<"ERROR: in ISR_Statistics::operator+=  in"<<m_name<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    return *this;
  }

  ISR_Statistics * job = ((ISR_Statistics*)(&ob));

  if (m_observables.size()==job->m_observables.size()) {
    for (size_t i=0; i<m_observables.size();++i) {
      (*m_observables[i])+=(*job->m_observables[i]);
    }
  }
  return *this;
}

#ifdef USING__ROOT
#include "My_Root.H"
#include "TH2D.h"

DECLARE_GETTER(Sprime_Y_Distribution_Getter,"SprimeY",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *
Sprime_Y_Distribution_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new Sprime_Y_Distribution(-10.0,0.0,200,-10.0,10.0,200);
}

void Sprime_Y_Distribution_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

Sprime_Y_Distribution::
Sprime_Y_Distribution(const double spmin,const double spmax,const size_t spbins,
		      const double ymin,const double ymax,const size_t ybins):
  m_ymin(ymin),
  m_ymax(ymax),
  m_ybins(ybins)
{ 
  m_xmin=spmin;
  m_xmax=spmax;
  m_nbins=spbins;
  (*MYROOT::myroot)(new TH2D((ATOOLS::ToString((long int)this)+"ME").c_str(),
			     "Sprime_Y_ME",m_nbins,m_xmin,m_xmax,m_ybins,m_ymin,m_ymax),
		    "Sprime_Y_ME");
  (*MYROOT::myroot)(new TH2D((ATOOLS::ToString((long int)this)+"PS").c_str(),
			     "Sprime_Y_PS",m_nbins,m_xmin,m_xmax,m_ybins,m_ymin,m_ymax),
		    "Sprime_Y_PS");
}

#define ANALYSE__Phase_Space_Handler

void Sprime_Y_Distribution::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
#ifndef ANALYSE__Phase_Space_Handler
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type()==btp::Signal_Process) {
      Blob_Data_Base *info=(*(*bit))["ISR_Info_lab"];
      if (info) {
	std::vector<double> isrinfo=info->Get<std::vector<double> >();
	double tau=isrinfo[PDF::iic::E_1]*isrinfo[PDF::iic::E_2];
	double y=0.5*log(isrinfo[PDF::iic::E_1]/isrinfo[PDF::iic::E_2]);
	((TH2D*)(*MYROOT::myroot)["Sprime_Y_ME"])->Fill(tau,y,weight);
	((TH2D*)(*MYROOT::myroot)["Sprime_Y_PS"])->Fill(tau,y,1.0);
      }
    }
  }
#endif
} 

Primitive_Observable_Base * Sprime_Y_Distribution::Copy() const 
{
  return new Sprime_Y_Distribution(m_xmin,m_xmax,m_nbins,m_ymin,m_ymax,m_ybins);
}
#endif

