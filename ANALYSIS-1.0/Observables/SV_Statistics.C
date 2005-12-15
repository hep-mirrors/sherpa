#include "Primitive_Observable_Base.H"
#include "Primitive_Analysis.H"
#include "Veto_Info.H"
#include "MyStrStream.H"
#include "Shell_Tools.H"
#include "Exception.H"
#include "Veto_Info.H"
#include <fstream>
#include <iostream>

#ifdef USING__Veto_Info

namespace ANALYSIS {

  class SV_Statistics: public Primitive_Observable_Base {  
  protected:

    std::vector<ATOOLS::Histogram*> m_observables;

    std::string m_stype;
    int         m_beam;

  public:

    SV_Statistics(const std::string &type,const std::string &inlistname);

    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base * Copy() const;

    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);

  };

}

using namespace ANALYSIS;
using namespace APACIC;

DECLARE_GETTER(SV_Statistics_Getter,"SVStats",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base * 
SV_Statistics_Getter::operator()(const String_Matrix &parameters) const
{
  if (parameters.size()<1 || parameters[0].size()<1) 
    THROW(fatal_error,"Statistics observable needs arguments.");
  std::string listname="FinalState", type=parameters[0][0];
  if (parameters.size()>1 && parameters[1].size()>0) listname=parameters[0][1];
  return new SV_Statistics(type,listname);
}

void SV_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"IS0|IS1|FS|IFS [list]"; 
}

using namespace ATOOLS;

SV_Statistics::SV_Statistics(const std::string &type,const std::string & listname):
  m_stype(type)
{
  if (m_stype=="FS") m_beam=-1;
  else if (m_stype=="IFS") m_beam=-2;
  else if (m_stype=="IS0") m_beam=0;
  else if (m_stype=="IS1") m_beam=1;
  else THROW(fatal_error,"No shower type specified.");
  m_name  = m_stype+"_VS_";
  m_type  = 0;
  m_listname    = listname;
  m_splitt_flag = true;
}

void SV_Statistics::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
  Particle_List * pl_cut = p_ana->GetParticleList(m_listname);
  if (pl_cut && pl_cut->size()==0) {
    for (size_t i=0;i<m_observables.size();++i) m_observables[i]->Insert(0.0,0.0,ncount);
    return;
  }
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    Blob_Data_Base *info(NULL);
    if (m_beam==-2) {
      if ((*bit)->Type()==btp::FS_Shower) info=(*(*bit))["IFS_VS"];
    }
    if (m_beam==-1) {
      if ((*bit)->Type()==btp::FS_Shower) info=(*(*bit))["FS_VS"];
    }
    else {
      if ((*bit)->Type()==btp::IS_Shower && 
	  (*bit)->Beam()==m_beam) info=(*(*bit))["IS_VS"];
    }
    if (info) {
      std::vector<int> sinfo=info->Get<std::vector<int> >();
      if (sinfo.size()>m_observables.size()) {
	for (size_t i(m_observables.size());i<sinfo.size();++i) {
	  m_observables.push_back(new ATOOLS::Histogram(0,0.0,10.0,10));
	  m_observables[i]->SetFills((long int)m_observables[0]->Fills());
	}
      }
      for (size_t i(0);i<sinfo.size();++i) {
	double ex(0.0);
	for (size_t j(svc::no_veto);j<svc::size;j=j<<1) {
	  m_observables[i]->Insert(ex++,sinfo[i]&j?weight:0.0,0);
	}
	m_observables[i]->Insert(0.0,0.0,ncount);
      }
      for (size_t i(sinfo.size());i<m_observables.size();++i)
	m_observables[i]->Insert(0.0,0.0,ncount);
      break;
    }
  }
}

Primitive_Observable_Base * SV_Statistics::Copy() const 
{
  return new SV_Statistics(m_stype,m_listname);
}

void SV_Statistics::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_observables.size();++i) {
    m_observables[i]->Finalize();
    if (scale!=1.) m_observables[i]->Scale(scale);
  }
}

void SV_Statistics::Output(const std::string & pname) 
{
  ATOOLS::MakeDir((pname).c_str(),448); 
  for (size_t i=0; i<m_observables.size();++i) 
    m_observables[i]->Output((pname+"/"+m_name+ToString(i)+".dat").c_str());
}

Primitive_Observable_Base & 
SV_Statistics::operator+=(const Primitive_Observable_Base & ob)
{
  SV_Statistics * job = ((SV_Statistics*)(&ob));
  for (size_t i(m_observables.size());i<job->m_observables.size();++i) {
    m_observables.push_back(new ATOOLS::Histogram(0,0.0,10.0,10));
    m_observables[i]->SetFills((long int)m_observables[0]->Fills());
  }
  for (size_t i(0);i<job->m_observables.size();++i) {
    (*m_observables[i])+=(*job->m_observables[i]);
  }
  size_t fills(job->m_observables.size()>0 ? 
	       (long int)job->m_observables[0]->Fills() : 0);
  for (size_t i(job->m_observables.size());i<m_observables.size();++i) {
    m_observables[i]->SetFills(fills);
  }
  return *this;
}

#endif
