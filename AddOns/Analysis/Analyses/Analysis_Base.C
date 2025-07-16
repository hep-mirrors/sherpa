#include "AddOns/Analysis/Analyses/Analysis_Base.H"

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Analysis_Base::Analysis_Base(const std::string &listname):
  Primitive_Observable_Base(1,0.,1.,100)
{
  m_listname=listname;
  m_name=m_listname;
}

Analysis_Base::~Analysis_Base()
{
  while (m_dists.size()) {
    delete m_dists.back();
    m_dists.pop_back();
  }
  while (m_histos.size()) {
    delete m_histos.back();
    m_histos.pop_back();
  }
}

void Analysis_Base::AddZeroPoint(const double &ntrial,const int &mode)
{
  for (size_t i(0);i<m_dists.size();++i)
    FillDist(i,0.0,0.0,0.0,ntrial,mode);
  for (size_t i(0);i<m_histos.size();++i)
    FillHisto(i,0.0,0.0,ntrial,mode);
}

void Analysis_Base::Evaluate(const ATOOLS::Particle_List &list,
			     double weight, double ncount)
{
  Evaluate(weight,ncount,0);
}

void Analysis_Base::EvaluateNLOcontrib(double weight,double ncount)
{
  Evaluate(weight,ncount,1);
}

void Analysis_Base::FillHisto
(const size_t &i,const double &x,const double &weight,
 const double &ntrial,const int &mode)
{
  if (mode==1) m_histos[i]->InsertMCB(x,weight,ntrial);
  else m_histos[i]->Insert(x,weight,ntrial);
}

void Analysis_Base::FillDist
(const size_t &i,const double &x,const double &y,
 const double &weight,const double &ntrial,const int &mode)
{
  if (mode==1) m_dists[i]->FillMCB(x,y,weight,ntrial);
  else m_dists[i]->Fill(x,y,weight,ntrial);
}

void Analysis_Base::EvaluateNLOevt()
{
  for (size_t i(0);i<m_dists.size();++i) m_dists[i]->FinishMCB();
  for (size_t i(0);i<m_histos.size();++i) m_histos[i]->FinishMCB();
}

Analysis_Object &Analysis_Base::operator+=
(const Analysis_Object &obj)
{
  const Analysis_Base *vob((const Analysis_Base*)&obj);
  for (size_t i(0);i<m_dists.size();++i) 
    *m_dists[i]+=*vob->m_dists[i];
  for (size_t i(0);i<m_histos.size();++i) 
    *m_histos[i]+=*vob->m_histos[i];
  return *this;
}

void Analysis_Base::EndEvaluation(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    m_dists[i]->EndEvaluation(scale);
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->MPISync();
    m_histos[i]->Finalize();
    m_histos[i]->Scale(scale);
  }
}

void Analysis_Base::Restore(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    m_dists[i]->Restore(scale);
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void Analysis_Base::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_dists.size();++i) {
    m_dists[i]->SetName(m_name+"_f"+ToString(i)+"_"+
			m_dists[i]->Name()+".dat");
    m_dists[i]->Output(pname);
  }
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->Output
      (pname+"/"+m_name+"_h"+ToString(i)+"_"+
       m_histos[i]->Name()+".dat");
  }
  msg_Debugging()<<"}\n";
}
