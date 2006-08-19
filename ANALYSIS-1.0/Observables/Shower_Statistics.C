#include "Shower_Statistics.H"

using namespace ANALYSIS;

DECLARE_GETTER(Shower_Statistics_Getter,"PSStats",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base * 
Shower_Statistics_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new Shower_Statistics(listname);
}

void Shower_Statistics_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"
#include "ISR_Info.H"
#include "MyStrStream.H"
#include "Shell_Tools.H"

#include <fstream>

using namespace ATOOLS;

Shower_Statistics::Shower_Statistics(const std::string & listname, int type)
{
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-200.,400.,200));
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-10000.,10000.,200));
  for (size_t i=0;i<2;++i) m_observables.push_back(new ATOOLS::Histogram(0,-10.,10.,200));
  m_observables.push_back(new ATOOLS::Histogram(0,0.,20000.,200));
  m_observables.push_back(new ATOOLS::Histogram(10,0.01,100.,200));
  m_observables.push_back(new ATOOLS::Histogram(0,0.,1.,200));
  m_observables.push_back(new ATOOLS::Histogram(10,0.01,100.,200));
  m_name  = "Shower_Statistics_";
  m_type  = type;
  m_listname    = listname;
  m_splitt_flag = true;
}

void Shower_Statistics::Evaluate(const Blob_List &  blobs,double weight,int ncount)
{
  Particle_List * pl_cut = p_ana->GetParticleList(m_listname);
  if (pl_cut && pl_cut->size()==0) {
    for (size_t i=0;i<m_observables.size();++i) m_observables[i]->Insert(1.,0.,ncount);
    return;
  }
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    if ((*bit)->Type()==btp::Signal_Process) {
      Blob_Data_Base * info=(*(*bit))["Shower_Info_cms"];
      if (info) {
	std::vector<double> showerinfo=info->Get<std::vector<double> >();
	m_observables[0]->Insert(showerinfo[PDF::iic::E_1],weight,ncount);
	m_observables[2]->Insert(showerinfo[PDF::iic::t_1],weight,ncount);
	m_observables[4]->Insert(showerinfo[PDF::iic::Em_1],weight,ncount);
	m_observables[6]->Insert(showerinfo[PDF::iic::mu_1],weight,ncount);
	m_observables[7]->Insert(showerinfo[PDF::iic::mu_1]/
				 showerinfo[PDF::iic::mu_2],weight,ncount);
      }
      info=(*(*bit))["Shower_Info_lab"];
      if (info) {
	std::vector<double> showerinfo=info->Get<std::vector<double> >();
	m_observables[1]->Insert(showerinfo[PDF::iic::E_1],weight,ncount);
	m_observables[3]->Insert(showerinfo[PDF::iic::t_1],weight,ncount);
	m_observables[5]->Insert(showerinfo[PDF::iic::Em_1],weight,ncount);
	m_observables[8]->Insert(showerinfo[PDF::iic::z_1],weight,ncount);
	m_observables[9]->Insert(showerinfo[PDF::iic::z_1]/
				 showerinfo[PDF::iic::z_2],weight,ncount);
      }
    }
  }
}

Primitive_Observable_Base * Shower_Statistics::Copy() const 
{
  return new Shower_Statistics(m_listname,m_type);
}

void Shower_Statistics::EndEvaluation(double scale) 
{
  for (size_t i=0; i<m_observables.size();++i) {
    m_observables[i]->Finalize();
    if (scale!=1.) m_observables[i]->Scale(scale);
    m_observables[i]->Output();
  }
}

void Shower_Statistics::Output(const std::string & pname) {
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

Primitive_Observable_Base & Shower_Statistics::operator+=(const Primitive_Observable_Base & ob)
{
  if (m_xmin!=ob.Xmin() || m_xmax!=ob.Xmax() || m_nbins!=ob.Nbins()) {
    msg.Error()<<"ERROR: in Shower_Statistics::operator+=  in"<<m_name<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    return *this;
  }

  Shower_Statistics * job = ((Shower_Statistics*)(&ob));

  if (m_observables.size()==job->m_observables.size()) {
    for (size_t i=0; i<m_observables.size();++i) {
      (*m_observables[i])+=(*job->m_observables[i]);
    }
  }
  return *this;
}

