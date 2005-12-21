#include "Statistics_Observable.H"

using namespace ANALYSIS;

DECLARE_GETTER(Statistics_Observable_Getter,"Statistics",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base * 
Statistics_Observable_Getter::operator()(const String_Matrix &parameters) const
{
  std::string listname="Analysed";
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new Statistics_Observable(listname);
}

void Statistics_Observable_Getter::PrintInfo(std::ostream &str,
					     const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"
#include "Shell_Tools.H"
#include "Run_Parameter.H"

#include <fstream>

using namespace ATOOLS;

Statistics_Observable::Statistics_Observable(const std::string & listname, int mode)
{
  m_name  = "Statistics_Observable";
  m_type  = mode;
  m_listname    = listname;
  m_splitt_flag = false;
}

void Statistics_Observable::Evaluate(const Blob_List &  blobs,double value, int ncount)
{
  Particle_List * pl     = p_ana->GetParticleList("FinalState");
  Particle_List * pl_cut = p_ana->GetParticleList(m_listname);

  std::string key="unknown";
  unsigned long blcount=0;
  for (Blob_List::const_iterator bit=blobs.begin();bit!=blobs.end();++bit) {
    ++blcount;
    if ((*bit)->Type()==btp::Signal_Process) key  = (*bit)->TypeSpec();
  }

  Statistics_Map::iterator cit=m_signal_process_statistics.find(key);
  if (cit==m_signal_process_statistics.end()) {
    m_signal_process_statistics[key]=Statistics_Data();
    cit=m_signal_process_statistics.find(key);
  }

  double xsecweight=(*p_ana)["XS_Weight"]->Get<double>();
  double sudweight=(*p_ana)["Sud_Weight"]->Get<double>();
  int xsecntrials=(*p_ana)["XS_NumberOfTrials"]->Get<int>();
  m_nevt+=ncount;
  cit->second.nevt+=ncount;
  cit->second.xsnevt+=1;
  cit->second.nblobssum+=blcount;
  cit->second.nplsum+=pl->size();
  cit->second.weightsum+=value;
  cit->second.xssum+=xsecweight;
  cit->second.sudsum+=sudweight;

  if (pl_cut) {
    if (pl_cut->size()>0) {
      cit->second.nevt_cut+=ncount;
      cit->second.nblobssum_cut+=blcount;
      cit->second.nplsum_cut+=pl_cut->size();
      cit->second.weightsum_cut+=value;
    }
  }
}

void Statistics_Observable::EndEvaluation(double scale) {
}

void Statistics_Observable::Output(const std::string & pname) {
  int  mode_dir = 448;
  ATOOLS::MakeDir((pname).c_str(),mode_dir); 
  std::ofstream file((pname+std::string("/")+m_name).c_str());

  double  sum(0.0);
  for (Statistics_Map::const_iterator cit=m_signal_process_statistics.begin();cit!=m_signal_process_statistics.end();++cit) {
    file<<cit->first<<" : "<<cit->second<<std::endl;
    sum+=cit->second.xssum*cit->second.sudsum/sqr(double(cit->second.xsnevt));
  }
  file<<"total xs: "<<sum<<"\n";
  file.close();
}

Primitive_Observable_Base * Statistics_Observable::Copy() const 
{
  return new Statistics_Observable(m_listname,m_type);
}

std::ostream & ANALYSIS::operator<<(std::ostream & str, const Statistics_Data & sd) {
  str<<sd.nevt<<","<<double(sd.nblobssum)/double(sd.nevt)<<","
     <<double(sd.nplsum)/double(sd.nevt)<<","
     <<sd.weightsum/double(sd.nevt)<<", "<<sd.xssum/double(sd.xsnevt)
     <<","<<sd.sudsum/double(sd.xsnevt)<<",    ";
  str<<sd.nevt_cut<<","<<double(sd.nblobssum_cut)/double(sd.nevt_cut)<<","
     <<double(sd.nplsum_cut)/double(sd.nevt_cut)<<","
     <<sd.weightsum/double(sd.nevt_cut);
  return str;
}
