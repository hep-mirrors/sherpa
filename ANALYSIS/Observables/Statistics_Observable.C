#ifndef Statistics_Observable_H
#define Statistics_Observable_H

#include "Primitive_Observable_Base.H"
#include "MyStrStream.H"
#include <iostream>
#include <iomanip>

namespace ANALYSIS {

  struct Statistics_Data {

    long unsigned int m_xsn, m_cutxsn, m_n, m_cutn;
    double m_wsum, m_cutwsum, m_xssum, m_cutxssum;
    double m_sudsum, m_cutsudsum;
    double m_sum, m_cutsum, m_sum2, m_cutsum2;

    Statistics_Data():
      m_xsn(0), m_cutxsn(0), m_n(0), m_cutn(0), m_wsum(0.0), m_cutwsum(0.0),
      m_xssum(0.0), m_cutxssum(0.0), m_sudsum(0.0), m_cutsudsum(0.0),
      m_sum(0.0), m_cutsum(0.0), m_sum2(0.0), m_cutsum2(0.0) {} 

    void Output(const std::string &pname);

  };// end of struct Statistics_Data

  typedef std::map<std::string,Statistics_Data> Statistics_Map;

  class Statistics_Observable: public Primitive_Observable_Base {  
  protected:

    unsigned long int m_n;
    Statistics_Map m_stats;

  public:

    Statistics_Observable(const std::string &inlistname,int mode=0);

    void Evaluate(const ATOOLS::Blob_List &blobs,
		  double weight=1.0,int ncount=1);

    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &obj);
    Primitive_Observable_Base * Copy() const;

  };// end of class Statistics_Observable

  std::ostream & operator<<(std::ostream & str,const Statistics_Data &sd);

}// end of namespace ANALYSIS

#endif

using namespace ANALYSIS;

DECLARE_GETTER(Statistics_Observable_Getter,"Statistics",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base * 
Statistics_Observable_Getter::operator()
  (const Argument_Matrix &parameters) const
{
  std::string listname="FinalState";
  if (parameters.size()>0 && parameters[0].size()>0) 
    listname=parameters[0][0];
  return new Statistics_Observable(listname);
}

void Statistics_Observable_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list]"; 
}

#include "Primitive_Analysis.H"
#include "Shell_Tools.H"
#include "Run_Parameter.H"

#include <fstream>

using namespace ATOOLS;

Statistics_Observable::Statistics_Observable
(const std::string & listname, int mode)
{
  m_listname=listname;
  m_name="Statistics_Observable_"+m_listname;
  m_type=mode;
  m_splitt_flag=false;
}

Primitive_Observable_Base &Statistics_Observable::operator+=
(const Primitive_Observable_Base &obj)
{
  const Statistics_Observable *vob((const Statistics_Observable*)&obj);
  m_n+=vob->m_n;
  m_stats.insert(vob->m_stats.begin(),vob->m_stats.end());
  return *this;
}

void Statistics_Observable::Evaluate
(const Blob_List &  blobs,double value, int ncount)
{
  Particle_List *pl_cut(p_ana->GetParticleList(m_listname));

  Blob *sp(blobs.FindFirst(btp::Signal_Process));
  std::string key(sp!=NULL?sp->TypeSpec():"<unknown>");

  Statistics_Map::iterator cit(m_stats.find(key));
  if (cit==m_stats.end()) {
    m_stats[key]=Statistics_Data();
    cit=m_stats.find(key);
  }

  Blob_Data_Base *data((*p_ana)["XS_Weight"]);
  double xsweight(data?data->Get<double>():1.0);
  data=(*p_ana)["Sud_Weight"];
  double sudweight(data?data->Get<double>():1.0);
  data=(*p_ana)["Process_Weight"];
  double procweight(data?data->Get<double>():1.0);
  xsweight/=procweight;
  sudweight*=xsweight;
  m_n+=ncount;
  cit->second.m_n+=ncount;
  cit->second.m_xsn+=ncount;
  cit->second.m_wsum+=value;
  cit->second.m_xssum+=xsweight;
  cit->second.m_sum+=sudweight;
  cit->second.m_sum2+=sqr(sudweight);

  if (pl_cut && pl_cut->size()>0) {
    cit->second.m_cutn+=ncount;
    cit->second.m_cutxsn+=ncount;
    cit->second.m_cutwsum+=value;
    cit->second.m_cutxssum+=xsweight;
    cit->second.m_cutsum+=sudweight;
    cit->second.m_cutsum2+=sqr(sudweight);
  }
}

void Statistics_Observable::EndEvaluation(double scale) 
{
}

void Statistics_Observable::Restore(double scale) 
{
}

void Statistics_Observable::Output(const std::string & pname) 
{
  ATOOLS::MakeDir(pname); 
  std::ofstream file((pname+"/"+m_name).c_str());
  file<<std::setw(30)<<std::left<<"# Process Name"<<std::right
      <<" "<<std::setw(9)<<"# evts"
      <<" "<<std::setw(15)<<"<w>"
      <<" "<<std::setw(15)<<"<sigma>"
      <<" "<<std::setw(15)<<"<Delta>"
      <<"   "<<std::setw(9)<<"# evts cut"
      <<" "<<std::setw(15)<<"<w> cut"
      <<" "<<std::setw(15)<<"<sigma> cut"
      <<" "<<std::setw(15)<<"<Delta> cut\n";
  double sig(0.0), dsig(0.0);
  double csig(0.0), cdsig(0.0);
  for (Statistics_Map::const_iterator cit(m_stats.begin());
       cit!=m_stats.end();++cit) {
    file<<std::setw(30)<<std::left<<cit->first<<std::right
	<<" "<<cit->second<<std::endl;
    double wm(cit->second.m_sum/cit->second.m_xsn);
    double w2m(cit->second.m_sum2/cit->second.m_xsn);
    sig+=wm;
    if (cit->second.m_xsn>1) 
      dsig+=(w2m-sqr(wm))/(cit->second.m_xsn-1.0);
    double cwm(cit->second.m_cutsum/cit->second.m_cutxsn);
    double cw2m(cit->second.m_cutsum2/cit->second.m_cutxsn);
    csig+=cwm;
    if (cit->second.m_cutxsn>1)
      cdsig+=(cw2m-sqr(cwm))/(cit->second.m_cutxsn-1.0);
  }
  dsig=sqrt(dsig);
  cdsig=sqrt(cdsig);
  file<<"\n# Total XS      : "<<sig<<" +- ( "
      <<dsig<<" = "<<dsig/sig*100.0<<" % )\n";
  file<<"# Total XS (cut): "<<csig<<" +- ( "
      <<cdsig<<" = "<<cdsig/csig*100.0<<" % )\n\n";
  if (m_listname=="FinalState")
    rpa.gen.SetVariable("TOTAL_CROSS_SECTION",ToString(sig));
  file.close();
}

Primitive_Observable_Base * Statistics_Observable::Copy() const 
{
  return new Statistics_Observable(m_listname,m_type);
}

std::ostream & ANALYSIS::operator<<
  (std::ostream & str,const Statistics_Data &sd) 
{
  str<<std::setw(9)<<sd.m_n
     <<" "<<std::setw(15)<<sd.m_wsum/double(sd.m_xsn)
     <<" "<<std::setw(15)<<sd.m_xssum/double(sd.m_xsn)
     <<" "<<std::setw(15)<<sd.m_sum/sd.m_xssum
     <<"   "<<std::setw(9)<<sd.m_cutn
     <<" "<<std::setw(15)<<sd.m_cutwsum/double(sd.m_cutxsn)
     <<" "<<std::setw(15)<<sd.m_cutxssum/double(sd.m_cutxsn)
     <<" "<<std::setw(15)<<sd.m_cutsum/sd.m_cutxssum;
  return str;
}
