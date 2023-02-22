#include "PHASIC++/Main/Helicity_Integrator.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace PHASIC;
using namespace ATOOLS;

std::ostream &PHASIC::operator<<(std::ostream &str,const hls::scheme &s)
{
  switch (s) {
  case hls::unknown: return str<<"<unknown>";
  case hls::sum: return str<<"sum";
  case hls::sample: return str<<"sample";
  }
  return str<<"<error>";
}

std::map<std::string,std::string> hls::HelicitySchemeTags()
{
  std::map<std::string,std::string> tags;
  tags["UNKNOWN"]=ToString((int)hls::unknown);
  tags["SUM"]=ToString((int)hls::sum);
  tags["SAMPLE"]=ToString((int)hls::sample);
  return tags;
}
 
Helicity_Integrator::Helicity_Integrator():
  m_iter(1), m_on(1) {}

Helicity_Integrator::~Helicity_Integrator() 
{
}

bool Helicity_Integrator::CheckChirs(const Int_Vector &chirs)
{
  size_t p(0), m(0), msv(0);
  Int_Vector q(94,0);
  for (size_t i(0);i<chirs.size();++i) {
    int ch(i<m_nin?-chirs[i]:chirs[i]);
    Flavour fl(i<m_nin?m_flavs[i].Bar():m_flavs[i]);
    if (fl.IsQuark() && !fl.IsMassive()) q[fl.Kfcode()]+=ch;
    if (ch>0) ++p;
    else if (ch<0) ++m;
    else THROW(fatal_error,"Invalid helicities");
    if (fl.IsMassive() || !fl.Strong()) ++msv;
  }
  if (msv) return true;
  for (size_t i(0);i<q.size();++i) 
    if (q[i]!=0) return false;
  return p>1 && m>1;
}

void Helicity_Integrator::Construct(Int_Vector chirs,const size_t i)
{
  if (i==m_chirs.size()) {
    if (CheckChirs(chirs)) {
      size_t id(MakeId(chirs));
      msg_Debugging()<<"adding helicity configuration "
		     <<chirs<<" -> "<<id<<"\n";
      m_weights[id]=1.0;
      ++m_valid;
    }
    return;
  }
  if (chirs[i]!=0) {
    Construct(chirs,i+1);
  }
  else {
    for (int ch(-1);ch<=1;ch+=2) {
      chirs[i]=ch;
      Construct(chirs,i+1);
    }
  }
}

bool Helicity_Integrator::Construct
(const Flavour_Vector &flavs,const int nin)
{
  DEBUG_FUNC(flavs);
  m_nin=nin;
  m_flavs=flavs;
  m_chirs.resize(m_flavs.size());
  m_valid=0;
  m_weights.resize(1<<m_flavs.size(),0.0);
  m_asum.resize(m_weights.size(),0.0);
  m_sum.resize(m_weights.size(),0.0);
  m_sum2.resize(m_weights.size(),0.0);
  m_n.resize(m_weights.size(),0);
  Construct(std::vector<int>(m_chirs.size(),0),0);
  double sum(0.0), asum(0.0);
  for (size_t i(0);i<m_weights.size();++i) sum+=m_weights[i];
  for (size_t i(0);i<m_weights.size();++i) {
    m_weights[i]/=sum;
    m_asum[i]=asum+=m_weights[i];
  }
  m_weight=m_valid;
  msg_Debugging()<<"found "<<m_valid<<" configurations\n";
  return true;
}

void Helicity_Integrator::WriteOut(const std::string &pid)
{
  My_Out_File file(pid+"/HW_"+ToString(m_chirs.size()));
  file.Open();
  file->precision(14);
  msg_Debugging()<<METHOD<<"(): Write {\n";
  for (size_t i(0);i<m_weights.size();++i) {
    *file<<m_weights[i]<<" "<<m_sum[i]<<" "<<m_sum2[i]<<" "<<m_n[i]<<"\n";
    msg_Debugging()<<"  "<<MakeId(i)<<" -> "<<m_weights[i]<<"\n";
  }
  msg_Debugging()<<"}\n";
}

void Helicity_Integrator::ReadIn(const std::string &pid)
{
  My_In_File file(pid+"/HW_"+ToString(m_chirs.size()));
  if (!file.Open()) return;
  file->precision(14);
  msg_Debugging()<<METHOD<<"(): Read {\n";
  double sum(0.0);
  for (size_t i(0);i<m_weights.size();++i) {
    *file>>m_weights[i]>>m_sum[i]>>m_sum2[i]>>m_n[i];
    m_asum[i]=sum+=m_weights[i];
    msg_Debugging()<<"  "<<MakeId(i)<<" -> "<<m_weights[i]<<" ("<<sum<<")\n";
  }
  msg_Debugging()<<"}\n";
}

bool Helicity_Integrator::GeneratePoint()
{
  if (!m_on) return true;
  size_t l(0), r(m_asum.size()-1), i((l+r)/2);
  double disc(ran->Get()), a(m_asum[i]);
  while (r-l>1) {
    if (disc<a) r=i;
    else l=i;
    i=(l+r)/2;
    a=m_asum[i];
  }
  while (r>0 && m_weights[r]==0.0) --r;
  if (disc<m_asum[l]) --r;
  msg_Debugging()<<"selected "<<r<<" -> "<<MakeId(r)<<" from l="
		 <<m_asum[l]<<" < d="<<disc<<" < r="<<m_asum[r]<<"\n";
  m_chirs=MakeId(m_id=r);
  m_new=true;
  return true;
}

double Helicity_Integrator::Weight()
{
  if (!m_on) return 1.0;
  if (m_id>m_weights.size()) THROW(fatal_error,"Invalid identifier");
  return 1.0/(m_valid*m_weights[m_id])*m_weight;
}

void Helicity_Integrator::AddPoint(const double &weight)
{
  if (!m_new) return;
  m_new=false;
  m_sum[m_id]+=weight;
  m_sum2[m_id]+=sqr(weight);
  ++m_n[m_id];
}

void Helicity_Integrator::MPISync()
{
#ifdef USING__MPI
  mpi->Allreduce(&m_sum[0],m_sum.size(),MPI_DOUBLE,MPI_SUM);
  mpi->Allreduce(&m_sum2[0],m_sum2.size(),MPI_DOUBLE,MPI_SUM);
  mpi->Allreduce(&m_n[0],m_n.size(),MPI_DOUBLE,MPI_SUM);
#endif
}

void Helicity_Integrator::Optimize()
{
  double norm(0.0), oldnorm(0.0);
  for (size_t i(0);i<m_weights.size();++i)
    if (m_n[i]==0) return;
  ++m_iter;
  for (size_t i(0);i<m_weights.size();++i) {
    double alpha(m_weights[i]);
    oldnorm+=alpha;
    alpha=m_sum[i]?sqrt(sqrt(alpha)*m_sum2[i]/m_sum[i]):0.;
    m_weights[i]=alpha;
    norm+=alpha;
  }
  norm/=oldnorm;
  oldnorm=0.0;
  for (size_t i(0);i<m_weights.size();++i) {
    if (m_sum2[i]!=0.0) m_weights[i]/=norm;
    m_asum[i]=oldnorm+=m_weights[i];
  }
  if (!IsEqual(oldnorm,1.0)) 
    THROW(fatal_error,"Summation does not agree.");
}

size_t Helicity_Integrator::MakeId(const Int_Vector &ids) const
{
  if (ids.size()!=m_chirs.size()) 
    THROW(fatal_error,"Invalid particle number");
  size_t id(0);
  for (size_t i(0);i<ids.size();++i) 
    if (ids[i]>0) id+=1<<i;
#ifdef DEBUG__CDBCF
  msg_Debugging()<<METHOD<<ids<<" -> "<<id<<"\n";
#endif
  return id;
}

Int_Vector Helicity_Integrator::MakeId(const size_t &id) const
{
  size_t ic(id);
  Int_Vector ids(m_chirs.size(),-1);
  for (size_t i(0);i<ids.size();++i) {
    size_t c(1<<i);
    if (ic&c) {
      ids[i]=1;
      ic-=c;
    }
  }
  if (ic!=0) THROW(fatal_error,"Invalid particle number");
  return ids;
}

