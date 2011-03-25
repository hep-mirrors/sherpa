#include "ATOOLS/Phys/RB_Data.H"

#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

RB_Key::RB_Key(const int type,const ATOOLS::Flavour &flij,
               const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj):
  m_type(type), m_flij(flij), m_fli(fli), m_flj(flj)
{
  if (m_type&1) std::swap(m_flij,m_fli);
}

RB_Key::RB_Key(const NLO_subevt *sub): m_type(0)
{
  if (sub->m_i<2) m_type|=1;
  if (sub->m_j<2) THROW(fatal_error,"Emitted particle in initial state.");
  if (sub->m_k<2) m_type|=2;
  for (size_t i(0);i<sub->m_n;++i)
    if (sub->p_id[i]&(1<<sub->m_i)) {
      m_flij=sub->p_fl[i];
      break;
    }
  m_fli=sub->p_real->p_fl[sub->m_i];
  m_flj=sub->p_real->p_fl[sub->m_j];
}

bool RB_Key::operator<(const RB_Key &rb) const
{
  if (m_type<rb.m_type) return true;
  if (m_type>rb.m_type) return false;
  if (m_flij<rb.m_flij) return true;
  else if (m_flij==rb.m_flij) {
    if (m_fli<rb.m_fli) return true;
    else if (m_fli==rb.m_fli)
      return m_flj<rb.m_flj;
    return false;
  }
  return false;
}

RB_Data::RB_Data():
  m_sum(0.0), m_max(0.0), m_n(0.0),
  m_rbint(0.0), m_xssum(0.0),
  m_ktres(0.0), m_bmax(0.0),
  p_wh(NULL)
{
  InitWeightHisto();
}

RB_Data::~RB_Data()
{
  if (p_wh!=NULL) delete p_wh;
}

RB_Data &RB_Data::operator=(const RB_Data &rb)
{
  m_sum=rb.m_sum;
  m_max=rb.m_max;
  m_n=rb.m_n;
  m_rbint=rb.m_rbint;
  m_xssum=rb.m_xssum;
  m_ktres=rb.m_ktres;
  m_bmax=rb.m_bmax;
  *p_wh=*rb.p_wh;
  return *this;
}

void RB_Data::AddPoint(const double &rb,const double &xs,const double &b)
{
  if (msg_LevelIsDebugging()) {
    size_t n(m_n+1);
    msg_Out()<<om::bold<<om::blink<<om::brown<<"RB:   "<<rb<<"   in "
        <<n<<(n%10==1?"st":(n%10==2?"nd":(n%10==3?"rd":"th")))<<" event.\n";
  }
  if (IsBad(rb)) {
    msg_Error()<<METHOD<<"(){\n  RB ratio is bad : "<<rb
        <<"  Point not considered.\n}\n";
    return;
  }
  m_n+=1.0;
  m_sum+=rb;
  m_max=Max(m_max,rb);
  m_bmax=Max(m_bmax,b);
  m_rbint+=rb*xs;
  m_xssum+=xs;
  if (p_wh) {
    if(rb!=0.0) p_wh->Insert(rb,xs);
    else p_wh->SetFills(p_wh->Fills()+1);
  }
}

void RB_Data::InitWeightHisto() 
{
  if (p_wh) delete p_wh;
  p_wh = new Histogram(10,1.0e-4,1.0e6,1000);
}

void RB_Data::ReadInHisto(const std::string &fname)
{
  if (!FileExists(fname+".rbh")) return;
  if (p_wh) delete p_wh;
  p_wh = new Histogram(fname+".rbh");
}

void RB_Data::WriteOutHisto(const std::string &fname) const
{
  if (p_wh) p_wh->Output(fname+".rbh");	
}

double RB_Data::GetMaxEps(const double &eps)
{
  if (p_wh==NULL) return m_max;
  double xstot(0.0), cutxs(0.0);
  for (int i(p_wh->Nbin()+1);i>0;--i)
    xstot+=p_wh->Value(i)*pow(10.0,p_wh->Xmin()+(i-0.5)*p_wh->BinSize());
  for (int i(p_wh->Nbin()+1);i>0;--i) {
    cutxs+=p_wh->Value(i)*pow(10.0,p_wh->Xmin()+(i-0.5)*p_wh->BinSize());
    if (cutxs>xstot*eps)
      return Min(pow(10.0,p_wh->Xmin()+i*p_wh->BinSize()),m_max);
  }
  return m_max;
}

namespace ATOOLS {

  std::ostream &operator<<(std::ostream &str,const RB_Key &rbk)
  {
    return str<<rbk.m_type<<" "<<((long int)rbk.m_flij)
	      <<" "<<((long int)rbk.m_fli)
	      <<" "<<((long int)rbk.m_flj);
  }

  std::istream &operator>>(std::istream &str,RB_Key &rbk)
  {
    long int ij, i, j;
    str>>rbk.m_type>>ij>>i>>j;
    rbk.m_flij=Flavour((kf_code)(abs(ij)),ij<0);
    rbk.m_fli=Flavour((kf_code)(abs(i)),i<0);
    rbk.m_flj=Flavour((kf_code)(abs(j)),j<0);
    return str;
  }

  std::ostream &operator<<(std::ostream &str,const RB_Data &rbd)
  {
    return str<<rbd.m_max<<" "<<rbd.m_sum<<" "<<rbd.m_n
	      <<" "<<rbd.m_rbint<<" "<<rbd.m_xssum
	      <<" "<<rbd.m_ktres<<" "<<rbd.m_bmax;
  }

  std::istream &operator>>(std::istream &str,RB_Data &rbd)
  {
    return str>>rbd.m_max>>rbd.m_sum>>rbd.m_n
	      >>rbd.m_rbint>>rbd.m_xssum
	      >>rbd.m_ktres>>rbd.m_bmax;
  }

}
