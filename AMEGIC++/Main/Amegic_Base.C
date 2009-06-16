#include "AMEGIC++/Main/Amegic_Base.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace PHASIC;

Amegic_Base::Amegic_Base() : p_pinfo(0)
{}

Amegic_Base::~Amegic_Base()
{
  if (p_pinfo) delete p_pinfo;
}

Pol_Info AMEGIC::Amegic_Base::ExtractPolInfo(const PHASIC::Subprocess_Info &spi)
{
  Pol_Info pl(spi.m_fl); 
  if (spi.m_id!="") {
    pl.pol_type='d';
    pl.type[0]=spi.m_id[0];
    return pl;
  }
  std::string pn(spi.m_pol);
  if (pn=="") return pl;
  char pc(' '), pp(' ');
  double pd(0.0), angle(0.0);
  size_t lh(pn.find('l'));
  if(lh!=std::string::npos) {
    pc='l';
    pp='+';
    std::string ha = pn.substr(lh+1,pn.length()-lh);
    angle=ToType<double>(ha);
  }
  else {
    pp=pn[pn.length()-1];
    pc=pn[pn.length()-2];
    pn.erase(pn.length()-2,2);
  }
  pd=ToType<double>(pn);	  
#ifdef Explicit_Pols
  int t1=mt::p_m, t2=mt::p_p;
  if (pc=='l') { t1=mt::p_l0; t2=mt::p_l1; }
  if (pc!=' ') pl.pol_type = pc;
  pl.angle    = angle;
  int type;
  switch(pp){
  case '-' : type = t1;     break;
  case '+' : type = t2;     break;
  case '0' : type = mt::p_l;break;
  default  : type = t1;
  }
  if(!spi.m_fl.IsTensor()){
    int tf[3] = {t1,t2,mt::p_l};
    if (ATOOLS::IsZero(pd-1.)) {
      pl.type[0]=type;
      pl.factor[0]=pl.num;
      pl.num=1;
    }
    else{
      for (int j=0;j<pl.num;j++) {
	pl.type[j]=tf[j];
	if(pl.type[j]==type)  pl.factor[j] = 1.+pd*(pl.num-1.);
	else  pl.factor[j] = 1.-pd;
      }
    }
  }
#endif
  return pl;
}

void AMEGIC::Amegic_Base::TranslateDecay(Process_Tags &info,const PHASIC::Subprocess_Info &spi)
{
  ATOOLS::Flavour_Vector ffl(spi.m_ps.size());
  std::vector<Pol_Info> fpl(spi.m_ps.size());
  for (size_t i(0);i<spi.m_ps.size();++i) {
    ffl[i]=spi.m_ps[i].m_fl;
    fpl[i]=ExtractPolInfo(spi.m_ps[i]);
  }
  Process_Tags *dec(info.FindDM(spi.m_id));
  dec->AddSubList(spi.m_ps.size(),&ffl.front(),&fpl.front());
  dec->m_maxqcdjets = spi.m_nmax;
  dec->m_zerowidth = spi.m_osf;
  for (size_t i(0);i<spi.m_ps.size();++i)
    if (spi.m_ps[i].m_id!="") TranslateDecay(info,spi.m_ps[i]);
}

Process_Tags *AMEGIC::Amegic_Base::Translate(const PHASIC::Process_Info &pi)
{
  Subprocess_Info ii(pi.m_ii), fi(pi.m_fi);
  ATOOLS::Flavour_Vector ffl(fi.m_ps.size());
  std::vector<Pol_Info> fpl(fi.m_ps.size());
  for (size_t i(0);i<fi.m_ps.size();++i) {
    ffl[i]=fi.m_ps[i].m_fl;
    fpl[i]=ExtractPolInfo(fi.m_ps[i]);
  }
  Process_Tags *info(new Process_Tags(NULL,NULL));
  info->AddSubList(fi.m_ps.size(),&ffl.front(),&fpl.front());
  info->m_maxqcdjets = fi.m_nmax;
  info->m_zerowidth = fi.m_osf;
  for (size_t i(0);i<fi.m_ps.size();++i)
    if (fi.m_ps[i].m_id!="") TranslateDecay(*info,fi.m_ps[i]);
  if (!info->CheckCompleteness()) THROW(fatal_error,"Missing decay");
  return info;
}


