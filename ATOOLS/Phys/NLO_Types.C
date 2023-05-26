#include "ATOOLS/Phys/NLO_Types.H"

#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const nlo_type::code &c)
{
  std::string out="";
  if (c&nlo_type::born) out+="B";
  if (c&nlo_type::loop) out+="V";
  if (c&nlo_type::vsub) out+="I";
  if (c&nlo_type::real) out+="R";
  if (c&nlo_type::rsub) out+="S";
  if (c&nlo_type::rvirt) out+="E";
  if (c&nlo_type::realreal) out+="W";
  return str<<out;
}

std::istream &ATOOLS::operator>>(std::istream &str,nlo_type::code &c)
{
  std::string tag;
  str>>tag;
  c=nlo_type::lo;
  if (tag.find('B')!=std::string::npos) c|=nlo_type::born;
  if (tag.find('V')!=std::string::npos) c|=nlo_type::loop;
  if (tag.find('I')!=std::string::npos) c|=nlo_type::vsub;
  if (tag.find('R')!=std::string::npos) c|=nlo_type::real;
  if (tag.find('S')!=std::string::npos) c|=nlo_type::rsub;
  if (tag.find('E')!=std::string::npos) c|=nlo_type::rvirt;
  if (tag.find('W')!=std::string::npos) c|=nlo_type::realreal;
  return str;
}

std::ostream & ATOOLS::operator<<(std::ostream & s,
                                  const ATOOLS::asscontrib::type & at)
{
  if (at==asscontrib::none)  s<<"none";
  if (at&asscontrib::EW)     s<<"EW";
  if (at&asscontrib::LO1)    s<<"LO1";
  if (at&asscontrib::LO2)    s<<"LO2";
  if (at&asscontrib::LO3)    s<<"LO3";
  return s;
}

std::istream &ATOOLS::operator>>(std::istream &s,ATOOLS::asscontrib::type &at)
{
  std::string tag;
  getline(s,tag);
  at=asscontrib::none;
  if (tag.find("EW")!=std::string::npos)  at|=asscontrib::EW;
  if (tag.find("LO1")!=std::string::npos) at|=asscontrib::LO1;
  if (tag.find("LO2")!=std::string::npos) at|=asscontrib::LO2;
  if (tag.find("LO3")!=std::string::npos) at|=asscontrib::LO3;
  return s;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const nlo_mode::code &c)
{
  if      (c==nlo_mode::none)       return str<<"none";
  else if (c==nlo_mode::fixedorder) return str<<"fixedorder";
  else if (c==nlo_mode::powheg)     return str<<"powheg";
  else if (c==nlo_mode::mcatnlo)    return str<<"mcatnlo";
  else if (c==nlo_mode::yfs)    return str<<"YFS";
  return str<<"unknown";
}

std::istream &ATOOLS::operator>>(std::istream &str,nlo_mode::code &c)
{
  std::string tag;
  str>>tag;
  c=nlo_mode::none;
  if      (tag.find("None")!=std::string::npos)        c=nlo_mode::none;
  else if (tag.find("Fixed_Order")!=std::string::npos) c=nlo_mode::fixedorder;
  else if (tag.find("1")!=std::string::npos)           c=nlo_mode::fixedorder;
  else if (tag.find("MC@NLO")!=std::string::npos)      c=nlo_mode::mcatnlo;
  else if (tag.find("3")!=std::string::npos)           c=nlo_mode::mcatnlo;
  else if (tag.find("YFS")!=std::string::npos)         c=nlo_mode::yfs;
  else if (tag.find("4")!=std::string::npos)           c=nlo_mode::yfs;
  else                                                 c=nlo_mode::unknown;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const vtype::type &v)
{
  std::string out="";
  if (v&vtype::devidedByBorn) out+="B";
  if (v&vtype::includesI)     out+="I";
  if (v&vtype::needsCoupling) out+="C";
  return str<<out;
}

std::istream &ATOOLS::operator>>(std::istream &str,vtype::type &v)
{
  std::string tag;
  str>>tag;
  v=vtype::none;
  if (tag.find('B')!=std::string::npos) v|=vtype::devidedByBorn;
  if (tag.find('I')!=std::string::npos) v|=vtype::includesI;
  if (tag.find('C')!=std::string::npos) v|=vtype::needsCoupling;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const cs_itype::type &c)
{
  std::string out="";
  if (c&cs_itype::I) out+="I";
  if (c&cs_itype::K) out+="K";
  if (c&cs_itype::P) out+="P";
  return str<<out;
}

std::istream &ATOOLS::operator>>(std::istream &str,cs_itype::type &c)
{
  std::string tag;
  str>>tag;
  c=cs_itype::none;
  if (tag.find('I')!=std::string::npos) c|=cs_itype::I;
  if (tag.find('K')!=std::string::npos) c|=cs_itype::K;
  if (tag.find('P')!=std::string::npos) c|=cs_itype::P;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,const cs_kcontrib::type &c)
{
  std::string out="";
  if (c&cs_kcontrib::Kb)  out+="B";
  if (c&cs_kcontrib::KFS) out+="S";
  if (c&cs_kcontrib::t)   out+="G";
  if (c&cs_kcontrib::Kt)  out+="T";
  return str<<out;
}

std::istream &ATOOLS::operator>>(std::istream &str,cs_kcontrib::type &c)
{
  std::string tag;
  str>>tag;
  c=cs_kcontrib::none;
  if (tag.find('B')!=std::string::npos) c|=cs_kcontrib::Kb;
  if (tag.find('S')!=std::string::npos) c|=cs_kcontrib::KFS;
  if (tag.find('G')!=std::string::npos) c|=cs_kcontrib::t;
  if (tag.find('T')!=std::string::npos) c|=cs_kcontrib::Kt;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const sbt::subtype &st)
{
  if      (st==sbt::none) return ostr<<"NONE";
  else if (st==sbt::qcd)  return ostr<<"QCD";
  else if (st==sbt::qed)  return ostr<<"QED";
  else if (st&sbt::qed && st&sbt::qed)  return ostr<<"QCD|QED";
  return ostr<<"UNKNOWN";
}

std::istream &ATOOLS::operator>>(std::istream &str,sbt::subtype &st)
{
  std::string tag;
  str>>tag;
  st=sbt::none;
  if (tag.find("QCD")!=std::string::npos) st|=sbt::qcd;
  if (tag.find("QED")!=std::string::npos) st|=sbt::qed;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const subscheme::code &ss)
{
  if      (ss==subscheme::CS)   return ostr<<"CS";
  else if (ss==subscheme::Dire) return ostr<<"Dire";
  else if (ss==subscheme::CSS)  return ostr<<"CSS";
  return ostr<<"UNKNOWN";
}

std::istream &ATOOLS::operator>>(std::istream &str,subscheme::code &ss)
{
  std::string tag;
  str>>tag;
  ss=subscheme::CS;
  if (tag.find("1")!=std::string::npos)    ss=subscheme::Dire;
  if (tag.find("Dire")!=std::string::npos) ss=subscheme::Dire;
  if (tag.find("2")!=std::string::npos)    ss=subscheme::CSS;
  if (tag.find("CSS")!=std::string::npos)  ss=subscheme::CSS;
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const dpt::dipoletype &dt)
{
  if      (dt==dpt::none) return ostr<<"NONE";
  else if (dt==dpt::f_f)  return ostr<<"FF";
  else if (dt==dpt::f_fm) return ostr<<"FFm";
  else if (dt==dpt::f_i)  return ostr<<"FI";
  else if (dt==dpt::f_im) return ostr<<"FIm";
  else if (dt==dpt::i_f)  return ostr<<"IF";
  else if (dt==dpt::i_fm) return ostr<<"IFm";
  else if (dt==dpt::i_i)  return ostr<<"II";
  return ostr<<"UNKNOWN";
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const spt::splittingtype &st)
{
  if      (st==spt::none)  return ostr<<"NONE";
  else if (st==spt::q2qg)  return ostr<<"q->qg";
  else if (st==spt::q2gq)  return ostr<<"q->gq";
  else if (st==spt::g2qq)  return ostr<<"g->qq";
  else if (st==spt::g2gg)  return ostr<<"g->gg";
  else if (st==spt::s2sg)  return ostr<<"s->sg";
  else if (st==spt::s2gs)  return ostr<<"s->gs";
  else if (st==spt::G2Gg)  return ostr<<"G->Gg";
  else if (st==spt::G2gG)  return ostr<<"G->gG";
  else if (st==spt::V2Vg)  return ostr<<"V->Vg";
  else if (st==spt::V2gV)  return ostr<<"V->gV";
  return ostr<<"UNKNOWN";
}

std::ostream &ATOOLS::operator<<(std::ostream &ostr,const ist::itype &it)
{
  if      (it==ist::none)  return ostr<<"NONE";
  else if (it==ist::q)     return ostr<<"q";
  else if (it==ist::g)     return ostr<<"g";
  else if (it==ist::Q)     return ostr<<"Q";
  else if (it==ist::V)     return ostr<<"V";
  else if (it==ist::sQ)    return ostr<<"sQ";
  else if (it==ist::sG)    return ostr<<"sG";
  return ostr<<"UNKNOWN";
}

