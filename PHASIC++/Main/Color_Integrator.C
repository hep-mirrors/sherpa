#include "PHASIC++/Main/Color_Integrator.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Smart_Pointer.C"
#include <iomanip>
#include <limits>

#include <set>

using namespace PHASIC;
using namespace ATOOLS;

namespace ATOOLS { template class SP(Color_Integrator); }

std::ostream &PHASIC::operator<<(std::ostream &str,const cls::scheme &s)
{
  switch (s) {
  case cls::unknown: return str<<"<unknown>";
  case cls::sum: return str<<"sum";
  case cls::sample: return str<<"sample";
  }
  return str<<"<error>";
}

std::map<std::string,std::string> cls::ColorSchemeTags()
{
  std::map<std::string,std::string> tags;
  tags["UNKNOWN"]=ToString((int)cls::unknown);
  tags["SUM"]=ToString((int)cls::sum);
  tags["SAMPLE"]=ToString((int)cls::sample);
  return tags;
}
 
std::ostream &PHASIC::operator<<(std::ostream &ostr,const Representation &v)
{
  if (v.Act())
    switch (v.Type()) {
    case -1: return ostr<<"|"<<v.J()<<">("<<v.Id()<<")";
    case 0: return ostr<<"|"<<v.J()<<">["<<v.Id()<<"]<"<<v.I()<<"|";
    case 1: return ostr<<"("<<v.Id()<<")<"<<v.I()<<"|";
    }
  else
    switch (v.Type()) {
    case -1: return ostr<<"|"<<v.J()<<">{"<<v.Id()<<"}";
    case 0: return ostr<<"|"<<v.J()<<">{"<<v.Id()<<"}<"<<v.I()<<"|";
    case 1: return ostr<<"{"<<v.Id()<<"}<"<<v.I()<<"|";
    }
  return ostr<<"<error>";
}

Representation::Representation(const size_t &id,
			       const int &type,const int &act): 
  m_id(id), m_i(0), m_j(0), m_type(type), m_act(act)
{
  m_ids=ID(m_id);
}

Color_Integrator::Color_Integrator():
  m_lastconf(0), m_alphamode(0), 
  m_check(false), m_on(true), 
  m_otfcc(false), m_fincc(true), m_nogen(true), m_won(true),
  m_n(0), m_nv(0), m_over(0.0) {}

Color_Integrator::~Color_Integrator()
{
  while (m_ids.size()>0) {
    delete m_ids.back();
    m_ids.pop_back();
  }
}

double Color_Integrator::Factorial(const double &n) const
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

bool Color_Integrator::ConstructRepresentations
(const Idx_Vector &ids,const Int_Vector &types,const Int_Vector &acts)
{
  m_weight=1.0;
  m_otfcc=ids.size()>10;
  if (ids.size()!=types.size()) THROW(fatal_error,"Internal error.");
  m_pairs=0;
  m_ids.resize(ids.size());
  int fermions(0);
  for (size_t i(0);i<ids.size();++i) {
    m_ids[i] = new Representation(ids[i],types[i],acts[i]);
    if (types[i]>=0 && acts[i]>0) m_weight*=3.0;
    if (types[i]>0) m_pairs+=1;
    fermions+=types[i];
  }
  if (fermions!=0) THROW(fatal_error,"Invalid number of fermions.");
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): Weight = "<<m_weight<<"\n";
#endif
  m_weight*=m_weight;
  return true;
}

size_t Color_Integrator::GenerateIndex()
{
  double rn(3.0*ran->Get());
  for (double disc(1.0);disc<=3.0;++disc)
    if (disc>=rn) return (size_t)disc;
  return std::string::npos;
}

bool Color_Integrator::GenerateColours()
{
  Idx_Vector iids, jids;
  for (size_t i(0);i<m_ids.size();++i)
    // collect indices
    if (m_ids[i]->Act()) { 
      if (m_ids[i]->Type()>=0) iids.push_back(i);
      if (m_ids[i]->Type()<=0) jids.push_back(i);
    }
  size_t nr(0), ng(0), nb(0);
  for (size_t i(0);i<iids.size();++i) {
    // select partner
    size_t j(Min(jids.size()-1,(size_t)(ran->Get()*jids.size())));
    // set colours
    size_t idx(GenerateIndex());
    m_ids[iids[i]]->SetI(idx);
    m_ids[jids[j]]->SetJ(idx);
    if (idx==1) ++nr;
    else if (idx==2) ++ng;
    else if (idx==3) ++nb;
    else THROW(fatal_error,"Internal error");
    // remove partner from list
    for (Idx_Vector::iterator jit(jids.begin());
	 jit!=jids.end();++jit) if (*jit==jids[j]) {
	jids.erase(jit);
	break;
      }
  }
  m_weight=pow(3.0,iids.size())*Factorial(iids.size())/
    (Factorial(nr)*Factorial(ng)*Factorial(nb));
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): w = "<<m_weight<<"\n";
#endif
  return true;
}

void Color_Integrator::SetPoint(const Int_Vector &ci,const Int_Vector &cj)
{
  if (ci.size()!=m_ids.size() || cj.size()!=m_ids.size()) 
    THROW(fatal_error,"Invalid number of colours");
  for (size_t i(0);i<m_ids.size();++i) {
    m_ids[i]->SetI(ci[i]);
    m_ids[i]->SetJ(cj[i]);
  }
  size_t nr(0), ng(0), nb(0);
  for (size_t k(0);k<m_ids.size();++k) {
    int idx(m_ids[k]->I());
    if (idx==1) ++nr;
    else if (idx==2) ++ng;
    else if (idx==3) ++nb;
  }
  if (nr==0) ++nr;
  if (ng==0) ++ng;
  if (nb==0) ++nb;
  size_t niids(0);
  for (size_t i(0);i<m_ids.size();++i)
    if (m_ids[i]->Act() && m_ids[i]->Type()>=0) ++niids;
  m_weight=pow(3.0,niids)*Factorial(niids)/
    (Factorial(nr)*Factorial(ng)*Factorial(nb));
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): w = "<<m_weight<<"\n";
#endif
  m_fincc=true;
  m_cweight=m_weight;
  m_valid=m_nogen?true:GenerateOrders();
}

void Color_Integrator::SetPoint(const Cluster_Amplitude *const ampl)
{
  if (ampl->Legs().size()!=m_ids.size()) 
    THROW(fatal_error,"Invalid number of colours");
  CI_Map cmap(ampl->ColorMap());
  cmap[0]=0;
  for (size_t i(0);i<m_ids.size();++i) {
    ColorID cc(ampl->Leg(i)->Col());
    CI_Map::iterator iit(cmap.find(cc.m_i)), jit(cmap.find(cc.m_j));
    if (iit!=cmap.end()) m_ids[i]->SetI(iit->second);
    else m_ids[i]->SetI(cc.m_i);
    if (jit!=cmap.end()) m_ids[i]->SetJ(jit->second);
    else m_ids[i]->SetJ(cc.m_j);
  }
#ifdef DEBUG__BG
  msg_Debugging()<<I()<<"\n";
  msg_Debugging()<<J()<<"\n";
#endif
  int nr(0), ng(0), nb(0), niids(0);
  if (ampl->NewCol()==1) --nr;
  if (ampl->NewCol()==2) --ng;
  if (ampl->NewCol()==3) --nb;
  for (size_t k(0);k<m_ids.size();++k) {
    int idx(m_ids[k]->I());
    if (idx==1) ++nr;
    else if (idx==2) ++ng;
    else if (idx==3) ++nb;
  }
  if (nr==0) ++nr;
  if (ng==0) ++ng;
  if (nb==0) ++nb;
  for (size_t i(0);i<m_ids.size();++i)
    if (m_ids[i]->Act() && m_ids[i]->Type()>=0) ++niids;
  m_weight=pow(3.0,niids)*Factorial(niids)/
    (Factorial(nr)*Factorial(ng)*Factorial(nb));
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): w = "<<m_weight<<"\n";
#endif
  m_fincc=true;
  m_cweight=m_weight;
  m_valid=m_nogen?true:GenerateOrders();
}

int Color_Integrator::ConstructConfigurations
(Idx_Vector ids,Idx_Vector perm,bool sing,
 double weight,Idx_Vector &nexti,bool one,size_t depth)
{
  if (perm.size()==m_ids.size()) {
    ++nexti[depth];
    // last step of permutation 
    if (m_ids[perm.front()]->Type()==0) {
      // pure gluonic -> last i must match first j
      if (m_ids[perm.back()]->I()!=
	  m_ids[perm.front()]->J()) return 0;
    }
    else if (2*m_pairs<perm.size()) {
      // quarks -> 1/NC weight for each indirect pair
      // and for each singlet gluon decaying into quarks
      size_t dpairs(1);
      for (size_t i(1);i<perm.size();++i){
	if (m_ids[perm[i-1]]->Type()>0 && 
	    m_ids[perm[i-1]]->Type()==-m_ids[perm[i]]->Type()) {
	  ++dpairs;
	  if (m_ids[perm[i-1]]->Id()==
	      m_ids[perm[i]]->Id()) weight/=-3.0;
	}
      }
    }
    // get particle indices for permutation
    for (size_t i(0);i<perm.size();++i) 
      perm[i]=m_ids[perm[i]]->Id();
    // add permutation and weight
    m_orders.push_back(perm);
    if (sing) weight/=-3.0;
    m_weights.push_back(weight);
#ifdef DEBUG__BG
    if (msg_LevelIsDebugging()) {
      msg_Out()<<"permutation "<<m_orders.size()<<" -> "<<perm<<" -> ";
      for (size_t i(0);i<perm.size();++i) msg_Out()<<*m_ids[perm[i]];
      msg_Out()<<" ("<<weight<<")\n";
    }
#endif
    return 1;
  }
  bool newstr(false);
  Idx_Vector tids(1,perm.back());
  if (m_ids[perm.back()]->Type()<0) {
    newstr=true;
    tids.pop_back();
    // find start for next string 
    // -> quark or singlet gluon
    Idx_Vector sids(0);
    for (size_t i(0);i<ids.size();++i) {
      switch (m_ids[ids[i]]->Type()) {
      case 0: sids.push_back(ids[i]); break;
      case 1: tids.push_back(ids[i]); break;
      }
    }
    if (tids.empty()) {
      // if new string starts with gluon, 
      // all remaining gluons are singlets
      /*
        // pick randomized any to start
        size_t cg(Max(sids.size()-1,
	              (size_t)(sids.size()*ran->Get())));
        tids.push_back(sids[cg]);
      */
      tids.push_back(sids.front());
      // broadcast that now all gluons 
      // must be in singlet state
      sing=true;
    }
  }
  int nc(0);
  if (newstr) perm.push_back(0);
  for (size_t l(0);l<tids.size();++l) {
    perm.back()=tids[l];
    size_t last(m_ids[perm.back()]->I());
    Idx_Vector pids;
    if (sing) {
      // test for singlet
      if (m_ids[perm.back()]->J()!=last) {
	++nexti[depth];
	return 0;
      }
      else {
	// take only one gluon ordering -> no 1/k!
	for (size_t i(0);i<ids.size();++i) 
	  if (ids[i]!=tids[l]) {
	    pids.push_back(ids[i]);
	    break;
	  }
	// add 1/NC weight
	weight/=-3.0;
      }
    }
    else {
      // find all matching partons
      for (size_t i(0);i<ids.size();++i) {
	if (m_ids[ids[i]]->Type()<=0 && 
	    m_ids[ids[i]]->J()==last) {
	  pids.push_back(ids[i]);
	  // for ew particles consider only one ordering
	  if (last==0) break;
	}
      }
    }
    Idx_Type &i(nexti[depth+1]);
    if (newstr && ids.size()==1) {
      // last parton has been used to end the string
      // permutation is finished
      // correct for last 1/NC weight -> added in last step
      weight*=-3.0;
      Idx_Vector nids;
      if (i==0) {
	int cnc(ConstructConfigurations
		(nids,perm,sing,weight,nexti,one,depth+1));
	if (cnc<0) return -1;
	nc+=cnc;
	if (one && nc>0) return nc;
      }
    }
    else {
      // partons left
      perm.push_back(0);
      Idx_Vector nids(newstr?ids.size()-2:ids.size()-1);
      while (i<pids.size()) {
	// loop over all possible next partons
	size_t shift(0);
	for (size_t j(0);j<ids.size();++j) {
	  // create vector of remaining indices
	  if (j-shift<nids.size()) nids[j-shift]=ids[j];
	  if ((newstr && ids[j]==tids[l]) || 
	      ids[j]==pids[i]) ++shift;
	}
	perm.back()=pids[i];
	// iterate
	int cnc(ConstructConfigurations
		(nids,perm,sing,weight,nexti,one,depth+1));
	if (cnc<0) return -1;
	nc+=cnc;
	if (one && nc>0) return nc;
      }  
    }
    i=0;
    perm.pop_back();
    ++nexti[depth];
  }
  return nc;
}

void Color_Integrator::InitConstruction
(Idx_Vector &ids,Idx_Vector &perm,Idx_Vector &nexti)
{
  perm.resize(1);
  ids.resize(m_ids.size()-1);
  nexti.resize(m_ids.size(),0);
  size_t fid(0);
  for (;fid<m_ids.size();++fid) 
    // find first fermion
    if (m_ids[fid]->Type()>0) break;
  // if no quark is present take any gluon
  if (fid==m_ids.size()) --fid;
  for (size_t i(0);i<=ids.size();++i) {
    // reorder amplitude, starting with a quark
    // the rest is ordered automatically, when 
    // searching for matching colour indices
    if (i>fid) ids[i-fid-1]=i; 
    else if (m_ids.size()-fid-1+i<ids.size())
      ids[m_ids.size()-fid-1+i]=i;  
    nexti[i]=0;
  }
  perm.back()=fid;
}

int Color_Integrator::ConstructConfigurations()
{
  if (m_otfcc) {
    bool one(NextOrder());
    m_fincc=true;
    return one;
  }
  m_orders.clear();
  m_weights.clear();
  // initialize construction
  InitConstruction(m_lastids,m_lastperm,m_nexti);
  // permute
  int nc(ConstructConfigurations
	 (m_lastids,m_lastperm,false,1.0,m_nexti,false,0));
  if (nc<0) return -1;
  return nc;
}

bool Color_Integrator::NextOrder()
{
  if (m_fincc) {
    // initialize construction
    InitConstruction(m_lastids,m_lastperm,m_nexti);
    m_fincc=false;
  }
  m_orders.clear();
  m_weights.clear();
  // permute
  int nc(ConstructConfigurations
	 (m_lastids,m_lastperm,false,1.0,m_nexti,true,0));
  if (nc>0) {
    if (nc>1) THROW(fatal_error,"Internal error");
    return true;
  }
  m_fincc=true;
  return false;
}

bool Color_Integrator::TrivialCheck()
{
  int sumr(0), sumg(0), sumb(0);
  for (size_t i(0);i<m_ids.size();++i) {
    sumr+=(m_ids[i]->I()==1)-(m_ids[i]->J()==1);
    sumg+=(m_ids[i]->I()==2)-(m_ids[i]->J()==2);
    sumb+=(m_ids[i]->I()==3)-(m_ids[i]->J()==3);
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"sum red = "<<sumr<<", sum green = "
		 <<sumg<<", sum blue = "<<sumb<<"\n";
#endif
  return sumr==0 && sumg==0 && sumb==0;
}

void Color_Integrator::SetDecayIds(const std::vector<size_t> &ids,
				   const Int_Vector &types,
				   const Int_Vector &acts)
{
  m_decids.resize(ids.size());
  for (size_t i(0);i<ids.size();++i)
    m_decids[i] = new Representation(ids[i],types[i],acts[i]);
}

bool Color_Integrator::CheckDecays()
{
  for (size_t j(0);j<m_decids.size();++j) {
    int sumr(0), sumg(0), sumb(0);
    const Int_Vector &ids(m_decids[j]->Ids());
    for (size_t i(0);i<ids.size();++i) {
      const Representation *cr(m_ids[ids[i]]);
      sumr+=(cr->I()==1)-(cr->J()==1);
      sumg+=(cr->I()==2)-(cr->J()==2);
      sumb+=(cr->I()==3)-(cr->J()==3);
    }
#ifdef DEBUG__BG
    msg_Debugging()<<"Decay "<<j<<" "<<*m_decids[j]
		   <<": sum red = "<<sumr<<", sum green = "
		   <<sumg<<", sum blue = "<<sumb<<"\n";
#endif
    if (m_decids[j]->Act()) {
      if (m_decids[j]->Type()==1) {
	if (!((sumr==1 && sumg==0 && sumb==0) ||
	      (sumr==0 && sumg==1 && sumb==0) ||
	      (sumr==0 && sumg==0 && sumb==1))) return false;
      }
      else if (m_decids[j]->Type()==-1) {
	if (!((sumr==-1 && sumg==0 && sumb==0) ||
	      (sumr==0 && sumg==-1 && sumb==0) ||
	      (sumr==0 && sumg==0 && sumb==-1))) return false;
      }
      else {
	if (!((sumr==1 && sumg==-1 && sumb==0) ||
	      (sumr==1 && sumg==0 && sumb==-1) ||
	      (sumr==-1 && sumg==1 && sumb==0) ||
	      (sumr==-1 && sumg==0 && sumb==1) ||
	      (sumr==0 && sumg==1 && sumb==-1) ||
	      (sumr==0 && sumg==-1 && sumb==1) ||
	      (sumr==0 && sumg==0 && sumb==0))) return false;
      }
    }
    else {
      if (!(sumr==0 && sumg==0 && sumb==0)) return false;
    }
  }
  return true;
}

bool Color_Integrator::CheckPermutation(const Idx_Vector &perm)
{
  std::set<size_t> all;
  for (size_t i(0);i<m_ids.size();++i) all.insert(m_ids[i]->Id());
  std::set<size_t> checked;
  for (size_t i(0);i<perm.size();++i) {
    // test for doubled indices
    if (checked.find(perm[i])!=checked.end()) {
      msg_Error()<<METHOD<<"(): Permutation "<<perm<<" contains index "
		 <<perm[i]<<" twice. Abort."<<std::endl;
      return false;
    }
    checked.insert(perm[i]);
    std::set<size_t>::iterator ait(all.find(perm[i]));
    // check for invalid index
    if (ait==all.end()) {
      msg_Error()<<METHOD<<"(): Permutation "<<perm
		 <<" contains invalid index "<<perm[i]
		 <<". Abort."<<std::endl;
      return false;      
    }
    all.erase(ait);
  } 
  // check whether all indices occur
  if (all.size()>0) {
    msg_Error()<<METHOD<<"(): Permutation "<<perm
	       <<" does not contain all indices. Abort."<<std::endl;
    return false;      
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"checked "<<perm<<" -> ok\n";
#endif
  return true;
}

bool Color_Integrator::GenerateOrders()
{
#ifdef DEBUG__BG
  if (msg_LevelIsDebugging()) {
    msg_Debugging()<<" --- colors --- \n";
    for (size_t i(0);i<m_ids.size();++i)
      msg_Debugging()<<i<<" -> "<<*m_ids[i]<<"\n";
  }
#endif
  if (!TrivialCheck()) return false;
#ifdef DEBUG__BG
  msg_Debugging()<<"color sums agree\n";
#endif
  if (ConstructConfigurations()==0) return false;
  if (m_check)
    for (size_t i(0);i<m_orders.size();++i)
      if (!CheckPermutation(m_orders[i])) return false;
  return true;
}

bool Color_Integrator::GenerateType(const size_t &type,
				    const bool orders)
{
  if (type>=m_ids.size()-1) return false;
  Idx_Vector perm(m_ids.size());
  for (size_t i(0);i<perm.size();++i) perm[i]=i;
  for (size_t i(1);i<=type;++i) 
    std::swap<Idx_Type>(perm[i],perm[i+1]);
  m_weight=1.0;
  for (size_t i(0);i<m_ids.size();++i) {
    m_weight*=3.0;
    m_ids[perm[i]]->SetI(i);
    m_ids[perm[i]]->SetJ(i+1);
  }
  m_ids[perm.front()]->SetI(m_ids[perm.back()]->J());
  m_cweight=m_weight*=m_weight;
  if (orders) return GenerateOrders();
  return true;
}

size_t Color_Integrator::IdentifyType(const Idx_Vector &perm) const
{
  size_t zero(0), one(0);
  for (;zero<perm.size();++zero) if (perm[zero]==0) break;
  Idx_Vector rp(perm.size());
  for (size_t i(0);i<perm.size();++i)
    rp[i]=i+zero<rp.size()?perm[i+zero]:perm[i+zero-rp.size()];
  for (;one<perm.size();++one) if (rp[one]==1) break;
  return one-1;
}

bool Color_Integrator::LookUp()
{
  if (m_over==0.0) return false;
  if (m_over>1.0) {
    m_over-=1.0;
    return true;
  }
  double rn(ran->Get());
  if (rn>=m_over) {
    m_orders.clear();
    m_weights.clear();
    m_over=0.0;
    return false;
  }
  m_over=0.0;
  return true;
}

int Color_Integrator::Generate()
{
  double weight(0.0);
  if (m_otfcc) {
    while (NextOrder()) {
      size_t type(IdentifyType(m_orders.front()));
      weight+=m_alpha[type];
    }
    m_fincc=true;
  }
  else {
    for (size_t i(0);i<m_orders.size();++i) {
      size_t type(IdentifyType(m_orders[i]));
      weight+=m_alpha[type];
    }
  }
  double rn(ran->Get());
  double cmax(m_alphamode>1?m_max:m_mean/m_weight*m_cmax);
  m_over=Max(0.0,weight/cmax-1.0);
  msg_Debugging()<<METHOD<<"(): amode = "<<m_alphamode<<", rn = "
		 <<rn<<", w = "<<weight<<"/"<<cmax<<" = "<<(weight/cmax)
		 <<", m_over = "<<m_over<<"\n";
  if (m_over==0.0 && weight<rn*cmax) {
    m_orders.clear();
    m_weights.clear();
    if (m_alphamode>1) return 0;
    return -1;
  }
  if (m_alphamode==1) m_cweight=m_mean/weight;
  else m_cweight=m_weight*m_max/weight;
  return 1;
}

bool Color_Integrator::GeneratePoint()
{
  if (!m_on) return m_valid=true;
  m_fincc=true;
  m_valid=false;
  if (m_alpha.empty() || m_alphamode==0) {
    GenerateColours();
    m_cweight=m_weight;
    if (!CheckDecays()) return false;
    return m_valid=m_nogen?true:GenerateOrders();
  }
  if (LookUp()) return m_valid=true;
  while (true) {
    GenerateColours();
    if (!GenerateOrders()) {
      if (m_alphamode>1) return false;
      continue;
    }
    switch (Generate()) {
    case 1: return m_valid=true;
    case 0: return false;
    }
  }
  THROW(fatal_error,"Internal error");
  return false;
}

bool Color_Integrator::Initialize()
{
  return true;
}

void Color_Integrator::SetI(const Int_Vector &i)
{
  for (size_t k(0);k<m_ids.size();++k) 
    m_ids[k]->SetI(i[k]);
}

void Color_Integrator::SetJ(const Int_Vector &j)
{
  for (size_t k(0);k<m_ids.size();++k) 
    m_ids[k]->SetJ(j[k]);
}

Int_Vector Color_Integrator::I() const
{
  Int_Vector is(m_ids.size());
  for (size_t i(0);i<m_ids.size();++i) 
    is[i]=m_ids[i]->I();
  return is;
}

Int_Vector Color_Integrator::J() const
{
  Int_Vector js(m_ids.size());
  for (size_t i(0);i<m_ids.size();++i) 
    js[i]=m_ids[i]->J();
  return js;
}

void Color_Integrator::SetAlpha(const Double_Vector &alpha)
{
  m_alpha=alpha;
  double sum(0.0);
  double min(std::numeric_limits<double>::max());
  double max(0.0);
  for (size_t i(0);i<m_alpha.size();++i) {
    sum+=m_alpha[i];
    min=Min(min,m_alpha[i]);
    max=Max(max,m_alpha[i]);
  }
  m_max=sum*Factorial(m_ids.size()-2);
  m_mean=m_max*pow(3.0,m_ids.size());
  double aexp(0.0);
  Default_Reader reader;
  if (reader.Read(aexp,"CI_ALPHA_EXP",0.0)) {
    msg_Info()<<METHOD<<"(): Set \\alpha exp "<<aexp<<"."<<std::endl;
  }
  m_cmax=pow(max/min,aexp);
  msg_Tracking()<<METHOD<<"(): m_max = "<<sum<<"*"
		<<Factorial(m_ids.size()-2)<<" = "<<m_max
		<<", m_cmax = "<<m_cmax<<"\n";
}
