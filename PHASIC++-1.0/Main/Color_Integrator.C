#include "Color_Integrator.H"

#include "Random.H"
#include "Exception.H"
#include "Message.H"
#include "Data_Reader.H"
#include "STL_Tools.H"
#include <iomanip>

#include <set>

using namespace PHASIC;
using namespace ATOOLS;

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

Color_Integrator::Color_Integrator():
  m_lastconf(0), m_alphamode(0), 
  m_check(false), m_iterate(false), m_on(true), 
  m_otfcc(false), m_fincc(true),
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
  m_confs.clear();
  m_cweights.clear();
  m_asum.clear();
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
  msg_Debugging()<<METHOD<<"(): Weight = "<<m_weight<<"\n";
  m_weight*=m_weight;
  return true;
}

size_t Color_Integrator::GenerateIndex()
{
  double rn(3.0*ran.Get());
  for (double disc(1.0);disc<=3.0;++disc)
    if (disc>=rn) return (size_t)disc;
  return std::string::npos;
}

bool Color_Integrator::DiceColours()
{
  for (size_t i(0);i<m_ids.size();++i) {
    if (!m_ids[i]->Act()) continue;
    switch (m_ids[i]->Type()) {
    case -1:
      m_ids[i]->SetJ(GenerateIndex());
      break;
    case 0:
      m_ids[i]->SetI(GenerateIndex());
      m_ids[i]->SetJ(GenerateIndex());
      break;
    case 1:
      m_ids[i]->SetI(GenerateIndex());
      break;
    }
    if (m_ids[i]->I()==std::string::npos ||
	m_ids[i]->J()==std::string::npos) return false;
  }
  return true;
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
    if (msg.LevelIsDebugging()) {
      msg.Out()<<"permutation "<<m_orders.size()<<" -> "<<perm<<" -> ";
      for (size_t i(0);i<perm.size();++i) msg.Out()<<*m_ids[perm[i]];
      msg.Out()<<" ("<<weight<<")\n";
    }
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
	              (size_t)(sids.size()*ran.Get())));
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
    if (newstr && ids.size()==1) {
      // last parton has been used to end the string
      // permutation is finished
      // correct for last 1/NC weight -> added in last step
      weight*=-3.0;
      Idx_Vector nids;
      int cnc(ConstructConfigurations
	      (nids,perm,sing,weight,nexti,one,depth+1));
      if (cnc<0) return -1;
      nc+=cnc;
      if (one && nc>0) return nc;
    }
    else {
      // partons left
      perm.push_back(0);
      Idx_Vector nids(newstr?ids.size()-2:ids.size()-1);
      Idx_Type &i(nexti[depth+1]);
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
      i=0;
      perm.pop_back();
    }
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
  msg_Debugging()<<"sum red = "<<sumr<<", sum green = "
		 <<sumg<<", sum blue = "<<sumb<<"\n";
  return sumr==0 && sumg==0 && sumb==0;
}

bool Color_Integrator::CheckPermutation(const Idx_Vector &perm)
{
  std::set<size_t> all;
  for (size_t i(0);i<m_ids.size();++i) all.insert(m_ids[i]->Id());
  std::set<size_t> checked;
  for (size_t i(0);i<perm.size();++i) {
    // test for doubled indices
    if (checked.find(perm[i])!=checked.end()) {
      msg.Error()<<METHOD<<"(): Permutation "<<perm<<" contains index "
		 <<perm[i]<<" twice. Abort."<<std::endl;
      return false;
    }
    checked.insert(perm[i]);
    std::set<size_t>::iterator ait(all.find(perm[i]));
    // check for invalid index
    if (ait==all.end()) {
      msg.Error()<<METHOD<<"(): Permutation "<<perm
		 <<" contains invalid index "<<perm[i]
		 <<". Abort."<<std::endl;
      return false;      
    }
    all.erase(ait);
  } 
  // check whether all indices occur
  if (all.size()>0) {
    msg.Error()<<METHOD<<"(): Permutation "<<perm
	       <<" does not contain all indices. Abort."<<std::endl;
    return false;      
  }
  msg_Debugging()<<"checked "<<perm<<" -> ok\n";
  return true;
}

bool Color_Integrator::GenerateOrders()
{
  if (msg.LevelIsDebugging()) {
    msg_Debugging()<<" --- colors --- \n";
    for (size_t i(0);i<m_ids.size();++i)
      msg_Debugging()<<i<<" -> "<<*m_ids[i]<<"\n";
  }
  if (!TrivialCheck()) return false;
  msg_Debugging()<<"color sums agree\n";
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
  for (size_t i(0);i<m_ids.size();++i) {
    m_ids[perm[i]]->SetI(i);
    m_ids[perm[i]]->SetJ(i+1);
  }
  m_ids[perm.front()]->SetI(m_ids[perm.back()]->J());
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
  double rn(ran.Get());
  if (rn>=m_over) {
    m_orders.clear();
    m_weights.clear();
    m_over=0.0;
    return false;
  }
  m_over=0.0;
  return true;
}

int Color_Integrator::Dice()
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
  double rn(ran.Get());
  double cmax(m_alphamode>1?m_max:m_cmax);
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

bool Color_Integrator::GeneratePoint(const bool orders)
{
  if (!m_on) return m_valid=true;
  m_fincc=true;
  m_valid=false;
  if (m_alpha.empty() || m_alphamode==0) {
    DicePoint();
    m_cweight=m_weight;
    if (m_confs.empty() || orders) 
      return m_valid=GenerateOrders();
    return m_valid=true;
  }
  if (LookUp()) return m_valid=true;
  while (true) {
    DicePoint();
    if (!GenerateOrders()) {
      if (m_alphamode>1) return false;
      continue;
    }
    switch (Dice()) {
    case 1: return m_valid=true;
    case 0: return false;
    }
  }
  THROW(fatal_error,"Internal error");
  return false;
}

void Color_Integrator::DicePoint()
{
  if (m_confs.empty()) {
    if (!DiceColours()) THROW(fatal_error,"Cannot dice colors");
    return;
  }
  Idx_Vector conf;
  if (m_iterate) {
    if (m_lastconf>=m_confs.size()) 
      THROW(fatal_error,"Index outy of bounds");
    conf=m_confs[m_lastconf++];      
    msg_Debugging()<<"selected "<<m_lastconf-1<<" "
		   <<" => "<<conf<<"\n";
  }
  else {
    size_t l(0), r(m_asum.size()-1), i((l+r)/2);
    double disc(ran.Get()), a(m_asum[i]);
    while (r-l>1) {
      if (disc<a) r=i;
      else l=i;
      i=(l+r)/2;
      a=m_asum[i];
    }
    if (disc<m_asum[l]) r=l;
    conf=m_confs[r];
    msg_Debugging()<<"selected "<<r<<" from l="<<m_asum[l]
		   <<"("<<l<<") < d="<<disc<<" < r="<<m_asum[r]
		   <<"("<<r<<") => "<<conf<<"\n";
  }
  if (m_confs.empty()) THROW(fatal_error,"Internal error");
  for (size_t i(0);i<m_ids.size();++i) {
    m_ids[i]->SetI(conf[2*i]);
    m_ids[i]->SetJ(conf[2*i+1]);
  }
}

bool Color_Integrator::SetConfiguration(const size_t &id)
{
  if (id>=m_confs.size()) THROW(fatal_error,"Index out of bounds");
  Idx_Vector &conf(m_confs[id]);
  msg_Debugging()<<"selected "<<id<<" => "<<conf<<"\n";
  for (size_t i(0);i<m_ids.size();++i) {
    m_ids[i]->SetI(conf[2*i]);
    m_ids[i]->SetJ(conf[2*i+1]);
  }
  return GenerateOrders();
}

bool Color_Integrator::AddConfiguration(const size_t &l)
{
  if (l==m_ids.size()) {
    m_orders.clear();
    m_weights.clear();
    msg_Debugging()<<" --- colors --- \n";
    if (msg.LevelIsDebugging()) {
      for (size_t i(0);i<m_ids.size();++i)
	msg_Debugging()<<i<<" -> "<<*m_ids[i]<<"\n";
    }
    ++m_n;
    if (!TrivialCheck()) return true;
    msg_Debugging()<<"color sums agree\n";
    m_fincc=true;
    if (!NextOrder()) return true;
    if (m_check)
      for (size_t i(0);i<m_orders.size();++i) {
	if (!CheckPermutation(m_orders[i])) return false;
      }
    m_confs.push_back(Idx_Vector(2*m_ids.size()));
    for (size_t i(0);i<m_ids.size();++i) {
      m_confs.back()[2*i]=m_ids[i]->I();
      m_confs.back()[2*i+1]=m_ids[i]->J();
    }
    m_cweights.push_back(1.0);
    m_asum.push_back(1.0);
    msg_Debugging()<<"adding "<<m_confs.back()<<"\n";
    ++m_nv;
    if (int(m_nv)%100==0 && msg.Modifiable())
      msg_Info()<<std::setw(12)<<m_nv<<" ("<<std::setw(12)<<m_n<<")"
		<<mm_left(27)<<std::flush;
    return true;
  }
  if (!m_ids[l]->Act()) {
    m_ids[l]->SetI(0);
    m_ids[l]->SetJ(0);
    AddConfiguration(l+1);
  }
  else {
    switch (m_ids[l]->Type()) {
    case -1: {
      for (size_t cj(1);cj<=3;++cj) {
	m_ids[l]->SetJ(cj);
	AddConfiguration(l+1);
      }
      return true;
    }
    case 0: {
      for (size_t ci(1);ci<=3;++ci) {
	m_ids[l]->SetI(ci);
	for (size_t cj(1);cj<=3;++cj) {
	  m_ids[l]->SetJ(cj);
	  AddConfiguration(l+1);
	}
      }
      return true;
    }
    case 1: {
      for (size_t ci(1);ci<=3;++ci) {
	m_ids[l]->SetI(ci);
	AddConfiguration(l+1);
      }
      return true;
    }
    }
  }
  return false;
}

bool Color_Integrator::Initialize()
{
  m_lastconf=0;
  if (!m_confs.empty() || m_otfcc) return true;
  m_nv=m_n=0;
  msg_Info()<<METHOD<<"(): Determining configurations ... "<<std::flush;
  if (!AddConfiguration(0)) return false;
  msg_Info()<<m_nv<<" of "<<m_n<<" ( "
	    <<((size_t)(m_nv*1000)/(size_t)(m_n))/10.0<<"% ) "
	    <<"done                           "<<std::endl;
  double sum(0.0), psum(0.0);
  for (size_t i(0);i<m_asum.size();++i)
    sum+=m_asum[i];
  for (size_t i(0);i<m_asum.size();++i)
    m_asum[i]=(psum+=m_asum[i])/sum;
  m_cweight=m_weight*=m_nv/m_n;
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
  for (size_t i(0);i<m_alpha.size();++i) {
    sum+=m_alpha[i];
    min=Min(min,m_alpha[i]);
  }
  m_max=sum*Factorial(m_ids.size()-2);
  m_mean=m_max*pow(3.0,m_ids.size());
  double aexp(0.0);
  Data_Reader read(" ",";","!","=");
  if (!read.ReadFromFile(aexp,"CI_ALPHA_EXP")) aexp=0.0;
  else msg_Info()<<METHOD<<"(): Set \\alpha exp "<<aexp<<".\n";
  m_cmax=pow(min,aexp)*pow(m_max,1.0-aexp);
  msg_Tracking()<<METHOD<<"(): m_max = "<<sum<<"*"
		<<Factorial(m_ids.size()-2)<<" = "<<m_max
		<<", m_cmax = "<<m_cmax<<"\n";
}
