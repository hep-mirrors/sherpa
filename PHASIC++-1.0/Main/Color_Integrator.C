#include "Color_Integrator.H"

#include "Random.H"
#include "Exception.H"
#include "Message.H"
#include <iomanip>

#include <set>

using namespace PHASIC;
using namespace ATOOLS;

std::ostream &PHASIC::operator<<
  (std::ostream &ostr,const std::vector<size_t> &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

std::ostream &PHASIC::operator<<
  (std::ostream &ostr,const std::vector<std::vector<size_t> > &v)
{
  ostr<<"{";
  for (int i(0);i<(int)v.size()-1;++i) ostr<<v[i]<<",";
  if (v.size()>0) ostr<<v.back();
  return ostr<<"}";
}

std::ostream &PHASIC::operator<<(std::ostream &ostr,const Representation &v)
{
  switch (v.Type()) {
  case -1: return ostr<<"|"<<v.J()<<">("<<v.Id()<<")";
  case 0: return ostr<<"|"<<v.J()<<">["<<v.Id()<<"]<"<<v.I()<<"|";
  case 1: return ostr<<"("<<v.Id()<<")<"<<v.I()<<"|";
  }
  return ostr<<"<error>";
}

Color_Integrator::Color_Integrator():
  m_check(false), m_n(0), m_nv(0) {}

Color_Integrator::~Color_Integrator()
{
  while (m_ids.size()>0) {
    delete m_ids.back();
    m_ids.pop_back();
  }
}

bool Color_Integrator::ConstructRepresentations
(const std::vector<size_t> &ids,const std::vector<int> &types)
{
  m_weight=1.0;
  m_confs.clear();
  m_cweights.clear();
  m_asum.clear();
  if (ids.size()!=types.size()) THROW(fatal_error,"Internal error.");
  m_pairs=0;
  m_ids.resize(ids.size());
  int quarks(0);
  for (size_t i(0);i<ids.size();++i) {
    m_ids[i] = new Representation(ids[i],types[i]);
    if (types[i]>=0) {
      m_weight*=3.0;
      if (types[i]>0) m_pairs+=1;
    }
    quarks+=types[i];
  }
  if (quarks!=0) THROW(fatal_error,"Invalid number of quarks.");
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
(std::vector<size_t> ids,std::vector<size_t> perm,
 bool sing,double weight)
{
  if (perm.size()==m_ids.size()) {
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
	if (m_ids[perm[i-1]]->Type()==1 && 
	    m_ids[perm[i]]->Type()==-1) {
	  ++dpairs;
	  if (m_ids[perm[i-1]]->Id()==
	      m_ids[perm[i]]->Id()) weight/=-3.0;
	}
      }
      weight/=pow(-3.0,m_pairs-dpairs);
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
  std::vector<size_t> tids(1,perm.back());
  if (m_ids[perm.back()]->Type()==-1) {
    newstr=true;
    tids.pop_back();
    // find start for next string 
    // -> quark or singlet gluon
    std::vector<size_t> sids(0);
    for (size_t i(0);i<ids.size();++i) {
      switch (m_ids[ids[i]]->Type()) {
      case 1: tids.push_back(ids[i]); break;
      case 0: sids.push_back(ids[i]); break;
      }
    }
    if (tids.empty()) {
      // if new string starts with gluon, 
      // all remaining gluons are singlets
      // pick randomized any to start
      size_t cg(Max(sids.size()-1,
		    (size_t)(sids.size()*ran.Get())));
      tids.push_back(sids[cg]);
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
    std::vector<size_t> pids;
    if (sing) {
      // test for singlet
      if (m_ids[perm.back()]->J()!=last) return 0;
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
	if (m_ids[ids[i]]->Type()!=1 && 
	    m_ids[ids[i]]->J()==last) 
	  pids.push_back(ids[i]);
      }
    }
    if (newstr && ids.size()==1) {
      // last parton has been used to end the string
      // permutation is finished
      // correct for last 1/NC weight -> added in last step
      weight*=-3.0;
      std::vector<size_t> nids;
      int cnc(ConstructConfigurations(nids,perm,sing,weight));
      if (cnc<0) return -1;
      nc+=cnc;
    }
    else {
      // partons left
      perm.push_back(0);
      std::vector<size_t> nids(newstr?ids.size()-2:ids.size()-1);
      for (size_t i(0);i<pids.size();++i) {
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
	int cnc(ConstructConfigurations(nids,perm,sing,weight));
	if (cnc<0) return -1;
	nc+=cnc;
      }  
      perm.pop_back();
    }
  }
  return nc;
}

int Color_Integrator::ConstructConfigurations()
{
  std::vector<size_t> ids(m_ids.size()-1), perm(1);
  size_t fid(0);
  for (;fid<m_ids.size();++fid) 
    // find first quark
    if (m_ids[fid]->Type()==1) break;
  // if no quark is present take any gluon
  if (fid==m_ids.size()) --fid;
  for (size_t i(0);i<=ids.size();++i) {
    // reorder amplitude, starting with a quark
    // the rest is ordered automatically, when 
    // searching for matching colour indices
    if (i>fid) ids[i-fid-1]=i; 
    else if (m_ids.size()-fid-1+i<ids.size())
      ids[m_ids.size()-fid-1+i]=i;  
  }
  perm.back()=fid;
  msg_Debugging()<<"start with ids="<<ids<<" perm="<<perm<<"\n";
  // permute
  int nc(ConstructConfigurations(ids,perm,false,1.0));
  if (nc<0) return -1;
  return nc;
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

bool Color_Integrator::CheckPermutation(const std::vector<size_t> &perm)
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

bool Color_Integrator::GeneratePoint()
{
  if (!m_confs.empty()) {
    size_t l(0), r(m_asum.size()-1), i((l+r)/2);
    double disc(ran.Get()), a(m_asum[i]);
    while (r-l>1) {
      if (disc<a) r=i;
      else l=i;
      i=(l+r)/2;
      a=m_asum[i];
    }
    std::vector<size_t> &conf(m_confs[r]);
    msg_Debugging()<<"selected "<<r<<" from l="
		   <<m_asum[l]<<" < d="<<disc<<" < r="<<m_asum[r]
		   <<" => "<<conf<<"\n";
    for (size_t i(0);i<m_ids.size();++i) {
      m_ids[i]->SetI(conf[2*i]);
      m_ids[i]->SetJ(conf[2*i+1]);
    }
    return true;
  }
  m_orders.clear();
  m_weights.clear();
  msg_Debugging()<<" --- colors --- \n";
  if (!DiceColours()) {
    msg.Error()<<METHOD<<"(): Cannot dice colors. Abort."<<std::endl;
    return false;
  }
  if (msg.LevelIsDebugging()) {
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
    if (ConstructConfigurations()==0) return true;
    if (m_check)
      for (size_t i(0);i<m_orders.size();++i) 
	if (!CheckPermutation(m_orders[i])) return false;
    m_confs.push_back(std::vector<size_t>(2*m_ids.size()));
    for (size_t i(0);i<m_ids.size();++i) {
      m_confs.back()[2*i]=m_ids[i]->I();
      m_confs.back()[2*i+1]=m_ids[i]->J();
    }
    m_cweights.push_back(1.0);
    m_asum.push_back(1.0);
    msg_Debugging()<<"adding "<<m_confs.back()<<"\n";
    ++m_nv;
    if (int(m_nv)%100==0)
      msg_Info()<<std::setw(12)<<m_nv<<" ("<<std::setw(12)<<m_n<<")"
		<<mm_left(27)<<std::flush;
    return true;
  }
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
  return false;
}

bool Color_Integrator::Initialize()
{
  if (!m_confs.empty()) return true;
  m_nv=m_n=0;
  msg_Info()<<METHOD<<"(): Determining configurations ... "<<std::flush;
  if (!AddConfiguration(0)) return false;
  msg_Info()<<m_nv<<" of "<<m_n<<" ( "
	    <<((size_t)(m_nv*1000)/(size_t)(m_n))/10.0<<"% ) "
	    <<"done                           "<<std::endl;
  double sum(0.0), psum(0.0);
  for (size_t i(0);i<m_asum.size();++i)
    sum+=m_asum[i];
  for (size_t i(0);i<m_asum.size();++i) {
    m_asum[i]=(psum+=m_asum[i])/sum;
  msg_Debugging()<<m_asum[i]<<" "<<psum<<" "<<sum<<"\n";
  }
  m_weight*=m_nv/m_n;
  return true;
}

std::vector<int> Color_Integrator::I() const
{
  std::vector<int> is(m_ids.size());
  for (size_t i(0);i<m_ids.size();++i) 
    is[i]=m_ids[i]->I();
  return is;
}

std::vector<int> Color_Integrator::J() const
{
  std::vector<int> js(m_ids.size());
  for (size_t i(0);i<m_ids.size();++i) 
    js[i]=m_ids[i]->J();
  return js;
}
