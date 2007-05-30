#include "N_Parton_CDBG_T.H"

#include "Message.H"
#include "Random.H"
#include "MyStrStream.H"
#include "Color.H"
#include "STL_Tools.H"
#include <iomanip>

using namespace EXTRAXS;
using namespace ATOOLS;

N_Parton_CDBG_T::N_Parton_CDBG_T(const size_t nin,const size_t nout,
				 const std::vector<Flavour> &flavs,
				 const std::vector<std::string> &models):
  m_nin(nin), m_nout(nout), m_flavs(flavs), m_mode(0), m_tests(0)
{ 
  OrderFlavours();
  std::vector<int> types(m_nin+m_nout,0);
  m_ampl.Construct(m_flavs,models);
  m_ng=0;
  Int_Vector ic(m_flavs.size(),0), jc(m_flavs.size(),0);
  for (size_t lc(0), k(0);k<m_flavs.size();++k) {
    if (m_flavs[k].IsGluon()) {
      jc[k]=(ic[k]=++lc)-1;
      ++m_ng;
    }
    else if (m_flavs[k].IsAnti()) jc[k]=lc; 
    else ic[k]=++lc;
  }
  if (m_flavs.front().IsGluon()) jc.front()=m_flavs.size();
  m_ampl.SetColors(ic,jc);
  m_moms.resize(m_nin+m_nout);
  std::map<Flavour,size_t> fc;
  for (size_t i(2);i<flavs.size();++i) {
    std::map<Flavour,size_t>::iterator fit(fc.find(flavs[i]));
    if (fit==fc.end()) {
      fc[flavs[i]]=0;
      fit=fc.find(flavs[i]);
    }
    ++fit->second;
  }
  m_sf=1.0;
  m_norm=pow(2.0,m_ng);
  msg_Debugging()<<METHOD<<"(): Construct symmetry factor {\n";
  msg_Debugging()<<"  Final state:\n";
  for (std::map<Flavour,size_t>::const_iterator fit(fc.begin());
       fit!=fc.end();++fit) {
    msg_Debugging()<<"  "<<std::setw(2)<<fit->second<<" "
		   <<std::setw(15)<<fit->first<<" -> "
		   <<std::setw(12)<<Factorial(fit->second)<<"\n";
    m_sf*=Factorial(fit->second);
  }
  msg_Debugging()<<"  Initial state:\n";
  for (size_t i(0);i<2;++i) {
    if (flavs[i].IsGluon()) {
      msg_Debugging()<<"     "<<std::setw(15)<<flavs[i]<<" -> 2*8\n";
      m_sf*=2*8;
    }
    else if (flavs[i].IsQuark()) {
      msg_Debugging()<<"     "<<std::setw(15)<<flavs[i]<<" -> 2*3\n";
      m_sf*=2*3;
    }
  }
  msg_Debugging()<<"} -> "<<m_sf<<" / "<<m_norm<<"\n";
  Construct();
  m_cress.resize(m_ampl.Results().size(),
		 std::vector<Complex>(m_cords.size()));
}

N_Parton_CDBG_T::~N_Parton_CDBG_T()
{
}

double N_Parton_CDBG_T::Factorial(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

class Order_KF {
public:
  bool operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
  { return a.Kfcode()<b.Kfcode(); }
};// end of class Order_KF

void N_Parton_CDBG_T::OrderFlavours()
{
  m_flavs[0]=m_flavs[0].Bar();
  m_flavs[1]=m_flavs[1].Bar();
  Flavour_Vector ofl(m_flavs);
  msg_Debugging()<<METHOD<<"(): Old: ";
  for (size_t i(0);i<m_flavs.size();++i) msg_Debugging()<<" "<<m_flavs[i];
  msg_Debugging()<<"\n";
  std::sort(m_flavs.begin(),m_flavs.end(),Order_KF());
  bool sl(false);
  for (size_t i(1);i<m_flavs.size();++i)
    if (m_flavs[i]==m_flavs.front().Bar()) {
      for (size_t i(1);i<m_flavs.size()-1;++i) m_flavs[i]=m_flavs[i+1];
      m_flavs.back()=m_flavs.front().Bar();
      if (m_flavs.front().IsAnti())
	std::swap<Flavour>(m_flavs.front(),m_flavs.back());
      sl=true;
      break;
    }
  for (size_t i(1);i<m_flavs.size()-(sl?2:1);++i)
    if (m_flavs[i].IsQuark()) {
      for (size_t j(i+1);j<m_flavs.size()-(sl?1:0);++j)
	if (m_flavs[i]==m_flavs[j].Bar()) {
	  std::swap<Flavour>(m_flavs[j],m_flavs[i+1]);
	  if (m_flavs[i+1].IsAnti())
	    std::swap<Flavour>(m_flavs[i],m_flavs[i+1]);
	  ++i;
	  break;
	}
    }
  msg_Debugging()<<METHOD<<"(): New: ";
  for (size_t i(0);i<m_flavs.size();++i) msg_Debugging()<<" "<<m_flavs[i];
  msg_Debugging()<<"\n";
  Int_Vector fd(ofl.size(),0);
  m_imap.resize(m_flavs.size(),-1);
  for (size_t i(0);i<m_flavs.size();++i)
    for (size_t j(0);j<ofl.size();++j)
      if (m_flavs[i]==ofl[j]&&fd[j]==0) {
	m_imap[j]=i;
	fd[j]=1;
	break;
      }
  msg_Debugging()<<METHOD<<"(): Map: "<<m_imap<<"\n";
}

void N_Parton_CDBG_T::Construct(Int_Vector ic,Int_Vector left)
{
  if (left.empty()) {
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): Add permutation "<<ic<<".\n";
#endif
    size_t id(m_omap.size());
    m_omap[ic]=id;
    return;
  }
  for (size_t i(0);i<left.size();++i) {
    std::vector<int> nl(left.size()-1);
    int k(0);
    for (size_t j(0);j<left.size();++j,++k) {
      if (j!=i) nl[k]=left[j];
      else --k;
    }
    ic.push_back(left[i]);
    Construct(ic,nl);
    ic.pop_back();
  }  
}

void N_Parton_CDBG_T::ReOrder(Int_Vector &ci,Int_Vector &cj)
{
  Int_Vector ncj(cj.size(),-1);
  for (size_t i(0);i<ci.size();++i)
    if (ci[i]!=(int)i) {
      for (size_t j(0);j<cj.size();++j) 
	if (cj[j]==ci[i]) {
	  ncj[j]=ci[i]=i;
	  break;
	}
    }
  for (size_t j(0);j<ncj.size();++j)
    if (ncj[j]<0) ncj[j]=cj[j];
  cj=ncj;
}

Complex N_Parton_CDBG_T::ColorFactor(const std::vector<int> &cor,
				     const std::vector<int> &coa)
{
  Expression ex;
  ex.back() = CNumber::New(Complex(1.0,0.0));
  if (!AddColors(&ex,false,cor)) return Complex(0.0,0.0);
  if (!AddColors(&ex,true,coa)) return Complex(0.0,0.0);
  ex.Evaluate();
  msg_Debugging()<<METHOD<<"(): "<<cor<<"(*)"<<coa
		 <<" -> "<<ex.Result()<<"\n";
  return ex.Result();
}

bool N_Parton_CDBG_T::AddColors(Expression *const ex,const bool &b,
				const std::vector<int> &co)
{
  size_t os(ex->size());
  Flavour pf(m_flavs[co.front()]);
  bool pg(pf.IsGluon());
  if (!pg&&!m_flavs[co.back()].IsAnti()) return false;
  for (size_t ai(b?201:101), bi(b?401:301), i(pg?0:1);
       i<m_flavs.size()-(pg?0:1);++i) {
    size_t id(co[i]);
    if (m_flavs[id].IsGluon()) {
      if (pf.IsAnti()) return false;
      ex->push_back(Fundamental::New(id,b?ai+1:ai,b?ai:ai+1));
      ++ai;
    }
    else if (m_flavs[id].IsAnti()) {
      if (pf.IsAnti()) return false;
      ex->push_back(Fundamental::New(bi,b?ai+1:ai,b?ai:ai+1));
      ex->push_back(Fundamental::New(bi,id,id));
      ++ai;
      ++bi;
    }
    else {
      if (!pf.IsAnti()) return false;
      if (b) ((Fundamental*)ex->back())->SetJ(id);
      else ((Fundamental*)ex->back())->SetI(id);
    }
    pf=m_flavs[id];
  }
  if (pg) {
    if (b) ((Fundamental*)ex->back())->SetI(201);
    else ((Fundamental*)ex->back())->SetJ(101);
  }
  else {
    if (b) {
      ((Fundamental*)(*ex)[os])->SetJ(co.front());
      if (pf.IsQuark()) ((Fundamental*)*----ex->end())->SetI(co.back());
      else ((Fundamental*)ex->back())->SetI(co.back());
    }
    else {
      ((Fundamental*)(*ex)[os])->SetI(co.front());
      if (pf.IsQuark()) ((Fundamental*)*----ex->end())->SetJ(co.back());
      else ((Fundamental*)ex->back())->SetJ(co.back());
    }
  }
  return true;
}

bool N_Parton_CDBG_T::CheckOrder(const Int_Vector &co)
{
  Flavour pf(m_flavs[co.back()]);
  for (size_t i(0);i<m_flavs.size()-1;++i) {
    if (m_flavs[co[i]].IsQuark()&&
	m_flavs[co[i]].IsAnti()==pf.IsAnti()) return false;
    pf=m_flavs[co[i]];
  }
  msg_Debugging()<<METHOD<<"(): Found "<<co<<"\n";
  return true;
}

bool N_Parton_CDBG_T::CheckOrder(const Int_Vector &cor,
				 const Int_Vector &coa)
{
  for (size_t i(0);i<m_flavs.size();++i)
    if (m_flavs[cor[i]]!=m_flavs[coa[i]]) return false;
  return true;
}

void N_Parton_CDBG_T::Construct()
{
  m_omap.clear();
  m_colfs.clear();
  m_evals.clear();
  Int_Vector ic(1,0), left(m_nin+m_nout-1);
  for (size_t i(1);i<m_nin+m_nout;++i) left[i-1]=i;
  Construct(ic,left);
  std::map<std::pair<Int_Vector,Int_Vector>,Complex> ress;
  std::map<Int_Vector,size_t>::iterator mit(m_omap.begin());
  for (;mit!=m_omap.end();++mit) {
    m_cords.push_back(mit->first);
    m_colfs.push_back(std::vector<Complex>(m_omap.size(),Complex(0.0,0.0)));
    m_evals.push_back(std::vector<int>(m_omap.size(),0));
    if (!CheckOrder(m_cords.back())) {
      m_cords.back().clear();
      continue;
    }
    for (std::map<Int_Vector,size_t>::iterator nit(m_omap.begin());
	 nit!=m_omap.end();++nit) {
      Int_Vector ci(mit->first), cj(nit->first);
      if (!CheckOrder(ci,cj)) continue;
      if (m_ng>(int)(m_nin+m_nout-4)) ReOrder(ci,cj);
      std::pair<Int_Vector,Int_Vector> fid(ci,cj);
      std::map<std::pair<Int_Vector,Int_Vector>,Complex>::const_iterator
	fit(ress.find(fid));
      if (fit!=ress.end()) {
	m_colfs[mit->second][nit->second]=fit->second;
      }
      else {
 	m_colfs[mit->second][nit->second]=ColorFactor(ci,cj);
	ress[fid]=m_colfs[mit->second][nit->second];
      }
    }
  }
  for (size_t j(0);j<m_colfs.size();++j) {
    m_evals[j][j]=1;
    for (size_t k(j+1);k<m_colfs[j].size();++k)
      if (m_colfs[j][k]+m_colfs[k][j]!=Complex(0.0,0.0)) 
	m_evals[k][j]=m_evals[j][k]=1;
  }
}

double N_Parton_CDBG_T::Differential(const std::vector<Vec4D> &momenta)
{
  m_moms[m_imap[0]]=-1.0*momenta[0];
  m_moms[m_imap[1]]=-1.0*momenta[1];
  for (size_t j(2);j<m_moms.size();++j) m_moms[m_imap[j]]=momenta[j];
  m_ampl.SetMomenta(m_moms);
  for (size_t i(0);i<m_cords.size();++i) {
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): Calculate permutation "<<i
		   <<" -> "<<m_cords[i]<<".\n";
#endif
    if (m_cords[i].empty()) continue;
    Int_Vector ic(m_flavs.size(),0), jc(m_flavs.size(),0);
    for (size_t lc(0), k(0);k<m_flavs.size();++k) {
      size_t id(m_cords[i][k]);
      if (m_flavs[id].IsGluon()) jc[id]=(ic[id]=++lc)-1;
      else if (m_flavs[id].IsAnti()) jc[id]=lc; 
      else ic[id]=++lc;
    }
    if (m_flavs.front().IsGluon()) jc.front()=m_flavs.size();
    m_ampl.SetColors(ic,jc);
    m_ampl.EvaluateAll();
    for (size_t j(0);j<m_ampl.Results().size();++j) {
      m_cress[j][i]=m_ampl.Results()[j];
#ifdef DEBUG__BG
      if (msg.LevelIsDebugging()) {
	msg.Out()<<"A["<<j<<"]("<<m_cords[i].front()
		 <<(m_ampl.Chiralities()[j][m_cords[i].front()]>0?'+':'-');
	for (size_t k(1);k<m_cords[i].size();++k)
	  msg.Out()<<","<<m_cords[i][k]
		   <<(m_ampl.Chiralities()[j][m_cords[i][k]]>0?'+':'-');
	msg.Out()<<") -> "<<m_cress[j][i]<<" "<<std::abs(m_cress[j][i])<<"\n";
      }
#endif
    }
  }
  double csum(0.0);
  for (size_t i(0);i<m_cress.size();++i) {
    Complex ccs(0.0,0.0);
    for (size_t j(0);j<m_cress[i].size();++j) {
      if (m_cords[j].empty()) continue;
      ccs+=m_cress[i][j]*std::conj(m_cress[i][j])*m_colfs[j][j];
      for (size_t k(j+1);k<m_cress[i].size();++k) {
	if (m_evals[j][k])
	  ccs+=(m_cress[i][j]*std::conj(m_cress[i][k])).real()*
	    (m_colfs[j][k]+m_colfs[k][j]);
      }
    }
    csum+=ccs.real();
  }
  csum*=m_norm/m_sf;
  return csum;
}

bool N_Parton_CDBG_T::GaugeTest(std::vector<Vec4D> momenta)
{
  m_moms[m_imap[0]]=-1.0*momenta[0];
  m_moms[m_imap[1]]=-1.0*momenta[1];
  for (size_t j(2);j<m_moms.size();++j) m_moms[m_imap[j]]=momenta[j];
  if (!m_ampl.GaugeTest(m_moms)) return false;
  return true;
}

