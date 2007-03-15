#include "N_Gluon_BG.H"

#include "Message.H"
#include "Random.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include "Color.H"
#include "STL_Tools.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

N_Gluon_BG::N_Gluon_BG(const size_t nin,const size_t nout,
			   const std::vector<Flavour> &flavs):
  m_nin(nin), m_nout(nout), m_mode(0), m_tests(0)
{ 
  std::vector<int> types(m_nin+m_nout,0);
  m_ampl.Construct(flavs,types);
  Int_Vector ic(m_nin+m_nout), jc(m_nin+m_nout);
  for (size_t j(1);j<m_nin+m_nout;++j) ic[j]=jc[j-1]=j;
  ic.front()=jc.back()=0;
  m_moms.resize(m_nin+m_nout);
  std::vector<size_t> ids(m_nin+m_nout,0);
  for (size_t i(0);i<ids.size();++i) {
    ids[i]=i;
    if (flavs[i].IsGluon()) types[i]=0;
    else if (flavs[i].IsAnti()) types[i]=-1;
    else types[i]=1;
  }
  m_as=(*MODEL::as)(sqr(rpa.gen.Ecms()));
  m_sf=Factorial(m_nout);
  Construct();
  m_cress.resize(m_ampl.Results().size(),
		 std::vector<Complex>(m_cords.size()));
}

N_Gluon_BG::~N_Gluon_BG()
{
}

double N_Gluon_BG::Factorial(const double &n) 
{
  if (n<=0.0) return 1.0;
  return n*Factorial(n-1.0);
}

bool N_Gluon_BG::GenerateConfiguration()
{
  return true;
}

void N_Gluon_BG::Construct(Int_Vector ic,Int_Vector left)
{
  if (left.empty()) {
    size_t n(ic.front());
    for (size_t i(1);i<ic.size();++i) ic[i-1]=ic[i];
    ic.back()=n;
    msg_Debugging()<<METHOD<<"(): Add permutation "<<ic<<".\n";
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

void N_Gluon_BG::ReOrder(Int_Vector &ci,Int_Vector &cj)
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

void N_Gluon_BG::Construct()
{
  m_omap.clear();
  m_colfs.clear();
  m_evals.clear();
  Int_Vector ic(2), left(m_nin+m_nout-2);
  ic[0]=m_nin+m_nout-1;
  ic[1]=0;
  for (size_t i(2);i<m_nin+m_nout;++i) left[i-2]=i-1;
  Construct(ic,left);
  std::map<std::pair<Int_Vector,Int_Vector>,Complex> ress;
  std::map<Int_Vector,size_t>::iterator mit(m_omap.begin());
  for (;mit!=m_omap.end();++mit) {
    m_cords.push_back(mit->first);
    m_colfs.push_back(std::vector<Complex>(m_omap.size(),Complex(0.0,0.0)));
    m_evals.push_back(std::vector<int>(m_omap.size(),0));
    for (std::map<Int_Vector,size_t>::iterator nit(m_omap.begin());
	 nit!=m_omap.end();++nit) {
      Int_Vector ci(mit->first), cj(nit->first);
      ReOrder(ci,cj);
      std::pair<Int_Vector,Int_Vector> fid(ci,cj);
      if (ress.find(fid)!=ress.end()) {
	m_colfs[mit->second][nit->second]=ress.find(fid)->second;
      }
      else {
	Expression ex;
	size_t ai(101);
	ex.back() = Adjoint::New(mit->first.front(),
				 *++mit->first.begin(),ai);
	for (size_t i(2);i<mit->first.size()-2;++i,++ai)
	  ex.push_back(Adjoint::New(ai,mit->first[i],ai+1));
	ex.push_back(Adjoint::New(ai,*----mit->first.end(),
				  mit->first.back()));
	size_t bi(201);
	ex.push_back(Adjoint::New(nit->first.front(),
				  *++nit->first.begin(),bi));
	for (size_t i(2);i<nit->first.size()-2;++i,++bi)
	  ex.push_back(Adjoint::New(bi,nit->first[i],bi+1));
	ex.push_back(Adjoint::New(bi,*----nit->first.end(),
				  nit->first.back()));
	ex.Evaluate();
	m_colfs[mit->second][nit->second]=ex.Result();
	ress[fid]=ex.Result();
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

double N_Gluon_BG::Differential(const std::vector<Vec4D> &momenta)
{
  PROFILE_HERE;
  m_moms[0]=-1.0*momenta[0];
  m_moms[1]=-1.0*momenta[1];
  for (size_t j(2);j<m_moms.size();++j) m_moms[j]=momenta[j];
  m_ampl.SetMomenta(m_moms);
  for (size_t i(0);i<m_cords.size();++i) {
    msg_Debugging()<<METHOD<<"(): Calculate permutation "<<i
		   <<" -> "<<m_cords[i]<<".\n";
    m_ampl.EvaluateAll(m_cords[i]);
    for (size_t j(0);j<m_ampl.Results().size();++j) {
      m_cress[j][i]=m_ampl.Results()[j];
      if (msg.LevelIsDebugging()) {
	msg.Out()<<"A["<<j<<"]("<<m_cords[i].front()
		 <<(m_ampl.Chiralities()[j][m_cords[i].front()]>0?'+':'-');
	for (size_t k(1);k<m_cords[i].size();++k)
	  msg.Out()<<","<<m_cords[i][k]
		   <<(m_ampl.Chiralities()[j][m_cords[i][k]]>0?'+':'-');
	msg.Out()<<") -> "<<m_cress[j][i]<<" "<<std::abs(m_cress[j][i])<<"\n";
      }
    }
  }
  double csum(0.0);
  for (size_t i(0);i<m_cress.size();++i) {
    Complex ccs(0.0,0.0);
    for (size_t j(0);j<m_cress[i].size();++j) {
      ccs+=m_cress[i][j]*std::conj(m_cress[i][j])*m_colfs[j][j];
      for (size_t k(j+1);k<m_cress[i].size();++k) {
	if (m_evals[j][k])
	  ccs+=(m_cress[i][j]*std::conj(m_cress[i][k])).real()*
	    (m_colfs[j][k]+m_colfs[k][j]);
      }
    }
    csum+=ccs.real();
  }
  csum*=pow(2.0,m_nout);
  // average over initial chiralities and colors 
  csum/=4.0*64.0*m_sf;
  // add couplings
  csum*=pow(4.0*M_PI*m_as,m_nout);
  return csum;
}

bool N_Gluon_BG::GaugeTest(std::vector<Vec4D> momenta)
{
  PROFILE_HERE;
  momenta[0]=-1.0*momenta[0];
  momenta[1]=-1.0*momenta[1];
  if (!m_ampl.GaugeTest(momenta)) return false;
  return true;
}

