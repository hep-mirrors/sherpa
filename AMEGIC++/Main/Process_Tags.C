#include "AMEGIC++/Main/Process_Tags.H"
#include "ATOOLS/Math/MathTools.H"
#include "AMEGIC++/Main/Point.H"
#include <algorithm>
#include <iostream>

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Check_External_Flavours AMEGIC::CF;

Process_Tags::Process_Tags(ATOOLS::Flavour *fl,Pol_Info *pl)
{
  if (fl) p_fl = new Flavour(*fl);
  else p_fl=NULL;
  if (pl) p_pl = new Pol_Info(*pl);
  else p_pl=NULL;
  vector<Process_Tags*> dummy;
  m_sublist.push_back(dummy);
  m_maxqcdjets = 0;
  m_zerowidth  = 0;
}

Process_Tags::Process_Tags(Process_Tags *pi)
{
  if (pi->Flav()) p_fl = new Flavour(*(pi->Flav()));
  else p_fl=NULL;
  if (pi->Pol()) p_pl = new Pol_Info(*(pi->Pol()));
  else p_pl=NULL;
  vector<Process_Tags*> dummy;
  m_sublist.push_back(dummy);
  m_sublist[0].clear();  
  for (int i=0;i<pi->Nout();i++) m_sublist[0].push_back(new Process_Tags(pi->m_sublist[0][i]));
  m_maxqcdjets = pi->m_maxqcdjets;
  m_zerowidth = pi->m_zerowidth;
}

Process_Tags::~Process_Tags()
{
  for (size_t j=1;j<m_sublist.size();j++) {
    for (size_t i=0;i<m_sublist[j].size();i++) 
      if(m_sublist[0][i]->p_fl->Size()>1) delete m_sublist[j][i]; 
    m_sublist[j].clear();
  }
  for (size_t i=0;i<m_sublist[0].size();i++) 
    if(m_sublist[0][i]) delete m_sublist[0][i];
  m_sublist[0].clear();
  m_sublist.clear();
  if (p_fl) delete p_fl;
  if (p_pl) delete p_pl;
}

void Process_Tags::AddSubList(int n,ATOOLS::Flavour* fl,Pol_Info* pl)
{
  if (m_sublist[0].size()>0) m_sublist[0].clear();
  for (int i=0;i<n;i++) m_sublist[0].push_back(new Process_Tags(&fl[i],&pl[i]));
}

void Process_Tags::ResetSubList(int n,ATOOLS::Flavour* fl,Pol_Info* pl)
{
  if (m_sublist[0].size()!= (size_t) n) {
    std::cout<<" Process_Tags::ResetSubList : wrong particle number: "<<n<<" vs. "<<m_sublist[0].size()<<std::endl;
    abort();
  }
  for (int i=0;i<n;i++) {
    *(m_sublist[0][i]->p_fl)=fl[i];;
    *(m_sublist[0][i]->p_pl)=pl[i];
  }
}

int Process_Tags::Nout()
{
  return m_sublist[0].size();
}

std::string Process_Tags::PNID() const
{
  size_t pn(0);
  std::string id;
  for (size_t i=0;i<m_sublist[0].size();i++) 
    if (m_sublist[0][i]->Nout()>0) {
      if (pn>0) id+=ToString(pn);
      id+="["+m_sublist[0][i]->PNID()+"]";
      pn=0;
    }
    else ++pn;
  return id+ToString(pn);
}

int Process_Tags::TotalNout()
{
  if (m_sublist[0].size()==0) return 1;
  int n=0;
  for (size_t i=0;i<m_sublist[0].size();i++) n+=m_sublist[0][i]->TotalNout();
  return n;
}

int Process_Tags::Nmax(int nin) 
{
  int k=m_sublist[0].size();
  for (size_t i=0;i<m_sublist[0].size();i++) k=Max(m_sublist[0][i]->Nmax(1),k);
  return k+nin;
}

int Process_Tags::Ndecays()
{
  int k=0;
  if (p_fl && m_sublist[0].size()>0) k=1;
  for (size_t i=0;i<m_sublist[0].size();i++) k+= m_sublist[0][i]->Ndecays();
  return k;
}

int Process_Tags::GetDPOffset(int n)
{
  int m=0;
  return GetDPOffset(n,m);
}

int Process_Tags::GetDPOffset(int &n,int &m)
{
  if (n<0) return 0;
  if (m_sublist[0].size()==0) {
    n++;m++;
    return 0;
  }
  if (n==0) return m+TotalNout()-Nout();
  for (size_t i=0;i<m_sublist[0].size();i++) {
    int a = m_sublist[0][i]->GetDPOffset(--n,m);
    if (a>0) return a;
  }
  return 0;
}

Process_Tags* Process_Tags::GetDecay(int &n)
{
  if (n<0 || m_sublist[0].size()==0) {
    n++;
    return NULL;
  }
  if (n==0) return this;
  Process_Tags* pi=NULL;
  for (size_t i=0;i<m_sublist[0].size();i++) {
    pi = m_sublist[0][i]->GetDecay(--n);
    if (pi) return pi;
  }
  return NULL;
}

void Process_Tags::GetFlavList(ATOOLS::Flavour* fl, int n)
{
  for (size_t i=0;i<m_sublist[n].size();i++) fl[i]=*(m_sublist[n][i]->p_fl);
}

size_t Process_Tags::GetStableFlavList(ATOOLS::Flavour* fl, int n)
{
  size_t cnt=0;
  for (size_t i=0;i<m_sublist[n].size();i++) 
    if (m_sublist[n][i]->Nout()==0) {  
      fl[cnt]=*(m_sublist[n][i]->p_fl);
      cnt++;
    }
  return cnt;
}

void Process_Tags::GetPolList(Pol_Info* pl)
{
  for (size_t i=0;i<m_sublist[0].size();i++) pl[i]=*(m_sublist[0][i]->p_pl);
}

Process_Tags* Process_Tags::GetSubProcess(int n)
{
  int dn=1;
  return GetSubProcess(n,dn);
}

Process_Tags* Process_Tags::GetSubProcess(int n, int& dn)
{
  Process_Tags* pi = new Process_Tags(p_fl,p_pl);
  pi->m_maxqcdjets = m_maxqcdjets;
  pi->m_zerowidth  = m_zerowidth;
  if (m_sublist[0].size()==0) return pi;

  int cn=n/dn;
  for (size_t j=0;j<m_sublist[0].size();j++) {
    int k=0;
    if (m_sublist.size()>1) k=1;
    if (n>=0) pi->m_sublist[0].push_back(m_sublist[cn%(m_sublist.size()-k)+k][j]->GetSubProcess(n,dn));
    else pi->m_sublist[0].push_back(m_sublist[0][j]->GetSubProcess(n,dn));
  }
  if (m_sublist.size()>2)dn*=m_sublist.size()-1;

  return pi;
}


int Process_Tags::GetTotalFlavList(ATOOLS::Flavour* fl, int n)
{
//    cout<<"TF: "<<n<<" "<<m_sublist.size()-1<<endl;
  if (m_sublist[0].size()==0) {
    fl[0]=*p_fl;
    return 1;
  }
  size_t i=0, dn=1;
  if (m_sublist.size()>2)dn*=m_sublist.size()-1;
//   cout<<"y: "<<m_sublist.size()<<" "<<dn<<endl;
  for (size_t j=0;j<m_sublist[0].size();j++) {
    int k=0;
    if (m_sublist.size()>1) k=1;
 //    cout<<"x: ";
//     if (p_fl) cout<<*p_fl;
//     cout<<" "<<n%(m_sublist.size()-k)+k<<" "<<dn<<" "<<n/dn<<endl;
    if (n>=0) dn*=m_sublist[n%(m_sublist.size()-k)+k][j]->GetTotalFlavList(&fl[i],n/dn);
    else m_sublist[0][j]->GetTotalFlavList(&fl[i]);
    i+=m_sublist[0][j]->TotalNout();
  }
//   cout<<"y: "<<m_sublist.size()<<endl;
//   if (m_sublist.size()>2)dn*=m_sublist.size()-1;
  return dn;
}

size_t Process_Tags::GetOnshellFlavList(ATOOLS::Flavour_Vector &fl, vector<Process_Tags*> &decaylist, bool first)
{
  if (m_sublist[0].size()==0) {
    fl.push_back(*p_fl);
    decaylist.push_back(NULL);
    return 1;
  }
  if (!first&&m_zerowidth==1) {
    fl.push_back(*p_fl);
    decaylist.push_back(this);
    return 1;
  }
  size_t i=0;
  for (size_t j=0;j<m_sublist[0].size();j++) {
    i+=m_sublist[0][j]->GetOnshellFlavList(fl,decaylist,0);
  }
  return i;
}

int Process_Tags::OSDecays()
{
  int k=m_zerowidth;
  for (size_t i=0;i<m_sublist[0].size();i++) k+= m_sublist[0][i]->OSDecays();
  return k;
}

void Process_Tags::GetTotalPolList(Pol_Info* pl)
{
  if (m_sublist[0].size()==0) {
    pl[0]=*p_pl;
    return;
  }
  size_t i=0;
  for (size_t j=0;j<m_sublist[0].size();j++) {
    m_sublist[0][j]->GetTotalPolList(&pl[i]);
    i+=m_sublist[0][j]->TotalNout();
  }
}

Process_Tags* Process_Tags::FindDM(std::string c)
{
  if (p_pl) {
    if (p_pl->pol_type=='d' && p_pl->type[0]==c[0]) return this;
  }
  for (size_t i=0;i<m_sublist[0].size();i++) {
    Process_Tags* pi=m_sublist[0][i]->FindDM(c);
    if (pi) return pi;
  }
  return NULL;
}

bool Process_Tags::CheckCompleteness()
{
  if (m_sublist[0].size()==0) {
    if (p_pl && p_pl->pol_type!='d') return true;
    return false;
  }
  for (size_t i=0;i<m_sublist[0].size();i++) {
    if (!m_sublist[0][i]->CheckCompleteness()) return false;
  }  
  return true;
}

void Process_Tags::GetOSConditions(vector<pair<string, double> >& osc,int &cnt)
{
  if (m_sublist[0].size()==0) {
    cnt++;
    return;
  }

  for (size_t i=0;i<m_sublist[0].size();i++) {
    if (m_sublist[0][i]->m_zerowidth) {
      int nt=m_sublist[0][i]->TotalNout();
      string str="";
      for (int j=cnt;j<nt+cnt;j++) str+=ToString(j);
      osc.push_back(pair<string, double>(str,sqr(m_sublist[0][i]->Flav()->Mass())));
    }
    m_sublist[0][i]->GetOSConditions(osc,cnt);
  }
  return;
}

void Process_Tags::Print()
{
  if (p_fl==0) cout<<" Final State:";
  for (size_t i=0;i<m_sublist[0].size();i++) {
    cout<<" "<<*(m_sublist[0][i]->p_fl);
    if (m_sublist[0][i]->m_sublist[0].size()>0) {
      if (m_sublist[0][i]->m_zerowidth) cout<<"|";
      cout<<"(->";
      m_sublist[0][i]->Print();
	cout<<")";
    }
  }
  if (p_fl==0) cout<<endl;
}

void Process_Tags::FullPrint()
{
  if (p_fl==0) cout<<" Final State:";
  if (m_sublist.size()==1) {
    for (size_t i=0;i<m_sublist[0].size();i++) {
      cout<<" "<<*(m_sublist[0][i]->p_fl);
      if (m_sublist[0][i]->m_sublist[0].size()>0) {
	if (m_sublist[0][i]->m_zerowidth) cout<<"|";
	cout<<"(->";
	m_sublist[0][i]->FullPrint();
	cout<<")";
      }
    }
  }
  else {
    for (size_t j=1;j<m_sublist.size();j++) {
      for (size_t i=0;i<m_sublist[j].size();i++) {
	cout<<" "<<*(m_sublist[j][i]->p_fl);
	if (m_sublist[j][i]->m_sublist[0].size()>0) {
	  if (m_sublist[j][i]->m_zerowidth) cout<<"|";
	  cout<<"(->";
	  m_sublist[j][i]->FullPrint();
	  cout<<")";
	}
      }
      if (j<m_sublist.size()-1) cout<<" |";
    }
  }
  if (p_fl==0) cout<<endl;
}

void Process_Tags::Expand()
{
  if (p_fl==0) {
    for (size_t i=0;i<m_sublist[0].size();i++) m_sublist[0][i]->Expand();
    return;
  }
  
  if (m_sublist.size()>1) return;
  kf_code container=0;
  int mode=0; 
  for (size_t i=0;i<m_sublist[0].size();i++) {
    if(m_sublist[0][i]->p_fl->Size()>1) {
      if (mode>0 && container!=m_sublist[0][i]->p_fl->Kfcode()) mode=-10;
      else mode++;
      container = m_sublist[0][i]->p_fl->Kfcode();
    }
  }
  if (container>0) {
    size_t nout=m_sublist[0].size();
    Flavour* flout = new Flavour[nout];
    Pol_Info* plout = new Pol_Info[nout];
    GetFlavList(flout);
    GetPolList(plout);
    size_t  * flindex = new size_t[nout];
    for (size_t i=0;i<nout;i++) flindex[i] = 0;
    bool flag = 1;
    for (;;) {
      if (!flag) break;
      for (size_t i=0;i<nout;++i) {
	if (m_sublist[0][i]->p_fl->Size() != 1) {
	  flout[i] = (*(m_sublist[0][i]->p_fl))[flindex[i]]; 
	}
      }
      if (CF.ValidProcess(1,p_fl,nout,flout)) {
	vector<Process_Tags*> dummy;
	for (size_t i=0;i<nout;i++) {
	  if (m_sublist[0][i]->p_fl->Size()==1) dummy.push_back(m_sublist[0][i]);
	  else {
	    dummy.push_back(new Process_Tags(&flout[i],&plout[i]));
	    dummy[dummy.size()-1]->m_zerowidth=m_sublist[0][i]->m_zerowidth;
	  }
	}
	m_sublist.push_back(dummy);
      }
      for (int i=nout-1;i>=0;--i) {
	if (m_sublist[0][i]->p_fl->Size()-1>flindex[i]) {
	  ++flindex[i];
	  break;
	}
	else {
	  if (i==0) flag = 0;
	  else {
	    size_t maxi = 0;
	    for (size_t j=0;j<i;++j) if (flindex[j]>maxi) maxi=flindex[j];
	    if (mode<0) flindex[i] = 0;
	    else flindex[i] = maxi+1;
	  }
	}
	if (i==0) break;
      }
    }
    delete [] plout;
    delete [] flout;
    delete [] flindex;
  }

  for (size_t j=0;j<m_sublist[0].size();j++) m_sublist[0][j]->Expand();  
}

int Process_Tags::NProcs()
{
  int cnt=1;
  for (size_t i=0;i<m_sublist[0].size();i++) {
    cnt*=m_sublist[0][i]->NProcs();
  }
  if (m_sublist.size()>2) cnt*=m_sublist.size()-1;
  return cnt;
}

void Process_Tags::MergePointList(Point** plist,Point* np, int nin) 
{
  int nd=0;
  int ep=nin;
  MergePointList(plist,np,nd,nin,ep);
}

Point* Process_Tags::MergePointList(Point** plist,Point* np,int &nd, int nin, int &ep)
{
  Point* hp;
  hp = np->CopyList(plist[nd++]);
  
  for (size_t i=0;i<m_sublist[0].size();i++) {
    for (size_t j=0;j<2*(m_sublist[0].size()+nin)-3;j++) {
      if (np[j].b==1 && np[j].number<99) {
	if ((int)i == np[j].number-nin) {
	  if (m_sublist[0][i]->m_sublist[0].size()>0) {
	    Point* hhp = m_sublist[0][i]->MergePointList(plist,hp+1,nd,1,ep);
	    np[j] = hp[1];
	    np[j].number = 100;
	    hp = hhp;
 	    np[j].t = 10+m_sublist[0][i]->m_maxqcdjets;
	    np[j].zwf = m_sublist[0][i]->m_zerowidth;
	  }
	  else {
	    np[j].number = ep;
	    np[j].b = 2;
	    ep++;
	  }
	}
      }
    }
  }
  return hp;
}

bool Check_External_Flavours::ValidProcess(int _nin,Flavour * _in,
					   int _nout,Flavour * _out) {
  // first : sum over all invariants and compare
  for (int i=0;i<_nin;i++)  { if (_in[i].Size()>1)  return 1; }
  for (int i=0;i<_nout;i++) { if (_out[i].Size()>1) return 1; }
  
  int    chin  = 0, chout  = 0;
  int    sin  = 0, sout  = 0;
  double bin  = 0, bout  = 0;
  for (int i=0;i<_nin;i++) {
    chin   += _in[i].IntCharge();
    sin   += _in[i].IntSpin();
    bin   += BaryonNumber(_in[i]);
  }
  for (int i=0;i<_nout;i++) {
    chout  += _out[i].IntCharge();
    sout  += _out[i].IntSpin();
    bout  += BaryonNumber(_out[i]);
  }
  sin = sin%2; sout = sout%2;
  
  if (chin  != chout) return 0;    // electric charge violation
  if (sin  != sout) return 0;    // spin/fermion number violation
  if (bin  != bout) return 0;    // baryon number violation
  return 1;
}

bool Check_External_Flavours::PureGluonic(int _nin,Flavour * _in,
					  int _nout,Flavour * _out) {
  for (int i=0;i<_nin;i++)  { if (!_in[i].IsGluon())  return 0; }
  for (int i=0;i<_nout;i++) { if (!_out[i].IsGluon()) return 0; }
  return 1;
}


///////////////////////////// MHV ///////////////////////////////////
bool Check_External_Flavours::MHVCalculable(const PHASIC::Process_Info& pi) {
  if (pi.m_fi.m_ps.size()!=pi.m_fi.NExternal()) return 0;
  vector<ATOOLS::Flavour> flin;
  vector<ATOOLS::Flavour> flout;
  pi.m_ii.GetExternal(flin);
  pi.m_fi.GetExternal(flout);

    int n_gl(0);
    int n_q(0);
    int n_l(0);
    for (size_t i=0;i<flin.size();i++)  {
      if (flin[i].IsMassive()) return 0;
      if (!flin[i].IsGluon()) {
	if (!flin[i].IsQuark()) {
	  if (!flin[i].IsLepton()) return 0;
	  else n_l++;
	}
	else n_q++;
      }
      else n_gl++;
    }
    for (size_t i=0;i<flout.size();i++) { 
      if (flout[i].IsMassive()) return 0;
      if (!flout[i].IsGluon()) {
	if (!flout[i].IsQuark()) {
	  if (!flout[i].IsLepton()) return 0;
	  else n_l++;
	}
	else n_q++;
      }
      else n_gl++;
    }
    if (flin.size()+flout.size()>9 || n_l>0) return 0;
    if (n_q<=2) return 1;
    // !!! TEMPORARY BUGFIX FOR BLACKHAT !!!
    // if (n_q==4 && (pi.m_oew==0||pi.m_maxoew==0) ) return 1;
    return 0;
//     if (n_q>4 || flin.size()+flout.size()>9 || n_l>2 || (n_l>0 && n_q!=2)) return 0;
    return 1;
}
/////////////////////////////////////////////////////////////////////


class Order_FVST {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_FVST::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a < b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->IsScalar() && !b->m_sublist[0][i]->p_fl->IsScalar()) return 0;
      if (a->m_sublist[0][i]->p_fl->IsVector() && !b->m_sublist[0][i]->p_fl->IsScalar() && !b->m_sublist[0][i]->p_fl->IsVector()) return 0;
      if (a->m_sublist[0][i]->p_fl->IsFermion() && !b->m_sublist[0][i]->p_fl->IsFermion() && 
	  !b->m_sublist[0][i]->p_fl->IsScalar() && !b->m_sublist[0][i]->p_fl->IsVector()) return 0;      
   }
    return 0;
  }
  if (a->p_fl->IsFermion() && !b->p_fl->IsFermion()) return 1;
  if (a->p_fl->IsVector() && !b->p_fl->IsFermion() && !b->p_fl->IsVector()) return 1;
  if (a->p_fl->IsScalar() && !b->p_fl->IsScalar() && 
      !b->p_fl->IsFermion() && !b->p_fl->IsVector()) return 1;
  return 0;
}

class Order_SVFT {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_SVFT::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a < b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->IsFermion() && !b->m_sublist[0][i]->p_fl->IsFermion()) return 0;
      if (a->m_sublist[0][i]->p_fl->IsVector() && !b->m_sublist[0][i]->p_fl->IsFermion() && !b->m_sublist[0][i]->p_fl->IsVector()) return 0;
      if (a->m_sublist[0][i]->p_fl->IsScalar() && !b->m_sublist[0][i]->p_fl->IsScalar() && 
	  !b->m_sublist[0][i]->p_fl->IsFermion() && !b->m_sublist[0][i]->p_fl->IsVector()) return 0;
    }
    return 0;
  }
  if (a->p_fl->IsScalar() && !b->p_fl->IsScalar()) return 1;
  if (a->p_fl->IsVector() && !b->p_fl->IsScalar() && !b->p_fl->IsVector()) return 1;
  if (a->p_fl->IsFermion() && !b->p_fl->IsFermion() && 
      !b->p_fl->IsScalar() && !b->p_fl->IsVector()) return 1;
  return 0;
}

class Order_Mass {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_Mass::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a > b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->Mass() < b->m_sublist[0][i]->p_fl->Mass()) return 0;
    }
    return 0;
  }
  if (a->p_fl->Mass() <= b->p_fl->Mass()) return 0;
  return 1;
}

class Order_InvMass {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_InvMass::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a < b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->Mass() > b->m_sublist[0][i]->p_fl->Mass()) return 0;
    }
    return 0;
  }
  if (a->p_fl->Mass() < b->p_fl->Mass()) return 1;
  return 0;
}


class Order_Kfc {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_Kfc::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a < b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->Kfcode() > b->m_sublist[0][i]->p_fl->Kfcode()) return 0;
    }
    return 0;
  }
  if (a->p_fl->Kfcode() < b->p_fl->Kfcode()) return 1;
  return 0;
}


class Order_Anti {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_Anti::operator()(const Process_Tags * a, const Process_Tags * b) {
  //    if "a < b" return 1  else 0;
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (!(a->m_sublist[0][i]->p_fl->IsFermion() && b->m_sublist[0][i]->p_fl->IsFermion())) return 0;
      if ((a->m_sublist[0][i]->p_fl->IsAnti() && !b->m_sublist[0][i]->p_fl->IsAnti())) return 0;
    }
    return 0;
  }
  if ((a->p_fl->IsFermion() && b->p_fl->IsFermion())
      && (!a->p_fl->IsAnti() && b->p_fl->IsAnti())) return 1;
  return 0;
}


class Order_Coupling {
public:
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_Coupling::operator()(const Process_Tags * a, const Process_Tags * b) {
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) {
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
      if (a->m_sublist[0][i]->p_fl->Strong() && !b->m_sublist[0][i]->p_fl->Strong()) return 0;
    }
    return 0;
  }
  if (!a->p_fl->Strong() && b->p_fl->Strong()) return 1;
  return 0;
}

static FMMap s_fmm;

class Order_FlavMulti {
//   FMMap* p_fmm;
public:
//   Order_FlavMulti(FMMap* fmm) {p_fmm=fmm;}
  int operator()(const Process_Tags * a, const Process_Tags * b);
};
int Order_FlavMulti::operator()(const Process_Tags * a, const Process_Tags * b) {
  if (*a->p_fl==*b->p_fl && (a->m_sublist[0].size()>0 || b->m_sublist[0].size()>0)) {
    if (a->m_sublist[0].size()>b->m_sublist[0].size()) return 1;
    if (a->m_sublist[0].size()<b->m_sublist[0].size()) return 0;
    for (size_t i=0;i<a->m_sublist[0].size();++i) 
      if (operator()(a->m_sublist[0][i],b->m_sublist[0][i])) return 1;
    return 0;
  }
  if (s_fmm[int(a->p_fl->Kfcode())]==0 || s_fmm[int(b->p_fl->Kfcode())]==0) return 0;
  if (s_fmm[int(a->p_fl->Kfcode())]>s_fmm[int(b->p_fl->Kfcode())]) return 1;
  return 0;
}


//
// Note: all order operator have to return 0 if 
//       two elements are equal!
//       Otherwise the order will change even for
//       equal elements.
//


void Process_Tags::Reshuffle(Process_Tags *cpi)
{
  int n=m_sublist[0].size();
  if (n==0) return;
  for (int i=0;i<n;++i) m_sublist[0][i]->Reshuffle();
  s_fmm.clear();

  int nin = 0;
  if (cpi) nin = cpi->Nout();
  if (nin>0) {
    for (int i=0;i<nin;i++) {
      Flavour *hfl=cpi->m_sublist[0][i]->Flav();
      if (s_fmm.find(int(hfl->Kfcode()))==s_fmm.end()) 
	s_fmm[int(hfl->Kfcode())]=0;
      if (hfl->IsFermion()) s_fmm[int(hfl->Kfcode())]++;
//         cout<<i<<"I: "<<*hfl<<" "<<int(hfl->Kfcode())<<" "<<s_fmm[int(hfl->Kfcode())]<<endl;
    }
  }

  Flavour *flav = new Flavour[n];
  GetFlavList(flav);
  Flavour heaviest(kf_photon);
  for (int i=0;i<n;++i) {
    if (flav[i].Mass()>heaviest.Mass()) heaviest=flav[i];
    else if (flav[i].Mass()==heaviest.Mass() &&
	     !flav[i].IsAnti()) heaviest=flav[i];
  }

  for (int i=0;i<n;++i) {
    if (s_fmm.find(int(flav[i].Kfcode()))==s_fmm.end()) 
      s_fmm[int(flav[i].Kfcode())]=0;
    if (flav[i].IsFermion() && !(flav[i].IsMassive())) s_fmm[int(flav[i].Kfcode())]++;
  } 

  delete [] flav;

  std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_Kfc());
  std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_Anti());
  std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_SVFT());
  if (s_fmm.size()>0) std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_FlavMulti());
  if (heaviest.IsAnti())  std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_InvMass());
  else   std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_Mass());
  std::stable_sort(m_sublist[0].begin(),m_sublist[0].end(),Order_Coupling());
}

string Process_Tags::GenerateName()
{
  string name;

  for (size_t i=0;i<m_sublist[0].size();i++) {
    name+=m_sublist[0][i]->p_fl->IDName();
    if (m_sublist[0][i]->p_pl->pol_type=='c' && m_sublist[0][i]->p_pl->num==1) {
      if (m_sublist[0][i]->p_pl->type[0]==-1) name += string("m");
      if (m_sublist[0][i]->p_pl->type[0]==+1) name += string("p");
      if (m_sublist[0][i]->p_pl->type[0]==0)  name += string("z");
    }
    else if (m_sublist[0][i]->p_pl->pol_type=='h' && m_sublist[0][i]->p_pl->num==1) {
      if (m_sublist[0][i]->p_pl->type[0]==-1) name += string("m");
      if (m_sublist[0][i]->p_pl->type[0]==+1) name += string("p");
    }
    if (m_sublist[0][i]->m_sublist[0].size()>0) {
      name+=string("[");
      name+=m_sublist[0][i]->GenerateName();
      name+=string("]");
    }
    if (i<m_sublist[0].size()-1) name+=string("__");
  }

  return name;
}
