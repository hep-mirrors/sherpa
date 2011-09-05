#include "PHASIC++/Selectors/Cut_Data.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


std::ostream & PHASIC::operator<<(std::ostream & s , Cut_Data & cd)
{
  s<<" Cut Data : "<<cd.ncut<<" \n"<<std::endl;
  for (int i=0;i<cd.ncut;++i) {
    for (int j=0;j<cd.ncut;++j) s<<cd.scut[i][j]<<"  ";
    s<<std::endl;
  }
  return s;
}


Cut_Data::Cut_Data() {
  energymin = 0;
  etmin = 0;
  cosmin    = 0;
  cosmax    = 0;
  scut      = 0;
  ncut      = 0;
}

Cut_Data::~Cut_Data() {
  if (!scut) return;
  for (short int i=0;i<ncut;i++) {
    delete[] cosmin[i];
    delete[] cosmax[i];
    delete[] cosmin_save[i];
    delete[] cosmax_save[i];
    delete[] scut[i];
    delete[] scut_save[i];
  }
  delete[] cosmin;
  delete[] cosmax;
  delete[] cosmin_save;
  delete[] cosmax_save;
  delete[] scut;
  delete[] scut_save;
  delete[] energymin;
  delete[] energymin_save;
  delete[] etmin;
}

void Cut_Data::Init(int _nin,const Flavour_Vector &_fl) {
  if (energymin != 0) return;
  smin = 0.;
  nin            = _nin;
  ncut           = _fl.size();
  fl             = &_fl.front();
  energymin      = new double[ncut];
  energymin_save = new double[ncut];
  etmin          = new double[ncut];
  cosmin         = new double*[ncut];
  cosmax         = new double*[ncut];
  cosmin_save    = new double*[ncut];
  cosmax_save    = new double*[ncut];
  scut           = new double*[ncut];
  scut_save      = new double*[ncut];

  for (int i=0;i<ncut;i++) {
    cosmin[i]      = new double[ncut];
    cosmax[i]      = new double[ncut];
    cosmin_save[i] = new double[ncut];
    cosmax_save[i] = new double[ncut];
    scut[i]        = new double[ncut];
    scut_save[i]   = new double[ncut];
    energymin[i]   = Max(0.,fl[i].SelMass());
    if (fl[i].IsKK()) energymin[i] = 0.;
    smin += energymin_save[i] = energymin[i];
    etmin[i]       = 0.;
  }
  smin = sqr(smin);

  for (int i=0;i<ncut;i++) {
    for (int j=i;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = cosmin_save[i][j] = -1.;
      cosmax[i][j] = cosmax[j][i] = cosmax_save[i][j] = 1.;
      /*      double sc =
	+sqr(fl[i].SelMass())+sqr(fl[j].SelMass())
	+2.*energymin[i]*energymin[j]
	-2.*sqrt(dabs(sqr(energymin[i])-sqr(fl[i].SelMass())))
	*sqrt(dabs(sqr(energymin[j])-sqr(fl[j].SelMass())))
	*cosmax[i][j];*/
      scut[i][j] = scut[j][i] = scut_save[i][j] =
	(i<nin)^(j<nin)?0.0:1.e-12*sqr(rpa->gen.Ecms());
	//Max(sc,1.e-12*sqr(rpa->gen.Ecms()));
    }
  }  
}

void Cut_Data::Complete()
{
  for (int i=0;i<ncut;i++) {
    for (int j=i+1;j<ncut;j++) {
      if ((i<nin)^(j<nin)) continue;
      scut[i][j] =  
// 	Max(scut[i][j],2.*energymin[i]*energymin[j]*(1.-cosmax[i][j])+sqr(fl[i].SelMass())+sqr(fl[j].SelMass()));
	Max(scut[i][j],2.*energymin[i]*energymin[j]-2.*sqrt(sqr(energymin[i])-sqr(fl[i].SelMass()))
	    *sqrt(sqr(energymin[j])-sqr(fl[j].SelMass()))*cosmax[i][j]
	    +sqr(fl[i].SelMass())+sqr(fl[j].SelMass()));
      scut[i][j] = scut[j][i] = 
	Max(scut[i][j],sqr(fl[i].SelMass()+fl[j].SelMass()));
//       std::cout<<i<<","<<j<<": "<<scut[i][j]<<std::endl;
    }
  } 

  MyStrStream strs;
  for (int i=0;i<ncut;i++) {
    energymin_save[i] = energymin[i];
    for (int j=i+1;j<ncut;j++) {
      cosmin_save[i][j] = cosmin[i][j];
      cosmax_save[i][j] = cosmax[i][j];
      scut_save[i][j]   = scut[i][j];
    }
    if (i>=2) strs<<GetIndexID(i);
  }
  std::string str;
  strs>>str;
  smin = 0.;
  double etmm = 0.; 
  double e1=0.,e2=0.;
  for (int i=2;i<ncut;i++) {
    if (etmin[i]>etmm) etmm = etmin[i];
    smin += etmin[i];
    e1 += energymin[i];
    e2 += energymin[i]*cosmax[0][i];
  }
  smin = Max(sqr(smin),sqr(e1)-sqr(e2));
  smin = Max(smin,sqr(2.*etmm));
  smin = Max(Getscut(str),smin);

  msg_Tracking()<<"Cut_Data::Complete(): s_{min} = "<<smin<<endl;
  m_smin_map.clear();
}

void Cut_Data::Reset(bool update)
{
  for (int i=0;i<ncut;i++) {
    energymin[i] = energymin_save[i];
    for (int j=i+1;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = cosmin_save[i][j];
      cosmax[i][j] = cosmax[j][i] = cosmax_save[i][j];
      scut[i][j]   = scut[j][i]   = scut_save[i][j];
    }
  }
  if (update) {
    map<string,double>::iterator it;
    for (it=m_smin_map.begin();it!=m_smin_map.end();++it) it->second = -1.;
  }
}

char Cut_Data::GetIndexID(int id)
{
  char c = id;
  c<10 ? c+=48 : c+=55;
  return c;
}

double Cut_Data::Getscut
(std::vector<char> pl,std::vector<char> pr,int n,int k,int li)
{
  if (n==k) {
    std::string idl, idr;
    for (size_t i(0);i<pl.size();++i) if (pl[i]!=' ') idl+=pl[i];
    for (size_t i(0);i<pr.size();++i) if (pr[i]!=' ') idr+=pr[i];
    double ml(sqrt(Getscut(idl))), mr(sqrt(Getscut(idr)));
#ifdef DEBUG__Cut_Data
    msg_Debugging()<<"m_{"<<idl<<"} + m_{"<<idr<<"} = "
		   <<ml<<" + "<<mr<<" = "<<ml+mr<<"\n";
#endif
    return sqr(ml+mr);
  }
  msg_Indent();
  double sc(0.0);
  for (size_t i(li+1);i<pl.size();++i) {
    std::swap<char>(pl[i],pr[i]);
    sc=Max(sc,Getscut(pl,pr,n,k+1,i));
    std::swap<char>(pl[i],pr[i]);
  }
  return sc;
}

double Cut_Data::Getscut(string str)
{
  map<string,double>::iterator it = m_smin_map.find(str);
  if (it!=m_smin_map.end())
    if (it->second>=0.) return it->second;

  int length = str.length();
  int *legs = new int[length];
  std::vector<char> pr(length);
  for (int i=0;i<length;i++) {
    pr[i] = str[i];
    if (pr[i]<58) legs[i]=pr[i]-48;
    else legs[i]=pr[i]-55;
  }
  double sc = 0.;

  if (length==1) {
    m_smin_map[str] = sc = sqr(fl[legs[0]].SelMass());
    delete[] legs;
    return sc;
  }

  if (length==2) {
    m_smin_map[str] = sc = scut[legs[0]][legs[1]];
    delete[] legs;
    return sc;
  }

  string help("0"), help2("");
  for (int i=0;i<length;i++) {
    help[0] = GetIndexID(legs[i]);
    sc += Getscut(help);
  }
  sc *= 2.-(double)length;
  help = string("00");
  for (int i=0;i<length;i++) {
    for (int j=i+1;j<length;j++) {
      help[0] = GetIndexID(legs[i]);
      help[1] = GetIndexID(legs[j]);
      sc += Getscut(help);
    }
  }

  std::vector<char> pl(length,' ');
  for (int i(1);i<=length/2;++i) sc=Max(sc,Getscut(pl,pr,i,0,-1));
  
  m_smin_map[str] = sc;
  delete[] legs;
  return sc;
}

void Cut_Data::Setscut(std::string str,double d)
{
  m_smin_map[str]=d;
}

void Cut_Data::Update(double sprime,double y) 
{
  // reset cuts to lab values
  Reset(false);
  // boost from lab to cms
  double chy(cosh(y)), shy(sinh(y));
  Poincare cms[2]={Poincare(Vec4D(chy,0.0,0.0,shy)),
		   Poincare(Vec4D(chy,0.0,0.0,-shy))};
  for (int a=0;a<2;++a) {
    for (int i=2;i<ncut;i++) {
      if (cosmax[a][i]<1.0 && !fl[i].IsMassive()) {
	Vec4D help(1.0,sqrt(1.0-sqr(cosmax[a][i])),0.0,cosmax[a][i]);
	cms[a].Boost(help);
	cosmax[a][i]=cosmax[i][a]=help[3]/help[0];
      } 
      else cosmax[a][i] = cosmax[i][a] = 1.0;
      if (cosmin[a][i]>-1.0 && !fl[i].IsMassive()) {
	Vec4D help(1.0,sqrt(1.0-sqr(cosmin[a][i])),0.0,cosmin[a][i]);
	cms[a].Boost(help);
	cosmin[a][i]=cosmin[i][a]=help[3]/help[0];
      } 
      else cosmin[a][i] = cosmin[i][a] = -1.0;
      double ct=sqrt(1.0-(sqr(etmin[i])-sqr(fl[i].SelMass()))/
		     (sprime/4.0-sqr(fl[i].SelMass())));
      if (etmin[i]<fl[i].SelMass()) ct=1.0;
      cosmax[i][a]=cosmax[a][i]=Min(cosmax[a][i],ct);
      cosmin[i][a]=cosmin[a][i]=Max(cosmin[a][i],-ct);
    }
  }
}
