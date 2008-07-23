#include "Cut_Data.H"
#include "Run_Parameter.H"
#include "MyStrStream.H"
#include "Poincare.H"
#include "Message.H"

using namespace ATOOLS;
using namespace std;


std::ostream & ATOOLS::operator<<(std::ostream & s , Cut_Data & cd)
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
  energymax = 0;
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
    delete[] cosmax_save[i];
    delete[] scut[i];
    delete[] scut_save[i];
  }
  delete[] cosmin;
  delete[] cosmax;
  delete[] cosmax_save;
  delete[] scut;
  delete[] scut_save;
  delete[] energymin;
  delete[] energymin_save;
  delete[] energymax;
  delete[] etmin;
}

void Cut_Data::Init(int _ncut,Flavour * _fl) {
  if (energymin != 0) return;
  smin = 0.;
  ncut           = _ncut;
  fl             = _fl;
  double E       = rpa.gen.Ecms();
  energymin      = new double[ncut];
  energymin_save = new double[ncut];
  energymax      = new double[ncut];
  etmin          = new double[ncut];
  cosmin         = new double*[ncut];
  cosmax         = new double*[ncut];
  cosmax_save    = new double*[ncut];
  scut           = new double*[ncut];
  scut_save      = new double*[ncut];

  for (int i=0;i<ncut;i++) {
    cosmin[i]      = new double[ncut];
    cosmax[i]      = new double[ncut];
    cosmax_save[i] = new double[ncut];
    scut[i]        = new double[ncut];
    scut_save[i]   = new double[ncut];
    energymin[i]   = Max(0.,fl[i].SelMass());
    if (fl[i].IsKK()) energymin[i] = 0.;
    smin += energymin_save[i] = energymin[i];
    energymax[i]   = E;
    etmin[i]       = 0.;
  }
  smin = sqr(smin);

  for (int i=0;i<ncut;i++) {
    for (int j=i;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.;
      cosmax[i][j] = cosmax[j][i] = cosmax_save[i][j] = 1.;
      /*      double sc =
	+sqr(fl[i].SelMass())+sqr(fl[j].SelMass())
	+2.*energymin[i]*energymin[j]
	-2.*sqrt(dabs(sqr(energymin[i])-sqr(fl[i].SelMass())))
	*sqrt(dabs(sqr(energymin[j])-sqr(fl[j].SelMass())))
	*cosmax[i][j];*/
      scut[i][j] = scut[j][i] = scut_save[i][j] = 1.e-12*sqr(rpa.gen.Ecms());
	//Max(sc,1.e-12*sqr(rpa.gen.Ecms()));
    }
  }  
}

void Cut_Data::Complete()
{
  for (int i=0;i<ncut;i++) {
    for (int j=i+1;j<ncut;j++) {
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

/*void Cut_Data::Init(Cut_Data * _cuts,Flavour * _fl) {
  ncut      = _cuts->ncut;
  fl        = _fl;
  if (!energymin) {
    energymin           = new double[ncut];
    energymax           = new double[ncut];
    cosmin              = new double*[ncut];
    cosmax              = new double*[ncut];
    scut                = new double*[ncut];
    scut_save           = new double*[ncut];
    for (int i=0;i<ncut;i++) {
      energymin[i]      = _cuts->energymin[i];
      energymax[i]      = _cuts->energymax[i];
      cosmin[i]         = new double[ncut];
      cosmax[i]         = new double[ncut];
      scut[i]           = new double[ncut];
      scut_save[i]      = new double[ncut];
      for (int j=i;j<ncut;j++) {
	cosmin[i][j]    = cosmin[j][i]    = _cuts->cosmin[i][j];
	cosmax[i][j]    = cosmax[j][i]    = _cuts->cosmax[i][j];
	scut[i][j]      = scut[j][i]      = _cuts->scut[i][j];
	scut_save[i][j] = scut_save[j][i] = scut[i][j];
      }
    }
  }
  else {
    for (int i=0;i<ncut;i++) {
      energymin[i]      = Min(energymin[i],_cuts->energymin[i]);
      energymax[i]      = Max(energymax[i],_cuts->energymax[i]);
      for (int j=i;j<ncut;j++) {
	cosmin[i][j]    = cosmin[j][i]    = Min(cosmin[i][j],_cuts->cosmin[i][j]);
	cosmax[i][j]    = cosmax[j][i]    = Max(cosmax[i][j],_cuts->cosmax[i][j]);
	scut[i][j]      = scut[j][i]      = Min(scut[i][j],_cuts->scut[i][j]);
	scut_save[i][j] = scut_save[j][i] = scut[i][j];
      }
    }
  }
  }*/




void Cut_Data::Update(double sprime,double y) 
{
  // boost cosmax from LAB to CMS
  double E1     = std::exp(y);  
  double E2     = std::exp(-y);  
  Poincare Forward(Vec4D(E1+E2,0.,0.,E1-E2));
  Poincare Backward(Vec4D(E1+E2,0.,0.,E2-E1));
  Vec4D help;
  for (int i=2;i<ncut;i++) {


    if (cosmax[0][i]<1.) {
      if (!fl[i].IsMassive() || y>=0.) {
	help = Vec4D(1.,sqrt(1.-sqr(cosmax[0][i])),0.,cosmax[0][i]);
	Forward.Boost(help);
	cosmax[0][i] = cosmax[i][0] = help[3]/help[0];
      } 
      else cosmax[0][i] = cosmax[i][0] = 1.;  // No better estimate for massive particles
      if (!fl[i].IsMassive() || y<=0.) {
	help = Vec4D(1.,sqrt(1.-sqr(cosmax[1][i])),0.,cosmax[1][i]);
	Backward.Boost(help);
	cosmax[1][i] = cosmax[i][1] = help[3]/help[0];
      }
      else cosmax[1][i] = cosmax[i][1] = 1.;  // No better estimate for massive particles
    }

    double ct= sqrt(1.0-(sqr(etmin[i])-sqr(fl[i].SelMass()))/(sprime/4.0-sqr(fl[i].SelMass())));
    if (etmin[i]<fl[i].SelMass()) ct=1.;
    cosmax[i][0] = cosmax[0][i] = Min(cosmax[0][i],ct);
    cosmax[i][1] = cosmax[1][i] = Min(cosmax[1][i],ct);
    cosmin[i][1] = cosmin[i][0] = cosmin[1][i] = cosmin[0][i] = -ct;
  }
  
  //for (int i=2;i<ncut;i++) 
    //cosmax[0][i] = cosmax[1][i] =1.;
  /*  double E = sqrt(sprime);
  for (int i=0;i<ncut;i++) {
    energymin[i]   = Max(0.,fl[i].SelMass());
    energymax[i]   = E;
    for (int j=i+1;j<ncut;j++) {
      cosmin[i][j] = cosmin[j][i] = -1.;
      cosmax[i][j] = cosmax[j][i] =  1.;
    }
  }
  for (int i=0;i<ncut;i++) {
    for (int j=i+1;j<ncut;j++) {
      scut[i][j] = scut[j][i] = Max(scut_save[i][j],1.e-12*sprime);
    }
    }  */
}







