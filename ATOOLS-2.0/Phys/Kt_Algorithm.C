#include "Kt_Algorithm.H"
#include "Particle_List.H"
#include "Message.H"

#include <iomanip>
#include <algorithm>

#ifdef PROFILE__Analysis_Phase
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace ATOOLS;
using namespace std;


Kt_Algorithm::Kt_Algorithm(double rmin) : 
  m_matrixsize(0), p_ktij(NULL), p_imap(NULL), p_kis(NULL), p_jets(NULL),
  p_kts(NULL)
{
  std::cout<<" Init Kt_Algorithm "<<std::endl;
}

Kt_Algorithm::~Kt_Algorithm()
{
  if (p_ktij) {
    for (int i=0;i<m_matrixsize;++i) {
      delete [] p_ktij[i];
    }
    delete [] p_ktij;
    delete [] p_imap;
    delete [] p_kis;
  }
}


void Kt_Algorithm::AddToKtlist(double kt2) {
  if (p_kts) {
    if (p_kts->size()>0) if (p_kts->back()==kt2) std::cout<<" WARNING something fishy in AddToKtlist "<<std::endl;
    p_kts->push_back(kt2);
  }
}

void Kt_Algorithm::AddToJetlist(const Vec4D & mom, bool bf) {
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf::jet),mom));
    else    p_jets->push_back(new Particle(p_jets->size(),Flavour(kf::bjet),mom));
  }
}



bool Kt_Algorithm::ConstructJets(const Particle_List * pl, Particle_List * jets, std::vector<double> * kts, double rmin) {
  //cout<<" in Kt_Algorithm::ConstructJets(("<<pl->size()<<"),("<<jets->size()<<"),"<<kts<<","<<rmin<<") "<<endl;
  // assume empty containers
  p_jets = jets;
  p_kts  = kts;
  m_r2min = sqr(rmin);

  // create vector list from particle list
  int n=0;
  //  cout<<" create vector "<<endl;
  Vec4D * moms = new Vec4D[pl->size()];
  bool * bflag = new bool[pl->size()];
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    if (!(*it)->Flav().IsLepton()) {
      moms[n]  = ((*it)->Momentum()); 
      bflag[n] = (((*it)->Flav()).Kfcode()==kf::b)&& !m_bflag;
      ++n;
    }
  }

  //cout<<" cluster vector "<<n<<endl;
  // cluster
  Ktmin(moms,bflag,n);

  //cout<<" finish "<<endl;
  delete [] moms;
  delete [] bflag;
  // finalize (sort and release used containers)

  SortPT();
  //  cout<<" done "<<endl;

  p_jets=0;
  p_kts =0;

  return true;
}


void Kt_Algorithm::Init(int size) 
{
  PROFILE_HERE;
  if (size>m_matrixsize ) {
    if (p_ktij) {
      for (int i=0;i<m_matrixsize;++i) delete [] p_ktij[i];
      delete [] p_ktij;
      delete [] p_imap;
      delete [] p_kis;
    }
    m_matrixsize = size;
    p_kis  = new double[size];
    p_imap = new int[size];
    p_ktij = new double*[size];
    for (int i=0;i<size;++i) p_ktij[i] = new double[size];
  }
  for (int i=0;i<size;++i) p_imap[i] = i;
}


double Kt_Algorithm::Ktmin(Vec4D * p, bool * bf, int n)
{
  if (n==0) return 0.;
  if (n==1) {
    AddToJetlist(p[0],bf[0]);
    double dmin=Kt2(p[0]);
    AddToKtlist(dmin);
    return dmin;;
  }

  Init(n);


  //cal first matrix
  int ii=0, jj=0;
  double dmin=Kt2(p[0]);
  double rmin=R2(p[0],p[1]);
  //cout<<"Ktmin start: "  << dmin<<" "<<rmin<<endl;
  {
  PROFILE_LOCAL(" first loop ");

  for (int i=0;i<n;++i) {
    double di = p_ktij[i][i] = Kt2(p[i]);
    if (di<dmin) { dmin=di; ii=i; jj=i;}
    for (int j=0;j<i;++j) {
      double dj  = p_ktij[j][j]; 
      double rij = p_ktij[i][j] = R2(p[i],p[j]);
      double dij = Min(di,dj)* rij /m_r2min;
      if (dij<dmin) {dmin=dij; ii=i; jj=j;}
      if (rij<rmin) rmin=rij;
    }
  }
  }
  //cout<<"Ktmin after: "  << dmin<<" "<<rmin<<"("<<ii<<","<<jj<<")"<<endl;

  // recalc matrix
  while (n>0) {
    PROFILE_LOCAL(" main loop ");

//     cout<<" winner is "<<ii<<","<<jj<<std::endl;
//     cout<<" current list is :"<<std::endl;
//     for (int i=0;i<n;++i) {
//       cout<<i<<" "<<p_imap[i]<<" : "<<p[p_imap[i]]<<endl;
//     }
//     cout.precision(6);
//     std::cout<<resetiosflags(ios::scientific);
//     for (int i=0; i<n;++i) {
//       for (int j=0; j<n;++j) {
// 	if (j<=i) cout<<setw(12)<<p_ktij[p_imap[i]][p_imap[j]];
// 	else      cout<<"            ";
//       }
//       cout<<endl;
//     }
// ----------------------------------------

    if (ii!=jj) {
      // combine precluster
      p[p_imap[jj]]+=p[p_imap[ii]];
      bf[p_imap[jj]] = bf[p_imap[jj]]||bf[p_imap[ii]];      
      AddToKtlist(dmin);
    }
    else {
      // add to jet list
      AddToJetlist(p[p_imap[jj]],bf[p_imap[jj]]);
      AddToKtlist(dmin);
    }

    --n;

    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];

    if (n==1) break;
    // update map (remove precluster)
    {
      PROFILE_LOCAL(" update loop ");

    
    // update matrix (only what is necessary)
    int jjx=p_imap[jj];
    p_ktij[jjx][jjx] = Kt2(p[jjx]);
    for (int j=0;j<jj;++j)   p_ktij[jjx][p_imap[j]] = R2(p[jjx],p[p_imap[j]]);
    for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = R2(p[p_imap[i]],p[jjx]);
    }
//     cout.precision(6);
//     std::cout<<resetiosflags(ios::scientific);
//     for (int i=0; i<n;++i) {
//       for (int j=0; j<n;++j) {
// 	if (j<=i) cout<<setw(12)<<p_ktij[p_imap[i]][p_imap[j]];
// 	else      cout<<"            ";
//       }
//       cout<<endl;
//     }
    // redetermin rmin and dmin
    ii=0, jj=0;
    {
      PROFILE_LOCAL(" second loop ");

    dmin=p_ktij[p_imap[0]][p_imap[0]];
    rmin=p_ktij[p_imap[1]][p_imap[0]];
    for (int i=0;i<n;++i) {
      int ix=p_imap[i];
      double di = p_ktij[ix][ix];
      if (di<dmin) { dmin=di; ii=i; jj=i;}
      for (int j=0;j<i;++j) {
	int jx=p_imap[j];
	double dj  = p_ktij[jx][jx];
	double rij = p_ktij[ix][jx];
	double dij = Min(di,dj)* rij/m_r2min;
	if (dij<dmin) {
// 	  cout<<" check: "<<i<<","<<j<<" dij="<<dij<<"   di="<<di<<"  dj="<<dj<<"  r2="<<rij<<endl;
	  dmin=dij; ii=i; jj=j;}
	if (rij<rmin) rmin=rij;
      }
    }
    }
  }

  // add remaining preclusters to jetlist
  for (int i=0;i<n;++i) {
    AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
    AddToKtlist(p_ktij[p_imap[i]][p_imap[i]]);
  }
  return dmin;
}

class Order_E {
public:
  int operator()(const Particle * a, const Particle * b) {
    if (a->Momentum()[0] > b->Momentum()[0]) return 1;
    return 0;
  }
};

class Order_PT {
public:
  int operator()(const Particle * a, const Particle * b) {
    if (Kt_Algorithm::Kt2(a->Momentum()) > Kt_Algorithm::Kt2(b->Momentum())) return 1;
    return 0;
  }
};

void Kt_Algorithm::SortE()
{
  if (p_jets) {
    std::sort(p_jets->begin(), p_jets->end(),Order_E());
  }
}
void Kt_Algorithm::SortPT()
{
  if (p_jets) {
    std::sort(p_jets->begin(), p_jets->end(),Order_PT());
  }
}

double Kt_Algorithm::CosDPhi12(const Vec4D & p1,const Vec4D & p2) const
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return (p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2);
}

double Kt_Algorithm::DCos12(const Vec4D & p1,const Vec4D & p2) const
{
  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Kt_Algorithm::DEta12(const Vec4D & p1, const Vec4D & p2) const
{
  // eta1,2 = -log(tan(theta_1,2)/2)   
  double c1=p1[3]/Vec3D(p1).Abs();
  double c2=p2[3]/Vec3D(p2).Abs();
  return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
}

double Kt_Algorithm::DPhi12(const Vec4D & p1,const Vec4D & p2) const
{
  double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
  double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
  return acos((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2));
}

Jet_Algorithm_Base::~Jet_Algorithm_Base()
{
}
