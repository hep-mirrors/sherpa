#include "Durham_Algorithm.H"
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


Durham_Algorithm::Durham_Algorithm(double rmin) : 
  m_matrixsize(0), p_imap(NULL), p_jets(NULL)
{

}

Durham_Algorithm::~Durham_Algorithm()
{
  if (p_imap) {
//     for (int i=0;i<m_matrixsize;++i) {
//       delete [] p_ktij[i];
//     }
//     delete [] p_ktij;
    delete [] p_imap;
    //    delete [] p_kis;
  }
}


void Durham_Algorithm::AddToKtlist(double kt2) {
   if (p_kts) {
     if (p_kts->size()>0 && p_kts->back()==kt2) std::cout<<"WARNING something fishy in AddToKtlist "<<std::endl;
     p_kts->push_back(sqr(kt2));
   }
}

void Durham_Algorithm::AddToJetlist(const Vec4D & mom, bool bf) {
  if (p_jets) {
    if(!bf) p_jets->push_back(new Particle(p_jets->size(),Flavour(kf::jet),mom));
    else    p_jets->push_back(new Particle(p_jets->size(),Flavour(kf::bjet),mom));
  }
}



bool Durham_Algorithm::ConstructJets(const Particle_List * pl, Particle_List * jets, 
				     std::vector<double> * kts, double ycut) {
  p_jets = jets;
  p_kts  = kts;
  m_ycut = ycut;
  

  Vec4D momsum=Vec4D(0.,0.,0.,0.);
  // create vector list from particle list
  int n=0;
  Vec4D * moms = new Vec4D[pl->size()];
  bool * bflag = new bool[pl->size()];
  for (Particle_List::const_iterator it=pl->begin(); it!=pl->end();++it) {
    momsum+=(*it)->Momentum();
    if (!(*it)->Flav().IsLepton()) {
      moms[n]  = ((*it)->Momentum()); 
      bflag[n] = (((*it)->Flav()).Kfcode()==kf::b)&& !m_bflag;
      ++n;
    }
  }
  m_sprime = momsum.Abs2();

  // cluster
  Ymin(moms,bflag,n);
  delete [] moms;
  delete [] bflag;
  // finalize (sort and release used containers)

  SortE();

  p_jets=0;
  p_kts =0;

  return true;
}


void Durham_Algorithm::Init(int size) 
{
  PROFILE_HERE;
  if (size>m_matrixsize ) {
    if (p_imap) {
      delete [] p_imap;
    }
    m_matrixsize = size;
    p_imap = new int[size];
  }
  for (int i=0;i<size;++i) p_imap[i] = i;
}


void Durham_Algorithm::Ymin(Vec4D * p, bool * bf, int n)
{
  if (n==0) return;
  if (n==1) {
    AddToJetlist(p[0],bf[0]);
    return;// dmin;
  }

  Init(n);

  int hit = 0;

  while (n>1) {
    int ii=0, jj=0;
    double ymin = 1.;//m_ycut;  
    for (int i=1;i<n;++i) {
      for (int j=0;j<i;++j) {
	double y = Y12(p[p_imap[i]],p[p_imap[j]]);
	if (y<ymin) { ymin=y; ii=i; jj=j;}
      }
    }
    if (!hit && ymin>=m_ycut) {
      hit = 1;
      for (int i=0;i<n;++i) AddToJetlist(p[p_imap[i]],bf[p_imap[i]]);
    }
    // combine precluster
    AddToKtlist(ymin);
    p[p_imap[jj]]+=p[p_imap[ii]];
    bf[p_imap[jj]] = bf[p_imap[jj]]||bf[p_imap[ii]];      
    --n;
    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];
  }
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
    if (Durham_Algorithm::Kt2(a->Momentum()) > Durham_Algorithm::Kt2(b->Momentum())) return 1;
    return 0;
  }
};

void Durham_Algorithm::SortE()
{
  if (p_jets) {
    std::sort(p_jets->begin(), p_jets->end(),Order_E());
  }
}
void Durham_Algorithm::SortPT()
{
  if (p_jets) {
    std::sort(p_jets->begin(), p_jets->end(),Order_PT());
  }
}

double Durham_Algorithm::DCos12(const Vec4D & p1,const Vec4D & p2) const
{
  return Vec3D(p1)*Vec3D(p2)/(Vec3D(p1).Abs()*Vec3D(p2).Abs());
}

double Durham_Algorithm::Y12(const Vec4D & p1, const Vec4D & p2) const
{
  return 2.*sqr(Min(p1[0],p2[0]))*(1.-DCos12(p1,p2))/m_sprime;
}

