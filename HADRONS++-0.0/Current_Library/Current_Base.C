#include "Current_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE HADRONS::Current_Base
#define PARAMETER_TYPE HADRONS::Flavour_Info
#include "Getter_Function.C"

using namespace ATOOLS;
using namespace HADRONS;
using namespace std;

Current_Base::Current_Base(const ATOOLS::Flavour* flavs, int n, int* decayindices, string name) :
    m_n(n), m_name(name)
{
  p_decayindices = new int[n];
  for(int i=0; i<m_n; i++) {
    p_decayindices[i] = decayindices[i];
  }
  p_flavs = new ATOOLS::Flavour[n];
  p_moms = new ATOOLS::Vec4D[n];
  p_masses = new double[n];
  p_masses2 = new double[n];
  p_intspins = new int[n];
  m_spin_combinations = 1;
  for(int currentindex=0; currentindex<m_n; currentindex++) {
    int decayindex = p_decayindices[currentindex];
    p_flavs[currentindex]   = flavs[decayindex];
    p_masses[currentindex]  = p_flavs[currentindex].PSMass();
    p_masses2[currentindex] = p_masses[currentindex]*p_masses[currentindex];
    p_intspins[currentindex]   = p_flavs[currentindex].IntSpin();
    m_spin_combinations *= (p_intspins[currentindex]+1);
  }
  p_results = new ComplexVec4D[m_spin_combinations];
  msg.Tracking()<<"  Initialized "<<m_name<<" current with "
                <<m_spin_combinations<<" spin combinations"<<endl;
  for(int i=0; i<m_n; i++) {
    msg.Debugging()<<"    flavs["<<i<<"]="<<p_flavs[i]<<endl;
    msg.Debugging()<<"    decayindices["<<i<<"]="<<p_decayindices[i]<<endl;
  }
}


Current_Base::~Current_Base()
{
  if (p_results!=NULL) delete[] p_results; p_results=NULL;
  if (p_flavs!=NULL) delete[] p_flavs; p_flavs=NULL;
  if (p_moms!=NULL) delete[] p_moms; p_moms=NULL;
  if (p_masses!=NULL) delete[] p_masses; p_masses=NULL;
  if (p_masses2!=NULL) delete[] p_masses2; p_masses2=NULL;
  if (p_intspins!=NULL) delete[] p_intspins; p_intspins=NULL;
  if (p_decayindices!=NULL) delete[] p_decayindices; p_decayindices=NULL;
}


const ComplexVec4D& Current_Base::operator[] (int i) const
{
    return p_results[i];
}


void Current_Base::SetMomsAndk0n ( const Vec4D * moms, const int k0n )
{
  for(int i=0;i<m_n;i++) {
    p_moms[i]=moms[p_decayindices[i]];
  }
  m_k0n=k0n;
}


void HADRONS::ContractCurrents( HADRONS::Current_Base* c1,
                                HADRONS::Current_Base* c2,
                                const ATOOLS::Vec4D * moms,
                                const int k0n,
                                double mefactor,
                                vector<Complex>* ampls,
                                vector<pair<int,int> >* sctindices)
{
  ampls->clear();
  sctindices->clear();
  c1->SetMomsAndk0n(moms,k0n);
  c1->Calc();
  c2->SetMomsAndk0n(moms,k0n);
  c2->Calc();
  msg.Debugging()<<(*c1)<<endl;
  msg.Debugging()<<(*c2)<<endl;
  msg.Debugging()<<"  mefactor="<<mefactor<<endl;
  for(int i=0; i<c1->m_spin_combinations; i++) {
    for(int j=0; j<c2->m_spin_combinations; j++) {
      ampls->push_back( mefactor*(c1->p_results[i]*c2->p_results[j]) );
      msg.Debugging()<<"  currents ampls pushing ["<<i<<","<<j<<"]:"
                     <<mefactor*(c1->p_results[i]*c2->p_results[j])<<endl;
    }
  }

  for(int i=c2->m_n-1; i>=0; i--) {
    if(c2->p_intspins[i]!=0) {
      sctindices->push_back( pair<int,int>(c2->p_decayindices[i],c2->p_intspins[i]) );
      msg.Debugging()<<"  currents sctindices pushing: "<<c2->p_decayindices[i]
                     <<" ["<<c2->p_intspins[i]<<"]"<<endl;
    }
  }
  for(int i=c1->m_n-1; i>=0; i--) {
    if(c1->p_intspins[i]!=0) {
      sctindices->push_back( pair<int,int>(c1->p_decayindices[i],c1->p_intspins[i]) );
      msg.Debugging()<<"  currents sctindices pushing: "<<c1->p_decayindices[i]
                     <<" ["<<c1->p_intspins[i]<<"]"<<endl;
    }
  }
    // example for n·n spin 1/2 particles
    // indexpicture in SCT:    c1[n] ... c1[0]  c2[n] ... c2[0]
    //                           +         +      +         +
    //                           +         +      +         -
    //                           +         +      -         +
    //                           +         +      -         -
    //                           +         -      +         +
    //                           +         -      +         -
    //                           +         -      -         +
    //                           +         -      -         -
    //                           -         +      +         +
    //                           -         +      +         -
    //                           -         +      -         +
    //                           -         +      -         -
    //                           -         -      +         +
    //                           -         -      +         -
    //                           -         -      -         +
    //                           -         -      -         -

}


ostream& operator<< (ostream& s, const HADRONS::Current_Base& cb)
{
  s<<cb.m_name<<" current with "<<cb.m_spin_combinations<<" spin combinations:"<<endl;
  for(int i=0; i<cb.m_spin_combinations; i++) {
    s<<"  "<<cb.p_results[i]<<endl;
  }
  return s;
  
}
