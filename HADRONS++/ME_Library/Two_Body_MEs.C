#include "HADRONS++/ME_Library/Two_Body_MEs.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void Baryon_Diquark_Quark::SetModelParameters( GeneralModel _md )
{
  m_Vxx_decay = _md("Vxx_decay",1.0);
  m_Vxx_production = _md("Vxx_production",1.0);
  m_GF = _md("GF",1.0);
}

void Baryon_Diquark_Quark::operator()(
                      const Vec4D             * p,
                      METOOLS::Spin_Amplitudes * amps)
{
  vector<pair<int,int> > spins(3);
  for(int h0=0; h0<m_flavs[p_i[0]].IntSpin()+1;++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for( int h1=0; h1<m_flavs[p_i[1]].IntSpin()+1; h1++ ) {
      spins[1] = make_pair(p_i[1],h1);
      for( int h2=0; h2<m_flavs[p_i[2]].IntSpin()+1; h2++ ) {
	spins[2] = make_pair(p_i[2],h2);
	amps->Add(Complex(1.0,0.0),spins);
      }
    }
  }
}

bool Baryon_Diquark_Quark::SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q, int n_g)
{
  int pos = m_anti ? 2 : 1;
  outparts[p_i[2]-1]->SetFlow(pos,-1);
  outparts[p_i[1]-1]->SetFlow(3-pos,outparts[p_i[2]-1]->GetFlow(pos));
  return true;
}

DEFINE_ME_GETTER(Baryon_Diquark_Quark,Baryon_Diquark_Quark_Getter,"Baryon_Diquark_Quark")

void Baryon_Diquark_Quark_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}
