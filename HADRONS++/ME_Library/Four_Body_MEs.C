#include "HADRONS++/ME_Library/Four_Body_MEs.H"
#include "ATOOLS/Org/Message.H"
#include "HADRONS++/Main/Tools.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Math/Random.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

void QQ_QQQQ_Spectator::SetModelParameters( GeneralModel _md )
{
  int decayerquark(0);
  if(m_flavs[p_i[0]].Kfcode()==kf_B || m_flavs[p_i[0]].Kfcode()==kf_B_plus ||
     m_flavs[p_i[0]].Kfcode()==kf_B_s || m_flavs[p_i[0]].Kfcode()==kf_B_c)
    decayerquark=-5;
  if(m_flavs[p_i[0]].IsAnti()) decayerquark=-decayerquark;
  decayerquark = int(_md("decayerquark",decayerquark));
  kf_code decayerkfc = (kf_code) abs(decayerquark);
  m_decayer = Flavour(decayerkfc, decayerquark<0);

  m_Vxx_decay=1.0;
  if(m_decayer.Kfcode()==kf_b) {
    if(m_flavs[p_i[1]].Kfcode()==kf_c) {
      m_Vxx_decay=Tools::Vcb;
    }
    else if(m_flavs[p_i[1]].Kfcode()==kf_u) {
      m_Vxx_decay=Tools::Vub;
    }
  }
  m_Vxx_decay = _md("Vxx_decay",m_Vxx_decay);

  m_Vxx_production=1.0;
  if((m_flavs[p_i[2]].Kfcode()==kf_c && m_flavs[p_i[3]].Kfcode()==kf_s) ||
     (m_flavs[p_i[2]].Kfcode()==kf_s && m_flavs[p_i[3]].Kfcode()==kf_c))
    m_Vxx_production=Tools::Vcs;
  else if((m_flavs[p_i[2]].Kfcode()==kf_c && m_flavs[p_i[3]].Kfcode()==kf_d) ||
          (m_flavs[p_i[2]].Kfcode()==kf_d && m_flavs[p_i[3]].Kfcode()==kf_c))
    m_Vxx_production=Tools::Vcd;
  else if((m_flavs[p_i[2]].Kfcode()==kf_u && m_flavs[p_i[3]].Kfcode()==kf_s) ||
          (m_flavs[p_i[2]].Kfcode()==kf_s && m_flavs[p_i[3]].Kfcode()==kf_u))
    m_Vxx_production=Tools::Vus;
  else if((m_flavs[p_i[2]].Kfcode()==kf_u && m_flavs[p_i[3]].Kfcode()==kf_d) ||
          (m_flavs[p_i[2]].Kfcode()==kf_d && m_flavs[p_i[3]].Kfcode()==kf_u))
    m_Vxx_production=Tools::Vud;
  
  m_Vxx_production = _md("Vxx_production",m_Vxx_production);
  m_GF = _md("GF",8.24748e-6);
  
  m_cR_decay   = _md("v_decay",1.)+_md("a_decay",-1.);
  m_cL_decay   = _md("v_decay",1.)-_md("a_decay",-1.);
  m_cR_production   = _md("v_production",1.)+_md("a_production",-1.);
  m_cL_production   = _md("v_production",1.)-_md("a_production",-1.);

  m_colourflip_ratio = _md("colourflip_ratio",0.0);
}

void QQ_QQQQ_Spectator::operator()(
                      const Vec4D             * p,
                      METOOLS::Spin_Amplitudes * amps)
{
  double factor = m_GF*m_Vxx_decay*m_Vxx_production;
  Flavour partonflavs[] = {m_decayer, m_flavs[p_i[1]], m_flavs[p_i[2]], m_flavs[p_i[3]]};
  Vec4D partonmoms[] = {p[p_i[1]]+p[p_i[2]]+p[p_i[3]], 
                        p[p_i[1]], p[p_i[2]], p[p_i[3]]};
  XYZFunc F(4,partonmoms,partonflavs,Tools::k0,m_anti);
  
  vector<pair<int,int> > spins(5);
  spins[0] = make_pair(p_i[0],0);
  for( int h1=0; h1<2; h1++ ) { // direct quark line
    spins[1] = make_pair(p_i[1],h1);
    for( int h2=0; h2<2; h2++ ) { // anti quark from W
      spins[2] = make_pair(p_i[2],h2);
      for( int h3=0; h3<2; h3++ ) { // quark from W
        spins[3] = make_pair(p_i[3],h3);
        for( int h4=0; h4<2; h4++ ) { // spectator quark
          spins[4] = make_pair(p_i[4],h4);
          int h0 = 1-h4; // flip for decayer

          Complex amp = factor * F.Z(0,h0, 1,h1, 3,h3, 2,h2,
                                     m_cR_decay, m_cL_decay, m_cR_production, m_cL_production);
          amps->Add(amp,spins);
        }
      }
    }
  }
}

bool QQ_QQQQ_Spectator::SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q, int n_g)
{
  int pos = m_anti ? 2 : 1;
  if(ran->Get()<m_colourflip_ratio) { // colourflip
    outparts[p_i[4]-1]->SetFlow(pos,-1);
    outparts[p_i[2]-1]->SetFlow(3-pos,outparts[p_i[4]-1]->GetFlow(pos));
    
    outparts[p_i[3]-1]->SetFlow(pos,-1);
    outparts[p_i[1]-1]->SetFlow(3-pos,outparts[p_i[3]-1]->GetFlow(pos));
  }
  else {
    outparts[p_i[4]-1]->SetFlow(pos,-1);
    outparts[p_i[1]-1]->SetFlow(3-pos,outparts[p_i[4]-1]->GetFlow(pos));

    outparts[p_i[3]-1]->SetFlow(pos,-1);
    outparts[p_i[2]-1]->SetFlow(3-pos,outparts[p_i[3]-1]->GetFlow(pos));
  }
  return true;
}

DEFINE_ME_GETTER(QQ_QQQQ_Spectator,QQ_QQQQ_Spectator_Getter,"QQ_QQQQ_Spectator")

void QQ_QQQQ_Spectator_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B^{+} \\rightarrow \\bar{c} u u \\bar{d} $ \n\n"
    <<"Order: 0 = Scalar ($B^{+}$), 1 = quark from decay line ($\\bar{c}$), "
    <<"2 = anti quark from W ($\\bar{d}$), 3 = quark from W ($u$), "
    <<"4 = spectator quark ($u$) \n\n"
    <<"\\[ \\mathcal{M} =  \\]"
    <<endl;
}



void Baryon_Diquark_3Quarks::SetModelParameters( GeneralModel _md )
{
  m_Vxx_decay = _md("Vxx_decay",1.0);
  m_Vxx_production = _md("Vxx_production",1.0);
  m_GF = _md("GF",8.033333333e-6);
}

void Baryon_Diquark_3Quarks::operator()(
                      const Vec4D             * p,
                      METOOLS::Spin_Amplitudes * amps)
{
  vector<pair<int,int> > spins(5);
  for(int h0=0; h0<m_flavs[p_i[0]].IntSpin()+1;++h0) {
    spins[0] = make_pair(p_i[0],h0);
    for( int h1=0; h1<m_flavs[p_i[1]].IntSpin()+1; h1++ ) {
      spins[1] = make_pair(p_i[1],h1);
      for( int h2=0; h2<m_flavs[p_i[2]].IntSpin()+1; h2++ ) {
	spins[2] = make_pair(p_i[2],h2);
	for( int h3=0; h3<m_flavs[p_i[3]].IntSpin()+1; h3++ ) {
	  spins[3] = make_pair(p_i[3],h3);
	  for( int h4=0; h4<m_flavs[p_i[4]].IntSpin()+1; h4++ ) {
	    spins[4] = make_pair(p_i[4],h4);
	    amps->Add(Complex(1.0,0.0),spins);
	  }
	}
      }
    }
  }
}

bool Baryon_Diquark_3Quarks::SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q, int n_g)
{
  int pos = m_anti ? 2 : 1;
  outparts[p_i[2]-1]->SetFlow(pos,-1);
  outparts[p_i[1]-1]->SetFlow(3-pos,outparts[p_i[2]-1]->GetFlow(pos));
  
  outparts[p_i[3]-1]->SetFlow(pos,-1);
  outparts[p_i[4]-1]->SetFlow(3-pos,outparts[p_i[3]-1]->GetFlow(pos));
  return true;
}

DEFINE_ME_GETTER(Baryon_Diquark_3Quarks,Baryon_Diquark_3Quarks_Getter,"Baryon_Diquark_3Quarks")

void Baryon_Diquark_3Quarks_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<endl;
}
