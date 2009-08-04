#include "HADRONS++/Current_Library/VA_B_B3.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"


using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

#include "HADRONS++/Current_Library/VA_B_B3_HO.C"
#include "HADRONS++/Current_Library/VA_B_B3_ST.C"

void VA_B_B3::SetModelParameters( struct GeneralModel model )
{
  m_V_CKM    = model("V_CKM",1.0);
  switch(int(model("mode",-1.0)+0.5)) {
  case 0:
  case 10:
  case 11:
    m_unnatural = false;
    break;
  case 1:
  case 12:
  case 13:
  case 14:
  case 15:
    m_unnatural = true;
    break;
  default:
    int temp = int(model("unnatural",-1.0)+0.5);
    if(temp==0 || temp==1) m_unnatural=bool(temp);
    else THROW(fatal_error, "\"unnatural\" not specified in "+m_name);
    break;
  }
  m_v        = model("v", 1.0);
  m_a        = model("a", -1.0);

  switch( int(model("FORM_FACTOR", 1)+0.5) ) {
    case 1:
      p_ff = new VA_B_B3_FFs::HO(model,p_masses,"HONR");
      msg_Tracking()<<"    Using HONR form factor model for "<<m_name<<std::endl;
      break;
    case 2:
      p_ff = new VA_B_B3_FFs::HO(model,p_masses,"HOSR");
      msg_Tracking()<<"    Using HOSR form factor model for "<<m_name<<std::endl;
      break;
    case 3:
      p_ff = new VA_B_B3_FFs::ST(model,p_masses,"STNR");
      msg_Tracking()<<"    Using STNR form factor model for "<<m_name<<std::endl;
      break;
    case 4:
      p_ff = new VA_B_B3_FFs::ST(model,p_masses,"STSR");
      msg_Tracking()<<"    Using STSR form factor model for "<<m_name<<std::endl;
      break;
    default:
      msg_Error()<<METHOD<<": You chose a form factor model which does not "
          <<"exist for current "<<m_name<<". Aborting."<<std::endl;
      abort();
  }
  /*
  Vec4D p0(Vec4D(p_masses[0],0,0,0));
  Vec4D p1(Vec4D(p_masses[1],0,0,0));
  p_ff->CalcFFs(p0, p1);
  cout<<setw(18)<<m_flavs[p_i[0]]<<" --> "<<setw(18)<<m_flavs[p_i[1]];
  cout<<"\tF1 = "<<p_ff->V1();
  cout<<"\tF2 = "<<p_ff->V2();
  cout<<"\tF3 = "<<p_ff->V3();
  cout<<"\tF4 = "<<p_ff->V4();
  cout<<"\tG1 = "<<p_ff->A1();
  cout<<"\tG2 = "<<p_ff->A2();
  cout<<"\tG3 = "<<p_ff->A3();
  cout<<"\tG4 = "<<p_ff->A4()<<endl;
  */
}

VA_B_B3::~VA_B_B3() {
  delete p_ff;
}

void VA_B_B3::Calc(const ATOOLS::Vec4D * moms)
{
  p_ff->CalcFFs(moms[p_i[0]], moms[p_i[1]]);

  Complex cR1 = m_v*p_ff->V1()+m_a*p_ff->A1();
  Complex cL1 = m_unnatural ?
    -m_v*p_ff->V1()+m_a*p_ff->A1() :
    m_v*p_ff->V1()-m_a*p_ff->A1();
  Complex cR2 = m_v*p_ff->V2()+m_a*p_ff->A2();
  Complex cL2 = m_unnatural ?
    -m_v*p_ff->V2()+m_a*p_ff->A2() :
    m_v*p_ff->V2()-m_a*p_ff->A2();
  Complex cR3 = m_v*p_ff->V3()+m_a*p_ff->A3();
  Complex cL3 = m_unnatural ?
    -m_v*p_ff->V3()+m_a*p_ff->A3() :
    m_v*p_ff->V3()-m_a*p_ff->A3();
  Complex cR4 = m_v*p_ff->V4()+m_a*p_ff->A4();
  Complex cL4 = m_unnatural ?
    -m_v*p_ff->V4()+m_a*p_ff->A4() :
    m_v*p_ff->V4()-m_a*p_ff->A4();

  XYZFunc F(2, moms, m_flavs, Tools::k0, m_anti);
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      // 1 is "the barred spinor" in the current, 0 is the not-barred one
      // (other way around than VA_F_F, maybe unify?)
      Vec4C amp(Vec4D(0.0,0.0,0.0,0.0));
      Vec4D w = moms[p_i[0]]/p_masses[0];
      amp += F.L31(1,h1, w, 0,h0, cR1, cL1);
//       amp += moms[p_i[0]]/p_masses[0]*w*F.Y31(1,h1, 0,h0, cR2, cL2);
//       amp += moms[p_i[1]]/p_masses[1]*w*F.Y31(1,h1, 0,h0, cR3, cL3);
//       amp += F.Y31(1,h1, 0,h0, cR4, cL4);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( m_V_CKM*amp ,spins );
    }
  }
}

DEFINE_CURRENT_GETTER(VA_B_B3,VA_B_B3_Getter,"VA_B_B3")

void VA_B_B3_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<std::endl;
}
