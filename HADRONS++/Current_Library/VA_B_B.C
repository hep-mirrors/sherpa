#include "HADRONS++/Current_Library/VA_B_B.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"


using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

#include "HADRONS++/Current_Library/VA_B_B_HO.C"
#include "HADRONS++/Current_Library/VA_B_B_ST.C"

VA_B_B_FFs::FormFactor_Base::~FormFactor_Base()
{
}

void VA_B_B::SetModelParameters( struct GeneralModel model )
{
  m_V_CKM    = model("V_CKM",1.0);
  switch(int(model("mode",-1.0)+0.5)) {
    case 0:
    case 2:
    case 10:
    case 11:
    case 12:
      m_unnatural = false;
      break;
    case 1:
    case 13:
    case 14:
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
      p_ff = new VA_B_B_FFs::HO(model,p_masses,"HONR");
      msg_Tracking()<<"    Using HONR form factor model for "<<m_name<<std::endl;
      break;
    case 2:
      p_ff = new VA_B_B_FFs::HO(model,p_masses,"HOSR");
      msg_Tracking()<<"    Using HOSR form factor model for "<<m_name<<std::endl;
      break;
    case 3:
      p_ff = new VA_B_B_FFs::ST(model,p_masses,"STNR");
      msg_Tracking()<<"    Using STNR form factor model for "<<m_name<<std::endl;
      break;
    case 4:
      p_ff = new VA_B_B_FFs::ST(model,p_masses,"STSR");
      msg_Tracking()<<"    Using STSR form factor model for "<<m_name<<std::endl;
      break;
    default:
      msg_Error()<<METHOD<<": You chose a form factor model which does not "
          <<"exist for current "<<m_name<<". Aborting."<<std::endl;
      Abort();
  }
  /*
  Vec4D p0(Vec4D(p_masses[0],0,0,0));
  Vec4D p1(Vec4D(p_masses[1],0,0,0));
  p_ff->CalcFFs(p0, p1);
  cout<<setw(18)<<m_flavs[p_i[0]]<<" --> "<<setw(18)<<m_flavs[p_i[1]];
  cout<<"\tF1 = "<<p_ff->V1();
  cout<<"\tF2 = "<<p_ff->V2();
  cout<<"\tF3 = "<<p_ff->V3();
  cout<<"\tG1 = "<<p_ff->A1();
  cout<<"\tG2 = "<<p_ff->A2();
  cout<<"\tG3 = "<<p_ff->A3()<<endl;
  */
}

VA_B_B::~VA_B_B() {
  delete p_ff;
}

void VA_B_B::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  p_ff->CalcFFs(moms[p_i[0]], moms[p_i[1]]);
  
  Complex cR1 = m_v*p_ff->V1()+m_a*p_ff->A1();
  Complex cL1 = m_unnatural ? -m_v*p_ff->V1()+m_a*p_ff->A1() : m_v*p_ff->V1()-m_a*p_ff->A1();
  Complex cR2 = m_v*p_ff->V2()+m_a*p_ff->A2();
  Complex cL2 = m_unnatural ? -m_v*p_ff->V2()+m_a*p_ff->A2() : m_v*p_ff->V2()-m_a*p_ff->A2();
  Complex cR3 = m_v*p_ff->V3()+m_a*p_ff->A3();
  Complex cL3 = m_unnatural ? -m_v*p_ff->V3()+m_a*p_ff->A3() : m_v*p_ff->V3()-m_a*p_ff->A3();

  XYZFunc F(moms, m_flavs, m_anti, p_i);
  for(int h0=0; h0<2; h0++) {
    for(int h1=0; h1<2; h1++) {
      // 1 is "the barred spinor" in the current, 0 is the not-barred one
      // (other way around than VA_F_F, maybe unify?)
      Vec4C amp(Vec4D(0.0, 0.0, 0.0, 0.0));
      amp += F.L(1,h1, 0,h0, cR1, cL1);
//       amp += moms[p_i[0]]/p_masses[0]*F.Y(1,h1, 0,h0, cR2, cL2);
//       amp += moms[p_i[1]]/p_masses[1]*F.Y(1,h1, 0,h0, cR3, cL3);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert( m_V_CKM*amp ,spins );
    }
  }
}

DEFINE_CURRENT_GETTER(VA_B_B,"VA_B_B")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_B_B>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<std::endl;
}
