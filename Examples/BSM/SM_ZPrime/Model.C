#include "MODEL/SM/Model.H"

// define a kf_code for the Zprime
#define kf_Zp 32

namespace MODEL {

  class Standard_Model_Zprime: public Standard_Model {
  private:

    void FixZprimeParameters();  // <-- new: sets Zprime couplings
    void ParticleZprimeInit();   // <-- new: sets Zprime particle properties
    void InitZprimeVertices();   // <-- new: initialises model vertices

  public :

    Standard_Model_Zprime();
    bool ModelInit();
    void InitVertices();

  };

}

#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "MODEL/Main/Single_Vertex.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(Standard_Model_Zprime,"SMZprime",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,Standard_Model_Zprime>::
operator()(const Model_Arguments &args) const
{
  return new Standard_Model_Zprime();
}

void Getter<Model_Base,Model_Arguments,Standard_Model_Zprime>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"The Standard Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- EW_SCHEME (EW input scheme, see documentation)\n"
     <<setw(width+7)<<" "<<"- EW_REN_SCHEME (EW renormalisation scheme, see documentation)\n"
     <<setw(width+7)<<" "<<"- WIDTH_SCHEME (Fixed or CMS, see documentation)\n"
     <<setw(width+7)<<" "<<"- ALPHAS(MZ) (strong coupling at MZ)\n"
     <<setw(width+7)<<" "<<"- ORDER_ALPHAS (0,1,2 -> 1, 2, 3-loop running)\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- ALPHAQED_DEFAULT_SCALE (scale for alpha_QED default)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+7)<<" "<<"- VEV (Higgs vev)\n"
     <<setw(width+7)<<" "<<"- CKM_ORDER (0,1,2,3 - order of CKM expansion in Cabibbo angle)\n"
     <<setw(width+7)<<" "<<"- CKM_CABIBBO (Cabibbo angle in Wolfenstein parameterization)\n"
     <<setw(width+7)<<" "<<"- CKM_A (Wolfenstein A)\n"
     <<setw(width+7)<<" "<<"- CKM_RHO (Wolfenstein Rho)\n"
     <<setw(width+7)<<" "<<"- CKM_ETA (Wolfenstein Eta)\n"
     <<setw(width+7)<<" "<<"- CKM_ELEMENT[<i>][<j>] (explicit value for element, supersedes parametrisation)\n"
     <<setw(width+7)<<" "<<"- Zprime mass/width via MASS[32] & WIDTH[32]\n"
     <<setw(width+7)<<" "<<"- multiplicative coupling parameter Zp_cpl_L\n"
     <<setw(width+7)<<" "<<"- multiplicative coupling parameter Zp_cpl_R\n"
     <<setw(width+4)<<" "<<"}";
  str<<"Infrared continuation of alphaS:\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"- AS_FORM (values 0,1,2,3,10, see documentation)\n"
     <<setw(width+7)<<" "<<"- Q2_AS (corresponding infrared parameter, see documentation)\n"
     <<setw(width+4)<<" "<<"}";
}

Standard_Model_Zprime::Standard_Model_Zprime() :
  Standard_Model()
{
  m_name="SM+Zprime";
  ParticleZprimeInit();
  RegisterDefaults();
  ReadParticleData();
}

bool Standard_Model_Zprime::ModelInit()
{
  FixZprimeParameters();
  Standard_Model::ModelInit();
}


void Standard_Model_Zprime::InitVertices()
{
  Standard_Model::InitVertices();
  InitZprimeVertices();
}

void Standard_Model_Zprime::FixZprimeParameters() {
  auto s = Settings::GetMainSettings()["Zprime"];
  s["Zp_cpl_L"].SetDefault(1);
  s["Zp_cpl_R"].SetDefault(1);
  p_constants->insert(make_pair("Zp_cpl_L",s["Zp_cpl_L"].Get<double>()));
  p_constants->insert(make_pair("Zp_cpl_R",s["Zp_cpl_R"].Get<double>()));
}

void Standard_Model_Zprime::ParticleZprimeInit()
{
  // add Zprime
  // kf_code,mass,radius,width,3*charge,strong,spin,majorana,take,stable,massive,idname,antiname,texname,antitexname
  s_kftable[kf_Zp] = new Particle_Info(kf_Zp,1000.,0,10.,0,0,2,-1,1,0,1,"Zprime","Zprime","Z^{\\prime}","Z^{\\prime}");
}


void Standard_Model_Zprime::InitZprimeVertices()
{
  // set up constants for the model
  Kabbala I("i",Complex(0.,1.)), rt2("\\sqrt(2)",sqrt(2.0));
  Kabbala sintW("\\sin\\theta_W",sqrt(ComplexConstant("csin2_thetaW")));
  Kabbala costW("\\cos\\theta_W",sqrt(ComplexConstant("ccos2_thetaW")));
  // coupling constants
  Kabbala g1("g_1",sqrt(4.*M_PI*ScalarConstant("alpha_QED")));
  Kabbala g2("g_1/\\cos\\theta_W", g1.Value()/costW.Value());

<<<<<<< HEAD
  Kabbala PL("P_L",ScalarConstant("Zp_cpl_L"));
  Kabbala PR("P_R",ScalarConstant("Zp_cpl_R"));

=======
>>>>>>> bcdf56de60e1167a5c903794ce67615ad83b97c7
  // the parameter specifying the LR model
  // - sqrt(2.) will describe a totally LR-symm model
  // - sqrt(2./3.) describes an E6-inspired model
  Kabbala alphaLR("\\alpha_{LR}",sqrt(2./3.));

  // create FFV vertices with Z' if it's on
  if (Flavour(kf_Zp).IsOn()) {
    for (short int i=1;i<17;++i) {
      // parse through all fermions that couple to Z' and create vertices
      if (i==7) i=11;
      Flavour flav((kf_code)i);
      if (!flav.IsOn()) continue;
      Kabbala B = Kabbala(string("B_{")+flav.TexName()+string("}"),
                          flav.IsQuark()?(flav.IsAnti()?-1./3.:1./3.):0.);
      Kabbala L = Kabbala(string("L_{")+ flav.TexName()+string("}"),
                          flav.IsLepton()?(flav.IsAnti()?-1:1):0.);
      Kabbala Y3R = Kabbala(string("YR_{")+flav.TexName()+string("}"),
                            flav.IsoWeak());

      // create the vertex for that particular fermion and a Z'.
      // Right-handed neutrinos will not take part in any interaction.
      Kabbala kcpl0;
      if (flav.Kfcode()==kf_nue || flav.Kfcode()==kf_numu ||
          flav.Kfcode()==kf_nutau)
        kcpl0 = Kabbala("0.0", 0.);
      else
        kcpl0 = -I * g2 * (Y3R * alphaLR + (L-B)/(alphaLR*2));
      Kabbala kcpl1 = -I * g2 * (L-B) / (alphaLR*2);
      m_v.push_back(Single_Vertex());
      m_v.back().AddParticle(flav.Bar());
      m_v.back().AddParticle(flav);
      m_v.back().AddParticle(Flavour(kf_Zp));
      m_v.back().Color.push_back
        (i>10?Color_Function(cf::None):
         Color_Function(cf::D,1,2));
      m_v.back().Color.push_back
        (i>10?Color_Function(cf::None):
         Color_Function(cf::D,1,2));
      m_v.back().Lorentz.push_back("FFVL");
      m_v.back().Lorentz.push_back("FFVR");
      m_v.back().cpl.push_back(kcpl0);
      m_v.back().cpl.push_back(kcpl1);
      m_v.back().order[1]=1;
    }
  }
}
