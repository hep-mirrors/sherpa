#include "Four_Jet_Angles.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                   \
  NAME::operator()(const String_Matrix &parameters) const               \
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(ostream &str,const size_t width) const      \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Four_Jet_Angle_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Four_Jet_Angle_Base(*old) { }                            \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

Four_Jet_Angle_Base::Four_Jet_Angle_Base(const String_Matrix & parameters,
                                                           const string & obsname):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<4) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
                 <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_listname = parameters[0].size()>4?parameters[0][4]:"Analysed";
    m_name="";
    if (m_listname!="Analysed") m_name=m_listname+string("_");
    m_name+="4jet_"+obsname;
    m_name+=".dat";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    if (m_name=="SherpaDefault") {
      m_name="";
      if (m_listname!="Analysed") m_name=m_listname+string("_");
      m_name+="4jet_"+obsname;
    }
    if (p_histo->Title()=="SherpaDefault") {
      string title = "4jet "+obsname;
      p_histo->SetTitle(title);
    }
  }
}

Four_Jet_Angle_Base::Four_Jet_Angle_Base(const Four_Jet_Angle_Base * old) :
    Primitive_Observable_Base(*old) {}

void Four_Jet_Angle_Base::Evaluate(const Blob_List & blobs,double value, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);

  // sort in Durham (and NOT in Final_Selector!)

  vector<Vec3D> moms;
 
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    moms.push_back((*pit)->Momentum());
  }

  if (moms.size()!=4) {
    p_histo->Insert(0.,0.,ncount);
    return;
  }

//   for (size_t i=0;i<4;++i)
//     cout<<" p["<<i<<"]="<<moms[i]<<endl;

  double cos_chi=Calc(moms);
  if (p_histo->Xmin()==0.) cos_chi=dabs(cos_chi);
//   cout<<Name()<<" cos_chi="<<cos_chi<<" ("<<p_histo->Xmin()<<")"<<endl;
  p_histo->Insert(cos_chi,value,ncount);
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Bengtsson_Zerwas_Angle,
			 Bengtsson_Zerwas_Angle_Getter,"BZAngle")

double Bengtsson_Zerwas_Angle::Calc(const vector<Vec3D> & moms)
{
  Vec3D p12=cross(moms[0],moms[1]);
  Vec3D p34=cross(moms[2],moms[3]);

  return (p12*p34)/(p12.Abs()*p34.Abs());
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Nachtmann_Reiter_Angle,
			 Nachtmann_Reiter_Angle_Getter,"NRAngle")

double Nachtmann_Reiter_Angle::Calc(const vector<Vec3D> & moms)
{
  Vec3D p1_2=moms[0]-moms[1];
  Vec3D p3_4=moms[2]-moms[3];

  return (p1_2*p3_4)/(p1_2.Abs()*p3_4.Abs());
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Koerner_Schierholz_Willrodt_Angle,
			 Koerner_Schierholz_Willrodt_Angle_Getter,"KSWAngle")

double Koerner_Schierholz_Willrodt_Angle::Calc(const vector<Vec3D> & moms)
{
  Vec3D p14=cross(moms[0],moms[3]);
  Vec3D p23=cross(moms[1],moms[2]);
  double c1423= (p14*p23)/(p14.Abs()*p23.Abs());
  Vec3D p13=cross(moms[0],moms[2]);
  Vec3D p24=cross(moms[1],moms[3]);
  double c1324= (p13*p24)/(p13.Abs()*p24.Abs());
  //  cout<<" phi1="<<acos(c1423)<<"  phi2="<<acos(c1324)<<endl;
  
  double phi_ksw=0.5*(acos(c1423)+acos(c1324));
  return cos(phi_ksw);
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Alpha34_Angle,Alpha34_Angle_Getter,"A34Angle")

double Alpha34_Angle::Calc(const vector<Vec3D> & moms)
{
  return (moms[2]*moms[3])/(moms[2].Abs()*moms[3].Abs());
}
