#include "Six_Particle_Observables.H"
#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                           \
  NAME::operator()(const String_Matrix &parameters) const               \
{ return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const      \
{ str<<"kf1 kf2 kf3 kf4 kf5 kf6 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Six_Particle_Observable_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Six_Particle_Observable_Base(*old) { }                          \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

Six_Particle_Observable_Base::Six_Particle_Observable_Base(const String_Matrix & parameters,
                                                           const std::string & obsname):
    Primitive_Observable_Base(parameters), m_flavs(6)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<10) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    for (short unsigned int i=0;i<6;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      m_flavs[i]=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
      if (kf<0) m_flavs[i]=m_flavs[i].Bar();
    }

    m_xmin  = ATOOLS::ToType<double>(parameters[0][6]);
    m_xmax  = ATOOLS::ToType<double>(parameters[0][7]);
    int nbins = ATOOLS::ToType<int>(parameters[0][8]);
    m_type  = HistogramType(parameters[0][9]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    MyStrStream str;
    str<<obsname<<m_flavs[0]<<m_flavs[1]<<m_flavs[2]<<m_flavs[3]<<m_flavs[4]<<m_flavs[5]<<".dat";
    str>>m_name;

    m_listname = parameters[0].size()>10?parameters[0][10]:"Analysed";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      for (short unsigned int j=0;j<6;++j) {
        if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
          int kf=ATOOLS::ToType<int>(parameters[i][1]);
          m_flavs[j] = ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
        }
      }
    }

    if (m_name=="SherpaDefault") {
      MyStrStream str;
      str<<obsname<<m_flavs[0]<<m_flavs[1]<<m_flavs[2]<<m_flavs[3]<<m_flavs[4]<<m_flavs[5];
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "";
      MyStrStream str;
      str<<obsname<<" of "<<m_flavs[0]<<", "<<m_flavs[1]<<", "<<m_flavs[2]<<", "
         <<m_flavs[3]<<", "<<m_flavs[4]<<", "<<m_flavs[5];
      str>>title;
      p_histo->SetTitle(title);
    }
  }
}

Six_Particle_Observable_Base::Six_Particle_Observable_Base(const Six_Particle_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_flavs=old->m_flavs;
}

void Six_Particle_Observable_Base::Evaluate(double value, double weight,
					     int ncount) {
  p_histo->Insert(value,weight,ncount); 
}

 
void Six_Particle_Observable_Base::Evaluate(int nout, const Vec4D* moms,
					     const Flavour* flavs,
					     double weight, int ncount) 
{
  for (int i=0;i<nout;i++) { 
    if (flavs[i]==m_flavs[0]) {
      for (int j=0;j<nout;j++) { 
        if (flavs[j]==m_flavs[1] && i!=j) {
          for (int k=0;k<nout;k++) { 
            if (flavs[k]==m_flavs[2] && k!=j && k!=i) {
              for (int l=0;l<nout;l++) { 
                if (flavs[l]==m_flavs[3] && l!=k && l!=j && l!=i) {
                  for (int m=0;m<nout;m++) { 
                    if (flavs[m]==m_flavs[4] && m!=l && m!=k && m!=j && m!=i) {
                      for (int n=0;n<nout;n++) { 
                        if (flavs[n]==m_flavs[5] && n!=m && n!=l && n!=k && n!=j && n!=i) {
                          Evaluate(moms[i],moms[j],moms[k],moms[l],moms[m],moms[n],weight,ncount);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        } 
      }
    }
  }
}


void Six_Particle_Observable_Base::Evaluate(const Particle_List& plist,
					     double weight, int ncount) {
  for(Particle_List::const_iterator plit1=plist.begin();
      plit1!=plist.end(); ++plit1) {
    if((*plit1)->Flav()==m_flavs[0]) {
      for(Particle_List::const_iterator plit2=plist.begin();
          plit2!=plist.end(); ++plit2) {
        if((*plit2)->Flav()==m_flavs[1] && plit1!=plit2) {
          for(Particle_List::const_iterator plit3=plist.begin();
              plit3!=plist.end(); ++plit3) {
            if((*plit3)->Flav()==m_flavs[2] && plit3!=plit2 && plit3!=plit1) {
              for(Particle_List::const_iterator plit4=plist.begin();
                  plit4!=plist.end(); ++plit4) {
                if((*plit4)->Flav()==m_flavs[3] &&
                    plit4!=plit3 && plit4!=plit2 && plit4!=plit1) {
                  for(Particle_List::const_iterator plit5=plist.begin();
                      plit5!=plist.end(); ++plit5) {
                    if((*plit5)->Flav()==m_flavs[4] &&
                        plit5 != plit4 && plit5!=plit3 && plit5!=plit2 && plit5!=plit1) {
                      for(Particle_List::const_iterator plit6=plist.begin();
                          plit6!=plist.end(); ++plit6) {
                        if((*plit6)->Flav()==m_flavs[5] &&
                            plit6 != plit5 && plit6 != plit4 && plit6!=plit3 && plit6!=plit2 && plit6!=plit1) {
                          Evaluate(
                              (*plit1)->Momentum(),(*plit2)->Momentum(),
                              (*plit3)->Momentum(),(*plit4)->Momentum(),
                              (*plit5)->Momentum(),(*plit6)->Momentum(),
                              weight, ncount);
                          return;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  p_histo->Insert(0.0,0.0,ncount);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Six_Particle_PlaneAngle,
			 Six_Particle_PlaneAngle_Getter,"PlaneAngle6")

void Six_Particle_PlaneAngle::Evaluate(
    const Vec4D & mom1,
    const Vec4D & mom2,
    const Vec4D & mom3,
    const Vec4D & mom4,
    const Vec4D & mom5,
    const Vec4D & mom6,
    double weight, int ncount)
{
  Vec4D inter1  = mom1+mom2;
  Vec4D inter2  = mom4+mom5;
  Vec3D normal1 = cross(Vec3D(inter1),Vec3D(mom3));
  Vec3D normal2 = cross(Vec3D(inter2),Vec3D(mom6));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs()); 
  p_histo->Insert(costh,weight,ncount); 
}
