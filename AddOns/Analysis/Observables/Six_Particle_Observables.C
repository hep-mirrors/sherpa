#include "AddOns/Analysis/Observables/Six_Particle_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<10) return NULL;
    std::vector<ATOOLS::Flavour> f(6);
    for (short unsigned int i=0;i<6;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0].size()>10?parameters[0][10]:finalstate_list;
    return new Class(f,HistogramType(parameters[0][9]),
         ATOOLS::ToType<double>(parameters[0][6]),
         ATOOLS::ToType<double>(parameters[0][7]),
         ATOOLS::ToType<int>(parameters[0][8]),list);
  }
  else if (parameters.size()<10) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::vector<ATOOLS::Flavour> f(6);
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    for (short unsigned int j=0;j<6;++j) {
      if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
        int kf=ATOOLS::ToType<int>(parameters[i][1]);
        f[j]=ATOOLS::Flavour((kf_code)abs(kf));
        if (kf<0) f[j]=f[j].Bar();
      }
    }
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(f,HistogramType(scale),min,max,bins,list);
}                 

#define DEFINE_GETTER_METHOD(CLASS,NAME)        \
  Primitive_Observable_Base *         \
  NAME::operator()(const Argument_Matrix &parameters) const   \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)         \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const  \
  { str<<"kf1 kf2 kf3 kf4 kf5 kf6 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)      \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix); \
  DEFINE_GETTER_METHOD(CLASS,NAME)          \
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

Six_Particle_Observable_Base::Six_Particle_Observable_Base
(const std::vector<Flavour>& flavs, int type, double xmin, double xmax,
 int nbins, const std::string& listname, const std::string& name)
  : Primitive_Observable_Base(type,xmin,xmax,nbins), f_special(false) {

  if(flavs.size()<6) {
    msg_Error()<<"Error in Six_Particle_Observable_Base:"<<std::endl
         <<"   No six flavours specified, try to copy flavours."
         <<std::endl;
    msg_Error()<<"number of flavours is: "<<flavs.size()<<std::endl;
  }
  MyStrStream str;
  str<<name<<flavs[0]<<flavs[1]<<flavs[2]<<flavs[3]<<flavs[4]<<flavs[5]<<".dat";
  str>>m_name;
  Flavour fl;
  for(size_t i=0; i<6; i++) {
    if(i<flavs.size()) fl=flavs[i];
    m_flavs.push_back(fl);
  }
  m_listname = listname;
  m_blobtype = std::string("");
  m_blobdisc = false;
  if(xmin>=0.0) f_special=true;

}

void Six_Particle_Observable_Base::Evaluate(double value, double weight,
               double ncount) {
  p_histo->Insert(value,weight,ncount); 
}

 
void Six_Particle_Observable_Base::Evaluate(int nout, const Vec4D* moms,
               const Flavour* flavs,
               double weight, double ncount) 
{
  for (int i=0;i<nout;i++) { 
    if (flavs[i]==m_flavs[0]) {
      for (int j=0;j<nout;j++) { 
        if (flavs[j]==m_flavs[0] && i!=j) {
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
               double weight, double ncount) {
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
                if((*plit4)->Flav()==m_flavs[3] && plit4!=plit3 && plit4!=plit2
                                                && plit4!=plit1) {
                  for(Particle_List::const_iterator plit5=plist.begin();
                      plit5!=plist.end(); ++plit5) {
                    if((*plit5)->Flav()==m_flavs[4] && plit5!=plit4 &&
                       plit5!=plit3 && plit5!=plit2 && plit5!=plit1) {
                      for(Particle_List::const_iterator plit6=plist.begin();
                          plit6!=plist.end(); ++plit6) {
                        if((*plit6)->Flav()==m_flavs[5] && plit6!=plit5 &&
                           plit6!=plit4 && plit6!=plit3 && plit6!=plit2 &&
                           plit6!=plit1) {
                          Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),
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


DEFINE_OBSERVABLE_GETTER(Six_Particle_PT,
       Six_Particle_PT_Getter,"PT6")

void Six_Particle_PT::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                               const Vec4D& mom3,const Vec4D& mom4,
                               const Vec4D& mom5,const Vec4D& mom6,
                               double weight, double ncount)
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]+mom3[1]+mom4[1]+mom5[1]+mom6[1]) +
                   sqr(mom1[2]+mom2[2]+mom3[2]+mom4[2]+mom5[2]+mom6[2]));
  p_histo->Insert(pt,weight,ncount);
}

Six_Particle_PT::Six_Particle_PT(const std::vector<Flavour>& flavs,
           int type,double xmin,double xmax,int nbins,
           const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"PT") {}

Primitive_Observable_Base* Six_Particle_PT::Copy() const
{
  return new Six_Particle_PT(m_flavs,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Six_Particle_ET,
       Six_Particle_ET_Getter,"ET6")

void Six_Particle_ET::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                               const Vec4D& mom3,const Vec4D& mom4,
                               const Vec4D& mom5,const Vec4D& mom6,
                               double weight, double ncount)
{
  double pt2 = sqr(mom1[1]+mom2[1]+mom3[1]+mom4[1]+mom5[1]+mom6[1]) +
               sqr(mom1[2]+mom2[2]+mom3[2]+mom4[2]+mom5[2]+mom6[2]);
  double p2  = sqr(mom1[3]+mom2[3]+mom3[3]+mom4[3]+mom5[3]+mom6[3])+pt2;
  double et  = (mom1[0]+mom2[0]+mom3[0]+mom4[0]+mom5[0]+mom6[0])*sqrt(pt2/p2);
  p_histo->Insert(et,weight,ncount);
}

Six_Particle_ET::Six_Particle_ET(const std::vector<Flavour>& flavs,
           int type,double xmin,double xmax,int nbins,
           const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"ET") {}

Primitive_Observable_Base* Six_Particle_ET::Copy() const
{
  return new Six_Particle_ET(m_flavs,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Six_Particle_Mass,
       Six_Particle_Mass_Getter,"6Mass")

void Six_Particle_Mass::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                                 const Vec4D& mom3,const Vec4D& mom4,
                                 const Vec4D& mom5,const Vec4D& mom6,
                                 double weight, double ncount)
{
  Vec4D  p = mom1+mom2+mom3+mom4+mom5+mom6;
  p_histo->Insert(p.Mass(),weight,ncount);
}

Six_Particle_Mass::Six_Particle_Mass(const std::vector<Flavour>& flavs,
                            int type,double xmin,double xmax,int nbins,
                            const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"6Mass") {}

Primitive_Observable_Base* Six_Particle_Mass::Copy() const
{
  return new Six_Particle_Mass(m_flavs,m_type,m_xmin,
                               m_xmax,m_nbins,m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partontriplett_DeltaPhi,
       Two_Partontriplett_DeltaPhi_Getter,"DPhi3_2")

void Two_Partontriplett_DeltaPhi::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                                           const Vec4D& mom3,const Vec4D& mom4,
                                           const Vec4D& mom5,const Vec4D& mom6,
                                           double weight, double ncount) {
  Vec4D vecA(mom1); vecA+=mom2; vecA+=mom3;
  Vec4D vecB(mom4); vecB+=mom5; vecB+=mom6;
  double dphi(vecA.Phi()-vecB.Phi());
  if (dphi >= 0.) {
    if (dphi < M_PI)
      p_histo->Insert(dphi, weight, ncount);
    else
      p_histo->Insert(2.*M_PI-dphi, weight, ncount);
  }
  else {
    if (dphi > -M_PI)
      p_histo->Insert(-dphi, weight, ncount);
    else
      p_histo->Insert(2.*M_PI+dphi, weight, ncount);
  }
}

Two_Partontriplett_DeltaPhi::Two_Partontriplett_DeltaPhi(
             const std::vector<Flavour>& flavs,
             int type, double xmin, double xmax,
             int nbins,
             const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,
                                 nbins,listname,"DPhi3_2") {}

Primitive_Observable_Base* Two_Partontriplett_DeltaPhi::Copy() const {
  return new Two_Partontriplett_DeltaPhi(m_flavs, m_type, m_xmin, m_xmax,
                                         m_nbins, m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partontriplett_DeltaEta,
       Two_Partontriplett_DeltaEta_Getter,"DEta3_2")

void Two_Partontriplett_DeltaEta::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                                           const Vec4D& mom3,const Vec4D& mom4,
                                           const Vec4D& mom5,const Vec4D& mom6,
                                           double weight, double ncount) {
  Vec4D vecA(mom1); vecA+=mom2; vecA+=mom3;
  Vec4D vecB(mom4); vecB+=mom5; vecB+=mom6;
  double deta(abs(vecA.Eta()-vecB.Eta()));
  p_histo->Insert(deta, weight, ncount);
}

Two_Partontriplett_DeltaEta::Two_Partontriplett_DeltaEta(
             const std::vector<Flavour>& flavs,
             int type, double xmin, double xmax,
             int nbins,
             const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,
                                 nbins,listname,"DEta3_2") {}

Primitive_Observable_Base* Two_Partontriplett_DeltaEta::Copy() const {
  return new Two_Partontriplett_DeltaEta(m_flavs, m_type, m_xmin, m_xmax,
                                         m_nbins, m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partontriplett_DR,
                         Two_Partontriplett_DR_Getter,"DR3_2")

void Two_Partontriplett_DR::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
                                     const Vec4D& mom3,const Vec4D& mom4,
                                     const Vec4D& mom5,const Vec4D& mom6,
                                     double weight, double ncount) {
  Vec4D vecA(mom1); vecA+=mom2; vecA+=mom3;
  Vec4D vecB(mom4); vecB+=mom5; vecB+=mom6;
  double pt1=sqrt(vecA[1]*vecA[1]+vecA[2]*vecA[2]);
  double pt2=sqrt(vecB[1]*vecB[1]+vecB[2]*vecB[2]);
  double dphi=acos((vecA[1]*vecB[1]+vecA[2]*vecB[2])/(pt1*pt2));
  double c1=vecA[3]/Vec3D(vecA).Abs();
  double c2=vecB[3]/Vec3D(vecB).Abs();
  double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
  double dR= sqrt(sqr(deta) + sqr(dphi));
  p_histo->Insert(dR, weight, ncount);
}

Two_Partontriplett_DR::Two_Partontriplett_DR(
             const std::vector<Flavour>& flavs,
             int type, double xmin, double xmax,
             int nbins,
             const std::string & listname)
  : Six_Particle_Observable_Base(flavs,type,xmin,xmax,
                                 nbins,listname,"DR3_2") {}

Primitive_Observable_Base* Two_Partontriplett_DR::Copy() const {
  return new Two_Partontriplett_DR(m_flavs, m_type, m_xmin, m_xmax,
                                   m_nbins, m_listname);
}

