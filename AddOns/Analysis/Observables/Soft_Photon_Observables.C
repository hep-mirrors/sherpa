#include "AddOns/Analysis/Observables/Soft_Photon_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7 || parameters[0].size()>15) return NULL;
    std::vector<ATOOLS::Flavour> f(parameters[0].size()-5);
    for (size_t i=0;i<f.size();++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0][parameters[0].size()-1];
    return new Class(f,HistogramType(parameters[0][parameters[0].size()-2]),
         ATOOLS::ToType<double>(parameters[0][parameters[0].size()-5]),
         ATOOLS::ToType<double>(parameters[0][parameters[0].size()-4]),
         ATOOLS::ToType<int>(parameters[0][parameters[0].size()-3]),list);
  }
  else return NULL;
//   else if (parameters.size()<10) return NULL;
//   double min=0.0, max=1.0;
//   size_t bins=100;
//   std::vector<ATOOLS::Flavour> f(10);
//   std::string list=finalstate_list, scale="Lin";
//   for (size_t i=0;i<parameters.size();++i) {
//     if (parameters[i].size()<2) continue;
//     for (short unsigned int j=0;j<6;++j) {
//       if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
//         int kf=ATOOLS::ToType<int>(parameters[i][1]);
//         f[j]=ATOOLS::Flavour((kf_code)abs(kf));
//         if (kf<0) f[j]=f[j].Bar();
//       }
//     }
//     if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
//     else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
//     else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
//     else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
//     else if (parameters[i][0]=="LIST") list=parameters[i][1];
//   }
//   return new Class(f,HistogramType(scale),min,max,bins,list);
}

#define DEFINE_GETTER_METHOD(CLASS,NAME)        \
  Primitive_Observable_Base *         \
  NAME::operator()(const Argument_Matrix &parameters) const   \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)         \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const  \
  { str<<"kf1 ... kfn min max bins Lin|LinErr|Log|LogErr list , 1<n<11"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)      \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix); \
  DEFINE_GETTER_METHOD(CLASS,NAME)          \
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

Soft_Photon_Observable_Base::Soft_Photon_Observable_Base
(const std::vector<Flavour>& flavs, int type, double xmin, double xmax,
 int nbins, const std::string& listname, const std::string& name)
  : Primitive_Observable_Base(type,xmin,xmax,nbins), f_special(false)
{
  if(flavs.size()<2) {
    msg_Error()<<"Error in Soft_Photon_Observable_Base:"<<std::endl
         <<"   Less than two flavours specified, system undefined."
         <<std::endl;
    msg_Error()<<"number of flavours is: "<<flavs.size()<<std::endl;
  }
  m_name = name + "_";
  for (size_t i=0;i<flavs.size();++i) {
    m_name+=flavs[i].ShellName();
    m_flavs.push_back(flavs[i]);
  }
  m_name += ".dat";
  m_listname = listname;
  m_blobtype = std::string("YFS-type QED Corrections to ME");
  m_blobdisc = true;
  if(xmin>=0.0) f_special=true;

}

void Soft_Photon_Observable_Base::Evaluate(double value, double weight,
                                           double ncount)
{
  p_histo->Insert(value,weight,ncount); 
}

void Soft_Photon_Observable_Base::Evaluate(int nout, const Vec4D* moms,
                                           const Flavour* flavs,
                                           double weight, double ncount)
{
  msg_Out()<<"Flavour array method"<<endl;
  msg_Out()<<"I don't do anything"<<endl;
}

void Soft_Photon_Observable_Base::Evaluate(const ATOOLS::Particle_List& plist,
                                           double weight, double ncount)
{
  msg_Out()<<"Particle list method"<<endl;
  msg_Out()<<"I don't do anything"<<endl;
}

void Soft_Photon_Observable_Base::Evaluate(const ATOOLS::Blob_List& blobs,
                                           double weight, double ncount)
{
  Blob * QEDblob = NULL;
  for (size_t i=0;i<blobs.size();++i) {
    if (blobs[i]->Type()==btp::Shower &&
        blobs[i]->TypeSpec()==m_blobtype) {
      QEDblob=blobs[i];
      break;
    }
  }
  if (!QEDblob) {
    Evaluate(0.,weight,ncount);
    return;
  }
  Particle_Vector parts;
  for (int i=0;i<QEDblob->NOutP();++i) {
    for (size_t j=0;j<m_flavs.size();++j) {
      if (QEDblob->OutParticle(i)->Flav()==m_flavs[j] ||
          QEDblob->OutParticle(i)->Info()=='S') {
        parts.push_back(QEDblob->OutParticle(i));
        j=m_flavs.size();
      }
    }
  }
  Evaluate(parts,weight,ncount);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DEFINE_OBSERVABLE_GETTER(Soft_Photon_Energy, Soft_Photon_Energy_Getter,
                         "SoftPhotonEnergy")

void Soft_Photon_Energy::Evaluate(const ATOOLS::Particle_Vector& parts,
                                  double weight, double ncount)
{
  Vec4D decayer(0.,0.,0.,0.);
  for (size_t i=0;i<parts.size();++i) {
    decayer+=parts[i]->Momentum();
  }
  Poincare decframe(decayer);
  double Ephots(0.);
  for (size_t i=0;i<parts.size();++i) {
    if (parts[i]->Flav().Kfcode()==kf_photon) {
      Vec4D k=parts[i]->Momentum();
      decframe.Boost(k);
      Ephots+=k[0];
    }
  }
  p_histo->Insert(Ephots,weight,ncount);
}

Soft_Photon_Energy::Soft_Photon_Energy(const std::vector<Flavour>& flavs,
                                       int type,double xmin,double xmax,
                                       int nbins,
                                       const std::string & listname)
  : Soft_Photon_Observable_Base(flavs,type,xmin,xmax,nbins,listname,
                                "SoftPhotonEnergy") {}

Primitive_Observable_Base* Soft_Photon_Energy::Copy() const
{
  return new Soft_Photon_Energy(m_flavs,m_type,m_xmin,m_xmax,m_nbins,
                                m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Soft_Photon_Angle,Soft_Photon_Angle_Getter,
                         "SoftPhotonAngle")

void Soft_Photon_Angle::Evaluate(const ATOOLS::Particle_Vector& parts,
                                 double weight, double ncount)
{
  Vec4D multipole(0.,0.,0.,0.);
  Vec4D axpart(0.,0.,0.,1.);
  Vec4D axis(0.,0.,0.,1.);
  int chargesum(0);
  // FS charged momentum sum
  for (size_t i=0;i<parts.size();++i) {
    if (parts[i]->Flav().IntCharge()!=0) {
      multipole+=axpart=parts[i]->Momentum();
      chargesum+=parts[i]->Flav();
    }
  }
  // neutral resonance -> nothing further
  // charged resonance -> add to multipole, change definition for theta=0
  if (chargesum!=0) {
    // add charged resonance to multipole
    Vec4D resonance(0.,0.,0.,0.);
    for (size_t i=0;i<parts.size();++i) resonance+=parts[i]->Momentum();
    multipole+=resonance;
    // resonance at theta=0
    axpart=resonance;
  }
  Poincare multipoleboost(multipole);
  multipoleboost.Boost(axpart);
  Poincare rotation(axpart,axis);
  for (size_t i=0;i<parts.size();++i) {
    if (parts[i]->Flav().IsPhoton()) {
      Vec4D k=parts[i]->Momentum();
      multipoleboost.Boost(k);
      rotation.Rotate(k);
      double theta = acos((Vec3D(k)*Vec3D(axis))/(Vec3D(k).Abs()));
      p_histo->Insert(theta,weight,ncount);
    }
  }
}

Soft_Photon_Angle::Soft_Photon_Angle(const std::vector<Flavour>& flavs,
                                     int type, double xmin, double xmax,
                                     int nbins,
                                     const std::string & listname)
  : Soft_Photon_Observable_Base(flavs,type,xmin,xmax,nbins,listname,
                                "SoftPhotonAngle") {}

Primitive_Observable_Base* Soft_Photon_Angle::Copy() const
{
  return new Soft_Photon_Angle(m_flavs,m_type,m_xmin,m_xmax,m_nbins,
                               m_listname);
}

