#include "AddOns/Analysis/Observables/Jet_Mass_and_Broadening.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Observables/Event_Shapes_EE.H"

namespace ANALYSIS {
  std::ostream& operator<<( std::ostream& ostr, const JetMass_Broadening_Data & data) {
    ostr<<"JetMass_Broadening_Data : "<<data.heavyjetmass<<","<<data.lightjetmass
	<<","<<data.widejetbroadening<<","<<data.narrowjetbroadening;
    return ostr;
  }
}

using namespace ANALYSIS;
using namespace ATOOLS;

#include "ATOOLS/Org/MyStrStream.H"

DECLARE_GETTER(JetMass_Broadening_Calculator_Getter,"JBCalc",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
JetMass_Broadening_Calculator_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string listname=finalstate_list;
  if (parameters.size()>0 && parameters[0].size()>0) listname=parameters[0][0];
  return new JetMass_Broadening_Calculator(listname);
}

void JetMass_Broadening_Calculator_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"[list] -> EEShapes"; 
}

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list] -> JBCalc"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)


JetMass_Broadening_Calculator::JetMass_Broadening_Calculator(const std::string & listname)
  :  m_inkey("EvtShapeData"), m_outkey(listname+"_JetMass_Broadening") 
{
  m_name = listname+"_JetMass_Broadening_Calculator";
  m_listname = listname;
}

void JetMass_Broadening_Calculator::Evaluate(const Blob_List & ,double weight, double ncount) {
  Particle_List * pl = p_ana->GetParticleList(m_listname);
  if (pl==NULL) {
    msg_Out()<<"WARNING in JetMass_Broadening_Calculator::Evaluate : particle list "
	     <<m_listname<<" not found "<<std::endl;
    return;
  }
  Vec3D thrustaxis;
  Blob_Data_Base * data = (*p_ana)[m_inkey];
  if (data) {
    thrustaxis = data->Get<Event_Shape_EE_Data>().thrustaxis;
  }
  else {
    msg_Out()<<"WARNING in JetMass_Broadening_Calculator::Evaluate : data container "
	     <<m_inkey<<" not found "<<std::endl;
    return;
  }

  double mh=0., ml=0., bw=0., bn=0.;
  Vec4D totalsum;
  Vec4D sum1;
  Vec4D sum2;
  double totalb=0., b1=0., b2=0.;
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    Vec4D mom = (*pit)->Momentum();
    Vec3D v   = Vec3D(mom);
    totalsum+=mom;
    totalb += v.Abs();
    if (thrustaxis*v>0.) {
      sum1 += mom;
      b1   += cross(v,thrustaxis).Abs();
    }
    else {
      sum2 += mom;
      b2   += cross(v,thrustaxis).Abs();
    }
  }
  if (pl->size()>0) {
    double evis2 = totalsum.Abs2();
    double m1    = sum1.Abs2()/evis2;
    double m2    = sum2.Abs2()/evis2;
    totalb*=2;
    b1/=totalb;
    b2/=totalb;
    if (m1<m2) {
      mh=m2; ml=m1;
    }
    else {
      mh=m1; ml=m2;
    }
    if (b1<b2) {
      bw=b2; bn=b1;
    }
    else {
      bw=b1; bn=b2;
    }
  }
  
  p_ana->AddData(m_outkey,
		 new Blob_Data<JetMass_Broadening_Data>(JetMass_Broadening_Data(mh,ml,bw,bn)));
}

Analysis_Object * JetMass_Broadening_Calculator::GetCopy() const
{
  return new JetMass_Broadening_Calculator(m_listname);
}



// ----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Heavy_Jet_Mass,Heavy_Jet_Mass_Getter,"MH")


Heavy_Jet_Mass::Heavy_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Heavy_Jet_Mass.dat");
}

void Heavy_Jet_Mass::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().heavyjetmass,weight,ncount);
  }
}
 
Primitive_Observable_Base * Heavy_Jet_Mass::Copy() const
{
  return new Heavy_Jet_Mass(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Light_Jet_Mass,Light_Jet_Mass_Getter,"ML")


Light_Jet_Mass::Light_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Light_Jet_Mass.dat");
}

void Light_Jet_Mass::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().lightjetmass,weight,ncount);
  }
}
 
Primitive_Observable_Base * Light_Jet_Mass::Copy() const
{
  return new Light_Jet_Mass(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Jet_Mass_Difference,
			 Jet_Mass_Difference_Getter,"MD")


Jet_Mass_Difference::Jet_Mass_Difference(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Jet_Mass_Difference.dat");
}

void Jet_Mass_Difference::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().heavyjetmass
		    -data->Get<JetMass_Broadening_Data>().lightjetmass,
		    weight,ncount);
  }
}
 
Primitive_Observable_Base * Jet_Mass_Difference::Copy() const
{
  return new Jet_Mass_Difference(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Single_Inclusive_Jet_Mass,
			 Single_Inclusive_Jet_Mass_Getter,"MS")


Single_Inclusive_Jet_Mass::Single_Inclusive_Jet_Mass(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Single_Inclusive_Jet_Mass.dat");
}

void Single_Inclusive_Jet_Mass::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().heavyjetmass,weight,ncount);
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().lightjetmass,weight,ncount);
  }
}
 
Primitive_Observable_Base * Single_Inclusive_Jet_Mass::Copy() const
{
  return new Single_Inclusive_Jet_Mass(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


// ----------------------------------------------------------------------

DEFINE_OBSERVABLE_GETTER(Wide_Jet_Broadening,Wide_Jet_Broadening_Getter,"BW")

Wide_Jet_Broadening::Wide_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Wide_Jet_Broadening.dat");
}

void Wide_Jet_Broadening::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().widejetbroadening,weight,ncount);
  }
}
 
Primitive_Observable_Base * Wide_Jet_Broadening::Copy() const
{
  return new Wide_Jet_Broadening(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Narrow_Jet_Broadening,Narrow_Jet_Broadening_Getter,"BN")

Narrow_Jet_Broadening::Narrow_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Narrow_Jet_Broadening.dat");
}

void Narrow_Jet_Broadening::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().narrowjetbroadening,weight,ncount);
  }
}
 
Primitive_Observable_Base * Narrow_Jet_Broadening::Copy() const
{
  return new Narrow_Jet_Broadening(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Total_Jet_Broadening,Total_Jet_Broadening_Getter,"BT")

Total_Jet_Broadening::Total_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Total_Jet_Broadening.dat");
}

void Total_Jet_Broadening::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().widejetbroadening+
		    data->Get<JetMass_Broadening_Data>().narrowjetbroadening,
		    weight,ncount);
  }
}
 
Primitive_Observable_Base * Total_Jet_Broadening::Copy() const
{
  return new Total_Jet_Broadening(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


DEFINE_OBSERVABLE_GETTER(Jet_Broadening_Difference,Jet_Broadening_Difference_Getter,"BD")

Jet_Broadening_Difference::Jet_Broadening_Difference(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Jet_Broadening_Difference.dat");
}

void Jet_Broadening_Difference::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().widejetbroadening-
		    data->Get<JetMass_Broadening_Data>().narrowjetbroadening,
		    weight,ncount);
  }
}
 
Primitive_Observable_Base * Jet_Broadening_Difference::Copy() const
{
  return new Jet_Broadening_Difference(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_OBSERVABLE_GETTER(Single_Inclusive_Jet_Broadening,Single_Inclusive_Jet_Broadening_Getter,"BS")


Single_Inclusive_Jet_Broadening::Single_Inclusive_Jet_Broadening(int type, double xmin, double xmax, int nbin, std::string listname)
  : Primitive_Observable_Base(type,xmin,xmax,nbin), m_key(listname+"_JetMass_Broadening")
{
  m_listname = listname;
  m_name = std::string("Single_Inclusive_Jet_Broadening.dat");
}

void Single_Inclusive_Jet_Broadening::Evaluate(const ATOOLS::Blob_List & bl, double weight, double ncount)
{
  Blob_Data_Base * data = (*p_ana)[m_key];
  if (data) {
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().widejetbroadening,weight,ncount);
    p_histo->Insert(data->Get<JetMass_Broadening_Data>().narrowjetbroadening,weight,ncount);
  }
}
 
Primitive_Observable_Base * Single_Inclusive_Jet_Broadening::Copy() const
{
  return new Single_Inclusive_Jet_Broadening(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

namespace ATOOLS {

template <> Blob_Data<JetMass_Broadening_Data>::~Blob_Data() { }
template class Blob_Data<JetMass_Broadening_Data>;

}
