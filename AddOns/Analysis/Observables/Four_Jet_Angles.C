#include "AddOns/Analysis/Observables/Four_Jet_Angles.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

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
  ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,CLASS>::operator()(const Argument_Matrix &parameters) const \
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void ATOOLS::Getter<Primitive_Observable_Base,Argument_Matrix,NAME>::PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(CLASS)

using namespace ATOOLS;

Four_Jet_Angle_Base::Four_Jet_Angle_Base(unsigned int type,double xmin,double xmax,int nbins,
					 const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_listname=lname;
  m_name  = std::string("4jet_");
  if (lname!=finalstate_list) m_name=lname+std::string("_")+m_name;
}

void Four_Jet_Angle_Base::Evaluate(const Blob_List & blobs,double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);

  // sort in Durham (and NOT in Final_Selector!)

  std::vector<Vec3D> moms;
 
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    moms.push_back((*pit)->Momentum());
  }

  if (moms.size()!=4) {
    p_histo->Insert(0.,0.,ncount);
    return;
  }

//   for (size_t i=0;i<4;++i)
//     std::cout<<" p["<<i<<"]="<<moms[i]<<std::endl;

  double cos_chi=Calc(moms);
  if (p_histo->Xmin()==0.) cos_chi=dabs(cos_chi);
//   std::cout<<Name()<<" cos_chi="<<cos_chi<<" ("<<p_histo->Xmin()<<")"<<std::endl;
  p_histo->Insert(cos_chi,value,ncount);
}

void Four_Jet_Angle_Base::EvaluateNLOcontrib(double value, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);

  // sort in Durham (and NOT in Final_Selector!)

  std::vector<Vec3D> moms;
 
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    moms.push_back((*pit)->Momentum());
  }

  if (moms.size()!=4) {
    p_histo->InsertMCB(0.,0.,ncount);
    return;
  }


  double cos_chi=Calc(moms);
  if (p_histo->Xmin()==0.) cos_chi=dabs(cos_chi);
  p_histo->InsertMCB(cos_chi,value,ncount);
}
 
void Four_Jet_Angle_Base::EvaluateNLOevt() 
{
  p_histo->FinishMCB();
}


// ======================================================================

DEFINE_OBSERVABLE_GETTER(Bengtsson_Zerwas_Angle,
			 Bengtsson_Zerwas_Angle_Getter,"BZAngle")

Bengtsson_Zerwas_Angle::Bengtsson_Zerwas_Angle(unsigned int type,double xmin,double xmax,int nbins,
					       const std::string & lname) :
  Four_Jet_Angle_Base(type,xmin,xmax,nbins,lname)
{
  m_name+=std::string("BZ.dat");
}

double Bengtsson_Zerwas_Angle::Calc(const std::vector<ATOOLS::Vec3D> & moms)
{
  Vec3D p12=cross(moms[0],moms[1]);
  Vec3D p34=cross(moms[2],moms[3]);

  return (p12*p34)/(p12.Abs()*p34.Abs());
}

Primitive_Observable_Base * Bengtsson_Zerwas_Angle::Copy() const 
{
  return new Bengtsson_Zerwas_Angle(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Nachtmann_Reiter_Angle,
			 Nachtmann_Reiter_Angle_Getter,"NRAngle")

Nachtmann_Reiter_Angle::Nachtmann_Reiter_Angle(unsigned int type,double xmin,double xmax,int nbins,
					       const std::string & lname) :
  Four_Jet_Angle_Base(type,xmin,xmax,nbins,lname)
{
  m_name+=std::string("mNR.dat"); // modified Nachtmann Reiter
}

double Nachtmann_Reiter_Angle::Calc(const std::vector<ATOOLS::Vec3D> & moms)
{
  Vec3D p1_2=moms[0]-moms[1];
  Vec3D p3_4=moms[2]-moms[3];

  return (p1_2*p3_4)/(p1_2.Abs()*p3_4.Abs());
}

Primitive_Observable_Base * Nachtmann_Reiter_Angle::Copy() const 
{
  return new Nachtmann_Reiter_Angle(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Koerner_Schierholz_Willrodt_Angle,
			 Koerner_Schierholz_Willrodt_Angle_Getter,"KSWAngle")

Koerner_Schierholz_Willrodt_Angle::Koerner_Schierholz_Willrodt_Angle(unsigned int type,double xmin,double xmax,int nbins,
					       const std::string & lname) :
  Four_Jet_Angle_Base(type,xmin,xmax,nbins,lname)
{
  m_name+=std::string("KSW.dat");
}

double Koerner_Schierholz_Willrodt_Angle::Calc(const std::vector<ATOOLS::Vec3D> & moms)
{
  Vec3D p14=cross(moms[0],moms[3]);
  Vec3D p23=cross(moms[1],moms[2]);
  double c1423= (p14*p23)/(p14.Abs()*p23.Abs());
  Vec3D p13=cross(moms[0],moms[2]);
  Vec3D p24=cross(moms[1],moms[3]);
  double c1324= (p13*p24)/(p13.Abs()*p24.Abs());
  //  std::cout<<" phi1="<<acos(c1423)<<"  phi2="<<acos(c1324)<<std::endl;
  
  double phi_ksw=0.5*(acos(c1423)+acos(c1324));
  return cos(phi_ksw);
}

Primitive_Observable_Base * Koerner_Schierholz_Willrodt_Angle::Copy() const 
{
  return new Koerner_Schierholz_Willrodt_Angle(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

// ======================================================================

DEFINE_OBSERVABLE_GETTER(Alpha34_Angle,Alpha34_Angle_Getter,"A34Angle")

Alpha34_Angle::Alpha34_Angle(unsigned int type,double xmin,double xmax,int nbins,
					       const std::string & lname) :
  Four_Jet_Angle_Base(type,xmin,xmax,nbins,lname)
{
  m_name+=std::string("a34.dat");
}

double Alpha34_Angle::Calc(const std::vector<ATOOLS::Vec3D> & moms)
{
  return (moms[2]*moms[3])/(moms[2].Abs()*moms[3].Abs());
}

Primitive_Observable_Base * Alpha34_Angle::Copy() const 
{
  return new Alpha34_Angle(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
