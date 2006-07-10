#include "Four_Particle_Observables.H"
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
{ str<<"kf1 kf2 kf3 kf4 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)                           \
  CLASS::CLASS(const String_Matrix & parameters) :                      \
      Four_Particle_Observable_Base(parameters, TAG) { }               \
                                                                        \
  CLASS::CLASS(const CLASS * old) :                                     \
      Four_Particle_Observable_Base(*old) { }                          \
                                                                        \
  Primitive_Observable_Base * CLASS::Copy() const { return new CLASS(this); }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)                                             \
  DEFINE_CONSTRUCTOR_METHODS(CLASS,TAG)

Four_Particle_Observable_Base::Four_Particle_Observable_Base(const String_Matrix & parameters,
                                                           const std::string & obsname):
    Primitive_Observable_Base(parameters), f_special(false), m_flavs(4)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<8) {
      msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
          <<"observable "<<obsname<<" in Analysis.dat";
      abort();
    }
    for (short unsigned int i=0;i<4;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      m_flavs[i]=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
      if (kf<0) m_flavs[i]=m_flavs[i].Bar();
    }

    m_xmin  = ATOOLS::ToType<double>(parameters[0][4]);
    m_xmax  = ATOOLS::ToType<double>(parameters[0][5]);
    int nbins = ATOOLS::ToType<int>(parameters[0][6]);
    m_type  = HistogramType(parameters[0][7]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    MyStrStream str;
    str<<obsname<<m_flavs[0]<<m_flavs[1]<<m_flavs[2]<<m_flavs[3]<<".dat";
    str>>m_name;

    m_listname = parameters[0].size()>8?parameters[0][8]:"Analysed";
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      for (short unsigned int j=0;j<4;++j) {
        if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
          int kf=ATOOLS::ToType<int>(parameters[i][1]);
          m_flavs[j] = ATOOLS::Flavour((ATOOLS::kf::code)abs(kf),kf<0);
        }
      }
    }

    if (m_name=="SherpaDefault") {
      MyStrStream str;
      str<<obsname<<m_flavs[0]<<m_flavs[1]<<m_flavs[2]<<m_flavs[3];
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = obsname;
      MyStrStream str;
      str<<" of "<<m_flavs[0]<<", "<<m_flavs[1]<<", "<<m_flavs[2]<<", "<<m_flavs[3];
      str>>title;
      p_histo->SetTitle(title);
    }
  }
  if(m_xmin>=0.0) f_special=true;
}

Four_Particle_Observable_Base::Four_Particle_Observable_Base(const Four_Particle_Observable_Base * old) :
    Primitive_Observable_Base(*old)
{
  m_flavs=old->m_flavs; f_special=old->f_special;
}

void Four_Particle_Observable_Base::Evaluate(double value, double weight,
					     int ncount) {
  p_histo->Insert(value,weight,ncount); 
}

 
void Four_Particle_Observable_Base::Evaluate(int nout, const Vec4D* moms,
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
		  Evaluate(moms[i],moms[j],moms[k],moms[l],weight,ncount);
		}
	      }
	    }
	  }
	} 
      }
    }
  }
}


void Four_Particle_Observable_Base::Evaluate(const Particle_List& plist,
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
		  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),
			   (*plit3)->Momentum(),(*plit4)->Momentum(),
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
  p_histo->Insert(0.0,0.0,ncount);
}

// ============================================================================

#define DEFINE_PRINT_METHOD2(NAME)                                      \
  void NAME::PrintInfo(std::ostream &str,const size_t width) const      \
{ str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER2(CLASS,NAME,TAG)                       \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD2(NAME)

DEFINE_OBSERVABLE_GETTER2(Di_Mass,
                          Di_Mass_Getter,"Di_Mass")

Di_Mass::Di_Mass(const String_Matrix & parameters):
    Primitive_Observable_Base(parameters)
{
  if (parameters.size()==1) {
    if (parameters[0].size()==4 || parameters[0].size()==5) {
      m_xmin  = ATOOLS::ToType<double>(parameters[0][0]);
      m_xmax  = ATOOLS::ToType<double>(parameters[0][1]);
      int nbins = ATOOLS::ToType<int>(parameters[0][2]);
      m_type  = HistogramType(parameters[0][3]);
      p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);
      m_listname=parameters[0].size()>4?parameters[0][4]:"Analysed";

      m_name  = std::string("4jet_");
      if (m_listname!="Analysed") m_name=m_listname+std::string("_")+m_name;
      m_name += "Di_Mass.dat";
    }
    else {
      msg.Error()<<"Error in "<<METHOD<<": Wrong number of parameters for "
          <<"observable Di_Mass in Analysis.dat";
      abort();
    }
  }
  else {
    if (m_listname=="") m_listname="Analysed";
    if (m_name=="SherpaDefault") {
      m_name  = std::string("4jet_");
      if (m_listname!="Analysed") m_name=m_listname+std::string("_")+m_name;
      m_name += "Di_Mass";
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "DiMass";
      p_histo->SetTitle(title);
    }
  }
}

void Di_Mass::Evaluate(const ATOOLS::Blob_List & blobs,double weight, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);

  if (pl->size()!=4) {
    p_histo->Insert(0.,0.,2*ncount);
    return;
  }

  std::vector<Vec4D> moms;
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    moms.push_back((*pit)->Momentum());
  }

  double m1a = (moms[0]+moms[1]).Abs2();
  double m1b = (moms[2]+moms[3]).Abs2();
  double d1 = dabs(m1a-m1b);

  double m2a = (moms[0]+moms[2]).Abs2();
  double m2b = (moms[1]+moms[3]).Abs2();
  double d2 = dabs(m2a-m2b);

  double m3a = (moms[0]+moms[3]).Abs2();
  double m3b = (moms[1]+moms[2]).Abs2();
  double d3 = dabs(m3a-m3b);

  if (d1<d2 && d1<d3) {
    p_histo->Insert(sqrt(m1a),weight,ncount);
    p_histo->Insert(sqrt(m1b),weight,ncount);
  }
  else if (d2<d3) {
    p_histo->Insert(sqrt(m2a),weight,ncount);
    p_histo->Insert(sqrt(m2b),weight,ncount);
  }
  else {
    p_histo->Insert(sqrt(m3a),weight,ncount);
    p_histo->Insert(sqrt(m3b),weight,ncount);
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Four_Particle_PlaneAngle,
			 Four_Particle_PlaneAngle_Getter,"PlaneAngle")

void Four_Particle_PlaneAngle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
					const Vec4D & mom3,const Vec4D & mom4,
					double weight, int ncount)
{
  Vec3D normal1 = cross(Vec3D(mom1),Vec3D(mom2));
  Vec3D normal2 = cross(Vec3D(mom3),Vec3D(mom4));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs()); 
  p_histo->Insert(costh,weight,ncount); 
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Four_Particle_PT,
			 Four_Particle_PT_Getter,"PT4")

void Four_Particle_PT::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				const Vec4D& mom3,const Vec4D& mom4,
				double weight, int ncount)
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]+mom3[1]+mom4[1]) +
		   sqr(mom1[2]+mom2[2]+mom3[2]+mom4[2]));
  p_histo->Insert(pt,weight,ncount);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partonpair_PTdiff,
			 Two_Partonpair_PTdiff_Getter,"PTdiff4")

void Two_Partonpair_PTdiff::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				     const Vec4D& mom3,const Vec4D& mom4,
				     double weight, int ncount) {
  Vec4D vecA(mom1); vecA+=mom2;
  Vec4D vecB(mom3); vecB+=mom4;
  double ptdiff=sqrt(sqr(vecA[1])+sqr(vecA[2]));
  ptdiff-=sqrt(sqr(vecB[1])+sqr(vecB[2]));
  if(f_special) ptdiff=dabs(ptdiff);
  p_histo->Insert(ptdiff, weight, ncount);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partonpair_Theta,
			 Two_Partonpair_Theta_Getter,"Theta4")

void Two_Partonpair_Theta::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				    const Vec4D& mom3,const Vec4D& mom4,
				    double weight, int ncount) {
  Vec4D vecA(mom1); vecA+=mom2;
  Vec4D vecB(mom3); vecB+=mom4;
  //
  if(!f_special) {
    Vec4D plab(vecA); plab+=vecB;
    //removing the z boost effect of the considered 4 particle system
    plab[1]=plab[2]=0.0;
    if(plab.Abs2()<=0.0) {
      p_histo->Insert(-M_PI/100.0, weight, ncount);
      msg.Error()<<__PRETTY_FUNCTION__<<":\n   Warning:"
		 <<" Not able to boost the system. Insert theta=-pi/100.\n"
		 <<std::endl;
      return;
    }
    Poincare fly(plab);
    fly.Boost(vecA);
    fly.Boost(vecB);
  }
  Vec3D vec1(vecA);
  Vec3D vec2(vecB);
  double theta=acos((vec1*vec2)/(vec1.Abs()*vec2.Abs()));
  p_histo->Insert(theta, weight, ncount);
}
