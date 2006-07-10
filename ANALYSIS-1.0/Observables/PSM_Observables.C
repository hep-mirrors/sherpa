#include "PSM_Observables.H"
#include "MyStrStream.H"
#include "Primitive_Analysis.H"
#include "Shell_Tools.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_GETTER_METHOD(CLASS,NAME)                                \
  Primitive_Observable_Base *                                   \
  NAME::operator()(const String_Matrix &parameters) const               \
  { return new CLASS(parameters); }

#define DEFINE_PRINT_METHOD(NAME)                                       \
  void NAME::PrintInfo(ostream &str,const size_t width) const   \
  { str<<"min max bins pn0 pn1 pn2 pn3 Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)                        \
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);     \
  DEFINE_GETTER_METHOD(CLASS,NAME)                                      \
  DEFINE_PRINT_METHOD(NAME)

DEFINE_OBSERVABLE_GETTER(PSM_Observable,PSM_Observable_Getter,"PSM")

PSM_Observable::PSM_Observable(const String_Matrix & parameters):
  Primitive_Observable_Base(parameters), m_pnb(4)
{
  if (parameters.size()==1) {
    if (parameters[0].size()<8) {
    msg.Error()<<"Error in "<<METHOD<<": Not enough parameters for "
        <<"observable PSM_Observable in Analysis.dat";
      abort();
    }
    m_xmin  = ToType<double>(parameters[0][0]);
    m_xmax  = ToType<double>(parameters[0][1]);
    int nbins = ToType<int>(parameters[0][2]);
    m_type  = HistogramType(parameters[0][3]);
    p_histo = new Histogram(m_type,m_xmin,m_xmax,nbins);

    m_pnb.clear();
    for(int i=3;i<7;i++) {
      m_pnb.push_back(ToType<int>(parameters[0][i]));
    }

    m_listname = parameters[0].size()>8?parameters[0][8]:"Analysed";
    MyStrStream str;
    if (m_listname!="Analysed") str<<m_listname<<string("_");
    str<<"PSM_";
    for (size_t i=0;i<m_pnb.size();i++) str<<m_pnb[i];
    str>>m_name;
  }
  else {
    for (size_t i=0;i<parameters.size();++i) {
      if (parameters[i].size()<2) continue;
      for (short unsigned int j=0;j<5;++j) {
        if (parameters[i][0]==std::string("PN")+ToString(j)) {
          m_pnb[i] = ToType<int>(parameters[i][1]);
        }
      }
    }

    if (m_listname=="") m_listname="Analysed";
    if (m_name=="SherpaDefault") {
      MyStrStream str;
      if (m_listname!="Analysed") str<<m_listname<<string("_");
      str<<"PSM_";
      for (size_t i=0;i<m_pnb.size();i++) str<<m_pnb[i];
      str>>m_name;
    }
    if (p_histo->Title()=="SherpaDefault") {
      std::string title = "";
      MyStrStream str;
      str<<"PSM ";
      for (size_t i=0;i<m_pnb.size();i++) str<<m_pnb[i];
      str>>title;
      p_histo->SetTitle(title);
    }
  }
}

PSM_Observable::PSM_Observable(const PSM_Observable * old) :
Primitive_Observable_Base(*old)
{
  m_pnb = old->m_pnb;
}

Primitive_Observable_Base * PSM_Observable::Copy() const { return new PSM_Observable(this); }

void PSM_Observable::Evaluate(const Particle_List & pl,double weight, int ncount)
{
  std::vector<Vec4D> moms;
  Vec4D smom(0.,0.,0.,0.);
  for (Particle_List::const_iterator it=pl.begin();it!=pl.end();++it) {
    smom+=(*it)->Momentum();
  }
  moms.push_back(Vec4D(0.5*(smom[0]+smom[3]),0.,0.,0.5*(smom[0]+smom[3])));
  moms.push_back(Vec4D(0.5*(smom[0]-smom[3]),0.,0.,-0.5*(smom[0]-smom[3])));
  for (Particle_List::const_iterator it=pl.begin();it!=pl.end();++it) {
    moms.push_back((*it)->Momentum());
  }
  
  Vec4D ms=Vec4D(0.,0.,0.,0.);
  if (m_pnb.size()>0) {
    for (size_t i=0;i<moms.size();i++){
      int hit=0;
      for(size_t j=0;j<m_pnb.size();j++) {
	if (m_pnb[j]==(int)i) hit = 1;
      }
      if (hit) {
	if (i<2) ms -= moms[i];
	else ms += moms[i];
      }
    } 
    double st=ms.Abs2();
    if (st<0.) st=-sqrt(-st);
    else st=sqrt(st);
    p_histo->Insert(st,weight,ncount);
  }
  else {
    ms = moms[0]+moms[1];
    double y = 0.5 * log( (ms[0]+ms[3])/(ms[0]-ms[3]) );
    p_histo->Insert(y,weight,ncount);
  }
}

void PSM_Observable::Evaluate(const Blob_List & blobs,double value, int ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);
  Evaluate(*pl,value, ncount);
}

// void PSM_Observable::EndEvaluation(double scale) {
//     p_histo->Finalize();
// }

// void PSM_Observable::Output(const std::string & pname) {
//   int  mode_dir = 448;
//   ATOOLS::MakeDir((pname).c_str(),mode_dir); 
//   p_histo->Output((pname + std::string("/") + m_name+std::string(".dat")).c_str());
// }

/*
Primitive_Observable_Base & PSM_Observable::operator+=(const Primitive_Observable_Base & ob)
{
 PSM_Observable * cob = ((PSM_Observable*)(&ob));
 if (p_histo) {
   //    if (cob->p_histo->Type()>99) (*p_histo)+=(*(static_cast<Histogram *>(cob->p_histo)));
   //                            else (*p_histo)+=(*(static_cast<Root_Histogram *>(cob->p_histo)));

    //(*p_histo)+=(cob->p_histo);
  }
  else {
    msg.Out()<<" warning "<<Name()<<" has not overloaded the operator+="<<std::endl;
  }

  return *this;
}
*/

// void PSM_Observable::Reset()
// {
//   p_histo->Reset();
// }
