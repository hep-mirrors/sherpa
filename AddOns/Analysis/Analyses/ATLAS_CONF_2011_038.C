#include "AddOns/Analysis/Analyses/Analysis_Base.H"

namespace ANALYSIS {

  class ATLAS_CONF_2011_038_Analysis: public Analysis_Base {  
  private:

    std::string m_jetslist;

  public:

    ATLAS_CONF_2011_038_Analysis(const std::string &listname);

    void Evaluate(double weight,double ncount,int mode);
    Primitive_Observable_Base *Copy() const;

  };// end of class ATLAS_CONF_2011_038_Analysis

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

ATLAS_CONF_2011_038_Analysis::ATLAS_CONF_2011_038_Analysis
(const std::string &listname): Analysis_Base(listname)
{
  m_name="ATLAS_CONF_2011_038__"+m_listname;
  m_dists.resize(48,NULL);
  m_dists[0] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_A1",1);
  m_dists[1] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_A2",1);
  m_dists[2] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_A3",1);
  m_dists[3] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_A4",1);
  m_dists[4] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_A5",1);
  m_dists[5] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_B1",1);
  m_dists[6] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_B2",1);
  m_dists[7] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_B3",1);
  m_dists[8] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_B4",1);
  m_dists[9] = new Normalized_Observable(4,50.0,500.0,45,"GapFrac_vs_pT_B5",1);
  m_dists[10] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A1",1);
  m_dists[11] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A2",1);
  m_dists[12] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A3",1);
  m_dists[13] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A4",1);
  m_dists[14] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A5",1);
  m_dists[15] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A6",1);
  m_dists[16] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_A7",1);
  m_dists[17] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B1",1);
  m_dists[18] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B2",1);
  m_dists[19] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B3",1);
  m_dists[20] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B4",1);
  m_dists[21] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B5",1);
  m_dists[22] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B6",1);
  m_dists[23] = new Normalized_Observable(4,0.0,6.0,12,"GapFrac_vs_DY_B7",1);
  m_dists[24] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_A1",1);
  m_dists[25] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_A2",1);
  m_dists[26] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_A3",1);
  m_dists[27] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_A4",1);
  m_dists[28] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_A5",1);
  m_dists[29] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_B1",1);
  m_dists[30] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_B2",1);
  m_dists[31] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_B3",1);
  m_dists[32] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_B4",1);
  m_dists[33] = new Normalized_Observable(4,50.0,500.0,45,"NJets_vs_pT_B5",1);
  m_dists[34] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A1",1);
  m_dists[35] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A2",1);
  m_dists[36] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A3",1);
  m_dists[37] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A4",1);
  m_dists[38] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A5",1);
  m_dists[39] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A6",1);
  m_dists[40] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_A7",1);
  m_dists[41] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B1",1);
  m_dists[42] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B2",1);
  m_dists[43] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B3",1);
  m_dists[44] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B4",1);
  m_dists[45] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B5",1);
  m_dists[46] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B6",1);
  m_dists[47] = new Normalized_Observable(4,0.0,6.0,12,"NJets_vs_DY_B7",1);
}

class ATLAS_CONF_2011_038_Order_PT {
public:
  bool operator()(const Particle *a,const Particle *b)
  { return a->Momentum().PPerp2()>b->Momentum().PPerp2(); }
};// end of class ATLAS_CONF_2011_038_Order_PT

class ATLAS_CONF_2011_038_Order_Y {
public:
  bool operator()(const Particle *a,const Particle *b)
  { return a->Momentum().Y()<b->Momentum().Y(); }
};// end of class ATLAS_CONF_2011_038_Order_Y

void ATLAS_CONF_2011_038_Analysis::Evaluate(double weight, double ncount,int mode)
{
  Particle_List jets(*p_ana->GetParticleList(m_listname));
  if (jets.size()<2) AddZero(ncount,mode);
  // Selection A
  {
    std::stable_sort(jets.begin(),jets.end(),ATLAS_CONF_2011_038_Order_PT());
    double y1(jets[0]->Momentum().Y()), y2(jets[1]->Momentum().Y());
    double ymin(Min(y1,y2)), ymax(Max(y1,y2));
    double njingap(0.0);
    for (size_t i(2);i<jets.size();++i) {
      double y(jets[i]->Momentum().Y());
      if (y>ymin && y<ymax) ++njingap;
    }
    double ptbarA((jets[0]->Momentum().PPerp()+jets[1]->Momentum().PPerp())/2.0);
    double gap(njingap?0.0:1.0), dy(ymax-ymin);
    FillDist(0,ptbarA,gap,(dy>1.0&&dy<2.0)?weight:0.0,ncount,mode);
    FillDist(1,ptbarA,gap,(dy>2.0&&dy<3.0)?weight:0.0,ncount,mode);
    FillDist(2,ptbarA,gap,(dy>3.0&&dy<4.0)?weight:0.0,ncount,mode);
    FillDist(3,ptbarA,gap,(dy>4.0&&dy<5.0)?weight:0.0,ncount,mode);
    FillDist(4,ptbarA,gap,(dy>5.0&&dy<6.0)?weight:0.0,ncount,mode);
    FillDist(10,dy,gap,(ptbarA>70.0&&ptbarA<90.0)?weight:0.0,ncount,mode);
    FillDist(11,dy,gap,(ptbarA>90.0&&ptbarA<120.0)?weight:0.0,ncount,mode);
    FillDist(12,dy,gap,(ptbarA>120.0&&ptbarA<150.0)?weight:0.0,ncount,mode);
    FillDist(13,dy,gap,(ptbarA>150.0&&ptbarA<180.0)?weight:0.0,ncount,mode);
    FillDist(14,dy,gap,(ptbarA>180.0&&ptbarA<210.0)?weight:0.0,ncount,mode);
    FillDist(15,dy,gap,(ptbarA>210.0&&ptbarA<240.0)?weight:0.0,ncount,mode);
    FillDist(16,dy,gap,(ptbarA>240.0&&ptbarA<270.0)?weight:0.0,ncount,mode);
    FillDist(24,ptbarA,njingap,(dy>1.0&&dy<2.0)?weight:0.0,ncount,mode);
    FillDist(25,ptbarA,njingap,(dy>2.0&&dy<3.0)?weight:0.0,ncount,mode);
    FillDist(26,ptbarA,njingap,(dy>3.0&&dy<4.0)?weight:0.0,ncount,mode);
    FillDist(27,ptbarA,njingap,(dy>4.0&&dy<5.0)?weight:0.0,ncount,mode);
    FillDist(28,ptbarA,njingap,(dy>5.0&&dy<6.0)?weight:0.0,ncount,mode);
    FillDist(34,dy,njingap,(ptbarA>70.0&&ptbarA<90.0)?weight:0.0,ncount,mode);
    FillDist(35,dy,njingap,(ptbarA>90.0&&ptbarA<120.0)?weight:0.0,ncount,mode);
    FillDist(36,dy,njingap,(ptbarA>120.0&&ptbarA<150.0)?weight:0.0,ncount,mode);
    FillDist(37,dy,njingap,(ptbarA>150.0&&ptbarA<180.0)?weight:0.0,ncount,mode);
    FillDist(38,dy,njingap,(ptbarA>180.0&&ptbarA<210.0)?weight:0.0,ncount,mode);
    FillDist(39,dy,njingap,(ptbarA>210.0&&ptbarA<240.0)?weight:0.0,ncount,mode);
    FillDist(40,dy,njingap,(ptbarA>240.0&&ptbarA<270.0)?weight:0.0,ncount,mode);
  }
  // Selection B
  {
    std::stable_sort(jets.begin(),jets.end(),ATLAS_CONF_2011_038_Order_Y());
    double ymin(jets.front()->Momentum().Y()), ymax(jets.back()->Momentum().Y());
    double njingap(jets.size()-2);
    double ptbarB((jets.front()->Momentum().PPerp()+
		   jets.back()->Momentum().PPerp())/2.0);
    double gap(njingap?0.0:1.0), dy(ymax-ymin);
    FillDist(5,ptbarB,gap,(dy>1.0&&dy<2.0)?weight:0.0,ncount,mode);
    FillDist(6,ptbarB,gap,(dy>2.0&&dy<3.0)?weight:0.0,ncount,mode);
    FillDist(7,ptbarB,gap,(dy>3.0&&dy<4.0)?weight:0.0,ncount,mode);
    FillDist(8,ptbarB,gap,(dy>4.0&&dy<5.0)?weight:0.0,ncount,mode);
    FillDist(9,ptbarB,gap,(dy>5.0&&dy<6.0)?weight:0.0,ncount,mode);
    FillDist(17,dy,gap,(ptbarB>70.0&&ptbarB<90.0)?weight:0.0,ncount,mode);
    FillDist(18,dy,gap,(ptbarB>90.0&&ptbarB<120.0)?weight:0.0,ncount,mode);
    FillDist(19,dy,gap,(ptbarB>120.0&&ptbarB<150.0)?weight:0.0,ncount,mode);
    FillDist(20,dy,gap,(ptbarB>150.0&&ptbarB<180.0)?weight:0.0,ncount,mode);
    FillDist(21,dy,gap,(ptbarB>180.0&&ptbarB<210.0)?weight:0.0,ncount,mode);
    FillDist(22,dy,gap,(ptbarB>210.0&&ptbarB<240.0)?weight:0.0,ncount,mode);
    FillDist(23,dy,gap,(ptbarB>240.0&&ptbarB<270.0)?weight:0.0,ncount,mode);
    FillDist(29,ptbarB,njingap,(dy>1.0&&dy<2.0)?weight:0.0,ncount,mode);
    FillDist(30,ptbarB,njingap,(dy>2.0&&dy<3.0)?weight:0.0,ncount,mode);
    FillDist(31,ptbarB,njingap,(dy>3.0&&dy<4.0)?weight:0.0,ncount,mode);
    FillDist(32,ptbarB,njingap,(dy>4.0&&dy<5.0)?weight:0.0,ncount,mode);
    FillDist(33,ptbarB,njingap,(dy>5.0&&dy<6.0)?weight:0.0,ncount,mode);
    FillDist(41,dy,njingap,(ptbarB>70.0&&ptbarB<90.0)?weight:0.0,ncount,mode);
    FillDist(42,dy,njingap,(ptbarB>90.0&&ptbarB<120.0)?weight:0.0,ncount,mode);
    FillDist(43,dy,njingap,(ptbarB>120.0&&ptbarB<150.0)?weight:0.0,ncount,mode);
    FillDist(44,dy,njingap,(ptbarB>150.0&&ptbarB<180.0)?weight:0.0,ncount,mode);
    FillDist(45,dy,njingap,(ptbarB>180.0&&ptbarB<210.0)?weight:0.0,ncount,mode);
    FillDist(46,dy,njingap,(ptbarB>210.0&&ptbarB<240.0)?weight:0.0,ncount,mode);
    FillDist(47,dy,njingap,(ptbarB>240.0&&ptbarB<270.0)?weight:0.0,ncount,mode);
  }
}

Primitive_Observable_Base *ATLAS_CONF_2011_038_Analysis::Copy() const 
{
  return new ATLAS_CONF_2011_038_Analysis(m_listname);
}

DECLARE_ND_GETTER(ATLAS_CONF_2011_038_Analysis_Getter,"ATLAS-CONF-2011-038",
		  Primitive_Observable_Base,Argument_Matrix,true);

Primitive_Observable_Base *
ATLAS_CONF_2011_038_Analysis_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new ATLAS_CONF_2011_038_Analysis(parameters[0][0]);
}

void ATLAS_CONF_2011_038_Analysis_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list";
}

