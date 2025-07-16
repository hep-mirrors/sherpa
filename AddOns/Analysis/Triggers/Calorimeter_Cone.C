#include "AddOns/Analysis/Triggers/Calorimeter_Cone.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Math/MathTools.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;


void Calorimeter_Cone::SetAnalysis(Primitive_Analysis  * ana)
{
  // give analysis
  Primitive_Detector * detector =  
    dynamic_cast<Primitive_Detector *>(ana->GetObject("Full Detector"));

  // get calorimeter
  if (detector) {
    p_calorimeter= 
      dynamic_cast<Primitive_Calorimeter *>(detector->GetElement("Hadronic Calorimeter"));

    if (p_calorimeter)
      p_calorimeter->GetDimensions(m_neta,m_nphi,m_mineta,m_maxeta);
    else {
      msg_Out()<<"WARNING  Calorimeter_Cone::Calorimeter_Cone no ""Hadronic Calorimeter"" "<<std::endl;
    }
  }
  else {
    msg_Out()<<"WARNING  Calorimeter_Cone::Calorimeter_Cone no ""Full Detector"" "<<std::endl;
  }
  if (p_jetno!=NULL) {
    for (int i=0; i<m_neta;++i) delete [] p_jetno[i];
    delete [] p_jetno;
  }
  p_jetno = new int*[m_neta];
  for (int i=0; i<m_neta;++i) p_jetno[i] = new int[m_nphi];
  m_delta_eta = (m_maxeta-m_mineta)/double(m_neta);
  m_delta_phi = 2.*M_PI/double(m_nphi);
}

Calorimeter_Cone::Calorimeter_Cone(const double Etcut,const double etamin, 
				   const double etamax,double sep) : 
  Jet_Algorithm_Base(NULL), m_dR(sep), m_dR2(sep*sep), m_Etcut(Etcut), m_Etstop(1.5), m_etamode(1), p_jetno(NULL)
{
  m_minetajet = etamin;
  m_maxetajet = etamax;
}

Calorimeter_Cone::~Calorimeter_Cone() 
{
  if (p_jetno) {
    for (int i=0; i<m_neta;++i) delete [] p_jetno[i];
    p_jetno=NULL;
  }
}

void Calorimeter_Cone::Test()
{
  Particle_List * pl = new Particle_List;
  pl->push_back(new Particle(1,Flavour(kf_p_plus),Vec4D(sqrt(500.),20.,10.,0.)));
  pl->push_back(new Particle(1,Flavour(kf_p_plus),Vec4D(50.,40.,0.,30.)));
  pl->push_back(new Particle(1,Flavour(kf_p_plus),Vec4D(10.,5.,3.,4.)));
  pl->push_back(new Particle(1,Flavour(kf_p_plus),Vec4D(50.,-40.,0.,30.)));
  pl->push_back(new Particle(1,Flavour(kf_p_plus),Vec4D(50.,20.,-20.,30.)));

  p_calorimeter->Fill(pl);
  p_calorimeter->Print();
  CalcJets();
  abort();
}

void Calorimeter_Cone::CalcJets()
{
  m_dneta     = int(m_neta*m_dR/(m_maxeta-m_mineta));
  m_dnphi     = int(m_nphi*m_dR/(2.*M_PI));
  
  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) {
      p_jetno[i][j]=0;
    }
  }
  m_jets.clear();
  double maxet, jetet, eta;
  double costheta, sintheta, cosphi, sinphi;
  Vec4D  jetmom;
  for (;;) {  
    // find highest tower
    maxet = 0;
    int ii=-1,jj=-1;
    for (int i=0; i<m_neta; ++i) {
      if (m_etamode==0) {
	if (m_mineta+i*m_delta_eta<m_minetajet ||
	    m_mineta+i*m_delta_eta>m_maxetajet) continue;
      }
      for (int j=0; j<m_nphi; ++j) {
	if (p_jetno[i][j]>0)                continue;
	if (p_calorimeter->Cell(i,j)<maxet) continue;
	maxet = p_calorimeter->Cell(i,j);
	ii = i; jj = j;
      }
    }
    if (ii==-1) break;
    if (maxet<m_Etstop) break;

    // add jet:
    jetet  = 0.;
    jetmom = Vec4D(0.,0.,0.,0.);
    for (int i=ii-m_dneta;i<=ii+m_dneta;++i) {
      if (i<0) i=0;
      if (i>=m_neta) break; 
      for (int jp=jj-m_dnphi;jp<=jj+m_dnphi;++jp) {
	int j=jp;
	if (j<0) j+=m_nphi;
	else if (j>=m_nphi) j-=m_nphi;

	double dr2 = sqr(m_delta_phi*(jp-jj))+sqr(m_delta_eta*(i-ii));
	if (dr2>m_dR2)       continue;
	if (p_jetno[i][j]>0) continue;
	
	p_jetno[i][j] = m_jets.size()+1;
	// add to jet
	double pt  = p_calorimeter->Cell(i,j);
	if (pt>0.) {
	  p_calorimeter->GetCosSinTheta(i,costheta,sintheta);
	  p_calorimeter->GetCosSinPhi(j,cosphi,sinphi);
	  double px  = pt/sintheta;
	  jetmom[0] += px;
	  jetmom[1] += pt*sinphi;
	  jetmom[2] += pt*cosphi;
	  jetmom[3] += px*costheta;
	  jetet     += pt;
	}
      }
    }
    if (jetet<m_Etcut) break;
    if (m_etamode==1) {
      eta = jetmom.Eta();
      if (!(eta>m_minetajet && eta<m_maxetajet)) continue; 
    }
    if (jetet>m_Etcut) {
      m_jets.push_back(Jet_Data(ii,jj,m_jets.size()+1,jetmom,jetet));
    }
  }
  SortPT();
}


bool Calorimeter_Cone::ConstructJets(const Particle_List *, Particle_List * jets,
				     std::vector<double> * kt2, double rmin)
{
  if (rmin!=-1.) {
    m_dR  = rmin;
    m_dR2 = rmin*rmin;
  }
  /*
  m_dneta     = int(m_neta*rmin/(m_maxeta-m_mineta));
  m_dnphi     = int(m_nphi*rmin/(2.*M_PI));
  */
  CalcJets();
  
  if (jets) {
    int i=1;
    for (std::vector<Jet_Data>::iterator it=m_jets.begin();it!=m_jets.end();++it,++i) {
      jets->push_back(new Particle(i,Flavour(kf_jet),it->mom));
      kt2->push_back(it->mom.PPerp2());
    }
  }
  return true;
}

void Calorimeter_Cone::SortPT()
{
  std::sort(m_jets.begin(),m_jets.end(),Order_PT_JData());
}

void Calorimeter_Cone::FillShape(int jetno,ATOOLS::Histogram * histo,
				 double weight,double ncount)
{
  if (jetno>(int)m_jets.size()) return;

  double dR, phi, et = m_jets[jetno-1].et;
  
  Vec4D  mom         = m_jets[jetno-1].mom;
  int    number      = m_jets[jetno-1].orig;
  double y           = p_calorimeter->PseudoRapidityNAzimuthalAngle(mom,phi);
  if (phi<0) phi    += 2.*M_PI;

  int ipos=m_jets[jetno-1].i, jpos=m_jets[jetno-1].j;

  /*
    int itest = int((y-m_mineta)/m_delta_eta);
    int jtest = int(phi/m_delta_phi);
    
    std::cout<<"Check this out : "<<ipos<<"/"<<jpos<<" -> "<<mom<<","<<et<<","
    <<itest<<"/"<<jtest<<"("<<number<<")"<<std::endl;
    p_calorimeter->Print();
  */
  
  double cells = 0.;

  //determine number of entries
  for (int i=ipos-m_dneta;i<=ipos+m_dneta;++i) {
    if (i<0) i=0;
    if (i>=m_neta) break; 
    for (int jp=jpos-m_dnphi;jp<=jpos+m_dnphi;++jp) {
      int j=jp;
      if (j<0) j+=m_nphi;
      else if (j>=m_nphi) j-=m_nphi;
      if (p_calorimeter->Cell(i,j)==0.) continue; 
      if (p_jetno[i][j]!=number) continue;
      cells+=1.;
    }
  }

  for (int i=ipos-m_dneta;i<=ipos+m_dneta;++i) {
    if (i<0) i=0;
    if (i>=m_neta) break; 
    for (int jp=jpos-m_dnphi;jp<=jpos+m_dnphi;++jp) {
      int j=jp;
      if (j<0) j+=m_nphi;
      else if (j>=m_nphi) j-=m_nphi;
      if (p_calorimeter->Cell(i,j)==0.) continue; 
      if (p_jetno[i][j]!=number) continue;
      double dphi = sqr(Min(dabs(j*m_delta_phi-phi),dabs(j*m_delta_phi-phi-2.*M_PI)));
      double deta = sqr(m_mineta+m_delta_eta*i-y);
      dR    = sqrt(dphi+deta);
      histo->Insert(dR,p_calorimeter->Cell(i,j)*weight/et*cells,ncount);
    }
  }
}
