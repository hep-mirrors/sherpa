#include "Calorimeter_Cone.H"
#include "Particle_List.H"
#include "MathTools.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Calorimeter_Cone::Calorimeter_Cone(const double dR,const double Etcut,
				   Primitive_Calorimeter * const calorimeter) : 
  m_dR(dR), m_dR2(dR*dR), m_Etcut(Etcut), m_Etstop(1.5), m_etamode(1), 
  p_calorimeter(calorimeter)
{
  p_calorimeter->GetDimensions(m_neta,m_nphi,m_mineta,m_maxeta);
  m_minetajet = m_mineta;
  m_maxetajet = m_maxeta;
  m_dneta     = int(m_neta*m_dR/(m_maxeta-m_mineta));
  m_dnphi     = int(m_nphi*m_dR/(2.*M_PI));
  m_delta_eta = (m_maxeta-m_mineta)/double(m_neta);
  m_delta_phi = 2.*M_PI/double(m_nphi);

  p_jetno = new int*[m_neta];
  for (int i=0; i<m_neta;++i) p_jetno[i] = new int[m_nphi];
}

Calorimeter_Cone::~Calorimeter_Cone() 
{
  if (p_jetno) {
    for (int i=0; i<m_neta;++i) delete [] p_jetno[i];
    p_jetno=NULL;
  }
}


void  Calorimeter_Cone::CalcJets()
{
  for (int i=0; i<m_neta; ++i) {
    for (int j=0; j<m_nphi; ++j) {
      p_jetno[i][j]=0;
    }
  }
  m_jets.clear();
  double maxet, jetet;
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
	if (p_jetno[i][j]>0)                 continue;
	if (p_calorimeter->Cell(i,j)<maxet) continue;
	maxet = p_calorimeter->Cell(i,j);
	ii = i; jj = j;
      }
    }
    if (ii==-1) break;
    if (maxet<m_Etstop) break;

    // add jet:
    jetet = 0.;

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
	p_calorimeter->GetCosSinTheta(i,costheta,sintheta);
	p_calorimeter->GetCosSinPhi(j,cosphi,sinphi);
	double px  = pt/sintheta;
	jetmom[0] += px;
	jetmom[1] += pt*cosphi;
	jetmom[2] += pt*sinphi;
	jetmom[3] += px*costheta;
	jetet     += pt;
      }
    }
    if (jetet<m_Etcut) break;
    m_jets.push_back(Jet_Data(ii,jj,jetmom,jetet));
  }
}


bool  Calorimeter_Cone::ConstructJets(Particle_List * jets,std::vector<double> * kt2)
{
  CalcJets();
  int i=1;
  double eta;
  for (std::vector<Jet_Data>::iterator it=m_jets.begin();it!=m_jets.end();++it,++i) {
    if (m_etamode==1) {
      eta = it->mom.Eta();
      if (eta>m_minetajet && eta<m_maxetajet) 
	jets->push_back(new Particle(i,Flavour(kf::jet),it->mom));
    }
    else jets->push_back(new Particle(i,Flavour(kf::jet),it->mom));
  }    

  SortPT(jets);
  for (Particle_Iterator pit=jets->begin();pit!=jets->end();++pit) {
    kt2->push_back((*pit)->Momentum().PPerp2());
  }
  return true;
}
