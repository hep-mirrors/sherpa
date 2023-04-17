#include "SHRiMPS/Event_Generation/Collinear_Emission_Generator.H"

#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Collinear_Emission_Generator::Collinear_Emission_Generator() :
  m_kt2min(MBpars.GetShowerLinkParameters().CEKT2min)
{
  p_alphaS = new Strong_Coupling(static_cast<Running_AlphaS *>
				 (s_model->GetScalarFunction(string("alpha_S"))),
				 asform::frozen,
				 m_kt2min);
}

Collinear_Emission_Generator::~Collinear_Emission_Generator() {
  delete p_alphaS;
}

int Collinear_Emission_Generator::GenerateEmissions(Blob_List * blobs) {
  if (!InitNewRound(blobs->FindLast(btp::Hard_Collision))) return 0;
  AddEmissions();
  blobs->push_back(MakeBlob());
  return 1;
}

bool Collinear_Emission_Generator::InitNewRound(ATOOLS::Blob * hard_blob) {
  if (!hard_blob || hard_blob->GetInParticles().size() != 2) return false;
  m_sladder   = (hard_blob->InParticle(0)->Momentum()+hard_blob->InParticle(1)->Momentum()).Abs2();
  m_inparts   = hard_blob->GetOutParticles();
  m_beamparts = hard_blob->GetInParticles();
  m_outparts.clear();
  m_kt2starts.clear();
  for (size_t i=0; i<m_inparts.size(); ++i ) {
    Particle * part(new Particle(-1,m_inparts[i]->Flav(),m_inparts[i]->Momentum(),m_inparts[i]->Info()));
    //if (dabs(part->Momentum().Abs2())>1.e-6) msg_Out()<<(*part)<<"\n";
    part->SetNumber(0);
    part->SetFlow(1,m_inparts[i]->GetFlow(1));
    part->SetFlow(2,m_inparts[i]->GetFlow(2));
    m_outparts.push_back(part);
    m_kt2starts[part->Number()] = m_sladder;
  }
  hard_blob->UnsetStatus(blob_status::needs_showers);
  return true;
}

void Collinear_Emission_Generator::AddEmissions() {
  bool cont;
  do {
    cont = false;
    size_t Npart(m_outparts.size());
    for (size_t i=0; i<Npart; ++i ) {
      Particle * part(m_outparts[i]), * spec(FindSpectator(part));
      if (!spec) continue;
      m_kt2max = m_kt2starts[part->Number()];
      bool   isgluon(part->Flav().Kfcode() == kf_gluon);
      double kt2,z;
      GetKt2(isgluon,kt2,z);
      //msg_Out()<<METHOD<<"(kt2(max) = "<<m_kt2max<<" --> kt2 = "<<kt2<<") for "
      //       <<"s = "<<(part->Momentum()+spec->Momentum()).Abs2()<<"\n";
      if (kt2 > 0.) {
	cont = true;
	Vec4D splitmom = part->Momentum(), specmom = spec->Momentum(), gluonmom;
	if(FixKinematics(kt2,z,splitmom,gluonmom,specmom)) {
	  part->SetMomentum(splitmom);
	  spec->SetMomentum(specmom);
	  Particle * gluon(new Particle(-1, Flavour(kf_gluon), gluonmom, 'F'));
	  gluon->SetNumber(0);
	  bool trip;
	  if (isgluon){
	    if (ran->Get() < 0.) trip = true;
	    else trip = false;
	  }
	  else {
	    if (part->Flav().IsAnti()) trip = false;
	    else trip = true;
	  }
	  if (trip) {
	    gluon->SetFlow(1,part->GetFlow(1));
	    gluon->SetFlow(2,-1);
	    part->SetFlow(1,gluon->GetFlow(2));
	  }
	  else {
	    gluon->SetFlow(2,part->GetFlow(2));
	    gluon->SetFlow(1,-1);
	    part->SetFlow(2,gluon->GetFlow(1));
	  }
	  m_outparts.push_back(gluon);
	  m_kt2starts[gluon->Number()] = kt2;
	  /*msg_Out()<<"------------------ particles after splitting:\n";
	    totmom = Vec4D(0.,0.,0.,0.);
	    for (PVIt piter2 = m_outparts.begin(); piter2 != m_outparts.end(); ++piter2){
	    msg_Out()<<**piter2<<endl;
	    totmom += (*piter2)->Momentum();
	    }
	    msg_Out()<<" --> total momentum: "<<totmom<<endl;
	    msg_Out()<<"------------------------------------------\n";*/
	}
      }
      m_kt2starts[part->Number()] = kt2;
    }
  } while (cont);
}

void Collinear_Emission_Generator::GetKt2(bool isgluon, double &kt2, double &z){
  if (m_kt2max < m_kt2min) {
    kt2 = 0.;
    return;
  }
  double pref,asmax((*p_alphaS)(m_kt2min));
  pref = asmax/(2.*M_PI);
  if (isgluon) pref *= 1./2.;
  else pref *= 2./3.;
  double kt2max,weight;
  kt2 = m_kt2max;
  do {
    kt2max = kt2;
    double R(ran->Get());
    kt2 = kt2max*pow(R,1./pref);
    if (kt2 < m_kt2min) {
      kt2 = 0.;
      return;
    }
    z = GetZ(isgluon);
    if (z < 0.5*(1.-sqrt(1.-m_kt2min/kt2)) || z > 0.5*(1.+sqrt(1.-m_kt2min/kt2))) {
      weight = 0.;
    }
    else {
      weight = (*p_alphaS)(kt2)/asmax;
    }
  } while (ran->Get() > weight);
  //msg_Out()<<"   --> final parameters: "<<kt2<<"  "<<z<<endl;
  return;
}

double Collinear_Emission_Generator::GetZ(bool isgluon) {
  double R(ran->Get());
  if (isgluon) {
    double z,wt;
    do {
      z = sqrt(R)/2.;
      wt = 1.-z;
    } while (ran->Get() > wt);
    if (ran->Get() > 0.5) z = 1.-z;
    return z;
  }
  else {
    return 1.-sqrt(1.-R);
  }
}

Particle * Collinear_Emission_Generator::FindSpectator(Particle * part) {
  for (PVIt piter = m_outparts.begin(); piter != m_outparts.end(); ++piter){
    if ( ((*piter)->GetFlow(1) != 0 && (*piter)->GetFlow(1) == part->GetFlow(2)) ||
	 ((*piter)->GetFlow(2) != 0 && (*piter)->GetFlow(2) == part->GetFlow(1))) {
      return *piter;
    }
  }
  //msg_Out()<<METHOD<<": Did not find colour partner in final state, no splitting allowed"<<endl;
  return NULL;
}

bool Collinear_Emission_Generator::
FixKinematics(double kt2, double z, Vec4D & split, Vec4D & gluon , Vec4D & spec){
  Vec4D cms = split+spec;
  if (4.*kt2>cms.Abs2() || split[0] == 0. || spec[0] == 0.) return false;
  //msg_Out()<<METHOD<<"(cms = "<<cms<<"\n"
  //	   <<"from "<<split<<"("<<split.Abs2()<<") + "<<spec<<"("<<spec.Abs2()<<")\n";
  Poincare boost(cms);
  boost.Boost(split);
  boost.Boost(spec);
  Poincare rotat(split,Vec4D(1.,0.,0.,1.));
  rotat.Rotate(split);
  rotat.Rotate(spec);
  //msg_Out()<<split<<" ("<<split.Abs2()<<")+ "<<spec<<" ("<<spec.Abs2()<<")\n";
  
  double phi(ran->Get()*2*M_PI);
  Vec4D kt(0., sqrt(kt2)*cos(phi), sqrt(kt2)*sin(phi), 0.);
  Vec4D pj = (1.-z)*split + kt2/(2.*(1.-z)*split*spec)*spec - kt;
  gluon = (z*split + kt2/(2.*z*split*spec)*spec + kt);
  Vec4D pk((1.-kt2/(2.*z*(1.-z)*split*spec))*spec);
  if (split+spec != pj+gluon+pk ||
      dabs(pj.Abs2())>1.e-6 || dabs(gluon.Abs2())>1.e-6 || dabs(pk.Abs2())>1.e-6){
    msg_Debugging()<<METHOD<<": 4-momentum not conserved:\n"<<endl
		   <<split<<" ("<<split.Abs2()<<") + "<<spec<<" ("<<spec.Abs2()<<")\n"
		   <<" = "<<split+spec<<" vs.\n"
		   <<pj<<" ("<<pj.Abs2()<<") + "<<gluon<<" ("<<gluon.Abs2()<<") + "<<pk<<" ("<<pk.Abs2()<<"\n"
		   <<" = "<<(pj+gluon+pk)<<"\n"
		   <<"for kt2 = "<<kt2<<" and z = "<<z<<".\n";
    THROW(fatal_error,"Four-momentum non-conservation.")
  }
  rotat.RotateBack(pj);
  rotat.RotateBack(gluon);
  rotat.RotateBack(pk);
  boost.BoostBack(pj);
  boost.BoostBack(gluon);
  boost.BoostBack(pk);
  //msg_Out()<<METHOD<<": kinematics works out for kt^2 = "<<kt2<<", z = "<<z<<":\n"
  //	   <<"   "<<split<<" --> "<<pj<<" ("<<pj.Abs2()<<"), \n"
  //	   <<"   "<<spec<<" --> "<<pk<<" ("<<pk.Abs2()<<"), \n"
  //	   <<"   gluon = "<<gluon<<" ("<<gluon.Abs2()<<")\n";
  //pk = split+spec-pj-gluon;
  split = pj;
  spec = pk;
  if (split[0] < 0. || spec[0] < 0. || gluon[0] < 0.) return false;
  return true;
}

Blob * Collinear_Emission_Generator::MakeBlob() {
  Blob * blob = new Blob();
  blob->SetId();
  blob->AddData("Weight",new Blob_Data<double>(1.));
  blob->AddData("Factorisation_Scale",new Blob_Data<double>(1.));
  blob->AddData("Renormalization_Scale",new Blob_Data<double>(1.));
  blob->SetType(btp::Shower);
  blob->SetTypeSpec("ShrimpsCollinearEmissions");
  blob->SetStatus(blob_status::needs_hadronization | blob_status::needs_beams);
  for (size_t i = 0; i < m_beamparts.size(); ++i) {
    Particle * part(new Particle(-1,m_beamparts[i]->Flav(),m_beamparts[i]->Momentum(),
				 m_beamparts[i]->Info()));
    part->SetNumber(0);
    part->SetBeam(m_beamparts[i]->Beam());
    part->SetFlow(1,m_beamparts[i]->GetFlow(1));
    part->SetFlow(2,m_beamparts[i]->GetFlow(2));
    blob->AddToInParticles(part);
    blob->AddToOutParticles(m_beamparts[i]);
  }
  for (size_t i = 0; i < m_inparts.size(); ++i)   blob->AddToInParticles(m_inparts[i]);
  for (size_t i = 0; i < m_outparts.size(); ++i)  blob->AddToOutParticles(m_outparts[i]);
  return blob;
}

void Collinear_Emission_Generator::CleanUp(){
  while (!m_outparts.empty()) { delete m_outparts.back(); m_outparts.pop_back(); }
}
