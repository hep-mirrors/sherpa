#include "REMNANTS/Main/Hadron_Remnant.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;

Hadron_Remnant::Hadron_Remnant(PDF::PDF_Base * pdf,const unsigned int beam):
  Remnant_Base(rtp::hadron,beam),
  p_pdf(pdf), p_partons(&(p_pdf->Partons())), m_beamflav(pdf->Bunch()),
  p_valence(nullptr), p_remnant(nullptr), p_recoiler(nullptr), p_spectator(nullptr), m_ff(Form_Factor()),
  m_valence(false), m_alpha(0.), m_gamma(1.), m_beta(-1.5),  m_invb(1./(m_beta+1)), m_LambdaQCD(0.25)
{
  m_scale2 = Max(4.0,p_pdf->Q2Min());
  ConstructConstituentFlavours();
}

void Hadron_Remnant::ConstructConstituentFlavours() {
  if (m_constituents.size()>0) return;
  int hadint=(m_beamflav.Kfcode()-(m_beamflav.Kfcode()/10000)*10000)/10;
  if ((hadint>100)&&(hadint<1000)) {
    m_constituents.push_back(Flavour((kf_code)(hadint)/100));
    m_constituents.push_back(Flavour((kf_code)((hadint-(hadint/100)*100)/10)));
    m_constituents.push_back(Flavour((kf_code)(hadint-(hadint/10)*10)));
  }
  else if ((hadint>10)&&(hadint<100)) {
    m_constituents.push_back(Flavour((kf_code)(hadint)/10));
    m_constituents.push_back(Flavour((kf_code)(hadint-(hadint/10)*10)));
  }
  else THROW(critical_error,"Cannot determine constituents.");
  if (m_beamflav.IsAnti()) {
    for(auto& flit : m_constituents) flit = flit.Bar();
  }
}

bool Hadron_Remnant::IsValence(Particle * part) {
  // only one valence parton.
  if (m_valence) return false;
  // assume valence(q) = pdf(q) - pdf(qbar)
  Flavour flav = part->Flav();
  for (const auto& flit : m_constituents) {
    if (flav==flit) {
      Vec4D   mom  = part->Momentum();
      double x = mom[0]/m_residualE;
      p_pdf->Calculate(x,sqr(flav.Mass())+m_scale2);
      double val = p_pdf->GetXPDF(flav)-p_pdf->GetXPDF(flav.Bar());
      double tot = p_pdf->GetXPDF(flav);
      m_valence = (val/tot > ran->Get());
      if (m_valence) p_valence = part;
      return m_valence;
    }
  }
  return false;
}

void Hadron_Remnant::MakeSpectator(Particle * parton) {
  // If a shower initiator is a sea-quark or antiquark, a corresponding
  // antiflavour has to be added to the spectators.
  p_spectator = nullptr;
  if (IsValence(parton)) return;
  Flavour flav = parton->Flav();
  if (!flav.IsQuark()) return;
  p_spectator = MakeParticle(flav.Bar());
  p_spectator->SetFlow((flav.Bar().IsAnti()?2:1),-1);
  p_spectator->SetPosition(parton->XProd());
  p_colours->AddColour(m_beam,(flav.Bar().IsAnti()?1:0),p_spectator);
  m_spectators.push_front(p_spectator);
}

Particle * Hadron_Remnant::MakeParticle(const Flavour & flav) {
  Particle * part = new Particle(-1,flav,Vec4D(0.,0.,0.,0.),'B');
  part->SetNumber();
  part->SetBeam(m_beam);
  part->SetPosition(m_position+m_ff());
  return part;
}

bool Hadron_Remnant::FillBlob(ParticleMomMap *ktmap,const bool & copy) {
  // Add remnants, diquark and quark, if necessary.
  if (!p_valence || !p_remnant) MakeRemnants();
  // Possibly adjust final pending colours with extra gluons - in prinicple one may have
  // to check that they are not singlets ....
  CompensateColours();
  msg_Debugging() << METHOD << ": Filling blob with remnants, extracted = "
                  << m_extracted << ", \n and spectators = " << m_spectators
                  << "\n";
  // Assume all remnant bases already produced a beam blob = p_beamblob
  MakeLongitudinalMomenta(ktmap,copy);
  bool colourconserved = p_beamblob->CheckColour(true);
  if (!colourconserved) {
    msg_Error()<<"Error in "<<METHOD<<" for \n"<<(*p_beamblob)<<"\n";
    p_colours->Output();
    return false;
  }
  return true;
}

void Hadron_Remnant::CompensateColours() {
  while (!p_colours->Colours(m_beam,0).empty() &&
         !p_colours->Colours(m_beam,1).empty() &&
	 p_colours->Colours(m_beam,0)!=p_colours->Colours(m_beam,1)) {
    Particle * gluon = MakeParticle(Flavour(kf_gluon));
    for (size_t i=0;i<2;i++) gluon->SetFlow(i+1,p_colours->NextColour(m_beam,i));
    m_spectators.push_back(gluon);
  }
}

bool Hadron_Remnant::MakeRemnants() {
  // If no valence quark has been extracted to date, a quark-diquark
  // pair must be constructed.  the idea is to pick one of the three flavours
  // at random for the quark, add it to the spectators, then construct the
  // "conjugate" diquark and add it as well to the spectators
  Flavour valflav;
  size_t  index;
  if (!p_valence) {
    int random = int(ran->Get()*m_constituents.size());
    FlavourList::iterator flit=m_constituents.begin();
    for (size_t i=0;i<random;i++) flit++;
    valflav    = *flit;
    p_valence  = MakeParticle(valflav);
    index      = ((valflav.IsQuark() && !valflav.IsAnti()) ||
		 (valflav.IsDiQuark() && valflav.IsAnti()))?0:1;
    p_valence->SetFlow(index+1,p_colours->NextColour(m_beam,index));
    m_spectators.push_back(p_valence);
  }
  else {
    valflav = p_valence->Flav();
    index      = ((valflav.IsQuark() && !valflav.IsAnti()) ||
		  (valflav.IsDiQuark() && valflav.IsAnti()))?0:1;
  }
  p_remnant    = p_recoiler = MakeParticle(RemnantFlavour(valflav));
  p_remnant->SetFlow(2-index,p_colours->NextColour(m_beam,1-index));
  m_spectators.push_front(p_recoiler);
  return true;
}

Flavour Hadron_Remnant::RemnantFlavour(const Flavour & flav) {
  // Counter taken to make sure only two flavours are used
  // to construct diquark - either qq'_0 or qq_1.
  bool taken = false;
  std::vector<int> kfs;
  for (const auto& flit : m_constituents) {
    if (taken && flav==flit) continue;
    kfs.push_back(((flit.IsAnti() && !m_beamflav.IsAnti())?-1:1)*flit.Kfcode());
    taken = true;
  }
  int kfcode = 1 + (kfs.size()==2 && kfs[0]==kfs[1]?2:0);
  for (size_t i=0;i<kfs.size();i++) kfcode += kfs[i]*pow(10,kfs.size()+1-i);
  return m_beamflav.IsAnti()?Flavour(kfcode).Bar():Flavour(kfcode);
}

void Hadron_Remnant::MakeLongitudinalMomenta(ParticleMomMap *ktmap,const bool & copy) {
  // Calculate the total momentum that so far has been extracted through
  // the shower initiators and use it to determine the still available
  // momentum; the latter will be successively reduced until the
  // rest is taken by the diquark.
  Vec4D availMom = p_beam->OutMomentum();
  for (auto pmit : m_extracted) {
    availMom -= pmit->Momentum();
    if (copy) {
      auto pcopy = new Particle(*pmit);
      pcopy->SetNumber();
      pcopy->SetBeam(m_beam);
      p_beamblob->AddToOutParticles(pcopy);
    }
    else p_beamblob->AddToOutParticles(pmit);
    (*ktmap)[pmit] = Vec4D();
  }
  msg_Debugging() << METHOD << ": Longitudinal momentum left for remnants = " << availMom
            << "\n";
  double remnant_masses = 0.;
  for (Particle  const * pit : m_spectators) {
    remnant_masses += Max(pit->Flav().HadMass(), m_LambdaQCD);
  }
  if (remnant_masses > m_residualE)
    msg_Debugging() << METHOD << ": Warning, HadMasses of remnants = "
                << remnant_masses << " vs. residual energy = " << m_residualE << "\n";
  for (auto part : m_spectators) {
    if (availMom[0] < 0) msg_Error() << METHOD << ": Negative Energy in Remnants! \n";
    if (part==m_spectators.back()) part->SetMomentum(availMom);
    else {
      part->SetMomentum(SelectZ(part->Flav(),availMom[0], remnant_masses)*availMom);
      availMom -= part->Momentum();
      remnant_masses -= Max(part->Flav().HadMass(), m_LambdaQCD);
    }
    msg_Debugging() << METHOD << ": set momentum for "<<part->Flav()<<" to "
              << part->Momentum() << "\n";
    if (copy) {
      auto pcopy = new Particle(*part);
      pcopy->SetNumber();
      pcopy->SetBeam(m_beam);
      p_beamblob->AddToOutParticles(pcopy);
    }
    else p_beamblob->AddToOutParticles(part);
    (*ktmap)[part] = Vec4D();
  }
}

double Hadron_Remnant::SelectZ(const ATOOLS::Flavour &flav, double restmom,
                               double remnant_masses) const {
  double zmin = Max(flav.HadMass(), m_LambdaQCD) / restmom;
  double zmax = zmin + (restmom - remnant_masses) / restmom;
  double z;
  if (zmax < zmin) {
    msg_Debugging() << METHOD << ": Error, zmin, zmax = " << zmin <<", "<<zmax << "\n";
    return 0;
  }
  // This assumes that there is exactly one diquark in the spectators, which
  // will be the valence quark, taking most of the remaining momentum
  if (!flav.IsDiQuark()) {
    // Assume functional from of z^beta with beta = -1.5 (default)
    // Maybe beta_gluon != beta_quark, but leave it for the time being
    if (m_beta!=-1) { 
      double rand = ran->Get();
      z = pow(rand*pow(zmax,m_beta+1.)+(1.-rand)*pow(zmin,m_beta+1.),m_invb);
    }
    else
      z = zmin * pow(zmax/zmin,ran->Get());
  }
  else {
    // If di-quark assume form peaking at 1, something like
    // exp(-gamma/z)*(1-z)^alpha -> realised by hit-or-miss
    double wtmax = pow((1.-zmin),m_alpha)*exp(-m_gamma/zmax);
    double wt;
    do {
      z  = zmin + (zmax-zmin)*ran->Get();
      wt = pow((1.-z),m_alpha)*exp(-m_gamma/z);
    } while (wt < wtmax*ran->Get());
  }
  return z;
}

void Hadron_Remnant::Reset(const bool & resc,const bool & DIS) {
  Remnant_Base::Reset();
  while (!m_spectators.empty()) {
    Particle * part = m_spectators.front();
    if (part->ProductionBlob())
      part->ProductionBlob()->RemoveOutParticle(part);
    if (part->DecayBlob())
      part->DecayBlob()->RemoveInParticle(part);
    delete part;
    m_spectators.pop_front();
  }
  /////// TODO: Have to fix this!!!!!
  m_spectators.clear();
  if (resc)
    msg_Out()<<METHOD<<"(resc = "<<resc<<"): "
	     <<p_beam->InMomentum()<<" - "<<p_beam->OutMomentum()<<"\n";
  m_residualE = ( resc ? (p_beam->InMomentum()-p_beam->OutMomentum())[0] :
		  p_beam->OutMomentum()[0] );
  m_valence   = false;
  p_valence   = p_remnant = p_recoiler = nullptr; 
}

bool Hadron_Remnant::TestExtract(const Flavour &flav,const Vec4D &mom) {
  DEBUG_FUNC("");
  // Is flavour element of flavours allowed by PDF?
  if (p_partons->find(flav)==p_partons->end()) {
    msg_Error()<<METHOD<<": flavour "<<flav<<" not found.\n";
    return false;
  }
  if (mom[0] < flav.HadMass()) {
    msg_Debugging() << METHOD << ": parton too soft, mass = " << flav.HadMass()
                    << " and energy = " << mom[0] << "\n";
    return false;
  }
  // Still enough energy for parton and its remnant quark?
  // We're checking for the remnant masses as well, so the stretching onto the
  // mass-shell works out later.
  double required_energy = EstimateRequiredEnergy()
      + mom[0] + Max(flav.HadMass(), m_LambdaQCD);
  if (required_energy>m_residualE) {
      msg_Debugging()<<METHOD<<": not enough energy for remnants, E_{required} "
                <<mom[0]<<" > E_{in beam} = "<<m_residualE<<".\n";
    return false;
  }
  // Still enough energy?  And in range?
  double x = mom[0]/m_residualE;
  if (x<p_pdf->XMin() || x>p_pdf->XMax()) {
    msg_Error()<<METHOD<<": out of limits, x = "<<x<<".\n";
    return false;
  }
  msg_Debugging()<<flav<<" with mom = "<<mom<<" can be extracted.\n";
  return true;
}

void Hadron_Remnant::Output() {
  msg_Out()<<METHOD<<"("<<m_beam<<", "<<m_beamflav<<").\n"
	   <<"   Constituents are [ ";
  for (const auto& flit : m_constituents)  msg_Out()<<flit<<" ";
  msg_Out()<<"]\n"
	   <<"   Partons are { ";
  for (const auto& flit : *p_partons) {
    msg_Out()<<" "<<flit;
  }
  msg_Out()<<"}.\n";
}

double Hadron_Remnant::EstimateRequiredEnergy() const
{
  double masses = 0.;
  for (Particle const * pit : m_spectators) {
    masses += Max(pit->Flav().HadMass(), m_LambdaQCD);
  }
  if (!m_valence) {
    masses += 2 * Flavour(kf_d).HadMass();
  }
  // Adding lambda_QCD for the gluon that might be added later in CompensateColours()
  return masses + m_LambdaQCD;
}
