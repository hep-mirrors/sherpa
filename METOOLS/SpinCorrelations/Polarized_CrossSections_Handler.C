#include "METOOLS/SpinCorrelations/Polarized_CrossSections_Handler.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/SpinCorrelations/PolWeight_Map.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace METOOLS;

Polarized_CrossSections_Handler::Polarized_CrossSections_Handler()
  : m_old_refmom(ATOOLS::Vec4D()), m_new_refmom(ATOOLS::Vec4D()), m_helicitybasis(false),
  m_customweights(std::map<std::string, std::string>())
 {
  // set polarization settings
  InitPolSettings();

  // Determine old and new reference momenta, if constant for all particles (spin basis is not helicity basis)
   InitRefMoms();
}

void Polarized_CrossSections_Handler::InitPolSettings() {
  auto& s = ATOOLS::Settings::GetMainSettings();
  auto pol = s["HARD_DECAYS"]["Pol_Cross_Section"];

  // set settings for polarization definition:
  // spin basis (default is helicity basis according to common definition
  // of polarization in terms of helicity) and reference system (default: Lab)
  // it is possible to specify more than one reference system but only one spin basis per simulation run
  m_spinbasis = pol["Spin_Basis"].SetDefault("Helicity").Get<std::string>();
  m_refsystem = pol["Reference_System"].SetDefault(std::vector<std::string>(1, "Lab"))
    .GetVector<std::string>();

  // read in specified desired custom weights:
  // maximum number of custom weights in 12 (Weight0-Weight10 + Weight)
  // but can be increased by Number_Of_Custom_Weights setting in YAML-File
  // set defaults
  int number_of_custom_weights = pol["Number_Of_Custom_Weights"].SetDefault(10).Get<int>()+1;
  std::vector<std::string> weights(number_of_custom_weights, "no weight");
  pol["Weight"].SetDefault("no weight");
  for (size_t i = 0; i < weights.size(); ++i){
    pol["Weight" + ATOOLS::ToString(i)].SetDefault(weights[i]);
  }
  // only add weights to m_customweights (which transport custom_weights to PolWeightsMap) if there are some one
  if (pol["Weight"].Get<std::string>()!="no weight") {
    m_customweights.emplace("Weight", pol["Weight"].Get<std::string>());
  }
  for (size_t i(0);i<weights.size();++i) {
    if (pol["Weight" + ATOOLS::ToString(i)].Get<std::string>()!="no weight"){
      m_customweights.emplace("Weight" + ATOOLS::ToString(i), pol["Weight" + ATOOLS::ToString(i)].
        Get<std::string>());
    }
  }

  // this setting is only necessary for investigation of identical intermediate particles which decays via different
  // decay channels, and where the different decay channels should specify which of the intermediate particles should be
  // taken as polarized, the decay channel of the polarized particle is set here.
  m_singlepol_channel = pol["Single_Polarized_Channel"].SetDefault("no channel").Get<std::string>();
}

void Polarized_CrossSections_Handler::InitRefMoms() {
  // Determining currently used reference vector from COMIX_DEFAULT_GAUGE-Setting
  static const double invsqrttwo(1.0/sqrt(2.0));
  ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  int n = s["COMIX_DEFAULT_GAUGE"].Get<int>();
  m_old_refmom = ATOOLS::Vec4D(1.0, 0.0, 1.0, 0.0);
  /*
  switch(n) {
    case 1: m_old_refmom=ATOOLS::Vec4D(1.0, 0.0, invsqrttwo, invsqrttwo); break;
    case 2: m_old_refmom=ATOOLS::Vec4D(1.0, invsqrttwo, 0.0, invsqrttwo); break;
    case 3: m_old_refmom=ATOOLS::Vec4D(1.0, invsqrttwo, invsqrttwo, 0.0); break;
  }
  */
  // PREPARATION FOR DEFINITION OF NEW SPINBASIS
  // setting reference vector for defining the spin basis for the polarization vectors in the new Amplitude2_Tensor
  // from user input
  // user input can only be constant reference vectors (not momentum dependent as in case of helicity basis)
  if (m_spinbasis!="Helicity" && m_spinbasis!="ComixDefault"){
    std::string spinbasis_temp(m_spinbasis);
    std::replace(spinbasis_temp.begin(),spinbasis_temp.end(),',',' ');
    auto new_ref_mom = ATOOLS::ToVector<double>(spinbasis_temp);
    if (new_ref_mom.size() > 4 || new_ref_mom.size() < 4){
      THROW(invalid_input, "Reference vector for spinbasis definition must have four components!")
    }
    for (size_t i(0); i<new_ref_mom.size(); ++i){
      m_new_refmom[i] = new_ref_mom[i];
    }
  }
  else if (m_spinbasis=="Helicity"){
    m_helicitybasis = true;
  }
    // setting new reference vector if ComixDefault is chosen for the spin basis
  else if (m_spinbasis=="ComixDefault"){
    m_new_refmom = ATOOLS::Vec4D(m_old_refmom);
  }
}

std::vector<METOOLS::PolWeight_Map*> Polarized_CrossSections_Handler::Treat(ATOOLS::Blob* signalblob,
                                                                            const METOOLS::Amplitude2_Tensor* prod_amps,
                                                                            const std::vector<METOOLS::Decay_Matrix>& decay_matrices)
const{
  if (decay_matrices.empty()){
    THROW(not_implemented, "Polarization for final or initial state particles is not supported yet.")
  }
  int num_particles= prod_amps->NumberParticles();
  if (num_particles!=decay_matrices.size()){
    std::cout << "number of decay matrices in DecayMatrices: " << decay_matrices.size() << std::endl;
    std::cout << "number of particles described by the Amplitude2_Tensor: " << num_particles << std::endl;
    THROW(fatal_error, "DecayMatrices vector does not contain the same number of decay matrices as particles are "
                       "described by the Amplitude2_Tensor pol_amps")
  }
  std::map<int, METOOLS::Polarization_Vector> default_polarization_vectors =
    std::map<int,METOOLS::Polarization_Vector>();
  std::map<int, SpinorType> default_spinors = std::map<int, SpinorType>();
  // Calculate polarization fractions for each desired reference system
  std::vector<METOOLS::PolWeight_Map*> polweights;
  for (size_t h(0); h<m_refsystem.size(); ++h){
    ATOOLS::Vec4D beta = ATOOLS::Vec4D();
    if (m_refsystem[h]!="Lab" && m_refsystem[h]!="RestFrames") {
      beta = Beta(signalblob, prod_amps, m_refsystem[h]);
    }
    polweights.push_back(Calculation(signalblob, prod_amps, decay_matrices, default_polarization_vectors,
                                     default_spinors, beta, m_refsystem[h]));
  }
  return polweights;
}

ATOOLS::Vec4D Polarized_CrossSections_Handler::Beta(const ATOOLS::Blob* signalblob,
                                                    const METOOLS::Amplitude2_Tensor* prod_amps,
                                                    std::string refsystem) const {
  // PREPARATION FOR DEFINITION OF NEW REFERENCE FRAME
  // Polarization is by default defined in laboratory system
  // possible user input keywords: Lab, COM, PPFr, RestFrames
  // Selecting particles in which center of mass frame the polarization vectors should be defined
  std::vector<int> particles;
  ATOOLS::Vec4D beta = ATOOLS::Vec4D(0.0, 0.0, 0.0, 0.0);
  if (refsystem!="Lab" && refsystem!="RestFrames"){
    // TODO: Do only particles have a DecayBlob() which should also be decayed in hard decays?
    // Center-of-Mass frame of the hard decaying particles
    if (refsystem=="COM"){
      for (size_t i(0); i<signalblob->NOutP(); ++i){
        // select hard decaying particles from all outcoming particles of hard process blob
        if (signalblob->ConstOutParticle(i)->DecayBlob()){
          particles.push_back(signalblob->ConstOutParticle(i)->Number());
          beta += (*(signalblob->ConstOutParticle(i)->DecayBlob()))["p_onshell"]->Get<ATOOLS::Vec4D>();
        }
      }
      if (particles.size() != prod_amps->NumberParticles()){
        THROW(fatal_error, "Internal error: Not all particles with hard decay blobs are used to define the COM "
                           "reference frame!")
      }
    }
      // Parton-Parton-frame
    else if (refsystem=="PPFr"){
      for (size_t i(0); i<signalblob->NInP(); ++i){
        particles.push_back(signalblob->ConstInParticle(i)->Number());
        beta += signalblob->ConstInParticle(i)->Momentum();
      }
      if (particles.size() != 2){
        std::cout << particles.size() << " Particles in initial state" << std::endl;
        THROW(fatal_error, "Internal error: More or less than two initial state particles")
      }
    }
      // TODO: Do we need other reference frames than all possible restframes which can be defined from hard process
      //       particles?
      // selecting the center of mass frame defining particles from the user input particle numbers,
      // only particles from hard process possible
      // particle numbers are used because with that particles can be identified without problems, with pdg-codes
      // one would run into problems if one particle is in the initial and the final states, furthermore
      // containers must be handled separately
    else {
      // TODO: Exception if user does not either give a valid key word nor valid particle numbers
      // Syntax is now z.B. 1.0 0.0 0.0 1.0 and not comma separated because "," separate different refsystems
      // only exception is, if only one weight is given, then comma separation is also possible
      std::replace(refsystem.begin(),refsystem.end(),',',' ');
      auto particle_numbers = ATOOLS::ToVector<int>(refsystem);
      for (size_t i(0); i<particle_numbers.size(); ++i){
        if (particle_numbers[i]+1 > signalblob->NInP()+signalblob->NOutP()){
          THROW(invalid_input, "Particle number inputs for defining reference system for polarization "
                               "definition are not valid, particles numbers must be bigger than 0 and "
                               "only particles from hard process can be considered for system definition")
        }
        // +1 due to discrepany between particle numbering in YAML-File and internal particle numbering
        // (former starts at 0, latter at 1)
        if (particle_numbers[i]+1 <= signalblob->NInP()){
          for (size_t j(0); j<signalblob->NInP(); ++j){
            if (signalblob->ConstInParticle(j)->Number() == particle_numbers[i]+1){
              particles.push_back(signalblob->ConstInParticle(j)->Number());
              beta+=signalblob->ConstInParticle(j)->Momentum();
              break;
            }
          }
        }
        else{
          for (size_t j(0); j<signalblob->NOutP(); ++j){
            if (signalblob->ConstOutParticle(j)->Number() == particle_numbers[i]+1){
              particles.push_back(signalblob->ConstOutParticle(j)->Number());
              if (signalblob->ConstOutParticle(j)->DecayBlob()){
                beta += (*(signalblob->ConstOutParticle(j)->DecayBlob()))["p_onshell"]->Get<ATOOLS::Vec4D>();
              }
              else beta += signalblob->ConstOutParticle(j)->Momentum();
              break;
            }
          }
        }
      }
      if (particles.size() != particle_numbers.size()){
        std::cout << "Not all given particle numbers describe particles from the hard process. " <<
                  particle_numbers.size() << " particle numbers are specified in YAML, but only "
                  << particles.size() << " particle numbers are found in hard process! Continue with particles found: "
                  << std::endl;
        for (size_t m(0); m < particles.size(); ++m){
          std::cout << particles[m] - 1 << std::endl;
        }
      }
    }

    if (particles.empty() && !(refsystem == "RestFrames" || refsystem=="Lab")){
      THROW(invalid_input, "Given reference frame " + refsystem + " is not supported/valid.")
    }
  }
  return beta;
}

PolWeight_Map* Polarized_CrossSections_Handler::Calculation(ATOOLS::Blob* signalblob,
                                                  const METOOLS::Amplitude2_Tensor* prod_amps,
                                                  const std::vector<METOOLS::Decay_Matrix>& decay_matrices,
                                                  std::map<int, METOOLS::Polarization_Vector>& default_polarization_vectors,
                                                  std::map<int, SpinorType>& default_spinors,
                                                  ATOOLS::Vec4D beta, std::string refsystem) const {

  METOOLS::Amplitude2_Tensor* pol_amps = new METOOLS::Amplitude2_Tensor(*prod_amps);
  std::vector<METOOLS::Decay_Matrix> trafo_decay_matrices = std::vector<METOOLS::Decay_Matrix>();
  ATOOLS::Vec4D new_ref_mom(m_new_refmom);

  if (!(m_spinbasis=="ComixDefault" && refsystem == "Lab")) {
    // Determination of transformation coefficients and transformation of decay matrices
    METOOLS::Amplitude2_Tensor *tmp_amps(pol_amps);

    // -- CALCULATION OF DEFAULT AND DESIRED POLARIZATION OBJECTS --
    // contains the transformation coefficient vectors for the several particles in the amplitude2_tensor
    std::vector<std::vector<std::vector<Complex>>> coeff_vec;
    std::vector<std::vector<std::vector<Complex>>> conj_coeff_vec;
    do {
      if (tmp_amps->SpinDegreesofFreedom() != 3) {
        THROW(fatal_error, "basis transformation for polarization definition is currently only implemented for "
                           "massive vector bosons")
      }

      // Determination of necessary particle information
      ATOOLS::Vec4D mom(tmp_amps->CurrentParticle().Momentum());
      // is current particle a final state particle?
      bool out(false);
      bool found(false);
      for (size_t i(0); i < signalblob->NOutP(); ++i) {
        // TODO: Are there cases where p_onshell does not exist or where its use here is wrong?
        if (*(signalblob->OutParticle(i)) == tmp_amps->CurrentParticle()){
          out = true;
          found = true;
          // use on shell particle momenta for the calculation of the default and desired polarization objects for
          // intermediate particles
          if (signalblob->ConstOutParticle(i)->DecayBlob()) {
            mom = (*(signalblob->ConstOutParticle(i)->DecayBlob()))["p_onshell"]->Get<ATOOLS::Vec4D>();
          }
        }
      }
      if (!found)
        THROW(fatal_error, "A particle in Amplitude2_tensor is not an outgoing particle of the signal blob.")
      bool anti = tmp_amps->CurrentParticle().RefFlav().IsAnti();
      double spin = tmp_amps->CurrentParticle().RefFlav().Spin();
      // from comparison with literature: seems that W- has the same polarization vectors than W+
      // (not complex conjugate)
      // TODO: Check with Stefan (Comix) and Giovanni
      if (tmp_amps->CurrentParticle().RefFlav().IDName()=="W-"){
        anti = false;
      }
      // TODO: particles with mom.PSpat()=0 in laboratory frame?
      ATOOLS::Vec4D direction = ATOOLS::Vec4D(1., ATOOLS::Vec3D(mom) / mom.PSpat());

      // CALCULATION OF DEFAULT POLARIZATION OBJECTS, if this is the first call of this function
      if (spin == 1) {
        if (default_polarization_vectors.find(tmp_amps->CurrentParticle().Number()) ==
        default_polarization_vectors.end()) {
          default_polarization_vectors.emplace(tmp_amps->CurrentParticle().Number(),
                                                 METOOLS::Polarization_Vector(mom, m_old_refmom));
        }
      }

      // for refsystem=RestFrames polarization of a particle is defined in its own rest frame, therefore beta is
      // different for each particle
      if (refsystem=="RestFrames") {
        beta = ATOOLS::Vec4D(mom);
      }
      // Boosting momentum to the desired reference frame
      ATOOLS::Poincare momboost(beta);
      if (beta != ATOOLS::Vec4D(0.0, 0.0, 0.0, 0.0)) {
        momboost.Boost(mom);
      }

      // DETERMINATION OF NEW REFERENCE VECTOR
      // For getting physical polarisation in helicity basis a special reference vector must be used
      // Reference vector according to Alnefjord et al. 2021
      if (m_helicitybasis) {
        // more general if-condition than only regarding RestFrames refsystem, since if user chooses
        // the center of mass of one intermediate particle as reference frame the same problems in defining a spin
        // axis occur compared to the RestFrames reference system
        // in this case, the spin axis is defined as the particle's direction of flight in the laboratory frame too
        if (mom.PSpat() < 1e-12) new_ref_mom = ATOOLS::Vec4D(direction[0], -direction[1], -direction[2], -direction[3]);
        else new_ref_mom = ATOOLS::Vec4D(1., -ATOOLS::Vec3D(mom) / mom.PSpat());
      }
      if (new_ref_mom.Abs2() > 1e-8) {
        THROW(fatal_error, "new reference vector has to be light-like")
      }

      // CALCULATION OF POLARIZATION OBJECTS IN THE NEW POLARIZATION BASIS
      std::vector<std::vector<Complex>> coeff_in;
      std::vector<std::vector<Complex>> coeff_out;
      if (spin==1){
        // calculate polarisation vectors with new reference momentum
        METOOLS::Polarization_Vector new_polarisation(mom, new_ref_mom);
        // transformation back to laboratory system where matrix elements are calculated
        if (beta != ATOOLS::Vec4D(0.0, 0.0, 0.0, 0.0)) {
          for (size_t i(0); i < 3; ++i) {
            momboost.BoostBack(new_polarisation[i]);
          }
          momboost.BoostBack(mom);
        }
        // -- DETERMINATION OF THE TRANSFORMATION COEFFICIENTS BETWEEN THE DEFAULT AND THE DESIRED POLARIZATION
        // DEFINITION ---
        // determine transformation coefficients under consideration that constructor of polarization vector only
        // generate ingoing polarization vectors for particles, for outgoing ones or antiparticles (or equivalently the
        // complex conjugate ones in the squared amplitude tensor) the transformation coefficients should be the complex
        // conjugate of the ingoing coefficients
        coeff_in = new_polarisation.BasisTrafo(
            default_polarization_vectors.find(tmp_amps->CurrentParticle().Number())->second);
        coeff_out = std::vector<std::vector<Complex>>(coeff_in);
        // determine complex conjugate coefficients
        for (size_t i(0); i < coeff_in.size(); ++i) {
          for (size_t j(0); j < coeff_in[i].size(); ++j) {
            coeff_out[i][j] = conj(coeff_in[i][j]);
          }
        }
        if ((out && !anti) || (!out && anti)) {
          coeff_vec.push_back(coeff_out);
          conj_coeff_vec.push_back(coeff_in);
        } else {
          coeff_vec.push_back(coeff_in);
          conj_coeff_vec.push_back(coeff_out);
        }
      }
      // --- TRANSFORMATION OF THE DECAY MATRICES ---
      for (size_t k(0); k<decay_matrices.size(); ++k){
        if (decay_matrices[k].Particle()->Number()==tmp_amps->CurrentParticle().Number()){
          METOOLS::Decay_Matrix decay_matrix(decay_matrices[k]);
          if (anti) decay_matrix.PolBasisTrafo(coeff_out, coeff_in);
          else decay_matrix.PolBasisTrafo(coeff_in, coeff_out);
          trafo_decay_matrices.push_back(decay_matrix);
          break;
        }
      }
      tmp_amps = tmp_amps->Next()[0];
    }
    while (tmp_amps->IsP_Next());
    // Transformation of production tensor
    pol_amps->PolBasisTrafo(coeff_vec, conj_coeff_vec, 0);
  }
  else trafo_decay_matrices = std::vector<METOOLS::Decay_Matrix>(decay_matrices);
  if (trafo_decay_matrices.size()!=decay_matrices.size()){
    THROW(fatal_error, "Vector of decay matrices does not contain the same hard decaying particles as the given "
                       "production tensor!")
  }

  // Multiplication of production tensor and decay matrices to get the actual polarized matrix elements
  if (!decay_matrices.empty()){
    for (size_t k(0); k<trafo_decay_matrices.size(); ++k){
      METOOLS::Decay_Matrix* pointer_to_Decay_Matrix = new METOOLS::Decay_Matrix(trafo_decay_matrices[k]);
      pol_amps->Multiply(pointer_to_Decay_Matrix);
      delete pointer_to_Decay_Matrix;
    }
  }
  // Labeling the polarized matrix elements and storing results in p_polweights
  PolWeight_Map* polWeightMap = new METOOLS::PolWeight_Map(pol_amps,m_customweights,
                                                           m_singlepol_channel);
  Tests((*signalblob)["ATensor"]->Get<METOOLS::Amplitude2_Tensor*>(), pol_amps);
  delete pol_amps;
  return polWeightMap;
}

void Polarized_CrossSections_Handler::Tests(const METOOLS::Amplitude2_Tensor* amps,
                                            const METOOLS::Amplitude2_Tensor* trafo_pol_amps) const {
  // test, whether unpolarized cross-section determined by spin correlation algorithm is identical to the
  // unpolarized cross-section one receives after changing of bases for polarization definition and calculating of
  // polarized cross-sections
  METOOLS::Amplitude2_Tensor* tmpTensor = new METOOLS::Amplitude2_Tensor(*amps);
  METOOLS::Amplitude2_Tensor* tmpTrafoTensor = new METOOLS::Amplitude2_Tensor(*trafo_pol_amps);
  Complex unpol = tmpTensor->Sum();
  Complex trafo_unpol = tmpTrafoTensor->Sum();
  if ((trafo_unpol-unpol).real()>ATOOLS::dabs(trafo_unpol.real())*1e-8 || ATOOLS::dabs(trafo_unpol.imag())>1e-8 ||
      ATOOLS::dabs(unpol.imag())>1e-8){
    std::cout<<"Polarization_Warning in"<< METHOD <<
               "Testing consistency between unpolarized cross section before and after transformation"
               " to another bases failed..." << std::endl;
    msg_Out() << "Unpolarized cross section resulting after spin correlation algorithm: "
              << std::setprecision(20) << unpol << std::endl;
    msg_Out() << "Unpolarized cross section after transformation:: " << std::setprecision(20) << trafo_unpol
              << std::endl;
  }
  delete tmpTensor;
  delete tmpTrafoTensor;
}

