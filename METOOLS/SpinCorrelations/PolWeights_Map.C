#include "PolWeights_Map.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHASIC++/Decays/Decay_Channel.H"
#include "ATOOLS/Org/Message.H"


#include <cmath>
#include <iostream>
#include <utility>
#include <cstring>

using namespace METOOLS;

PolWeights_Map::PolWeights_Map():
  m_unpolcrosssec(Complex(0., 0.)), m_massive_vb(true), m_custom_weights(std::map<std::string, std::string>()),
  m_singlepol_channel("no channel"), m_trans_mode(1), m_interference_weights(std::set<std::string>()),
  p_all_weights(NULL), m_pol_checks(false) {
}

PolWeights_Map::PolWeights_Map(const METOOLS::Amplitude2_Tensor* amps, int trans_mode,
                               std::map<std::string, std::string> custom_weights, std::string singlepol_channel,
                               bool pol_checks)
  : m_custom_weights(custom_weights), m_trans_mode(trans_mode), m_singlepol_channel(singlepol_channel),
    m_interference_weights(std::set<std::string>()), m_pol_checks(pol_checks){
  // set attribute values
  // unpolarized cross section
  METOOLS::Amplitude2_Tensor* tmpTensor = new METOOLS::Amplitude2_Tensor(*amps);
  m_unpolcrosssec = tmpTensor->Sum();
  // Test, whether unpolarized result is real
  if (m_pol_checks && m_unpolcrosssec.imag() > 1e-8){
    std::cout<<"Polarization_Warning in "<< METHOD <<
             ": unpolarized result is not real" << std::endl;
    msg_Out() << "imaginary part of the unpolarized result: " << std::endl;
    msg_Out() << m_unpolcrosssec.imag() << std::endl;
  }
  delete tmpTensor;

  m_massive_vb = false;
  p_all_weights = new PolWeights_Map();

  // FILL MAP
  // Add all basic polarization labels
  LabelAndSeparate(amps, "start");
  if (m_pol_checks){
    Tests();
    erase("polsum");
  }

  // First calculation of partially unpolarized weights and adding to polweights map to enable transverse calculation of
  // base AND partially unpolarized weights
  std::vector<std::string> finished_weights = Unpol(amps);

  // Calculation of single-polarized channel weights for identical particles -> private feature, not in manual yet
  // Add single channel weights
  if (m_singlepol_channel != "no channel") AddSinglePolWeights(amps);

  // Calculation of transverse weights for massive vector bosons
  if (m_massive_vb) {
    if (m_trans_mode == 0 || m_trans_mode == 2) Transverse(false);
    if (m_trans_mode == 1 || m_trans_mode == 2) Transverse(true);
    if (m_trans_mode > 2) THROW(not_implemented, "Given coherent_weights_mode not implemented")
  }

  // Add user specified custom weights
  if (!m_custom_weights.empty()) AddCustomWeights(amps, finished_weights);
}

PolWeights_Map::~PolWeights_Map() {
  delete p_all_weights;
}

void PolWeights_Map::LabelAndSeparate(const METOOLS::Amplitude2_Tensor* amps, const std::string& mode,
                                      const std::string& prefix, bool nonzero_weight, std::string spin_label){
  if (mode != "start" && mode != "pol" && mode != "int") THROW(fatal_error, "Invalid mode for PolWeights_Map::Init")

  std::vector<Amplitude2_Tensor*> next_amps;
  if (amps->IsP_Next()){
    next_amps = amps->Next();
    int m_nhel = std::sqrt(next_amps.size());
    if (m_nhel>3) THROW(not_implemented, "Particles with spin bigger than 1 are not implemented for polarized cross "
                                         "sections yet")
    if (m_nhel == 3) m_massive_vb = true;

    // list of strings for the different possible polarization
    std::vector<std::string> spin_strings;
    // TODO: LABELING OF TRANSVERSE POLARIZED MATRIXELEMENTS IS SWITCHED HERE TO GET THE RIGHT LABELING IN THE EVENT
    //       OUTPUT UNTIL SWITCHED ORDERING ISSUE IN MATRIXELEMENT GENERATORS IS FIXED,
    //       POSSIBLE TESTS: DECAY ANGLE OF + -  AND - - POLARIZED VECTOR BOSONS IN VECTOR BOSON PRODUCTION PROCESSES
    spin_strings.push_back("-");
    spin_strings.push_back("+");
    if (m_nhel == 3) spin_strings.push_back("0");

    // recursive generation of the spin label as well as separation between on- and off-diagonal entries
    // (new_mode = "pol": terms where all particles are in a definite polarization state or new_mode="int":
    // terms describing interference between different polarizations)
    // spin labels have the form:
    // <particle1>.<polarization_index_in_matrix_element><polarization_index_in_complex_conjugate_matrix_element>_
    // <particle2>.<polarization_index_in_matrix_element><polarization_index_in_complex_conjugate_matrix_element> ...
    // order of the intermediate particles in the label is according to the order of the particles in the
    // Amplitude2_Tensor
    for (size_t i(0); i<m_nhel*m_nhel; ++i) {
      std::string new_mode;
      if (i % (m_nhel + 1) == 0 && mode!="int") new_mode = "pol";
      else new_mode = "int";

      if (mode=="start") LabelAndSeparate(next_amps[i], new_mode, prefix, nonzero_weight,
                                          spin_label + amps->CurrentParticle().RefFlav().IDName() + "."
                                          + spin_strings[i % m_nhel] + spin_strings[i / m_nhel]);
      else LabelAndSeparate(next_amps[i], new_mode, prefix, nonzero_weight,
                            spin_label + "_" + amps->CurrentParticle().RefFlav().IDName() + "."
                            + spin_strings[i % m_nhel] + spin_strings[i / m_nhel]);
    }
  }
  else{
    // add founded polarization fraction to the map, key = spin label
    p_all_weights->emplace(prefix + spin_label, double(nonzero_weight) * (amps->Value() / m_unpolcrosssec));
    if (mode=="pol") {
      // add polarization fraction with shortened spin label for polarized contributions (i.e. only one polarization
      // index per particle)
      emplace(prefix + ShortName(spin_label), double(nonzero_weight) * (amps->Value() / m_unpolcrosssec));
      // calculate sum of polarizations for later tests: look whether it already exists in the map; if not,
      // generate a new one with the current polarization fraction; if yes, find the corresponding value, add
      // the current fraction and write it to the map
      if (m_pol_checks){
        auto it = find(prefix + "polsum");
        if (it==end()) emplace(prefix + "polsum", double(nonzero_weight) * (amps->Value() / m_unpolcrosssec));
        else{
          Complex tmp = it->second;
          tmp += amps->Value() / m_unpolcrosssec;
          (*this)[prefix + "polsum"] = tmp;
        }
      }
    }
      // add interference term to the map, if it does not exist, otherwise adding the current interference term to
      // one already existing in the map
    else if (mode=="int"){
      m_interference_weights.emplace(prefix+spin_label);
      auto it1 = find(prefix+"int");
      if (it1==end()) emplace(prefix+"int", double(nonzero_weight)*(amps->Value() / m_unpolcrosssec));
      else{
        Complex tmp = it1->second;
        tmp += double(nonzero_weight)*(amps->Value() / m_unpolcrosssec);
        (*this)[prefix + "int"] = tmp;
      }
    }
    // for Amplitude2_Tensor with only one entry (all intermediate particles considered as unpolarized)
    else if (mode=="start" && amps->Value()!=Complex(-1,0))
      emplace(prefix, double(nonzero_weight)*(amps->Value() / m_unpolcrosssec));
    else THROW(fatal_error, "No Tensor")
  }
}

std::set<std::string> PolWeights_Map::ListofKeys() const{
  std::set<std::string> keys;
  for (auto const& element : *this) {
    keys.emplace(element.first);
  }
  return keys;
}

std::string PolWeights_Map::ShortName(std::string name) const{
  std::string prefix;
  std::string temp_name(name);
  std::replace(temp_name.begin(),temp_name.end(),'_',' ');
  auto label_parts = ATOOLS::ToVector<std::string>(temp_name);

  // determine possible existing prefix
  if (label_parts[0]=="dc" || label_parts[0].substr(0,6)=="Weight"){
    prefix = label_parts[0] + "_";
    label_parts.erase(label_parts.begin());
  }
  std::string new_spin_label(prefix);

  // replace all doubly appearing polarization indices by one single index
  for (size_t k(0); k<label_parts.size(); ++k){
    // for partially unpolarized particles, coint, int, polsum are not appearing after the first underscore
    if (label_parts[k] == "coint" || label_parts[k] == "int" || label_parts[k] == "polsum") return name;
    std::replace(label_parts[k].begin(),label_parts[k].end(),'.',' ');
    auto label_parts2 = ATOOLS::ToVector<std::string>(label_parts[k]);
      if (label_parts2[1][0]==label_parts2[1][1]) new_spin_label += label_parts2[0] +  "." + label_parts2[1][0];
      else new_spin_label += label_parts2[0] +  "." + label_parts2[1];
    if (k<label_parts.size()-1) new_spin_label += "_";
  }
  return new_spin_label;
}

std::set<std::string> PolWeights_Map::TransverseKeys(std::set<std::string> keys, int level, bool coherent) const{
  std::set<std::string> tmp_keys;
  int num_particles(level);

  // look at each key in the key list and check whether the level particle is a vector boson with + or - helicity
  // if yes, generate the corresponding transverse polarized label (where +/- -> T/t)
  while (!keys.empty()){
    std::string prefix;
    auto it = keys.begin();
    std::string temp_label = (*it);
    std::replace(temp_label.begin(),temp_label.end(),'_',' ');
    auto label_parts = ATOOLS::ToVector<std::string>(temp_label);
    // determine possible prefix
    if (label_parts[0] == "dc" || label_parts[0].substr(0, 6) == "Weight") {
      prefix = label_parts[0] + "_";
      label_parts.erase(label_parts.begin());
    }

    num_particles = label_parts.size();
    std::replace(label_parts[level].begin(), label_parts[level].end(),'.',' ');
    auto parts = ATOOLS::ToVector<std::string>(label_parts[level]);
    std::string new_string(prefix);

    if ((parts[0]=="W+" || parts[0]=="W-" || parts[0]=="Z") && (parts[1]=="++" || parts[1]=="--")){
      // generate the spin label for the new transverse weight
      // string of particles at smaller levels
      for (size_t j(0); j<level; ++j) {
        new_string += label_parts[j] + "_";
      }
      // string of particle at the current level
      new_string += parts[0];
      if (coherent) new_string += ".T";
      else new_string += ".t";
      // if not the last particle is currently considered the string of the particles with higher level
      // will be added to the new spin label
      if (level+1 < label_parts.size()) new_string += "_";
      for (size_t j(level+1); j<label_parts.size(); ++j){
        new_string += label_parts[j];
      }
    }
      // add weights where the current particle is unpolarized or in a definite polarization state (longitudinally
      // polarized massive VB or different particle species) since particles at higher levels can be transverse
      // polarized VBs
      // with that, the interference terms in p_all_weights are excluded
    else if (parts[1].size()==1 || parts[1][0]==parts[1][1]) new_string = (*it);
    if (!new_string.empty() && new_string!=prefix) tmp_keys.emplace(new_string);
    keys.erase(it);
  }
  return level+1 != num_particles ? TransverseKeys(tmp_keys, level + 1, coherent) : tmp_keys;
}

std::vector<std::string> PolWeights_Map::ExpandLabels(const std::vector<std::string>& transverse_labels, int level, int num_particles) const {
  std::vector<std::string> tmp_keys;
  if (level==0 && transverse_labels.size() != 1) THROW(fatal_error, "PolWeights_Map::ExpandLabels() can only determine the "
                                                        "polarization fractions to add for one single transverse "
                                                        "polarization combination per call!");
  for (const auto & current_string : transverse_labels) {
    std::string prefix;
    std::string temp_label = current_string;
    std::replace(temp_label.begin(), temp_label.end(), '_', ' ');
    auto label_parts = ATOOLS::ToVector<std::string>(temp_label);

    // determine spin label prefix
    if (label_parts[0] == "dc" || label_parts[0].substr(0, 6) == "Weight") {
      prefix = label_parts[0] + "_";
      label_parts.erase(label_parts.begin());
    }

    if (level == 0) num_particles = label_parts.size();
    std::replace(label_parts[level].begin(), label_parts[level].end(), '.', ' ');
    auto parts = ATOOLS::ToVector<std::string>(label_parts[level]);

    // spin labels of polarization fractions to add are generated recursively by iterating through the particles
    // described by the spin labels, at each call, the T or t of the particle at level level in the spin label is
    // replace by ++, -- and for the coherent transverse signal definition also: +-,-+
    // next call of the method than runs on the intermediate spin labels to add obtained in the call before such that
    // each +, - combination is covered
    if ((parts[0] == "W+" || parts[0] == "W-" || parts[0] == "Z") && (parts[1] == "T" || parts[1] == "t")) {
      std::vector<std::string> expanded_labels;
      int number_contributing_weights;
      if (parts[1] == "T") number_contributing_weights = 4;
      else number_contributing_weights = 2;
      // add current_string of the particles at lower levels
      for (size_t k(0); k < number_contributing_weights; ++k) {
        expanded_labels.push_back(prefix);
      }
      for (size_t j(0); j < level; ++j) {
        for (size_t k(0); k < number_contributing_weights; ++k) {
          expanded_labels[k] += label_parts[j] + "_";
        }
      }
      // expand transverse weight (resulting in two to four new spin labels to add per transversely polarized particle)
      expanded_labels[0] += parts[0] + ".++";
      expanded_labels[1] += parts[0] + ".--";
      // interference terms to add for coherent transverse signal definition
      if (parts[1] == "T") {
        expanded_labels[2] += parts[0] + ".+-";
        expanded_labels[3] += parts[0] + ".-+";
      }
      // if not the last particle is currently considered the current_string of the particles with higher levels
      // will be added to the expanded transverse_labels
      if (level + 1 < label_parts.size()) {
        for (size_t k(0); k < number_contributing_weights; ++k) {
          expanded_labels[k] += "_";
        }
      }
      for (size_t j(level + 1); j < label_parts.size(); ++j) {
        for (size_t k(0); k < number_contributing_weights; ++k) {
          expanded_labels[k] += label_parts[j];
        }
      }
      for (const auto & expanded_string : expanded_labels){
        tmp_keys.push_back(expanded_string);
      }
    }
    else{
      // at a new level, all current_strings look the same since they are all based on the same weight, which is passed
      // to the method during its initial call
      // if no vector boson and T or t is found this holds for all spin labels in transverse_labels
      tmp_keys=transverse_labels;
      break;
    }
  }
  return level+1 != num_particles ? ExpandLabels(tmp_keys, level + 1, num_particles) : tmp_keys;
}

void PolWeights_Map::Transverse(bool coherent) {
  // To calculate spin labels describing all possible polarization combinations that includes at least on transverse
  // polarized particle, the following steps are done:

  // 1. Determine all spin labels with transverse VBs
  // runs on p_all_weights since the spin labels in *this are already reduced: generated spin labels will be the basis
  // for the determination of all polarization fractions which need to be added to get the fraction which corresponds
  // to the generated (transverse) spin label; this could (in the case of the coherent transverse polarization
  // definition) also include interference weights, such that spin labels generated from the output of TransverseKeys()
  // need to be contained in p_all_weights (with two instead of one polarization index per particle)
  std::set<std::string> transverse_keys = TransverseKeys(p_all_weights->ListofKeys(), 0, coherent);
  std::set<std::string> interference_weights(m_interference_weights);

  for (const auto & current_string : transverse_keys){
    // 2. Determine contributing polarization weights for each transverse weight and add them together
    std::vector<std::string> weights_to_add = ExpandLabels(std::vector<std::string>(1,
                                                                                    current_string), 0);
    Complex new_weight(0);
    for (const auto & current_weight : weights_to_add){
      new_weight += p_all_weights->find(current_weight)->second;
      // delete all interference contributions which are already included into the transverse definition
      // this is done to know which interference contributions do not contribute to any transverse weight such that they
      // can be totaled to a new (reduced) interference contribution if all transverse weights are determined
      interference_weights.erase(current_weight);
    }
    // 3. Shorten (transverse) spin label and add it to *this
    std::string short_name = ShortName(current_string);
    emplace(short_name, new_weight);
  }
  // 4. calculate remaining interference; separate interference weight for each prefix
  // calculating interference by totaling remaining interference terms instead of subtraction to enable consistency
  // checks
  if (coherent){
    while (!interference_weights.empty()){
      std::string prefix;
      std::string prefix_label = *interference_weights.begin();
      std::string unpol_label;
      Complex new_int_weight(0);
      std::replace(prefix_label.begin(), prefix_label.end(), '_', ' ');
      auto label_parts = ATOOLS::ToVector<std::string>(prefix_label);
      if (label_parts[0] == "dc" || label_parts[0].substr(0, 6) == "Weight") {
        // determine possible prefix
        prefix = label_parts[0];
        label_parts.erase(label_parts.begin());
        // determine label of particles which are considered as unpolarized (if existent)
        // is necessary to name the corresponding coint weight consistent with the other polarization fractions
        for (auto & label_part : label_parts){
          std::replace(label_part.begin(), label_part.end(), '.', ' ');
          std::vector<std::string> parts = ATOOLS::ToVector<std::string>(label_part);
          if (parts[1]=="U") unpol_label += "_" + parts[0] + "." + parts[1];
        }
      }
      // iterate through all interference spin labels and add all ones which have the same prefix; erase all
      // interference weights which were already totaled
      auto it = interference_weights.begin();
      while (it != interference_weights.end()){
        std::string temp_label = *it;
        std::replace(temp_label.begin(), temp_label.end(), '_', ' ');
        std::string potential_prefix = ATOOLS::ToVector<std::string>(temp_label)[0];
        // potential_prefix is not empty for base weights
        if (potential_prefix == prefix || (prefix.empty() && !(potential_prefix == "dc" ||
        potential_prefix.substr(0, 6) == "Weight"))) {
          new_int_weight += p_all_weights->find(*it)->second;
          it = interference_weights.erase(it);
        }
        else it++;
      }
      if (!unpol_label.empty()) unpol_label += "_";
      if (!prefix.empty() && unpol_label.empty()) prefix += "_";
      emplace(prefix + unpol_label + "coint", new_int_weight);
    }
  }
}

// valid input: comma separated of particles numbers according to Sherpa's particle ordering describing which particles
// should be considered as unpolarized in the custom weights
std::vector<std::string> PolWeights_Map::Unpol(const Amplitude2_Tensor *amps,
                                               const std::vector<int>& unpol_particle_numbers, bool non_zero_weight){
  std::vector<std::string> finished_weights;
  if (unpol_particle_numbers.empty()){
    for (auto  &w: m_custom_weights) {
      // add custom spin label as specified in the SHERPA run card since produced spin labels without that would not
      // be unambiguous in all cases
      std::string prefix = w.first + "_";
      std::string current_weight(w.second);
      std::replace(current_weight.begin(),current_weight.end(),',',' ');
      std::vector<int> particle_numbers;
      // test whether custom_weights contains particle numbers
      try{
        particle_numbers = ATOOLS::ToVector<int>(current_weight);
      }
      catch(const ATOOLS::Exception& error){
        continue;
      }
      Amplitude2_Tensor* tmp_amps = new Amplitude2_Tensor(*amps);
      int number_particles = amps->NumberParticles();
      int particle_counter(0);
      // if particle numbers are given, the unpolarized weights for the specified particles are calculated
      if (!particle_numbers.empty()) {
        for (int particle_number : particle_numbers){
          // Sherpa's numbering within the code starts with one while the numbering in the run card starts with zero
          std::pair<int, const ATOOLS::Particle*> particle = amps->Search(particle_number + 1);
          // Contract polarization indices of particles which should be considered as unpolarized with a decay matrix
          // filled with ones
          if (particle.second!=NULL) {
            Amplitude2_Matrix* D = new Amplitude2_Matrix(particle.second, 1);
            tmp_amps->Contract(D);
            delete D;
            // avoid underscore at the end of the spin label for the weight where all particles are unpolarized
            ++particle_counter;
            if (particle_counter==number_particles) prefix += particle.second->RefFlav().IDName() + ".U";
            else prefix += particle.second->RefFlav().IDName() + ".U" + "_";
          }
          else THROW(invalid_input, "Particle with given particle number not found!")
        }
        // Label remaining Amplitude2_Tensor entries and calculate new interference contribution
        LabelAndSeparate(tmp_amps, "start", prefix);
        // Consistency checks with polsum and int
        if (m_pol_checks){
          Tests(prefix);
          // delete polsum from *this since it can be calculated from the other polarized cross sections output
          erase(prefix + "polsum");
        }
        finished_weights.push_back(w.first);
      }
      delete tmp_amps;
    }
  }
  else{
    // partially unpolarized weights for identical particles where particles decaying via a certain decay channel are
    // considered as polarized / unpolarized (private feature, not documented in Sherpa manual)
    Amplitude2_Tensor* tmp_amps = new Amplitude2_Tensor(*amps);
    std::string prefix="dc_";
    for (int unpol_particle_number : unpol_particle_numbers){
      std::pair<int, const ATOOLS::Particle*> particle = amps->Search(unpol_particle_number + 1);
      if (particle.second!=NULL) {
        Amplitude2_Matrix* D = new Amplitude2_Matrix(particle.second, 1);
        tmp_amps->Contract(D);
        delete D;
        // underscore: one intermediate particle is always considered as polarized for this special case (= the one
        // decaying into the user specified decay channel)
        prefix += particle.second->RefFlav().IDName() + ".U" + "_";
      }
      else THROW(invalid_input, "Particle with given particle number not found!")
    }
    LabelAndSeparate(tmp_amps, "start", prefix, non_zero_weight);
    // Consistency checks with polsum and int
    // Tests only make sense if not all fractions are set to zero due to unwanted decay channel combinations appearing
    if (m_pol_checks){
      if (non_zero_weight) Tests(prefix);
      // delete polsum from *this since it can be calculated from the other polarized cross sections output
      erase(prefix + "polsum");
    }
    delete tmp_amps;
  }
  return finished_weights;
}

// valid input: comma separated spin labels which should be added
// custom weights specified by spin labels are named after the corresponding setting in YAML-File (Weight, Weight1, ...
// Weightn)
void PolWeights_Map::AddCustomWeights(const METOOLS::Amplitude2_Tensor* amps,
                                      const std::vector<std::string>& finished_custom_weights){
  for (auto  &w: m_custom_weights) {
    // ignore all custom weights which were already calculated (= partially unpolarized polarization fractions)
    bool next_weight(false);
    for (const auto & finished_custom_weight : finished_custom_weights){
      if (finished_custom_weight==w.first){
        next_weight=true;
        break;
      }
    }
    if (next_weight) continue;

    std::string current_weight(w.second);
    std::replace(current_weight.begin(),current_weight.end(),',',' ');
    auto weights_to_add = ATOOLS::ToVector<std::string>(current_weight);
    Complex new_weight(0.0);
    // Searching user specified weights in PolWeights_Map containing "basic" polarization weights directly from
    // matrix element and in case of VB also transverse weights
    // If the user specified weights are found they are added together (all one which are specified comma separated
    // under one weight setting in YAML-file)
    for (size_t j(0); j<weights_to_add.size(); ++j){
      auto it = find(weights_to_add[j]);
      if (it != end()){
        if (weights_to_add.size()>1){
          new_weight += it->second;
        }
      }
      else{
        // Allows to also add and print out interference weights
        auto it2 = p_all_weights->find(weights_to_add[j]);
        if (it2 != end()) new_weight += it2->second;
          // if one spin label can not be found in PolWeightMap and no particle numbers are used instead
          // the whole custom weight is ignored
        else{
          std::cout << weights_to_add[j] << ", which should be added to a new custom weight, does not exist in "
                                            "PolWeightsMap, ignore the whole custom weight"
                    << std::endl;
          break;
        }
      }
    }
    emplace(w.first, new_weight);
  }
}

// Currently private feature
// TODO: currently only for processes where all intermediate particle are of the same type
// TODO: currently only single polarized cross sections (takes one decay channel which characterizes the particle
//       which should be considered as the only polarized
void PolWeights_Map::AddSinglePolWeights(const METOOLS::Amplitude2_Tensor* amps) {
  METOOLS::Amplitude2_Tensor* tmp_amps = (METOOLS::Amplitude2_Tensor*) amps;
  std::vector<int> unpol_particle_numbers;
  int found_channel(0);
  std::string tmp_string(m_singlepol_channel);
  std::replace(tmp_string.begin(),tmp_string.end(),',',' ');
  int pdg_code = ATOOLS::ToVector<int>(tmp_string)[0];
  do {
    // current particle matches the desired decay channel
    if ((*(tmp_amps->CurrentParticle().OriginalPart()->DecayBlob()))["dc"]->Get<PHASIC::Decay_Channel*>()
          ->IDCode()==m_singlepol_channel){
      if (found_channel==1) unpol_particle_numbers.push_back(tmp_amps->CurrentParticle().Number()-1);
      ++found_channel;
    }
    // current particle is not the searched one, should be unpolarized in the final polarization weights
    else unpol_particle_numbers.push_back(tmp_amps->CurrentParticle().Number()-1);
    tmp_amps = tmp_amps->Next()[0];
  } while (tmp_amps->IsP_Next());
  // Calculate single polarized weights by the help of Unpol method
  if (found_channel==1) Unpol(amps, unpol_particle_numbers, true);
  else if (found_channel==0) {
    unpol_particle_numbers.pop_back();
    Unpol(amps, unpol_particle_numbers, false);
  }
  else if (found_channel==2) Unpol(amps, unpol_particle_numbers, false);
}

void PolWeights_Map::Tests(std::string prefix) {
  // Consistency checks, when PolWeights_Map is finished
  Complex interference = find(prefix + "int") -> second;
  Complex polsum = find(prefix + "polsum") -> second;
  // interference, polsum real?
  if (interference.imag()>1e-8){
    std::cout<<"Polarization_Warning in "<< METHOD <<
             ": Imaginary parts of amplitude2_tensor does not sum up to zero" << std::endl;
    msg_Out() << "imaginary part of interference term: " << std::endl;
    msg_Out() << interference.imag() << std::endl;
  }
  if (polsum.imag()>1e-8){
    std::cout<<"Polarization_Warning in "<< METHOD <<
             ": Sum of polarizations is not real!" << std::endl;
    msg_Out() << "imaginary part of polarization sum: " << std::endl;
    msg_Out() << polsum.imag() << std::endl;
  }
  // Consistency check polsum + int = unpol, unpol real
  if ((m_unpolcrosssec * (interference + polsum)  - m_unpolcrosssec).real() >
      fabs(m_unpolcrosssec.real())*1e-8 ||
      (m_unpolcrosssec * (interference + polsum)).imag() > 1e-8 || m_unpolcrosssec.imag() > 1e-8) {
    std::cout << "Polarization_Warning in " << METHOD <<
              ": Testing consistency between polarization sum + interference and unpolarized result failed"
              << std::endl;
    msg_Out() << "Polarization sum plus interference:" << m_unpolcrosssec * (interference + polsum)
              << std::endl;
    msg_Out() << "Unpolarized result" << m_unpolcrosssec << std::endl;
  }
}
