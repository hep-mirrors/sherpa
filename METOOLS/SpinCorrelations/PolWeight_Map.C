#include "PolWeight_Map.H"
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

using namespace METOOLS;

PolWeight_Map::PolWeight_Map():
  m_unpolcrosssec(Complex(0., 0.)), m_massive_vb(true), m_custom_weights(std::map<std::string, std::string>()),
  m_singlepol_channel("no channel") {
}

PolWeight_Map::PolWeight_Map(const METOOLS::Amplitude2_Tensor* amps, std::map<std::string, std::string> custom_weights,
                             std::string singlepol_channel)
  :m_custom_weights(custom_weights),
   m_singlepol_channel(singlepol_channel){
  // set attribute values
  // unpolarized cross section
  METOOLS::Amplitude2_Tensor* tmpTensor = new METOOLS::Amplitude2_Tensor(*amps);
  m_unpolcrosssec = tmpTensor->Sum();
  // Test, whether unpolarized result is real
  if (m_unpolcrosssec.imag() > 1e-8){
    std::cout<<"Polarization_Warning in"<< METHOD <<
               ": unpolarized result is not real" << std::endl;
    msg_Out() << "imaginary part of the unpolarized result: " << std::endl;
    msg_Out() << m_unpolcrosssec.imag() << std::endl;
  }
  delete tmpTensor;

  m_massive_vb = false;

  // generate map
  Init(amps);
}

PolWeight_Map::~PolWeight_Map() {
}

void PolWeight_Map::Init(const METOOLS::Amplitude2_Tensor* amps, std::string mode, std::string spin_label) {
  if (mode != "start" && mode != "pol" && mode != "int"){
    THROW(fatal_error, "Invalid mode for PolWeight_Map::Init")
  }
  std::vector<Amplitude2_Tensor*> next_amps;
  // searching for and labeling polarized fractions and interference fractions
  if (amps->IsP_Next()){
    next_amps = amps->Next();
    int m_nhel = std::sqrt(next_amps.size());
    if (m_nhel>3){
      THROW(not_implemented,
            "Particles with spin bigger than 1 are not implemented for polarized cross sections yet")
    }
    if (m_nhel==3){
      m_massive_vb = true;
    }

    // lists of strings for the different possible polarization
    std::vector<std::string> spin_strings;
    // TODO: LABELING OF TRANSVERSE POLARIZED MATRIXELEMENTS IS SWITCHED HERE TO GET THE RIGHT LABELING IN THE EVENT
    //       OUTPUT UNTIL SWITCHED ORDERING ISSUE IN MATRIXELEMENT GENERATORS IS FIXED,
    //       POSSIBLE TESTS: DECAY ANGLE OF + -  AND - - POLARIZED VECTOR BOSONS IN VECTOR BOSON PRODUCTION PROCESSES
    spin_strings.push_back("-");
    spin_strings.push_back("+");
    if(m_nhel == 3) {
      spin_strings.push_back("0");
    }
    for (size_t i(0); i<m_nhel*m_nhel; ++i) {
      // polarizations parts
      if(i % (m_nhel + 1) == 0 && mode=="start"){
        Init(next_amps[i], "pol", spin_label + amps->CurrentParticle().RefFlav().IDName() + "."
                                  + spin_strings[i / (m_nhel + 1)]);
      }
      else if (i % (m_nhel + 1) == 0 && mode!="int" && mode!="start"){
        Init(next_amps[i], "pol", spin_label + "_" + amps->CurrentParticle().RefFlav().IDName() + "."
                                  + spin_strings[i / (m_nhel + 1)]);
      }
        // interference parts
      else{
        Init(next_amps[i], "int");
      }
    }
  }
  else{
    // add founded polarisation fraction to the map, key = spinlabel
    if (mode=="pol"){
      emplace(spin_label, amps->Value() / m_unpolcrosssec);
      auto it = find("polsum");
      // calculate sum of polarizations for later tests: looks whether it already exists in the map, if not,
      // generate a new one with the current polarization fraction, if yes, find the corresponding value, add
      // the current fraction and write it to the map
      if (it==end()){
        emplace("polsum", amps->Value() / m_unpolcrosssec);
      }
      else{
        Complex tmp = it->second;
        tmp += amps->Value() / m_unpolcrosssec;
        (*this)["polsum"] = tmp;
      }
    }
      // add interference term to the map, if it does not exist, otherwise adding the current interference term to
      // one already existing in the map
    else if (mode=="int"){
      auto it1 = find("int");
      if (it1==end()){
        emplace("int", amps->Value() / m_unpolcrosssec);
      }
      else{
        Complex tmp = it1->second;
        tmp += amps->Value() / m_unpolcrosssec;
        (*this)["int"] = tmp;
      }
    }
    else{
      THROW(fatal_error, "No Tensor")
    }
  }

    if (mode=="start") {
        // Tests, when PolWeight_Map is finished
        Complex interference = find("int") -> second;
        Complex polsum = find("polsum") -> second;
        if (interference.imag()>1e-8){
            std::cout<<"Polarization_Warning in"<< METHOD <<
            ": Imaginary parts of amplitude2_tensor does not sum up to zero" << std::endl;
            msg_Out() << "imaginary part of interference term: " << std::endl;
            msg_Out() << interference.imag() << std::endl;
        }
        if (polsum.imag()>1e-8){
          std::cout<<"Polarization_Warning in"<< METHOD <<
                   ": Sum of polarizations is not real!" << std::endl;
          msg_Out() << "imaginary part of polarization sum: " << std::endl;
          msg_Out() << polsum.imag() << std::endl;
        }
        if ((m_unpolcrosssec * (interference + polsum)  - m_unpolcrosssec).real() >
        fabs(m_unpolcrosssec.real())*1e-8 ||
        (m_unpolcrosssec * (interference + polsum)).imag() > 1e-8 || m_unpolcrosssec.imag() > 1e-8){
      std::cout<<"Polarization_Warning in"<< METHOD <<
               ": Testing consistency between polarisation sum + interference and unpolarized result failed" << std::endl;
      msg_Out() << "Polarisation sum plus interference:" << m_unpolcrosssec * (interference + polsum)
                << std::endl;
      msg_Out() << "Unpolarized result" << m_unpolcrosssec << std::endl;
    }

    // Calculation of transverse weights
    // VB as well as custom weights first generate a separate PolWeightsMap since it always needs the keys of the
    // original weights
    if (m_massive_vb){
      PolWeight_Map TransverseWeights = Transverse(0);
      insert(TransverseWeights.begin(), TransverseWeights.end());
    }

    std::vector<PolWeight_Map> temp_maps;
    // Add user specified custom weights
    if (!m_custom_weights.empty()){
      temp_maps.push_back(AddCustomWeights(amps));
    }
    // Add single channel weights
    if (m_singlepol_channel!="no channel"){
      temp_maps.push_back(AddSinglePolWeights(amps));
    }
    // Add additional weights to overall PolWeightsMap only if all additional weights are calculated since
    // their calculation may expect an overall PolWeightsMap containing only basis polarization weights
    for (size_t i(0); i<temp_maps.size(); ++i){
      insert(temp_maps[i].begin(), temp_maps[i].end());
    }
  }
}

std::set<std::string> PolWeight_Map::ListofKeys() const{
  std::set<std::string> keys;
  for (auto const& element : *this) {
    keys.emplace(element.first);
  }
  return keys;
}

PolWeight_Map PolWeight_Map::Transverse(int level) const{
  // Generating of transverse weights is done by going through all the existing polarization weights which come
  // directly out of the matrix element and for each particle adding + - and - helicity contribution for each
  // polarization weight which are the same beside the helicity of this particle
  // The first particle in the key names runs on the key list of the original weight map; all other particles
  // only run on the extra polarization map which than containing only partially transverse weights and the weights
  // describing longitudinal contributions of particles at former levels
  PolWeight_Map temp_map;
  std::set<std::string> keys = ListofKeys();
  keys.erase("int");
  keys.erase("polsum");
  int num_particles(level);
  // look at each key in the key list and check whether the level particle is a vector boson with + or - helicity
  // if yes, generate the corresponding "partner" weight name were + is replaced by - helicity or vice versa
  while (!keys.empty()){
    auto it = keys.begin();
    std::string temp_label = (*it);
    std::replace(temp_label.begin(),temp_label.end(),'_',' ');
    auto label_parts = ATOOLS::ToVector<std::string>(temp_label);
    num_particles = label_parts.size();
    std::replace(label_parts[level].begin(), label_parts[level].end(),'.',' ');
    auto parts = ATOOLS::ToVector<std::string>(label_parts[level]);
    if (parts[0]=="W+" || parts[0]=="W-" || parts[0]=="Z"){
      if(parts[1]=="+" || parts[1]=="-"){
        // generate partner string to current key and the weight name for the new transverse weight
        // string of particles at smaller levels
        std::string other_string = "";
        std::string new_string = "";
        for (size_t j(0); j<level; ++j) {
          other_string += label_parts[j];
          other_string += "_";
          new_string += label_parts[j];
          new_string += "_";
        }
        // string of particle at the current level
        other_string += parts[0];
        new_string += parts[0];
        new_string += ".T";
        if (parts[1]=="+"){
          other_string += ".-";
        }
        else{
          other_string += ".+";
        }
        // if not the last particle is currently considered the string of the particles with higher level
        // will be added to searched partner string and new weightname
        if (level+1 < label_parts.size()){
          new_string += "_";
          other_string += "_";
        }
        for (size_t j(level+1); j<label_parts.size(); ++j){
          other_string += label_parts[j];
          new_string += label_parts[j];
        }
        // Searching the two weight names of + - and - helicity of current VB in the weights map, adding them
        // together and store the new weight in a temporal map
        Complex new_weight = find(other_string)->second + find((*it))->second;
        temp_map.emplace(new_string, new_weight);
        keys.erase(keys.find(other_string));
      }
      else{
        // Adding the weights where the current particle is longitudinal since VB at higher levels can have + or
        // - contribution
        temp_map.emplace((*it), find((*it))->second);
      }
    }
      // current particle is not a vector boson; since all particles at one level have the same flavour, one can
      // directly to the next level in that case
    else{
      temp_map = (*this);
      break;
    }
    keys.erase(it);
  }
  return level+1 != num_particles ? temp_map.Transverse(level + 1) : temp_map;
}

// valid input: comma separated weight names which should be added or particles numbers according to Sherpa ordering
// describing which particles should be considered as unpolarized in the custom weights
// custom weights specified by weight names are named after the corresponding setting in YAML-File (Weight, Weight1, ...
// Weightn)
PolWeight_Map PolWeight_Map::AddCustomWeights(const METOOLS::Amplitude2_Tensor* amps) const{
  PolWeight_Map tmp_map;
  for (auto  &w: m_custom_weights) {
    std::string current_weight(w.second);
    std::replace(current_weight.begin(),current_weight.end(),',',' ');
    auto weights_to_add = ATOOLS::ToVector<std::string>(current_weight);
    Complex new_weight(0.0);
    bool undefined_weight(false), test_particle_numbers(false);
    // Searching user specified weights in PolWeight_Map containing "basic" polarization weights directly from
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
        // if already the first weight could not be found it could be possible that particle numbers are given
        // instead of weight names
      else if (j==0){
        test_particle_numbers = true;
        break;
      }
        // if one weight name can not be found in PolWeightMap and no particle numbers are used instead
        // the whole custom weight is ignored
      else{
        std::cout << "Weight " << weights_to_add[j] << ", which should be added to a new custom weight, does not"
                                                       " exist in PolWeightsMap; ignore the whole custom weight"
                  << std::endl;
        undefined_weight = true;
        break;
      }
    }
    // if the given custom weight contains only valid weight names, calculated custom weight is write to the
    // temporal map which storing the custom weights
    if (!undefined_weight && !test_particle_numbers){
      tmp_map.emplace(w.first, new_weight);
    }
    // if it seems that no weight names are used test whether particle numbers are used for specifing partially
    // unpolarized weights
    //TODO: More user friendly error for wrong weight names then failed to parse error in ATOOLS::ToVector
    if (test_particle_numbers){
      auto particle_numbers = ATOOLS::ToVector<int>(current_weight);
      // if particle numbers are given the unpolarized weights for the specified particles are calculated
      if (!particle_numbers.empty()){
        PolWeight_Map unpol_weights = Unpol(amps, particle_numbers);
        tmp_map.insert(unpol_weights.begin(), unpol_weights.end());
      }
      else{
        std::cout << "Weight " << current_weight << ", which should be added to a new custom weight, does "
                                                    "not exist in PolWeightsMap and also does not describe"
                                                    " hard process particles with possibly given particle "
                                                    "numbers; ignore the whole custom weight"
                  << std::endl;
      }
    }
  }
  return tmp_map;
}

PolWeight_Map PolWeight_Map::Unpol(const METOOLS::Amplitude2_Tensor* amps, std::vector<int> particle_numbers) const {
  // principe of this method is similar to Transverse method with the difference that here only polarizations of
  // the particles are added which were specified as should be considered unpolarized by the user
  PolWeight_Map temp_map;
  std::set<std::string> keys = ListofKeys();
  keys.erase("int");
  keys.erase("polsum");
  // place in weightnames which describes particle with particle number particle_numbers[0]
  std::vector<std::string> spin_strings;
  // searching for particle with the desired particle number
  // +1 due to difference in YAML-particle numbering and Sherpa internal particle numbering
  //TODO: Discrepancy between internal particle numbering in Sherpa and YAML-particle numbering?
  std::pair<int, ATOOLS::Particle> particle = amps->Search(particle_numbers[0]+1);
  int place(particle.first);
  if (particle.second.Number()!=-1){
    int m_nhel=particle.second.RefFlav().IntSpin()+1;
    if (m_nhel==3 && ATOOLS::IsZero(particle.second.RefFlav().Mass())) m_nhel=2;
    // TODO: LABELING OF TRANSVERSE POLARIZED MATRIXELEMENTS IS SWITCHED HERE TO GET THE RIGHT LABELING IN THE EVENT
    //       OUTPUT UNTIL SWITCHED ORDERING ISSUE IN MATRIXELEMENT GENERATORS IS FIXED,
    //       POSSIBLE TESTS: DECAY ANGLE OF + -  AND - - POLARIZED VECTOR BOSONS IN VECTOR BOSON PRODUCTION PROCESSES
    spin_strings.push_back("-");
    spin_strings.push_back("+");
    if(m_nhel == 3) {
      spin_strings.push_back("0");
      spin_strings.push_back("T");
    }
  }
    // ignore particle number if the corresponding particle can not be found in the Amplitude2_Tensor
  else{
    std::cout << "Particle number " << particle_numbers[0] << " does not match any hard decaying particle! "
                                                              "This particle number will be ignored in the "
                                                              "following!" << std::endl;
    particle_numbers.erase(particle_numbers.begin());
    return particle_numbers.empty() ? (*this) : Unpol(amps, particle_numbers);
  }
  do{
    // generating the labels of all weights which should be added (only add weights which are identical beside the
    // helicity of the current particle which should be unpolarized
    // in addition to that, the weight name of the new weight is generated
    auto it = keys.begin();
    std::string temp_label = (*it);
    std::replace(temp_label.begin(),temp_label.end(),'_',' ');
    auto label_parts = ATOOLS::ToVector<std::string>(temp_label);
    std::replace(label_parts[place-1].begin(), label_parts[place-1].end(),'.',' ');
    auto parts = ATOOLS::ToVector<std::string>(label_parts[place-1]);
    std::string other_string = "";
    std::string new_string = "";
    // weightname stays the same before the unpolarized particle
    for (size_t j(0); j<place-1; ++j) {
      other_string += label_parts[j];
      other_string += "_";
      new_string += label_parts[j];
      new_string += "_";
    }
    other_string += parts[0];
    new_string += parts[0];
    new_string += ".u";
    if (place < label_parts.size()) {
      new_string += "_";
    }
    // weightname stays the same after the unpolarized particle
    for (size_t k(place);k < label_parts.size(); ++k) {
      new_string += label_parts[k];
    }
    if (temp_map.find(new_string)!=temp_map.end()){
      keys.erase(it);
      continue;
    }
    Complex new_weight(0.0, 0.0);
    size_t added_trans_weights(0);
    // To avoid double adding transverse weight (once through unpolarizedparticle.T and one through
    // unpolarizedparticle.+ and unpolarizedparticle.-) to unpolarizedparticle.u
    for (size_t j(0); j<spin_strings.size(); ++j) {
      if (spin_strings[j]=="T" && added_trans_weights==2){
        break;
      }
      else if (spin_strings[j]=="T" && added_trans_weights!=0 && added_trans_weights!=2){
        std::cout << "added_trans_weights" << added_trans_weights << std::endl;
        THROW(fatal_error, "Internal error while calculating unpolarized polarization weights")
      }
      // finishing generation of the weight names to add
      std::string tmp_string(other_string);
      tmp_string += "." + spin_strings[j];
      if (place < label_parts.size()) {
        tmp_string += "_";
      }
      for (size_t k(place);k < label_parts.size(); ++k) {
        tmp_string += label_parts[k];
      }
      // calculating new weight and adding it to temporal map
      auto it2 = find(tmp_string);
      // in connection with transverse weights, weight names like W+.T_Z.+ were generated above which can not be
      // found in the overall PolWeightMap, but this is not a problem, since for e.g. an unpolarized Z only
      // W+.T_Z.T and W+.T_Z.0 needs to be added to get W+_T.Z.u
      if (it2!=end()){
        if (spin_strings[j]=="+" || spin_strings[j]=="-"){
          added_trans_weights += 1;
        }
        new_weight += it2->second;
      }
      if (keys.find(tmp_string)!=keys.end()){
        keys.erase(keys.find(tmp_string));
      }
    }
    temp_map.emplace(new_string, new_weight);
  } while (!keys.empty());
  particle_numbers.erase(particle_numbers.begin());
  return particle_numbers.empty() ? temp_map : temp_map.Unpol(amps, particle_numbers);
}

//TODO: Shall single polarized weights be part of the official implementation? If yes, add manual entry for that.
PolWeight_Map PolWeight_Map::AddSinglePolWeights(const METOOLS::Amplitude2_Tensor* amps) const{
  METOOLS::Amplitude2_Tensor* tmp_amps = (METOOLS::Amplitude2_Tensor*) amps;
  PolWeight_Map tmpMap;
  std::vector<int> unpol_particle_numbers;
  bool found_channel(false);
  std::string tmp_string(m_singlepol_channel);
  std::replace(tmp_string.begin(),tmp_string.end(),',',' ');
  int pdg_code = ATOOLS::ToVector<int>(tmp_string)[0];
  std::vector<std::string> particle_strings;
  do {
    // Generating keys for the new polarization weights
    if (tmp_amps->CurrentParticle().RefFlav().operator long()==pdg_code && particle_strings.empty()){
      std::string particle_string = tmp_amps->CurrentParticle().RefFlav().IDName();
      int m_nhel = std::sqrt(tmp_amps->Next().size());
      particle_strings.push_back(particle_string + ".+");
      particle_strings.push_back(particle_string +".-");
      if(m_nhel == 3) {
        particle_strings.push_back(particle_string +".0");
        particle_strings.push_back(particle_string +".T");
      }
    }
    // current particle matches the desired decay channel
    if ((*(tmp_amps->CurrentParticle().OriginalPart()->DecayBlob()))["dc"]->Get<PHASIC::Decay_Channel*>()
            ->IDCode()==m_singlepol_channel){
      // two hard decaying particles have the same decay channel, polarized particle is not unique, set single
      // polarized weights to zero for this event
      if(found_channel){
        for (size_t i(0); i<particle_strings.size(); ++i){
          tmpMap.emplace(particle_strings[i], Complex(0.0, 0.0));
        }
        return tmpMap;
      }
      else{
        found_channel = true;
      }
    }
    else{
      // current particle is not the searched one, should be unpolarized in the final polarization weights
      unpol_particle_numbers.push_back(tmp_amps->CurrentParticle().Number()-1);
    }
    tmp_amps = tmp_amps->Next()[0];
  } while (tmp_amps->IsP_Next());

  if (particle_strings.empty()){
    std::cout << "The decaying particle specified in the channel for single polarized cross sections is not a "
                 "hard decaying particle! Process further without calculating single polarized cross sections!" <<
              std::endl;
  }
  // Calculate single polarized weights by the help of Unpol method
  if (found_channel){
    PolWeight_Map tmptmpMap = Unpol(amps, unpol_particle_numbers);
    if (tmptmpMap.size() != particle_strings.size()){
      THROW(fatal_error, "Internal error while calculating single polarized cross sections!")
    }
    // Adjust weight names
    size_t index(0);
    for (auto  &e: tmptmpMap) {
      tmpMap.emplace(particle_strings[index], e.second);
      ++index;
    }
  }
    // all hard decaying particles decays into different channels than the desired one, set the single polarized weights
    // to zero for this event
  else{
    for (size_t i(0); i<particle_strings.size(); ++i){
      tmpMap.emplace(particle_strings[i], Complex(0.0, 0.0));
    }
  }
  return tmpMap;
}
