#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Poincare.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Amplitude2_Tensor::Amplitude2_Tensor(const std::vector<ATOOLS::Particle*>& parts,
                                     size_t level) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (level>parts.size()) THROW(fatal_error, "Internal error 1");

  if (level==parts.size()) {
    m_value=Complex(1.0,0.0);
  }
  else {
    p_part=parts[level];

    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      (*p_next)[i]=new Amplitude2_Tensor(parts, level+1);
    }
  }
}


Amplitude2_Tensor::Amplitude2_Tensor(const vector<ATOOLS::Particle*>& parts,
                                     size_t level,
                                     const std::vector<Spin_Amplitudes*>& diagrams,
                                     vector<int>& spin_i,
                                     vector<int>& spin_j) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (level>parts.size()) THROW(fatal_error, "Internal error 1");

  if (level==parts.size()) {
    m_value=Complex(0.0, 0.0);
    for (size_t i(0); i<diagrams.size(); ++i) {
      for (size_t j(0); j<diagrams.size(); ++j) {
        m_value+=diagrams[i]->Get(spin_i)*
          conj(diagrams[j]->Get(spin_j));
      }
    }
  }
  else {
    p_part=parts[level];

    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      spin_i[level]=(i%m_nhel);
      spin_j[level]=(i/m_nhel);
      (*p_next)[i]=new Amplitude2_Tensor(parts, level+1,
                                         diagrams, spin_i, spin_j);
    }
  }
}


Amplitude2_Tensor::Amplitude2_Tensor(const vector<ATOOLS::Particle*>& parts,
                                     const vector<int>& permutation,
                                     size_t level,
                                     const vector<Spin_Amplitudes>& diagrams,
                                     vector<int>& spin_i, vector<int>& spin_j) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (level>parts.size()) THROW(fatal_error, "Internal error 1");

  if (level==parts.size() || parts[level]->RefFlav().IsStable()) {
    m_value=ContractRemaining(parts,permutation,level,diagrams,
                              spin_i,spin_j, 1.0);
  }
  else {
    p_part=parts[level];

    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      spin_i[level]=(i%m_nhel);
      spin_j[level]=(i/m_nhel);
      (*p_next)[i]=new Amplitude2_Tensor(parts, permutation, level+1,
                                         diagrams, spin_i, spin_j);
    }
  }
}


Amplitude2_Tensor::Amplitude2_Tensor(const Amplitude2_Tensor& other)
{
  m_value=other.m_value;
  m_nhel=other.m_nhel;
  p_part=other.p_part;

  if (other.p_next) {
    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      (*p_next)[i]=new Amplitude2_Tensor(*(other.p_next->at(i)));
    }
  }
  else p_next=NULL;
}

Complex Amplitude2_Tensor::ContractRemaining
(const std::vector<ATOOLS::Particle*>& parts,
 const vector<int>& permutation,
 size_t level,
 const vector<Spin_Amplitudes>& diagrams,
 vector<int>& spin_i, vector<int>& spin_j, double factor) const
{
  if (level>parts.size()) THROW(fatal_error, "Internal error 1");

  Complex ret(0.0, 0.0);

  if (level==parts.size()) {
    vector<int> spin_i_perm(spin_i.size()), spin_j_perm(spin_j.size());
    for (size_t p=0; p<spin_i.size(); ++p) {
      spin_i_perm[p]=spin_i[permutation[p]];
      spin_j_perm[p]=spin_j[permutation[p]];
    }
    for (size_t i(0); i<diagrams.size(); ++i) {
      for (size_t j(0); j<diagrams.size(); ++j) {
        ret+=diagrams[i].Get(spin_i_perm)*
          conj(diagrams[j].Get(spin_j_perm))*factor;
      }
    }
  }
  else {
    int nlambda=parts[level]->RefFlav().IntSpin()+1;
    if (nlambda==3 && IsZero(parts[level]->RefFlav().Mass())) nlambda=2;
    double newfactor=factor/double(nlambda);
    for (size_t i=0; i<nlambda; ++i) {
      spin_i[level]=i;
      spin_j[level]=i;
      ret+=ContractRemaining(parts, permutation, level+1,
                             diagrams, spin_i, spin_j, newfactor);
    }
  }
  return ret;
}

Amplitude2_Tensor::~Amplitude2_Tensor()
{
  if (p_next) {
    for (size_t i=0; i<p_next->size(); ++i) {
      if ((*p_next)[i]) {
        delete (*p_next)[i];
        (*p_next)[i]=NULL;
      }
    }
    delete p_next;
  }
}

void Amplitude2_Tensor::Contract(const Amplitude2_Matrix* D) {
  const Particle* part=D->Particle();
  DEBUG_FUNC(*part);
  DEBUG_VAR(Trace());
  if (part==p_part) {
    if (p_next) {
      DEBUG_INFO("found. summing hels.");
      (*p_next)[0]->Multiply((*D)[0]);
      for (size_t i=1; i<p_next->size(); ++i)
        (*p_next)[0]->Add((*p_next)[i], (*D)[i]);

      DEBUG_INFO("deleting all but remaining.");
      for (size_t i=1; i<p_next->size(); ++i) delete (*p_next)[i];
      Amplitude2_Tensor* tmp=(*p_next)[0];

      DEBUG_INFO("setting the remaining as this.");
      p_part=tmp->p_part;
      tmp->p_part=NULL;
      m_value=tmp->m_value;
      m_nhel=tmp->m_nhel;
      tmp->m_nhel=0;
      if (tmp->p_next) {
        p_next->clear();
        p_next->insert(p_next->end(), tmp->p_next->begin(), tmp->p_next->end());
        tmp->p_next->clear();
      }
      else {
        delete p_next;
        p_next=NULL;
      }
      delete tmp;
    }
    else THROW(fatal_error, "Particle not found");
  }
  else {
    DEBUG_INFO("not here. looking further down the tree.");
    if (p_next) {
      for (size_t i(0);i<p_next->size();++i) {
        (*p_next)[i]->Contract(D);
      }
      DEBUG_INFO("finished");
    }
    else THROW(fatal_error, "Particle not found");
  }
  DEBUG_VAR(Trace());
}      

Amplitude2_Matrix Amplitude2_Tensor::ReduceToMatrix(const Particle* left) const
{
  if (!p_part || !p_next) THROW(fatal_error, "Internal1");

  Amplitude2_Matrix sigma(left);
  if (p_part==left) {
    for (size_t i(0); i<p_next->size(); ++i) {
      sigma[i]=(*p_next)[i]->Trace();
    }
  }
  else {
    sigma.assign(sigma.size(), Complex(0.0,0.0));
    // contract with delta
    // have to normalise delta_ij?
    Complex OneOverN=Complex(1.0/double(m_nhel), 0.0);
    for (size_t i(0); i<m_nhel; ++i) {
      sigma.Add((*p_next)[i*m_nhel+i]->ReduceToMatrix(left),OneOverN);
    }
  }
  return sigma;
}




void Amplitude2_Tensor::Add(const Amplitude2_Tensor* amp, const Complex& factor)
{
  if (p_part!=amp->p_part) THROW(fatal_error,"Particles don't match.");
  if (p_next) {
    if (p_next->size() != amp->p_next->size()) THROW(fatal_error, "Internal1.");
    for (size_t i(0);i<p_next->size();++i) {
      (*p_next)[i]->Add((*amp->p_next)[i], factor);
    }
  }
  else {
    if (m_value==Complex(-1.0,0.0) || amp->m_value==Complex(-1.0,0.0))
      THROW(fatal_error, "Internal2.");
    if (amp->p_next) THROW(fatal_error, "Internal3.");
    m_value+=factor*amp->m_value;
  }
}

void Amplitude2_Tensor::Multiply(const Complex& factor)
{
  if (p_next) {
    for (size_t i(0);i<p_next->size();++i) {
      (*p_next)[i]->Multiply(factor);
    }
  }
  else m_value*=factor;
}

void Amplitude2_Tensor::Multiply(const Amplitude2_Matrix* D)
{
  const Particle* part=D->Particle();
  if (part==p_part) {
    if (p_next) {
      if (p_next->size() != D->size()){
        THROW(fatal_error, "InternalError");
      }
      for (size_t i=0; i<p_next->size(); ++i)
        (*p_next)[i]->Multiply((*D)[i]);
    }
    else THROW(fatal_error, "Particle not found");
  }
  else {
    if (p_next) {
      for (size_t i(0);i<p_next->size();++i) {
        (*p_next)[i]->Multiply(D);
      }
    }
    else THROW(fatal_error, "Particle not found");
  }
}

Complex Amplitude2_Tensor::Sum()
{
  while (p_next) {
    Amplitude2_Matrix* Ones = new Amplitude2_Matrix(p_part);
    fill(Ones->begin(), Ones->end(), Complex(1.0, 0.0));
    Contract(Ones);
    delete Ones;
  }
  return m_value;
}

Complex Amplitude2_Tensor::Trace() const
{
  if (!p_part) {
    //if (m_value<0.0) THROW(fatal_error, "Internal.");
    return m_value;
  }
  else {
    size_t pos(0);
    Complex val(0.,0.);
    for (size_t i=0; i<m_nhel; ++i) {
      val += (*p_next)[pos]->Trace();
      pos += m_nhel+1;
    }
    return val;
  }
}

bool Amplitude2_Tensor::Contains(const ATOOLS::Particle* part) const
{
  if (p_part==part) {
    return true;
  }
  else {
    if (p_next) {
      for (size_t i(0);i<p_next->size();++i) {
        if ((*p_next)[i]->Contains(part)) return true;
      }
    }
  }
  return false;
}

void Amplitude2_Tensor::UpdateParticlePointers(const std::map<Particle*,Particle*>& pmap)
{
  if (p_part) {
    std::map<Particle*,Particle*>::const_iterator pit(pmap.find(p_part));
    if (pit!=pmap.end()) p_part=pit->second;
    else THROW(fatal_error, "Could not update particle pointer.");
  }
  if (p_next) {
    for (size_t i(0);i<p_next->size();++i) {
      (*p_next)[i]->UpdateParticlePointers(pmap);
    }
  }
}

void Amplitude2_Tensor::Print(std::ostream& ostr, string label) const
{
  if (m_value!=Complex(-1.0,0.0)) {
    ostr<<"  "<<label<<": "<<m_value<<endl;
  }
  else if (p_next) {
    for (size_t i=0; i<p_next->size(); ++i) {
      (*p_next)[i]->Print(ostr,
          label+" "+ToString(p_part->Flav())+"["+ToString(i)+"]");
    }
  }
  else {
    ostr<<"  nothing here yet, ";
  }
}

namespace METOOLS {
  std::ostream& operator<<(std::ostream& ostr, const Amplitude2_Tensor& t) {
    t.Print(ostr, "");
    return ostr;
  }
}

bool Amplitude2_Tensor::SortCrit(const pair<Particle*, size_t>& p1,
                                        const pair<Particle*, size_t>& p2)
{
  return p1.first->RefFlav().IsStable()<p2.first->RefFlav().IsStable();
}

void Amplitude2_Tensor::PolBasisTrafo(const std::vector<std::vector<std::vector<Complex> > >& coeff,
                                      const std::vector<std::vector<std::vector<Complex> > >& conj_coeff, int level,
                                      std::vector<std::vector<std::vector<Complex> > > coeff_tmp,
                                      std::vector<std::vector<std::vector<Complex> > > conj_coeff_tmp,
                                      Amplitude2_Tensor* old_amps) {
  // TODO: Tests of properties of coefficient matrices like unitarity
  // TODO: calculate conj_coeff here such that no input necessary?
  // TODO: Tests of some fundamental properties of transformed Amplitude2_Tensor like all polarized entries are real,
  //       unpol identical, real & overall interferences zero ...
  // copy old amplitude2_tensor since each component of the original one will be overridden but depends one all the
  // entries of old tensor
  if (level == 0) {
    old_amps = new METOOLS::Amplitude2_Tensor(*this);
    coeff_tmp = std::vector<std::vector<std::vector<Complex> > >(NumberParticles(),
                                                               std::vector<std::vector<Complex> >(1,
                                                                 std::vector<Complex>(1, Complex())));
    conj_coeff_tmp = std::vector<std::vector<std::vector<Complex> > >(NumberParticles(),
                                                                    std::vector<std::vector<Complex> >(1,
                                                                      std::vector<Complex>(1, Complex())));
  }
  // runs through the original, in this context new amplitude2_tensor and selects for each entry the transformation
  // coefficients according to the helicities of the entry under consideration
  if (level < coeff.size()) {
    if (p_next) {
      for (size_t i(0); i < p_next->size(); ++i) {
        // coeff_tmp and conj_coeff_tmp are necessary since coeff, conj_coeff must be unchanged for the next loop step
        coeff_tmp[level][0]=coeff[level][i - (i / m_nhel) * m_nhel];
        conj_coeff_tmp[level][0]=conj_coeff[level][i / m_nhel];
        (*p_next)[i]->PolBasisTrafo(coeff, conj_coeff, level+1, coeff_tmp, conj_coeff_tmp, old_amps);
      }
    }
    else {
      THROW(fatal_error, "size of coefficient vector does not match number of particles in Amplitude2_Tensor")
    }
  }
    // if entering the end of the orignal amplitude tensor, its new value is calculated by running through the
    // old (copied) amplitude tensor
    // recursive calling of this function by its own now determine the corresponding coefficients for the
    // transformation of the single (old) entries from the coefficient matrices filtered above (see next else if)
    // and multiply it with them
    // if this is done for all entries of the old amplitude2_tensor all entries of the tensor are summed up to
    // build the entry of the new tensor considered above
  else if (level == coeff.size()) {
    if (p_next) {THROW(fatal_error,
                       "size of coefficient vector does not match number of particles in Amplitude2_Tensor")}
    if (old_amps->IsP_Next()) {
      Amplitude2_Tensor* old_amps_tmp = new Amplitude2_Tensor(*old_amps);
      old_amps_tmp->PolBasisTrafo(coeff_tmp, conj_coeff_tmp, level+1,
                                  std::vector<std::vector<std::vector<Complex> > >(old_amps->NumberParticles(),
                                                                                 std::vector<std::vector<Complex> >(1,
                                                                                   std::vector<Complex>(1, Complex()))),
                                  std::vector<std::vector<std::vector<Complex> > >(old_amps->NumberParticles(),
                                                                                 std::vector<std::vector<Complex> >(1,
                                                                                   std::vector<Complex>(1, Complex()))));
      m_value = old_amps_tmp->Sum();
      delete old_amps_tmp;
    }
    else {
      THROW(fatal_error, "Internal error")
    }
  }
    // determine transformation coefficients for the old entries to build up the new one considered above
  else if (level > coeff.size() && level < coeff.size() * 2 + 1) {
    if (p_next) {
      for (size_t j(0); j < p_next->size(); ++j) {
        coeff_tmp[level - coeff.size() - 1][0][0]=coeff[level - coeff.size() - 1][0][j - (j / m_nhel) * m_nhel];
        conj_coeff_tmp[level - coeff.size() - 1][0][0]=conj_coeff[level - coeff.size() - 1][0][j / m_nhel];
        (*p_next)[j]->PolBasisTrafo(coeff, conj_coeff, level + 1, coeff_tmp, conj_coeff_tmp);
      }
    }
    else {
      THROW(fatal_error, "size of coefficient vector does not match number of particles in Amplitude2_Tensor")
    }
  }
    // Multiplication of transformation coefficients with the old tensor entry
  else if (level == coeff.size() * 2 + 1) {
    if (p_next) {THROW(fatal_error,
                       "size of coefficient vector does not match number of particles in Amplitude2_Tensor")}
    for (size_t k(0); k < coeff.size(); ++k) {
      m_value *= coeff_tmp[k][0][0] * conj_coeff_tmp[k][0][0];
    }
  }
  else {
    THROW(fatal_error, "Internal error")
  }
  if (level==0){
    delete old_amps;
  }
}

bool Amplitude2_Tensor::IsP_Next() const {
  if (p_next){
    return true;
  }
  return false;
}

int Amplitude2_Tensor::NumberParticles(int num) const {
  if (p_next) {
    num += 1;
    return (*p_next)[0]->NumberParticles(num);
  }
  else{
    return num;
  }
}

std::pair<const int, const ATOOLS::Particle> Amplitude2_Tensor::Search(const int part_number, int level) const {
  if (p_part->Number()==part_number) {
    return std::make_pair(level, *p_part);
  }
  else {
    if (p_next) {
      return (*p_next)[0]->Search(part_number, level+1);
    }
  }
  return std::make_pair(level, ATOOLS::Particle());
}

namespace ATOOLS {

  template <> Blob_Data<Amplitude2_Tensor_SP>::~Blob_Data() {}
  template class Blob_Data<Amplitude2_Tensor_SP>;
  template Amplitude2_Tensor_SP& Blob_Data_Base::Get<Amplitude2_Tensor_SP>();

}
