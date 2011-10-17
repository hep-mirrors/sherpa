#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Amplitude2_Tensor::Amplitude2_Tensor(std::deque<ATOOLS::Particle*> parts) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (parts.size()>0) {
    p_part=parts.front();
    parts.pop_front();

    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      (*p_next)[i]=new Amplitude2_Tensor(parts);
    }
  }
  else {
    m_value=Complex(1.0,0.0);
  }
}

Amplitude2_Tensor::Amplitude2_Tensor(std::deque<Particle*> parts,
                                     const PHASIC::DiagColVec& diagrams,
                                     const vector<vector<Complex> >& cmatrix,
                                     vector<int> spin_i,
                                     vector<int> spin_j) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (parts.size()>0) {
    p_part=parts.front();
    parts.pop_front();
  
    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      vector<int> spin_i_branch(spin_i);
      spin_i_branch.push_back(i%m_nhel);
      vector<int> spin_j_branch(spin_j);
      spin_j_branch.push_back(i/m_nhel);
      (*p_next)[i]=new Amplitude2_Tensor(parts,diagrams,cmatrix,
                                      spin_i_branch, spin_j_branch);
    }
  }
  else {
    m_value=Complex(0.0,0.0);
    for (size_t i(0); i<diagrams.size(); ++i) {
      Spin_Amplitudes* Ai=diagrams[i].first;
      for (size_t j(0); j<diagrams.size(); ++j) {
        Spin_Amplitudes* Aj=diagrams[j].first;
        m_value+=cmatrix[i][j]*Ai->Get(spin_i)*conj(Aj->Get(spin_j));
      }
    }
  }
}

Amplitude2_Tensor::Amplitude2_Tensor(std::deque<ATOOLS::Particle*> parts,
                                     const vector<Spin_Amplitudes>& diagrams,
                                     const std::vector<std::vector<Complex> >& cmatrix,
                                     vector<int> spin_i,
                                     vector<int> spin_j) :
  p_next(NULL), m_value(-1.0,0.0), p_part(NULL), m_nhel(0)
{
  if (parts.size()>0) {
    p_part=parts.front();
    parts.pop_front();
  
    m_nhel=p_part->RefFlav().IntSpin()+1;
    if (m_nhel==3 && IsZero(p_part->RefFlav().Mass())) m_nhel=2;

    p_next=new vector<Amplitude2_Tensor*>(m_nhel*m_nhel);
    for (size_t i=0; i<p_next->size(); ++i) {
      vector<int> spin_i_branch(spin_i);
      spin_i_branch.push_back(i%m_nhel);
      vector<int> spin_j_branch(spin_j);
      spin_j_branch.push_back(i/m_nhel);
      (*p_next)[i]=new Amplitude2_Tensor(parts,diagrams,cmatrix,
                                      spin_i_branch, spin_j_branch);
    }
  }
  else {
    m_value=Complex(0.0,0.0);
    for (size_t i(0); i<diagrams.size(); ++i) {
      Spin_Amplitudes Ai=diagrams[i];
      for (size_t j(0); j<diagrams.size(); ++j) {
        Spin_Amplitudes Aj=diagrams[j];
        m_value+=cmatrix[i][j]*Ai.Get(spin_i)*conj(Aj.Get(spin_j));
      }
    }
  }
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

namespace ATOOLS {
  template <> Blob_Data<METOOLS::Amplitude2_Tensor*>::~Blob_Data() {}
  template class Blob_Data<METOOLS::Amplitude2_Tensor*>;
  template METOOLS::Amplitude2_Tensor*&
  Blob_Data_Base::Get<METOOLS::Amplitude2_Tensor*>();
}
