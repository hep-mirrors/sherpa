#include "Spin_Correlation_Tensor.H"
#include "Message.H"
#include "Blob.H"
#include "Smart_Pointer.H"
#include "prof.hh"
#include <typeinfo>


  //#define SCT_Debug

namespace ATOOLS{

  // INITIALIZATION OF THE STATIC MEMBERS OF SPIN_CORRELATION_TENSOR

  int Spin_Correlation_Tensor::m_k0_n(-1);
  std::set<int> Spin_Correlation_Tensor::m_possible_particles;
  scmode::code Spin_Correlation_Tensor::m_mode(scmode::None);
  bool Spin_Correlation_Tensor::m_created(false);


  // ACCESS METHODS FOR THE STATIC MEMBERS

  void Spin_Correlation_Tensor::AddPossibleParticles( std::set<kf_code> * flavs )
  {
    for (std::set<kf_code>::iterator iter=flavs->begin(); iter!=flavs->end(); ++iter)
      m_possible_particles.insert(int(*iter));  
  }

  void Spin_Correlation_Tensor::AddPossibleParticle( kf_code flav )
  {
    m_possible_particles.insert(int(flav));
  }
  
  bool Spin_Correlation_Tensor::PossibleParticle( kf_code flav)
  {
    for (std::set<int>::iterator iter=m_possible_particles.begin(); 
	 iter!=m_possible_particles.end(); ++iter)
      if (int(flav) == (*iter)) return true;
    return false;
  }
   
  void Spin_Correlation_Tensor::PrintPossibleParticles()
  {
    for (std::set<int>::iterator iter=m_possible_particles.begin(); 
	 iter!=m_possible_particles.end(); ++iter)
      std::cout<<(*iter)<<" ";  
    std::cout<<std::endl;
  }
  
  void Spin_Correlation_Tensor::SetMode(scmode::code i)
  {
    if ( m_created && (i != m_mode) ) {
      msg_Error()<<"Warning in Spin_Correlation_Tensor::SetMode:"<<std::endl;
      msg_Error()<<" New spin-correlation mode "<<i<<" was selected but"<<std::endl;
      msg_Error()<<" some Spin_Correlation_Tensor objects were already created."<<std::endl;
      msg_Error()<<" Old mode was "<<m_mode<<std::endl;
    }
    m_mode = i;
  }

  scmode::code Spin_Correlation_Tensor::Mode()  { return m_mode; }


  // CONSTRUCTOR AND DESTRUCTOR FOR AN SCT-TREE STRUCTURE

  Spin_Correlation_Tensor::Spin_Correlation_Tensor(
      std::vector<std::pair<int,int> >* particles, 
      std::vector<Complex>* Amplitudes, 
      size_t pPos, size_t aPos1, size_t aPos2,
      size_t add):
    p_next(NULL)
  {
#ifdef SCT_Debug
    if (m_mode == scmode::None) {
      PRINT_INFO("m_mode == None; why do you want to create an SCT structure ?");
    }
#endif
    m_created = true;

    if( particles->size() ) {
      if (pPos<particles->size()) {
        // Store particle index and create follow-up nodes.
        m_particle  = (*particles)[pPos].first;
        int nrspins = (*particles)[pPos].second + 1;      // 2*spin + 1

        if (m_mode == scmode::Diagonal) {
          // only store diagonals
          p_next = new std::vector<Spin_Correlation_Tensor*>(nrspins);
          for( int leg=0; leg<nrspins; ++leg ) {
            (*p_next)[leg] = 
              new Spin_Correlation_Tensor(particles,Amplitudes,pPos+1,
                  aPos1+add*leg, aPos2+add*leg,add*nrspins);
          }	 
        }
        if (m_mode == scmode::Full) { 
          // store off-diagonals, too
          p_next = new std::vector<Spin_Correlation_Tensor*>(nrspins*nrspins);
          for( int leg1=0; leg1<nrspins; ++leg1 ) {
            for( int leg2=0; leg2<nrspins; ++leg2 ) {
//              std::cout<<om::blue<<m_particle<<" ("<<nrspins<<") : "<<aPos1+add*leg1<<" "<<aPos2+add*leg2<<"    "<<old_add<<" "<<add<<om::reset<<std::endl;
              (*p_next)[nrspins*leg1+leg2] = 
                new Spin_Correlation_Tensor(particles,Amplitudes,pPos+1,
                    aPos1+add*leg1, aPos2+add*leg2,add*nrspins);
            }
          }
        }
      } else {
        // Set the values MM*
        m_particle = -1;
//        std::cout<<aPos1<<" "<<aPos2<<"   "<<Amplitudes->size()<<std::endl;
        if( aPos1>=Amplitudes->size() || aPos2>=Amplitudes->size() ) {
          msg_Error()<<"ERROR in Spin_Correlation_Tensor constructor : "<<std::endl
            <<"     Tried to access an element of the amplitude tensor that does not exist."<<std::endl
            <<"     Don't know, what to do. Will abort."<<std::endl;
          abort();
        }
        m_value = (*Amplitudes)[aPos1] * conj((*Amplitudes)[aPos2]);
      }
    }
    else {
      m_particle  = -1;
      m_value = ( Amplitudes->size() )? (*Amplitudes)[0] : 1.;
    }
  }

  Spin_Correlation_Tensor::~Spin_Correlation_Tensor() 
  {   
    if (p_next!=NULL) {
      for (std::vector<Spin_Correlation_Tensor*>::iterator sit=p_next->begin();
          sit!=p_next->end();++sit) {
        delete (*sit);
      }
      delete p_next;
    }
  }


  // ACCESS METHODS FOR THE SCT-TREE STRUCTURE

  Complex Spin_Correlation_Tensor::Trace( Spin_Density_Matrix * sigma0 )
  {
#ifdef SCT_DEBUG
    if (m_mode == scmode::None) {
      PRINT_INFO("Warning! Spin_Correlation_Tensor::Trace(Spin_Density_Matrix) was called"
		 <<" but m_mode == None");
    }
#endif
    if (m_particle==-1) return m_value;

    if (m_particle != 0) { 
      Complex val(0., 0.);
      if (m_mode == scmode::Diagonal) {
	for (size_t i=0; i<p_next->size(); ++i) val += (*p_next)[i]->Trace();
      }
      else if (m_mode == scmode::Full) {
	size_t max(GetIdxRange());
	size_t pos(0);
	Complex val(0., 0.);
	for (size_t i=0; i<max; ++i) {
	    val += (*p_next)[pos]->Trace();
	    pos += max+1;
	}
      }
      return val;
    } 
    
    // soft-contraction with sigma0.
#ifdef SCT_Debug
    if ((m_mode == scmode::Full) && (p_next->size() != sigma0->NrEntries())) {
      PRINT_INFO("Error in Spin_Correlation_Tensor::Trace():");
      PRINT_INFO("The size of legs does not coincide with the size of the given "
		 <<"spin-density matrix. Abort the run.");
      abort();
    }
#endif
    Complex sum (0.,0.);
    if (m_mode == scmode::Diagonal)
      for (size_t i=0; i<p_next->size(); ++i)
	sum += (*sigma0)[i*p_next->size()+i] * (*p_next)[i]->Trace();
    else if (m_mode == scmode::Full)
      for (size_t i=0; i<p_next->size(); ++i)
	sum += (*sigma0)[i] * (*p_next)[i]->Trace();
    return sum;  
  }

  Complex Spin_Correlation_Tensor::Trace()
  {
    if (m_particle==-1) return m_value;
    else {
      size_t max(GetIdxRange());
      size_t pos(0);

      Complex val(0., 0.);
      if (m_mode == scmode::Full)
	for (size_t i=0; i<max; ++i) {
	  val += (*p_next)[pos]->Trace();
	  pos += max+1;
	}
      if (m_mode == scmode::Diagonal)
	for (size_t i=0; i<p_next->size(); ++i) val += (*p_next)[i]->Trace();
	  
      return val;
    }
  }

  // return SDM by using a "soft contraction" over mother SDM
  Spin_Density_Matrix Spin_Correlation_Tensor::GetSigma(int i, Spin_Density_Matrix * sigma0 )
  {
#ifdef SCT_Debug
    if (m_mode == 0)
      PRINT_INFO("And how do you want to GetSigma from nothing? m_mode==0 !");
    if (m_particle==-1 ) {
      PRINT_INFO("A terrible accident happened. The requested index "<<i<<" couldn't be found."
          <<"Let's return nothing, then.");
      return Spin_Density_Matrix();
    }
#endif
    size_t max(GetIdxRange());
    Spin_Density_Matrix anSDM(max);

    if (m_particle == i) {
      if (m_mode == scmode::Diagonal)
	for (size_t idx=0; idx<p_next->size(); ++idx)
	  anSDM[idx*p_next->size()+idx] = (*p_next)[idx]->Trace();
      if (m_mode == scmode::Full)
	for (size_t idx=0; idx<p_next->size(); ++idx)
	  anSDM[idx] = (*p_next)[idx]->Trace();
      return anSDM;
    } 
    else {
      if( m_particle != 0 ) {
	if (m_mode == scmode::Diagonal)
	  for (size_t idx=0; idx<p_next->size(); ++idx)
	    anSDM += (*p_next)[idx]->GetSigma(i);
	if (m_mode == scmode::Full)
	  for (size_t idx=0; idx<max; ++idx)
	    anSDM += (*p_next)[idx*max+idx]->GetSigma(i);
        return anSDM;
      }
      else {
#ifdef SCT_Debug
        if (p_next->size() != sigma0->NrEntries()) {
          PRINT_INFO("Error in Spin_Correlation_Tensor::Trace():");
          PRINT_INFO("The size of legs does not coincide with the size of the given "
              <<"spin-density matrix. Abort the run.");
          abort();
        }
#endif
        Spin_Density_Matrix loc_sdm(max);
	if (m_mode == scmode::Full)
	  for (size_t idx=0; idx<p_next->size(); ++idx) {
	    loc_sdm = (*p_next)[idx]->GetSigma(i);
	    loc_sdm *= (*sigma0)[idx];
	    anSDM += loc_sdm; 
	  }
	if (m_mode == scmode::Diagonal)
	  for (size_t idx=0; idx<p_next->size(); ++idx) {
	    loc_sdm = (*p_next)[idx]->GetSigma(i);
	    loc_sdm *= (*sigma0)[idx*p_next->size()+idx];
	    anSDM += loc_sdm;
	  }
      }
    }
    return anSDM;
  }
  
  Spin_Density_Matrix Spin_Correlation_Tensor::GetSigma(int i) 
  {
    if (m_particle==-1) {
      msg_Error()<<"Spin_Correlation_Tensor::GetSigma: The requested index "<<i<<" couldn't be found."
		 <<"Let's return nothing, then."<<std::endl;
      return Spin_Density_Matrix();
    }
    size_t max(GetIdxRange());
    Spin_Density_Matrix anSDM(max);

    if (m_particle == i) {

      if (m_mode==scmode::Diagonal)
	for (size_t idx=0; idx<p_next->size(); ++idx)
	  anSDM[idx*p_next->size()+idx] = (*p_next)[idx]->Trace();

      if (m_mode==scmode::Full)
	for (size_t idx=0; idx<p_next->size(); ++idx) anSDM[idx] = (*p_next)[idx]->Trace();
      
      return anSDM;
    }

    if (m_mode == scmode::Diagonal)
      for (size_t idx=0; idx<p_next->size(); ++idx)
	anSDM += (*p_next)[idx]->GetSigma(i);
    if (m_mode == scmode::Full)
      for (size_t idx=0; idx<max; ++idx)
	anSDM += (*p_next)[idx*max+idx]->GetSigma(i);

    return anSDM;
  }

  void Spin_Correlation_Tensor::Contract(int i, Spin_Density_Matrix* SDM) 
  {
    if (m_particle == -1) {
      msg_Error()<<"Spin_Correlation_Tensor::Contract: The index ("<<i<<") could not be found!"<<std::endl;
      return;
    }
    if (m_particle == i) {
      // Do matrix multiplication
      if( SDM ) {
        if (m_mode == scmode::Full)
          for (size_t idx=0; idx<p_next->size(); ++idx) 
            *(*p_next)[idx] *= (*SDM)[idx];      
        if (m_mode == scmode::Diagonal)
          for (size_t idx=0; idx<p_next->size(); ++idx)
            *(*p_next)[idx] *= (*SDM)[idx*p_next->size()+idx];
      }
      else { // multiply with unit-matrix; this is obsolete for m_mode == scmode::Diagonal
        if (m_mode == scmode::Full) {
          Spin_Density_Matrix * unitmatrix = new Spin_Density_Matrix(GetIdxRange());
          unitmatrix->SetUnitMatrix();
          for (size_t idx=0; idx<p_next->size(); ++idx) {
            *(*p_next)[idx] *= (*unitmatrix)[idx];      
          }
          delete unitmatrix;
        }
      }

      // An sum it up and delete obsolete branches
      for (size_t idx=1; idx<p_next->size();++idx) {
        *(*p_next)[0] += *(*p_next)[idx];
        delete (*p_next)[idx];
      }      

      // Copy A1 to A and soft_delete A1
      Spin_Correlation_Tensor* dummy = (*p_next)[0];
      m_particle = dummy->m_particle;
      m_value = dummy->m_value;
      delete p_next;
      p_next = dummy->p_next;
      dummy->Soft_Delete();
    } else {
      for (size_t idx=0; idx<p_next->size(); ++idx)
        (*p_next)[idx]->Contract(i, SDM);      
    }
  }


  // INTERNAL HELPER METHODS

  void Spin_Correlation_Tensor::Soft_Delete()
  {
    p_next = NULL;
    delete this;
  }
	            
  Spin_Correlation_Tensor& Spin_Correlation_Tensor::operator+=(Spin_Correlation_Tensor& SCT) 
  {
    if (m_particle==-1) 
      m_value+=SCT.m_value;
    else 
      for (size_t idx=0; idx<p_next->size(); ++idx) *(*p_next)[idx] += *SCT(idx);
    return *this;
  } 
  
  Spin_Correlation_Tensor& Spin_Correlation_Tensor::operator*=(Complex& c)
  {
    if (m_particle==-1)  m_value*=c;
    else for (size_t idx=0; idx<p_next->size(); ++idx) *(*p_next)[idx] *= c;
    return *this;
  }

  Spin_Correlation_Tensor* Spin_Correlation_Tensor::operator()(const size_t &number)
  {
      return (*p_next)[number]; 
  } 

  std::ostream& operator<<(std::ostream &ostr, Spin_Correlation_Tensor &sct)
  {
    if (sct.m_particle==-1) ostr<<sct.m_value;
    else {
      ostr<<sct.m_particle<<":{";
      for (size_t idx=0; idx<sct.p_next->size()-1; ++idx)
        ostr<<*sct(idx)<<",";   
      ostr<<*sct(sct.p_next->size()-1)<<"}";   
    }
    return ostr;
  }

  size_t Spin_Correlation_Tensor::GetDepth( size_t i )
  {
    if( m_particle == -1 ) return i;
    else return (*p_next)[0]->GetDepth(i+1);
  }

  size_t Spin_Correlation_Tensor::GetIdxRange()
  { if (m_mode == scmode::Diagonal) return p_next->size();
    if (p_next == NULL) return 0;
    switch (p_next->size()) 
      {
      case  1: return 1;
      case  4: return 2;
      case  9: return 3;
      case 16: return 4;
      case 25: return 5;
      case 36: return 6;
      default: {
	msg_Error()<<"Error in Spin_Correlation_Tensor::GetIdxRange():"<<std::endl;
	msg_Error()<<"Can't assign size of legs ("<<p_next->size()<<") to a reasonable "
		   <<"index range. Return 0"<<std::endl;
	return 0;}
      }
  }


template SP(Spin_Correlation_Tensor) &Blob_Data_Base::Get<SP(Spin_Correlation_Tensor) >();
template <> Blob_Data<SP(Spin_Correlation_Tensor) >::~Blob_Data() { }
template class Blob_Data<SP(Spin_Correlation_Tensor) >;

} //end of namespace ATOOLS
