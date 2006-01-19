#include "Spin_Correlation_Tensor.H"
#include "Message.H"
#include "Blob.H"
#include "Smart_Pointer.H"

INSTANTIATE_SMART_POINTER(Spin_Correlation_Tensor);

#define SCT_Debug

namespace ATOOLS{

  int Spin_Correlation_Tensor::m_k0_n(-1);

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

  Spin_Correlation_Tensor::Spin_Correlation_Tensor(
      std::vector<std::pair<int,int> >* particles, 
      std::vector<Complex>* Amplitudes, 
      size_t pPos, size_t aPos1, size_t aPos2):
    p_next(NULL)
  {
    if( particles->size() ) {
      if (pPos<particles->size()) {
        // Store particle index and create follow-up nodes.
        m_particle  = (*particles)[pPos].first;
        int nrspins = (*particles)[pPos].second + 1;      // 2*spin + 1
        p_next = new std::vector<Spin_Correlation_Tensor*>(nrspins*nrspins);
        size_t add = pow( nrspins, pPos );
        for( int leg1=0; leg1<nrspins; ++leg1 ) {
          for( int leg2=0; leg2<nrspins; ++leg2 ) {
            (*p_next)[nrspins*leg1+leg2] = 
              new Spin_Correlation_Tensor(particles,Amplitudes,pPos+1,aPos1+add*leg1, aPos2+add*leg2);
          }
        }
      } else {
        // Set the values MM*
        m_particle = -1;
        m_value = (*Amplitudes)[aPos1] * conj((*Amplitudes)[aPos2]);
      }
    }
    else {
      m_particle  = -1;
      m_value = ( Amplitudes->size() )? (*Amplitudes)[0] : 1.;
    }
  }

  Spin_Correlation_Tensor::Spin_Correlation_Tensor(double x, std::vector<int>* creation_list, 
						   size_t pos) :
    p_next(NULL)
  {
    if (pos<creation_list->size()) {
      m_particle=(*creation_list)[pos];
      
      // Create the follow-up nodes
      p_next = new std::vector<Spin_Correlation_Tensor*>(4);
      (*p_next)[0] = new Spin_Correlation_Tensor(x/2, creation_list, pos+1);
      (*p_next)[1] = new Spin_Correlation_Tensor(0., creation_list, pos+1);
      (*p_next)[2] = new Spin_Correlation_Tensor(0., creation_list, pos+1);
      (*p_next)[3] = new Spin_Correlation_Tensor(x/2, creation_list, pos+1);
    }
    else {
      m_particle = -1;
      m_value    =Complex(x, 0.);
    };
  }

  Complex Spin_Correlation_Tensor::Trace( Spin_Density_Matrix * sigma0 )
  {
    if (m_particle==-1) return m_value;
    if (m_particle != 0) {
      size_t max(GetIdxRange());
      size_t pos(0);
      Complex val(0., 0.);
      for (size_t i=0; i<max; ++i) {
	val += (*p_next)[pos]->Trace();
	pos += max+1;
      }
      return val;
    } // return (*p_next)[0]->Trace(sigma0) + (*p_next)[3]->Trace(sigma0); 

    // soft-contraction with sigma0.
#ifdef SCT_Debug
    if (p_next->size() != sigma0->NrEntries()) {
      PRINT_INFO("Error in Spin_Correlation_Tensor::Trace():");
      PRINT_INFO("The size of legs does not coincide with the size of the given "
		 <<"spin-density matrix. Abort the run.");
      abort();
    }
#endif
    Complex sum (0.,0.);
    for (size_t i=0; i<p_next->size(); ++i)
      sum += (*sigma0)[i] * (*p_next)[0]->Trace();
    return sum;  
  }

  // return SDM by using a "soft contraction" over mother SDM
  Spin_Density_Matrix Spin_Correlation_Tensor::GetSigma(int i, Spin_Density_Matrix * sigma0 ) 
  {
    if (m_particle==-1 ) {
      PRINT_INFO("A terrible accident happened. The requested index couldn't be found."
          <<"Let's return nothing, then.");
      return Spin_Density_Matrix();
    }
    size_t max(GetIdxRange());
    Spin_Density_Matrix anSDM(max);

    if (m_particle == i) {
      for (size_t idx=0; idx<p_next->size(); ++idx)
        anSDM[idx] = (*p_next)[idx]->Trace();
      return anSDM;
    } 
    else {
      if( m_particle != 0 ) {
        for (size_t idx=0; idx<p_next->size(); ++idx)
          anSDM += (*p_next)[idx]->GetSigma(i);
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
        for (size_t idx=0; idx<p_next->size(); ++idx) {
          loc_sdm = (*p_next)[idx]->GetSigma(i);
          loc_sdm *= (*sigma0)[idx];
          anSDM += loc_sdm;
        }
      }
    }
    return anSDM;
  } // returns same SDM as using "hard" contraction over mother sigma

  Complex Spin_Correlation_Tensor::Trace()
  {
    if (m_particle==-1) return m_value;
    else {
      size_t max(GetIdxRange());
      size_t pos(0);
      Complex val(0., 0.);
      for (size_t i=0; i<max; ++i) {
	val += (*p_next)[pos]->Trace();
	pos += max+1;
      }
      return val;
    } //return (*p_next)[0]->Trace() + (*p_next)[3]->Trace();		     
  }

  Spin_Density_Matrix Spin_Correlation_Tensor::GetSigma(int i) 
  {
    if (m_particle==-1) {
      PRINT_INFO("A terrible accident happened. The requested index couldn't be found."
          <<"Let's return nothing, then.");
      return Spin_Density_Matrix();
    }
    size_t max(GetIdxRange());
    Spin_Density_Matrix anSDM(max);

    if (m_particle == i) {
      for (size_t idx=0; idx<p_next->size(); ++idx)
        anSDM[idx] = (*p_next)[idx]->Trace();
      return anSDM;
    }

    for (size_t idx=0; idx<p_next->size(); ++idx)
      anSDM += (*p_next)[idx]->GetSigma(i);
    return anSDM;
  }

  void Spin_Correlation_Tensor::Contract(int i, Spin_Density_Matrix* SDM) 
  {
    if (m_particle == -1) {
      PRINT_INFO("The index to be contracted could not be found!");
      return;
    }
    if (m_particle == i) {
      // Do matrix multiplication
      for (size_t idx=0; idx<p_next->size(); ++idx)
        *(*p_next)[idx] *= (*SDM)[idx];      

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

  void Spin_Correlation_Tensor::Soft_Delete()
  {
    p_next = NULL;
    delete this;
  }
	            
  Spin_Correlation_Tensor& Spin_Correlation_Tensor::operator+=(Spin_Correlation_Tensor& SCT) 
  {
    if (m_particle==-1) {
      m_value+=SCT.m_value;
    } else {
      *(*p_next)[0] += *SCT(0);    
      *(*p_next)[1] += *SCT(1);  
      *(*p_next)[2] += *SCT(2);  
      *(*p_next)[3] += *SCT(3);
    } 
    return *this;
  }
  
  Spin_Correlation_Tensor& Spin_Correlation_Tensor::operator*=(Complex& c)
  {
    if (m_particle==-1) {
      m_value*=c;
    } else {
      *(*p_next)[0] *= c;   
      *(*p_next)[1] *= c;    
      *(*p_next)[2] *= c;    
      *(*p_next)[3] *= c;
    }
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
      ostr<<*sct(0)<<",";   
      ostr<<*sct(1)<<",";      
      ostr<<*sct(2)<<",";   
      ostr<<*sct(3)<<"}";
    }
    return ostr;
  }

  size_t Spin_Correlation_Tensor::GetDepth( size_t i )
  {
    if( m_particle == -1 ) return i;
    else return (*p_next)[0]->GetDepth(i+1);
  }

  size_t Spin_Correlation_Tensor::GetIdxRange()
  {
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
	PRINT_INFO("Error in Spin_Correlation_Tensor::GetIdxRange():");
	PRINT_INFO("Can't assign size of legs ("<<p_next->size()<<") to a reasonable "
		   <<"index range. Return 0");
	return 0;}
      }
  }

template class Blob_Data<SP(Spin_Correlation_Tensor) >;
template SP(Spin_Correlation_Tensor) &Blob_Data_Base::Get<SP(Spin_Correlation_Tensor) >();
template <> Blob_Data<SP(Spin_Correlation_Tensor) >::~Blob_Data() { }

} //end of namespace ATOOLS
