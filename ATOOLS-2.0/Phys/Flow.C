// former Flow.cxx of HepMC
#include "Flow.H"
#include "Parton.H"
//#include "Blob.H"


namespace APHYTOOLS {
//  Counter          * count = 0;
  long Flow::qcd_counter = 600;
}


using namespace APHYTOOLS;

Flow::Flow( Parton* parton_owner ) : 
  m_parton_owner(parton_owner) 
{}

Flow::Flow( const Flow& inflow ) : 
  m_parton_owner(inflow.m_parton_owner)
{
  // copies both the m_icode AND the m_parton_owner
  *this = inflow;
}

Flow::~Flow() {
  m_icode.clear();
}

void Flow::print( std::ostream& ostr ) const {
  ostr << "Flow(" << m_parton_owner << "): " << *this << std::endl;
}
    
std::set<Parton*> Flow::connected_partners( int code, int code_index, 
						 int num_indices ) const {
  // Returns all flow partners which have "code" in any  of the 
  //  num_indices beginning with index code_index.
  //  m_parton_owner is included in the result.
  //  Return is by value since the set should never be very big.
  // EXAMPLE: if you want to find all flow partners that have the same
  //   code in indices 2,3,4 as parton p has in index 2, you would use:
  //   set<Parton*> result = 
  //             p->Flow().connected_partners(p->Flow().icode(2),2,3);
  //
  std::set<Parton*> output;
  for ( int i = code_index; i!=code_index+num_indices; ++i ) {
    if ( icode(i)==code ) {
      output.insert(m_parton_owner);
      connected_partners( &output, code, code_index, num_indices );
      break;
    } 
  }
  return output;
}


void Flow::connected_partners( std::set<Parton*>* output, int code, 
			       int code_index, int num_indices ) const
{
  /*
  // protected: for recursive use by Flow::connected_partners()
  //
  if ( !m_parton_owner ) return; // nothing to do
  // look for connected partners joined to this m_parton_owner
  // through its end_vertex
  if ( m_parton_owner->end_vertex() ) {
    for ( Blob::parton_iterator p 
	    = m_parton_owner->end_vertex()->partons_begin(family);
	  p != m_parton_owner->end_vertex()->partons_end(family);
	  ++p ) {
      // if the parton has the correct flow code and is not yet in 
      // the set, then we recursively call connected_partners
      for ( int index = code_index; index!=code_index+num_indices; 
	    ++index ){
	if ( (*p)->Flow(index)==code &&
	     output->insert(*p).second ) {
	  (*p)->Flow().connected_partners( output, code,
					   code_index, 
					   num_indices );
	}
      }
    }
  }
  // same for production_vertex
  if ( m_parton_owner->production_vertex() ) {
    for ( Blob::parton_iterator p 
	    = m_parton_owner->production_vertex()->
	    partons_begin( family );
	  p != m_parton_owner->production_vertex()->
	    partons_end( family ); ++p ) {
      // if the parton has the correct flow code and is not yet in 
      // the set, then we recursively call connected_partners
      for ( int index = code_index; index!=code_index+num_indices; 
	    ++index ){
	if ( (*p)->Flow(index)==code &&
	     output->insert(*p).second ) {
	  (*p)->Flow().connected_partners( output, code,
					   code_index, 
					   num_indices );
	}
      }
    }
  }
  */
}

std::set<Parton*> Flow::dangling_connected_partners( int code, 
							  int code_index, int num_indices ) const {
  std::set<Parton*> output;
  std::set<Parton*> visited_partons;
  for ( int i = code_index; i!=code_index+num_indices; ++i ) {
    if ( icode(i)==code ) {
      visited_partons.insert(m_parton_owner);
      dangling_connected_partners( &output, &visited_partons, code,
				   code_index, num_indices );
      break;
    }
  }
  return output;
}

void Flow::dangling_connected_partners( std::set<Parton*>* output, 
					std::set<Parton*>* 
					visited_partons,
					int code, int code_index, 
					int num_indices ) const 
{
  /*
  // protected: for recursive use by Flow::dangling_connected_partners
  //
  if ( !m_parton_owner ) return; // nothing to do
  int count_partners = 0;
  // look for connected partners joined to this m_parton_owner
  // through its end_vertex
  if ( m_parton_owner->end_vertex() ) {
    for ( Blob::parton_iterator p 
	    = m_parton_owner->end_vertex()->partons_begin(family);
	  p != m_parton_owner->end_vertex()->partons_end(family);
	  ++p ) {
      // if the parton has the correct flow code and is not yet in 
      // the set, then we recursively call connected_partners
      for ( int index = code_index; index!=code_index+num_indices; 
	    ++index ){
	if ( (*p)->Flow(index)==code ) {
	  if ( *p!=m_parton_owner ) ++count_partners;
	  if ( visited_partons->insert(*p).second ) {
	    (*p)->Flow().dangling_connected_partners( output, 
						      visited_partons, code,
						      code_index, num_indices );

	  }
	}		
      }
    }
  }
  // same for production_vertex
  if ( m_parton_owner->production_vertex() ) {
    for ( Blob::parton_iterator p = m_parton_owner->
	    production_vertex()->
	    partons_begin( family );
	  p != m_parton_owner->production_vertex()->
	    partons_end( family ); 
	  ++p ) {
      // if the parton has the correct flow code and is not yet in 
      // the set, then we recursively call connected_partners
      for ( int index = code_index; index!=code_index+num_indices; 
	    ++index ){
	if ( (*p)->Flow(index)==code ) {
	  if ( *p!=m_parton_owner ) ++count_partners;
	  if ( visited_partons->insert(*p).second ) {
	    (*p)->Flow().dangling_connected_partners( output, 
						      visited_partons, code,
						      code_index, num_indices );

	  }
	}
      }
    }
  }
  if ( count_partners <= 1 ) output->insert( m_parton_owner );
  */
}
	
/////////////
// Friends //
/////////////

std::ostream& APHYTOOLS::operator<<( std::ostream& ostr, const Flow& f ) {
  ostr << f.m_icode.size();
  for ( std::map<int,int>::const_iterator i = f.m_icode.begin();
	i != f.m_icode.end(); ++i ) {
    ostr << " " << (*i).first << " " << (*i).second;
  }
  return ostr;
}


















