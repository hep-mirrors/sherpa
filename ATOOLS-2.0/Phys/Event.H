// former GenEvent.H from HepMC
#ifndef EVENT_H
#define EVENT_H

#include <set>
#include <vector>
#include <algorithm>
#include <iostream>

#include "Blob.H"
#include "Particle.H"

namespace ATOOLS {

  class Event {

  public:
    Event( int signal_process_id = 0, int event_number = 0,
	   Blob* signal_vertex = 0,
	   const std::vector<long int>& randomstates
	   = std::vector<long int>() );
    Event( const Event& inevent );          // deep copy
    Event& operator=( const Event& inevent ); // deep.
    virtual ~Event(); //deletes all vertices/particles in this evt
    
    void print( std::ostream& ostr = std::cout ) const; // dumps to ostr

    ////////////////////
    // access methods //
    ////////////////////

    int signal_process_id() const;
    int event_number() const;
    Blob* signal_process_vertex() const;

    std::vector<long int> random_states() const;

    void set_signal_process_id( int id );
    void set_event_number( int eventno );
    void set_signal_process_vertex( Blob* );
    void set_random_states( const std::vector<long int>& randomstates );

    // following in temporary ... don't use it in your code!
    static unsigned int counter(); //num Event objects in memory
   
  public:
    ///////////////////////////////
    // vertex_iterator           //
    ///////////////////////////////
	
    typedef std::set<Blob*>::const_iterator vertex_const_iterator;
    typedef std::set<Blob*>::iterator vertex_iterator;
    int     vertices_size() const;
    bool    vertices_empty() const;
    bool    add_vertex( Blob* vtx );    // adds to evt and adopts
    bool    remove_vertex( Blob* vtx ); // erases vtx from evt, 
    // does not delete
    vertex_iterator            vertices_begin();
    vertex_iterator            vertices_end();
    vertex_const_iterator      vertices_begin() const;
    vertex_const_iterator      vertices_end() const;

  public:
    ///////////////////////////////
    // particle_iterator         //
    ///////////////////////////////
    class particle_iterator :
      public std::forward_iterator<Blob*,std::ptrdiff_t>{
      // Iterates over all particles in this event
    public:
      particle_iterator();
      particle_iterator( const Event& evt_root );
      particle_iterator( const particle_iterator& );
      virtual             ~particle_iterator();
      particle_iterator&  operator=( const particle_iterator& );
      Particle*        operator*(void) const;
      particle_iterator&  operator++(void);  //Pre-fix increment 
      particle_iterator   operator++(int);   //Post-fix increment
      bool                operator==( const particle_iterator& ) const;
      bool                operator!=( const particle_iterator& ) const;
    protected:
      Particle*        advance_to_first_();
    private:
      const Event*           m_event;
      Event::vertex_iterator m_vertex_iterator;
      Blob::edge_iterator  m_edge;// points to return particle
    };	
    particle_iterator       particles_begin() const;
    particle_iterator       particles_end() const;

    ////////////////////////////////////////////////

  protected:
    void delete_all_vertices();

  private: // data members
    int                   m_signal_process_id;
    int                   m_event_number;  
    Blob*            m_signal_process_vertex;
    std::set<Blob*>  m_vertices;// container of vertices
    // is used by default for hit and miss
    std::vector<long int> m_random_states; // container of rndm num 
    // generator states
    static unsigned int   s_counter;
  };

  ///////////////////////////
  // INLINE Access Methods //
  ///////////////////////////

  inline int Event::signal_process_id() const 
  { return m_signal_process_id; }

  inline int Event::event_number() const 
  { return m_event_number; }
 
  inline Blob* Event::signal_process_vertex() const {
    // returns a (mutable) pointer to the signal process vertex
    return m_signal_process_vertex;
  }  

  inline std::vector<long int> Event::random_states() const 
  { return m_random_states; }

  inline void Event::set_signal_process_id( int id )
  { m_signal_process_id = id; }

  inline void Event::set_event_number( int eventno )
  { m_event_number = eventno; }

  inline void Event::set_signal_process_vertex( Blob* vtx ) {
    m_signal_process_vertex = vtx;
    if ( m_signal_process_vertex ) add_vertex( m_signal_process_vertex );
  }

  inline void Event::set_random_states( const std::vector<long int>&
					randomstates )
  { m_random_states = randomstates; }

  //////////////////////////////
  // INLINES vertex_iterator  //
  //////////////////////////////

  inline int Event::vertices_size() const {
    return (int)m_vertices.size();
  }
  inline bool Event::vertices_empty() const {
    return (bool)m_vertices.empty();
  }
  inline bool Event::add_vertex( Blob* vtx ) {
    // returns true if successful - generally will only return false
    // if the inserted vertex is already included in the event.
    if ( !vtx ) return 0;
    return m_vertices.insert( vtx ).second;
  }
  inline bool Event::remove_vertex( Blob* vtx ) {
    // this removes vtx from the event but does NOT delete it.
    // returns True is an entry vtx existed in the table and was erased
    if ( m_signal_process_vertex == vtx ) m_signal_process_vertex = 0;
    return (bool)m_vertices.erase( vtx );
  }
  inline Event::vertex_iterator Event::vertices_begin() { 
    return m_vertices.begin(); 
  }
  inline Event::vertex_iterator Event::vertices_end() { 
    return m_vertices.end(); 
  }
  inline Event::vertex_const_iterator Event::vertices_begin() const { 
    return m_vertices.begin(); 
  }
  inline Event::vertex_const_iterator Event::vertices_end() const { 
    return m_vertices.end(); 
  }

  ///////////////////////////////
  // INLINES particle iterator //
  ///////////////////////////////

  inline bool Event::particle_iterator::operator==( 
						   const Event::particle_iterator& a ) const {
    return **this == *a; 
  }
  inline bool Event::particle_iterator::operator!=( 
						   const Event::particle_iterator& a ) const {
    return !(**this == *a); 
  }
  inline Event::particle_iterator 
  Event::particles_begin() const {
    return Event::particle_iterator( *this );
  }
  inline Event::particle_iterator 
  Event::particles_end() const {
    return Event::particle_iterator();
  }

}

#endif




