#ifndef ATOOLS_Org_List_Algorithms_H
#define ATOOLS_Org_List_Algorithms_H

namespace ATOOLS {

    // To create a list from an iterator, use: (i.e. for a list of particles);
    // #include <algorithm>
    //     list<Particle*> thelist;
    //     copy( evt->particles_begin(), evt->particles_end(), 
    //           back_inserter(thelist) );
    // to create a list subject to a condition (predicate) use:
    //     list<GenParticle*> thelist;
    //     copy_if( evt->particles_begin(), evt->particles_end(), 
    //                     back_inserter(thelist), is_photon() );
    // where is_photon() is a predicate like:
    //     class is_photon {
    //       public:
    //         bool operator() ( const Particle* p ) {
    //             if ( p && p->pdg_id() == 22 ) return 1;
    //             return 0;
    //         }
    //     };
    // which the user defines herself, since it was somehow dropped from
    // the set of algorithms provided by the standard library:
    template <class InIter, class OutIter, class Predicate>
    void copy_if( InIter first, InIter last, OutIter out,
		  Predicate pred ) {
	for ( ; first != last; ++first ) { if ( pred(*first) ) out = *first; }
    }


  /*!
    \file
    \brief contains the copy_if template (somehow  missing in STL) 
  */

  /*!
    \fn copy_if ( InIter first, InIter last, OutIter out, Predicate pred )
    \brief copy_if is provided since it is somehow missing in STL
  */
}
#endif
