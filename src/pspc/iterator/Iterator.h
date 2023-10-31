#ifndef PSPC_ITERATOR_H
#define PSPC_ITERATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/FSArray.h>
#include <util/global.h>                  

namespace Pscf {
namespace Pspc
{

   template <int D>
   class System;

   using namespace Util;

   /**
   * Base class for iterative solvers for SCF equations.
   *
   * \ingroup Pspc_Iterator_Module
   */
   template <int D>
   class Iterator : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      Iterator();

      /**
      * Constructor.
      * 
      * \param system parent System object
      */
      Iterator(System<D>& system);

      /**
      * Destructor.
      */
      ~Iterator();

      /**
      * Iterate to solution.
      *
      * \param isContinuation true iff a continuation within a sweep
      * \return error code: 0 for success, 1 for failure.
      */
      virtual int solve(bool isContinuation) = 0;
      
      /**
      * Log output timing results 
      */
      virtual void outputTimers(std::ostream& out) = 0;
      
      /**
      * Clear timers 
      */
      virtual void clearTimers() = 0;

      /**
      * Does this iterator use a symmetry-adapted Fourier basis?
      */
      bool isSymmetric() const
      {  return (isSymmetric_); }

      /**
      * Return true iff unit cell has any flexible lattice parameters.
      */
      bool isFlexible() const
      {  return (isFlexible_); }

      /**
      * Get the array indicating which lattice parameters are flexible.
      *
      * This array should be nParameters long, where the i-th entry is a 
      * boolean indicating whether parameter i is flexible. 
      */
      FSArray<bool,6> flexibleParams() const
      {  return flexibleParams_; }

      /**
      * Get the number of flexible lattice parameters.
      */
      int nFlexibleParams() const;

   protected:

      /**
      * Does this iterator use a symmetry-adapted basis?
      */
      bool isSymmetric_;
    
      /** 
      * Are any lattice parameters flexible during iteration?
      */
      bool isFlexible_;

      /**
      * Get parent system by const reference.
      */
      System<D> const & system() const
      {  return *sysPtr_; }

      /**
      * Get parent system by non-const reference.
      */
      System<D>& system() 
      {  return *sysPtr_; }

      /**
      * Set the array indicating which lattice parameters are flexible.
      *
      * \param flexParams array of indices of flexible lattice parameters
      */ 
      void setFlexibleParams(FSArray<bool,6> const & flexParams);

      /**
      * Array of indices of the lattice parameters that are flexible.
      */
      FSArray<bool,6> flexibleParams_;

   private:

      /// Pointer to the associated system object.
      System<D>* sysPtr_;
      
   };

   // Inline member functions

   // Default constructor
   template <int D>
   inline Iterator<D>::Iterator()
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(0)
   {  setClassName("Iterator"); }

   // Constructor
   template <int D>
   Iterator<D>::Iterator(System<D>& system)
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(&system)
   {  setClassName("Iterator"); }

   // Destructor
   template <int D>
   Iterator<D>::~Iterator()
   {}

   // Get the number of flexible lattice parameters
   template <int D>
   int Iterator<D>::nFlexibleParams() const
   {
      UTIL_CHECK(flexibleParams_.size()==system().unitCell().nParameter());
      int nFlexParams = 0;
      for (int i = 0; i < flexibleParams_.size(); i++) {
         if (flexibleParams_[i]) nFlexParams++;
      }
      return nFlexParams;
   }

   // Set the array indicating which lattice parameters are flexible.
   template <int D>
   void Iterator<D>::setFlexibleParams(FSArray<bool,6> const & flexParams)
   {  
      flexibleParams_ = flexParams; 
      if (nFlexibleParams() == 0) {
         isFlexible_ = false;
      } else {
         isFlexible_ = true;
      }
   }

} // namespace Pspc
} // namespace Pscf
#endif
