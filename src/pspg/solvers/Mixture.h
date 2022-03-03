#ifndef PSPG_MIXTURE_H
#define PSPG_MIXTURE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"
#include "Solvent.h"
#include <pscf/solvers/MixtureTmpl.h>
#include <pscf/inter/Interaction.h>
#include <util/containers/DArray.h>

namespace Pscf { 
   template <int D> class Mesh; 
}
 
namespace Pscf {
namespace Pspg
{


   /**
   * Solver for a mixture of polymers and solvents.
   *
   * A Mixture contains a list of Polymer and Solvent objects. Each
   * such object can solve the single-molecule statistical mechanics 
   * problem for an ideal gas of the associated species in a set of
   * specified chemical potential fields, and thereby compute 
   * concentrations and single-molecule partition functions. A
   * Mixture is thus both a chemistry descriptor and an ideal-gas 
   * solver.
   *
   * A Mixture is associated with a Mesh<D> object, which models a
   * spatial discretization mesh. 
   *
   * \ingroup Pspg_Solvers_Module
   */
   template <int D>
   class Mixture : public MixtureTmpl< Polymer<D>, Solvent<D> >
   {

   public:

      // Public typedefs

      /**
      * Monomer chemical potential field type.
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Monomer concentration or volume fraction field type.
      */
      typedef typename Propagator<D>::CField CField;

      // Public member functions

      /**
      * Constructor.
      */
      Mixture();

      /**
      * Destructor.
      */
      ~Mixture();

      /**
      * Read all parameters and initialize.
      *
      * This function reads in a complete description of
      * the chemical composition and structure of all species,
      * as well as the target contour length step size ds.
      *
      * \param in input parameter stream
      */
      void readParameters(std::istream& in);

      /**
      * Create an association with the mesh and allocate memory.
      * 
      * The Mesh<D> object must have already been initialized, 
      * e.g., by reading its parameters from a file, so that the
      * mesh dimensions are known on entry.
      *
      * \param mesh associated Mesh<D> object (stores address).
      */
      void setMesh(Mesh<D> const & mesh);

      /**
      * Set unit cell parameters used in solver.
      * 
      * \param unitCell UnitCell<D> object that contains Bravais lattice.
      * \param waveList WaveList<D> object that contains wavevector data
      */
      void setupUnitCell(UnitCell<D> const & unitCell, 
                         WaveList<D> const & wavelist);

      /**
      * Reset statistical segment length for one monomer type.
      * 
      * This function resets the kuhn or statistical segment length value
      * for a monomer type, and updates the associcated value in every 
      * block of that monomer type.
      *
      * \param monomerId  monomer type id
      * \param kuhn  new value for the statistical segment length
      */
      void setKuhn(int monomerId, double kuhn);

      /**
      * Compute concentrations.
      *
      * This function calls the compute function of every molecular
      * species, and then adds the resulting block concentration
      * fields for blocks of each type to compute a total monomer
      * concentration (or volume fraction) for each monomer type.
      * Upon return, values are set for volume fraction and chemical 
      * potential (mu) members of each species, and for the 
      * concentration fields for each Block and Solvent. The total
      * concentration for each monomer type is returned in the
      * cFields output parameter. 
      *
      * The arrays wFields and cFields must each have size nMonomer(),
      * and contain fields that are indexed by monomer type index. 
      *
      * \param wFields array of chemical potential fields (input)
      * \param cFields array of monomer concentration fields (output)
      */
      void 
      compute(DArray<WField> const & wFields, DArray<CField>& cFields);

      /**
      * Get monomer reference volume.
      */
      void computeStress(WaveList<D>& wavelist);

      /**
      * Get derivative of free energy w/ respect to a unit cell parameter.
      *
      * Get the pre-computed derivative with respect to unit cell 
      * parameter number n of the free energy per monomer (i.e., of the 
      * product of the free energy density and the monomer reference 
      * volume). The returned value is precomputed by the computeStress()
      * function.
      *
      * \param n  index of unit cell parameter
      */
      double stress(int n) const;


      /**
      * Get monomer reference volume.
      */
      double vMonomer() const;

      /**
      * Determine if the mixture is canonical. Otherwise, it is
      * mixed or open.
      */
      bool isCanonical();
	  
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::nMonomer;
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::nPolymer;
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::nSolvent;
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::polymer;
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::monomer;
      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::solvent;

   protected:

      using MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::setClassName;
      using ParamComposite::read;
      using ParamComposite::readOptional;

   private:

      /// Derivatives of free energy w/ respect to cell parameters.
      FArray<double, 6> stress_;
    
      /// Monomer reference volume (set to 1.0 by default).
      double vMonomer_;

      /// Optimal contour length step size.
      double ds_;

      /// Number of unit cell parameters.
      int nParams_;

      /// Pointer to associated Mesh<D> object.
      Mesh<D> const * meshPtr_;

      /// Return associated domain by reference.
      Mesh<D> const & mesh() const;

      /// Has stress been computed for current w fields?
      bool hasStress_;

   };

   // Inline member function

   /*
   * Get monomer reference volume (public).
   */
   template <int D>
   inline double Mixture<D>::vMonomer() const
   {  return vMonomer_; }

   /*
   * Get Mesh<D> by constant reference (private).
   */
   template <int D>
   inline Mesh<D> const & Mixture<D>::mesh() const
   {   
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   // Stress with respect to unit cell parameter n.
   template <int D>
   inline double Mixture<D>::stress(int n) const
   {
      UTIL_CHECK(hasStress_);  
      return stress_[n]; 
   }

   #ifndef PSPG_MIXTURE_TPP
   // Suppress implicit instantiation
   extern template class Mixture<1>;
   extern template class Mixture<2>;
   extern template class Mixture<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
//#include "Mixture.tpp"
#endif
