#ifndef PSPG_MIXTURE_TPP
#define PSPG_MIXTURE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include <pspg/math/GpuResources.h>

#include <cmath>

namespace Pscf { 
namespace Pspg
{

   template <int D>
   Mixture<D>::Mixture()
    : vMonomer_(1.0),
      ds_(-1.0),
      meshPtr_(0),
      hasStress_(false)
   {  setClassName("Mixture"); }

   template <int D>
   Mixture<D>::~Mixture()
   {}

   template <int D>
   void Mixture<D>::readParameters(std::istream& in)
   {
      MixtureTmpl< Pscf::Pspg::Polymer<D>, Pscf::Pspg::Solvent<D> >::readParameters(in);
      vMonomer_ = 1.0; // Default value
      readOptional(in, "vMonomer", vMonomer_);
      read(in, "ds", ds_);

      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);
   }

   template <int D>
   void Mixture<D>::setMesh(Mesh<D> const& mesh)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);

      meshPtr_ = &mesh;

      // Set discretization for all blocks
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).setDiscretization(ds_, mesh);
         }
      }

      // Set spatial discretization for solvents
      if (nSolvent() > 0) {
         for (int i = 0; i < nSolvent(); ++i) {
            solvent(i).setDiscretization(mesh);
         }
      }
   }

   template <int D>
   void Mixture<D>::setupUnitCell(const UnitCell<D>& unitCell, const WaveList<D>& wavelist)
   {
      nParams_ = unitCell.nParameter();
      for (int i = 0; i < nPolymer(); ++i) {
         polymer(i).setupUnitCell(unitCell, wavelist);
      }

      hasStress_ = false;
   }

      /*
   * Reset statistical segment length for one monomer type.
   */
   template <int D>
   void Mixture<D>::setKuhn(int monomerId, double kuhn)
   {
      // Set new Kuhn length for relevant Monomer object
      monomer(monomerId).setKuhn(kuhn);

      // Update kuhn length for all blocks of this monomer type
      for (int i = 0; i < nPolymer(); ++i) {
         for (int j =  0; j < polymer(i).nBlock(); ++j) {
            if (monomerId == polymer(i).block(j).monomerId()) {
               polymer(i).block(j).setKuhn(kuhn);
            }
         }
      }
      hasStress_ = false;
   }

   /*
   * Compute concentrations (but not total free energy).
   */
   template <int D>
   void Mixture<D>::compute(DArray<Mixture<D>::WField> const & wFields, 
                            DArray<Mixture<D>::CField>& cFields)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(mesh().size() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);
      UTIL_CHECK(wFields.capacity() == nMonomer());
      UTIL_CHECK(cFields.capacity() == nMonomer());

      int nMesh = mesh().size();
      int nm = nMonomer();
      int i, j;

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nMesh, nBlocks, nThreads);

      // Clear all monomer concentration fields
      for (i = 0; i < nm; ++i) {
         UTIL_CHECK(cFields[i].capacity() == nMesh);
         UTIL_CHECK(wFields[i].capacity() == nMesh);
         assignUniformReal<<<nBlocks, nThreads>>>(cFields[i].cDField(), 0.0, nMesh);
      }


      int monomerId;
      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields);
      }

      // Accumulate monomer concentration fields
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nm);
            CField& monomerField = cFields[monomerId];
            CField& blockField = polymer(i).block(j).cField();
            UTIL_CHECK(blockField.capacity()==nMesh);
            pointWiseAdd<<<nBlocks, nThreads>>>(monomerField.cDField(), blockField.cDField(), nMesh);
         }
      }
      
      // Process solvent species
      for (i = 0; i < nSolvent(); i++) {
         monomerId = solvent(i).monomerId();
         UTIL_CHECK(monomerId >= 0);
         UTIL_CHECK(monomerId < nm);

         // Compute solvent concentration
         solvent(i).compute(wFields[monomerId]);

         // Add solvent contribution to relevant monomer concentration
         CField& monomerField = cFields[monomerId];
         CField const & solventField = solvent(i).concField();
         UTIL_CHECK(solventField.capacity() == nMesh);
         pointWiseAdd<<<nBlocks, nThreads>>>(monomerField.cDField(), solventField.cDField(), nMesh);


      }

      hasStress_ = false;
   }

   /*  
   * Compute Total Stress.
   */  
   template <int D>
   void Mixture<D>::computeStress(WaveList<D>& wavelist)
   {   
      int i, j;

      // Compute stress for each polymer.
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).computeStress(wavelist);
      } 

      // Accumulate total stress 
      for (i = 0; i < nParams_; ++i) {
         stress_[i] = 0.0;
         for (j = 0; j < nPolymer(); ++j) {
            stress_[i] += polymer(j).stress(i);
         }   
      }

      // Note: Solvent does not contribute to derivatives of f_Helmholtz
      // with respect to unit cell parameters at fixed volume fractions.

      hasStress_ = true;
   }

   template <int D>
   bool Mixture<D>::isCanonical()
   {
      // Check ensemble of all polymers
      for (int i = 0; i < nPolymer(); ++i) {
         if (polymer(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Check ensemble of all solvents
      for (int i = 0; i < nSolvent(); ++i) {
         if (solvent(i).ensemble() == Species::Open) {
            return false;
         }
      }
      // Returns true if false was never returned
      return true;
   }

} // namespace Pspg
} // namespace Pscf
#endif
