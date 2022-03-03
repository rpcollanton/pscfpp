#ifndef PSPG_SYSTEM_H
#define PSPG_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspg/field/FieldIo.h>            // member
#include <pspg/iterator/Iterator.h>        // member
#include <pspg/solvers/Mixture.h>          // member
#include <pspg/solvers/Solvent.h>        // member
#include <pspg/field/Domain.h>             // member
#include <pspg/solvers/WaveList.h>         // member
#include <pspg/field/RDField.h>            // typedef
#include <pspg/field/RDFieldDft.h>         // typedef

#include <pscf/crystal/Basis.h>            // member
#include <pscf/mesh/Mesh.h>                // member
#include <pscf/crystal/UnitCell.h>         // member
#include <pscf/homogeneous/Mixture.h>      // member
#include <pscf/inter/ChiInteraction.h>     // member

#include <util/param/ParamComposite.h>     // base class
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/Array.h>         // function parameter

namespace Pscf {
namespace Pspg
{
   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;

   using namespace Util;

   /**
   * Main class in SCFT simulation of one system.
   *
   * \ingroup Pscf_Pspg_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      /// Base class for WField and CField
      typedef RDField<D> Field;

      /// Monomer chemical potential field type.
      typedef typename Propagator<D>::WField WField;

      /// Monomer concentration / volume fraction field type.
      typedef typename Propagator<D>::CField CField;

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      /// \name Lifetime (Actions)
      //@{

      /**
      * Process command line options.
      */
      void setOptions(int argc, char **argv);

      /** 
      * Set number of blocks and number of threads
      */
      void setGpuResources (int nBlocks, int nThreads);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      */
      void readParam();

      /**
      * Read body of input parameters block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read command script.
      * 
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      */
      void readCommands();

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * iterator().solve(). Resulting values are returned by the 
      * freeEnergy() and pressure() accessor functions.
      */
      void computeFreeEnergy();

      /**
      * Output thermodynamic properties to a file. 
      *
      * This function outputs Helmholtz free energy per monomer,
      * pressure (in units of kT per monomer volume), and the
      * volume fraction and chemical potential of each species.
      *
      * \param out output stream 
      */
      void outputThermo(std::ostream& out);

      //@}
      /// \name Chemical Potential Fields (W Fields)
      //@{

      /**
      * Get an array of chemical potential fields, in a basis.
      *
      * This function returns an array in which each element is an
      * array containing the coefficients of the chemical potential
      * field (w field) in a symmetry-adapted basis for one monomer 
      * type. The array capacity is the number of monomer types.
      */
      DArray<RDField <D> >& wFields();

      /**
      * Get chemical potential field for one monomer type, in a basis.
      *
      * This function returns an array containing coefficients of 
      * the chemical potential field (w field) in a symmetry-adapted
      * basis for a specified monomer type.
      *
      * \param monomerId integer monomer type index
      */
      RDField<D>& wField(int monomerId);

      /**
      * Get array of chemical potential fields, on an r-space grid.
      *
      * This function returns an array in which each element is a
      * WField object containing values of the chemical potential field 
      * (w field) on a regular grid for one monomer type. The array 
      * capacity is the number of monomer types.
      */
      DArray<WField>& wFieldsRGrid();

      /**
      * Get the chemical potential field for one monomer type, on a grid.
      *
      * \param monomerId integer monomer type index
      */
      WField& wFieldRGrid(int monomerId);

      /**
      * Get array of chemical potential fields, in Fourier space.
      *
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RDFieldDft<D> >& wFieldsKGrid();

      /**
      * Get the chemical potential field for one monomer, in Fourier space.
      *
      * \param monomerId integer monomer type index
      */
      RDFieldDft<D>& wFieldKGrid(int monomerId);

      //@{
      /// \name Monomer Concentration / Volume Fraction Fields (C Fields)
      //@{
      
      /**
      * Get an array of all monomer concentration fields, in a basis
      *
      * This function returns an array in which each element is an
      * array containing the coefficients of the monomer concentration
      * field (cfield) for one monomer type in a symmetry-adapted basis.
      * The array capacity is equal to the number of monomer types.
      */
      DArray<RDField <D> >& cFields();

      /**
      * Get the concentration field for one monomer type, in a basis.
      *
      * This function returns an array containing the coefficients of 
      * the monomer concentration / volume fraction field (c field) 
      * for a specific monomer type. 
      *
      * \param monomerId integer monomer type index
      */
      RDField<D>& cField(int monomerId);

      /**
      * Get array of all concentration fields (c fields), on a grid.
      *
      * This function returns an array in which each element is the
      * monomer concentration field for one monomer type on a regular 
      * grid (an r-grid). 
      */
      DArray<CField>& cFieldsRGrid();

      /**
      * Get the concentration (c field) for one monomer type, on a grid.
      *
      * \param monomerId integer monomer type index
      */
      CField& cFieldRGrid(int monomerId);

      /**
      * Get all monomer concentration fields, in Fourier space (k-grid).
      *
      * This function returns an arrray in which each element is the
      * discrete Fourier transform (DFT) of the concentration field
      * (c field) for on monomer type.
      */
      DArray<RDFieldDft<D> >& cFieldsKGrid();

      /**
      * Get the c field for one monomer type, in Fourier space (k-grid).
      *
      * This function returns the discrete Fourier transform (DFT) of the 
      * concentration field (c field) for monomer type index monomerId.
      *
      * \param monomerId integer monomer type index
      */
      RDFieldDft<D>& cFieldKGrid(int monomerId);

      //@}
      /// \name Accessors (access objects by reference)
      //@{

      /**
      * Get Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get spatial discretization mesh by reference.
      */
      Mesh<D> & mesh();

      /**
      * Get crystal unitCell (i.e., lattice type and parameters) by reference.
      */
      UnitCell<D> & unitCell();

      /**
      * Get Domain by const reference.
      */
      Domain<D> const & domain() const;

      /**
      * Get interaction (i.e., excess free energy model) by reference.
      */
      ChiInteraction& interaction();
	  
      /**
      * Get the iterator by reference.
      */
      Iterator<D>& iterator();

      /**
      * Get basis object by reference.
      */
      Basis<D> & basis();

      /**
      * Get container for wavevector data.
      */
      WaveList<D>& wavelist();

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D> & fieldIo();

      /**
      * Get associated FFT objecti by reference.
      */
      FFT<D> & fft();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      //@}
      /// \name Accessors (return values)
      //@{
      
      /** 
      * Get the group name string.
      */  
      std::string groupName();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure x monomer volume kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double pressure() const;

      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * \param fields  array of new w (chemical potential) fields
      */  
      void setWRGrid(DArray<RDField<D>> const & fields);

      /**
      * Set an association to a new unit cell.
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set new parameters of the associated unit cell.
      *
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      /** 
      * Have monomer chemical potential fields (w fields) been set?
      *
      * A true value is returned if and only if values have been set on a 
      * real space grid. The READ_W_BASIS command must immediately convert 
      * from a basis to a grid to satisfy this requirement.
      */
      bool hasWFields() const;

      /** 
      * Have monomer concentration fields (c fields) been computed?
      *
      * A true value is returned if and only if monomer concentration fields
      * have been computed by solving the modified diffusion equation for the
      * current w fields, and values are known on a grid (cFieldsRGrid).
      */  
      bool hasCFields() const;

      /** 
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file with name "filename",
      * which must contain chemical potential fields in symmetry-adapted 
      * basis format, stores these fields in the system wFields array,
      * converts these fields to real-space grid format and stores the
      * result in the wFieldsRGrid array. On exit hasWFields is set true 
      * and hasCFields is false.
      *
      * \param filename name of input w-field basis file
      */
      void readWBasis(const std::string & filename);

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields (wFields
      * and wFieldRGrid). This requires solution of the modified 
      * diffusion equation for all polymers, computation of Boltzmann
      * weights for all solvents, computation of molecular partition
      * functions for all species, and computation of concentration
      * fields for blocks and solvents, and computation of overall 
      * concentrations for all monomer types. This function does not 
      * compute the canonical (Helmholtz) free energy or grand-canonical 
      * free energy (i.e., pressure). Upon return, the flag hasCFields 
      * is set true.
      *
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \param needStress true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Iteratively solve a SCFT problem.
      * 
      * This function calls the iterator to attempt to solve the SCFT
      * problem for the current mixture and system parameters, using
      * the current chemical potential fields (wFields and wFieldRGrid) 
      * and current unit cell parameter values as initial guesses.  
      * Upon exist, hasCFields is set true whether or not convergence 
      * is obtained to within the desired tolerance.  The Helmholtz free 
      * energy and pressure are computed if and only if convergence is
      * obtained. 
      *
      * \pre The hasWFields flag must be true on entry.
      * \return returns 0 for successful convergence, 1 for failure.
      */
      int iterate();

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the 
      * parameter file to solve the SCF problem at a sequence of points
      * along a line in parameter space. The nature of this sequence of
      * points is determined by implementation of a subclass of Sweep
      * and the parameters passed to the sweep object in the parameter 
      * file.  The Iterator that is initialized in the parameter file 
      * is called at each state point.
      *
      * An Exception is thrown if this is called when no Sweep has been 
      * created (i.e., if hasSweep_ == false).
      */
      void sweep();

      /**
      * Convert a field from symmetry-adapted basis to r-grid format.
      *
      * This function uses the arrays that stored monomer concentration 
      * fields for temporary storage, and thus corrupts any previously
      * stored values. As a result, flag hasCFields is false on output.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /** 
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * This function uses the arrays that stored monomer concentration 
      * fields for temporary storage, and thus corrupts any previously
      * stored values. As a result, flag hasCFields is false on return.
      *
      * \param inFileName name of input file
      * \param outFileName name of output file
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Write chemical potential fields in symmetry adapted basis format.
      *
      * \param filename name of output file
      */
      void writeWBasis(const std::string & filename);

      /**
      * Write concentrations in symmetry-adapted basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename);

      /**
      * Write last contour slice of the propagator in real space grid format.
      *
      * \param filename name of output file
      */
      void writePropagatorRGrid(const std::string & filename, int polymerID, int blockID);

      //@}

      #if 0
      // Additional functions for field-theoretic Monte-Carlo
      
      RDField<D>& pressureField();

      RDField<D>& compositionField();
      #endif

   private:

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
      * Domain object (crystallography and mesh).
      */
      Domain<D> domain_;

      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      /**
      * Homogeneous mixture, for reference.
      */
      Homogeneous::Mixture homogeneous_;

      /**
      * Pointer to Interaction (excess free energy model).
      */
      ChiInteraction* interactionPtr_;
	  
      /**
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;
     
      /**
      * Pointer to iterator factory object
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Container for wavevector data.   
      */ 
      WaveList<D>* wavelistPtr_;

      /**
      * Pointer to an Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RDField<D> > wFields_;

      /**
      * Array of chemical potential fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<WField> wFieldsRGrid_;

      /**
      * work space for chemical potential fields
      *
      */
      DArray<RDFieldDft<D> > wFieldsKGrid_;

      /**
      * Array of concentration fields for monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray< RDField<D> > cFields_;

      DArray<CField> cFieldsRGrid_;

      DArray<RDFieldDft<D> > cFieldsKGrid_;

      /**
      * Work array (size = # of grid points).
      */
      DArray<double> f_;

      /**
      * Work array (size = # of monomer types).
      */
      DArray<double> c_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Pressure times monomer volume / kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_; 

      /**
      * Has memory been allocated for fields?
      */
      bool isAllocated_;

      /**
      * Have W fields been set?
      *
      * True iff wFieldsRGrid_ has been set.
      */
      bool hasWFields_;

      /**
      * Have C fields been computed by solving the MDE ?
      *
      * True iff cFieldsRGrid_ has been set.
      */
      bool hasCFields_;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep_;

      /** 
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RDField<D> > tmpFields_;

      /** 
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
       DArray<CField> tmpFieldsRGrid_;

      /** 
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      DArray<RDFieldDft<D> > tmpFieldsKGrid_;


      IntVec<D> kMeshDimensions_;

      RDField<D> workArray;

      cudaReal* d_kernelWorkSpace_;

      cudaReal* kernelWorkSpace_;

      /**
      * Allocate memory for fields (private)
      */
      void allocate();

      /**
      * Initialize Homogeneous::Mixture object (private).
      */
      void initHomogeneous();

      /**
      * Read a filename string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      #if 0
      // Additional member variables for field-theoretic Monte Carlo
      
      RDField<D> compositionField_; //rField version

      RDFieldDft<D> compositionKField_; //kField

      RDField<D> pressureField_;

      // Free energy of the new configuration due to random change
      double fHelmholtzOld_;
      #endif

   };

   // Inline member functions

   template <int D>
   inline Domain<D> const & System<D>::domain() const
   { return domain_; }

   template <int D>
   inline UnitCell<D> & System<D>::unitCell()
   { return domain_.unitCell(); }

   // Get the Mesh<D> object.
   template <int D>
   inline Mesh<D> & System<D>::mesh()
   { return domain_.mesh(); }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> & System<D>::basis()
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> & System<D>::fft()
   {  return domain_.fft(); }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D> & System<D>::fieldIo()
   {  return domain_.fieldIo(); }

   // Get the groupName string.
   template <int D>
   inline std::string System<D>::groupName()
   { return domain_.groupName(); }

   // Get the associated Mixture<D> object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline ChiInteraction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Iterator.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the FileMaster.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the Homogeneous::Mixture object.
   template <int D>
   inline 
   Homogeneous::Mixture& System<D>::homogeneous()
   {  return homogeneous_; }

   // Get the WaveList<D> object.
   template <int D>
   inline WaveList<D>& System<D>::wavelist()
   {  return *wavelistPtr_; }

   // Get all monomer chemical potential (w field), in a basis.
   template <int D>
   inline
   DArray<RDField<D> >& System<D>::wFields()
   {  return wFields_; }

   // Get a single monomer chemical potential (w field), in a basis.
   template <int D>
   inline
   RDField<D>& System<D>::wField(int id)
   {  return wFields_[id]; }

   // Get all monomer excess chemical potential fields, on a grid.
   template <int D>
   inline 
   DArray< typename System<D>::WField >& System<D>::wFieldsRGrid()
   {  return wFieldsRGrid_; }

   // Get a single monomer hemical potential field, on a grid.
   template <int D>
   inline 
   typename System<D>::WField& System<D>::wFieldRGrid(int id)
   {  return wFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RDFieldDft<D> >& System<D>::wFieldsKGrid()
   {  return wFieldsKGrid_; }

   template <int D>
   inline
   RDFieldDft<D>& System<D>::wFieldKGrid(int id)
   {  return wFieldsKGrid_[id]; }

   // Get all monomer concentration fields, in a basis.
   template <int D>
   inline
   DArray<RDField<D> >& System<D>::cFields()
   {  return cFields_; }

   // Get a single monomer concentration field, in a basis.
   template <int D>
   inline
   RDField<D>& System<D>::cField(int id)
   {  return cFields_[id]; }

   // Get all monomer concentration fields, on a grid.
   template <int D>
   inline
   DArray< typename System<D>::CField >& System<D>::cFieldsRGrid()
   {  return cFieldsRGrid_; }

   // Get a single monomer concentration field, on a grid.
   template <int D>
   inline typename System<D>::CField& System<D>::cFieldRGrid(int id)
   {  return cFieldsRGrid_[id]; }

   template <int D>
   inline
   DArray<RDFieldDft<D> >& System<D>::cFieldsKGrid()
   {  return cFieldsKGrid_; }

   template <int D>
   inline
   RDFieldDft<D>& System<D>::cFieldKGrid(int id)
   {  return cFieldsKGrid_[id]; }

   // Get precomputed Helmoltz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {  return fHelmholtz_; }

   // Get precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {  return pressure_; }

   // Have the w fields been set?
   template <int D>
   inline bool System<D>::hasWFields() const
   {  return hasWFields_; }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   #if 0
   // Additional functions for field-theoretic Monte-Carlo
   
   template <int D>
   inline RDField<D>& System<D>::pressureField()
   { return pressureField_;}

   template <int D>
   inline RDField<D>& System<D>::compositionField()
   { return compositionField_;}
   #endif

   #ifndef PSPG_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Pspg
} // namespace Pscf
//#include "System.tpp"
#endif
