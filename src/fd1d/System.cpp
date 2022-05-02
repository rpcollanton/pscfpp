/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"

#include <fd1d/iterator/Iterator.h>
#include <fd1d/iterator/IteratorFactory.h>
#include <fd1d/sweep/Sweep.h>
#include <fd1d/sweep/SweepFactory.h>
#include <fd1d/iterator/NrIterator.h>
#include <fd1d/misc/HomogeneousComparison.h>
#include <fd1d/misc/FieldIo.h>

#include <pscf/inter/Interaction.h>
#include <pscf/inter/ChiInteraction.h>
#include <pscf/homogeneous/Clump.h>

#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/param/BracketPolicy.h>

#include <string>
#include <unistd.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   /*
   * Constructor.
   */
   System::System()
    : mixture_(),
      domain_(),
      fileMaster_(),
      homogeneous_(),
      interactionPtr_(0),
      iteratorPtr_(0),
      iteratorFactoryPtr_(0),
      sweepPtr_(0),
      sweepFactoryPtr_(0),
      wFields_(),
      cFields_(),
      f_(),
      c_(),
      fHelmholtz_(0.0),
      pressure_(0.0),
      hasMixture_(0),
      hasDomain_(0),
      hasFields_(0),
      hasSweep_(0)
   {  
      setClassName("System"); 

      interactionPtr_ = new ChiInteraction(); 
      iteratorFactoryPtr_ = new IteratorFactory(*this); 
      sweepFactoryPtr_ = new SweepFactory(*this);

      BracketPolicy::set(BracketPolicy::Optional);
   }

   /*
   * Destructor.
   */
   System::~System()
   {
      if (interactionPtr_) {
         delete interactionPtr_;
      }
      if (iteratorPtr_) {
         delete iteratorPtr_;
      }
      if (iteratorFactoryPtr_) {
         delete iteratorFactoryPtr_;
      }
      if (sweepPtr_) {
         delete sweepPtr_;
      }
      if (sweepFactoryPtr_) {
         delete sweepFactoryPtr_;
      }
   }

   /*
   * Process command line options.
   */
   void System::setOptions(int argc, char **argv)
   {
      bool eflag = false;  // echo
      bool pFlag = false;  // param file 
      bool cFlag = false;  // command file 
      bool iFlag = false;  // input prefix
      bool oFlag = false;  // output prefix
      char* pArg = 0;
      char* cArg = 0;
      char* iArg = 0;
      char* oArg = 0;
   
      // Read program arguments
      int c;
      opterr = 0;
      while ((c = getopt(argc, argv, "er:p:c:i:o:f")) != -1) {
         switch (c) {
         case 'e':
            eflag = true;
            break;
         case 'p': // parameter file
            pFlag = true;
            pArg  = optarg;
            break;
         case 'c': // command file
            cFlag = true;
            cArg  = optarg;
            break;
         case 'i': // input prefix
            iFlag = true;
            iArg  = optarg;
            break;
         case 'o': // output prefix
            iFlag = true;
            oArg  = optarg;
            break;
         case '?':
           Log::file() << "Unknown option -" << optopt << std::endl;
           UTIL_THROW("Invalid command line option");
         }
      }
   
      // Set flag to echo parameters as they are read.
      if (eflag) {
         Util::ParamComponent::setEcho(true);
      }

      // If option -p, set parameter file name
      if (pFlag) {
         fileMaster().setParamFileName(std::string(pArg));
      }

      // If option -c, set command file name
      if (cFlag) {
         fileMaster().setCommandFileName(std::string(cArg));
      }

      // If option -i, set path prefix for input files
      if (iFlag) {
         fileMaster().setInputPrefix(std::string(iArg));
      }

      // If option -o, set path prefix for output files
      if (oFlag) {
         fileMaster().setOutputPrefix(std::string(oArg));
      }

   }

   /*
   * Read parameters and initialize.
   */
   void System::readParameters(std::istream& in)
   {
      readParamComposite(in, mixture());
      hasMixture_ = true;

      // Initialize homogeneous object
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      int ns = mixture().nSolvent(); 
      homogeneous_.setNMolecule(np+ns);
      homogeneous_.setNMonomer(nm);
      initHomogeneous();

      interaction().setNMonomer(mixture().nMonomer());
      readParamComposite(in, interaction());

      readParamComposite(in, domain());
      hasDomain_ = true;
      allocateFields();

      // Initialize iterator
      //readParamComposite(in, iterator());

      // Instantiate and initialize an Iterator 
      std::string className;
      bool isEnd;
      iteratorPtr_ = iteratorFactoryPtr_->readObject(in, *this, 
                                                     className, isEnd);
      if (!iteratorPtr_) {
         std::string msg = "Unrecognized Iterator subclass name ";
         msg += className;
         UTIL_THROW(msg.c_str());
      }

      // Optionally instantiate a Sweep object
      readOptional<bool>(in, "hasSweep", hasSweep_);
      if (hasSweep_) {
         std::string className;
         bool isEnd;
         sweepPtr_ = 
            sweepFactoryPtr_->readObject(in, *this, className, isEnd);
         if (!sweepPtr_) {
            UTIL_THROW("Unrecognized Sweep subclass name");
         }
         sweepPtr_->setSystem(*this);
      }
   }

   /*
   * Read default parameter file.
   */
   void System::readParam(std::istream& in)
   {
      readBegin(in, className().c_str());  
      readParameters(in);  
      readEnd(in);  
   }

   /*
   * Read default parameter file.
   */
   void System::readParam()
   {  readParam(fileMaster().paramFile()); }

   /*
   * Read parameters and initialize.
   */
   void System::allocateFields()
   {
      // Preconditions
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasDomain_);

      // Allocate memory in mixture
      mixture().setDomain(domain());

      // Allocate wFields and cFields
      int nMonomer = mixture().nMonomer();
      wFields_.allocate(nMonomer);
      cFields_.allocate(nMonomer);
      int nx = domain().nx();
      for (int i = 0; i < nMonomer; ++i) {
         wField(i).allocate(nx);
         cField(i).allocate(nx);
      }
      hasFields_ = true;
   }

   /*
   * Read and execute commands from a specified command file.
   */
   void System::readCommands(std::istream &in)
   {
      UTIL_CHECK(hasMixture_);
      UTIL_CHECK(hasDomain_);

      std::string command;
      std::string filename;
      FieldIo fieldIo(*this);

      std::istream& inBuffer = in;

      bool readNext = true;
      while (readNext) {

         inBuffer >> command;
         Log::file() << command;

         if (command == "FINISH") {
            Log::file() << std::endl;
            readNext = false;
         } else
         if (command == "READ_W") {
            readEcho(inBuffer, filename);
            fieldIo.readFields(wFields(), filename);  
         } else
         if (command == "ITERATE") {
            Log::file() << std::endl;
            iterator().solve();
            outputThermo(Log::file());
         } else 
         if (command == "COMPARE_HOMOGENEOUS") {
            int mode;
            inBuffer >> mode;
            Log::file() << std::endl;
            Log::file() << "mode       = " << mode << std::endl;

            HomogeneousComparison comparison(*this);
            comparison.compute(mode);
            comparison.output(mode, Log::file());
         } else 
         if (command == "SWEEP") {
            UTIL_CHECK(hasSweep_);
            UTIL_CHECK(sweepPtr_);
            sweepPtr_->sweep();
         } else 
         if (command == "WRITE_W") {
            readEcho(inBuffer, filename);
            fieldIo.writeFields(wFields(), filename);  
         } else
         if (command == "WRITE_C") {
            readEcho(inBuffer, filename);
            fieldIo.writeFields(cFields(), filename);  
         } else
         if (command == "WRITE_BLOCK_C") {
            readEcho(inBuffer, filename);
            fieldIo.writeBlockCFields(filename);  
         } else
         if (command == "WRITE_VERTEX_Q") {
            readEcho(inBuffer, filename);
            int polymerId, vertexId;
            inBuffer >> polymerId;
            Log::file() << "polymerId = " 
                        << Int(polymerId, 5) << std::endl;
            inBuffer >> vertexId;
            Log::file() << "vertexId  = " 
                        << Int(vertexId, 5) << std::endl;
            fieldIo.writeVertexQ(polymerId, vertexId, filename);  
         } else
         if (command == "REMESH_W") {
            int nx;
            inBuffer >> nx;
            Log::file() << std::endl;
            Log::file() << "nx      = " << Int(nx, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo.remesh(wFields(), nx, filename);
         } else
         if (command == "EXTEND_W") {
            int m;
            inBuffer >> m;
            Log::file() << std::endl;
            Log::file() << "m       = " << Int(m, 20) << std::endl;
            inBuffer >> filename;
            Log::file() << "outfile = " << Str(filename, 20) << std::endl;
            fieldIo.extend(wFields(), m, filename);
         } else {
            Log::file() << "  Error: Unknown command  " << command << std::endl;
            readNext = false;
         }

      }
   }

   /*
   * Read and execute commands from the default command file.
   */
   void System::readCommands()
   {  
      if (fileMaster().commandFileName().empty()) {
         UTIL_THROW("Empty command file name");
      }
      readCommands(fileMaster().commandFile()); 
   }

   /*
   * Read a string (e.g., a filename) and echo to the log file.
   */
   void System::readEcho(std::istream& in, std::string& string) const
   {
      in >> string;
      Log::file() << "  " << Str(string, 20) << std::endl;
   }

   /*
   * Compute Helmoltz free energy and pressure
   */
   void System::computeFreeEnergy()
   {
      fHelmholtz_ = 0.0;
      int np = mixture().nPolymer();
      int ns = mixture().nSolvent();
      int nm = mixture().nMonomer();
      double phi, mu;
 
      // Polymeric ideal gas contributions to fHelhmoltz_
      if (np > 0) {
         Polymer* polymerPtr;
         double length;
         for (int i = 0; i < np; ++i) {
            polymerPtr = &mixture().polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            // Recall: mu = ln(phi/q)
            length = polymerPtr->length();
            fHelmholtz_ += phi*( mu - 1.0 )/length;
         }
      }

      // Solvent ideal gas contributions to fHelhmoltz_
      if (ns > 0) {
         Solvent* solventPtr;
         double size;
         for (int i = 0; i < ns; ++i) {
            solventPtr = &mixture().solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            // Recall: mu = ln(phi/q)
            size = solventPtr->size();
            fHelmholtz_ += phi*( mu - 1.0 )/size;
         }
      }

      // Apply Legendre transform subtraction
      for (int i = 0; i < nm; ++i) {
         fHelmholtz_ -= 
                  domain().innerProduct(wFields_[i], cFields_[i]);
      }

      // Add average interaction free energy density per monomer
      int nx = domain().nx();
      if (!f_.isAllocated()) f_.allocate(nx);
      if (!c_.isAllocated()) c_.allocate(nm);
      int j;
      for (int i = 0; i < nx; ++i) { 
         // Store c_[j] = local concentration of species j
         for (j = 0; j < nm; ++j) {
            c_[j] = cFields_[j][i];
         }
         // Compute f_[i] = excess free eenrgy at grid point i
         f_[i] = interaction().fHelmholtz(c_);
      }
      fHelmholtz_ += domain().spatialAverage(f_);

      // Compute pressure
      pressure_ = -fHelmholtz_;
      if (np > 0) {
         double length;
         Polymer* polymerPtr;
         for (int i = 0; i < np; ++i) {
            polymerPtr = & mixture().polymer(i);
            phi = polymerPtr->phi();
            mu = polymerPtr->mu();
            length = polymerPtr->length();
            pressure_ += phi*mu/length;
         }
      }
      if (ns > 0) {
         double size;
         Solvent* solventPtr;
         for (int i = 0; i < ns; ++i) {
            solventPtr = & mixture().solvent(i);
            phi = solventPtr->phi();
            mu = solventPtr->mu();
            size = solventPtr->size();
            pressure_ += phi*mu/size;
         }
      }

   }

   void System::initHomogeneous()
   {

      // Set number of molecular species and monomers
      int nm = mixture().nMonomer(); 
      int np = mixture().nPolymer(); 
      int ns = mixture().nSolvent(); 
      UTIL_CHECK(homogeneous_.nMolecule() == np + ns);
      UTIL_CHECK(homogeneous_.nMonomer() == nm);

      // Allocate c_ work array, if necessary
      if (c_.isAllocated()) {
         UTIL_CHECK(c_.capacity() == nm);
      } else {
         c_.allocate(nm);
      }

      int i;   // molecule index
      int j;   // monomer index
      int k;   // block or clump index
      int nb;  // number of blocks
      int nc;  // number of clumps
 
      // Loop over polymer molecule species
      if (np > 0) {
         for (i = 0; i < np; ++i) {
   
            // Initial array of clump sizes 
            for (j = 0; j < nm; ++j) {
               c_[j] = 0.0;
            }
   
            // Compute clump sizes for all monomer types.
            nb = mixture().polymer(i).nBlock(); 
            for (k = 0; k < nb; ++k) {
               Block& block = mixture().polymer(i).block(k);
               j = block.monomerId();
               c_[j] += block.length();
            }
    
            // Count the number of clumps of nonzero size
            nc = 0;
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  ++nc;
               }
            }
            homogeneous_.molecule(i).setNClump(nc);
    
            // Set clump properties for this Homogeneous::Molecule
            k = 0; // Clump index
            for (j = 0; j < nm; ++j) {
               if (c_[j] > 1.0E-8) {
                  homogeneous_.molecule(i).clump(k).setMonomerId(j);
                  homogeneous_.molecule(i).clump(k).setSize(c_[j]);
                  ++k;
               }
            }
            homogeneous_.molecule(i).computeSize();
   
         }
      }

      // Add solvent contributions
      if (np > 0) {
         double size;
         int monomerId;
         for (int is = 0; is < ns; ++is) {
            i = is + np;
            monomerId = mixture().solvent(is).monomerId();
            size = mixture().solvent(is).size();
            homogeneous_.molecule(i).setNClump(1);
            homogeneous_.molecule(i).clump(0).setMonomerId(monomerId);
            homogeneous_.molecule(i).clump(0).setSize(size);
            homogeneous_.molecule(i).computeSize();
         }
      }

   }

   void System::outputThermo(std::ostream& out)
   {
      out << std::endl;
      out << "fHelmholtz = " << Dbl(fHelmholtz(), 18, 11) << std::endl;
      out << "pressure   = " << Dbl(pressure(), 18, 11) << std::endl;
      out << std::endl;

      // Polymers
      int np = mixture().nPolymer();
      if (np > 0) {
         out << "Polymers:" << std::endl;
         out << "    i"
             << "        phi[i]      "
             << "        mu[i]       " 
             << std::endl;
         for (int i = 0; i < np; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture().polymer(i).phi(),18, 11)
                << "  " << Dbl(mixture().polymer(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

      // Solvents
      int ns = mixture().nSolvent();
      if (ns > 0) {
         out << "Solvents:" << std::endl;
         out << "    i"
             << "        phi[i]      "
             << "        mu[i]       " 
             << std::endl;
         for (int i = 0; i < ns; ++i) {
            out << Int(i, 5) 
                << "  " << Dbl(mixture().solvent(i).phi(),18, 11)
                << "  " << Dbl(mixture().solvent(i).mu(), 18, 11)  
                << std::endl;
         }
         out << std::endl;
      }

   }

} // namespace Fd1d
} // namespace Pscf
