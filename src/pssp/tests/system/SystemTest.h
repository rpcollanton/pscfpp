#ifndef PSSP_SYSTEM_TEST_H
#define PSSP_SYSTEM_TEST_H

#include <test/UnitTest.h>
#include <test/UnitTestRunner.h>

#include <pssp/System.h>
#include <pssp/iterator/AmIterator.h>
#include <pscf/mesh/MeshIterator.h>

#include <fstream>

using namespace Util;
using namespace Pscf;
using namespace Pscf::Pssp;

class SystemTest : public UnitTest
{

public:

   void setUp()
   {}

   void tearDown()
   {}

   void testConstructor1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;
   }

   void testReadParameters1D()
   {
      printMethod(TEST_FUNC);
      System<1> system;

      std::ifstream in;
      openInputFile("in/domainOn/System1D", in);
      system.readParam(in);
      in.close();
   }
 
   void testConversion1D_lam() 
   {   
      printMethod(TEST_FUNC);
 
      System<1> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System1D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_1d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_1d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      // Compare fields
      double err;
      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
            err = wFields_check[i][j] - system.wFields()[i][j];
            if ((std::abs(err) > 1.0E-8)){
               diff = false;
               break;
            } else {
               diff = true;
            }
         }   
         if (diff == false) break;
      }   
      TEST_ASSERT(diff);
   }  

   void testConversion2D_hex() 
   {   
      printMethod(TEST_FUNC);
 
      System<2> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System2D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_2d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_2d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      double err;
      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           err = wFields_check[i][j] - system.wFields()[i][j];
           if((std::abs(err) > 1.0E-8)){
               diff = false;
               break;
            }   
            else
               diff = true;
         }   
         if (diff==false) break;
      }   
      TEST_ASSERT(diff);
   }  

   void testConversion3D_bcc() 
   {   
      printMethod(TEST_FUNC);
 
      System<3> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System3D", in);
          
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/conv/Conversion_3d_step1", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){         
            wFields_check[i][j] = system.wFields() [i] [j];
         }
      }

      std::ifstream command_2;
      openInputFile("in/conv/Conversion_3d_step2", command_2);
      system.readCommands(command_2);
      command_2.close();

      // Compare fields
      double err;
      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           err = wFields_check[i][j] - system.wFields()[i][j];
           if((std::abs(err) > 1.0E-8)){
               diff = false;
               break;
            }
            else
               diff = true;
         } 
         if(diff==false) break;
      }
      TEST_ASSERT(diff);
   }   

   void testIterate1D_lam_rigid()
   {
      printMethod(TEST_FUNC);

      System<1> system;
      std::ifstream in;
      openInputFile("in/domainOff/System1D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_lam", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate1d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 5.07058e-08)){
               // The above is the minimum error in the omega field.
               // Occurs for the first star                 
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }    
         if(diff==false){
            break;
         }    
      }    
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.006583929) < 1.0E-8){
         //0.006583929 is the stress calculated 
         //for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }

   void testIterate1D_lam_flex()
   {
      printMethod(TEST_FUNC);
 
      System<1> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System1D", in);
    
      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_lam", command);
      system.readCommands(command);
      command.close();
 
      int nMonomer = system.mixture().nMonomer();
      DArray<RField<1> > wFields_check;
      DArray<RField<1> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){    
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate1d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) > 1.0E-8)){
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }   
            else
               diff = true;
         }   
         if(diff==false){
            break;
         }
      }   
      TEST_ASSERT(diff);
   }

   void testIterate2D_hex_rigid()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      std::ifstream in;
      openInputFile("in/domainOff/System2D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_hex", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate2d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.60828e-07)){
               // The above is the minimum error in the omega field.
               // Occurs for the first star            
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }    
         if(diff==false){
            break;
         }    
      }    
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.010633960) < 1.0E-8){
         //0.010633960 is the stress calculated 
         //for this omega field for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }

   void testIterate2D_hex_flex()
   {
      printMethod(TEST_FUNC);

      System<2> system;
      std::ifstream in;
      openInputFile("in/domainOn/System2D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_hex", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<2> > wFields_check;
      DArray<RField<2> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j];
         }
      }

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate2d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 2.58007e-07)){
               // The above is the maximum error in the omega field.
               // Occurs for the first star
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }
            else
               diff = true;
         }
         if(diff==false){

            break;
         }
      }
      TEST_ASSERT(diff);
   }

   void testIterate3D_bcc_rigid()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in;
      openInputFile("in/domainOff/System3D", in); 

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOff/ReadOmega_bcc", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }    

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j]; 
         }    
      }    

      std::ifstream command_2;
      openInputFile("in/domainOff/Iterate3d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >= 1.02291e-07)){
               // The above is the maximum error in the omega field.
               // Occurs for the second star.               
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }    
            else 
               diff = true;
         }
         if(diff==false){
            break;
         }
      }
      bool stress = false;
      if(std::abs(system.mixture().TStress[0] - 0.005242863) < 1.0E-8){
         //0.005242863 is the stress calculated for this omega field 
         //for no stress relaxation using Fortran
         stress = true;
      }

      TEST_ASSERT(stress);
      TEST_ASSERT(diff);
   }

   void testIterate3D_bcc_flex()
   {
      printMethod(TEST_FUNC);

      System<3> system;
      std::ifstream in; 
      openInputFile("in/domainOn/System3D", in);

      system.readParam(in);
      in.close();
      std::ifstream command;
      openInputFile("in/domainOn/ReadOmega_bcc", command);
      system.readCommands(command);
      command.close();

      int nMonomer = system.mixture().nMonomer();
      DArray<RField<3> > wFields_check;
      DArray<RField<3> > wFields;
      wFields_check.allocate(nMonomer);
      wFields.allocate(nMonomer);
      int ns = system.basis().nStar();
      for (int i = 0; i < nMonomer; ++i){
          wFields_check[i].allocate(ns);
      }   

      for (int i = 0; i < nMonomer; ++i){
         for (int j = 0; j < ns; ++j){
            wFields_check[i][j] = system.wFields() [i] [j];
         }   
      }   

      std::ifstream command_2;
      openInputFile("in/domainOn/Iterate3d", command_2);
      system.readCommands(command_2);
      command_2.close();

      bool diff = true;
      for (int j = 0; j < ns; ++j) {
         for (int i = 0; i < nMonomer; ++i) {
           if((std::abs(wFields_check[i][j] - system.wFields()[i][j]) >=  1.09288e-07)){ 
               // The above is the maximum error in the omega field.
               // Occurs for the second star.
               diff = false;
               std::cout <<"This is error for break:"<< 
                  (std::abs(wFields_check[i][j] - system.wFields()[i][j])) <<std::endl;
               std::cout <<"ns = "<< j << std::endl;
               break;
            }
            else
               diff = true;
         }
         if(diff==false){

            break;
         }
      }
      TEST_ASSERT(diff);

   }


};

TEST_BEGIN(SystemTest)
TEST_ADD(SystemTest, testConstructor1D)
TEST_ADD(SystemTest, testReadParameters1D)
TEST_ADD(SystemTest, testConversion1D_lam)
TEST_ADD(SystemTest, testConversion2D_hex)
TEST_ADD(SystemTest, testConversion3D_bcc)
TEST_ADD(SystemTest, testIterate1D_lam_rigid)
TEST_ADD(SystemTest, testIterate1D_lam_flex)
TEST_ADD(SystemTest, testIterate2D_hex_rigid)
TEST_ADD(SystemTest, testIterate2D_hex_flex)
TEST_ADD(SystemTest, testIterate3D_bcc_rigid)
TEST_ADD(SystemTest, testIterate3D_bcc_flex)

TEST_END(SystemTest)

#endif
