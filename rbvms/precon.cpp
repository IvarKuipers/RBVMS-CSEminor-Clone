// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "precon.hpp"
#include <chrono>

using namespace mfem;
using namespace RBVMS;

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{  
   auto Precon_start = std::chrono::high_resolution_clock::now();
   BlockOperator *jacobian = (BlockOperator *) &op;
   //std::cout << "Operator address for block " << 0 << ": " << &jacobian->GetBlock(0,0) << std::endl;
   //std::cout << "Operator address for block " << 1 << ": " << &jacobian->GetBlock(1,1) << std::endl;
   if (is_operator_set){
      
      for (int i = 0; i < prec.Size(); ++i)
      {
         //std::cout << "\nSetting old preconditioner as operator" << std::endl;
         SetDiagonalBlock(i, prec[i]);

         for (int j = i+1; j < prec.Size(); ++j)
         {
            SetBlock(j,i, const_cast<Operator*>(&jacobian->GetBlock(j,i)));
         }
      }

      return;
      }

   

   if (prec[0] == nullptr)
   {
      prec[0] = new HypreILU();//*Jpp);new HypreSmoother();//HypreILU()
      //std::cout << "Making a new precondtioner\n";
   }
   
   //prec[0]->SetOperator(jacobian->GetBlock(0,0));
   if (prec[1] == nullptr)
   {
      HypreParMatrix* Jpp = dynamic_cast<HypreParMatrix*>(&jacobian->GetBlock(1,1));
      HypreParMatrix *Jpp2 = const_cast<HypreParMatrix*>(Jpp);
    //  prec[1] = new HypreSmoother();//*Jpp);
      HypreILU *ilu = new HypreILU();
      //Parasails->SetType(HypreSmoother::FIR);
      //Parasails->SetParams(0.05, 1);
      //Parasails->SetFilter(0.01); 
      //Parasails->SetReuse(1);        
      //Parasails->SetLogging(1); 
      prec[1] = ilu;//*Jpp);
      //std::cout << "Making a new precondtioner\n";
   }
  // std::cout << "Setting Operator for block 1,1" << std::endl;
   
   for (int i = 0; i < prec.Size(); ++i)
   {
      //std::cout << "\nSetting new preconditioner as operator" << std::endl;
    
      prec[i]->SetOperator(jacobian->GetBlock(i,i));
      SetDiagonalBlock(i, prec[i]);

      for (int j = i+1; j < prec.Size(); ++j)
      {
         SetBlock(j,i, const_cast<Operator*>(&jacobian->GetBlock(j,i)));
      }
   }
   auto Precon_end = std::chrono::high_resolution_clock::now();
   auto Precon_duration = std::chrono::duration_cast<std::chrono::microseconds>(Precon_end - Precon_start).count();
   std::cout << std::endl <<"Setup time for new preconditioner:  " << Precon_duration/1000000.0 <<  std::endl;
   //Has to be set to true to work
   is_operator_set = true;
}

void JacobianPreconditioner::ResetOperatorSetup()
{
   is_operator_set = false;
}

// Destructor
JacobianPreconditioner::~JacobianPreconditioner()
{
   for (int i = 0; i < prec.Size(); ++i)
   {
      delete prec[i];
   }
}

