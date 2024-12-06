// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "precon.hpp"

using namespace mfem;
using namespace RBVMS;

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{
   BlockOperator *jacobian = (BlockOperator *) &op;

   if (prec[0] == nullptr)
   {
      prec[0] = new HypreILU();//*Jpp);new HypreSmoother();//HypreILU()
      std::cout << "Making a new precondtioner\n";
   }
   
   //prec[0]->SetOperator(jacobian->GetBlock(0,0));
   if (prec[1] == nullptr)
   {
      HypreParMatrix* Jpp = dynamic_cast<HypreParMatrix*>(&jacobian->GetBlock(1,1));
      HypreParMatrix *Jpp2 = const_cast<HypreParMatrix*>(Jpp);
    //  prec[1] = new HypreSmoother();//*Jpp);
      HypreSmoother *Parasails = new HypreSmoother(*Jpp);
      Parasails->SetType(HypreSmoother::l1Jacobi);
      //Parasails->SetParams(0.05, 1);
      //Parasails->SetFilter(0.01); 
      //Parasails->SetReuse(1);        
      //Parasails->SetLogging(1); 
      prec[1] = Parasails;//*Jpp);
      std::cout << "Making a new precondtioner\n";
   }
  // std::cout << "Setting Operator for block 1,1" << std::endl;
   
   for (int i = 0; i < prec.Size(); ++i)
   {
      //std::cout << "\nSetting preconditioner as operator" << std::endl;
      prec[i]->SetOperator(jacobian->GetBlock(i,i));
      SetDiagonalBlock(i, prec[i]);

      for (int j = i+1; j < prec.Size(); ++j)
      {
         SetBlock(j,i, const_cast<Operator*>(&jacobian->GetBlock(j,i)));
      }
   }
}

// Destructor
JacobianPreconditioner::~JacobianPreconditioner()
{
   for (int i = 0; i < prec.Size(); ++i)
   {
      delete prec[i];
   }
}

