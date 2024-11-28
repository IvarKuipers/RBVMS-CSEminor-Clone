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
   }
   if (prec[1] == nullptr)
   {
      HypreParMatrix* Jpp = dynamic_cast<HypreParMatrix*>(&jacobian->GetBlock(1,1));
      HypreParMatrix *Jpp2 = const_cast<HypreParMatrix*>(Jpp);
    //  prec[1] = new HypreSmoother();//*Jpp);
      HypreBoomerAMG *amg = nullptr;
      amg = new HypreBoomerAMG(*Jpp);

      amg->SetPrintLevel(1);
      //8 i1-Jacobi smoother. Try 6 (symmetric gauss-seidel) next, afterwards 16 (ILU)
      amg->SetRelaxType(8);
      // //Aggressive coarsening
      amg->SetStrengthThresh(0.1);
      // //6 HMIS coarsening, try 8 next, may be better for non-symmetric problems
      // amg->SetCoarsening(6);
      prec[1] = amg;//*Jpp);
   }

   for (int i = 0; i < prec.Size(); ++i)
   {
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

