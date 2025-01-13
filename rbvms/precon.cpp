// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "precon.hpp"

using namespace mfem;
using namespace RBVMS;

// Set the diagonal and off-diagonal operators
void JacobianPreconditioner::SetOperator(const Operator &op)
{
   double t0 = MPI_Wtime();
   BlockOperator *jacobian = (BlockOperator *) &op;

   if (prec[0] == nullptr)
   {
      prec[0] = new HypreILU();
   }

   std::cout << "Setting Operator for block 0,0" << std::endl;
   prec[0]->SetOperator(jacobian->GetBlock(0,0));
   if (prec[1] == nullptr)
   {
      HypreParMatrix* Jpp = dynamic_cast<HypreParMatrix*>(&jacobian->GetBlock(1,1));
      HypreParMatrix *Jpp2 = const_cast<HypreParMatrix*>(Jpp);

	  // Hier de Parasails teringzooi voor later, jahhhhh
	  
      //HypreParaSails *Parasails = new HypreParaSails(*Jpp);
      //Parasails->SetFilter(0.01);
      //Parasails-> SetSymmetry(2);   

      // HypreILU parameter tuning testing. 
      HypreILU *ilu = new HypreILU();  

      // Tune ILU parameters here
      //ilu->SetLevelOfFill(1);              // ILU(1)
      //ilu->SetDropThreshold(1e-4);         // Drop tolerance
      ilu->SetType(1);                     // RCM reordering
      //ilu->SetMaxIter(2);                  // Iterative application

      prec[1] = ilu;
      
   }
   //std::cout << "Setting Operator for block 1,1" << std::endl;

   for (int i = 0; i < prec.Size(); ++i)
   {
      prec[i]->SetOperator(jacobian->GetBlock(i,i));
      SetDiagonalBlock(i, prec[i]);

      for (int j = i+1; j < prec.Size(); ++j)
      {
         SetBlock(j,i, const_cast<Operator*>(&jacobian->GetBlock(j,i)));
      }
   }
   double t1 = MPI_Wtime(); 
   double setup_time = t1 - t0;
   std::cout << "Setup time for new preconditioner:  " << setup_time << std::endl;
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
