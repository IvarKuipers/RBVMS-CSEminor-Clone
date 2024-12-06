// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_PRECON_HPP
#define RBVMS_PRECON_HPP

#include "mfem.hpp"
#include "tau.hpp"

using namespace mfem;

namespace RBVMS
{

// Custom block preconditioner for the Jacobian
class JacobianPreconditioner : public
   BlockLowerTriangularPreconditioner
{
protected:
   Array<Solver *> prec;
   bool is_operator_set;

public:
   /// Constructor
   JacobianPreconditioner(Array<int> &offsets)
      : BlockLowerTriangularPreconditioner (offsets), prec(offsets.Size()-1), is_operator_set(false)
   { prec = nullptr;};

   //Reset the Boolean
   virtual void ResetOperatorSetup();
   /// Set the diagonal and off-diagonal operators
   virtual void SetOperator(const Operator &op);

   // Destructor
   virtual ~JacobianPreconditioner();
};

} // namespace RBVMS

#endif
