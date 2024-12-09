// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#ifndef RBVMS_EVOLUTION_HPP
#define RBVMS_EVOLUTION_HPP

#include "mfem.hpp"
#include "weakform.hpp"

using namespace std;
using namespace mfem;

namespace RBVMS
{

// Predefine class
class ParTimeDepBlockNonlinForm;

/** This class provide the correct interface between the time-dependent 
    block nonlinear form (defined below) and the MFEM::ODESolver.
*/
class Evolution : public TimeDependentOperator
{
private:
   ParTimeDepBlockNonlinForm &form;
   Solver &solver;
   Vector dudt;
   int count = 0;

public:
   /// Constructor
   Evolution(ParTimeDepBlockNonlinForm &form,
             Solver &solver);

   /// Stub for explicit solve of time dependent problem
   virtual void Mult(const Vector &x, Vector &k) const override { k = 0.0;};

   /// Solve time dependent problem
   virtual void ImplicitSolve(const real_t dt,
                              const Vector &x,
                              Vector &k) override;

   /// Get the CFL Number
   real_t GetCFL() const;

   /// Get the Force on each of the boundaries
   DenseMatrix GetForce();

   /// Destructor
   ~Evolution() {}
};


/** This class is a specialized ParBlockNonlinearForm that includes timestepping
    interpolation, vis:
       add(x0,dt,dx,x);   // x = x0 + dt*dx

    Both x and dx need to be passed to the FormIntegrator.
*/
class ParTimeDepBlockNonlinForm : public ParBlockNonlinearForm
{
private:
   RBVMS::IncNavStoIntegrator &integrator;

   /// Numerical parameters
   real_t dt;
   mutable real_t cfl;

   Array<int> strongBCBdr;
   Array<int> weakBCBdr;
   Array<int> outflowBdr;

   /// Solution & Residual vector
   mutable Vector x0, x;
   mutable BlockVector dxs;
   mutable BlockVector dxs_true;

   /// Conservative boundary forces
   mutable DenseMatrix bdrForce;
   mutable bool hasGrad;

public:
   /// Constructor
   ParTimeDepBlockNonlinForm(Array<ParFiniteElementSpace *> &pfes,
                             RBVMS::IncNavStoIntegrator &integrator);

   void SetStrongBC (Array<int> strong_bdr);
   void SetWeakBC   (Array<int> weak_bdr);
   void SetOutflowBC(Array<int> outflow_bdr);

   /// Set the solution of the previous time step @a x0
   /// and the timestep size @a dt of the current solve.
   void SetTimeAndSolution(const real_t t,
                           const real_t dt,
                           const Vector &x0);

   void ResetGradient();

   /// Get the CFL-Number
   real_t GetCFL() { return cfl;};

   /// Get the conservative boundary forces
   DenseMatrix& GetForce() { return bdrForce;};

   /// Specialized version of GetEnergy() for BlockVectors
   //real_t GetEnergyBlocked(const BlockVector &bx) const;

   /// Block T-Vector to Block T-Vector
   void Mult(const Vector &x, Vector &y) const;

   /// Specialized version of Mult() for BlockVectors
   /// Block L-Vector to Block L-Vector
   void MultBlocked(const BlockVector &bx,
                    const BlockVector &dbx,
                    BlockVector &by) const;

   virtual BlockOperator &GetGradient(const Vector &x) const;

   /// Return the local block gradient matrix for the given true-dof vector x
   const BlockOperator& GetLocalGradient(const Vector &x) const;

   /// Specialized version of GetGradient() for BlockVector
   void ComputeGradientBlocked(const BlockVector &bx,
                               const BlockVector &dbx) const;
};

} // namespace RBVMS

#endif

