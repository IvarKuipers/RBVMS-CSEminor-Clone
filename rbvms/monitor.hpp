// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------


#ifndef RBVMS_MONITOR_HPP
#define RBVMS_MONITOR_HPP

#include "mfem.hpp"
#include "tau.hpp"

using namespace mfem;

namespace RBVMS
{

/// This class help monitor the convergence of the linear Krylov solve.
class GeneralResidualMonitor : public IterativeSolverMonitor
{
private:
   const std::string prefix;
   int print_level;
   mutable real_t norm0;

public:
   /// Constructor
   GeneralResidualMonitor(MPI_Comm comm,
                          const std::string& prefix_,
                          int print_lvl);

   /// Print residual
   virtual void MonitorResidual(int it,
                                real_t norm,
                                const Vector &r,
                                bool final);
};

/// This class help monitor the convergence of the nonlinear Newton solve.
class SystemResidualMonitor : public IterativeSolverMonitor
{
private:
   const std::string prefix;
   int print_level, nvar;
   mutable Vector norm0;
   Array<int> &bOffsets;

public:
   /// Constructor
   SystemResidualMonitor(MPI_Comm comm,
                         const std::string& prefix_,
                         int print_lvl,
                         Array<int> &offsets);

   /// Calculate residual
   virtual void ComputeResiduals(const Vector &r, Vector &vnorm);

   /// Print residual
   virtual void MonitorResidual(int it,
                                real_t norm,
                                const Vector &r,
                                bool final);
};

} // namespace RBVMS

#endif
