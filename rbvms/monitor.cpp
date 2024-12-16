// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "monitor.hpp"
#include <chrono>

using namespace mfem;
using namespace RBVMS;

// Constructor
GeneralResidualMonitor::GeneralResidualMonitor(MPI_Comm comm,
                                               const std::string& prefix_,
                                               int print_lvl)
   : prefix(prefix_)
{
   if (Mpi::Root())
   {
      print_level = print_lvl;
   }
   else
   {
      print_level = -1;
   }
}

// Print residual
void GeneralResidualMonitor::MonitorResidual(int it,
                                             real_t norm,
                                             const Vector &r,
                                             bool final)
{
   static auto Krylov_start = std::chrono::high_resolution_clock::now();
   if (it == 0)
   {
      norm0 = norm;
      Krylov_start = std::chrono::high_resolution_clock::now();
   }

   if ((print_level > 0) && (( it%print_level == 0) || final))
   {
      auto Krylov_end = std::chrono::high_resolution_clock::now();
      auto Krylov_duration = std::chrono::duration_cast<std::chrono::milliseconds>(Krylov_end - Krylov_start).count();
      mfem::out<<prefix<<" iteration "<<std::setw(3)<<it
               <<": ||r|| = "
               <<std::setw(8)<<std::defaultfloat<<std::setprecision(3)
               <<norm
               <<", ||r||/||r_0|| = "
               <<std::setw(6)<<std::fixed<<std::setprecision(2)
               <<100*norm/norm0<<  " %" << ", Time spent: " << Krylov_duration/1000.0 << " seconds"<< "\n";
      Krylov_start = Krylov_end;
   }
}

// Constructor
SystemResidualMonitor::SystemResidualMonitor(MPI_Comm comm,
                                             const std::string& prefix_,
                                             int print_lvl,
                                             Array<int> &offsets)
   : prefix(prefix_), bOffsets(offsets), iteration_count(0)
{
   if (Mpi::Root())
   {
      print_level = print_lvl;
   }
   else
   {
      print_level = -1;
   }
   nvar = bOffsets.Size()-1;
   norm0.SetSize(nvar);
}

void SystemResidualMonitor::ComputeResiduals(const Vector &r, Vector &vnorm)
{
   for (int i = 0; i < nvar; ++i)
   {
      Vector r_i(r.GetData() + bOffsets[i], bOffsets[i + 1] - bOffsets[i]);
      vnorm[i] = sqrt(InnerProduct(MPI_COMM_WORLD, r_i, r_i));
   }
}

int SystemResidualMonitor::GetIterationCount()
{return iteration_count;}

void SystemResidualMonitor::ResetCounter()
{iteration_count = 0;}
// Print residual
void SystemResidualMonitor::MonitorResidual(int it,
                                            real_t norm,
                                            const Vector &r,
                                            bool final)
{
   static auto Newt_start = std::chrono::high_resolution_clock::now();
   static Vector vnorm(nvar);
   ComputeResiduals(r, vnorm);
   for (int i = 0; i < nvar; ++i)
   {
      if (it == 0) { norm0[i] = vnorm[i];  Newt_start = std::chrono::high_resolution_clock::now();}
   }

   bool print = (print_level > 0) &&  ((it%print_level == 0) || final);
   if (print)
   {
      auto Newt_end = std::chrono::high_resolution_clock::now();
      auto Newt_duration = std::chrono::duration_cast<std::chrono::milliseconds>(Newt_end - Newt_start).count();
      mfem::out << prefix << " iteration " << std::setw(3) << it <<"\n"
                << " ||r||  \t"<< "||r||/||r_0||" << ", Time spent: " << Newt_duration/1000.0 << " seconds"<< std::endl;
      for (int i = 0; i < nvar; ++i)
      {
         mfem::out<<std::setw(8)<<std::defaultfloat<<std::setprecision(4)
                  <<vnorm[i]<<"\t"
                  <<std::setw(8)<<std::fixed<<std::setprecision(2)
                  <<100*vnorm[i]/norm0[i] << " %\n";
      }
      Newt_start = Newt_end;
   }
   iteration_count++;
   mfem::out<<std::flush;
}
