// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "weakform.hpp"
#include <chrono>

using namespace mfem;
using namespace RBVMS;

// Constructor
IncNavStoIntegrator::IncNavStoIntegrator(Coefficient &mu,
                                         VectorCoefficient &force,
                                         VectorCoefficient &sol)
   : c_mu(mu), c_force(force), c_sol(sol)
{
   dim = force.GetVDim();
   u.SetSize(dim);
   dudt.SetSize(dim);
   f.SetSize(dim);
   res_m.SetSize(dim);
   up.SetSize(dim);
   traction.SetSize(dim);
   grad_u.SetSize(dim);
   hess_u.SetSize(dim, (dim*(dim+1))/2);
   grad_p.SetSize(dim);
   Gij.SetSize(dim);
   nor.SetSize(dim);
   hn.SetSize(dim);
   hmap.SetSize(dim,dim);

   if (dim == 2)
   {
      hmap(0,0) = 0;
      hmap(0,1) =  hmap(1,0) =  1;
      hmap(1,1) = 2;
   }
   else if (dim == 3)
   {
      hmap(0,0) = 0;
      hmap(0,1) = hmap(1,0) = 1;
      hmap(0,2) = hmap(2,0) = 2;
      hmap(1,1) = 3;
      hmap(1,2) = hmap(2,1) = 4;
      hmap(2,2) = 5;
   }
   else
   {
      mfem_error("Only implemented for 2D and 3D");
   }
}

// Compute RBVMS stabilisation parameters
void IncNavStoIntegrator::GetTau(real_t &tau_m, real_t &tau_c, real_t &cfl2,
                                 real_t &mu, Vector &u,
                                 ElementTransformation &T)
{
   real_t Cd = 24.0;
   real_t Ct = 1.0;

   // Metric tensor
   MultAtB(T.InverseJacobian(),T.InverseJacobian(),Gij);

   // Temporal part
   tau_m = Ct/(dt*dt);

   // Convective part
   tau_c = 0.0;
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_c += Gij(i,j)*u[i]*u[j];
      }
   }
   cfl2 = tau_c/tau_m;
   tau_m += tau_c;

   // Diffusive part
   for (int j = 0; j < dim; j++)
   {
      for (int i = 0; i < dim; i++)
      {
         tau_m  += Cd*Cd*Gij(i,j)*Gij(i,j)*mu*mu;
      }
   }

   // Momentum stabilisation parameter
   tau_m = 1.0/sqrt(tau_m);

   // Continuity stabilisation parameter
   tau_c = 1.0/(tau_m*Gij.Trace());
}

// Compute Weak Dirichlet stabilisation parameters
void IncNavStoIntegrator::GetTauB(real_t &tau_b, real_t &tau_n,
                                  real_t &mu, Vector &u,
                                  Vector &nor,
                                  FaceElementTransformations &Tr)
{
   real_t Cb = 12.0;

   Tr.Elem1->InverseJacobian().Mult(nor,hn);
   tau_b = Cb*mu*hn.Norml2();
   tau_n = 0.0;
}

// Get energy
real_t IncNavStoIntegrator::GetElementEnergy(
   const Array<const FiniteElement *>&el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate)
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::GetElementEnergy"
                 " has incorrect block finite element space size!");
   }
   int dof_u = el[0]->GetDof();
   Vector tau(3);

   sh_u.SetSize(dof_u);
   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   real_t energy = 0.0;

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);

      real_t w = ip.weight * Tr.Weight();

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);

      energy += w*(u*u)/2;
   }

   return energy;
}

// Assemble the element interior residual vectors
void IncNavStoIntegrator::AssembleElementVector(
   const Array<const FiniteElement *> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array<Vector *> &elvec,
   real_t &elem_cfl)
{
   if (el.Size() != 2)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " has finite element space of incorrect block number");
   }

   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int spaceDim = Tr.GetSpaceDim();
   bool hess = false;//(el[0]->GetDerivType() == (int) FiniteElement::HESS);
   if (dim != spaceDim)
   {
      mfem_error("IncNavStoIntegrator::AssembleElementVector"
                 " is not defined on manifold meshes");
   }
   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   shh_u.SetSize(dof_u, (dim*(dim+1))/2);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, cfl2;
   elem_cfl = 0.0;
   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();
      real_t mu = c_mu.Eval(Tr, ip);
      c_force.Eval(f, Tr, ip);

      // Compute shape and interpolate
      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
                               hess_u(i,hmap(j,i)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }
      real_t res_c = grad_u.Trace();

      // Compute stability params
      GetTau(tau_m, tau_c, cfl2, mu, u, Tr);
      elem_cfl = fmax(elem_cfl, cfl2);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      u += up;
      p -= tau_c*res_c;

      // Compute momentum weak residual
      flux.Diag(-p, dim);                         // Add pressure
      grad_u.Symmetrize();                        // Grad to strain
      flux.Add(2*mu,grad_u);                      // Add stress to flux
      AddMult_a_VVt(-1.0, u, flux);               // Add convection to flux
      AddMult_a_ABt(w, shg_u, flux, elv_u);       // Add flux term to rhs
      f -= dudt;                                  // Add Acceleration to force
      AddMult_a_VWt(-w, sh_u, f, elv_u);          // Add force + acc term to rhs

      // Compute continuity weak residual
      elvec[1]->Add(-w*res_c, sh_p);              // Add Galerkin term
      shg_p.Mult(up, sh_p);                       // PSPG help term
      elvec[1]->Add(w, sh_p);                     // Add PSPG term
   }

   elem_cfl = sqrt(elem_cfl);
}






// // Assemble the element interior gradient matrices
// void IncNavStoIntegrator::AssembleElementGrad(
//    const Array<const FiniteElement*> &el,
//    ElementTransformation &Tr,
//    const Array<const Vector *> &elsol,
//    const Array<const Vector *> &elrate,
//    const Array2D<DenseMatrix *> &elmats)
// {
//    double ForloopIntroductionTime = 0;
//    double MomentumVelocityBlockTime = 0;
//    double MomentumPressureBlockTime = 0;
//    double ContinuityVelocityBlockTime = 0;
//    double ContinuityPressureBlockTime = 0;

//    double ForloopIntroductionTimeSum = 0;
//    double MomentumVelocityBlockTimeSum = 0;
//    double MomentumPressureBlockTimeSum = 0;
//    double ContinuityVelocityBlockTimeSum = 0;
//    double ContinuityPressureBlockTimeSum = 0;
//    double MomentumVelocityBlockCounter1 = 0;
//    double MomentumVelocityBlockCounter2 = 0;
//    double MomentumVelocityBlockCounter3 = 0;

//    int ForloopIntroductionCounter = 0;
//    int MomentumVelocityBlockCounter = 0;
//    int MomentumPressureBlockCounter = 0;
//    int ContinuityVelocityBlockCounter = 0;
//    int ContinuityPressureBlockCounter = 0;

//    // Start Measuring Time
//    auto TimeStart1 = std::chrono::high_resolution_clock::now();

//    int dof_u = el[0]->GetDof();
//    int dof_p = el[1]->GetDof();

//    bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

//    elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
//    elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);

//    elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
//    elmats(0,1)->SetSize(dof_u*dim, dof_p);
//    elmats(1,0)->SetSize(dof_p, dof_u*dim);
//    elmats(1,1)->SetSize(dof_p, dof_p);

//    *elmats(0,0) = 0.0;
//    *elmats(0,1) = 0.0;
//    *elmats(1,0) = 0.0;
//    *elmats(1,1) = 0.0;

//    sh_u.SetSize(dof_u);
//    shg_u.SetSize(dof_u, dim);
//    ushg_u.SetSize(dof_u);
//    dupdu.SetSize(dof_u);
//    sh_p.SetSize(dof_p);
//    shg_p.SetSize(dof_p, dim);

//    int intorder = 2*el[0]->GetOrder();
//    const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
//    real_t tau_m, tau_c, cfl2;

//    // Timing
//    auto TimeEnd1 = std::chrono::high_resolution_clock::now();
   
//    // For-loop dimensions:
//    // ir.GetNPoints(), dim, dof_u, dof_p
   
//    for (int i = 0; i < ir.GetNPoints(); ++i)
//    {
//       // Start Measuring Time
//       auto TimeStart2 = std::chrono::high_resolution_clock::now();

//       const IntegrationPoint &ip = ir.IntPoint(i);
//       Tr.SetIntPoint(&ip);
//       real_t w = ip.weight * Tr.Weight();

//       real_t mu = c_mu.Eval(Tr, ip);

//       el[0]->CalcPhysShape(Tr, sh_u);
//       elf_u.MultTranspose(sh_u, u);
//       elf_du.MultTranspose(sh_u, dudt);

//       el[0]->CalcPhysDShape(Tr, shg_u);
//       MultAtB(elf_u, shg_u, grad_u);

//       shg_u.Mult(u, ushg_u);

//       el[1]->CalcPhysShape(Tr, sh_p);
//       real_t p = sh_p*(*elsol[1]);

//       el[1]->CalcPhysDShape(Tr, shg_p);
//       shg_p.MultTranspose(*elsol[1], grad_p);

//       // Compute strong residual
//       MultAtB(elf_u, shg_u, grad_u);
//       grad_u.Mult(u,res_m);   // Add convection
//       res_m += dudt;          // Add acceleration
//       res_m += grad_p;        // Add pressure
//       res_m -= f;             // Subtract force

//       if (hess)               // Add diffusion
//       {
//          el[0]->CalcPhysHessian(Tr,shh_u);
//          MultAtB(elf_u, shh_u, hess_u);
//          for (int i = 0; i < dim; ++i)
//          {
//             for (int j = 0; j < dim; ++j)
//             {
//                res_m[j] -= mu*(hess_u(j,hmap(i,i)) +
//                                hess_u(i,hmap(j,i)));
//             }
//          }
//       }
//       else                   // No diffusion in strong residual
//       {
//          shh_u = 0.0;
//          hess_u = 0.0;
//       }

//       // Compute stability params
//       GetTau(tau_m, tau_c, cfl2, mu, u, Tr);

//       // Small scale reconstruction
//       up.Set(-tau_m,res_m);
//       u += up;

//       // Compute small scale jacobian
//       for (int j_u = 0; j_u < dof_u; ++j_u)
//       {
//          dupdu(j_u) = -tau_m*(sh_u(j_u) + ushg_u(j_u)*dt);
//       }

//       // Recompute convective gradient
//       MultAtB(elf_u, shg_u, grad_u);

//       auto TimeEnd2 = std::chrono::high_resolution_clock::now();
      
//       auto TimeStart3 = std::chrono::high_resolution_clock::now();

//       // Momentum - Velocity block (w,u)
//       for (int i_u = 0; i_u < dof_u; ++i_u)
//       {
//          for (int j_u = 0; j_u < dof_u; ++j_u)
//          {
//             // Diffusion
//             real_t mat = 0.0;
//             for (int dim_u = 0; dim_u < dim; ++dim_u)
//             {
//                mat += shg_u(i_u,dim_u)*shg_u(j_u,dim_u);
//                MomentumVelocityBlockCounter1++;
//             }
//             mat *= mu*dt;

//             // Acceleration
//             mat += sh_u(i_u)*sh_u(j_u);

//             // Convection -- frozen convection
//             mat -= ushg_u(i_u)*sh_u(j_u)*dt;           // Galerkin
//             mat -= ushg_u(i_u)*dupdu(j_u);             // SUPG

//             mat *= w;
//             for (int dim_u = 0; dim_u < dim; ++dim_u)
//             {
//                (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
//             }

//             for (int i_dim = 0; i_dim < dim; ++i_dim)
//             {
//                for (int j_dim = 0; j_dim < dim; ++j_dim)
//                {
//                   (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u)
//                   += mu*shg_u(i_u,j_dim)*shg_u(j_u,i_dim)*w*dt;
//                   MomentumVelocityBlockCounter2++;
//                }
//             }
//             for (int i_dim = 0; i_dim < dim; ++i_dim)
//             {
//                for (int j_dim = 0; j_dim < dim; ++j_dim)
//                {
//                   (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u)
//                   += tau_c*shg_u(i_u,i_dim)*shg_u(j_u,j_dim)*w*dt;
//                   MomentumVelocityBlockCounter3++;
//                }
//             }

//          }
//       }
//       auto TimeEnd3 = std::chrono::high_resolution_clock::now();

//       auto TimeStart4 = std::chrono::high_resolution_clock::now();
//       // Momentum - Pressure block (w,p) & // Continuity - Velocity block (q,u)
//       for (int j_u = 0; j_u < dof_u; ++j_u)
//       {
//          for (int i_p = 0; i_p < dof_p; ++i_p)
//          {
//             for (int dim_u = 0; dim_u < dim; ++dim_u)
//             {
//                // Momentum - Pressure block (w,p)
//                (*elmats(0,1))(j_u + dof_u * dim_u, i_p)
//                += (shg_p(i_p,dim_u)*tau_m*ushg_u(j_u)
//                    - shg_u(j_u,dim_u)*sh_p(i_p))*w*dt;
               
//                // MomentumPressureBlockCounter++;

//                // // Continuity - Velocity block (q,u)
//                // (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
//                // -= sh_p(i_p)*shg_u(j_u,dim_u)*w*dt;
//                // (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
//                // += shg_p(i_p, dim_u)*dupdu(j_u)*w;
//             }
//          }
//       }
//       auto TimeEnd4 = std::chrono::high_resolution_clock::now();      

//       auto TimeStart5 = std::chrono::high_resolution_clock::now();
//       // Continuity - Velocity block (q,u)
//       for (int j_u = 0; j_u < dof_u; ++j_u)
//       {
//          for (int i_p = 0; i_p < dof_p; ++i_p)
//          {
//             for (int dim_u = 0; dim_u < dim; ++dim_u)
//             {
//                (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
//                -= sh_p(i_p)*shg_u(j_u,dim_u)*w*dt;
//                (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
//                += shg_p(i_p, dim_u)*dupdu(j_u)*w;
//                // ContinuityVelocityBlockCounter++;
//             }
//          }
//       }
//       auto TimeEnd5 = std::chrono::high_resolution_clock::now();

//       auto TimeStart6 = std::chrono::high_resolution_clock::now();
      
//       // Continuity - Pressure block (w,p)
//       AddMult_a_AAt(-w*tau_m*dt, shg_p, *elmats(1,1));

//       auto TimeEnd6 = std::chrono::high_resolution_clock::now();
      
//       bool TimeMeasurement = false;

//       if (TimeMeasurement){
//          ForloopIntroductionCounter++;
//          ForloopIntroductionTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd2 - TimeStart2).count();
//          ForloopIntroductionTimeSum += ForloopIntroductionTime;

//          MomentumVelocityBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd3 - TimeStart3).count();
//          MomentumVelocityBlockTimeSum += MomentumVelocityBlockTime;

//          MomentumPressureBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd4 - TimeStart4).count();
//          MomentumPressureBlockTimeSum += MomentumPressureBlockTime;

//          ContinuityVelocityBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd5 - TimeStart5).count();
//          ContinuityVelocityBlockTimeSum += ContinuityVelocityBlockTime;

//          ContinuityPressureBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd6 - TimeStart6).count();
//          ContinuityPressureBlockTimeSum += ContinuityPressureBlockTime;
//          ContinuityPressureBlockCounter++;
//       }
//    }

//    bool printMeasurements = true;
//    if (printMeasurements){
//       auto IntroductionTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd1 - TimeStart1).count();
//       auto TotalTimeEnd = std::chrono::high_resolution_clock::now();
//       auto TotalTime = std::chrono::duration_cast<std::chrono::microseconds>(TotalTimeEnd - TimeStart1).count();

//       // std::cout << "ir.GetNPoints() = " << ir.GetNPoints() << ", dim = " << dim << ", dof_u = " << dof_u << ", dof_p = " << dof_p << std::endl;
//       // std::cout << "---------------------- AEG  Profile -------------------" << std::endl;
//       // std::cout << "Introduction Elapsed time: " << IntroductionTime << " microseconds (1 call)" << std::endl;
//       // std::cout << "For-loop Introduction Elapsed time: " << ForloopIntroductionTimeSum << " microseconds (" << ForloopIntroductionCounter << "calls) (Mean time per call: " << ForloopIntroductionTimeSum / ForloopIntroductionCounter << std::endl;
//       // std::cout << "Momentum Velocity Block Elapsed time: " << MomentumVelocityBlockTimeSum << " microseconds (" << MomentumVelocityBlockCounter1 << " & " << MomentumVelocityBlockCounter2 << "calls) (Mean time per call: " << MomentumVelocityBlockTimeSum / MomentumVelocityBlockCounter1 << std::endl;   
//       // std::cout << "Momentum Pressure Block Elapsed time: " << MomentumPressureBlockTimeSum << " microseconds (" << MomentumPressureBlockCounter << "calls) (Mean time per call: " << MomentumPressureBlockTimeSum / MomentumPressureBlockCounter << std::endl;
//       // std::cout << "Continuity Velocity Elapsed time: " << ContinuityVelocityBlockTimeSum << " microseconds (" << ContinuityVelocityBlockCounter << "calls) (Mean time per call: " << ContinuityVelocityBlockTimeSum / ContinuityVelocityBlockCounter << std::endl;
//       // std::cout << "Continuity Pressure Elapsed time: " << ContinuityPressureBlockTimeSum << " microseconds (" << ContinuityPressureBlockCounter << "calls) (Mean time per call: " << ContinuityPressureBlockTimeSum / ContinuityPressureBlockCounter << std::endl;
//       // std::cout << "Total AEG Elapsed time: " << TotalTime << " microseconds" << std::endl;
//       // std::cout << "Sum of Elapsed time: " << IntroductionTime + ForloopIntroductionTimeSum + MomentumVelocityBlockTimeSum + MomentumPressureBlockTimeSum + ContinuityVelocityBlockTimeSum + ContinuityPressureBlockTimeSum<< " microseconds" << std::endl;
//       // std::cout << "--------------------- Next AEG call ---------------------\n" << std::endl;

//       std::cout << IntroductionTime << "\t" << ForloopIntroductionTimeSum << "\t" << MomentumVelocityBlockTimeSum << "\t" << MomentumPressureBlockTimeSum << "\t" << ContinuityVelocityBlockTimeSum << "\t" << ContinuityPressureBlockTimeSum << "\t" << TotalTime << std::endl;
   
//    }

// }




// Assemble the element interior gradient matrices

void IncNavStoIntegrator::AssembleElementGrad(
   const Array<const FiniteElement*> &el,
   ElementTransformation &Tr,
   const Array<const Vector *> &elsol,
   const Array<const Vector *> &elrate,
   const Array2D<DenseMatrix *> &elmats)
{
   double ForloopIntroductionTime = 0;
   double MomentumVelocityBlockTime = 0;
   double MomentumPressureBlockTime = 0;
   double ContinuityVelocityBlockTime = 0;
   double ContinuityPressureBlockTime = 0;

   double ForloopIntroductionTimeSum = 0;
   double MomentumVelocityBlockTimeSum = 0;
   double MomentumPressureBlockTimeSum = 0;
   double ContinuityVelocityBlockTimeSum = 0;
   double ContinuityPressureBlockTimeSum = 0;
   double MomentumVelocityBlockCounter1 = 0;
   double MomentumVelocityBlockCounter2 = 0;
   double MomentumVelocityBlockCounter3 = 0;

   int ForloopIntroductionCounter = 0;
   int MomentumVelocityBlockCounter = 0;
   int MomentumPressureBlockCounter = 0;
   int ContinuityVelocityBlockCounter = 0;
   int ContinuityPressureBlockCounter = 0;
   
   // Start Measuring Time
   auto TimeStart1 = std::chrono::high_resolution_clock::now();

   // Same
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   bool hess = false;// = (el[0]->GetDerivType() == (int) FiniteElement::HESS);

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elf_du.UseExternalData(elrate[0]->GetData(), dof_u, dim);
   // Same


   // New way of dereferencing elmats
   DenseMatrix &mat_wu = *elmats(0,0);
   DenseMatrix &mat_wp = *elmats(0,1);
   DenseMatrix &mat_qu = *elmats(1,0);
   DenseMatrix &mat_qp = *elmats(1,1);

   DenseMatrix mat_wu1(dof_u, dof_u);

   mat_wu.SetSize(dof_u*dim, dof_u*dim);
   mat_wp.SetSize(dof_u*dim, dof_p);
   mat_qu.SetSize(dof_p, dof_u*dim);
   mat_qp.SetSize(dof_p, dof_p);

   mat_wu1 = 0.0;
   mat_wu = 0.0;
   mat_wp = 0.0;
   mat_qu = 0.0;
   mat_qp = 0.0;
   // New way of dereferencing elmats


   // Same
   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);
   // Same


   // Introduce Transpose of shg_u
   DenseMatrix shg_uT(dim, dof_u);
   // Introduce Transpose of shg_u


   // Same
   int intorder = 2*el[0]->GetOrder()-2;
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);
   real_t tau_m, tau_c, cfl2;
   // Same

   for (int i = 0; i < ir.GetNPoints(); ++i)
   {
      // Same
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      real_t w = ip.weight * Tr.Weight();

      real_t mu = c_mu.Eval(Tr, ip);

      el[0]->CalcPhysShape(Tr, sh_u);
      elf_u.MultTranspose(sh_u, u);
      elf_du.MultTranspose(sh_u, dudt);

      el[0]->CalcPhysDShape(Tr, shg_u);
      MultAtB(elf_u, shg_u, grad_u);

      shg_u.Mult(u, ushg_u);

      el[1]->CalcPhysShape(Tr, sh_p);
      real_t p = sh_p*(*elsol[1]);

      el[1]->CalcPhysDShape(Tr, shg_p);
      shg_p.MultTranspose(*elsol[1], grad_p);

      // Compute strong residual
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Mult(u,res_m);   // Add convection
      res_m += dudt;          // Add acceleration
      res_m += grad_p;        // Add pressure
      res_m -= f;             // Subtract force

      if (hess)               // Add diffusion
      {
         el[0]->CalcPhysHessian(Tr,shh_u);
         MultAtB(elf_u, shh_u, hess_u);
         for (int i = 0; i < dim; ++i)
         {
            for (int j = 0; j < dim; ++j)
            {
               res_m[j] -= mu*(hess_u(j,hmap(i,i)) + hess_u(i,hmap(j,i)));
            }
         }
      }
      else                   // No diffusion in strong residual
      {
         shh_u = 0.0;
         hess_u = 0.0;
      }

      // Compute stability params
      GetTau(tau_m, tau_c, cfl2, mu, u, Tr);

      // Small scale reconstruction
      up.Set(-tau_m,res_m);
      u += up;

      // Compute small scale jacobian
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         dupdu(j_u) = -tau_m*(sh_u(j_u) + ushg_u(j_u)*dt);
      }

      // Recompute convective gradient
      MultAtB(elf_u, shg_u, grad_u);

      shg_uT.Transpose(shg_u);
      // Same


      // Momentum - Velocity block (w,u)
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         real_t tmp1 = sh_u(j_u);
         real_t tmp2 = dupdu(j_u);
         for (int i_u = 0; i_u < dof_u; ++i_u)
         {

            // Diffusion
            real_t mat = 0.0;
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               //  mat += shg_u(i_u,dim_u)*shg_u(j_u,dim_u);
               mat += shg_uT(dim_u, i_u)*shg_uT(dim_u, j_u);
            }
            mat *= mu*dt;

            // Acceleration
            mat += sh_u(i_u)*tmp1;

            // Convection -- frozen convection
            mat -= ushg_u(i_u)*(tmp1*dt+tmp2);   // tmp1 = Galerkin tmp2 = SUPG

            mat_wu1(i_u, j_u) += mat*w;
         }
      }

      for (int j_dim = 0; j_dim < dim; ++j_dim)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            int j_dof = j_u + j_dim*dof_u;
            real_t tmp1 = shg_u(j_u,j_dim)*w*dt;

            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               real_t tmp2 = shg_uT(i_dim,j_u)*w*dt;
               for (int i_u = 0; i_u < dof_u; ++i_u)
               {
                  mat_wu(i_u + i_dim*dof_u, j_dof)
                  += mu*shg_u(i_u,j_dim)*tmp2 + tau_c*shg_u(i_u,i_dim)*tmp1;
               }
            }
         }
      }

      // Momentum - Pressure block (w,p)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         real_t tmp1 = sh_p(i_p)*w*dt;
         for (int dim_u = 0; dim_u < dim; ++dim_u)
         {
            real_t tmp2 = shg_p(i_p,dim_u)*tau_m*w*dt;
            for (int j_u = 0; j_u < dof_u; ++j_u)
            {
               mat_wp(j_u + dof_u * dim_u, i_p)
               += (tmp2*ushg_u(j_u) - shg_u(j_u,dim_u)*tmp1);
            }
         }
      }

      // Continuity - Velocity block (q,u)
      for (int dim_u = 0; dim_u < dim; ++dim_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            int j_dof = j_u + dof_u * dim_u;
            real_t tmp1 = shg_u(j_u,dim_u)*w;
            real_t tmp2 = dupdu(j_u)*w;
            for (int i_p = 0; i_p < dof_p; ++i_p)
            {
               mat_qu(i_p, j_dof)
               += (-sh_p(i_p)*tmp1*dt + shg_p(i_p, dim_u)*tmp2);
            }
         }
      }

      // Continuity - Pressure block (w,p)
      AddMult_a_AAt(-w*tau_m*dt, shg_p, mat_qp);

      // bool TimeMeasurement = false;

      // if (TimeMeasurement){
      //    ForloopIntroductionCounter++;
      //    ForloopIntroductionTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd2 - TimeStart2).count();
      //    ForloopIntroductionTimeSum += ForloopIntroductionTime;

      //    MomentumVelocityBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd3 - TimeStart3).count();
      //    MomentumVelocityBlockTimeSum += MomentumVelocityBlockTime;

      //    MomentumPressureBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd4 - TimeStart4).count();
      //    MomentumPressureBlockTimeSum += MomentumPressureBlockTime;

      //    ContinuityVelocityBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd5 - TimeStart5).count();
      //    ContinuityVelocityBlockTimeSum += ContinuityVelocityBlockTime;

      //    ContinuityPressureBlockTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd6 - TimeStart6).count();
      //    ContinuityPressureBlockTimeSum += ContinuityPressureBlockTime;
      //    ContinuityPressureBlockCounter++;
      // }
   }

   for (int dim_u = 0; dim_u < dim; ++dim_u)
   {
      for (int j_u = 0; j_u < dof_u; ++j_u)
      {
         for (int i_u = 0; i_u < dof_u; ++i_u)
         {
            mat_wu(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat_wu1(i_u, j_u);
         }
      }
   }

   bool printMeasurements = true;
   if (printMeasurements){
      double IntroductionTime = 0;
      // auto IntroductionTime = std::chrono::duration_cast<std::chrono::microseconds>(TimeEnd1 - TimeStart1).count();
      auto TotalTimeEnd = std::chrono::high_resolution_clock::now();
      auto TotalTime = std::chrono::duration_cast<std::chrono::microseconds>(TotalTimeEnd - TimeStart1).count();

      // std::cout << "ir.GetNPoints() = " << ir.GetNPoints() << ", dim = " << dim << ", dof_u = " << dof_u << ", dof_p = " << dof_p << std::endl;
      // std::cout << "---------------------- AEG  Profile -------------------" << std::endl;
      // std::cout << "Introduction Elapsed time: " << IntroductionTime << " microseconds (1 call)" << std::endl;
      // std::cout << "For-loop Introduction Elapsed time: " << ForloopIntroductionTimeSum << " microseconds (" << ForloopIntroductionCounter << "calls) (Mean time per call: " << ForloopIntroductionTimeSum / ForloopIntroductionCounter << std::endl;
      // std::cout << "Momentum Velocity Block Elapsed time: " << MomentumVelocityBlockTimeSum << " microseconds (" << MomentumVelocityBlockCounter1 << " & " << MomentumVelocityBlockCounter2 << "calls) (Mean time per call: " << MomentumVelocityBlockTimeSum / MomentumVelocityBlockCounter1 << std::endl;   
      // std::cout << "Momentum Pressure Block Elapsed time: " << MomentumPressureBlockTimeSum << " microseconds (" << MomentumPressureBlockCounter << "calls) (Mean time per call: " << MomentumPressureBlockTimeSum / MomentumPressureBlockCounter << std::endl;
      // std::cout << "Continuity Velocity Elapsed time: " << ContinuityVelocityBlockTimeSum << " microseconds (" << ContinuityVelocityBlockCounter << "calls) (Mean time per call: " << ContinuityVelocityBlockTimeSum / ContinuityVelocityBlockCounter << std::endl;
      // std::cout << "Continuity Pressure Elapsed time: " << ContinuityPressureBlockTimeSum << " microseconds (" << ContinuityPressureBlockCounter << "calls) (Mean time per call: " << ContinuityPressureBlockTimeSum / ContinuityPressureBlockCounter << std::endl;
      // std::cout << "Total AEG Elapsed time: " << TotalTime << " microseconds" << std::endl;
      // std::cout << "Sum of Elapsed time: " << IntroductionTime + ForloopIntroductionTimeSum + MomentumVelocityBlockTimeSum + MomentumPressureBlockTimeSum + ContinuityVelocityBlockTimeSum + ContinuityPressureBlockTimeSum<< " microseconds" << std::endl;
      // std::cout << "--------------------- Next AEG call ---------------------\n" << std::endl;

      std::cout << IntroductionTime << "\t" << ForloopIntroductionTimeSum << "\t" << MomentumVelocityBlockTimeSum << "\t" << MomentumPressureBlockTimeSum << "\t" << ContinuityVelocityBlockTimeSum << "\t" << ContinuityPressureBlockTimeSum << "\t" << TotalTime << std::endl;
   }
}






















// Assemble the outflow boundary residual vectors
void IncNavStoIntegrator
::AssembleOutflowVector(const Array<const FiniteElement *> &el1,
                        const Array<const FiniteElement *> &el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   *elvec[0] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor);

      real_t w = ip.weight * Tr.Weight();// instead???* 0.5;//

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      real_t un = u*nor;
      AddMult_a_VWt(w*un, sh_u, u, elv_u);
   }
}

// Assemble the outflow boundary gradient matrices
void IncNavStoIntegrator
::AssembleOutflowGrad(const Array<const FiniteElement *>&el1,
                      const Array<const FiniteElement *>&el2,
                      FaceElementTransformations &Tr,
                      const Array<const Vector *> &elsol,
                      const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   *elmats(0,0) = 0.0;

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      CalcOrtho(Tr.Jacobian(), nor);

      real_t w = ip.weight * Tr.Weight(); //
      //real_t w = ip.weight * 0.5;// instead???

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      real_t un = u*nor;

      // Momentum - Velocity block (w,u)
      for (int i_u = 0; i_u < dof_u; ++i_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            real_t mat = sh_u(i_u)*sh_u(j_u)*un*w*dt;

            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
            }
         }
      }
   }
}

// Assemble the weak Dirichlet BC boundary residual vectors
void IncNavStoIntegrator
::AssembleWeakDirBCVector(const Array<const FiniteElement *> &el1,
                          const Array<const FiniteElement *> &el2,
                          FaceElementTransformations &Tr,
                          const Array<const Vector *> &elsol,
                          const Array<Vector *> &elvec)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();

   elvec[0]->SetSize(dof_u*dim);
   elvec[1]->SetSize(dof_p);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);
   elv_u.UseExternalData(elvec[0]->GetData(), dof_u, dim);

   sh_u.SetSize(dof_u);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();

      real_t mu = c_mu.Eval(*Tr.Elem1, eip);
      c_sol.Eval(up, *Tr.Elem1, eip);

      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight();

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      up -= u;
      up.Neg();
      real_t un = up*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      MultAtB(elf_u, shg_u, grad_u);
      grad_u.Symmetrize();           // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      real_t p = sh_p*(*elsol[1]);

      GetTauB(tau_b, tau_n, mu, u, nor,Tr);

      // Traction
      grad_u.Mult(nor, traction);
      traction *= -2*mu;             // Consistency
      traction.Add(p, nor);          // Pressure
      traction.Add(tau_b,up);        // Penalty
      traction.Add(tau_n*un,nor);    // Penalty -- normal
      AddMult_a_VWt(w, sh_u, traction, elv_u);

      // Dual consistency
      MultVWt(nor,up, flux);
      flux.Symmetrize();
      AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

      // Continuity
      elvec[1]->Add(w*un, sh_p);
   }
}

// Assemble the weak Dirichlet BC boundary gradient matrices
void IncNavStoIntegrator
::AssembleWeakDirBCGrad(const
                        Array<const FiniteElement *>&el1,
                        const Array<const FiniteElement *>&el2,
                        FaceElementTransformations &Tr,
                        const Array<const Vector *> &elsol,
                        const Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el1[0]->GetDof();
   int dof_p = el1[1]->GetDof();

   elf_u.UseExternalData(elsol[0]->GetData(), dof_u, dim);

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(dof_u*dim, dof_p);
   elmats(1,0)->SetSize(dof_p, dof_u*dim);
   elmats(1,1)->SetSize(dof_p, dof_p);

   *elmats(0,0) = 0.0;
   *elmats(0,1) = 0.0;
   *elmats(1,0) = 0.0;
   *elmats(1,1) = 0.0;

   sh_u.SetSize(dof_u);
   shg_u.SetSize(dof_u, dim);
   ushg_u.SetSize(dof_u);
   dupdu.SetSize(dof_u);
   sh_p.SetSize(dof_p);
   shg_p.SetSize(dof_p, dim);

   int intorder = 2*el1[0]->GetOrder();
   const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), intorder);
   real_t tau_b, tau_n;
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      // Set the integration point in the face and the neighboring element
      Tr.SetAllIntPoints(&ip);

      // Access the neighboring element's integration point
      const IntegrationPoint &eip = Tr.GetElement1IntPoint();
      real_t mu = c_mu.Eval(*Tr.Elem1, eip);
      CalcOrtho(Tr.Jacobian(), nor);
      nor /= nor.Norml2();

      real_t w = ip.weight * Tr.Weight(); //
      //real_t w = ip.weight * 0.5;// instead???

      el1[0]->CalcPhysShape(*Tr.Elem1, sh_u);
      elf_u.MultTranspose(sh_u, u);

      //  real_t un = u*nor;

      el1[0]->CalcPhysDShape(*Tr.Elem1, shg_u);
      //  MultAtB(elf_u, shg_u, grad_u);
      //   grad_u.Symmetrize();  // Grad to strain

      el1[1]->CalcPhysShape(*Tr.Elem1, sh_p);
      //  real_t p = sh_p*(*elsol[1]);

      GetTauB(tau_b, tau_n, mu, u, nor,Tr);

      /*
            // Traction
            grad_u.Mult(nor, traction);
            traction *= -2*mu;             // Consistency
            traction.Add(p, nor);          // Pressure
            traction.Add(lambda,up);       // Penalty
            traction.Add(lambda_n*un,nor); // Penalty -- normal
            AddMult_a_VWt(w, sh_u, traction, elv_u);

            // Dual consistency
            MultVWt(nor,up, flux);
            flux.Symmetrize();
            AddMult_a_ABt(-w*2*mu, shg_u, flux, elv_u);

            // Continuity
            elvec[1]->Add(w*un, sh_p);

      */

      // Momentum - Velocity block (w,u)
      for (int i_u = 0; i_u < dof_u; ++i_u)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {

            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
               // Consistency 1
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + i_dim*dof_u)
                    -= mu*sh_u(i_u)*nor(j_dim)*shg_u(j_u,j_dim)*w*dt;
               }
               // Consistency 2
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u)
                    -= mu*sh_u(i_u)*nor(j_dim)*shg_u(j_u,i_dim)*w*dt;
               }

               // Dual Consistency 1
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + i_dim*dof_u)
                    -= mu*shg_u(i_u,j_dim)*nor(j_dim)*sh_u(j_u)*w*dt;
               }
               // Dual Consistency 2
               for (int j_dim = 0; j_dim < dim; ++j_dim)
               {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u)
                    -= mu*shg_u(i_u,i_dim)*nor(j_dim)*sh_u(j_u)*w*dt;
               }
            }

            // Penalty
            real_t mat = sh_u(i_u)*sh_u(j_u)*tau_b*w*dt;

            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,0))(i_u + dim_u*dof_u, j_u + dim_u*dof_u) += mat;
            }
         }
      }

      // Momentum - Pressure block (w,p)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(0,1))(j_u + dof_u * dim_u, i_p)
               += sh_u(j_u)*nor(dim_u)*sh_p(i_p)*w*dt;
            }
         }
      }

      // Continuity - Velocity block (q,u)
      for (int i_p = 0; i_p < dof_p; ++i_p)
      {
         for (int j_u = 0; j_u < dof_u; ++j_u)
         {
            for (int dim_u = 0; dim_u < dim; ++dim_u)
            {
               (*elmats(1,0))(i_p, j_u + dof_u * dim_u)
               += sh_p(i_p)*sh_u(j_u)*nor(dim_u)*w*dt;
            }
         }
      }
   }
}
