// This file is part of the RBVMS application. For more information and source
// code availability visit https://idoakkerman.github.io/
//
//   _____  ______      ____  __  _____
//   |  __ \|  _ \ \    / /  \/  |/ ____|
//   | |__) | |_) \ \  / /| \  / | (___
//   |  _  /|  _ < \ \/ / | |\/| |\___ \ 
//   | | \ \| |_) | \  /  | |  | |____) |
//   |_|  \_\____/   \/   |_|  |_|_____/
//
//
// RBVMS is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license.
//------------------------------------------------------------------------------

#include "mfem.hpp"
#include "coefficients.hpp"
#include "weakform.hpp"
#include "evolution.hpp"
#include "precon.hpp"
#include "monitor.hpp"

#include <sys/stat.h>
#include <chrono>

using namespace std;
using namespace mfem;

extern void printInfo();
extern void line(int len);

int main(int argc, char *argv[])
{
   // 1. Initialize MPI and HYPRE and print info
   Mpi::Init(argc, argv);
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();
   printInfo();

   // 2. Parse command-line options.
   OptionsParser args(argc, argv);

   // Mesh and discretization parameters
   const char *mesh_file = "../../mfem/data/inline-quad.mesh";
   const char *ref_file  = "";
   int order = 1;
   int ref_levels = 0;
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_file, "-rf", "--ref-file",
                  "File with refinement data");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order isoparametric space.");

   // Problem parameters
   Array<int> strong_bdr;
   Array<int> weak_bdr;
   Array<int> outflow_bdr;

   real_t mu_param = 1.0;
   const char *lib_file = "libfun.so";

   args.AddOption(&weak_bdr, "-wbc", "--weak-bdr",
                  "List of boundaries where Dirichelet BCs are enforced weakly."
                  "\n\t - Default: strong Dirichelet BCs");
   args.AddOption(&outflow_bdr, "-out", "--outflow-bdr",
                  "List of outflow boundaries.");
   args.AddOption(&lib_file, "-l", "--lib",
                  "Library file for case specific function definitions:\n\t"
                  " - Initial condition\n\t"
                  " - Boundary condition\n\t"
                  " - Forcing\n\t"
                  " - Diffusion\n\t");
   args.AddOption(&mu_param, "-m", "--mu",
                  "Sets the diffusion parameters, should be positive.");

   // Time stepping params
   int ode_solver_type = 35;
   real_t t_final = 10.0;
   real_t dt = 0.01;
   real_t dt_max = 1.0;
   real_t dt_min = 0.0001;

   real_t cfl_target = 2.0;
   real_t dt_gain = -1.0;

   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  ODESolver::Types.c_str());
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--dt",
                  "Time step.");
   args.AddOption(&dt_min, "-dtmn", "--dt-min",
                  "Minimum time step size.");
   args.AddOption(&dt_max, "-dtmx", "--dt-max",
                  "Maximum time step size.");
   args.AddOption(&cfl_target, "-cfl", "--cfl-target",
                  "CFL target .");
   args.AddOption(&dt_gain, "-dtg", "--dt-gain",
                  "Gain coefficient for time step adjustment.");

   // Solver parameters
   double GMRES_RelTol = 1e-3;
   int    GMRES_MaxIter = 500;
   double Newton_RelTol = 1e-3;
   int    Newton_MaxIter = 10;
   int    PreconCounter = 0; 

   args.AddOption(&GMRES_RelTol, "-lt", "--linear-tolerance",
                  "Relative tolerance for the GMRES solver.");
   args.AddOption(&GMRES_MaxIter, "-li", "--linear-itermax",
                  "Maximum iteration count for the GMRES solver.");
   args.AddOption(&Newton_RelTol, "-nt", "--newton-tolerance",
                  "Relative tolerance for the Newton solver.");
   args.AddOption(&Newton_MaxIter, "-ni", "--newton-itermax",
                  "Maximum iteration count for the Newton solver.");

   // Solution input/output params
   bool restart = false;
   int restart_interval = -1;
   real_t dt_vis = 10*dt;
   args.AddOption(&restart, "-rs", "--restart", "-f", "--fresh",
                  "Restart from solution.");
   args.AddOption(&restart_interval, "-ri", "--restart-interval",
                  "Interval between archieved time steps. \n\t"
                  "For negative values output is skipped.");
   args.AddOption(&dt_vis, "-dtv", "--dt_vis",
                  "Time interval between visualization points.");

   // Parse parameters
   args.Parse();
   if (!args.Good())
   {
      if (Mpi::Root()) { args.PrintUsage(cout); }
      return 1;
   }
   if (Mpi::Root()) { args.PrintOptions(cout); }

   // 3. Read the mesh from the given mesh file.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // Refine mesh
   {
      if (mesh.NURBSext && (strlen(ref_file) != 0))
      {
         mesh.RefineNURBSFromFile(ref_file);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
      if (Mpi::Root()) { mesh.PrintInfo(); }
   }

   // Partition mesh
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();

   // Select the time integrator
   unique_ptr<ODESolver> ode_solver = ODESolver::Select(ode_solver_type);
   int nstate = ode_solver->GetState() ? ode_solver->GetState()->MaxSize() : 0;

   if (nstate > 1 && ( restart || restart_interval > 0 ))
   {
      mfem_error("RBVMS restart not available for this time integrator \n"
                 "Time integrator can have a maximum of one statevector.");
   }

   // 4. Define a finite element space on the mesh.
   Array<FiniteElementCollection *> fecs(2);
   fecs[0] = FECollection::NewH1(order, dim, pmesh.IsNURBS());
   fecs[1] = FECollection::NewH1(order, dim, pmesh.IsNURBS());

   Array<ParFiniteElementSpace *> spaces(2);
   spaces[0] = new ParFiniteElementSpace(&pmesh, fecs[0],
                                         dim); //, Ordering::byVDIM);
   spaces[1] = new ParFiniteElementSpace(&pmesh, fecs[1]);

   // Report the degree of freedoms used
   {
      Array<int> tdof(num_procs),udof(num_procs),pdof(num_procs);
      tdof = 0;
      tdof[myid] = spaces[0]->TrueVSize();
      MPI_Reduce(tdof.GetData(), udof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      tdof = 0;
      tdof[myid] = spaces[1]->TrueVSize();
      MPI_Reduce(tdof.GetData(), pdof.GetData(), num_procs,
                 MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

      int udof_t = spaces[0]->GlobalTrueVSize();
      int pdof_t = spaces[1]->GlobalTrueVSize();
      if (Mpi::Root())
      {
         mfem::out << "Number of finite element unknowns:\n";
         mfem::out << "\tVelocity = "<< udof_t << endl;
         mfem::out << "\tPressure = "<< pdof_t << endl;
         mfem::out << "Number of finite element unknowns per partition:\n";
         mfem::out <<  "\tVelocity = "; udof.Print(mfem::out, num_procs);
         mfem::out <<  "\tPressure = "; pdof.Print(mfem::out, num_procs);
      }
   }

   // Get vector offsets
   Array<int> bOffsets(3);
   bOffsets[0] = 0;
   bOffsets[1] = spaces[0]->TrueVSize();
   bOffsets[2] = spaces[1]->TrueVSize();
   bOffsets.PartialSum();

   // 5. Define the time stepping algorithm

   // Set up the preconditioner
   RBVMS::JacobianPreconditioner jac_prec(bOffsets);
   std::cout << "Preconditioner has been made\n";
   // Set up the Jacobian solver
   RBVMS::GeneralResidualMonitor j_monitor(MPI_COMM_WORLD,"\t\tFGMRES", 25);
   FGMRESSolver j_gmres(MPI_COMM_WORLD);
   j_gmres.iterative_mode = false;
   j_gmres.SetRelTol(GMRES_RelTol);
   j_gmres.SetAbsTol(1e-12);
   j_gmres.SetMaxIter(GMRES_MaxIter);
   j_gmres.SetPrintLevel(-1);
   j_gmres.SetMonitor(j_monitor);
   j_gmres.SetPreconditioner(jac_prec);

   // Set up the Newton solver
   RBVMS::SystemResidualMonitor newton_monitor(MPI_COMM_WORLD,
                                               "Newton", 1,
                                               bOffsets);
   NewtonSolver newton_solver(MPI_COMM_WORLD);
   newton_solver.iterative_mode = true;
   newton_solver.SetPrintLevel(-1);
   newton_solver.SetMonitor(newton_monitor);
   newton_solver.SetRelTol(Newton_RelTol);
   newton_solver.SetAbsTol(1e-12);
   newton_solver.SetMaxIter(Newton_MaxIter );
   newton_solver.SetSolver(j_gmres);

   // Define the physical parameters
   LibVectorCoefficient sol(dim, lib_file, "sol_u");
   LibCoefficient mu(lib_file, "mu", false, mu_param);
   LibVectorCoefficient force(dim, lib_file, "force");

   // Define weak form and evolution
   RBVMS::IncNavStoIntegrator integrator(mu, force, sol);
   RBVMS::ParTimeDepBlockNonlinForm form(spaces, integrator);
   RBVMS::Evolution evo(form, newton_solver);
   ode_solver->Init(evo);

   // Boundary conditions
   int amax = spaces[0]->GetMesh()->bdr_attributes.Max();
   Array<bool> bnd_flag(amax+1);
   bnd_flag = false;

   // Check weak boundaries
   if (Mpi::Root() && weak_bdr.Size() > 0 )
   {
      cout<<"Weak    = "; weak_bdr.Print();
   }
   for (int b = 0; b < weak_bdr.Size(); b++)
   {
      if ( weak_bdr[b] < 0 || weak_bdr[b] > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flag[weak_bdr[b]])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flag[weak_bdr[b]] = true;
   }

   // Check outflow boundaries
   if (Mpi::Root() && outflow_bdr.Size() > 0)
   {
      cout<<"Outflow = "; outflow_bdr.Print();
   }
   for (int b = 0; b < outflow_bdr.Size(); b++)
   {
      if ( outflow_bdr[b] < 0 || outflow_bdr[b] > amax )
      {
         mfem_error("Boundary out of range.");
      }
      if (bnd_flag[outflow_bdr[b]])
      {
         mfem_error("Boundary specified more then once.");
      }
      bnd_flag[outflow_bdr[b]] = true;
   }

   // Assign strong boundaries
   strong_bdr.SetSize(amax - weak_bdr.Size() - outflow_bdr.Size());
   strong_bdr = -9;
   for (int b = 1, s = 0; b < amax+1; b++)
   {
      if (!bnd_flag[b])
      {
         strong_bdr[s] = b;
         s++;
      }
   }
   if (Mpi::Root() && strong_bdr.Size() > 0)
   {
      cout<<"Strong  = "; strong_bdr.Print();
   }

   // Set boundaries in the weakform
   form.SetStrongBC (strong_bdr);
   form.SetWeakBC   (weak_bdr);
   form.SetOutflowBC(outflow_bdr);

   // 6. Define the solution vector, grid function and output
   BlockVector xp(bOffsets);
   BlockVector dxp(bOffsets);
   BlockVector xp0(bOffsets);
   BlockVector xpi(bOffsets);

   // Define the gridfunctions
   ParGridFunction x_u(spaces[0]);
   ParGridFunction x_p(spaces[1]);

   Array<ParGridFunction*> dx_u(nstate);
   Array<ParGridFunction*> dx_p(nstate);

   for (int i = 0; i < nstate; i++)
   {
      dx_u[i] = new ParGridFunction(spaces[0]);
      dx_p[i] = new ParGridFunction(spaces[1]);
   }

   // Define the visualisation output
   VisItDataCollection vdc("step", &pmesh);
   vdc.SetPrefixPath("solution");
   vdc.RegisterField("u", &x_u);
   vdc.RegisterField("p", &x_p);

   // Define the restart output
   VisItDataCollection rdc("step", &pmesh);
   rdc.SetPrefixPath("restart");
   rdc.SetPrecision(18);

   // Get the start vector(s) from file -- or from function
   real_t t;
   int si, ri, vi;
   struct stat info;
   if (restart && stat("restart/step.dat", &info) == 0)
   {
      // Read

      if (Mpi::Root())
      {
         real_t dtr;
         std::ifstream in("restart/step.dat", std::ifstream::in);
         in>>t>>si>>ri>>vi;
         in>>dtr;
         in.close();
         cout<<"Restarting from step "<<ri-1<<endl;
         if (dt_gain > 0) dt = dtr;
      }
      // Synchronize
      MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&si, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&ri, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&vi, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // Open data files
      rdc.Load(ri-1);

      x_u = *rdc.GetField("u");
      x_p = *rdc.GetField("p");

      x_u.GetTrueDofs(xp.GetBlock(0));
      x_p.GetTrueDofs(xp.GetBlock(1));

      if (nstate == 1)
      {
         *dx_u[0] = *rdc.GetField("du");
         *dx_p[0] = *rdc.GetField("dp");

         dx_u[0]->GetTrueDofs(dxp.GetBlock(0));
         dx_p[0]->GetTrueDofs(dxp.GetBlock(1));

         ode_solver->GetState()->Append(dxp);
      }
   }
   else
   {
      // Define initial condition from file
      t = 0.0; si = 0; ri = 1; vi = 1;
      LibVectorCoefficient sol(dim, lib_file, "sol_u");
      x_u.ProjectCoefficient(sol);
      x_p = 0.0;

      x_u.GetTrueDofs(xp.GetBlock(0));
      x_p.GetTrueDofs(xp.GetBlock(1));

      // Visualize initial condition
      vdc.SetCycle(0);
      vdc.SetTime(0.0);
      vdc.Save();

      // Define the restart writer
      rdc.RegisterField("u", &x_u);
      rdc.RegisterField("p", &x_p);
      if (nstate == 1)
      {
         rdc.RegisterField("du", dx_u[0]);
         rdc.RegisterField("dp", dx_p[0]);
      }
   }

   // 7. Actual time integration

   // Open output file
   std::ofstream os;
   if (Mpi::Root())
   {
      std::ostringstream filename;
      filename << "output_"<<std::setw(6)<<setfill('0')<<si<< ".dat";
      os.open(filename.str().c_str());
   }

   // Loop till final time reached
   while (t < t_final)
   {
      // Print header
      if (Mpi::Root())
      {
         line(80);
         cout<<std::defaultfloat<<std::setprecision(4)<<std::setw(10);
         cout<<" step = " << si << " : dt = " << dt << endl;
         cout<<" time = [" << t << ", " << t+dt <<"]"<< endl;
         line(80);
      }

      // Actual time step
      xp0 = xp;
      auto newton_start = std::chrono::high_resolution_clock::now();
      ode_solver->Step(xp, t, dt);

      PreconCounter++;
      
      auto newton_end = std::chrono::high_resolution_clock::now();
      auto newton_duration = std::chrono::duration_cast<std::chrono::milliseconds>(newton_end - newton_start).count();
      if (Mpi::Root())
      {std::cout << std::endl <<"Time taken for one Time step: " << newton_duration/1000.0 << " seconds"<<  std::endl;}
      si++;
      //Reset the preconditioner for the next time step
      // if (PreconCounter % 20 == 0){
      //    jac_prec.ResetOperatorSetup();
      // }
      // Postprocess solution
      real_t cfl = evo.GetCFL();
      DenseMatrix bdrForce = evo.GetForce();
      if (Mpi::Root())
      {
         // Print to file
         int nbdr = bdrForce.Height();
         os << std::setw(10);
         os << si<<"\t"<<t<<"\t"<<dt<<"\t"<<cfl<<"\t";
         for (int b=0; b<nbdr; ++b)
            for (int v=0; v<dim; ++v)
            {
               os<<bdrForce(b,v)<<"\t";
            }
         os<<"\n"<< std::flush;

         // Print force to screen in table
         auto pline = [](int len)
         {
            cout<<" +";
            for (int b=0; b<len; ++b) { cout<<"-"; }
            cout<<"+\n";
         };

         // Print boundary header
         cout<<"\n";
         pline(10+13*nbdr);
         cout<<" | Boundary | ";
         for (int b=0; b<nbdr; ++b)
         {
            cout<<std::setw(10)<<b+1<<" | ";
         }
         cout<<"\n";
         pline(10+13*nbdr);

         // Print actual forces
         char dimName[] = "xyz";
         for (int v=0; v<dim; ++v)
         {
            cout<<" | Force "<<dimName[v]<<"  | ";
            for (int b=0; b<nbdr; ++b)
            {
               cout<<std::defaultfloat<<std::setprecision(4)<<std::setw(10);
               cout<<bdrForce(b,v)<<" | ";
            }
            cout<<"\n";
         }
         pline(10+13*nbdr);
         cout<<"\n"<<std::flush;
      }

      // Write visualization files
      while (t >= dt_vis*vi)
      {
         // Interpolate solution
         real_t fac = (t-dt_vis*vi)/dt;
         add (fac, xp0,(1.0-fac), xp, xpi);

         // Report to screen
         if (Mpi::Root())
         {
            line(80);
            cout << "Visit output: " <<vi << endl;
            cout << "        Time: " <<t-dt<<" "<<t-fac*dt<<" "<<t<<endl;
            line(80);
         }

         // Copy solution in grid functions
         x_u.Distribute(xpi.GetBlock(0));
         x_p.Distribute(xpi.GetBlock(1));

         // Actually write to file
         vdc.SetCycle(vi);
         vdc.SetTime(dt_vis*vi);
         vdc.Save();
         vi++;
      }

      // Change time step
      real_t dt0 = dt;
      if ((dt_gain > 0))
      {
         dt *= pow(cfl_target/cfl, dt_gain);
         dt = min(dt, dt_max);
         dt = max(dt, dt_min);
      }

      // Print cfl and dt to screen
      if (Mpi::Root())
      {
         line(80);
         cout<<" cfl = "<<cfl<<endl;
         cout<<" dt  = "<<dt0<<" --> "<<dt<<endl;
         line(80);
      }

      // Write restart files
      if (restart_interval > 0 && si%restart_interval == 0)
      {
         // Report to screen
         if (Mpi::Root())
         {
            line(80);
            cout << "Restart output:" << ri << endl;
            line(80);
         }

         // Copy solution in grid functions
         x_u.Distribute(xp.GetBlock(0));
         x_p.Distribute(xp.GetBlock(1));

         if (nstate == 1)
         {
            ode_solver->GetState()->Get(0,dxp);
            dx_u[0]->Distribute(dxp.GetBlock(0));
            dx_p[0]->Distribute(dxp.GetBlock(1));
         }

         // Actually write to file
         rdc.SetCycle(ri);
         rdc.SetTime(t);
         rdc.Save();
         ri++;

         // print meta file
         if (Mpi::Root())
         {
            std::ofstream step("restart/step.dat", std::ifstream::out);
            step<<t<<"\t"<<si<<"\t"<<ri<<"\t"<<vi<<endl;
            step<<dt<<endl;
            step.close();
         }
      }

      if (Mpi::Root()) { cout<<endl<<endl; }
   }
   os.close();

   // 8. Free the used memory.
   for (int i = 0; i < fecs.Size(); ++i)
   {
      delete fecs[i];
   }
   for (int i = 0; i < spaces.Size(); ++i)
   {
      delete spaces[i];
   }

   return 0;
}

