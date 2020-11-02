#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_BCUtil.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include <AmrCoreAdv.H>
#include <AmrCoreAdv_F.H>

using namespace amrex;

// ---------------------------------------------------------------------------
// Function to fill physical domain boundary (fill ghost cells)
void 
AmrCoreAdv::FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc, Real cur_time)
{
    BL_PROFILE_VAR("AmrCoreAdv::FillDomainBoundary()", dombndry);
    int myproc = ParallelDescriptor::MyProc();
    // Print(myproc) << "rank= " << myproc << ", entered AmrCoreAdv::FillDomainBoundary()" << "\n";
    if (geom.isAllPeriodic()) return;
    if (phi.nGrow() == 0) return;

    AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

#if !(defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA) && defined(AMREX_GPU_PRAGMA_NO_HOST))
    if (Gpu::inLaunchRegion())
    {
#endif  
        GpuBndryFuncFab<dummy_gpu_fill_extdir> gpu_bndry_func(AmrCoreAdv::dummy_gpu_fill_extdir{});
        PhysBCFunct<GpuBndryFuncFab<dummy_gpu_fill_extdir> > physbcf
            (geom, bc, gpu_bndry_func);
        physbcf(phi, 0, phi.nComp(), phi.nGrowVect(), 0.0, 0);
        // Print(myproc) << "rank= " << myproc << ", reached GpuBndryFuncFab()" << "\n";
#if !(defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA) && defined(AMREX_GPU_PRAGMA_NO_HOST))
    }
    else
    {
        CpuBndryFuncFab cpu_bndry_func(outletBC_partialwall);;
        // CpuBndryFuncFab cpu_bndry_func(outletBC_hoextrap);;
        // CpuBndryFuncFab cpu_bndry_func(nullptr);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf(phi, 0, phi.nComp(), phi.nGrowVect(), cur_time, 0);
        // Print(myproc) << "rank= " << myproc << ", reached CpuBndryFuncFab()" << "\n";
    }
#endif
}

// This is a boundary condition where part of the boundary is wall 
// and part is outlet (first order extrapolation)
void 
AmrCoreAdv::outletBC_partialwall (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    int myproc = ParallelDescriptor::MyProc();
    const int ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, mach = 5;
    const Real gamma = 1.4;

    Real pfs = 0.0, roin = 0.0, pin = 0.0, uin = 0.0, vin = 0.0, ltube = 0.0, lwall = 0.0, loffset = 0.0, inflow_time = 1000.0;
    int ncellw = 5, tagprob;
    ParmParse pp("prob");
    pp.get("p2",pin);
    pp.get("ro2",roin);
    pp.get("u2",uin);
    pp.get("v2",vin);
    pp.get("pfs",pfs);
    pp.get("probtag",tagprob);

    pp.query("ncellw",ncellw);

    Real vmodin, ptot, alf, ssin, riem1;

    // Print() << "pin = " << pin << ", roin= " << roin << ", uin= " << uin << ", vin= " << vin
            // << "\n";

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const Box qbx(dest);
    const auto qlo = amrex::lbound(qbx);
    const auto qhi = amrex::ubound(qbx);

    const Box& domain = geom.Domain();
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    const int is = std::max(qlo.x, domlo.x);
    const int ie = std::min(qhi.x, domhi.x);
    const int ilo = domlo.x;
    const int ihi = domhi.x;

    Real xlo = geom.ProbLo(0);
    Real xhi = geom.ProbHi(0);
    Real dx  = geom.CellSize(0);

    Real lenwall = xlo + (ncellw - 0.5)*dx;
    // Print() << "lenwall= " << lenwall << "\n";
#if AMREX_SPACEDIM >= 2
    const int js = std::max(qlo.y, domlo.y);
    const int je = std::min(qhi.y, domhi.y);
    const int jlo = domlo.y;
    const int jhi = domhi.y;

    Real ylo = geom.ProbLo(1);
    Real yhi = geom.ProbHi(1);
    Real dy  = geom.CellSize(1);
#endif

    if(tagprob == 9 || tagprob == 10 || tagprob == 11){
      vmodin = std::sqrt(std::pow(uin,2.0) + std::pow(vin,2.0));
      ptot = pin + (0.5*roin*std::pow(vmodin,2.0));
      alf = 0.0;
      pp.query("alf",alf);
      ssin = std::sqrt(gamma*pin/roin);
      riem1 = vmodin + (2.0*ssin/(gamma-1));    
    }

    if(tagprob == 10 || tagprob == 11){
      // tube is 5 cells long in y-direction and wall is 2 cells aove the tube
      ltube = ylo + 5.0*dy;
      lwall = 2.0*dy;
      loffset = 5.0*dy;
      pp.query("ltube",ltube); 
      pp.query("lwall",lwall); 
      pp.query("loffset",loffset);
      inflow_time = 1000.0; 
      pp.query("inflow_time",inflow_time); 
      // Print(0) << "time= " << time << ", inflow_time= " << inflow_time << "\n";
    }

    Array4<Real> q(dest);
    BCRec const& bc = bcr[pre];

    // Print() << "ncomp= " << numcomp << "\n";

    if (lo.x < ilo) {
      const int imin = lo.x;
      const int imax = ilo-1;
      if (bc.lo(0) == BCType::ext_dir) {
        if(tagprob == 7){
          for(int k = lo.z; k <= hi.z; ++k){
            for(int j = lo.y; j <= hi.y; ++j){
              for(int i = imin; i <= imax; ++i){
                // first, set pressure
                dest(i,j,k,pre) = dest(ilo,j,k,pre);
                // next, set mass density
                dest(i,j,k,ro) = roin*(std::pow(dest(ilo,j,k,pre)/pin,1.0/gamma));
                // set momentum density
                dest(i,j,k,rou) = roin*uin;
                dest(i,j,k,rov) = roin*vin;
                Real temp = std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0);
                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1.0) + (0.5*temp/dest(i,j,k,ro));
                dest(i,j,k,mach) = std::sqrt(temp)/dest(i,j,k,ro);
              }
            }
          }
        }
      else if(tagprob == 8 || tagprob == 9){
          for(int k = lo.z; k <= hi.z; ++k){
            for(int j = lo.y; j <= hi.y; ++j){
              for(int i = imin; i <= imax; ++i){
                // vmod, sstmp, riem2 calculated at i = ilo
                Real vmod = std::sqrt(std::pow(dest(ilo,j,k,rou),2.0) + std::pow(dest(ilo,j,k,rov),2.0))/dest(ilo,j,k,ro);
                Real sstmp = std::sqrt(gamma*dest(ilo,j,k,pre)/dest(ilo,j,k,ro));
                Real riem2 = vmod - (2.0*sstmp/(gamma-1));
                // get mod of velocity at ghost cell uisng riem1 and riem2
                // riem1 is calculated based on inlet conditions, riem2 from i = ilo cell
                Real q = 0.5*(riem1 + riem2); // velmod at ghost cell
                Real ss = 0.25*(gamma-1)*(riem1 - riem2); // sound speed at ghost cell
                // First calculate mach number at ghost cells
                dest(i,j,k,mach) = q/ss;
                Real k1 = std::pow(1.0+(0.5*(gamma-1)*std::pow(dest(i,j,k,mach),2.0)), gamma/(gamma-1));
                dest(i,j,k,pre) = ptot/k1;
                dest(i,j,k,ro) = gamma*dest(i,j,k,pre)/std::pow(ss,2.0);
                dest(i,j,k,rou) = dest(i,j,k,ro)*q*std::cos(alf);
                dest(i,j,k,rov) = dest(i,j,k,ro)*q*std::sin(alf);
                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1) 
                              + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
              }
            }
          }
        } 
      else if(tagprob == 10 || tagprob == 11){
          for(int k = lo.z; k <= hi.z; ++k){
            for(int j = lo.y; j <= hi.y; ++j){
              for(int i = imin; i <= imax; ++i){
                Real y = ylo + (j + 0.5)*dy;
                Real dist_off = (y + 0.5*dy) - (loffset);
                Real dist_tube = (y + 0.5*dy) - (loffset + ltube);
                Real dist_w1 = (y + 0.5*dy) - (loffset - lwall);
                Real dist_w2 = (y + 0.5*dy) - (loffset + lwall);
                if(dist_off > 0.0 && dist_tube <= 0.0 && time <= inflow_time){
                  // vmod, sstmp, riem2 calculated at i = ilo
                  Real vmod = std::sqrt(std::pow(dest(ilo,j,k,rou),2.0) + std::pow(dest(ilo,j,k,rov),2.0))/dest(ilo,j,k,ro);
                  Real sstmp = std::sqrt(gamma*dest(ilo,j,k,pre)/dest(ilo,j,k,ro));
                  Real riem2 = vmod - (2.0*sstmp/(gamma-1));
                  // get mod of velocity at ghost cell uisng riem1 and riem2
                  // riem1 is calculated based on inlet conditions, riem2 from i = ilo cell
                  Real q = 0.5*(riem1 + riem2); // velmod at ghost cell
                  Real ss = 0.25*(gamma-1)*(riem1 - riem2); // sound speed at ghost cell
                  // First calculate mach number at ghost cells
                  dest(i,j,k,mach) = q/ss;
                  Real k1 = std::pow(1.0+(0.5*(gamma-1)*std::pow(dest(i,j,k,mach),2.0)), gamma/(gamma-1));
                  dest(i,j,k,pre) = ptot/k1;
                  dest(i,j,k,ro) = gamma*dest(i,j,k,pre)/std::pow(ss,2.0);
                  dest(i,j,k,rou) = dest(i,j,k,ro)*q*std::cos(alf);
                  dest(i,j,k,rov) = dest(i,j,k,ro)*q*std::sin(alf);
                  dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1) 
                              + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
                } 
                else if(dist_off > 0.0 && dist_tube <= 0.0 && time > inflow_time){
                  // dest(i,j,k,ro) = dest(ilo+(ilo-i)-1,j,k,ro);
                  // dest(i,j,k,rov) = dest(ilo+(ilo-i)-1,j,k,rov);
                  // dest(i,j,k,roE) = dest(ilo+(ilo-i)-1,j,k,roE);
                  // dest(i,j,k,pre) = dest(ilo+(ilo-i)-1,j,k,pre);
                  // dest(i,j,k,mach) = dest(ilo+(ilo-i)-1,j,k,mach);
                  // dest(i,j,k,rou) = -dest(ilo+(ilo-i)-1,j,k,rov); 
                  for(int n = 0; n < numcomp; ++n){
                    dest(i,j,k,n) = q(ilo,j,k,n);
                  }
                }
                else if( (dist_w1 > 0.0 && dist_off <= 0.0) || (dist_off > 0.0 && dist_w2 <= 0.0 && dist_tube > 0.0) ){
                  // impose no-penetration wall BC
                  dest(i,j,k,ro) = dest(ilo+(ilo-i)-1,j,k,ro);
                  dest(i,j,k,rov) = dest(ilo+(ilo-i)-1,j,k,rov);
                  dest(i,j,k,roE) = dest(ilo+(ilo-i)-1,j,k,roE);
                  dest(i,j,k,pre) = dest(ilo+(ilo-i)-1,j,k,pre);
                  dest(i,j,k,mach) = dest(ilo+(ilo-i)-1,j,k,mach);
                  dest(i,j,k,rou) = -dest(ilo+(ilo-i)-1,j,k,rov);                
                } else{
                  for(int n = 0; n < numcomp; ++n){
                    dest(i,j,k,n) = q(ilo,j,k,n);
                  }
                } 
              }
            }
          }        
        } 
      }
    }    

#if AMREX_SPACEDIM >= 2

    if (hi.y > jhi) {
      const int jmin = jhi+1;
      const int jmax = hi.y;

      if (bc.hi(1) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int i = lo.x; i <= hi.x; ++i) {
              for (int j = jmin; j <= jmax; ++j) {
                Real x = xlo + (i + 0.5)*dx;
                if(x <= lenwall){
                  // impose no-penetration wall BC upto a certain length lenwall
                  dest(i,j,k,ro) = dest(i,jhi-(j-jhi)+1,k,ro);
                  dest(i,j,k,rou) = dest(i,jhi-(j-jhi)+1,k,rou);
                  dest(i,j,k,roE) = dest(i,jhi-(j-jhi)+1,k,roE);
                  dest(i,j,k,pre) = dest(i,jhi-(j-jhi)+1,k,pre);
                  dest(i,j,k,mach) = dest(i,jhi-(j-jhi)+1,k,mach);
                  dest(i,j,k,rov) = -dest(i,jhi-(j-jhi)+1,k,rov);
                }else{
                  for(int n = 0; n < numcomp; ++n){
                    // first-order extrapolation in other places
                    dest(i,j,k,n) = dest(i,jhi,k,n);
                  }
                }
              }
            }

          }

        }
      }  
#endif  

}
