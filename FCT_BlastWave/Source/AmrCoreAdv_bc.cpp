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
AmrCoreAdv::FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc)
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
        // CpuBndryFuncFab cpu_bndry_func(outletBC_riemann);;
        // CpuBndryFuncFab cpu_bndry_func(outletBC_hoextrap);;
        CpuBndryFuncFab cpu_bndry_func(nullptr);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf(phi, 0, phi.nComp(), phi.nGrowVect(), 0.0, 0);
        // Print(myproc) << "rank= " << myproc << ", reached CpuBndryFuncFab()" << "\n";
    }
#endif
}


void 
AmrCoreAdv::outletBC_riemann (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    int myproc = ParallelDescriptor::MyProc();
    const int ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, mach = 5;
    const Real gamma = 1.4;

    Real pfs = 0.0;
    ParmParse pp("prob");
    pp.get("pfs",pfs);

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

#if AMREX_SPACEDIM >= 2
    const int js = std::max(qlo.y, domlo.y);
    const int je = std::min(qhi.y, domhi.y);
    const int jlo = domlo.y;
    const int jhi = domhi.y;
#endif

    Array4<Real> q(dest);
    BCRec const& bc = bcr[pre];

    // Print() << "ncomp= " << numcomp << "\n";

    if (lo.x < ilo) {
      const int imin = lo.x;
      const int imax = ilo-1;

      if (bc.lo(0) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            //  First, calculate Mach number at i = ilo
            Real vmod_ilo = std::sqrt( std::pow(dest(ilo,j,k,rou)/dest(ilo,j,k,ro),2.0) 
                  + std::pow(dest(ilo,j,k,rov)/dest(ilo,j,k,ro),2.0) );
            Real ss_ilo = std::sqrt(gamma*dest(ilo,j,k,pre)/dest(ilo,j,k,ro));
            dest(ilo,j,k,mach) = vmod_ilo/ss_ilo;

            if(dest(ilo,j,k,mach) > 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int i = imin; i <= imax; ++i) {
                for (int n = 0; n < numcomp; ++n){
                    dest(i,j,k,n) = dest(ilo,j,k,n);
                }
              }
            }else{
            // Subsonic outlet BC (based on Riemann invariants)
              for(int i = imin; i <= imax; ++i){
                 dest(i,j,k,pre) = pinf;
              }
            // y-component of velocity, 1st Riemann invariant and entropy are set 
            // based on the end cell value
              Real yvel_ilo = dest(ilo,j,k,rov)/dest(ilo,j,k,ro);
              Real ent_ilo = dest(ilo,j,k,pre)/(std::pow(dest(ilo,j,k,ro),gamma));
              for(int i = imin; i <= imax; ++i){
                dest(i,j,k,ro) = std::pow(dest(i,j,k,pre)/ent_ilo, 1.0/gamma);
              }
              Real riem1 = vmod_ilo + 2.0*ss_ilo/(gamma-1);
                for(int i = imin; i <= imax; ++i){
                    Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                    Real vmod_loc = riem1 - (2.0*ss_loc/(gamma-1));
                    if(vmod_loc < 0.0){
                    	vmod_loc = 0.0;
                      // Print(myproc) << "rank= " << myproc << ", vmod -ve = " << vmod_loc
                      // << "i = " << i << ", j= " << j << "\n";
                    } 
                    Real xvel = std::sqrt(std::pow(vmod_loc,2.0) - std::pow(yvel_ilo,2.0));
                	if(std::isnan(xvel)){
                		xvel = 0.0;
                	}                      
                    dest(i,j,k,rou) = dest(i,j,k,ro)*xvel;
                    dest(i,j,k,rov) = dest(i,j,k,ro)*yvel_ilo;
                    dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                    + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
                    dest(i,j,k,mach) = vmod_loc/ss_loc;
                  }

                }
                  // Print(myproc) << "rank= " << myproc << ", i = " << imax << ", j= " << j
                  // << ", pre = " << dest(imax,j,k,pre) << "\n";
              }
            }
          }
      }

        if (hi.x > ihi) {
            const int imin = ihi+1;
            const int imax = hi.x;

            if (bc.hi(0) == BCType::ext_dir) {
              // Print(myproc) << "rank= " << myproc << ", pressure xlobc= " << bc.lo(0) 
              // << ", pressure xhibc=" << bc.hi(0) << "\n";
              Real pinf = pfs;
              for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                  //  First, calculate Mach number at i = ihi
                  Real vmod_ihi = std::sqrt( std::pow(dest(ihi,j,k,rou)/dest(ihi,j,k,ro),2.0) 
                        + std::pow(dest(ihi,j,k,rov)/dest(ihi,j,k,ro),2.0) );
                  Real ss_ihi = std::sqrt(gamma*dest(ihi,j,k,pre)/dest(ihi,j,k,ro));
                  dest(ihi,j,k,mach) = vmod_ihi/ss_ihi;

                  if(dest(ihi,j,k,mach) >= 1.0){
                    // supersonic BC is only 1st order extrapolation for all quantities
                    for (int n = 0; n < numcomp; ++n){
                    	for (int i = imin; i <= imax; ++i) {
                          dest(i,j,k,n) = dest(ihi,j,k,n);
                      }
                    }
                  }else{

                    // Subsonic outlet BC (based on Riemann invariants)
                    for(int i = imin; i <= imax; ++i){
                      dest(i,j,k,pre) = pinf;
                    }
                    // y-component of velocity, 1st Riemann invariant and entropy are set 
                    // based on the end cell value
                    Real yvel_ihi = dest(ihi,j,k,rov)/dest(ihi,j,k,ro);
                    Real ent_ihi = dest(ihi,j,k,pre)/(std::pow(dest(ihi,j,k,ro),gamma));
                    for(int i = imin; i <= imax; ++i){
                      dest(i,j,k,ro) = std::pow(dest(i,j,k,pre)/ent_ihi, 1.0/gamma);
                    }
                    Real riem1 = vmod_ihi + (2.0*ss_ihi/(gamma-1));
                    for(int i = imin; i <= imax; ++i){
                      Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                      Real vmod_loc = riem1 - (2.0*ss_loc/(gamma-1));
                      if(vmod_loc < 0.0){
                      	vmod_loc = 0.0;
                        // Print(myproc) << "rank= " << myproc << ", vmod -ve = " << vmod_loc
                        // << "i = " << i << ", j= " << j << "\n";
                      } 
                      Real xvel = std::sqrt(std::pow(vmod_loc,2.0) - std::pow(yvel_ihi,2.0));
                	  if(std::isnan(xvel)){
                		xvel = 0.0;
                	  }         
                      dest(i,j,k,rou) = dest(i,j,k,ro)*xvel;
                      dest(i,j,k,rov) = dest(i,j,k,ro)*yvel_ihi;
                      dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                      + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
                      dest(i,j,k,mach) = vmod_loc/ss_loc;

                      // Print(myproc) << "subsonic BC, rank= " << myproc << ", i= " << i << ", j= " << j
                      // << ", ro= " << dest(i,j,k,ro) << ", rou= " << dest(i,j,k,rou)
                      // << ", rov= " << dest(i,j,k,rov) << "yvel_ihi= " << yvel_ihi <<
                      //  ", roE= " << dest(i,j,k,roE)
                      // << ", pre= " << dest(i,j,k,pre) << ", mach= " << dest(i,j,k,mach) << "\n"; 
                    }

                  }

                }
              }
            }
        }

#if AMREX_SPACEDIM >= 2

    if (lo.y < jlo) {
      const int jmin = lo.y;
      const int jmax = jlo-1;

      if (bc.lo(1) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int i = lo.x; i <= hi.x; ++i) {
            //  First, calculate Mach number at j = jlo
            Real vmod_jlo = std::sqrt(std::pow(dest(i,jlo,k,rou)/dest(i,jlo,k,ro),2.0) 
                  + std::pow(dest(i,jlo,k,rov)/dest(i,jlo,k,ro),2.0));
            Real ss_jlo = std::sqrt(gamma*dest(i,jlo,k,pre)/dest(i,jlo,k,ro));
            dest(i,jlo,k,mach) = vmod_jlo/ss_jlo;

            if(dest(i,jlo,k,mach) >= 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int j = jmin; j <= jmax; ++j) {
                for (int n = 0; n < numcomp; ++n){
                    dest(i,j,k,n) = dest(i,jlo,k,n);
                }
              }
            }else{
            // Subsonic outlet BC (based on Riemann invariants)
              for(int j = jmin; j <= jmax; ++j){
                 dest(i,j,k,pre) = pinf;
              }
            // y-component of velocity, 1st Riemann invariant and entropy are set 
            // based on the end cell value
              Real xvel_jlo = dest(i,jlo,k,rou)/dest(i,jlo,k,ro);
              Real ent_jlo = dest(i,jlo,k,pre)/(std::pow(dest(i,jlo,k,ro),gamma));
              for(int j = jmin; j <= jmax; ++j){
                dest(i,j,k,ro) = std::pow(dest(i,j,k,pre)/ent_jlo, 1.0/gamma);
              }
              Real riem1 = vmod_jlo + 2.0*ss_jlo/(gamma-1);

              for(int j = jmin; j <= jmax; ++j){
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                Real vmod_loc = riem1 - (2.0*ss_loc/(gamma-1));
                if(vmod_loc < 0.0){
                	vmod_loc = 0.0;
                  // Print(myproc) << "rank= " << myproc << ", vmod -ve = " << vmod_loc
                  //  << "i = " << i << ", j= " << j << "\n";
                } 
                Real yvel = std::sqrt(std::pow(vmod_loc,2.0) - std::pow(xvel_jlo,2.0));
                if(std::isnan(yvel)){
                	yvel = 0.0;
                }
                dest(i,j,k,rou) = dest(i,j,k,ro)*xvel_jlo;
                dest(i,j,k,rov) = dest(i,j,k,ro)*yvel;
                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;

                for(int n = 0; n < numcomp; ++n){
                  if(std::isnan(dest(i,j,k,n))){
                    Print(myproc) << "rank= " << myproc << ",Nan found, i= " << i << ",j= " << j
                    << ", n= " << n  << ", xvel_jlo= " << xvel_jlo << 
                    ", vmod_loc= " << vmod_loc << "yvel= " << yvel << ", rou= " << dest(i,j,k,rou) 
                    << ", vmod_jlo= " << vmod_jlo << ", riem1= " << riem1 << ", ss_jlo" << ss_jlo << "\n";
                  }
                }
              }

            }

          }
        }
      }

    if(lo.x < ilo){
      if(bc.lo(0) == BCType::reflect_even){
        const int imin = lo.x;
        const int imax = ilo-1;
        //apply symmetry BC at lower corners of domain based on BC setup
        for(int n = 0; n < numcomp; ++n){
          for(int k = lo.z; k <= hi.z; ++k){
            for(int j = jmin; j <= jmax; ++j){
              for(int i = imin; i <= imax; ++i){
                dest(i,j,k,n) = dest(ilo+(ilo-i)-1,j,k,n);
              }
            }
          }
        }
      }
    }  
  }

    if (hi.y > jhi) {
      const int jmin = jhi+1;
      const int jmax = hi.y;

      if (bc.hi(1) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int i = lo.x; i <= hi.x; ++i) {
            //  First, calculate Mach number at j = jlo
            Real vmod_jhi = std::sqrt( std::pow(dest(i,jhi,k,rou)/dest(i,jhi,k,ro),2.0) 
                  + std::pow(dest(i,jhi,k,rov)/dest(i,jhi,k,ro),2.0) );
            Real ss_jhi = std::sqrt(gamma*dest(i,jhi,k,pre)/dest(i,jhi,k,ro));
            dest(i,jhi,k,mach) = vmod_jhi/ss_jhi;

            if(dest(i,jhi,k,mach) >= 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int j = jmin; j <= jmax; ++j) {
                for (int n = 0; n < numcomp; ++n){
                    dest(i,j,k,n) = dest(i,jhi,k,n);
                }
              }
            }else{
            // Subsonic outlet BC (based on Riemann invariants)
              for(int j = jmin; j <= jmax; ++j){
                 dest(i,j,k,pre) = pfs;
              }
            // y-component of velocity, 1st Riemann invariant and entropy are set 
            // based on the end cell value
              Real xvel_jhi = dest(i,jhi,k,rou)/dest(i,jhi,k,ro);
              Real ent_jhi = dest(i,jhi,k,pre)/(std::pow(dest(i,jhi,k,ro),gamma));
              for(int j = jmin; j <= jmax; ++j){
                dest(i,j,k,ro) = std::pow(dest(i,j,k,pre)/ent_jhi, 1.0/gamma);
              }
              Real riem1 = vmod_jhi + 2.0*ss_jhi/(gamma-1);

              for(int j = jmin; j <= jmax; ++j){
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                Real vmod_loc = riem1 - (2.0*ss_loc/(gamma-1));
                if(vmod_loc < 0.0){
                	vmod_loc = 0.0;
                  // Print(myproc) << "rank= " << myproc << ", vmod -ve = " << vmod_loc
                  // << "i = " << i << ", j= " << j << "\n";
                } 
                Real yvel = std::sqrt(std::pow(vmod_loc,2.0) - std::pow(xvel_jhi,2.0));
                if(std::isnan(yvel)){
                		yvel = 0.0;
               	}  
                dest(i,j,k,rou) = dest(i,j,k,ro)*xvel_jhi;
                dest(i,j,k,rov) = dest(i,j,k,ro)*yvel;
                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }

            }

          }
        }
      }

    if(lo.x < ilo){
      if(bc.lo(0) == BCType::reflect_even){
        const int imin = lo.x;
        const int imax = ilo-1;
        //apply symmetry BC at lower corners of domain based on BC setup
        for(int n = 0; n < numcomp; ++n){
          for(int k = lo.z; k <= hi.z; ++k){
            for(int j = jmin; j <= jmax; ++j){
              for(int i = imin; i <= imax; ++i){
                dest(i,j,k,n) = dest(ilo+(ilo-i)-1,j,k,n);
              }
            }
          }
        }
      }
    }   
  }
#endif
    ParallelDescriptor::Barrier();   

}

void 
AmrCoreAdv::outletBC_hoextrap (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    int myproc = ParallelDescriptor::MyProc();
    const int ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, mach = 5;
    const Real gamma = 1.4;

    Real pfs = 0.0;
    ParmParse pp("prob");
    pp.get("pfs",pfs);

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

#if AMREX_SPACEDIM >= 2
    const int js = std::max(qlo.y, domlo.y);
    const int je = std::min(qhi.y, domhi.y);
    const int jlo = domlo.y;
    const int jhi = domhi.y;
#endif

    Array4<Real> q(dest);
    BCRec const& bc = bcr[pre];

    // Print() << "ncomp= " << numcomp << "\n";

    if (lo.x < ilo) {
      const int imin = lo.x;
      const int imax = ilo-1;

      if (bc.lo(0) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            //  First, calculate Mach number at i = ilo
            Real vmod_ilo = std::sqrt( std::pow(dest(ilo,j,k,rou)/dest(ilo,j,k,ro),2.0) 
                  + std::pow(dest(ilo,j,k,rov)/dest(ilo,j,k,ro),2.0) );
            Real ss_ilo = std::sqrt(gamma*dest(ilo,j,k,pre)/dest(ilo,j,k,ro));
            dest(ilo,j,k,mach) = vmod_ilo/ss_ilo;

            if(dest(ilo,j,k,mach) >= 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int i = imin; i <= imax; ++i) {
                  dest(i,j,k,pre) = dest(ilo,j,k,pre);

                  dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                  + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                  Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                    + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                  Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                  dest(i,j,k,mach) = vmod_loc/ss_loc;
              }
            }else{
            // Subsonic outlet BC (pressure set to freestream value)
              for(int i = imin; i <= imax; ++i){
                dest(i,j,k,pre) = pinf;

                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                  + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));

                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }
            }
          }
        }
      }
    }

    if (hi.x > ihi) {
      const int imin = ihi+1;
      const int imax = hi.x;

      if (bc.hi(0) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int j = lo.y; j <= hi.y; ++j) {
            //  First, calculate Mach number at i = ihi
            Real vmod_ihi = std::sqrt( std::pow(dest(ihi,j,k,rou)/dest(ihi,j,k,ro),2.0) 
                  + std::pow(dest(ihi,j,k,rov)/dest(ihi,j,k,ro),2.0) );
            Real ss_ihi = std::sqrt(gamma*dest(ihi,j,k,pre)/dest(ihi,j,k,ro));
            dest(ihi,j,k,mach) = vmod_ihi/ss_ihi;

            if(dest(ihi,j,k,mach) >= 1.0){
              for (int i = imin; i <= imax; ++i) {
                  dest(i,j,k,pre) = dest(ihi,j,k,pre);

                  dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                  + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                  Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                    + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                  Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                  dest(i,j,k,mach) = vmod_loc/ss_loc;
                }
              }else{
                // Subsonic outlet BC (based on Riemann invariants)
                for(int i = imin; i <= imax; ++i){
                  dest(i,j,k,pre) = pinf;

                  dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                  + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                  Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                    + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                  Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                  dest(i,j,k,mach) = vmod_loc/ss_loc;
                }

              }

            }
          }
        }
    }

#if AMREX_SPACEDIM >= 2

    if (lo.y < jlo) {
      const int jmin = lo.y;
      const int jmax = jlo-1;

      if (bc.lo(1) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int i = lo.x; i <= hi.x; ++i) {
            //  First, calculate Mach number at j = jlo
            Real vmod_jlo = std::sqrt( std::pow(dest(i,jlo,k,rou)/dest(i,jlo,k,ro),2.0) 
                  + std::pow(dest(i,jlo,k,rov)/dest(i,jlo,k,ro),2.0) );
            Real ss_jlo = std::sqrt(gamma*dest(i,jlo,k,pre)/dest(i,jlo,k,ro));
            dest(i,jlo,k,mach) = vmod_jlo/ss_jlo;

            if(dest(i,jlo,k,mach) > 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int j = jmin; j <= jmax; ++j) {
                dest(i,j,k,pre) = dest(i,jlo,k,pre);

                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                  + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }
            }else{
            // Subsonic outlet BC (based on Riemann invariants)
              for(int j = jmin; j <= jmax; ++j){
                dest(i,j,k,pre) = pinf;

                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                  + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }

            }

          }
        }
      } 
  }

    if (hi.y > jhi) {
      const int jmin = jhi+1;
      const int jmax = hi.y;

      if (bc.hi(1) == BCType::ext_dir) {
        Real pinf = pfs;
        for (int k = lo.z; k <= hi.z; ++k) {
          for (int i = lo.x; i <= hi.x; ++i) {
            //  First, calculate Mach number at j = jlo
            Real vmod_jhi = std::sqrt( std::pow(dest(i,jhi,k,rou)/dest(i,jhi,k,ro),2.0) 
                  + std::pow(dest(i,jhi,k,rov)/dest(i,jhi,k,ro),2.0) );
            Real ss_jhi = std::sqrt(gamma*dest(i,jhi,k,pre)/dest(i,jhi,k,ro));
            dest(i,jhi,k,mach) = vmod_jhi/ss_jhi;

            if(dest(i,jhi,k,mach) >= 1.0){
            // supersonic BC is only 1st order extrapolation for all quantities
              for (int j = jmin; j <= jmax; ++j) {
                dest(i,j,k,pre) = dest(i,jhi,k,pre);

                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                  + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }
            }else{
            // Subsonic outlet BC (based on Riemann invariants)
              for(int j = jmin; j <= jmax; ++j){
                dest(i,j,k,pre) = pfs;

                dest(i,j,k,roE) = dest(i,j,k,pre)/(gamma-1)
                + 0.5*((std::pow(dest(i,j,k,rou),2.0) + std::pow(dest(i,j,k,rov),2.0))/dest(i,j,k,ro));

                Real vmod_loc = std::sqrt( std::pow(dest(i,j,k,rou)/dest(i,j,k,ro),2.0)
                  + std::pow(dest(i,j,k,rov)/dest(i,j,k,ro),2.0) );
                Real ss_loc = std::sqrt(gamma*dest(i,j,k,pre)/dest(i,j,k,ro));
                dest(i,j,k,mach) = vmod_loc/ss_loc;
              }
            }

          }
        }
      }  
  }
#endif
    ParallelDescriptor::Barrier();   

}