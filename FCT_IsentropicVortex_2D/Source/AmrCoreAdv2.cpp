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

// Functions in this file (member functions of AmrCoreAdv)

// initVortex() : function to set the initial conditions of isentropic vortex

// WriteCheckpointFile() : function to write checkpoint output file

// ReadCheckpointFile () : function to read checkpoint output file

// GotoNextLine() : utility to skip to next line in Header

// WriteErrFile(), ErrFileName(), ErrFileMF(), ErrFileVarNames() : Writing plotfiles for the exact solution and error to disk

// PlotFileName(), PlotFileMF(), PlotFileVarNames(), WritePlotFile() : Writing plotfile of conserved variables to disk

// ---------------------------------------------------------------------------------------------
void 
AmrCoreAdv::initVortex (int &lev, Box const& bx, Array4<Real> const& a, const Geometry& geom)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);

   if(probtag == 1){
   	alpha = 0.25*pi;
   	Mach_fs = std::sqrt(2.0/gamma);
   	rofs = 1.0;
   	pfs = 1.0;	Tfs = 1.0;
   	R = 1.0; sigma = 1.0;
   	beta = 0.25*Mach_fs*5.0*std::sqrt(2.0)*std::exp(0.5)/pi;
   	ss_fs = std::sqrt(gamma*pfs/rofs);
   }
      ufs = Mach_fs*std::cos(alpha);
      vfs = Mach_fs*std::sin(alpha);
      // ufs = 0.0;
      // vfs = 0.0;

   for(int k = lo.z; k <= hi.z; ++k){
   		for(int j = lo.y; j <= hi.y; ++j){
   			for(int i = lo.x; i <= hi.x; ++i){
   				Real x = geom.ProbLo(0) + (i + 0.5)*geom.CellSize(0);
   				Real y = geom.ProbLo(1) + (j + 0.5)*geom.CellSize(1);

   				Real fxy = -0.5*( std::pow( (x-xcm)/R,2.0 ) + std::pow( (y-ycm)/R,2.0 ) )/(std::pow(sigma,2.0));
   				Real omega = beta*std::exp(fxy);
   				Real du = -omega*(y-ycm)/R;
   				Real dv = omega*(x-xcm)/R;
   				Real dT = -0.5*(gamma-1.0)*std::pow(omega,2.0)/std::pow(sigma,2.0);

   				a(i,j,k,ro) = std::pow(1+dT,c1)*rofs;
   				a(i,j,k,rou) = a(i,j,k,ro)*(ufs + du)*ss_fs;
   				a(i,j,k,rov) = a(i,j,k,ro)*(vfs + dv)*ss_fs;
   				a(i,j,k,pre) = pfs*std::pow(1.0 + dT, gamma/(gamma-1))/gamma;
   				a(i,j,k,roE) = a(i,j,k,pre)*c1 + 0.5*( ( std::pow(a(i,j,k,rou),2.0) 
   							 + std::pow(a(i,j,k,rov),2.0) )/a(i,j,k,ro) );
               a(i,j,k,ent) = a(i,j,k,pre)/std::pow(a(i,j,k,ro),gamma);
   			}
   		}
   }

}


void
AmrCoreAdv::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                     std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreAdv\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

}


void
AmrCoreAdv::ReadCheckpointFile ()
{

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int nc = 6;
        int ng = 4;
        phi_old[lev].define(grids[lev], dmap[lev], nc, ng);
        phi_new[lev].define(grids[lev], dmap[lev], nc, ng);
        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}

// utility to skip to next line in Header
void
AmrCoreAdv::GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}

// Writing plotfiles for the exact solution and error to disk
// get plotfile name
// write errorfile to disk
void
AmrCoreAdv::WriteErrFile () const
{
    BoxArray ba = exact.boxArray();
    const DistributionMapping& dm = exact.DistributionMap();
    MultiFab diff(ba,dm,1,0);
    // MultiFab::Copy(diff, exact, 1, 0, 1, 0);
    MultiFab::Copy(diff, exact, 3, 0, 1, 0);

    const std::string& errfilename = ErrFileName(istep[0]);
    const auto& mferr = ErrFileMF();
    const auto& errnames = ErrFileVarNames();

    amrex::Print() << "\nMax-norm of the error is " << diff.norm0() << "\n";

    amrex::Print() << "Writing errorfile " << errfilename << "\n";

    amrex::WriteMultiLevelPlotfile(errfilename, 1, mferr, errnames,
                   Geom(), t_new[0], istep, refRatio());
}
std::string
AmrCoreAdv::ErrFileName (int lev) const
{
    return amrex::Concatenate(err_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreAdv::ErrFileMF () const
{
    Vector<const MultiFab*> r;
    // r.push_back(&err);
    r.push_back(&exact);
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::ErrFileVarNames () const
{
    return {"roexact","rouexact","rovexact","roerr","rouerr","roverr"};
    // return {"roexact","roerr"};
}
//--------------------------------------------------------------------------
// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreAdv::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
    r.push_back(&phi_new[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::PlotFileVarNames () const
{
    return {"ro","rou","rov","roE","pre","ent"};
    // return {"ro"};
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                   Geom(), t_new[0], istep, refRatio());
}
// ---------------------------------------------------------------------------
// Function to fill physical domain boundary (fill ghost cells)
void 
AmrCoreAdv::FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc)
{
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
        CpuBndryFuncFab cpu_bndry_func(dummy_cpu_fill_extdir);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf(phi, 0, phi.nComp(), phi.nGrowVect(), 0.0, 0);
        // Print(myproc) << "rank= " << myproc << ", reached CpuBndryFuncFab()" << "\n";
    }
#endif
}

void 
AmrCoreAdv::dummy_cpu_fill_extdir (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
        // do something for external Dirichlet (BCType::ext_dir) if there are
}

// function to calculate pressure and entropy from conservative quantities (for testing purposes)
void 
AmrCoreAdv::CalcAuxillary (int lev, Box const& bx, Array4<Real> const& a, const Geometry& geom)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);
   const int nf = a.nComp();

   for(int k = lo.z; k <= hi.z; ++k){
      for(int j = lo.y; j <= hi.y; ++j){
         for(int i = lo.x; i <= hi.x; ++i){
            a(i,j,k,pre) = (gamma-1)*( a(i,j,k,roE)
            -   0.5*( ( pow(a(i,j,k,rou),2) + pow(a(i,j,k,rov),2) )/a(i,j,k,ro) ) );
            a(i,j,k,ent) = a(i,j,k,pre)/std::pow(a(i,j,k,ro),gamma);
         }
      }
   }

}
