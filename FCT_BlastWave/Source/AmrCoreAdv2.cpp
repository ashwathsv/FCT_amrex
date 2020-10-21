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

// initBlastWave() : function to set the initial conditions of blast wave problem

// WriteCheckpointFile() : function to write checkpoint output file

// ReadCheckpointFile () : function to read checkpoint output file

// GotoNextLine() : utility to skip to next line in Header

// WriteErrFile(), ErrFileName(), ErrFileMF(), ErrFileVarNames() : Writing plotfiles for the exact solution and error to disk

// PlotFileName(), PlotFileMF(), PlotFileVarNames(), WritePlotFile() : Writing plotfile of conserved variables to disk

// ---------------------------------------------------------------------------------------------
void 
AmrCoreAdv::initBlastWave (int &lev, Box const& bx, Array4<Real> const& a, const Geometry& geom)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);

   if(probtag != 2){
    xcm = geom.ProbLo(0);
    ycm = 0.5*(geom.ProbLo(1) + geom.ProbHi(1));
   }else{
    xcm = 0.5*(geom.ProbLo(0) + geom.ProbHi(0));
    ycm = 0.5*(geom.ProbLo(1) + geom.ProbHi(1));    
   }

   for(int k = lo.z; k <= hi.z; ++k){
   		for(int j = lo.y; j <= hi.y; ++j){
   			for(int i = lo.x; i <= hi.x; ++i){
   				Real x = geom.ProbLo(0) + (i + 0.5)*geom.CellSize(0);
   				Real y = geom.ProbLo(1) + (j + 0.5)*geom.CellSize(1);

   				Real dist = std::pow(x-xcm,2.0) + std::pow(y-ycm,2.0) - std::pow(rad_bw,2.0);
          Real du = 1E-10;
          Real dv = 1E-10;

          if(dist <= 0.0){
            a(i,j,k,ro) = ro2;
            a(i,j,k,rou) = ro2*(u2);
            a(i,j,k,rov) = ro2*(v2);
            // a(i,j,k,rou) = ro2*(u2+du);
            // a(i,j,k,rov) = ro2*(v2+dv);
            a(i,j,k,pre) = p2;
          }else{
            a(i,j,k,ro) = ro1;
            a(i,j,k,rou) = ro1*(u1);
            a(i,j,k,rov) = ro1*(v1);
            // a(i,j,k,rou) = ro1*(u1+du);
            // a(i,j,k,rov) = ro1*(v1+dv);
            a(i,j,k,pre) = p1;            
          }
            a(i,j,k,roE) = a(i,j,k,pre)*c1 + 0.5*( ( std::pow(a(i,j,k,rou),2.0) 
                   + std::pow(a(i,j,k,rov),2.0) )/a(i,j,k,ro) );
            Real ss = std::sqrt(gamma*a(i,j,k,pre)/a(i,j,k,ro));
            Real velmod = std::sqrt( std::pow(a(i,j,k,rou)/a(i,j,k,ro),2.0) 
                        + std::pow(a(i,j,k,rov)/a(i,j,k,ro),2.0) );
            a(i,j,k,mach) = velmod/ss;
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
    return {"ro","rou","rov","roE","pre","Mach"};
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
            Real ss = std::sqrt(gamma*a(i,j,k,pre)/a(i,j,k,ro));
            Real velmod = std::sqrt( std::pow(a(i,j,k,rou)/a(i,j,k,ro),2.0) 
                        + std::pow(a(i,j,k,rov)/a(i,j,k,ro),2.0) );
            a(i,j,k,mach) = velmod/ss;
         }
      }
   }

}

Real
AmrCoreAdv::get_gradp(int lev, Box const& validbox, Array4<Real const> const& a, int dir)
{
   Real gradpmax = 0.0;
   const auto lo = lbound(validbox);
   const auto hi = ubound(validbox);
   const int nf = a.nComp();
   Real gradp;
   const Geometry& geom1 = geom[lev];
   Real dx[2];
   dx[0] = geom1.CellSize(0);
   dx[1] = geom1.CellSize(1);
   int myproc = ParallelDescriptor::MyProc();

    if(dir == 0){
        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    gradp = 0.5*(a(i+1,j,k,pre) - a(i-1,j,k,pre))/dx[0];
                    if (fabs(gradp) > fabs(gradpmax)){
                        gradpmax = gradp;
                    }
                }
            }
        }
    }else if(dir == 1){
        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    gradp = 0.5*(a(i,j+1,k,pre) - a(i,j-1,k,pre))/dx[1];
                    if (fabs(gradp) > fabs(gradpmax)){
                        gradpmax = gradp;
                    }
                }
            }
        }        
    }
   // Print(ParallelDescriptor::MyProc()) << "rank= " << myproc << "gradromax= " << gradromax << "\n";
    return fabs(gradpmax);
}

void
AmrCoreAdv::GetProbeDets ()
{
    int myproc = ParallelDescriptor::MyProc();
    MultiFab& phi = phi_new[0];
    const Geometry& geom1 = geom[0];
    const Box& domain = geom1.Domain();
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    Print() << "domlo= " << domlo << ", domhi= " << domhi << "\n";

    proberank.resize(nprobes,-1);
    
    for(int n = 0; n < nprobes; ++n){
      Real xprobe = geom1.ProbLo(0) + (iprobe[n] + 0.5)*geom1.CellSize(0);
      Real yprobe = geom1.ProbLo(1) + (jprobe[n] + 0.5)*geom1.CellSize(1);
      if(xprobe < geom1.ProbLo(0) || xprobe > geom1.ProbHi(0) || 
         yprobe < geom1.ProbLo(1) || yprobe > geom1.ProbHi(1)  ){
          amrex::Error("Point for probe not found in domain, exiting...");
      }
      // find out which rank has the point to be probed
      for (MFIter mfi(phi); mfi.isValid(); ++mfi)
      {
        const Box& box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        //  Function to print results
        FArrayBox& mfab = phi[mfi];
        Array4<Real> const& a = mfab.array();
        for     (int k = lo.z; k <= hi.z; ++k) {
          for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
              if(i == iprobe[n] && j == jprobe[n]){
                Print(myproc) << "point is present in rank= " << myproc << ", lo.x= " << lo.x << 
                ", hi.x= " << hi.x << ", lo.y= " << lo.y << ", hi.y= " << hi.y << ", iprobe"
                << iprobe[n] << ", jprobe= " << jprobe[n] << "\n";
                proberank[n] = myproc;
              }
            }
          }
        }
      }
    }
    ParallelDescriptor::Barrier();
}



void
AmrCoreAdv::WriteProbeFile (int lev, Real cur_time, int stepnum)
{
    MultiFab& phi = phi_new[lev];
    const Geometry& geom1 = geom[lev];
    const Box& domain = geom1.Domain();
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    int myproc = ParallelDescriptor::MyProc();

    for(int n = 0; n < nprobes; ++n){
      std::string plotname = "Outputs_P";
      plotname = amrex::Concatenate(plotname, probtag, 1);
      plotname = plotname + "L";
      plotname = amrex::Concatenate(plotname, max_level, 1);
      plotname = plotname + "/varstimeseriespr";
      plotname = amrex::Concatenate(plotname, n, 2);
      std::string filename;
      filename = plotname + ".txt";

      // Print()<<"filename= " << filename << "\n";

      Real xprobe = geom1.ProbLo(0) + (iprobe[n] + 0.5)*geom1.CellSize(0);
      Real yprobe = geom1.ProbLo(1) + (jprobe[n] + 0.5)*geom1.CellSize(1);

      if(proberank[n] == myproc){
        std::ofstream ofs;
        // Print(myproc) << "rank= " << myproc << ", xprobe= " << xprobe << ", yprobe= " << yprobe << "\n";
        if(stepnum == 0 && cur_time == 0){
          ofs.open(filename, std::ofstream::out);
          Print(myproc,ofs) << "# time ro rou rov roE pre mach" << "\n";
          Print(myproc,ofs) << "# x = " << xprobe << "\n";
          Print(myproc,ofs) << "# y = " << yprobe << "\n";
        }else{
          ofs.open(filename, std::ofstream::app);
        }
        for (MFIter mfi(phi); mfi.isValid(); ++mfi){
          const Box& box = mfi.validbox();
          const auto lo = lbound(box);
          const auto hi = ubound(box);
          FArrayBox& mfab = phi[mfi];
          Array4<Real> const& a = mfab.array(); 
          Print(myproc, ofs).SetPrecision(8) << std::left << std::setw(12) << cur_time << "\t" 
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,ro)  << "\t" 
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,rou) << "\t"
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,rov) << "\t" 
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,roE) << "\t" 
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,pre) << "\t"
          << std::left << std::setw(12) << a(iprobe[n],jprobe[n],lo.z,mach) << "\n";               
        }
        ofs.close();

      }
      ParallelDescriptor::Barrier();
    }
}

