#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

class AmrCoreAdv
    : public amrex::AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrCoreAdv ();
    virtual ~AmrCoreAdv();

    // advance solution to final time
    void Evolve ();

    // initializes multilevel data
    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from phi_old and/or phi_new into another multifab
    void GetData (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // advance a level by dt
    // includes a recursive call for finer levels
    void timeStep (int lev, amrex::Real time, int iteration, int step);

    // advance a single level for a single time step, updates flux registers
    void Advance (int lev, amrex::Real time, amrex::Real dt_lev, int iteration, int ncycle);

    // a wrapper for EstTimeStep(0
    void ComputeDt ();

    // compute dt from CFL considerations
    Real EstTimeStep (int lev, bool local=false);

    // function to initialize the isentropic vortex
    void initBlastWave (int &lev, Box const& bx, Array4<Real> const& a, const Geometry& geom);

    void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc, Real cur_time);

    // Need to declare as static, else pass function pointer in AmrCoreAdv::FillDomainBoundary
    static void outletBC_partialwall (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp);

    struct dummy_gpu_fill_extdir
    {
        AMREX_GPU_DEVICE
        void operator() (const IntVect& iv, Array4<Real> const& dest,
                         const int dcomp, const int numcomp,
                         GeometryData const& geom, const Real time,
                         const BCRec* bcr, const int bcomp,
                         const int orig_comp) const
            {
                // do something for external Dirichlet (BCType::ext_dir) if there are
            }
    };

    void CalcAuxillary (int lev, Box const& bx, Array4<Real> const& a, const Geometry& geom);

    // Functions to write data from probes to file
    void GetProbeDets();

    Real get_gradp(int lev, Box const& validbox, Array4<Real const> const& a, int dir);
    void WriteProbeFile (int lev, Real cur_time, int stepnum);

    // get plotfile name
    std::string PlotFileName (int lev) const;
    std::string ErrFileName (int lev) const;

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> PlotFileMF () const;
    amrex::Vector<const amrex::MultiFab*> ErrFileMF () const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;
    amrex::Vector<std::string> ErrFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;
    void WriteErrFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    // utility to skip to next line in Header
    static void GotoNextLine (std::istream& is);
    ////////////////
    // private data members

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".
    amrex::Vector<amrex::MultiFab> phi_new;
    amrex::Vector<amrex::MultiFab> phi_old;
    amrex::MultiFab exact;
    // amrex::MultiFab err;

    // BCRec is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    amrex::Vector<BCRec> bcs;  // 1-component

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    amrex::Vector<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    // if >= 0 we restart from a checkpoint
    std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    amrex::Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    int do_reflux = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    std::string err_file {"err"};
    int plot_int = -1;

    // checkpoint prefix and frequency
    std::string chk_file {"chk"};
    int chk_int = -1;

    int ncomp  = 6;
    const int ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, mach = 5;
    int nghost = 4;
    int max_rk = 2;
    int probtag = -1;

    // Probes for measuring quantities varying with time
    Vector<int> iprobe, jprobe, proberank;
    int nprobes = 0;
    Real xcm = 0.0, ycm = 0.0, rad_bw = 0.01;

    Real p2, p1, ro2, ro1, u2, u1, v2, v1;

    Real diff1 = 1.0, pmin, romin;

    // some parameters
    const Real pi = 3.141592653589793238;
    const Real gamma = 1.4, c1 = 1/(gamma-1);
};

#endif
