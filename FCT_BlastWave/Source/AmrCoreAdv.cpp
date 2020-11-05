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

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        nsubsteps[lev] = 1;
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    bcs.resize(ncomp);

    // Print() << "probtag= " << probtag << "\n";

    //  Set periodic BCs in y-direction
    if(probtag == 1){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::reflect_even);
            bcs[n].setLo(1, BCType::foextrap);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::foextrap);            
        }
        bcs[rou].setLo(0, BCType::reflect_odd);
    }else if(probtag == 2){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::reflect_even);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::reflect_even);
            bcs[n].setHi(1, BCType::reflect_even);            
        }
        bcs[rou].setLo(0, BCType::reflect_odd);
        bcs[rou].setHi(0, BCType::reflect_odd);
        bcs[rov].setLo(1, BCType::reflect_odd);
        bcs[rov].setHi(1, BCType::reflect_odd);
    }else if(probtag == 3 || probtag == 4){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::reflect_even);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::reflect_even);            
        }
        bcs[rou].setLo(0, BCType::reflect_odd);
        bcs[rov].setLo(1, BCType::reflect_odd);
        bcs[rov].setHi(1, BCType::reflect_odd);        
    }else if(probtag == 5){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::reflect_even);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::foextrap);            
        }
        bcs[rou].setLo(0, BCType::reflect_odd);
        bcs[rov].setLo(1, BCType::reflect_odd);        
    }else if(probtag == 6){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::reflect_even);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::ext_dir);            
        }
        bcs[rou].setLo(0, BCType::reflect_odd);
        bcs[rov].setLo(1, BCType::reflect_odd);       
    }else if(probtag == 7){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::ext_dir);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::ext_dir);            
        }
        bcs[rov].setLo(1, BCType::reflect_odd);         
    }
    else if(probtag == 8){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::ext_dir);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::reflect_even);            
        }
        bcs[rov].setLo(1, BCType::reflect_odd); 
        bcs[rov].setHi(1, BCType::reflect_odd);        
    }
    else if(probtag == 9){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::ext_dir);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::ext_dir);            
        }
        bcs[rov].setLo(1, BCType::reflect_odd);        
    }
    else if(probtag == 10 || probtag == 12){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::ext_dir);
            bcs[n].setLo(1, BCType::reflect_even);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::foextrap);            
        }
        bcs[rov].setLo(1, BCType::reflect_odd);        
    }    
    else if(probtag == 11){
        for(int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::ext_dir);
            bcs[n].setLo(1, BCType::foextrap);
            bcs[n].setHi(0, BCType::foextrap);
            bcs[n].setHi(1, BCType::foextrap);            
        }     
    }  
    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);
}

AmrCoreAdv::~AmrCoreAdv ()
{
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
    ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
    pp.query("max_rk", max_rk);
    }

    {
    ParmParse pp("amr"); // Traditionally, these have prefix, amr.

    pp.query("regrid_int", regrid_int);
    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("chk_file", chk_file);
    pp.query("chk_int", chk_int);
    pp.query("restart",restart_chkfile);
    }

    {
    ParmParse pp("adv");

    pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
    pp.query("diff1", diff1);
    }

    {
        ParmParse pp("prob");
        
        pp.get("p2",p2);
        pp.get("p1",p1);
        pp.get("ro2",ro2);
        pp.get("ro1",ro1);
        pp.get("u2",u2);
        pp.get("u1",u1);
        pp.get("v2",v2);
        pp.get("v1",v1);
        pp.query("rad_blast",rad_bw);
        pp.get("probtag",probtag);              
    }

    {
        ParmParse pp("aux");
        
        pp.query("nprobes",nprobes);
        if(nprobes > 0){
            pp.getarr("iprobes", iprobe);
            pp.getarr("jprobes", jprobe);
            if(iprobe.size() != nprobes || jprobe.size() != nprobes){
                amrex::Error("Coorodinate vector of probes does not equal number of probes");
            }
        }               
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    int myproc = ParallelDescriptor::MyProc();
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        ParallelDescriptor::Barrier();
        AverageDown();
        ParallelDescriptor::Barrier();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

        if (plot_int > 0) {
            WritePlotFile();
        }

        if(nprobes > 0){
            AmrCoreAdv::GetProbeDets();
        }
        AmrCoreAdv::WriteProbeFile(0, 0.0, 0);
    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();

        if(nprobes > 0){
            AmrCoreAdv::GetProbeDets();
        }
    }

    // BoxArray ba = phi_new[0].boxArray();
    // const DistributionMapping& dm = phi_new[0].DistributionMap();
    // exact.define(ba,dm,6,0);
    // // err.define(ba,dm,3,0);
    // MultiFab::Copy(exact,phi_new[0],0,0,3,0);
    // MultiFab::Copy(exact,phi_new[0],0,3,3,0);
    // MultiFab::Subtract(exact, exact,0,3,3,0);
    ParallelDescriptor::Barrier();
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                    const DistributionMapping& dm)
{
    int myproc = ParallelDescriptor::MyProc();
    Print(myproc) << "rank= " << myproc << "entering MakeNewLevelFromCoarse()" << "\n";
    const int ncomp = phi_new[lev-1].nComp();
    const int nghost = phi_new[lev-1].nGrow();

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);

    phi_new[lev].FillBoundary();
    phi_new[lev].FillBoundary(geom[lev].periodicity());
    AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);

    phi_old[lev].FillBoundary();
    phi_old[lev].FillBoundary(geom[lev].periodicity());
    AmrCoreAdv::FillDomainBoundary(phi_old[lev],geom[lev],bcs,time);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
             const DistributionMapping& dm)
{
    int myproc = ParallelDescriptor::MyProc();
    Print(myproc) << "rank= " << myproc << "entering RemakeLevel()" << "\n";
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    new_state.FillBoundary();
    new_state.FillBoundary(geom[lev].periodicity());
    AmrCoreAdv::FillDomainBoundary(new_state,geom[lev],bcs,time);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    phi_new[lev].FillBoundary();
    phi_new[lev].FillBoundary(geom[lev].periodicity());
    AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);

    phi_old[lev].FillBoundary();
    phi_old[lev].FillBoundary(geom[lev].periodicity());
    AmrCoreAdv::FillDomainBoundary(phi_old[lev],geom[lev],bcs,time);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    int myproc = ParallelDescriptor::MyProc();
    Print(myproc) << "rank= " << myproc << "entering ClearLevel()" << "\n";
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                      const DistributionMapping& dm)
{
    int myproc = ParallelDescriptor::MyProc();
    // Print(myproc) << "rank= " << myproc << ", entered MakeNewLevelFromScratch" << "\n";
    // Real ro2, ro1, u2, u1, v2, v1, p, xw, yw;

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    int nc = phi_new[lev].nComp();

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    Real pminfrac = 0.01, rominfrac = 0.01;
    {
        ParmParse pp("prob");
        pp.query("pminfrac",pminfrac);
        pp.query("pminfrac",rominfrac);
    }

    if (lev > 0 && do_reflux) {
    flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real* prob_hi = geom[lev].ProbHi();
    Real cur_time = t_new[lev];

    MultiFab& state = phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

        FArrayBox& mfab = state[mfi];
        Array4<Real> const& a = mfab.array();

        AmrCoreAdv::initBlastWave(lev,box,a,geom[lev]);
    }
    // Print(myproc) << "rank= " << myproc << ", does phi_new have NaNs (before BC)? " << phi_new[lev].contains_nan() << "\n";
    ParallelDescriptor::Barrier();
    if (lev == 0) {
        FillPatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        phi_new[lev].FillBoundary();
        phi_new[lev].FillBoundary(geom[lev].periodicity());
        AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);
    }
    else {
        FillPatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        phi_new[lev].FillBoundary();
        phi_new[lev].FillBoundary(geom[lev].periodicity());
        AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);
    }
    pmin = pminfrac*phi_new[lev].min(pre);
    romin = rominfrac*phi_new[lev].min(ro);
    // Print(myproc) << "rank= " << myproc << ", does phi_new have NaNs (after BC)? " << phi_new[lev].contains_nan() << "\n";
    Print() << "min(ro)= " << phi_new[lev].min(ro,4) << ", min(pre)= " << phi_new[lev].min(pre,4)
    << ", min(roE)= " << phi_new[lev].min(roE,4) << ", min(mach)= " << phi_new[lev].min(mach,4) << "\n";
    ParallelDescriptor::Barrier();
    // Print(myproc) << "rank= " << myproc << ", reached end of MakeNewLevelFromScratch()" << "\n";
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    int myproc = ParallelDescriptor::MyProc();
    // Print(myproc) << "rank= " << myproc << ", lev= " << lev << ", entering ErrorEst()" << "\n";
    ParallelDescriptor::Barrier();
    static bool first = true;
    // static Vector<Real> tagfrac;
    static Real tagfrac = 0.1;

    // only do this during the first call to ErrorEst
    if (first)
    {
    first = false;
       // read in an array of "phierr", which is the tagging threshold
        // in this example, we tag values of "phi" which are greater than phierr
        // for that particular level
        // in subroutine state_error, you could use more elaborate tagging, such
        // as more advanced logical expressions, or gradients, etc.
    ParmParse pp("adv");
    pp.query("tagfrac",tagfrac);
    }

    // if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = phi_new[lev];
    int comp_lo = ro;
    int comp_hi = pre;
    Real maxgradpx = 0.0, maxgradpy = 0.0, gradpxtemp, gradpytemp;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

        for(MFIter mfi(state, true); mfi.isValid(); ++mfi){
            const Box& tmpbox  = mfi.validbox();
            Array4<Real const> const& a = state[mfi].const_array();

            gradpxtemp = get_gradp(lev, tmpbox, a, 0);
            if(fabs(gradpxtemp) > fabs(maxgradpx)){
                maxgradpx = gradpxtemp;
            }

            gradpytemp = get_gradp(lev, tmpbox, a, 1);
            if(fabs(gradpytemp) > fabs(maxgradpy)){
                maxgradpy = gradpytemp;
            }
        }
    // Get global maximum of pressure gradient (use this as criterion for refinement)
    ParallelDescriptor::ReduceRealMax(maxgradpx);
    ParallelDescriptor::ReduceRealMax(maxgradpy);

    for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
    {
        const Box& validbox  = mfi.validbox();

        TagBox&     tagfab  = tags[mfi];

        // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
        // So we are going to get a temporary integer array.
            // set itags initially to 'untagged' everywhere
            // we define itags over the tilebox region
        tagfab.get_itags(itags, validbox);

            // data pointer and index space
        int*        tptr    = itags.dataPtr();
        const int*  tlo     = validbox.loVect();
        const int*  thi     = validbox.hiVect();

        //Function to calculate the maximum pressure gradient and establish refinement criterion
            // tag cells for refinement
        state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
            BL_TO_FORTRAN_3D(state[mfi]),
            &tagval, &clearval,
            AMREX_ARLIM_3D(validbox.loVect()), AMREX_ARLIM_3D(validbox.hiVect()),
            AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &comp_lo, &comp_hi, &maxgradpx, &maxgradpy,
            &tagfrac,&rad_bw,&lev);
        //
        // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
        //
        tagfab.tags_and_untags(itags, validbox);
    }
    }
    // Print(myproc) << "rank= " << myproc << ", reached end of ErrorEst(b Barrier)" << "\n";
    ParallelDescriptor::Barrier();
    // Print(myproc) << "rank= " << myproc << ", reached end of ErrorEst()" << "\n";
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    int myproc = ParallelDescriptor::MyProc();
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
    amrex::average_down(phi_new[lev+1], phi_new[lev],
                            geom[lev+1], geom[lev],
                            0, phi_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(phi_new[crse_lev+1], phi_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, phi_new[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
    Vector<MultiFab*> smf;
    Vector<Real> stime;
    GetData(0, time, smf, stime);

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> physbc(geom[lev], bcs, bfunc);
    amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                    geom[lev], physbc, 0);
    }
    else
    {
    Vector<MultiFab*> cmf, fmf;
    Vector<Real> ctime, ftime;
    GetData(lev-1, time, cmf, ctime);
    GetData(lev  , time, fmf, ftime);

        BndryFuncArray bfunc(phifill);
        PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
        PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

    Interpolater* mapper = &cell_cons_interp;

    amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                  0, icomp, ncomp, geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

 //    if (cmf.size() != 1) {
    // amrex::Abort("FillCoarsePatch: how did this happen?");
 //    }

    BndryFuncArray bfunc(phifill);
    PhysBCFunct<BndryFuncArray> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<BndryFuncArray> fphysbc(geom[lev  ],bcs,bfunc);

    Interpolater* mapper = &cell_cons_interp;

    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                 cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                 mapper, bcs, 0);
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
    data.push_back(&phi_new[lev]);
    datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
    data.push_back(&phi_old[lev]);
    datatime.push_back(t_old[lev]);
    }
    else
    {
    data.push_back(&phi_old[lev]);
    data.push_back(&phi_new[lev]);
    datatime.push_back(t_old[lev]);
    datatime.push_back(t_new[lev]);
    }
}

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    int myproc = ParallelDescriptor::MyProc();
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration, step);

	cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;
        amrex::Print() << "Max mach= " << phi_new[0].norm0(mach) << "\n";
        // Print(myproc) << "rank= " << myproc << ", does phi_new have NaNs? " << phi_new[0].contains_nan() << "\n";

	// sync up time
	for (lev = 0; lev <= finest_level; ++lev) {
	    t_new[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();
        // WriteErrFile();
	}
    AmrCoreAdv::WriteProbeFile(0, cur_time, step);

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

	if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (chk_int > 0) {
       WriteCheckpointFile();
       // WriteErrFile();
    }

    if(plot_int > 0){
        WritePlotFile();
    }
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
    int myproc = ParallelDescriptor::MyProc();
    Vector<Real> dt_tmp(finest_level+1);

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        dt_tmp[lev] = EstTimeStep(lev, true);
    }
    ParallelDescriptor::ReduceRealMin(&dt_tmp[0], dt_tmp.size());

    constexpr Real change_max = 1.1;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
    dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
    n_factor *= nsubsteps[lev];
    dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const Real eps = 1.e-3*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
    dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
    dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
    
    for (int lev = 0; lev <= finest_level; ++lev) {
    dt[lev] = dt[finest_level];
    }
    // Print(myproc) << "rank= " << myproc << ", dt = " << dt[0] <<"\n";
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, bool local)
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    int myproc = ParallelDescriptor::MyProc();

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    MultiFab& S_new = phi_new[lev];
    int nc = phi_new[lev].nComp();
    int ng = nghost;
    Real dt_calc = 0.0;
    Real umax = 0.0, vmax = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
       FArrayBox uface[BL_SPACEDIM];

       for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
       {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
           const Box& bx = mfi.nodaltilebox(i);
           uface[i].resize(bx,1);
           }

            const Box& box = mfi.tilebox();
            const int* lo  = box.loVect();
            const int* hi  = box.hiVect();

           get_face_velocity_dt(&lev, &cur_time,
                     AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                     BL_TO_FORTRAN(uface[1]),
                     BL_TO_FORTRAN(uface[2])),
                     dx, prob_lo, &umax, &vmax, &nc,
                     BL_TO_FORTRAN_3D(S_new[mfi]), 
                     AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi));
            // Print(myproc) << "rank= " << myproc << ", umax= " << umax << ", vmax= " << vmax
            // << "dx= "<< geom[lev].CellSize(0) << ", dy= " << geom[lev].CellSize(1) << "\n";
        dt_calc = std::min(geom[lev].CellSize(0)/umax, geom[lev].CellSize(1)/vmax);
        dt_est = std::min(dt_est, dt_calc);
       }
    }

    if (!local) {
    ParallelDescriptor::ReduceRealMin(dt_est);
    }

    dt_est *= cfl;

    return dt_est;
}
// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStep (int lev, Real time, int iteration, int step)
{
    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        Print() << "lev= " << lev << ", istep= " << istep[lev] << 
        ", last_regrid_step= " << last_regrid_step[lev] << "\n";

        if (lev < max_level && istep[lev] > last_regrid_step[lev])
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
		int old_finest = finest_level;
		regrid(lev, time);

                // mark that we have regridded this level already
		for (int k = lev; k <= finest_level; ++k) {
		    last_regrid_step[k] = istep[k];
		}

                // if there are newly created levels, set the time step
		for (int k = old_finest+1; k <= finest_level; ++k) {
		    dt[k] = dt[k-1] / MaxRefRatio(k-1);
		}
        for (int k = 0; k <= finest_level; ++k){
            dt[k] = dt[finest_level];
        }
	    }
	}
    }

    if (Verbose()) {
	amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
	amrex::Print() << "ADVANCE with time = " << t_new[lev]
                       << " dt = " << dt[lev] << std::endl;
    }

    // advance a single level for a single time step, updates flux registers
    Advance(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
	amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
	   for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	   {
	       timeStep(lev+1, time+(i-1)*dt[lev+1], i, step);
	   }

        ParallelDescriptor::Barrier();
        phi_new[lev].FillBoundary();
        phi_new[lev].FillBoundary(geom[lev].periodicity());
        AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);

	   if (do_reflux)
	   {
            // update lev based on coarse-fine flux mismatch
	       flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);

            ParallelDescriptor::Barrier();

            for (MFIter mfi(phi_new[lev], true); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();
                FArrayBox& fab = phi_new[lev][mfi];
                Array4<Real> const& a = fab.array();
                CalcAuxillary (lev, box, a, geom[lev]);
            }//for(MFIter)

            ParallelDescriptor::Barrier();
            phi_new[lev].FillBoundary();
            AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);
            phi_new[lev].FillBoundary(geom[lev].periodicity());
	   }
        ParallelDescriptor::Barrier();
        AverageDownTo(lev); // average lev+1 down to lev
        for (MFIter mfi(phi_new[lev], true); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            FArrayBox& fab = phi_new[lev][mfi];
            Array4<Real> const& a = fab.array();
            CalcAuxillary (lev, box, a, geom[lev]);
        }//for(MFIter)
        if(lev == 0){
            phi_new[lev].FillBoundary();
            AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);
            phi_new[lev].FillBoundary(geom[lev].periodicity());
        }else{
            phi_new[lev].FillBoundary();
            phi_new[lev].FillBoundary(geom[lev].periodicity());
            AmrCoreAdv::FillDomainBoundary(phi_new[lev],geom[lev],bcs,time);
        }

    }

    ParallelDescriptor::Barrier();
}

// advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle)
{
    int myproc = ParallelDescriptor::MyProc();
    // Print(myproc) << "rank= " << myproc << "entered Advance()" << "\n";
    constexpr int num_grow = 4;
    int ngrow = num_grow;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    MultiFab& S_new = phi_new[lev];
    MultiFab& S_old = phi_old[lev];
    MultiFab& exact_sol = exact;
    const Geometry& geom1 = geom[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const Real* dx = geom1.CellSize();
    const Real* prob_lo = geom1.ProbLo();
    int nc = S_new.nComp();
    int rk_max = max_rk;

    MultiFab fluxes[BL_SPACEDIM], flux1[BL_SPACEDIM], uface1[BL_SPACEDIM];
    if (do_reflux)
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    BoxArray ba = grids[lev];
	    ba.surroundingNodes(i);
	    fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
        flux1[i].define(ba, dmap[lev], S_new.nComp(), num_grow);
        uface1[i].define(ba, dmap[lev], S_new.nComp(), num_grow);
	}
    }

    // State with ghost cells (partially convected variables in x, y directions)
    MultiFab Sconvx(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    MultiFab Sconvy(grids[lev], dmap[lev], S_new.nComp(), num_grow);

    // Print(myproc) << "rank= " << myproc << ", does S_new contains nan before FillPatch? " << S_new.contains_nan() 
    // << ", " << Sconvx.contains_nan() << ", " << Sconvy.contains_nan() << ", max p = " 
    // << S_new.norm0(pre) << "\n";

    ParallelDescriptor::Barrier();
    if(lev == 0){
        FillPatch(lev, time, Sconvx, 0, Sconvx.nComp());
        FillPatch(lev, time, Sconvy, 0, Sconvy.nComp());
        FillPatch(lev, time, S_new, 0, S_new.nComp());

        S_new.FillBoundary();
        S_new.FillBoundary(geom1.periodicity());
        AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
        Sconvx.FillBoundary();
        Sconvx.FillBoundary(geom1.periodicity());
        AmrCoreAdv::FillDomainBoundary(Sconvx,geom1,bcs,time);
        Sconvy.FillBoundary();
        AmrCoreAdv::FillDomainBoundary(Sconvy,geom1,bcs,time);
        Sconvy.FillBoundary(geom1.periodicity());
    }else{
        FillCoarsePatch(lev, time, Sconvx, 0, Sconvx.nComp());
        FillCoarsePatch(lev, time, Sconvy, 0, Sconvy.nComp());
        FillCoarsePatch(lev, time, S_new, 0, S_new.nComp());

        S_new.FillBoundary();
        S_new.FillBoundary(geom1.periodicity());
        AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
        Sconvx.FillBoundary();
        Sconvx.FillBoundary(geom1.periodicity());
        AmrCoreAdv::FillDomainBoundary(Sconvx,geom1,bcs,time);
        Sconvy.FillBoundary();
        AmrCoreAdv::FillDomainBoundary(Sconvy,geom1,bcs,time);
        Sconvy.FillBoundary(geom1.periodicity());       
    }
    // Print(myproc) << "rank= " << myproc << ", does S_new contains nan after FillPatch? " << S_new.contains_nan() 
    // << ", " << Sconvx.contains_nan() << ", " << Sconvy.contains_nan() << ", max p = " 
    // << S_new.norm0(pre) << "\n";

    ParallelDescriptor::Barrier();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	   FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

       for(int rk = 1 ; rk <= rk_max ; ++rk){
        for(int fct_step = 1 ; fct_step <= 2 ; ++fct_step){
            for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                FArrayBox& uvel = uface1[0][mfi];
                FArrayBox& vvel = uface1[1][mfi];
                FArrayBox& flxx = flux1[0][mfi];
                FArrayBox& flxy = flux1[1][mfi];

                FArrayBox& statex   = Sconvx[mfi];
                FArrayBox& statey   = Sconvy[mfi];
                FArrayBox& stateout = S_new[mfi];
                FArrayBox& stateold = S_old[mfi];

                const int* lo  = bx.loVect();
                const int* hi  = bx.hiVect();

                // Allocate fabs for fluxes and Godunov velocities.
                if(fct_step == 1 && rk == 1){
                for (int i = 0; i < BL_SPACEDIM ; i++) {
                        const Box& bxtmp = amrex::surroundingNodes(bx,i);
                        flux[i].resize(amrex::grow(bxtmp,ngrow),S_new.nComp()); // need to modify this to reduce memory
                        uface[i].resize(amrex::grow(bxtmp,ngrow),1);
                    }
                }

                // compute velocities on faces (prescribed function of space and time)
                if(fct_step == 1){
                        get_face_velocity(&lev, &ctr_time, 
                            AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
                            AMREX_D_DECL(BL_TO_FORTRAN(uvel),
                            BL_TO_FORTRAN(vvel),
                            BL_TO_FORTRAN(uface[2])),
                            dx, prob_lo, 
                            BL_TO_FORTRAN_3D(stateold), 
                            &nc);
                }//if(fct_step==1)

                //advection step (solve the conservation equations & compute new state)
                advect(&lev, &time, &rk, &rk_max, &fct_step, &nc, 
                    AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
                    BL_TO_FORTRAN_3D(stateold),
                    BL_TO_FORTRAN_3D(statex),
                    BL_TO_FORTRAN_3D(statey),
                    BL_TO_FORTRAN_3D(stateout),
                    AMREX_D_DECL(BL_TO_FORTRAN_3D(uvel),
                    BL_TO_FORTRAN_3D(vvel),
                    BL_TO_FORTRAN_3D(uface[2])),
                    AMREX_D_DECL(BL_TO_FORTRAN_3D(flxx),
                    BL_TO_FORTRAN_3D(flxy),
                    BL_TO_FORTRAN_3D(flux[2])),
                    dx, &dt_lev, &diff1, &pmin, &romin);
                if(do_reflux && rk == rk_max && fct_step == 2){
                    // for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi){
                        for(int i = 0; i < BL_SPACEDIM ; i++){
                            fluxes[i][mfi].copy<RunOn::Host>(flux1[i][mfi],mfi.nodaltilebox(i));
                        }
                    // }
                }
            }//for MFIter(mfi)
            ParallelDescriptor::Barrier();
            // if(S_new.contains_nan()){
                // Print() << "rank= " << myproc << ", does S_new contains nan before BC? " << S_new.contains_nan() 
                // << ", " << Sconvx.contains_nan() << ", " << Sconvy.contains_nan() << "\n";
            // }
            if (lev == 0){
                    S_new.FillBoundary();
                    AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
                    S_new.FillBoundary(geom1.periodicity());
                    Sconvx.FillBoundary();
                    AmrCoreAdv::FillDomainBoundary(Sconvx,geom1,bcs,time);
                    Sconvx.FillBoundary(geom1.periodicity());
                    Sconvy.FillBoundary();
                    AmrCoreAdv::FillDomainBoundary(Sconvy,geom1,bcs,time);
                    Sconvy.FillBoundary(geom1.periodicity());
            }
            else{
                    S_new.FillBoundary();
                    S_new.FillBoundary(geom1.periodicity());
                    AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
                    Sconvx.FillBoundary();
                    Sconvx.FillBoundary(geom1.periodicity());
                    AmrCoreAdv::FillDomainBoundary(Sconvx,geom1,bcs,time);
                    Sconvy.FillBoundary();
                    Sconvy.FillBoundary(geom1.periodicity());
                    AmrCoreAdv::FillDomainBoundary(Sconvy,geom1,bcs,time);
            }
            // Check if density/pressure/Mach number become negative
            if(S_new.min(ro,ngrow) < 0.0 || S_new.min(pre,ngrow) < 0.0 || S_new.min(mach,ngrow) < 0.0 || S_new.min(roE,ngrow) < 0.0){
                Print() << "End of FCT step= " << fct_step << ", RK= " << rk << "\n";
                Print() << "min ro= " << S_new.min(ro,ngrow) 
                        << ", min pre= " << S_new.min(pre,ngrow) << ", min roE= " << S_new.min(roE,ngrow)
                        << ", min mach= " << S_new.min(mach,ngrow) <<"\n";
                WritePlotFile();
                amrex::Error("Pressure/density is negative, aborting...");
            }
            if(S_new.contains_nan() || Sconvx.contains_nan() || Sconvy.contains_nan()){
                Print() << "End of FCT step= " << fct_step << ", RK= " << rk << "\n";
                Print() << "S_new contains nan after BC? " << S_new.contains_nan() 
                << ", " << Sconvx.contains_nan() << ", " << Sconvy.contains_nan() << "\n";
                amrex::Error("NaN value found in conserved variables, aborting...");
            }
        ParallelDescriptor::Barrier();
        }//for(fct_step)
    }//for(rk)

    if (lev == 0){
        S_new.FillBoundary();
        AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
        S_new.FillBoundary(geom1.periodicity());
    }
    else{
        S_new.FillBoundary();
        S_new.FillBoundary(geom1.periodicity());
        AmrCoreAdv::FillDomainBoundary(S_new,geom1,bcs,time);
    }    

    }
    ParallelDescriptor::Barrier();
    // increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) {
	if (flux_reg[lev+1]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev+1/lev flux register (index lev+1)
	        flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(), -1.0);
	    }
	}
	if (flux_reg[lev]) {
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
	        // update the lev/lev-1 flux register (index lev)
		flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(), 1.0);
	    }
	}
    }
    ParallelDescriptor::Barrier();
}
