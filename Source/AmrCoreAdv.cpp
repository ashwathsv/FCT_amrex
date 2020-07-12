
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
	nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    bcs.resize(ncomp);

    if (domdir == 1) {
        //  Set periodic BCs in y-direction
        for (int n = 0; n < ncomp; ++n){
            bcs[n].setLo(1, BCType::int_dir);
            bcs[n].setHi(1, BCType::int_dir);
        }
        // Density BCs
        bcs[ro].setLo(0,BCType::reflect_even); 
        bcs[ro].setHi(0,BCType::reflect_even);
        // x-momentum BCs
        bcs[rou].setLo(0,BCType::reflect_odd); 
        bcs[rou].setHi(0,BCType::reflect_odd);
        // y-momentum BCs
        bcs[rov].setLo(0,BCType::reflect_even);
        bcs[rov].setHi(0,BCType::reflect_even);
        // Energy BCs
        bcs[roE].setLo(0,BCType::reflect_even);
        bcs[roE].setHi(0,BCType::reflect_even);
        //Pressure BCs
        bcs[pre].setLo(0,BCType::reflect_even);
        bcs[pre].setHi(0,BCType::reflect_even);
    }
    else {
        // Set periodic BCs in x
        for (int n = 0; n < ncomp; ++n){
            bcs[n].setLo(0, BCType::int_dir);
            bcs[n].setHi(0, BCType::int_dir);
        }
        // Density BCs
        bcs[ro].setLo(1,BCType::reflect_even); 
        bcs[ro].setHi(1,BCType::reflect_even);
        // x-momentum BCs
        bcs[rou].setLo(1,BCType::reflect_odd); 
        bcs[rou].setHi(1,BCType::reflect_odd);
        // y-momentum BCs
        bcs[rov].setLo(1,BCType::reflect_even);
        bcs[rov].setHi(1,BCType::reflect_even);
        // Energy BCs
        bcs[roE].setLo(1,BCType::reflect_even);
        bcs[roE].setHi(1,BCType::reflect_even);
        //Pressure BCs
        bcs[pre].setLo(1,BCType::reflect_even);
        bcs[pre].setHi(1,BCType::reflect_even);
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

// advance solution to final time
void
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

	ComputeDt();

	int lev = 0;
	int iteration = 1;
	timeStep(lev, cur_time, iteration);

	cur_time += dt[0];

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0]  << std::endl;

	// sync up time
	for (lev = 0; lev <= finest_level; ++lev) {
	    t_new[lev] = cur_time;
	}

	if (plot_int > 0 && (step+1) % plot_int == 0) {
	    last_plot_file_step = step+1;
	    WritePlotFile();

        Real time = t_new[0];
        int stepnum = istep[0];
        for (int i = 0; i <= finest_level; ++i) {
            AmrCoreAdv::writetxtfile(i, time, stepnum);
        }
	}

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

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
	   WritePlotFile();
        Real time = t_new[0];
        int stepnum = istep[0];
        for (int i = 0; i <= finest_level; ++i) {
            AmrCoreAdv::writetxtfile(i, time, stepnum);
        }
    }
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
        // writetxtfile();
    }
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				    const DistributionMapping& dm)
{
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
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

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
    // int ncomp = 5;
    // const int nghost = 4;
    // int domdir = 1;
    // {    
    //     ParmParse pp;
    //     pp.query("ncomp",ncomp);
    //     pp.query("domdir",domdir);
    // }
    Real ro2ro1, v2v1, p2p1;
    {
        ParmParse pp("prob");
        pp.query("ro2ro1",ro2ro1);
        pp.query("p2p1",p2p1);
        if ( domdir == 1 ) {
            pp.query("v2v1",v2v1);
        }
        if ( domdir == 2 ) {
            pp.query("v2v1",v2v1);
        }        
    }

    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    int nc = phi_new[lev].nComp();

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

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
        // FArrayBox& fab = state[mfi];
        // Array4<Real> const& a = fab.array();

	initdata(&lev, &cur_time, &nc, &domdir, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
		 BL_TO_FORTRAN_3D(state[mfi]), AMREX_ZFILL(dx),
		 AMREX_ZFILL(prob_lo), AMREX_ZFILL(prob_hi), &ro2ro1, &v2v1, &p2p1);
    //  Function to print results
        // FArrayBox& fab = state[mfi];
        // Array4<Real> const& a = fab.array();
        // printout(box, a);

    }
    // if (not geom.isAllPeriodic()) {
    // GpuBndryFuncFab<MyExtBCFill> bf(MyExtBCFill{});
    // PhysBCFunct<GpuBndryFuncFab<MyExtBCFill> > physbcf(geom, bc, bf);
    // physbcf(mf, 0, mf.nComp(), mf.nGrowVector(), time, 0);
// }
    // phi_new[lev].FillBoundary();
    if (lev == 0) {
        FillPatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        amrex::FillDomainBoundary(phi_new[lev],geom[lev],bcs);
        phi_new[lev].FillBoundary(geom[lev].periodicity());
        // amrex::FillDomainBoundary (phi_new[lev], geom[lev], bcs);
    }
    else {
        FillPatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        phi_new[lev].FillBoundary();
    }

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        // FArrayBox& fab = state[mfi];
        // Array4<Real> const& a = fab.array();
    //  Function to print results
        FArrayBox& fab = state[mfi];
        Array4<Real> const& a = fab.array();
        // printout(box, a);

    }

    MultiFab& phi = phi_new[lev];
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        // FArrayBox& fab = state[mfi];
        // Array4<Real> const& a = fab.array();
        //  Function to print results
        FArrayBox& mfab = phi[mfi];
        Array4<Real> const& arr = mfab.array();
        WriteTxtFileInit(lev, box, arr, istep[lev], geom[lev]);

    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    static bool first = true;
    static Vector<Real> phierr;

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
	int n = pp.countval("phierr");
	if (n > 0) {
	    pp.getarr("phierr", phierr, 0, n);
	}
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = phi_new[lev];
    int nc = phi_new[lev].nComp();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

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

            // tag cells for refinement
	    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
			BL_TO_FORTRAN_3D(state[mfi]),
			&tagval, &clearval,
			AMREX_ARLIM_3D(validbox.loVect()), AMREX_ARLIM_3D(validbox.hiVect()),
			AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &phierr[lev], &nc, &domdir);
	    //
	    // Now update the tags in the TagBox in the tilebox region
            // to be equal to itags
	    //
	    tagfab.tags_and_untags(itags, validbox);
	}
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    pp.query("domdir", domdir);
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
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
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

    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

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


// advance a level by dt
// includes a recursive call for finer levels
void
AmrCoreAdv::timeStep (int lev, Real time, int iteration)
{
    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
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
	    timeStep(lev+1, time+(i-1)*dt[lev+1], i);
	}

	if (do_reflux)
	{
            // update lev based on coarse-fine flux mismatch
	    flux_reg[lev+1]->Reflux(phi_new[lev], 1.0, 0, 0, phi_new[lev].nComp(), geom[lev]);
	}

	AverageDownTo(lev); // average lev+1 down to lev
    }

}

// advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::Advance (int lev, Real time, Real dt_lev, int iteration, int ncycle)
{
    constexpr int num_grow = 4;
    int ngrow = num_grow;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    MultiFab& S_new = phi_new[lev];
    MultiFab& S_old = phi_old[lev];
    const Geometry& geom1 = geom[lev];

    const Real old_time = t_old[lev];
    const Real new_time = t_new[lev];
    const Real ctr_time = 0.5*(old_time+new_time);

    const Real* dx = geom1.CellSize();
    // Print() << "dx = " << geom1.CellSize(0) << ", dy = " << geom1.CellSize(1) << "\n";
    const Real* prob_lo = geom1.ProbLo();
    int ddir = domdir;
    int nc = S_new.nComp();
    int rk_max = max_rk;

    MultiFab fluxes[BL_SPACEDIM];
    if (do_reflux)
    {
	for (int i = 0; i < BL_SPACEDIM; ++i)
	{
	    BoxArray ba = grids[lev];
	    ba.surroundingNodes(i);
	    fluxes[i].define(ba, dmap[lev], S_new.nComp(), 0);
	}
    }

    // State with ghost cells
    MultiFab Sconvx(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    MultiFab Sconvy(grids[lev], dmap[lev], S_new.nComp(), num_grow);

    if(lev == 0){
        FillPatch(lev, time, Sconvx, 0, Sconvx.nComp());
        FillPatch(lev, time, Sconvy, 0, Sconvy.nComp());
        FillPatch(lev, time, S_new, 0, S_new.nComp());
    }else{
        FillCoarsePatch(lev, time, Sconvx, 0, Sconvx.nComp());
        FillCoarsePatch(lev, time, Sconvy, 0, Sconvy.nComp());
        FillCoarsePatch(lev, time, S_new, 0, S_new.nComp());       
    }

    // MultiFab S_old(grids[lev], dmap[lev], S_new.nComp(), num_grow);
    // FillPatch(lev, time, S_old, 0, S_old.nComp());
    

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	   FArrayBox flux[BL_SPACEDIM], uface[BL_SPACEDIM];

	   for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	   {
	        const Box& bx = mfi.tilebox();

	        FArrayBox& statex   = Sconvx[mfi];
            FArrayBox& statey   = Sconvy[mfi];
	        FArrayBox& stateout = S_new[mfi];
            FArrayBox& stateold = S_old[mfi];

            const int* lo  = bx.loVect();
            const int* hi  = bx.hiVect();

	       // Allocate fabs for fluxes and Godunov velocities.
	       for (int i = 0; i < BL_SPACEDIM ; i++) {
		      const Box& bxtmp = amrex::surroundingNodes(bx,i);
		      flux[i].resize(amrex::grow(bxtmp,ngrow),S_new.nComp()); // need to modify this to reduce memory
		      uface[i].resize(amrex::grow(bxtmp,ngrow),1);

                // Print() << "i= " << i << "flux_lim: lo->" << flux[i].smallEnd() << ": hi-> " 
                //         << flux[i].bigEnd() << "\n";
                // Print() << "i= " << i << "uface_lim: lo->" << uface[i].smallEnd() << ": hi-> " 
                //         << uface[i].bigEnd() << "\n";
	        }
            // Print() << "lo= " << *bx.loVect(3) << "\n";
            // compute velocities on faces (prescribed function of space and time)


            // compute new state (stateout) and fluxes.
            for(int rk = 1 ; rk <= rk_max ; ++rk){
                if(rk == 1){
                            get_face_velocity(&lev, &ctr_time, 
                            AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
                            AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                            BL_TO_FORTRAN(uface[1]),
                            BL_TO_FORTRAN(uface[2])),
                            dx, prob_lo, 
                            BL_TO_FORTRAN_3D(stateold), 
                            &nc, &ddir);                
                }else{
                            get_face_velocity(&lev, &ctr_time, 
                            AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
                            AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
                            BL_TO_FORTRAN(uface[1]),
                            BL_TO_FORTRAN(uface[2])),
                            dx, prob_lo, 
                            BL_TO_FORTRAN_3D(stateout), 
                            &nc, &ddir);                 
                }//if(rk==1)
            for(int fct_step = 1 ; fct_step <= 2 ; ++fct_step){
                // Array4<Real> const& a = stateout.array();
                Array4<Real> const& b = stateold.array();

                advect(&lev, &time, &rk, &rk_max, &fct_step, &ddir, &nc, 
                    AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
                    BL_TO_FORTRAN_3D(stateold),
		            BL_TO_FORTRAN_3D(statex),
		            BL_TO_FORTRAN_3D(statey),
                    BL_TO_FORTRAN_3D(stateout),
		            AMREX_D_DECL(BL_TO_FORTRAN_3D(uface[0]),
			        BL_TO_FORTRAN_3D(uface[1]),
			        BL_TO_FORTRAN_3D(uface[2])),
		            AMREX_D_DECL(BL_TO_FORTRAN_3D(flux[0]),
			        BL_TO_FORTRAN_3D(flux[1]),
			        BL_TO_FORTRAN_3D(flux[2])),
		            dx, &dt_lev);
                    // printout(bx, a);
                    // printout(bx,b); 
                    Array4<Real> const& a = stateout.array();
                
                    if (lev == 0){
                        // printout(bx, a);
                        amrex::FillDomainBoundary(S_new,geom1,bcs);
                        S_new.FillBoundary(geom1.periodicity());
                        amrex::FillDomainBoundary(Sconvx,geom1,bcs);
                        Sconvx.FillBoundary(geom1.periodicity());
                        amrex::FillDomainBoundary(Sconvy,geom1,bcs);
                        Sconvy.FillBoundary(geom1.periodicity());
                        // Print() << "Boundaries filled ------------------------ \n";
                        // printout(bx, a);
                        // FillPatch(lev, time, S_new, 0, S_new.nComp());
                    }
                    else{
                        S_new.FillBoundary();
                        Sconvx.FillBoundary();
                        Sconvy.FillBoundary();
                    }
                    Print() << "End of FCT step= " << fct_step << ", RK= " << rk << "\n";
                    // Print() << "-------------------------------" << "\n";
                    // Array4<Real> const& b = stateout.array();
                    // printout(bx, b);
                    // else{
                    //     FillCoarsePatch(lev, time, S_new, 0, S_new.nComp());
                    // }
                }//for(fct_step)           
            }//for(rk)

        if (lev == 0){
                    amrex::FillDomainBoundary(S_new,geom1,bcs);
                    S_new.FillBoundary(geom1.periodicity());
        }else{
            S_new.FillBoundary();
        }
        // if (lev == 0){
        //         FillPatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        // }


        // else{
        //         FillCoarsePatch(lev, time, phi_new[lev], 0, phi_new[lev].nComp());
        // }
	    if (do_reflux) {
		  for (int i = 0; i < BL_SPACEDIM ; i++) {
		      fluxes[i][mfi].copy<RunOn::Host>(flux[i],mfi.nodaltilebox(i));
                // Print() << "i= " << i << ", reached end of advance()" << "\n";
		  }
	    }
	   }//for(MFIter mfi)
    }

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
}

// a wrapper for EstTimeStep
void
AmrCoreAdv::ComputeDt ()
{
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
    Print() << "dt = " << dt[0] <<"\n";
}

// compute dt from CFL considerations
Real
AmrCoreAdv::EstTimeStep (int lev, bool local)
{
    BL_PROFILE("AmrCoreAdv::EstTimeStep()");

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    const Real cur_time = t_new[lev];
    MultiFab& S_new = phi_new[lev];
    int nc = phi_new[lev].nComp();
    int ng = nghost;
    int ddir = domdir;
    // dt_calc.resize(1, 0.0);
    Real dt_calc = 0.0;
    Real umax = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	   FArrayBox uface[BL_SPACEDIM];

	   for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	   {
	       for (int i = 0; i < BL_SPACEDIM ; i++) {
            // const Box& validbox  = mfi.validbox();
		   const Box& bx = mfi.nodaltilebox(i);
		   uface[i].resize(bx,1);
           Print() << "dir= " << i << ", bx= " << bx << "dx = " << geom[lev].CellSize(0) << "\n";
	       }

            const Box& box = mfi.tilebox();
            const int* lo  = box.loVect();
            const int* hi  = box.hiVect();

	       get_face_velocity_dt(&lev, &cur_time,
			         AMREX_D_DECL(BL_TO_FORTRAN(uface[0]),
				     BL_TO_FORTRAN(uface[1]),
				     BL_TO_FORTRAN(uface[2])),
			         dx, prob_lo, &ddir, &umax, &nc,
                     BL_TO_FORTRAN_3D(S_new[mfi]), 
                     AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi));
            Print() << "umax= " << umax << "\n";
            if(ddir == 1){
                dt_calc = geom[lev].CellSize(0)/umax;
            }else{
                dt_calc = geom[lev].CellSize(1)/umax;
            }

		  dt_est = std::min(dt_est, dt_calc);
	   }
    }

    if (!local) {
	ParallelDescriptor::ReduceRealMin(dt_est);
    }

    dt_est *= cfl;

    Print() << "dt_est = " << dt_est << "\n";

    return dt_est;
}

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
    return {"ro","rou","rov","roE","pre"};
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
        int ncomp = 5;
        int nghost = 4;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);
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

// function to printout results (for testing purposes)
void 
AmrCoreAdv::printout (Box const& bx, Array4<Real> const& a)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);
   const int nf = a.nComp();

   // for     (int k = lo.z; k <= hi.z; ++k) {
   //   for   (int j = lo.y -  nghost; j <= hi.y + nghost; ++j) {
   //     for (int i = lo.x - nghost; i <= hi.x + nghost; ++i) {
   //      if (domdir == 1){
   //          // if (j==1){
   //          Print() << "i= " << i << ", j= " << j << 
   //                  // ",ro=" << a(i,j,k,ro) << ", rou=" << a(i,j,k,rou) 
   //              // << ", rov=" << a(i,j,k,rov) << "roE= " << a(i,j,k,roE) << 
   //              ", pre= " << a(i,j,k,pre)
   //              << "\n";
   //          // }
   //      }
   //      else {
   //          if (i==1){
   //          Print() << "j= " << j << ",ro=" << a(i,j,k,ro) << ", rou=" << a(i,j,k,rou) 
   //          << ", rov=" << a(i,j,k,rov) << "roE= " << a(i,j,k,roE) << ", pre= " << a(i,j,k,pre)
   //          << "\n";
   //          }
   //      }
   //     }
   //   }
   // }

   for     (int k = lo.z; k <= hi.z; ++k) {
     // for   (int j = lo.y -  nghost; j <= hi.y + nghost; ++j) {
       for (int i = lo.x - nghost; i <= hi.x + nghost; ++i) {
        if (domdir == 1){
            // if (j==1){
            Print() << "i= " << i << ", pre_lo = " << a(i,lo.y,k,pre) << ", " << a(i,hi.y+1,k,pre)
                    << ", pre_hi = " << a(i,hi.y,k,pre) << ", " << a(i,lo.y-1,k,pre) << "\n";
            // }
        }
        // else {
        //     if (i==1){
        //     Print() << "j= " << j << ",ro=" << a(i,j,k,ro) << ", rou=" << a(i,j,k,rou) 
        //     << ", rov=" << a(i,j,k,rov) << "roE= " << a(i,j,k,roE) << ", pre= " << a(i,j,k,pre)
        //     << "\n";
        //     }
        // }
       }
     // }
   }
}

// function to write results to file (for testing purposes)
void 
AmrCoreAdv::WriteTxtFileInit (int lev, Box const& bx, Array4<Real> const& a, int stepnum, const Geometry& geom)
{
   const auto lo = lbound(bx);
   const auto hi = ubound(bx);
   const int nf = a.nComp();

   std::string plotname = amrex::Concatenate("plt",stepnum,5);
   
   plotname = plotname + "d";
   plotname = amrex::Concatenate(plotname, domdir, 1);
   plotname = plotname + "l";
   plotname = amrex::Concatenate(plotname, lev, 1);

   std::string filename = plotname + ".txt";

   std::ofstream ofs(filename, std::ofstream::out);
    if (domdir == 1){
        Print(ofs) << "# x ro rou rov roE pre" << "\n";
    } else {
        Print(ofs) << "# y ro rou rov roE pre" << "\n";
    }
  

   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y -  nghost; j <= hi.y + nghost; ++j) {
       for (int i = lo.x - nghost; i <= hi.x + nghost; ++i) {
        if (domdir == 1){
            if (j==1){
            Real len = geom.ProbHi(0) - geom.ProbLo(0);
            Real x = geom.ProbLo(0) + (i + 0.5)*geom.CellSize(0);
            Print(ofs) << x/len << "\t" << a(i,j,k,ro) << "\t" << a(i,j,k,rou) << "\t"
            << a(i,j,k,rov) << "\t" << a(i,j,k,roE) << "\t" << a(i,j,k,pre) << "\n";
            }
        }
        else {
            if (i==1){
            Real len = geom.ProbHi(1) - geom.ProbLo(1);
            Real y = geom.ProbLo(1) + (j + 0.5)*geom.CellSize(1);
            Print(ofs) << y/len << "\t\t" << a(i,j,k,ro) << "\t" << a(i,j,k,rou) << "\t"
            << a(i,j,k,rov) << "\t" << a(i,j,k,roE) << "\t" << a(i,j,k,pre) << "\n";
            }
        }
       }
     }
   }

   ofs.close();
}

void
AmrCoreAdv::writetxtfile (int lev, Real cur_time, int stepnum)
{
    MultiFab& phi = phi_new[lev];
    const Geometry& geom1 = geom[lev];
    std::string plotname = amrex::Concatenate("plt",stepnum,5);
    plotname = plotname + "d";
    plotname = amrex::Concatenate(plotname, domdir, 1);
    plotname = plotname + "l";
    plotname = amrex::Concatenate(plotname, lev, 1);

    std::string filename = plotname + ".txt";

    std::ofstream ofs(filename, std::ofstream::out);

    if (domdir == 1){
        Print(ofs) << "# x ro rou rov roE pre" << "\n";
    } else {
        Print(ofs) << "# y ro rou rov roE pre" << "\n";
    }
    Print(ofs) << "# Time = " << cur_time << "\n";

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const auto lo = lbound(box);
        const auto hi = ubound(box);
        // FArrayBox& fab = state[mfi];
        // Array4<Real> const& a = fab.array();
        //  Function to print results
        FArrayBox& mfab = phi[mfi];
        Array4<Real> const& a = mfab.array();
        const int nf = a.nComp();
    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (domdir == 1){
                    if (j==1){
                        Real len = geom1.ProbHi(0) - geom1.ProbLo(0);
                        Real x = geom1.ProbLo(0) + (i + 0.5)*geom1.CellSize(0);
                        Print(ofs).SetPrecision(8) << std::left << std::setw(12) << x/len << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,ro)  << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,rou) << "\t"
                        << std::left << std::setw(12) << a(i,j,k,rov) << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,roE) << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,pre) << "\n";
                    }
                }
                else {
                    if (i==1){
                        Real len = geom1.ProbHi(1) - geom1.ProbLo(1);
                        Real y = geom1.ProbLo(1) + (j + 0.5)*geom1.CellSize(1);
                        Print(ofs).SetPrecision(8) << std::left << std::setw(12) << y/len << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,ro)  << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,rou) << "\t"
                        << std::left << std::setw(12) << a(i,j,k,rov) << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,roE) << "\t" 
                        << std::left << std::setw(12) << a(i,j,k,pre) << "\n";
                    }
                }
            }
        }
    }
}
   ofs.close();
}
