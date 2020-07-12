#ifndef _AmrCoreAdv_F_H_
#define _AmrCoreAdv_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    void initdata(const int* level, const amrex_real* time, const int* nc,
    	  const int* domdir,
		  const int* lo, const int* hi,
		  BL_FORT_FAB_ARG_3D(state),
		  const amrex_real* dx, const amrex_real* problo, const amrex_real* probhi,
		  const amrex_real* ro2ro1, const amrex_real* v2v1, const amrex_real* p2p1);

    void get_face_velocity(const int* level, const amrex_real* time, 
    			const int* lo, const int* hi,
			   AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
				  BL_FORT_FAB_ARG(yvel),
				  BL_FORT_FAB_ARG(zvel)),
			   const amrex_real* dx, const amrex_real* problo,
			   BL_FORT_FAB_ARG_3D(state),
			   const int* nc, const int* ddir);

    void get_face_velocity_dt(const int* level, const amrex_real* time, 
			   AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
				  BL_FORT_FAB_ARG(yvel),
				  BL_FORT_FAB_ARG(zvel)),
			   // BL_FORT_FAB_ARG_3D(state), const int* nc,
			   const amrex_real* dx, const amrex_real* problo,
			   const int* domdir, amrex_real* umaxs,
			   const int* nc, BL_FORT_FAB_ARG_3D(phi),
			   const int* lo, const int* hi);

    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const amrex_real* dx, const amrex_real* problo,
		     const amrex_real* time, const amrex_real* phierr, const int* nc,
		     const int* domdir);

    void advect(const int* level, const amrex_real* time, const int* rk, const int* rk_max,
    	 const int* fct_step,
    	const int* ddir, const int* nc,
    	const int* lo, const int* hi,
    	BL_FORT_FAB_ARG_3D(stateold),
		BL_FORT_FAB_ARG_3D(statex),
		BL_FORT_FAB_ARG_3D(statey),
		BL_FORT_FAB_ARG_3D(stateout),
		AMREX_D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
		       const BL_FORT_FAB_ARG_3D(yvel),
		       const BL_FORT_FAB_ARG_3D(zvel)),
		AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx),
		       BL_FORT_FAB_ARG_3D(fy),
		       BL_FORT_FAB_ARG_3D(fz)),
		const amrex_real* dx, const amrex_real* dt);

    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const amrex_real* dx, const amrex_real* glo, 
		  const amrex_real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
