! $Id$

!

!  subroutines in the chosen set of physics modules.

!

module Equ
use Cdata
use Messages
use Boundcond
use Mpicomm
use Grid, only: calc_pencils_grid, get_grid_mn
implicit none
public :: pde, debug_imn_arrays, initialize_pencils
public :: impose_floors_ceilings
private
real, dimension(ny*nz,nx) :: dt1_max_array
contains
include 'pencil_init.inc'
subroutine pde(f,df,p)
use Chiral
use Chemistry
use Density
use Detonate, only: detonate_before_boundary
use Diagnostics
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: forcing_after_boundary
use General, only: ioptest, loptest
use GhostFold, only: fold_df, fold_df_3points
use Gpu
use Gravity
use Hydro
use Interstellar, only: interstellar_before_boundary, check_SN
use Magnetic
use Hypervisc_strict, only: hyperviscosity_strict
use Hyperresi_strict, only: hyperresistivity_strict
use NSCBC
use Particles_main
use Poisson
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_shock_profile, calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells, freeze_solid_cells,  dsolid_dt_integrate
use Special, only: special_before_boundary,special_after_boundary
use Sub
use Testfield
use Testflow
use Testscalar
use Viscosity, only: viscosity_after_boundary
use Grid, only: coarsegrid_interp
real, dimension (mx,my,mz,mfarray) :: f
real, dimension (mx,my,mz,mvar) :: df
type (pencil_case) :: p
intent(inout):: f
intent(out)  :: df,p
logical :: early_finalize
real, dimension(nx) :: pfreeze,pfreeze_int,pfreeze_ext
real, dimension(1)  :: mass_per_proc
integer :: iv
logical, dimension(npencils) :: lpenc_loc
headtt = headt .and. lfirst .and. lroot
if (headtt.or.ldebug) print*,'pde: ENTER'
if (headtt) call svn_id(  "$Id$")
ldiagnos   =lfirst.and.lout
l1davgfirst=lfirst.and.l1davg
l2davgfirst=lfirst.and.l2davg
l1dphiavg=lcylinder_in_a_box.and.l1davgfirst
lchemonly=.false.
if (ldiagnos   ) tdiagnos   =t
if (l1davgfirst) t1ddiagnos =t
if (l2davgfirst)  then
t2davgfirst=t
lpencil(i_rcyl_mn)=.true.
lpencil(i_z_mn)=.true.
endif
if (lfirst .and. lshift_datacube_x) then
call boundconds_x(f)
do  n=n1,n2; do m=m1,m2
f(:,m,n,:)=cshift(f(:,m,n,:),1,1)
enddo; enddo
endif
early_finalize=test_nonblocking.or.  leos_ionization.or.lradiation_ray.or.  lhyperviscosity_strict.or.lhyperresistivity_strict.or.  ltestscalar.or.ltestfield.or.ltestflow.or.  lparticles_spin.or.lsolid_cells.or.  lchemistry.or.lweno_transport .or. lbfield .or.  lvisc_smag .or.  lyinyang .or.  ncoarse>1
if (crash_file_dtmin_factor > 0.0) call output_crash_files(f)
call impose_floors_ceilings(f)
if (lparticles) call particles_boundconds(f)
if (lpointmasses) call boundconds_pointmasses
call calc_selfpotential(f)
if (.not. lchemistry) then
if (ldustdensity) call null_dust_vars(f)
if (ldustdensity .and. lmdvar .and. lfirst) call redist_mdbins(f)
endif
if (linterstellar) call interstellar_before_boundary(f)
if (ldensity.or.lboussinesq) call density_before_boundary(f)
if (lhydro.or.lhydro_kinematic) call hydro_before_boundary(f)
if (lmagnetic)     call magnetic_before_boundary(f)
call energy_before_boundary(f)
if (lshear)        call shear_before_boundary(f)
if (lchiral)       call chiral_before_boundary(f)
if (lspecial)      call special_before_boundary(f)
if (ltestflow)     call testflow_before_boundary(f)
if (ltestfield)    call testfield_before_boundary(f)
if (lparticles)    call particles_before_boundary(f)
if (ldetonate)     call detonate_before_boundary(f)
if (lparticles.and.lspecial) call particles_special_bfre_bdary(f)
if (lshock) call calc_shock_profile(f)
call boundconds_x(f)
if (nghost>0) then
if (ldebug) print*,'pde: before initiate_isendrcv_bdry'
call initiate_isendrcv_bdry(f)
if (early_finalize) then
call finalize_isendrcv_bdry(f)
if (lcoarse) call coarsegrid_interp(f)
call boundconds_y(f)
call boundconds_z(f)
endif
endif
if (lsolid_cells .and. lchemistry)  call chemspec_normalization_N2(f)
call update_solid_cells(f)
if (lhyperviscosity_strict)   call hyperviscosity_strict(f)
if (lhyperresistivity_strict) call hyperresistivity_strict(f)
if (ldynamical_diffusion) call set_dyndiff_coeff(f)
call fill_farray_pressure(f)
if (lfirst.and.ldt) then
if (dtmax/=0.0) then
if (lfractional_tstep_advance) then
dt1_max=1./(dt_incr*t)
else
dt1_max=1./dtmax
endif
else
dt1_max=0.0
endif
endif
if (leos_ionization.or.leos_temperature_ionization) call ioncalc(f)
if (lradiation_ray) call radtransfer(f)
if (lshock) call calc_shock_profile_simple(f)
call timing('pde','before hydro_after_boundary')
if (lhydro)                 call hydro_after_boundary(f)
if (lviscosity)             call viscosity_after_boundary(f)
if (lmagnetic)              call magnetic_after_boundary(f)
if (lenergy)                call energy_after_boundary(f)
if (lgrav)                  call gravity_after_boundary(f)
if (lforcing)               call forcing_after_boundary(f)
if (lpolymer)               call calc_polymer_after_boundary(f)
if (ltestscalar)            call testscalar_after_boundary(f)
if (ltestfield)             call testfield_after_boundary(f)
if (lpscalar)               call pscalar_after_boundary(f)
if (ldensity)               call density_after_boundary(f)
if (ltestflow)              call calc_ltestflow_nonlin_terms(f,df)
if (lspecial)               call special_after_boundary(f)
if (lchemistry .and. ldustdensity) then
call chemspec_normalization(f)
endif
if (lchemistry .and. ldensity) call calc_for_chem_mixture(f)
call timing('pde','after calc_for_chem_mixture')
if (lgpu) then
call rhs_gpu(f,itsub,early_finalize)
if (ldiagnos.or.l1davgfirst.or.l1dphiavg.or.l2davgfirst) then
call copy_farray_from_GPU(f)
call calc_all_module_diagnostics(f,p)
endif
else
call rhs_cpu(f,df,p,mass_per_proc,early_finalize)
endif
if (linterstellar) call check_SN(f)
if (lanelastic) call anelastic_after_mn(f,p,df,mass_per_proc)
call timing('pde','at the end of the mn_loop')
if (lsolid_cells) call dsolid_dt_integrate
if (igpotselfx/=0) then
call initiate_isendrcv_bdry(f,igpotselfx,igpotselfz)
call finalize_isendrcv_bdry(f,igpotselfx,igpotselfz)
call boundconds_x(f,igpotselfx,igpotselfz)
call boundconds_y(f,igpotselfx,igpotselfz)
call boundconds_z(f,igpotselfx,igpotselfz)
endif
if (lparticles) then
call particles_pde_blocks(f,df)
call particles_pde(f,df,p)
endif
if (lpointmasses) call pointmasses_pde(f,df)
if (lfold_df) then
if (lhydro .and. (.not. lpscalar) .and. (.not. lchemistry)) then
call fold_df(df,iux,iuz)
else
call fold_df(df,iux,mvar)
endif
endif
if (lfold_df_3points) then
call fold_df_3points(df,iux,mvar)
endif
do imn=1,ny*nz
n=nn(imn)
m=mm(imn)
if (any(lfreeze_varext).or.any(lfreeze_varint)) then
lpenc_loc=.false.
if (lcylinder_in_a_box.or.lcylindrical_coords) then
lpenc_loc(i_rcyl_mn)=.true.
else
lpenc_loc(i_r_mn)=.true.
endif
call calc_pencils_grid(f,p,lpenc_loc)
if (any(lfreeze_varint)) then
if (headtt) print*, 'pde: freezing variables for r < ', rfreeze_int,  ' : ', lfreeze_varint
if (lcylinder_in_a_box.or.lcylindrical_coords) then
if (wfreeze_int==0.0) then
where (p%rcyl_mn<=rfreeze_int)
pfreeze_int=0.0
elsewhere
pfreeze_int=1.0
endwhere
else
pfreeze_int=quintic_step(p%rcyl_mn,rfreeze_int,wfreeze_int,  SHIFT=fshift_int)
endif
else
if (wfreeze_int==0.0) then
where (p%r_mn<=rfreeze_int)
pfreeze_int=0.0
elsewhere
pfreeze_int=1.0
endwhere
else
pfreeze_int=quintic_step(p%r_mn,rfreeze_int,wfreeze_int,  SHIFT=fshift_int)
endif
endif
do iv=1,nvar
if (lfreeze_varint(iv)) df(l1:l2,m,n,iv)=pfreeze_int*df(l1:l2,m,n,iv)
enddo
endif
if (any(lfreeze_varext)) then
if (headtt) print*, 'pde: freezing variables for r > ', rfreeze_ext,  ' : ', lfreeze_varext
if (lcylinder_in_a_box.or.lcylindrical_coords) then
if (wfreeze_ext==0.0) then
where (p%rcyl_mn>=rfreeze_ext)
pfreeze_ext=0.0
elsewhere
pfreeze_ext=1.0
endwhere
else
pfreeze_ext=1.0-quintic_step(p%rcyl_mn,rfreeze_ext,wfreeze_ext,  SHIFT=fshift_ext)
endif
else
if (wfreeze_ext==0.0) then
where (p%r_mn>=rfreeze_ext)
pfreeze_ext=0.0
elsewhere
pfreeze_ext=1.0
endwhere
else
pfreeze_ext=1.0-quintic_step(p%r_mn,rfreeze_ext,wfreeze_ext,  SHIFT=fshift_ext)
endif
endif
do iv=1,nvar
if (lfreeze_varext(iv)) df(l1:l2,m,n,iv) = pfreeze_ext*df(l1:l2,m,n,iv)
enddo
endif
endif
if (any(lfreeze_varsquare)) then
if (headtt) print*, 'pde: freezing variables inside square : ',  lfreeze_varsquare
pfreeze=1.0-quintic_step(x(l1:l2),xfreeze_square,wfreeze,SHIFT=-1.0)* quintic_step(spread(y(m),1,nx),yfreeze_square,-wfreeze,SHIFT=-1.0)
do iv=1,nvar
if (lfreeze_varsquare(iv)) df(l1:l2,m,n,iv) = pfreeze*df(l1:l2,m,n,iv)
enddo
endif
if (lfrozen_bcs_x) then
if (.not. lperi(1)) then
if (lfirst_proc_x) where (lfrozen_bot_var_x(1:nvar)) df(l1,m,n,1:nvar) = 0.
if (llast_proc_x) where (lfrozen_top_var_x(1:nvar)) df(l2,m,n,1:nvar) = 0.
endif
endif
if (lfrozen_bcs_y) then
if (.not. lperi(2)) then
if (lfirst_proc_y .and. (m == m1)) then
do iv=1,nvar
if (lfrozen_bot_var_y(iv)) df(l1:l2,m,n,iv) = 0.
enddo
endif
if (llast_proc_y .and. (m == m2)) then
do iv=1,nvar
if (lfrozen_top_var_y(iv)) df(l1:l2,m,n,iv) = 0.
enddo
endif
endif
endif
if (lfrozen_bcs_z) then
if (.not. lperi(3)) then
if (lfirst_proc_z .and. (n == n1)) then
do iv=1,nvar
if (lfrozen_bot_var_z(iv)) df(l1:l2,m,n,iv) = 0.
enddo
endif
if (llast_proc_z .and. (n == n2)) then
do iv=1,nvar
if (lfrozen_top_var_z(iv)) df(l1:l2,m,n,iv) = 0.
enddo
endif
endif
endif
call freeze_solid_cells(df)
enddo
if (lnscbc) call nscbc_boundtreat(f,df)
if (ldiagnos) then
call diagnostic(fname,nname)
call diagnostic(fname_keep,nname,lcomplex=.true.)
endif
if (l1davgfirst) then
if (lwrite_xyaverages) call xyaverages_z
if (lwrite_xzaverages) call xzaverages_y
if (lwrite_yzaverages) call yzaverages_x
endif
if (l1dphiavg) call phizaverages_r
if (l2davgfirst) then
if (lwrite_yaverages)   call yaverages_xz
if (lwrite_zaverages)   call zaverages_xy
if (lwrite_phiaverages) call phiaverages_rz
endif
if (.not.l2davgfirst.and.ldiagnos.and.ldiagnos_need_zaverages) then
if (lwrite_zaverages) call zaverages_xy
endif
if (ldiagnos) then
if (lmagnetic) call calc_mfield
if (lhydro)    call calc_mflow
if (lpscalar)  call calc_mpscalar
endif
lwrite_prof=.false.
endsubroutine pde
subroutine calc_all_module_diagnostics(f,p)
use Density, only: calc_diagnostics_density
use Energy, only: calc_diagnostics_energy
use Hydro, only: calc_diagnostics_hydro
use Magnetic, only: calc_diagnostics_magnetic
use Forcing, only: calc_diagnostics_forcing
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
type (pencil_case)                ,intent(INOUT) :: p
integer :: nyz, imn
nyz=ny*nz
do imn=1,nyz
n=nn(imn)
m=mm(imn)
lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
if (lcoarse_mn) then
lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
if (ninds(0,m,n)<=0) cycle
endif
lfirstpoint=(imn==1)
llastpoint=(imn==nyz)
call calc_all_pencils(f,p)
call calc_diagnostics_density(f,p)
call calc_diagnostics_energy(f,p)
call calc_diagnostics_hydro(f,p)
call calc_diagnostics_magnetic(f,p)
if (lforcing_cont) call calc_diagnostics_forcing(p)
enddo
endsubroutine calc_all_module_diagnostics
subroutine calc_all_pencils(f,p)
use Diagnostics, only: calc_phiavg_profile
use BorderProfiles, only: calc_pencils_borderprofiles
use Chiral
use Chemistry
use Cosmicray
use CosmicrayFlux
use Density
use Dustvelocity
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: calc_pencils_forcing
use Gravity
use Heatflux
use Hydro
use Lorenz_gauge
use Magnetic
use NeutralDensity
use NeutralVelocity
use Particles_main
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_pencils_shock, calc_shock_profile,  calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells_pencil
use Special, only: calc_pencils_special, dspecial_dt
use Ascalar
use Testfield
use Testflow
use Testscalar
use Viscosity, only: calc_pencils_viscosity
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
type (pencil_case)                ,intent(INOUT) :: p
call update_solid_cells_pencil(f)
if (.not. lcartesian_coords .or. .not.all(lequidist)) call get_grid_mn
call calc_pencils_grid(f,p)
if ((l2davgfirst.and.lwrite_phiaverages )  .or.  (l1dphiavg  .and.lwrite_phizaverages))   call calc_phiavg_profile(p)
call calc_pencils_hydro(f,p)
call calc_pencils_density(f,p)
if (lpscalar)         call calc_pencils_pscalar(f,p)
if (lascalar)         call calc_pencils_ascalar(f,p)
call calc_pencils_eos(f,p)
if (lshock)           call calc_pencils_shock(f,p)
if (lchemistry)       call calc_pencils_chemistry(f,p)
call calc_pencils_energy(f,p)
if (lviscosity)       call calc_pencils_viscosity(f,p)
if (lforcing_cont)    call calc_pencils_forcing(f,p)
if (llorenz_gauge)    call calc_pencils_lorenz_gauge(f,p)
if (lmagnetic)        call calc_pencils_magnetic(f,p)
if (lpolymer)         call calc_pencils_polymer(f,p)
if (lgrav)            call calc_pencils_gravity(f,p)
if (lselfgravity)     call calc_pencils_selfgravity(f,p)
if (ldustvelocity)    call calc_pencils_dustvelocity(f,p)
if (ldustdensity)     call calc_pencils_dustdensity(f,p)
if (lneutralvelocity) call calc_pencils_neutralvelocity(f,p)
if (lneutraldensity)  call calc_pencils_neutraldensity(f,p)
if (lcosmicray)       call calc_pencils_cosmicray(f,p)
if (lcosmicrayflux)   call calc_pencils_cosmicrayflux(f,p)
if (lchiral)          call calc_pencils_chiral(f,p)
if (lradiation)       call calc_pencils_radiation(f,p)
if (lshear)           call calc_pencils_shear(f,p)
if (lborder_profiles) call calc_pencils_borderprofiles(f,p)
if (lpointmasses)     call calc_pencils_pointmasses(f,p)
if (lparticles)       call particles_calc_pencils(f,p)
if (lspecial)         call calc_pencils_special(f,p)
if (lheatflux)        call calc_pencils_heatflux(f,p)
endsubroutine calc_all_pencils
subroutine rhs_cpu(f,df,p,mass_per_proc,early_finalize)
use Diagnostics
use Chiral
use Chemistry
use Cosmicray
use CosmicrayFlux
use Density
use Dustvelocity
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: calc_pencils_forcing, calc_diagnostics_forcing
use General, only: notanumber
use GhostFold, only: fold_df, fold_df_3points
use Gravity
use Heatflux
use Hydro
use Lorenz_gauge
use Magnetic
use NeutralDensity
use NeutralVelocity
use Particles_main
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_pencils_shock, calc_shock_profile,  calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells_pencil, dsolid_dt
use Special, only: calc_pencils_special, dspecial_dt
use Sub, only: sum_mn
use Ascalar
use Testfield
use Testflow
use Testscalar
use Viscosity, only: calc_pencils_viscosity
use omp_lib
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
real, dimension (mx,my,mz,mvar)   ,intent(OUT  ) :: df
type (pencil_case)                ,intent(INOUT) :: p
real, dimension(1)                ,intent(INOUT) :: mass_per_proc
logical                           ,intent(IN   ) :: early_finalize
integer :: nyz
real, dimension (nx,3) :: df_iuu_pencil
integer :: num_omp_ranks, omp_rank
real :: dt1_poly_relax, dt1_preac
type (pencil_case), dimension(ny*nz) :: p_array
real, dimension (mx,my,mz,mfarray) :: f_temp
real, dimension (mx,my,mz,mvar) :: df_temp
real, dimension (nx,3) :: df_iuu_pencil_temp
nyz=ny*nz
call OMP_set_num_threads(1)
calc_pencils_loop: do imn=1,nyz
n=nn(imn)
m=mm(imn)
lfirstpoint=(imn==1)
llastpoint=(imn==nyz)
lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
if (lcoarse_mn) then
lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
if (ninds(0,m,n)<=0) cycle
endif
call timing('pde','before lanelastic',mnloop=.true.)
if (lanelastic) then
df_iuu_pencil = df(l1:l2,m,n,iuu:iuu+2)
df(l1:l2,m,n,iuu:iuu+2)=0.0
endif
if (.not.early_finalize.and.necessary(imn)) then
call finalize_isendrcv_bdry(f)
call boundconds_y(f)
call boundconds_z(f)
endif
call timing('pde','finished boundconds_z',mnloop=.true.)
if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
advec_cs2=0.0
maxadvec=0.
if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or.  (ltestfield_z.and.iuutest>0))  advec2=0.
if (ldensity.or.lviscosity.or.lmagnetic.or.lenergy)  advec2_hypermesh=0.0
maxdiffus=0.
maxdiffus2=0.
maxdiffus3=0.
maxsrc=0.
endif
call calc_all_pencils(f,p)
p_array(imn) = p
enddo calc_pencils_loop
call OMP_set_num_threads(4)
p = p_array(1)
call calc_modules(f,df,p,1,nyz)
call set_dt1_max(p)
!$omp parallel default(firstprivate) reduction(max:dt1_max) shared(f,df,fname,fnamer,fnamex,fnamey,fnamez,fnamexy,fnamexz,fnamerz,itype_name)

mn_loop: do imn=nyz-1,2,-1
!$omp critical

p = p_array(imn)
call calc_modules(f,df,p,imn,nyz)
call set_dt1_max(p)
!$omp end critical

enddo mn_loop
!$omp endparallel

p = p_array(ny*nz)
call calc_modules(f,df,p,ny*nz,nyz)
call set_dt1_max(p)
endsubroutine rhs_cpu
subroutine calc_rest(df_iuu_pencil, num_omp_ranks, omp_rank, dt1_advec,  dt1_diffus, dt1_src, dt1_reac, dt1_poly_relax, dt1_preac, f, df, p, mass_per_proc, early_finalize, imn,nyz,advec2_array,advec_cs2_array,maxadvec_array,advec2_hypermesh_array,  cdt_array,iproc_world_array,dt1_diffus_array,maxdiffus_array,cdtv_array,maxdiffus2_array,cdtv2_array,maxdiffus3_array, cdtv3_array,dt1_src_array,dt1_max_array,dt1_advec_array,dt1_poly_relax_array,dt1_reac_array,idiag_dtv_array,  idiag_dtdiffus_array,idiag_dtdiffus2_array,idiag_dtdiffus3_array,idiag_Rmesh_array,  idiag_Rmesh3_array,tini_array,pi_1_array,pi5_1_array)
use Diagnostics
use Chiral
use Chemistry
use Cosmicray
use CosmicrayFlux
use Density
use Dustvelocity
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: calc_pencils_forcing, calc_diagnostics_forcing
use General, only: notanumber
use GhostFold, only: fold_df, fold_df_3points
use Gravity
use Heatflux
use Hydro
use Lorenz_gauge
use Magnetic
use NeutralDensity
use NeutralVelocity
use Particles_main
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_pencils_shock, calc_shock_profile,  calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells_pencil, dsolid_dt
use Special, only: calc_pencils_special, dspecial_dt
use Sub, only: sum_mn
use Ascalar
use Testfield
use Testflow
use Testscalar
use Viscosity, only: calc_pencils_viscosity
use omp_lib
real, dimension (nx,3), intent(INOUT) :: df_iuu_pencil
integer, intent(INOUT) :: num_omp_ranks, omp_rank
real, dimension(nx), intent(INOUT) :: dt1_advec, dt1_diffus, dt1_src, dt1_reac
real, intent(INOUT) :: dt1_poly_relax, dt1_preac
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
real, dimension (mx,my,mz,mvar)   ,intent(INOUT  ) :: df
type (pencil_case)                ,intent(INOUT) :: p
real, dimension(1)                ,intent(INOUT) :: mass_per_proc
logical                           ,intent(IN   ) :: early_finalize
integer, intent(IN) :: imn
integer, intent(IN) :: nyz
real, dimension (ny*nz,nx), intent(INOUT) :: maxadvec_array, advec2_array, advec2_hypermesh_array, advec_cs2_array,  dt1_advec_array, dt1_diffus_array, maxdiffus_array, maxdiffus2_array, maxdiffus3_array,  dt1_src_array, dt1_max_array, dt1_reac_array
real, dimension (ny*nz,1), intent(INOUT) :: cdt_array, cdtv_array, cdtv2_array, cdtv3_array,  dt1_poly_relax_array, tini_array, pi_1_array, pi5_1_array
integer, dimension (ny*nz,1), intent(INOUT) :: idiag_dtv_array, idiag_Rmesh_array, idiag_Rmesh3_array,  iproc_world_array, idiag_dtdiffus_array, idiag_dtdiffus2_array, idiag_dtdiffus3_array
n=nn(imn)
m=mm(imn)
lfirstpoint=(imn==1)
llastpoint=(imn==nyz)
lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
if (lcoarse_mn) then
lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
if (ninds(0,m,n)<=0) return
endif
if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
advec2=advec2+advec_cs2
if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or.  (ltestfield_z.and.iuutest>0))  maxadvec=maxadvec+sqrt(advec2)
if (ldensity.or.lhydro.or.lmagnetic.or.lenergy) maxadvec=maxadvec+sqrt(advec2_hypermesh)
if (any(lfreeze_varint)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
where (p%rcyl_mn<=rfreeze_int)
maxadvec=0.0
maxdiffus=0.0
endwhere
else
where (p%r_mn<=rfreeze_int)
maxadvec=0.0
maxdiffus=0.0
endwhere
endif
endif
if (any(lfreeze_varext)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
where (p%rcyl_mn>=rfreeze_ext)
maxadvec=0.0
maxdiffus=0.0
endwhere
else
where (p%r_mn>=rfreeze_ext)
maxadvec=0.0
maxdiffus=0.0
endwhere
endif
endif
dt1_advec = maxadvec/cdt
if (notanumber(dt1_advec)) then
print*, 'pde: dt1_advec contains a NaN at iproc=', iproc_world
if (lenergy) print*, 'advec_cs2  =',advec_cs2
call fatal_error_local('pde','')
endif
dt1_diffus = maxdiffus/cdtv + maxdiffus2/cdtv2 + maxdiffus3/cdtv3
dt1_src = maxsrc/cdtsrc
dt1_max = max(dt1_max, sqrt(dt1_advec**2 + dt1_diffus**2 + dt1_src**2))
if (lchemistry .and. .not.llsode) then
dt1_reac = reac_chem/cdtc
dt1_max = max(dt1_max,dt1_reac)
endif
if (lpolymer) then
dt1_poly_relax = 1./(trelax_poly*cdt_poly)
dt1_max = max(dt1_max,dt1_poly_relax)
endif
if (any(lfreeze_varint)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
else
endif
endif
if (any(lfreeze_varext)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
else
endif
endif
if (ldiagnos) then
if (idiag_dtv/=0) call max_mn_name(maxadvec/cdt,idiag_dtv,l_dt=.true.)
if (idiag_dtdiffus/=0) call max_mn_name(maxdiffus/cdtv,idiag_dtdiffus,l_dt=.true.)
if (idiag_dtdiffus2/=0) call max_mn_name(maxdiffus2/cdtv2,idiag_dtdiffus2,l_dt=.true.)
if (idiag_dtdiffus3/=0) call max_mn_name(maxdiffus3/cdtv3,idiag_dtdiffus3,l_dt=.true.)
if (idiag_Rmesh/=0) call max_mn_name(pi_1*maxadvec/(maxdiffus+tini),idiag_Rmesh)
if (idiag_Rmesh3/=0) call max_mn_name(pi5_1*maxadvec/(maxdiffus3+tini),idiag_Rmesh3)
call max_mn_name(maxadvec,idiag_maxadvec)
endif
endif
if (lanelastic) then
call calc_pencils_density(f,p)
f(l1:l2,m,n,irhs) = p%rho*df(l1:l2,m,n,iuu)
f(l1:l2,m,n,irhs+1) = p%rho*df(l1:l2,m,n,iuu+1)
f(l1:l2,m,n,irhs+2) = p%rho*df(l1:l2,m,n,iuu+2)
df(l1:l2,m,n,iuu:iuu+2) = df_iuu_pencil(1:nx,1:3) + df(l1:l2,m,n,iuu:iuu+2)
call sum_mn(p%rho,mass_per_proc(1))
endif
call timing('pde','end of mn loop',mnloop=.true.)
headtt=.false.
endsubroutine calc_rest
subroutine calc_modules(f,df,p,imn,nyz)
use Diagnostics
use Chiral
use Chemistry
use Cosmicray
use CosmicrayFlux
use Density
use Dustvelocity
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: calc_pencils_forcing, calc_diagnostics_forcing
use General, only: notanumber
use GhostFold, only: fold_df, fold_df_3points
use Gravity
use Heatflux
use Hydro
use Lorenz_gauge
use Magnetic
use NeutralDensity
use NeutralVelocity
use Particles_main
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_pencils_shock, calc_shock_profile,  calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells_pencil, dsolid_dt
use Special, only: calc_pencils_special, dspecial_dt
use Sub, only: sum_mn
use Ascalar
use Testfield
use Testflow
use Testscalar
use Viscosity, only: calc_pencils_viscosity
use omp_lib
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
real, dimension (mx,my,mz,mvar)   ,intent(OUT  ) :: df
type (pencil_case)                ,intent(INOUT) :: p
integer, intent(IN) :: imn, nyz
n=nn(imn)
m=mm(imn)
lfirstpoint=(imn==1)
llastpoint=(imn==nyz)
lcoarse_mn=lcoarse.and.mexts(1)<=m.and.m<=mexts(2)
if (lcoarse_mn) then
lcoarse_mn=lcoarse_mn.and.ninds(0,m,n)>0
if (ninds(0,m,n)<=0) return
endif
if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
advec_cs2=0.0
maxadvec=0.
if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or.  (ltestfield_z.and.iuutest>0))  advec2=0.
if (ldensity.or.lviscosity.or.lmagnetic.or.lenergy)  advec2_hypermesh=0.0
maxdiffus=0.
maxdiffus2=0.
maxdiffus3=0.
maxsrc=0.
endif
call duu_dt(f,df,p)
call dlnrho_dt(f,df,p)
call denergy_dt(f,df,p)
if (lmagnetic) call daa_dt(f,df,p)
if (llorenz_gauge) call dlorenz_gauge_dt(f,df,p)
if (lpolymer) call dpoly_dt(f,df,p)
if (ltestscalar) call dcctest_dt(f,df,p)
if (ltestfield) call daatest_dt(f,df,p)
if (ltestflow) call duutest_dt(f,df,p)
if (lpscalar) call dlncc_dt(f,df,p)
if (lascalar) call dacc_dt(f,df,p)
if (ldustvelocity) call duud_dt(f,df,p)
if (ldustdensity) call dndmd_dt(f,df,p)
if (lneutraldensity) call dlnrhon_dt(f,df,p)
if (lneutralvelocity) call duun_dt(f,df,p)
if (lgrav) then
if (lhydro.or.ldustvelocity.or.lneutralvelocity)  call duu_dt_grav(f,df,p)
endif
if (lselfgravity) call duu_dt_selfgrav(f,df,p)
if (lcosmicray) call decr_dt(f,df,p)
if (lcosmicrayflux) call dfcr_dt(f,df,p)
if (lchiral) call dXY_chiral_dt(f,df,p)
if (lradiation_fld) call de_dt(f,df,p)
if (lchemistry) call dchemistry_dt(f,df,p)
if (lheatflux) call dheatflux_dt(f,df,p)
if (lforcing_cont) call calc_diagnostics_forcing(p)
if (lspecial) call dspecial_dt(f,df,p)
if (lradiation_ray.and.lenergy) then
call radiative_cooling(f,df,p)
call radiative_pressure(f,df,p)
endif
if (lsolid_cells) call dsolid_dt(f,df,p)
if (lshear) call shearing(f,df,p)
if (lparticles) call particles_pde_pencil(f,df,p)
if (lpointmasses) call pointmasses_pde_pencil(f,df,p)
if (ldiagnos) then
if (lhydro) call df_diagnos_hydro(df,p)
if (lmagnetic) call df_diagnos_magnetic(df,p)
endif
if (l2davgfirst) then
call phisum_mn_name_rz(p%rcyl_mn,idiag_rcylmphi)
call phisum_mn_name_rz(p%phi_mn,idiag_phimphi)
call phisum_mn_name_rz(p%z_mn,idiag_zmphi)
call phisum_mn_name_rz(p%r_mn,idiag_rmphi)
endif
if (ltime_integrals.and.llast.and..not.lgpu) then
if (lhydro)    call time_integrals_hydro(f,p)
if (lmagnetic) call time_integrals_magnetic(f,p)
endif
endsubroutine calc_modules
subroutine calc_dt(df_iuu_pencil, num_omp_ranks, omp_rank, dt1_advec,  dt1_diffus, dt1_src, dt1_reac, dt1_poly_relax, dt1_preac, f, df, p, mass_per_proc, early_finalize, imn, nyz)
use Diagnostics
use Chiral
use Chemistry
use Cosmicray
use CosmicrayFlux
use Density
use Dustvelocity
use Dustdensity
use Energy
use EquationOfState
use Forcing, only: calc_pencils_forcing, calc_diagnostics_forcing
use General, only: notanumber
use GhostFold, only: fold_df, fold_df_3points
use Gravity
use Heatflux
use Hydro
use Lorenz_gauge
use Magnetic
use NeutralDensity
use NeutralVelocity
use Particles_main
use Pscalar
use PointMasses
use Polymer
use Radiation
use Selfgravity
use Shear
use Shock, only: calc_pencils_shock, calc_shock_profile,  calc_shock_profile_simple
use Solid_Cells, only: update_solid_cells_pencil, dsolid_dt
use Special, only: calc_pencils_special, dspecial_dt
use Sub, only: sum_mn
use Ascalar
use Testfield
use Testflow
use Testscalar
use Viscosity, only: calc_pencils_viscosity
use omp_lib
real, dimension (nx,3), intent(INOUT) :: df_iuu_pencil
integer, intent(INOUT) :: num_omp_ranks, omp_rank
real, dimension(nx), intent(INOUT) :: dt1_advec, dt1_diffus, dt1_src, dt1_reac
real, intent(INOUT) :: dt1_poly_relax, dt1_preac
real, dimension (mx,my,mz,mfarray),intent(INOUT) :: f
real, dimension (mx,my,mz,mvar)   ,intent(INOUT  ) :: df
type (pencil_case)                ,intent(INOUT) :: p
real, dimension(1)                ,intent(INOUT) :: mass_per_proc
logical                           ,intent(IN   ) :: early_finalize
integer, intent(INOUT) :: nyz
integer, intent(IN) :: imn
call timing('pde','before lanelastic',mnloop=.true.)
if (lanelastic) then
df_iuu_pencil = df(l1:l2,m,n,iuu:iuu+2)
df(l1:l2,m,n,iuu:iuu+2)=0.0
endif
endsubroutine calc_dt
subroutine debug_imn_arrays
open(1,file=trim(directory)//'/imn_arrays.dat')
do imn=1,ny*nz
if (necessary(imn)) write(1,'(a)') '----necessary=.true.----'
write(1,'(4i6)') imn,mm(imn),nn(imn)
enddo
close(1)
endsubroutine debug_imn_arrays
subroutine output_crash_files(f)
use Snapshot
real, dimension(mx,my,mz,mfarray) :: f
integer, save :: icrash=0
character (len=10) :: filename
character (len=1) :: icrash_string
if ( (it>1) .and. lfirst .and. (dt<=crash_file_dtmin_factor*dtmin) ) then
write(icrash_string, fmt='(i1)') icrash
filename='crash'//icrash_string//'.dat'
call wsnap(filename,f,mvar_io,ENUM=.false.)
if (lroot) then
print*, 'Time-step is very low - writing '//trim(filename)
print*, '(it, itsub=', it, itsub, ')'
print*, '(t, dt=', t, dt, ')'
endif
icrash=icrash+1
icrash=mod(icrash,10)
endif
endsubroutine output_crash_files
subroutine set_dyndiff_coeff(f)
use Density,   only: dynamical_diffusion
use Energy,    only: dynamical_thermal_diffusion
use Magnetic,  only: dynamical_resistivity
use Sub,       only: find_max_fvec, find_rms_fvec
use Viscosity, only: dynamical_viscosity
real, dimension(mx,my,mz,mfarray), intent(in) :: f
real :: uc
if (ldyndiff_useumax) then
uc = find_max_fvec(f, iuu)
else
uc = find_rms_fvec(f, iuu)
endif
if (ldensity)                      call dynamical_diffusion(uc)
if (lmagnetic .and. .not. lbfield) call dynamical_resistivity(uc)
if (lenergy)                       call dynamical_thermal_diffusion(uc)
if (lviscosity)                    call dynamical_viscosity(uc)
endsubroutine set_dyndiff_coeff
subroutine impose_floors_ceilings(f)
use Cosmicray, only: impose_ecr_floor
use Density, only: impose_density_floor,impose_density_ceiling
use Dustdensity, only: impose_dustdensity_floor
use Energy, only: impose_energy_floor
use Hydro, only: impose_velocity_ceiling
real, dimension(mx,my,mz,mfarray), intent(inout) :: f
call impose_density_floor(f)
call impose_density_ceiling(f)
call impose_velocity_ceiling(f)
call impose_energy_floor(f)
call impose_dustdensity_floor(f)
call impose_ecr_floor(f)
endsubroutine impose_floors_ceilings
subroutine set_dt1_max(p)
use General, only: notanumber
type (pencil_case), intent(IN) :: p
real, dimension(nx) :: dt1_max_loc, dt1_advec, dt1_diffus, dt1_src
real :: dt1_preac
if (lfirst.and.ldt.and.(.not.ldt_paronly)) then
advec2=advec2+advec_cs2
if (lenergy.or.ldensity.or.lmagnetic.or.lradiation.or.lneutralvelocity.or.lcosmicray.or.  (ltestfield_z.and.iuutest>0))  maxadvec=maxadvec+sqrt(advec2)
if (ldensity.or.lhydro.or.lmagnetic.or.lenergy) maxadvec=maxadvec+sqrt(advec2_hypermesh)
dt1_advec = maxadvec/cdt
dt1_diffus = maxdiffus/cdtv + maxdiffus2/cdtv2 + maxdiffus3/cdtv3
dt1_src = maxsrc/cdtsrc
dt1_max_loc = sqrt(dt1_advec**2 + dt1_diffus**2 + dt1_src**2)
if (ldustdensity) dt1_max_loc = max(dt1_max_loc,reac_dust/cdtc)
if (lchemistry .and. .not.llsode) then
dt1_max_loc = max(dt1_max_loc,reac_chem/cdtc)
endif
if (lpolymer) dt1_max_loc = max(dt1_max_loc,1./(trelax_poly*cdt_poly))
if (any(lfreeze_varint)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
where (p%rcyl_mn<=rfreeze_int)
dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
endwhere
else
where (p%r_mn<=rfreeze_int)
dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
endwhere
endif
endif
if (any(lfreeze_varext)) then
if (lcylinder_in_a_box.or.lcylindrical_coords) then
where (p%rcyl_mn>=rfreeze_ext)
dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
endwhere
else
where (p%r_mn>=rfreeze_ext)
dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
endwhere
endif
endif
if (any(lfreeze_varsquare).and.y(m)>yfreeze_square) then
where (x(l1:l2)>xfreeze_square)
dt1_max_loc=0.; maxadvec=0.; maxdiffus=0.; maxdiffus2=0.; maxdiffus3=0.
endwhere
endif
dt1_max=max(dt1_max,dt1_max_loc)
if (notanumber(maxadvec)) then
print*, 'pde: maxadvec contains a NaN at iproc=', iproc_world
if (lenergy) print*, 'advec_cs2  =',advec_cs2
call fatal_error_local('pde','')
endif
dt1_max_array(imn,:) = dt1_max
endif
endsubroutine set_dt1_max
endmodule Equ