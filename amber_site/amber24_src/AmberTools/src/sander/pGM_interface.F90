#include "../include/dprec.fh"
#include "../include/assert.fh"

!-------------------------------------------------------------------------------
module pol_gauss_interface
  implicit none
  private

  public POL_GAUSS_readparm,POL_GAUSS_deallocate,pGM_NonBond_eval
#ifdef MPI
  public POL_GAUSS_bcast
#endif

contains
!-------------------------------------------------------------------------------
#ifdef MPI
subroutine POL_GAUSS_bcast(numatoms)
  use pol_gauss_mdin, only : pGM_mdin_bcast
  use pol_gauss_multipoles, only : pGM_MPOLE_bcast
  use pol_gauss_induced, only : pGM_INDUCED_bcast
  implicit none

  integer,intent(in) :: numatoms

  call pGM_mdin_bcast()
  call pGM_MPOLE_bcast()
  call pGM_INDUCED_bcast(numatoms)
end subroutine POL_GAUSS_bcast
#endif
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_readparm(nf,numatoms)
  use pol_gauss_mdin, only : ipgm
  use pol_gauss_multipoles, only : pGM_MPOLE_readparm
  use pol_gauss_induced, only : pGM_INDUCED_readparm
  implicit none

  integer,intent(in) :: nf,numatoms

  integer :: mpole_valid,polar_valid

  ! only proceed if ipgm has been turned on
  if ( ipgm == 0 ) return

  call pGM_check_parm_legal(nf)
  mpole_valid = pGM_MPOLE_readparm(nf,numatoms)
  polar_valid = pGM_INDUCED_readparm(nf,numatoms)
  if ( mpole_valid /= 1 .or. polar_valid /= 1 ) then
    write(6,*)'POL_GAUSS parameter reading failed '
    call mexit(6,1)
  end if

end subroutine POL_GAUSS_readparm
!-------------------------------------------------------------------------------
subroutine pGM_check_parm_legal(nf)
  use pol_gauss_mdin, only : ipgm
  implicit none
  integer,intent(in) :: nf

  integer :: iok,ionerr,indicator
  character(len=80) :: fmt,fmtin,dtype

  fmtin = '(I5)'
  dtype = 'POL_GAUSS_FORCEFIELD'
  ionerr = 1 ! not fatal if missing
  call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
  if ( iok == 0 )then !this data type found in prmtop
    read(nf,fmt)indicator
    if ( ipgm /= 1 )then
       write(6,*)'POL_GAUSS style prmtop, but pgm NOT set to one'
       call mexit(6,1)
    endif
    return ! pGM prmtop, pGM set to 1
  else
    if ( ipgm == 1 )then
       write(6,*)'NOT an POL_GAUSS style prmtop, but pgm IS set to one'
       call mexit(6,1)
    endif
    return ! NON pGM prmtop, pGM NOT set to 1
  endif
end subroutine pGM_check_parm_legal
!-------------------------------------------------------------------------------
subroutine POL_GAUSS_deallocate()
  use pol_gauss_multipoles, only : pGM_MPOLE_deallocate
  call pGM_MPOLE_deallocate()
  !call pGM_INDUCED_deallocate()
  !call pGM_RECIP_deallocate()
end subroutine POL_GAUSS_deallocate
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_eval(numatoms,crd,frc,sander_vir,x,nummols,molsiz,amass,ipairs, &
                       ntypes,iac,ico,cn1,cn2, &
                       evdw,eelt,diprms,dipiter)
  use pol_gauss_mdin, only : pol_gauss_verbose
  use pol_gauss_multipoles, only : pGM_MPOLE_local_to_global,start_multipoles,end_multipoles, &
                               global_multipole,MAXMP,coulomb_const_kcal_per_mole
  use pol_gauss_induced,only : pGM_INDUCED_eval
  use nblist, only: recip,adjust_imagcrds, map_coords
  use pol_gauss_recip, only : pGM_RECIP_allocate
  implicit none

  integer,intent(in) :: numatoms
  _REAL_ :: crd(3,*), amass(*)
  _REAL_,intent(inout) :: frc(3,*),sander_vir(4)
  _REAL_,intent(in) :: x(*)
  integer nummols,molsiz(*)
  integer, intent(in) :: ipairs(*)
  integer,intent(in) :: ntypes,iac(*),ico(*)
  _REAL_,intent(in) :: cn1(*),cn2(*)
  _REAL_,intent(out) :: evdw,eelt,dipiter,diprms

# include "box.h"
# include "extra.h"
# include "ew_erfc_spline.h"
#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
#endif

  integer ier
  _REAL_ :: pgm_virial(3,3), vdw_virial(3,3)

  evdw = 0.d0
  eelt = 0.d0

  pgm_virial = 0.0d0
  vdw_virial = 0.0d0

  ! allocate for reciprocal jobs
  if ( eedmeth .ne. 4 ) call pGM_RECIP_allocate(numatoms)

  ! update the imaged crds
  call map_coords(crd,numatoms,recip)
  call adjust_imagcrds(crd,numatoms)

  ! update the global moments
  call pGM_MPOLE_local_to_global(crd)

  ! do the induction iteration here
  call pGM_INDUCED_eval(numatoms,crd,x,ipairs,diprms,dipiter)

  ! compute final energy, forces, and virial
  call pGM_NonBond_ene_frc(numatoms,crd,x,nummols,molsiz,amass,ipairs,&
           ntypes,iac,ico,cn1,cn2,eelt,evdw,frc,pgm_virial,vdw_virial)

  if ( ntb > 0 ) then
    ! the factor 0.5 is due to the fact that sander uses molecular kinetic energy as the
    ! kinetic part of virial, which is actually half of the real virial
    pgm_virial = 0.5d0 * coulomb_const_kcal_per_mole * pgm_virial
    vdw_virial = 0.5d0 * vdw_virial

    sander_vir(1) = sander_vir(1) + pgm_virial(1,1) + vdw_virial(1,1)
    sander_vir(2) = sander_vir(2) + pgm_virial(2,2) + vdw_virial(2,2)
    sander_vir(3) = sander_vir(3) + pgm_virial(3,3) + vdw_virial(3,3)
    sander_vir(4) = 0.0d0
    ! RL: the total virial will be computed by runmd() after reducing
    !sander_vir(4) = sander_vir(4) + pgm_virial(1,1) + pgm_virial(2,2) + pgm_virial(3,3) +&
    !                                vdw_virial(1,1) + vdw_virial(2,2) + vdw_virial(3,3)

    if ( master .and. pol_gauss_verbose == 2 )then
      write(6,'(a,3(1x,g16.8))') &
                ' pGM nonbond vir = ',pgm_virial(1,1),pgm_virial(1,2),pgm_virial(1,3)
      write(6,'(a,3(1x,g16.8))') &
                ' pGM nonbond vir = ',pgm_virial(2,1),pgm_virial(2,2),pgm_virial(2,3)
      write(6,'(a,3(1x,g16.8))') &
                ' pGM nonbond vir = ',pgm_virial(3,1),pgm_virial(3,2),pgm_virial(3,3)
    end if
  end if

end subroutine pGM_NonBond_eval
!-------------------------------------------------------------------------------
end module pol_gauss_interface
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_perm_fields(numatoms,is_polarizable,crd,x,ipairs, &
                                   cart_dipole_field)
  use pol_gauss_direct, only : pGM_DIRECT_perm_field
  use pol_gauss_recip, only : pGM_RECIP_perm_field
  use pol_gauss_self, only : pGM_SELF_perm_field
  implicit none

  integer,intent(in) :: numatoms
  logical,intent(in) :: is_polarizable(*)
  _REAL_,intent(in) :: crd(3,*)
  _REAL_,intent(in) :: x(*)
  integer, intent(in) :: ipairs(*)
  _REAL_,intent(inout) :: cart_dipole_field(3,*)

# include "ew_erfc_spline.h"
#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
# include "ew_parallel.h"
#endif

#if 0
integer :: i
#endif

#ifdef MPI
  integer :: commsander_numtasks, commsander_mytaskid, ier

  commsander_numtasks = numtasks
  commsander_mytaskid = mytaskid
#endif

  ! All computed fields are actually grad phi
  call zero_array(cart_dipole_field,3*numatoms)

#ifdef MPI
  if ( i_do_recip )then
    call mpi_comm_size(recip_comm,numtasks,ier)
    call mpi_comm_rank(recip_comm,mytaskid,ier)
#endif

  if (eedmeth.ne.4) call pGM_RECIP_perm_field(numatoms,crd,cart_dipole_field,x)

#ifdef MPI
    numtasks = commsander_numtasks
    mytaskid = commsander_mytaskid
  end if
#endif

#ifdef MPI
  if ( i_do_direct ) then
    call mpi_comm_size(direct_comm,numtasks,ier)
    call mpi_comm_rank(direct_comm,mytaskid,ier)
#endif

    call pGM_DIRECT_perm_field(ipairs,x,cart_dipole_field)

#ifdef MPI
    numtasks = commsander_numtasks
    mytaskid = commsander_mytaskid
  end if
#endif

  if (eedmeth.ne.4) call pGM_SELF_perm_field(numatoms,cart_dipole_field)

#ifdef MPI
  call MPI_ALLREDUCE(MPI_IN_PLACE,cart_dipole_field,3*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
#endif


end subroutine pGM_NonBond_perm_fields
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_dip_dip_fields_short(numatoms,x,ind_dip,dip_field)
  use pol_gauss_direct, only : pGM_DIRECT_dip_dip_field_short
  implicit none

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(out) :: dip_field(3,*)

#ifdef MPI
# include "parallel.h"
# include "ew_parallel.h"
#endif

#ifdef MPI
  integer :: commsander_numtasks, commsander_mytaskid, ier

  commsander_numtasks = numtasks
  commsander_mytaskid = mytaskid
  
  if ( i_do_direct ) then
    call mpi_comm_size(direct_comm,numtasks,ier)
    call mpi_comm_rank(direct_comm,mytaskid,ier)
#endif

    call zero_array(dip_field,3*numatoms)

    ! All computed fields are actually grad phi
    call pGM_DIRECT_dip_dip_field_short(ind_dip,dip_field)

#ifdef MPI
    numtasks = commsander_numtasks
    mytaskid = commsander_mytaskid
  end if
#endif

end subroutine pGM_NonBond_dip_dip_fields_short
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_dip_dip_fields(numatoms,x,ind_dip,dip_field)
  use pol_gauss_recip, only : pGM_RECIP_dipole_field
  use pol_gauss_direct, only : pGM_DIRECT_dip_dip_field
  use pol_gauss_self, only : pGM_SELF_dipole_field
  implicit none

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: x(*)
  _REAL_,intent(in) :: ind_dip(3,*)
  _REAL_,intent(out) :: dip_field(3,*)

# include "ew_erfc_spline.h"
#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
# include "ew_parallel.h"
#endif

#ifdef MPI
  integer :: commsander_numtasks, commsander_mytaskid, ier

  commsander_numtasks = numtasks
  commsander_mytaskid = mytaskid
#endif

  call zero_array(dip_field,3*numatoms)
  
#ifdef MPI
  if ( i_do_recip )then
    call mpi_comm_size(recip_comm,numtasks,ier)
    call mpi_comm_rank(recip_comm,mytaskid,ier)
#endif

    if (eedmeth.ne.4) call pGM_RECIP_dipole_field(numatoms,x,ind_dip,dip_field)

#ifdef MPI
    numtasks = commsander_numtasks
    mytaskid = commsander_mytaskid
  end if
#endif

#ifdef MPI
  if ( i_do_direct ) then
    call mpi_comm_size(direct_comm,numtasks,ier)
    call mpi_comm_rank(direct_comm,mytaskid,ier)
#endif

    call pGM_DIRECT_dip_dip_field(ind_dip,dip_field)

#ifdef MPI
    numtasks = commsander_numtasks
    mytaskid = commsander_mytaskid
  end if
#endif

  if (eedmeth.ne.4) call pGM_SELF_dipole_field(numatoms,ind_dip,dip_field)

#ifdef MPI
  call MPI_ALLREDUCE(MPI_IN_PLACE,dip_field,3*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
#endif

end subroutine pGM_NonBond_dip_dip_fields
!-------------------------------------------------------------------------------
subroutine pGM_NonBond_ene_frc(numatoms,crd,x,nummols,molsiz,amass,ipairs, &
                       ntypes,iac,ico,cn1,cn2, &
                       ene_elec,ene_vdw,frc,pgm_virial,vdw_virial)
  use pol_gauss_recip, only : pGM_RECIP_ene_frc
  use pol_gauss_direct, only : pGM_DIRECT_ene_frc
  use pol_gauss_self, only : pGM_SELF_ene_frc
  use pol_gauss_induced, only : ind_dip
  use pol_gauss_mdin, only : pol_gauss_verbose,dipole_print,saved_nprint
  use pol_gauss_multipoles, only : global_multipole,coulomb_const_kcal_per_mole

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*),x(*),amass(*)
  integer nummols,molsiz(*)
  integer,intent(in) :: ipairs(*)
  integer,intent(in) :: ntypes,iac(*),ico(*)
  _REAL_,intent(in) :: cn1(*),cn2(*)
  _REAL_,intent(out) :: ene_elec,ene_vdw
  _REAL_,intent(inout) :: frc(3,numatoms),pgm_virial(3,3),vdw_virial(3,3)

# include "box.h"
# include "ew_erfc_spline.h"
# include "extra.h"
#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
# include "ew_parallel.h"
#endif

  _REAL_ :: e_rec_perm,e_rec_ind,e_dir_perm,e_dir_ind,e_adj_perm, &
            e_adj_ind,e_self_perm,e_self_ind,e_dir_vdw,e_adj_vdw,e_rec_vdw
  _REAL_, allocatable :: pgm_force(:,:)
  _REAL_, allocatable :: phi(:,:), phi_rec(:,:), phi_dir(:,:), phi_self(:,:)
  _REAL_, allocatable :: displacement(:,:)
  _REAL_ :: mass, massmol, xmol, ymol, zmol
  integer :: imol, num, isave, iatom, i
  _REAL_ :: xmolt, ymolt, zmolt, mom, momt
  _REAL_ :: xmolp, ymolp, zmolp, momp

  ! electric potentials and derivatives
  allocate(pgm_force(3,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_rec(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_dir(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(phi_self(10,numatoms), stat=ier)
  REQUIRE(ier==0)
  allocate(displacement(3,numatoms), stat=ier)
  REQUIRE(ier==0)

  ene_ind = 0.0d0
  ene_vdw = 0.0d0
  pgm_force = 0.0d0

  phi = 0.0d0
  phi_rec = 0.0d0
  phi_dir = 0.0d0
  phi_self = 0.0d0

#ifdef MPI
  call mpi_comm_size(recip_comm,numtasks,ier)
  call mpi_comm_rank(recip_comm,mytaskid,ier)
  master = mytaskid.eq.0
#endif

  if (eedmeth .ne. 4) call pGM_RECIP_ene_frc(numatoms,crd,x,ind_dip,pgm_virial,phi_rec)

#ifdef MPI
  call mpi_comm_rank(commsander,mytaskid,ier)
  call mpi_comm_size(commsander,numtasks,ier)
  master = mytaskid.eq.0
#endif

  ! get ready for virial calculation
  if ( ntb > 0 ) then
    i = 0
    do imol = 1,nummols
      massmol = 0.d0
      xmol = 0.d0
      ymol = 0.d0
      zmol = 0.d0
      num = molsiz(imol)
      isave = i
     
      ! get c.o.m. of molecule - TODO: use precalced masses
     
      do iatom = 1, num
        i = i + 1
        mass = amass(i)
        massmol = massmol + mass
        xmol = xmol + mass*crd(1,i)
        ymol = ymol + mass*crd(2,i)
        zmol = zmol + mass*crd(3,i)
      end do
      xmol = xmol / massmol
      ymol = ymol / massmol
      zmol = zmol / massmol

      i = isave
      do iatom = 1,num
        i = i + 1
        displacement(1,i) = crd(1,i) - xmol
        displacement(2,i) = crd(2,i) - ymol
        displacement(3,i) = crd(3,i) - zmol
      end do
    end do
  end if

  call pGM_DIRECT_ene_frc(ipairs,ntypes,iac,ico,cn1,cn2,crd,x,ind_dip, &
           ene_vdw,frc,pgm_virial,vdw_virial,phi_dir,numatoms,displacement)

  if (eedmeth .ne. 4) call pGM_SELF_ene_frc(numatoms,ind_dip,phi_self)

  phi = phi_rec + phi_dir + phi_self

  ! Now get full copy of phi's on all PEs 
  ! Haixin: To compute virial, also need other phi's, so do it here too
#ifdef MPI
  call mpi_barrier(commsander, ier)
  call mpi_allreduce(MPI_IN_PLACE,phi,10*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
  if ( ntb > 0 ) then
    !call mpi_allreduce(MPI_IN_PLACE,pgm_virial,9,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
    if ( eedmeth .ne. 4 ) then
    call mpi_allreduce(MPI_IN_PLACE,phi_rec,10*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
    ! Haixin: This is for the commented out virial procedure
    !call mpi_allreduce(MPI_IN_PLACE,phi_self,10*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
    end if
    !call mpi_allreduce(MPI_IN_PLACE,phi_dir,10*numatoms,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ier)
  end if
#endif

  call pGM_energy_force_virial(numatoms,crd,displacement,ene_elec,pgm_virial,&
               pgm_force,ind_dip,phi,phi_rec,phi_dir,phi_self)

  ene_elec = coulomb_const_kcal_per_mole * ene_elec
  frc = frc + coulomb_const_kcal_per_mole * pgm_force

  ! if requesting dipole printings, do it here
  saved_nprint = saved_nprint + 1
  if ( master .and. dipole_print /= 0 ) then
  if ( mod(saved_nprint,dipole_print) == 0 ) then
    write(100,*) ' pGM Current dipole moments @ Step # ', saved_nprint
    i = 0
    do imol = 1, nummols
      xmol = 0.d0; ymol = 0.d0; zmol = 0.d0
      xmolp = 0.d0; ymolp = 0.d0; zmolp = 0.d0
      xmolt = 0.d0; ymolt = 0.d0; zmolt = 0.d0
      num = molsiz(imol)
      do iatom = 1, num
        i = i + 1
        xmol = xmol + ind_dip(1,i)
        ymol = ymol + ind_dip(2,i)
        zmol = zmol + ind_dip(3,i)
        xmolp = xmolp +                global_multipole(2,i) + global_multipole(1,i) * crd(1,i)
        ymolp = ymolp +                global_multipole(3,i) + global_multipole(1,i) * crd(2,i)
        zmolp = zmolp +                global_multipole(4,i) + global_multipole(1,i) * crd(3,i)
        xmolt = xmolt + ind_dip(1,i) + global_multipole(2,i) + global_multipole(1,i) * crd(1,i)
        ymolt = ymolt + ind_dip(2,i) + global_multipole(3,i) + global_multipole(1,i) * crd(2,i)
        zmolt = zmolt + ind_dip(3,i) + global_multipole(4,i) + global_multipole(1,i) * crd(3,i)
      end do
      mom  = sqrt(xmol**2  + ymol**2  + zmol**2)
      momp = sqrt(xmolp**2 + ymolp**2 + zmolp**2)
      momt = sqrt(xmolt**2 + ymolt**2 + zmolt**2)
      write(100,*) ' Water #, induced, permanent, and total moments', imol, mom, momp, momt
    end do
    write(200,*) ' pGM Current induced dipole vectors @ Step # ', saved_nprint
    do iatom = 1, numatoms
      write(200,*) ' Atom #, induced moment', iatom, ind_dip(1:3, iatom)
    end do
  end if
  end if
  if ( master .and. pol_gauss_verbose == 2 )then
    write(6,'(a,/,4(1x,f14.4))') &
            ' e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm = ', &
              e_rec_perm,e_dir_perm,e_adj_perm,e_self_perm
    write(6,'(a,/,4(1x,f14.4))') &
            ' e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind = ', &
              e_rec_ind,e_dir_ind,e_adj_ind,e_self_ind
    write(6,'(a,/,3(1x,f14.4))') &
            ' e_dir_vdw,e_adj_vdw,e_rec_vdw = ', &
             e_dir_vdw,e_adj_vdw,e_rec_vdw
  endif

  deallocate(pgm_force, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_rec, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_dir, stat=ier)
  REQUIRE(ier==0)
  deallocate(phi_self, stat=ier)
  REQUIRE(ier==0)
  deallocate(displacement, stat=ier)
  REQUIRE(ier==0)

end subroutine pGM_NonBond_ene_frc
!-------------------------------------------------------------------------------
subroutine pGM_energy_force_virial(numatoms,crd,displacement,pgm_energy,pgm_virial,&
               pgm_force,induced_dipole,phi,phi_rec,phi_dir,phi_self)
  use pol_gauss_multipoles, only : global_multipole,covalent_ptr,covalent_atm,covalent_dip,&
          start_multipoles,end_multipoles
  use nblist, only: ucell,recip
  implicit none

  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,numatoms),displacement(3,numatoms)
  _REAL_,intent(out) :: pgm_energy,pgm_virial(3,3)
  _REAL_,intent(out) :: pgm_force(3,numatoms)
  _REAL_,intent(in) :: induced_dipole(3,numatoms)
  _REAL_,intent(inout) :: phi(10,numatoms),phi_rec(10,numatoms),phi_dir(10,numatoms),phi_self(10,numatoms)

# include "box.h"
# include "pol_gauss_mpole_index.h"
#ifdef MPI
  include 'mpif.h'
# include "extra.h"
# include "parallel.h"
# include "ew_parallel.h"
#endif

  integer :: i, j, k, atom_first, atom_last, covalent_atom
  _REAL_ :: distance, deltr_dot_field, deltr_dot_field2!, deltr_dot_dipole
  
  ! Haixin: This is for the commented out virial procedure
  !_REAL_ :: virial_force(3,numatoms)
  !virial_force = 0.0

  ! Haixin: convert grad of phi to electric field and its derivatives
  phi(2:10,:) = -phi(2:10,:)
  if ( ntb > 0 ) then
  phi_rec(2:10,:) = -phi_rec(2:10,:)
  ! Haixin: This is for the commented out virial procedure
  !phi_dir(2:10,:) = -phi_dir(2:10,:)
  !phi_self(2:10,:) = -phi_self(2:10,:)
  end if

  do i = start_multipoles, end_multipoles
    ! energy can be expressed as 1/2*q*phi and -1/2*mu*E
    pGM_energy = pGM_energy + 0.5d0*global_multipole(1,i)*phi(1,i)&
                            - 0.5d0*global_multipole(2,i)*phi(2,i)&
                            - 0.5d0*global_multipole(3,i)*phi(3,i)&
                            - 0.5d0*global_multipole(4,i)*phi(4,i)!&

    ! forces due to covalent monopoles
    pgm_force(1:3,i) = pgm_force(1:3,i) + global_multipole(1,i)*phi(2:4,i)

    ! forces due to total dipoles (covalent dipoles + induced dipoles)
    pgm_force(1,i) = pgm_force(1,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_200,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_110,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_101,i)
    pgm_force(2,i) = pgm_force(2,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_110,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_020,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_011,i)
    pgm_force(3,i) = pgm_force(3,i) + (global_multipole(2,i)+induced_dipole(1,i))*phi(Ind_101,i)&
                                    + (global_multipole(3,i)+induced_dipole(2,i))*phi(Ind_011,i)&
                                    + (global_multipole(4,i)+induced_dipole(3,i))*phi(Ind_002,i)

    ! forces due to changes in covalent dipoles 
    atom_last = covalent_ptr(i)
    atom_first = covalent_ptr(i-1) + 1
    do j = atom_first, atom_last
      covalent_atom = covalent_atm(j)

      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
      distance = sqrt(distance)
      if ( distance .gt. 3.0d0 ) then
        write(6,*) 'pGM Fatal Error: covalent bond too long!'
        call mexit(6,1)
      end if
      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi(2:4,i))
      do k = covalent_ptr(covalent_atom-1)+1, covalent_ptr(covalent_atom)
        if (covalent_atm(k) == i) exit
      end do
      deltr_dot_field2 = dot_product((crd(1:3,i) - crd(1:3,covalent_atom)),phi(2:4,covalent_atom))

      pgm_force(1,i) = pgm_force(1,i) - &
                       covalent_dip(j)*(&
                                                                phi(2,i)/distance - &
                       (crd(1,covalent_atom) - crd(1,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(2,covalent_atom)/distance - &
                       (crd(1,i) - crd(1,covalent_atom))*deltr_dot_field2/distance**3&
                       )
      pgm_force(2,i) = pgm_force(2,i) - &
                       covalent_dip(j)*(&
                                                                phi(3,i)/distance - &
                       (crd(2,covalent_atom) - crd(2,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(3,covalent_atom)/distance - &
                       (crd(2,i) - crd(2,covalent_atom))*deltr_dot_field2/distance**3&
                       )
      pgm_force(3,i) = pgm_force(3,i) - &
                       covalent_dip(j)*(&
                                                                phi(4,i)/distance - &
                       (crd(3,covalent_atom) - crd(3,i))*deltr_dot_field/distance**3&
                       ) + &
                       covalent_dip(k)*(&
                                                     phi(4,covalent_atom)/distance - &
                       (crd(3,i) - crd(3,covalent_atom))*deltr_dot_field2/distance**3&
                       )

    end do
  end do

  if ( ntb > 0 ) then
!  Haixin: The alternative virial procedure
!  phi = phi_dir + phi_self ! now the phi array is contaminated, this is only for virial calculation
!
!  do i = 1, numatoms
!
!    atom_last = covalent_ptr(i)
!    atom_first = covalent_ptr(i-1) + 1
!    do j = atom_first, atom_last
!      covalent_atom = covalent_atm(j)
!      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
!      distance = sqrt(distance)
!      if ( distance .gt. 3.0d0 ) then
!        write(6,*) 'pGM Fatal Error: covalent bond too long!'
!        call mexit(6,1)
!      endif
!      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi(2:4,i))
!      do k = covalent_ptr(covalent_atom-1)+1, covalent_ptr(covalent_atom)
!        if (covalent_atm(k) == i) exit
!      end do
!      deltr_dot_field2 = dot_product((crd(1:3,i) - crd(1:3,covalent_atom)),phi(2:4,covalent_atom))
!
!      virial_force(1,i) = virial_force(1,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(2,i)/distance - &
!                       (crd(1,covalent_atom) - crd(1,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(2,covalent_atom)/distance - &
!                       (crd(1,i) - crd(1,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!      virial_force(2,i) = virial_force(2,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(3,i)/distance - &
!                       (crd(2,covalent_atom) - crd(2,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(3,covalent_atom)/distance - &
!                       (crd(2,i) - crd(2,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!      virial_force(3,i) = virial_force(3,i) - &
!                       covalent_dip(j)*(&
!                                                                phi(4,i)/distance - &
!                       (crd(3,covalent_atom) - crd(3,i))*deltr_dot_field/distance**3&
!                       ) + &
!                       covalent_dip(k)*(&
!                                                     phi(4,covalent_atom)/distance - &
!                       (crd(3,i) - crd(3,covalent_atom))*deltr_dot_field2/distance**3&
!                       )
!
!    end do
!  end do
!
!  do i = 1, numatoms
!
!    pgm_virial(1,1) = pgm_virial(1,1) - &
!                      (virial_force(1,i)*crd(1,i) - phi_rec(2,i)*induced_dipole(1,i))
!    pgm_virial(1,2) = pgm_virial(1,2) - &
!                      (virial_force(1,i)*crd(2,i) - phi_rec(2,i)*induced_dipole(2,i))
!    pgm_virial(1,3) = pgm_virial(1,3) - &
!                      (virial_force(1,i)*crd(3,i) - phi_rec(2,i)*induced_dipole(3,i))
!    pgm_virial(2,1) = pgm_virial(2,1) - &
!                      (virial_force(2,i)*crd(1,i) - phi_rec(3,i)*induced_dipole(1,i))
!    pgm_virial(2,2) = pgm_virial(2,2) - &
!                      (virial_force(2,i)*crd(2,i) - phi_rec(3,i)*induced_dipole(2,i))
!    pgm_virial(2,3) = pgm_virial(2,3) - &
!                      (virial_force(2,i)*crd(3,i) - phi_rec(3,i)*induced_dipole(3,i))
!    pgm_virial(3,1) = pgm_virial(3,1) - &
!                      (virial_force(3,i)*crd(1,i) - phi_rec(4,i)*induced_dipole(1,i))
!    pgm_virial(3,2) = pgm_virial(3,2) - &
!                      (virial_force(3,i)*crd(2,i) - phi_rec(4,i)*induced_dipole(2,i))
!    pgm_virial(3,3) = pgm_virial(3,3) - &
!                      (virial_force(3,i)*crd(3,i) - phi_rec(4,i)*induced_dipole(3,i))
!    ! this is a different way of collecting part of direct and self related
!    ! virial, different from below which is one obvious way of implimentation of
!    ! the formula in the paper. If interested, it can be a good exercise to
!    ! prove this two ways of collecting are equal.
!
!    atom_last = covalent_ptr(i)
!    atom_first = covalent_ptr(i-1) + 1
!    do j = atom_first, atom_last
!      covalent_atom = covalent_atm(j)
!      distance = dot_product((crd(1:3,covalent_atom)-crd(1:3,i)),(crd(1:3,covalent_atom)-crd(1:3,i)))
!      distance = sqrt(distance)
!      if ( distance .gt. 3.0d0 ) then
!        write(6,*) 'pGM Fatal Error: covalent bond too long!'
!        call mexit(6,1)
!      endif
!      deltr_dot_field = dot_product((crd(1:3,covalent_atom) - crd(1:3,i)),phi_rec(2:4,i))
!
!      pgm_virial(1,1) = pgm_virial(1,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(1,2) = pgm_virial(1,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(1,3) = pgm_virial(1,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(1,covalent_atom)-crd(1,i))*(crd(3,covalent_atom)-crd(3,i))
!      pgm_virial(2,1) = pgm_virial(2,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(2,2) = pgm_virial(2,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(2,3) = pgm_virial(2,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(2,covalent_atom)-crd(2,i))*(crd(3,covalent_atom)-crd(3,i))
!      pgm_virial(3,1) = pgm_virial(3,1) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(1,covalent_atom)-crd(1,i))
!      pgm_virial(3,2) = pgm_virial(3,2) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(2,covalent_atom)-crd(2,i))
!      pgm_virial(3,3) = pgm_virial(3,3) - (-deltr_dot_field*covalent_dip(j)/distance**3)*&
!                                          (crd(3,covalent_atom)-crd(3,i))*(crd(3,covalent_atom)-crd(3,i))
!    end do
!
!  end do

    do i = start_multipoles, end_multipoles
    !do i = 1,numatoms

       pgm_virial(1,1) = pgm_virial(1,1) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(1,i)&
                         )

       pgm_virial(1,2) = pgm_virial(1,2) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(2,i)&
                         )

       pgm_virial(1,3) = pgm_virial(1,3) - &
                         (-phi_rec(2,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_200,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_110,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_101,i)*displacement(3,i)&
                         )

       pgm_virial(2,1) = pgm_virial(2,1) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(1,i)&
                         )

       pgm_virial(2,2) = pgm_virial(2,2) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(2,i)&
                         )

       pgm_virial(2,3) = pgm_virial(2,3) - &
                         (-phi_rec(3,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_110,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_020,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_011,i)*displacement(3,i)&
                         )

       pgm_virial(3,1) = pgm_virial(3,1) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(1,i)+global_multipole(2,i)+global_multipole(1,i)*displacement(1,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(1,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(1,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(1,i)&
                         )

       pgm_virial(3,2) = pgm_virial(3,2) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(2,i)+global_multipole(3,i)+global_multipole(1,i)*displacement(2,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(2,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(2,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(2,i)&
                         )

       pgm_virial(3,3) = pgm_virial(3,3) - &
                         (-phi_rec(4,i)*&
                          (induced_dipole(3,i)+global_multipole(4,i)+global_multipole(1,i)*displacement(3,i))-&
                          (induced_dipole(1,i)+global_multipole(2,i))*phi_rec(Ind_101,i)*displacement(3,i)-&
                          (induced_dipole(2,i)+global_multipole(3,i))*phi_rec(Ind_011,i)*displacement(3,i)-&
                          (induced_dipole(3,i)+global_multipole(4,i))*phi_rec(Ind_002,i)*displacement(3,i)&
                         )

    end do
  end if

end subroutine pGM_energy_force_virial
!-------------------------------------------------------------------------------
subroutine pGM_get_numlist(header,nf,num_list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf
  integer, intent(out) :: num_list

  integer :: iok,ionerr
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  fmtin = '(10I8)'
  dtype = header//'NUM_LIST'
  ionerr = 1 ! not fatal if missing
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  if ( iok == 0 )then !this data type found in prmtop
    read(nf,fmt)num_list
  else !either old style prmtop or data not found
    num_list = 0 ! upon return this will invalidate valence_term
  endif
end subroutine pGM_get_numlist
!------------------------------------------------------------------------
subroutine pGM_read_list_data(header,nf,dim1,num_list,list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf,dim1,num_list
  integer, intent(out) :: list(dim1,num_list)

  integer :: iok,ionerr,j,k
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  ionerr = 0 !fatal if missing
  fmtin = '(10I8)'
  dtype = header//'LIST'
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  read(nf,fmt)((list(j,k),j=1,dim1),k=1,num_list)
end subroutine pGM_read_list_data
!----------------------------------------------------------
subroutine pGM_read_real_list_data(header,nf,dim1,num_list,list)
  implicit none

  character(len=*), intent(in) :: header
  integer, intent(in) :: nf,dim1,num_list
  _REAL_, intent(out) :: list(dim1,num_list)

  integer :: iok,ionerr,j,k
  character(len=80) :: fmt
  character(len=80) :: fmtin,dtype

  ionerr = 0 !fatal if missing
  fmtin = '(5E16.8)'
  dtype = header//'LIST'
  call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
  read(nf,fmt)((list(j,k),j=1,dim1),k=1,num_list)
end subroutine pGM_read_real_list_data
!--------------------------------------------------------------------
subroutine pGM_read_real_scalar(flag,nf,scalar_value)
  implicit none

  character(len=*), intent(in) :: flag
  integer, intent(in) :: nf
  _REAL_, intent(out) :: scalar_value

  integer :: iok,ionerr
  character(len=80) :: fmt
  character(len=80) :: fmtin

  ionerr = 0 !fatal if missing
  fmtin = '(E16.8)'
  call nxtsec(nf,  6,  ionerr,fmtin,  flag,  fmt,  iok)
  read(nf,fmt)scalar_value
end subroutine pGM_read_real_scalar
!--------------------------------------------------------------------
