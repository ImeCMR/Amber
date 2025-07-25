module xray_interface_impl_cpu_module
   use pmemd_lib_mod
   use file_io_mod
   use xray_globals_module
   use xray_contracts_module
   use xray_bulk_mask_data_module, only: k_sol, b_sol
   use xray_interface_pre_init_data, only: num_scatter_types => n_scatter_types
   implicit none
   private

   public :: finalize
   public :: init
   public :: xray_get_derivative
   public :: xray_read_mdin
   public :: xray_read_parm
   public :: xray_write_md_state
   public :: xray_write_min_state
   public :: xray_write_options

   namelist /xray/ &
         pdb_infile, pdb_outfile, &
         fave_outfile, fmtz_outfile, &
         pdb_read_coordinates, &
         pdb_use_segid, &
         pdb_wrap_names, &
         spacegroup_name, &
         reflection_infile, &
         xray_weight_initial, &
         xray_weight_final, &
         target, &
         solvent_mask_adjustment, &
         solvent_mask_probe_radius, &
         ntwsf, &
         sf_outfile, &
         atom_selection_mask, &
         k_sol, b_sol,  &
         mask_update_period, scale_update_period, &
         ml_update_period, bulk_solvent_model
   
   !-------------------------------------------------------------------
contains

   subroutine xray_read_mdin(mdin_lun)
      use gbl_constants_mod, only : error_hdr
      implicit none
      integer, intent(in) :: mdin_lun

      character(len=512) :: line
      integer :: stat
      
      call xray_init_globals()
      if (.not.xray_active) then
        return
      end if

      ! Trap to catch users trying to run &xray in parallel
#ifdef MPI
      write(stdout, '(A)') 'Running simulations with an &xray namelist requires a serial &
                           &installation.'
      write(stdout, '(A)') 'Use the GPU extensions for best performance.'
      call mexit(stdout,1)
#endif
      rewind(mdin_lun)
      read(unit=mdin_lun,nml=xray,iostat=stat)

      if (stat /= 0) then
        backspace(mdin_lun)
        read(mdin_lun, fmt='(A)') line
        write(stdout, '(A)') 'Invalid line in &xray namelist '//': '//trim(line)
        call mexit(stdout,1)
      end if

      ! Check input
      if (ntwsf < 0) then
        write(stdout, '(A,A)') error_hdr, 'ntwsf must be >= 0'
        call mexit(stdout,1)
      end if
      
      if (xray_weight_initial == sentinel_xray_weight) then
         call check_requirement(xray_weight_final == sentinel_xray_weight, &
             & "`xray_weight_final` requires `xray_weight_initial` to be explicitly set")
         xray_weight_initial = default_xray_weight
      end if
   
      if (xray_weight_final == sentinel_xray_weight) then
         xray_weight_final = xray_weight_initial
      end if
      
      !write(unit=6,nml=xray)
   end subroutine xray_read_mdin

   subroutine xray_write_options()
      use xray_bulk_mask_data_module, only: k_sol, b_sol
      implicit none

      write(stdout,'(/,A)') 'X-ray Refinement Parameters:'
      write(stdout,'(5X,2A)') 'PDB InFile: ',trim(pdb_infile)
      if( pdb_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'PDB OutFile:',trim(pdb_outfile)
      if( fave_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'FCALC_AVE OutFile:',trim(fave_outfile)
      if( fmtz_outfile /= '' ) &
         write(stdout,'(5X,2A)') 'FMTZ OutFile:',trim(fmtz_outfile)
      write(stdout,'(5X,A,L1)') 'PDB Read Coordinates: ',pdb_read_coordinates
      write(stdout,'(5X,A,L1)') 'PDB Use SegID: ',pdb_use_segid
      write(stdout,'(5X,A,L1)') 'PDB Wrap Names: ',pdb_wrap_names
      write(stdout,'(5X,2A)') 'Spacegroup: ',trim(spacegroup_name)
      write(stdout,'(5X,2A)') 'Reflection InFile: ',trim(reflection_infile)
      write(stdout,'(5X,A,E10.3,A,E10.3)') 'X-ray weights: ', xray_weight_initial, ' ... ', xray_weight_final
      write(stdout,'(5X,A,A4)') 'Use target: ',target
      write(stdout,'(5X,A,I5)') 'Scale update Interval: ',scale_update_period
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask adjustment: ',solvent_mask_adjustment
      ! write(stdout,'(5X,A,F8.3)') 'Solvent mask probe radius: ',solvent_mask_probe_radius
      ! write(stdout,'(5X,2A)') 'Solvent Mask OutFile:',trim(solvent_mask_outfile)
      ! write(stdout,'(5X,2A)') 'Solvent Mask Reflection OutFile:',trim(solvent_mask_reflection_outfile)
      write(stdout,'(5X,A,I5)') 'Solvent Mask Update Interval: ',mask_update_period
      write(stdout,'(5X,2(A,F8.3))') 'Solvent scale:',k_sol,', B-factor:', b_sol
      ! write(stdout,'(5X,2(A,F8.3))') 'B-Factor Min:',bfactor_min,', Max: ',bfactor_max
      ! write(stdout,'(5X,A,I4)') 'B-factor Refinement Interval: ',bfactor_refinement_interval
      write(stdout,'(5X,2A)') 'Atom Selection Mask: ',trim(atom_selection_mask)
      return
   end subroutine xray_write_options

   ! Read X-ray data from the PRMTOP file, and also save pointers to global
   ! PRMTOP data.
   subroutine xray_read_parm(prmtop_lun, out_lun)
      use nextprmtop_section_mod, only: nxtsec, nxtsec_reset
      use prmtop_dat_mod, only: natom, nres
      use xray_interface_pre_init_data, only: n_scatter_coeffs, scatter_coefficients
      implicit none
      integer, intent(in) :: prmtop_lun, out_lun
      ! local
      character(len=32) :: fmt
      integer :: ierr
      logical :: master=.true.

      ! if (pdb_outfile /= '') then
         num_atoms = natom
         num_residues = nres

         allocate(atom_bfactor(natom), atom_occupancy(natom), &
               atom_selection(natom), residue_chainid(nres), residue_icode(nres), &
               atom_element(natom), atom_altloc(natom), residue_number(nres))

         call nxtsec_reset()
         call nxtsec(prmtop_lun,STDOUT,0,'(20I4)','RESIDUE_NUMBER',fmt,ierr)
         read(prmtop_lun,fmt) residue_number
         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','RESIDUE_CHAINID',fmt,ierr)
         read(prmtop_lun,fmt) residue_chainid

         call nxtsec(prmtop_lun,STDOUT,1,'*','RESIDUE_ICODE',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) residue_icode
         else
            residue_icode=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,1,'*','ATOM_ALTLOC',fmt,ierr)
         if (ierr==0) then
            read(prmtop_lun,fmt) atom_altloc
         else
            atom_altloc=' '
         end if

         call nxtsec(prmtop_lun,STDOUT,0,'(20A4)','ATOM_ELEMENT',fmt,ierr)
         read(prmtop_lun,fmt) atom_element
      ! end if

      if (reflection_infile == '') return

      call nxtsec(prmtop_lun,out_lun,0,'(I4)','XRAY_NUM_SCATTER_TYPES',fmt,ierr)
      if (fmt=='*') then
         write(stdout,'(A)') &
            'ERROR: XRAY_NUM_SCATTER_TYPES not found in PRMTOP file.'
         call mexit(stdout,1)
      end if
      read(prmtop_lun,fmt) num_scatter_types

      allocate(atom_scatter_type(natom),  &
            scatter_coefficients(2, n_scatter_coeffs, num_scatter_types))

      call nxtsec(prmtop_lun,out_lun,0,'(20I4)','XRAY_ATOM_SCATTER_TYPE_INDEX',fmt,ierr)
      read(prmtop_lun,fmt) atom_scatter_type
      call nxtsec(prmtop_lun,out_lun,0,'(F12.6)','XRAY_SCATTER_COEFFICIENTS',fmt,ierr)
      read(prmtop_lun,fmt) scatter_coefficients
      call nxtsec(prmtop_lun,out_lun,1,'*','XRAY_SYMMETRY_TYPE',fmt,ierr)
      if (ierr==-2) then
         if( master ) write(STDOUT,*) &
               'XRAY_SYMMETRY_TYPE not found in PRMTOP file; assuming P1'
         num_symmops = 1
         spacegroup_number = 1
         spacegroup_name = 'P 1'
         au_type = 1
      else
         stop 'ONLY P1 SUPPORTED FOR NOW'
         read(prmtop_lun,fmt) num_symmops, spacegroup_number, au_type, system
         call nxtsec(prmtop_lun,out_lun,1,'*', &
               'XRAY_SYMMETRY_OPERATORS',fmt,ierr)
         ! ...
      end if
   end subroutine xray_read_parm

   subroutine xray_read_pdb(filename)
      use prmtop_dat_mod, only: gbl_res_atms,gbl_labres,atm_igraph
      use inpcrd_dat_mod, only: atm_crd
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      character(len=4) :: name,resName,segID,element,altLoc,chainID,iCode
      integer :: serial,resSeq
      real(real_kind) :: xyz(3),occupancy,tempFactor
      character(len=80) :: line
      integer :: unit, iostat, iatom, ires, i, j, ndup, nmiss
      real(real_kind), parameter :: MISSING = -999.0_rk_
      logical :: master=.true.
      ! begin
      atom_occupancy(:)=MISSING
      call amopen(allocate_lun(unit),filename,'O','F','R')
      ndup=0
      iatom=1
      ires=1
      do
         read(unit,'(A)',iostat=iostat) line
         if (iostat/=0) exit
         if (line(1:6)=='END   ') exit
         if (line(1:6)=='ATOM  ' .or. line(1:6)=='HETATM') then
            read(line,'(6X,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)') &
                  serial,name,altLoc,resName,chainID,resSeq,iCode, &
                  xyz,occupancy,tempFactor,segID,element
            i = find_atom(name,resName,chainID,resSeq,iCode)
            if (i<0) then
               write(stdout,'(A)') 'Atom not found:'
               write(stdout,'(A)') trim(line)
               stop
            end if
            if (atom_occupancy(i) >= 0) then
               ndup=ndup+1
               if (ndup<10) then
                  if( master ) write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Duplicate ATOM:', &
                        name,resName,chainID(1:1),resSeq,iCode(1:1)
               end if
            end if
            if (pdb_read_coordinates) atm_crd(1:3,i) = xyz
            atom_bfactor(i) = tempFactor
            atom_occupancy(i) = occupancy
         end if
      end do
      nmiss = count(atom_occupancy==MISSING)
      if (nmiss>0) then
         if( master ) write(stdout,'(A,I4,A)') 'PDB: missing data for ',nmiss,' atoms.'
         j=0
         do i=1,num_atoms
            if (atom_occupancy(i)==MISSING) then
               atom_occupancy(i)=0
               ires = residue_number(i)
               j=j+1
               if (j<=10) then
                  if( master ) then
                     write(stdout,'(3(A,1X),A,I4,A)') 'PDB: Missing ATOM:', &
                             atm_igraph(i),gbl_labres(ires),residue_chainID(ires)(1:1),&
                             ires,residue_iCode(ires)(1:1)
                  end if
               end if
            end if
         end do
      end if
      if (nmiss==0 .and. ndup==0) then
         if( master ) write(stdout,'(A)') 'PDB: All atoms read successfully.'
      end if
      close(unit)
      return
   end subroutine xray_read_pdb

   function find_atom(name,resName,chainID,resSeq,iCode) result(atom_serial)
      use prmtop_dat_mod, only: gbl_res_atms,gbl_labres,atm_igraph
      implicit none
      integer :: atom_serial
      character(len=4), intent(in) :: name, resName, chainID, iCode
      integer, intent(in) :: resSeq
      ! locals
      character(len=4) :: lname
      integer, save :: ires = 1
      integer :: i,j
      lname = adjustl(name)
      ! first find the matching residue:
      do i=1,num_residues
         if (resSeq==residue_number(ires) &
               .and. chainID==residue_chainid(ires) &
               .and. iCode==residue_icode(ires) &
               .and. resName==gbl_labres(ires)) then
            ! then find the matching atom name:
            do j = gbl_res_atms(ires),gbl_res_atms(ires+1)-1
               if (lname==atm_igraph(j)) then
                  atom_serial = j
                  return
               end if
            end do
            ! Continue searching, just in case there is a residue
            ! that has been split into two parts.
         end if
         ires = ires + 1
      end do
      atom_serial = -1
      return
   end function find_atom

   subroutine xray_write_pdb(filename)
      use file_io_dat_mod, only: owrite
      use inpcrd_dat_mod, only: atm_crd
      use prmtop_dat_mod, only: &
            gbl_res_atms,gbl_labres,atm_igraph, &
            num_bonds=>nbona
      implicit none
      character(len=*), intent(in) :: filename
      ! locals
      integer :: unit, iatom, ires, ierr
      integer :: first1, last1, first2, last2, ibond
      integer :: iatom_p, ires_p
      character(len=4) :: name
      character(len=8) :: date
      character(len=10) :: time
      ! character(len=1) :: altloc
      ! character(len=4) :: segid
      character(len=3) :: resName
      logical          :: isStandardRes
      character(len=3), parameter :: standard_pdb_residues(28) = (/ &
            "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE", &
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL", &
            " DG"," DA"," DT"," DC","  G","  A","  U","  C" /)
      ! Amber modres types: "CYX","HID","HIE","HIP",
      character(len=*), parameter :: pdbfmt_MODRES = &
            '("MODRES",1X,A4,1X,A3,1X,A1,1X,I4,A1,1X,A3,2X,A41)'
      ! GMS: Fix for pgf90 compiler
      character(len=4) :: this_residue_chainid

      call amopen(allocate_lun(unit),filename,owrite,'F','R')
      call date_and_time(date,time)
      ! if (title/='') write(unit,'(2A)') 'REMARK  ', title
      ! if (title1/='') write(unit,'(2A)') 'REMARK  ', title1
      write(unit,'(12A)') 'REMARK  Written by Amber 24, PMEMD, ', &
            date(1:4),'.',date(5:6),'.',date(7:8),'  ', &
            time(1:2),':',time(3:4),':',time(5:6)

      ! Actually '(6A,3F9.3A9,3F7.2,1X,A11,I4)', with last value = Z
      write(unit,'(A6,3F9.3,3F7.2,1X,A11)') &
            'CRYST1', unit_cell%as_array(), spacegroup_name
#if 0
      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         if (gbl_labres(ires)=='HID') then
            write(unit,pdbfmt_MODRES) &
                  '----','HID',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HE2 ATOM REMOVED'
         else if (gbl_labres(ires)=='HIE') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIE',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 ATOM REMOVED'
         else if (gbl_labres(ires)=='HIP') then
            write(unit,pdbfmt_MODRES) &
                  '----','HIP',
            residue_chainid(ires)(1:1), &
                  residue_number(ires),residue_icode(ires)(1:1), &
                  'HIS','HD1 AND HE2 ATOMS REMOVED'
         end if
      end do
#endif

      do ires = 1,num_residues
         if (residue_chainid(ires)=='*') cycle
         do iatom = gbl_res_atms(ires), gbl_res_atms(ires+1)-1
            ! ***NOTE***
            ! This code only adds a leading space to give element-alignment
            ! where possible. It is impossible to follow the PDB version 3
            ! "remediated" format correctly, because it has no alignment rules.
            ! Instead, it assumes you have a complete database of all known
            ! residues, and any other residue names are a fatal error.
            name = atm_igraph(iatom)
            if (atom_element(iatom)(1:1)==' ' &
                  .and. name(1:1) == atom_element(iatom)(2:2)) then
               if (len_trim(name) < 4 .or. pdb_wrap_names) then
                  name = name(4:4)//name(1:3)
               end if
            end if
            resName=gbl_labres(ires)(1:3)
            resName=adjustr(resName)
            ! GMS: Fix for pgf90 compiler
            this_residue_chainid = residue_chainid(ires)
            ! DRR: PGI does not seem to like any() inside merge() intrinsic.
            isStandardRes = any(resName==standard_pdb_residues)
            ! don't overflow atom or residue numbers:
            iatom_p = mod( iatom, 100000 )
            ires_p = mod( residue_number(ires), 10000 )
            write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,6X,2A4)')&
                  merge('ATOM  ', 'HETATM', isStandardRes), &
                  iatom_p,name,atom_altloc(iatom)(1:1), &
                  resName,residue_chainid(ires)(1:1), &
                  ires_p,residue_icode(ires)(1:1), &
                  atm_crd(1:3,iatom), &
                  atom_occupancy(iatom), &
                  atom_bfactor(iatom), &
                  merge(this_residue_chainid,'    ',pdb_use_segid), &
                  atom_element(iatom)
         end do
      end do
      write(unit,'(A)') 'END'
      close(unit)
      return
   end subroutine xray_write_pdb

   subroutine init()

      use pbc_mod, only: pbc_box, pbc_alpha, pbc_beta, pbc_gamma
      use xray_pure_utils, only: pack_index
      use xray_unit_cell_module, only: unit_cell_t
      use findmask_mod, only: atommask
      use prmtop_dat_mod, only: natom,nres,atm_igraph,gbl_res_atms, gbl_labres, atm_isymbl, atm_atomicnumber
      use inpcrd_dat_mod, only: atm_crd
      use xray_interface_pre_init_data, only: scatter_coefficients
      use xray_interface2_module, only: init_interface2 => init
      use xray_debug_dump_module, only: xray_dump => dump  ! FIXME: remove this line in release
      implicit none
      ! local
      integer :: hkl_lun, i, j, k, nstlim = 1, NAT_for_mask1
      double precision :: resolution, fabs_solvent, phi_solvent
      complex(real_kind), allocatable, dimension(:) :: Fobs
      real(real_kind) :: phi
      integer :: has_f_solvent ! Keep it for compatibility with legacy input file format

      logical :: master=.true.
      ! following is local: copied into f_mask in this routine, after
      !     f_mask itself is allocated.  (could be simplified)
      if (pdb_infile /= '') call xray_read_pdb(trim(pdb_infile))

      if (reflection_infile == '') xray_active = .false.

      call unit_cell%init(pbc_box(1), pbc_box(2), pbc_box(3),  pbc_alpha, pbc_beta, pbc_gamma)

      if (.not.xray_active) then
         spacegroup_name = 'P 1'
         return
      end if

      write(stdout,'(A,3F9.3,3F7.2)') &
            'XRAY: UNIT CELL= ',pbc_box(1), pbc_box(2), pbc_box(3),  &
                                pbc_alpha, pbc_beta, pbc_gamma
      !--------------------------------------------------------------
      ! Read reflection data
      call amopen(allocate_lun(hkl_lun),reflection_infile,'O','F','R')
      read(hkl_lun,*,end=1,err=1) num_hkl, has_f_solvent

      allocate(hkl_index(3,num_hkl),abs_Fobs(num_hkl),sigFobs(num_hkl), &
            & test_flag(num_hkl))

      if (has_f_solvent > 0 ) then
         write(stdout,'(A)') 'f_solvent support is discontinued'
         call mexit(stdout, 1)
      endif

      if (fave_outfile /= '') then
         write(stdout,'(A)') 'fave_outfile is not yet implemented'
         call mexit(stdout, 1)
      endif

      !  each line contains h,k,l and two reals
      !  if target == "ls" or "ml", these are Fobs, sigFobs (for diffraction)
      !  if target == "vls",  these are Fobs, phiFobs (for cryoEM)

      do i = 1,num_hkl
         read(hkl_lun,*,end=1,err=1) &
            hkl_index(1:3,i),abs_Fobs(i),sigFobs(i),test_flag(i)
      end do

      if (atom_selection_mask/='') then
         call atommask(natom=natom,nres=nres,prnlev=0, &
               igraph=atm_igraph,isymbl=atm_isymbl,ipres=gbl_res_atms, &
               lbres=gbl_labres,crd=atm_crd, &
               maskstr=atom_selection_mask,mask=atom_selection)
         NAT_for_mask1 = sum(atom_selection)
         if( master ) write(6,'(a,i6,a,a)') 'Found ',NAT_for_mask1, &
              ' atoms in ', atom_selection_mask
         !  also ignore any atoms with zero occupancy:
         do i=1,natom
            if( atom_occupancy(i) == 0._rk_) atom_selection(i) = 0
         end do
         NAT_for_mask = sum(atom_selection)
         if( master .and. NAT_for_mask1 /= NAT_for_mask ) &
            write(6,'(a,i4,a)') 'Removing ',NAT_for_mask1 - NAT_for_mask, &
              ' additional atoms with zero occupancy'
      end if

      ! FIXME: account for phase in Fobs
      Fobs = abs_Fobs
      
      call init_interface2( &
         & target, bulk_solvent_model, &
         & hkl_index, Fobs, sigFobs, test_flag==1, &
         & unit_cell, scatter_coefficients, &
         & atom_bfactor, atom_occupancy, atom_scatter_type, &
         & atom_selection==1, atm_atomicnumber, &
         & mask_update_period, scale_update_period, &
         & ml_update_period, k_sol, b_sol, &
         & solvent_mask_adjustment, solvent_mask_probe_radius &
      & )
      
!      ! Dump xray init parameters to test/debug ! FIXME: remove it in release
!      call xray_dump( &
!          & "dump.txt", &
!          & hkl_index, Fobs, sigFobs, test_flag==1, &
!          & unit_cell_, scatter_coefficients, &
!          & atom_bfactor, atom_occupancy, atom_scatter_type, &
!          & atom_selection==1, atm_atomicnumber &
!      )

      return
      1 continue
      write(stdout,'(A)') 'Error reading HKL file.'
      call mexit(stdout,1)
   end subroutine init

   subroutine xray_init_globals()
      pdb_infile = ''
      pdb_outfile = ''
      fave_outfile = ''
      fmtz_outfile = ''
      sf_outfile = ''
      pdb_read_coordinates = .false.
      pdb_use_segid = .false.
      pdb_wrap_names = .false.
      target = 'ml'
      spacegroup_name = 'P 1'
      reflection_infile = ''
      xray_weight_initial = sentinel_xray_weight
      xray_weight_final = sentinel_xray_weight
      solvent_mask_adjustment = 1.11
      solvent_mask_probe_radius = 0.9
      solvent_mask_reflection_outfile = ''
      solvent_mask_outfile = ''
      atom_selection_mask = '!@H='
      mask_update_period = 100
      scale_update_period = 100
      ml_update_period = 100
      bulk_solvent_model = 'afonine-2013'
      return
   end subroutine xray_init_globals

   ! Write X-ray output files and deallocate. Bond info is included
   ! here only to check for places to insert TER in PDB output files.
   subroutine finalize()
      use xray_interface2_module, only: finalize2 => finalize
      implicit none
      ! local
      integer :: i
      real(real_kind) :: phi
      logical :: master = .true.

      if (.not.xray_active) return
      

      if (master .and. pdb_outfile /= '') then
         call xray_write_pdb(trim(pdb_outfile))
      end if
      if (master .and. fave_outfile /= '') then
         ! TODO: call xray_interface2::write_fave()
      endif

      if (master .and. fmtz_outfile /= '') then
         ! TODO: call xray_interface2::write_mtz_file()
      endif

      call finalize2()
      
      deallocate(atom_bfactor,atom_occupancy,atom_scatter_type, &
            atom_selection,residue_chainid,residue_icode, &
            atom_element,atom_altloc,residue_number, &
            hkl_index,abs_Fobs,sigFobs,test_flag)

   end subroutine finalize

   ! gets xray_energy and derivatives
   subroutine xray_get_derivative(xyz, force, current_step, xray_e)
      use mdin_ctrl_dat_mod, only: total_steps=> nstlim
      use xray_interface2_module, only: calc_force2 => calc_force, get_r_factors
      implicit none
      real(real_kind), intent(in) :: xyz(:, :)
      real(real_kind), intent(out) :: force(:, :)
      integer, intent(in) :: current_step
      real(real_kind), intent(out) :: xray_e
      ! local
      real(real_kind) :: xray_weight

      call check_precondition(size(xyz, 1) == 3)
      call check_precondition(size(xyz, 2) == size(force, 2))
      call check_precondition(size(force, 1) == 3)
!      call check_precondition(size(force, 2) == n_atom)

      if (.not. xray_active) then
         xray_e = 0
         xray_energy = 0
         return
      end if

      xray_weight = get_xray_weight(current_step, total_steps)

      call calc_force2(xyz, current_step, xray_weight, force, xray_e)
      xray_energy = xray_e
      call get_r_factors(r_work, r_free)

   end subroutine xray_get_derivative

   subroutine xray_write_md_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f14.4))') &
        'Exray  = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_md_state

   subroutine xray_write_min_state(unit)
      integer, intent(in) :: unit
      write(unit,'(3(1x,A,f13.4))') &
        'Exray   = ', xray_energy, ' Rwork   = ', r_work, ' Rfree      = ', r_free
   end subroutine xray_write_min_state


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Allocate the next unused logical UNIT number from a predefined range.
   ! Intended to be used with amopen, which aborts on error, because there
   ! is no corresponding deallocate function.
   function allocate_lun(lun_return) result(lun)
      integer :: lun
      integer, intent(out), optional :: lun_return
      integer, parameter :: FILE_UNIT_FIRST=201, FILE_UNIT_LAST=250
      integer, save :: unit = FILE_UNIT_FIRST-1
      ! locals
      integer :: i
      logical :: opened
      ! begin
      lun = FILE_UNIT_FIRST
      do i = FILE_UNIT_FIRST, FILE_UNIT_LAST
         unit = unit + 1
         if (unit > FILE_UNIT_LAST) unit = FILE_UNIT_FIRST
         ! This assumes that an allocated unit always results in an opened state.
         inquire(unit=unit,opened=opened)
         if (.not. opened) then
            lun = unit
            if (present(lun_return)) lun_return = unit
            return
         end if
      end do
      write(stdout,'(A)') &
              'ERROR in ALLOCATE_LUN(): ran out of available LUNs!'
      call mexit(stdout,2)
   end function allocate_lun
   
   function get_xray_weight(current_step, total_steps) result(result)
      integer, intent(in) :: current_step
      integer, intent(in) :: total_steps
      real(real_kind) :: result
      
      real(real_kind) :: weight_increment
      
      call check_precondition(current_step <= total_steps)
      
      if (total_steps > 1) then
         weight_increment = (xray_weight_final - xray_weight_initial) / (total_steps - 1)
      else
         weight_increment = 0
      end if
      
      result = xray_weight_initial + weight_increment * current_step
   end function get_xray_weight

end module xray_interface_impl_cpu_module
