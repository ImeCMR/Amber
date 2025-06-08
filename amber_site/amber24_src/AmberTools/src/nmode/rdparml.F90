!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program rdparml here]
subroutine rdparml(parm, amass)

   implicit double precision (a-h, o-z)
   double precision amass(*)
#include "lmsizes.h"
   character(len=80) fmt,parm
   character(len=80) fmtin,ifmt,afmt,rfmt,type
   character(len=4) dumc(ma)
   double precision dumf(ma)

   ifmt = '(12I6)'
   afmt = '(20A4)'
   rfmt = '(5E16.8)'
   
   !     ----- Read preliminary info from PARM file
   
   nf = 51
   call amopen(nf, parm, 'O', 'F', 'R')
   fmtin = afmt
   type = 'TITLE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt)  (dumc(i), i=1,20)
   
   fmtin = ifmt
   type = 'POINTERS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt)   natom, ntypes, nbonh,  mbona,  ntheth, mtheta,      &! 5
         nphih, mphia,  nhparm, nparm,  nnb,    nres,        &! 10
         nbona, ntheta, nphia,  numbnd, numang, nptra,       &! 15
         natyp, nphb, ifpert,idum,idum,idum,idum,idum, &
         idum,idum,idum,idum,idum
   nr3 = 3 * natom

   !     ----- Read rest of the parm file
   
   ntype = ntypes * ntypes
   nttyp = ntypes * (ntypes+1) / 2
   
   !     ----- read atomic and residue properties of the molecule(s)
   
   fmtin = afmt
   type = 'ATOM_NAME'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt) ( dumc(i), i = 1,natom )
   
   fmtin = rfmt
   type = 'CHARGE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (dumf(i), i=1,natom )
   
   fmtin = rfmt
   type = 'MASS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (amass(i),i=1,natom)
   
#if 0   /* we don't, for now, other sections are not needed  */

   fmtin = ifmt
   type = 'ATOM_TYPE_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=miac,miac+natom-1)
   
   fmtin = ifmt
   type = 'NUMBER_EXCLUDED_ATOMS'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=miblo,miblo+natom-1)
   
   fmtin = ifmt
   type = 'NONBONDED_PARM_INDEX'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=mico,mico+ntype-1)
   
   fmtin = afmt
   type = 'RESIDUE_LABEL'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt) ( ih(i), i = mlbres, mlbres + nres - 1 )
   
   fmtin = ifmt
   type = 'RESIDUE_POINTER'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=mipres,mipres+nres-1)
   
   !     ----- the last ipres array element should indicate the
   !     ----- beginning of a residue beyond the last one as if
   !     ----- it were present
   
   ix(mipres + nres) = natom + 1
   
   !     ----- read parameters such as bond lengths, force constants
   
   fmtin = rfmt
   type = 'BOND_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mrk,mrk+numbnd-1)
   
   fmtin = rfmt
   type = 'BOND_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mreq,mreq+numbnd-1)
   
   fmtin = rfmt
   type = 'ANGLE_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mtk,mtk+numang-1)
   
   fmtin = rfmt
   type = 'ANGLE_EQUIL_VALUE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mteq,mteq+numang-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_FORCE_CONSTANT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mpk,mpk+nptra-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PERIODICITY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mpn,mpn+nptra-1)
   
   fmtin = rfmt
   type = 'DIHEDRAL_PHASE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mphase,mphase+nptra-1)
   
   fmtin = rfmt
   type = 'SOLTY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=msolty,msolty+natyp-1)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mcn1,mcn1+nttyp-1)
   
   fmtin = rfmt
   type = 'LENNARD_JONES_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mcn2,mcn2+nttyp-1)
   
   !     ----- read bonding information
   
   fmtin = ifmt
   type = 'BONDS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ mibh-1), &
         ix(i+ mjbh-1), &
         ix(i+micbh-1), &
         i=1,nbonh)
   
   fmtin = ifmt
   type = 'BONDS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ miba-1), &
         ix(i+ mjba-1), &
         ix(i+micba-1), &
         i=1,nbona)
   
   fmtin = ifmt
   type = 'ANGLES_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ mith-1), &
         ix(i+ mjth-1), &
         ix(i+ mkth-1), &
         ix(i+micth-1), &
         i=1,ntheth)
   
   fmtin = ifmt
   type = 'ANGLES_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ mita-1), &
         ix(i+ mjta-1), &
         ix(i+ mkta-1), &
         ix(i+micta-1), &
         i=1,ntheta)
   
   fmtin = ifmt
   type = 'DIHEDRALS_INC_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ miph-1), &
         ix(i+ mjph-1), &
         ix(i+ mkph-1), &
         ix(i+ mlph-1), &
         ix(i+micph-1), &
         i=1,nphih)
   
   fmtin = ifmt
   type = 'DIHEDRALS_WITHOUT_HYDROGEN'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt)      (ix(i+ mipa-1), &
         ix(i+ mjpa-1), &
         ix(i+ mkpa-1), &
         ix(i+ mlpa-1), &
         ix(i+micpa-1), &
         i=1,nphia)
   
   fmtin = ifmt
   type = 'EXCLUDED_ATOMS_LIST'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i+minb-1),i=1,nnb)
   
   !     ----- read H-bond parameters
   
   fmtin = rfmt
   type = 'HBOND_ACOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=masol,masol+nphb-1)
   
   fmtin = rfmt
   type = 'HBOND_BCOEF'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mbsol,mbsol+nphb-1)
   
   fmtin = rfmt
   type = 'HBCUT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (x(i),i=mhbcut,mhbcut+nphb-1)
   
   !     ----- read the symbol, tree, join, irotat arrays:
   
   fmtin = afmt
   type = 'AMBER_ATOM_TYPE'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt) ( ih(i), i = msymbl, msymbl + natom - 1 )
   
   fmtin = afmt
   type = 'TREE_CHAIN_CLASSIFICATION'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read (nf,fmt) ( ih(i), i = mitree , mitree  + natom - 1 )
   
   fmtin = ifmt
   type = 'JOIN_ARRAY'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=mjoin,mjoin+natom-1)
   
   fmtin = ifmt
   type = 'IROTAT'
   call nxtsec(nf,  6,  0,fmtin,  type,  fmt,  iok)
   read(nf,fmt) (ix(i),i=mrotat,mrotat+natom-1)
   
   !   ---duplicate the multiple-term torsions:
   
   idum = nphih
   call dihdup(nphih,idum,ix(miph),ix(mjph),ix(mkph),ix(mlph), &
         ix(micph),x(mpn),maxdup)
   call dihdup(mphia,nphia,ix(mipa),ix(mjpa),ix(mkpa),ix(mlpa), &
         ix(micpa),x(mpn),maxdup)

   
   !     ----- Scale the charges by the dielectric constatnt
   
   if (dielc /= 1.0) then
      sqdiel  = sqrt (dielc)
      do 10 i = 1,natom
         cg(i) = cg(i) / sqdiel
      10 continue
   end if

#endif

   close (nf)
   return
   
end subroutine rdparml 

