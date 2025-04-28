!@ update by Bekk@14 hamza
! hamza_bekkali@um5.ac.ma
! Atomic position conversion utilities from input to internal to output format
!
!-----------------------------------------------------------------------
SUBROUTINE convert_tau( tau_format, nat, tau )
  !................................................... don't change the original code 
END SUBROUTINE convert_tau
!
!----------------------------------------------------------------------------
SUBROUTINE output_tau( print_lattice, print_final )
  ! ................................keep the original code
  
  !
  IF ( print_final  ) WRITE( stdout, '("End final coordinates")') 
  !--------------------------------------------------------------- add code 
  ! After writing "End final coordinates", CONVERT to POSCAR
  IF ( print_final ) CALL write_poscar(nat, atm, ityp, tau_out, cell_units, tau_format)
  !--------------------------------------------------------------- end add code 
  WRITE( stdout, '(/)' )
  !
  DEALLOCATE( tau_out )
  !
  RETURN
  !
END SUBROUTINE output_tau


SUBROUTINE output_tau_rescaled(rescale)
  
END SUBROUTINE

  !-----------------------------------------------------------------------------
  !   SAVE the final coordinates [cell and positions atomic]  as POSCAR.vasp to be easy to read by VESTA  
  !------------------------------------------------------------------------------
  !----------------------------------------------------------------------------
SUBROUTINE write_poscar(nat, atm, ityp, tau_out, cell_units, tau_format)
  !----------------------------------------------------------------------------
  !! Writes current structure to POSCAR_<step>.vasp in VASP format.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : ionode, stdout
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : at, alat
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nat
  CHARACTER(LEN=3), INTENT(IN) :: atm(*)
  INTEGER, INTENT(IN) :: ityp(nat)
  REAL(DP), INTENT(IN) :: tau_out(3,nat)
  CHARACTER(LEN=*), INTENT(IN) :: cell_units, tau_format
  !
  ! Local variables
  INTEGER, SAVE :: counter = 0
  INTEGER :: un, i, na, unique_count, current_type
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=3), ALLOCATABLE :: unique_atm(:)
  INTEGER, ALLOCATABLE :: atm_count(:)
  REAL(DP) :: cell(3,3)
  LOGICAL :: is_new

  IF (.NOT. ionode) RETURN  ! Only master process writes

  ! Increment counter and create filename
  counter = counter + 1
  WRITE(filename, '("POSCAR_",I0,".vasp")') counter

  ! Get cell parameters in Angstrom
  SELECT CASE(TRIM(cell_units))
  CASE('angstrom')
     cell = at * alat * bohr_radius_angs
  CASE('bohr')
     cell = at * alat * bohr_radius_angs
  CASE('alat')
     cell = at * alat * bohr_radius_angs
  END SELECT

  ! Identify unique atomic species and counts
  ALLOCATE(unique_atm(nat), atm_count(nat))
  unique_count = 0
  DO na = 1, nat
     is_new = .TRUE.
     DO i = 1, unique_count
        IF (TRIM(atm(ityp(na))) == TRIM(unique_atm(i))) THEN
           atm_count(i) = atm_count(i) + 1
           is_new = .FALSE.
           EXIT
        END IF
     END DO
     IF (is_new) THEN
        unique_count = unique_count + 1
        unique_atm(unique_count) = atm(ityp(na))
        atm_count(unique_count) = 1
     END IF
  END DO

  ! Write POSCAR file
  OPEN(NEWUNIT=un, FILE=filename, STATUS='unknown')
  WRITE(un, '("final coordinates from pwoscf Quantum ESPRESSO")')
  WRITE(un, '("1.0")')  ! Scale factor
  WRITE(un, '(3F20.10)') cell(1:3,1)  ! Lattice vectors
  WRITE(un, '(3F20.10)') cell(1:3,2)
  WRITE(un, '(3F20.10)') cell(1:3,3)
  WRITE(un, '(*(A3,1X))') (TRIM(unique_atm(i)), i=1,unique_count)  ! Species
  WRITE(un, '(*(I5,1X))') (atm_count(i), i=1,unique_count)  ! Counts

  IF (TRIM(tau_format) == 'crystal') THEN
     WRITE(un, '("Direct")')
  ELSE
     WRITE(un, '("Cartesian")')
  END IF

  WRITE(un, '(3F20.10)') (tau_out(1:3,na), na=1,nat)  ! Positions

  CLOSE(un)
  DEALLOCATE(unique_atm, atm_count)

END SUBROUTINE write_poscar
