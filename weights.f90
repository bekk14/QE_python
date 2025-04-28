module orbital_utils
implicit none
private

public :: orbitalsName, write_atom_weights,write_element_totals,get_unique_elements

logical, public :: verbosity = .false.

type ,public :: state_info
    integer :: state_index, atom_index, n, l, m, shell_index
    character(len=2) :: atom_symbol
    character(len=10) :: orbital
end type state_info



!public ::  write_global_totals
contains
subroutine orbitalsName(orbital_type, shell_index, orbital_name)
    ! Purpose: Returns the name of an atomic orbital based on its type and index.
    character(len=*), intent(in)  :: orbital_type  ! Input orbital type ('2S'/'3S'/'2P'/...)
    integer, intent(in)           :: shell_index   ! Index within orbital type (m-component index)
    character(len=32), intent(out) :: orbital_name ! Output orbital name

    character(len=20), parameter :: ORBITAL_NAMES(16) = [ &
        character(len=20) :: &
        's                  ', &  ! 1
        'p_z               ', &   ! 2
        'p_x               ', &   ! 3
        'p_y               ', &   ! 4
        'd_z2              ', &   ! 5
        'd_zx              ', &   ! 6
        'd_zy              ', &   ! 7
        'd_x2-y2           ', &   ! 8
        'd_xy              ', &   ! 9
        'f_z3              ', &   ! 10
        'f_zx2             ', &   ! 11
        'f_yz2             ', &   ! 12
        'f_zx2-y2          ', &   ! 13
        'f_xyz             ', &   ! 14
        'f_xx2-3y2         ', &   ! 15
        'f_y3x2-y2         ']     ! 16

    integer :: n_orbital
    character(len=1) :: suborbital

    orbital_name = 'Unknown orbital'
    read(orbital_type(1:1), '(I1)') n_orbital  ! Extract n (if needed)
    suborbital = orbital_type(2:2)             ! Extract orbital letter (S/P/D/F)

    select case (adjustl(trim(suborbital)))
        case ('S', 's')
            if (shell_index == 1) then
                orbital_name = trim(ORBITAL_NAMES(1))  ! s (only 1 state)
            else
                orbital_name = 'Invalid S index'
            end if

        case ('P', 'p')
            if (shell_index >= 1 .and. shell_index <= 3) then
                orbital_name = trim(ORBITAL_NAMES(1 + shell_index))  ! 2-4: p_z, p_x, p_y
            else
                orbital_name = 'Invalid P index'
            end if

        case ('D', 'd')
            if (shell_index >= 1 .and. shell_index <= 5) then
                orbital_name = trim(ORBITAL_NAMES(4 + shell_index))  ! 5-9: d_orbitals
            else
                orbital_name = 'Invalid D index'
            end if

        case ('F', 'f')
            if (shell_index >= 1 .and. shell_index <= 7) then
                orbital_name = trim(ORBITAL_NAMES(9 + shell_index))  ! 10-16: f_orbitals
            else
                orbital_name = 'Invalid F index'
            end if

        case default
            orbital_name = 'Unknown orbital'
    end select
end subroutine orbitalsName

subroutine write_atom_weights(dirname, atom_symbol, atom_id, orbital_names, weights, nkpt, nband, norb)
    character(len=*), intent(in) :: dirname, atom_symbol
    integer, intent(in)          :: atom_id, nkpt, nband, norb
    character(len=32), intent(in) :: orbital_names(norb)
    real*8, intent(in)             :: weights(norb, nkpt, nband)  ! Shape: (norb, nkpt, nband)
    
    character(len=256) :: filename
    integer :: u, ikpt, ibnd, iorb
    real*8 :: total, p_sum, d_sum, f_sum, s_val
    integer :: p_start, p_end, d_start, d_end, f_start, f_end
    logical :: s_exists, p_exists, d_exists, f_exists

    ! Create filename: e.g., 'out_pBands/atom_C_1_weights.dat'
    write(filename, '(a,"/atom_",a,"_",i0,"_wi.dat")') trim(dirname), trim(atom_symbol), atom_id

    open(newunit=u, file=trim(filename), status='replace')
    
    ! Detect orbital ranges
    s_exists = .false.
    p_start = 0; p_end = 0
    d_start = 0; d_end = 0
    f_start = 0; f_end = 0
    !
    do iorb = 1, norb
        if (trim(orbital_names(iorb)) == 's') then
            s_exists = .true.
        elseif (index(orbital_names(iorb), 'p_') == 1) then
            if(p_start == 0) p_start = iorb
            p_end = iorb
        elseif(index(orbital_names(iorb), 'd_') == 1) then
            if(d_start == 0) d_start = iorb
            d_end = iorb
        elseif(index(orbital_names(iorb), 'f_') == 1) then
            if(f_start == 0) f_start = iorb
            f_end = iorb
        endif
    enddo  

    p_exists = p_start > 0
    d_exists = d_start > 0
    f_exists = f_start > 0

    ! Write header
    if (verbose) then 
          write(u, '(a)', advance='no') '#ikpt  iband'
    else 
          write(u, '(a)', advance='no') '#ikpt  '
    endif 

    do iorb = 1, norb
      write(u, '(a16)', advance='no') trim(orbital_names(iorb))
    enddo

    if (s_exists) write(u, '(a16)', advance='no') 's_total'
    if (p_exists) write(u, '(a16)', advance='no') 'p_total'
    if (d_exists) write(u, '(a16)', advance='no') 'd_total'
    if (f_exists) write(u, '(a16)', advance='no') 'f_total'
    write(u, '(a16)') 'Total'

    !if 
    !write(u, '(a16,a16,a16,a16)') 's_total', 'p_total', 'd_total', 'f_total', 'Total'
    !write(u, *)  ! Newline after header
    
    ! Write data ! with sums 
    !do ikpt = 1, nkpt
    !  do ibnd = 1, nband
    do ibnd = 1 , nband 
        do ikpt =1 , nkpt 
        ! claculate sums 
        s_val = 0.0d0; p_sum = 0.0d0; d_sum = 0.0d0; f_sum = 0.0d0
                    ! Get s value
        if (s_exists) then
                do iorb = 1, norb
                    if (trim(orbital_names(iorb)) == 's') then
                        s_val = weights(iorb, ikpt, ibnd)
                        exit
                    endif
                enddo
        endif
        ! Calculate p/d/f sums
        if (p_exists) then
                do iorb = p_start, p_end
                    p_sum = p_sum + weights(iorb, ikpt, ibnd)
                enddo
        endif
        if (d_exists) then
                do iorb = d_start, d_end
                    d_sum = d_sum + weights(iorb, ikpt, ibnd)
                enddo
        endif
        if (f_exists) then
                do iorb = f_start, f_end
                    f_sum = f_sum + weights(iorb, ikpt, ibnd)
                enddo
        endif
        ! old method
        !do iorb = 1, norb
        !        if(iorb == 1) s_val = weights(iorb, ikpt, ibnd)
        !        if(iorb >= p_start .and. iorb <= p_end) p_sum = p_sum + weights(iorb, ikpt, ibnd)
        !        if(iorb >= d_start .and. iorb <= d_end) d_sum = d_sum + weights(iorb, ikpt, ibnd)
        !        if(iorb >= f_start .and. iorb <= f_end) f_sum = f_sum + weights(iorb, ikpt, ibnd)
        !enddo
        total = 0.0d0
        if (s_exists) total = total + s_val
        if (p_exists) total = total + p_sum
        if (d_exists) total = total + d_sum
        if (f_exists) total = total + f_sum
        !total = s_val + p_sum + d_sum + f_sum

        ! write data as old one 
        if (verbosity) then 
          write(u, '(i5, i5)', advance='no') ikpt, ibnd
        else 
         write(u, '(i5)', advance='no') ikpt
        endif 
        do iorb = 1, norb
          write(u, '(F16.10)', advance='no') weights(iorb, ikpt, ibnd)
        enddo
        !if 
        !write(u, '(4F16.10,F16.10)') s_val, p_sum, d_sum, f_sum, total
        !write(u, *)  ! Newline after each band
        if (s_exists) write(u, '(F16.10)', advance='no') s_val
        if (p_exists) write(u, '(F16.10)', advance='no') p_sum
        if (d_exists) write(u, '(F16.10)', advance='no') d_sum
        if (f_exists) write(u, '(F16.10)', advance='no') f_sum
        write(u, '(F16.10)') total
      enddo
      write(u, *)  ! Newline after each band      
    enddo
    close(u)
    print *, "Created file: ", trim(filename)

end subroutine write_atom_weights

subroutine write_element_totals(dirname, states, natoms, nshall, dataofwieghts, nkpt, nband)
    character(len=*), intent(in) :: dirname
    type(state_info), intent(in) :: states(*)
    integer, intent(in) :: natoms, nshall, nkpt, nband
    real*8, intent(in) :: dataofwieghts(natoms, nshall, nkpt, nband)
    
    character(len=32) :: name
    character(len=256) :: filename
    integer :: u, ikpt, ibnd, i, j, e, orb
    character(len=2), allocatable :: elements(:)
    integer :: nelements !, n_orb
    real*8, allocatable :: elem_totals(:,:,:) ! 3D orbitlas, kpts , nband 
    logical :: found
    
    ! Find unique elements
    allocate(elements(natoms))
    nelements = 0
    do i = 1, natoms
        found = .false.
        do e = 1, nelements
            if(elements(e) == states((i-1)*nshall+1)%atom_symbol) then
                found = .true.
                exit
            endif
        enddo
        if(.not. found) then
            nelements = nelements + 1
            elements(nelements) = states((i-1)*nshall+1)%atom_symbol
        endif
    enddo

        ! Process each element type
    do e = 1, nelements
        allocate(elem_totals(nshall, nkpt, nband))
        elem_totals = 0.0d0
            ! Sum weights for all atoms of this element
        do i = 1, natoms
            if(states((i-1)*nshall+1)%atom_symbol == elements(e)) then
                elem_totals = elem_totals + dataofwieghts(i,:,:,:)
            endif
        enddo

        ! Write element total file
        write(filename, '(a,"/atom_",a,"_tot.dat")') trim(dirname), trim(elements(e))
        open(newunit=u, file=trim(filename), status='replace')
        
        ! Write header
        if (verbosity)  then
            write(u, '(a)', advance='no') '# ikpt ibnd'
        else 
              write(u, '(a)', advance='no') '# ikpt '
        endif 
        write(u, '(a)', advance='no') '# ikpt ibnd'
        do orb = 1, nshall
            call orbitalsName(states(orb)%orbital, states(orb)%m, name)
            write(u, '(a16)', advance='no') trim(name)
        enddo
        write(u, '(a16)') 'Total'

        ! Write data
        !do ikpt = 1, nkpt
        !    do ibnd = 1, nband
        do ibnd = 1 , nband 
           do ikpt = 1, nkpt 
                if (verbosity)  then 
                write(u, '(i5,i5)', advance='no') ikpt, ibnd
                else 
                write(u, '(i5)', advance='no') ikpt
                endif 
                do orb = 1, nshall
                    write(u, '(F16.10)', advance='no') elem_totals(orb, ikpt, ibnd)
                enddo
                write(u, '(F16.10)') sum(elem_totals(:, ikpt, ibnd))
            enddo
            write(u,*)
        enddo
        
        close(u)
        deallocate(elem_totals)
    enddo
    
    deallocate(elements)
end subroutine write_element_totals

subroutine get_unique_elements(states, natoms, nshall, elements, nelements)
    type(state_info), intent(in) :: states(natoms*nshall)
    integer, intent(in) :: natoms, nshall
    character(len=2), allocatable, intent(out) :: elements(:)
    integer, intent(out) :: nelements
    integer :: i,j
    logical :: found
    
    allocate(elements(natoms))
    nelements = 0
    do i = 1, natoms
        found = .false.
        do j = 1, nelements
            if(elements(j) == states((i-1)*nshall+1)%atom_symbol) then
                found = .true.
                exit
            endif
        enddo
        if(.not. found) then
            nelements = nelements + 1
            elements(nelements) = states((i-1)*nshall+1)%atom_symbol
        endif
    enddo
end subroutine get_unique_elements

end module orbital_utils

program process_weights
    use orbital_utils
    !use iso_c_binding, only: c_sleep => sleep
    implicit none
    character(len=256) :: filename, dirname, arg, fmt_str
    character(len=1000) :: buffer, prev_line
    integer :: ios, natoms, ntype, nstates, nkpt, nband, nshall
    integer :: line_counter, current_state, stat, pos, max_l, dummy_int
    integer :: i, j, k, l, atom_idx, shall_idx, l_type, current_kpt, current_band
    integer ::  num_args 
    integer :: actual_kpt, actual_band, col, start_col, end_col, ikpt, ibnd, ierr, unit
    real*8 :: weight_value, weight
    logical :: ff_found, file_exists
    
    character(len=32) :: name
    integer :: prev_atom_index
    character(len=2) :: prev_atom_symbol
    integer :: current_shell
    type(state_info), allocatable :: states(:)
    ! 5D array: (atom, shall, orbital_type, kpt, band)
    !real, allocatable :: weights(:,:,:,:,:)
    ! 4D array: (atom, shall, orbital_type, kpt, band)
    real*8, allocatable :: dataofwieghts(:,:,:,:)
    integer :: ncols
    real :: dummy_reals(4)
    character(len=32), allocatable :: orbital_names(:)
    real*8, allocatable :: atom_weights(:,:,:) 

    print *, "                                                            ,,,  "
    print *, "                                                          /'^'\ "
    print *, "                                                         ( o o )"
    print *, "-------------------------------------------------------oOOO--(_)--OOOo------"
    print *, ""
    print *, "        ProjectionBands 2.01"
    print *, "          Fortran version  "
    print *, "        ---------------------"
    print *, "        (._.) Description:"
    print *, "              This script helps to analyze the data produced by code projwfc.x"
    print *, "              for band project, then extract the contribution of each state in "
    print *, "              each band using the weight of states."
    print *, "        (._.) Created by Bekk@ Hamza Lab. LaMCScI FSR/Morocco."
    print *, "        (._.) Coordinate: "
    print *, "                    .:. Dev: Bekkali Hamza"
    print *, "                    .:. Email: hamza_bekkali@um5.ac.ma"
    print *, "                        (C) 2022/2023"
    print *, "                                                        .oooO"
    print *, "                                                        (   )   Oooo. Bekk@ Hamza"
    print *, "---------------------------------------------------------\\ (----(   )-------"
    print *, "                                                          \\_)    ) /"
    print *, "                                                                (_/"
    print *, ""
    num_args = command_argument_count()
    ! Get input filename from command line
    !call get_command_argument(1, filename)
    !if (filename == '') then
    !    print *, "Error: No input file provided."
    !    stop
    !endif
    filename = ''
    !print *, num_args
    do i = 1, num_args
        call get_command_argument(i, arg)
        select case(trim(arg))
        case ('-v', '-V', '--verbose')
            verbosity = .true.
            if (verbosity) print *, "Verbose mode enabled"
        case default
            !print *,arg
            if (filename == '') then
                filename = trim(arg)
            else
                print *, "Error: Unknown option or extra argument: ", trim(arg)
                stop
            endif
        end select
    enddo
    call system('sleep 1')

    ! Open input file
    if (filename == '') then
        print *, "Error: No input file provided."
        print *, "Usage: ./program [-v] filename"
        stop
    endif
    
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        print *, "Error: File '", trim(filename), "' not found!"
        stop
    endif

    open(unit=10, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "Error opening file:", trim(filename)
        stop
    endif

    ! Read initial parameters --------------------------------------------------
    line_counter = 0
    read(10, '(a)', iostat=ios) buffer  ! skip first (possibly empty) line
    if (ios /= 0) then
        print *, "Error reading first line of file"
        stop
    endif
    line_counter = 1

    read(10, '(a)', iostat=ios) buffer  ! read second line
    if (ios /= 0) then
        print *, "Error reading second line of file"
        stop
    endif
    line_counter = 2
    ! Parse the second line which contains 8 integers
    read(buffer, *, iostat=ios) (dummy_int, i =1,6), natoms, ntype
    if (ios /= 0) then
        print *, "Error parsing parameters from line", line_counter
        print *, "Line content: ", trim(buffer)
        stop
    endif
    !print "(a,i0,a)", "Read from line ", line_counter, ":"
    ! Find F F marker ---------------------------------------------------------
    ff_found = .false.
    line_counter = 2
    ! Already read first two lines
    do while (.not. ff_found)
        prev_line = buffer
        read(10, '(a)', iostat=ios) buffer
        if (ios /= 0) exit
        line_counter = line_counter + 1
        ! Check for "F F" pattern with any spacing
        if (index(buffer, 'F') > 0) then
            if (scan(buffer, 'F', back=.true.) > index(buffer, 'F')) then
                ff_found = .true.
                ! Process previous line (line 35 in grep: 112 299 74)
                read(prev_line, *) nstates, nkpt, nband
                exit
            endif
        endif
    enddo
    if (.not. ff_found) then
       ! print *, "FATAL ERROR: marker not found!"
        print *, "Stopped searching at line", line_counter
        print *, "Last line read: ", trim(buffer)
        stop
    endif
    ! Validate and print critical parameters -----------------------------------
    nshall = nstates / natoms
    if (nshall * natoms /= nstates) then
        print *, "Error: nstates(", nstates, ") not divisible by natoms(", natoms, ")"
        stop
    endif
    print *, "Critical parameters:"
    print *, "--------------------------------------------------"
    print *, "  natoms  =", natoms 
    print *, "  ntype   =", ntype
    print *, "  nstates =", nstates
    print *, "  nkpt    =", nkpt
    print *, "  nband   =", nband
    print *, "  nshall  =", nshall
    print *, "--------------------------------------------------"
    print *, "Data to Read : ",nband*nkpt*nstates
    ! Read state headers and weight data -----------------------------------------
    allocate(states(nstates))
    max_l = 3  ! Maximum orbital type (s=0, p=1, d=2, f=3)
    ! Initialize 5D weights array: (atom, shall, orbital_type, kpt, band)
    !allocate(weights(natoms, nshall, max_l+1, nkpt, nband))
    allocate(dataofwieghts(natoms, nshall, nkpt, nband))
    !weights = 0.0
    dataofwieghts=0.0000000000
    ! Initialize previous values for shell index tracking
    prev_atom_index = -1
    prev_atom_symbol = '  '
    current_shell = 0
    rewind(10)
    do i = 1, line_counter
        read(10, '(a)') buffer
    end do
    ! Process each state
    do current_state = 1, nstates
        read(10, '(a)', iostat=ios) buffer
        if (ios /= 0) then
            print *, "Error reading line for state", current_state
            print *, trim(buffer)
            stop
        endif
        
        print *,"WFC(",current_state,")"
        read(buffer, *, iostat=ios) &
            dummy_int, &
            states(current_state)%atom_index, &
            states(current_state)%atom_symbol, &
            states(current_state)%orbital, &
            states(current_state)%n, &
            states(current_state)%l, &
            states(current_state)%m
        ! Determine shell index
        if (states(current_state)%atom_index == prev_atom_index .and. &
            states(current_state)%atom_symbol == prev_atom_symbol) then
            current_shell = current_shell + 1
        else
            current_shell = 1
        end if
        states(current_state)%shell_index = current_shell

        if (ios /= 0) then
            print *, "ERROR: Failed to parse state header:"
            print *, "Line content: [", trim(buffer), "]"
            stop
        endif
    
    
       if (verbosity) then 
        print *, "[:....................................."
        print "(a,i3,a,1x,a2,a,1x,a10,a,1x,a,a10,a,i2)", &
            " ", current_state, &
            "| ", trim(states(current_state)%atom_symbol), &
            states(current_state)%atom_index, &
            "|  ", trim(states(current_state)%orbital), &
            "| ", states(current_state)%n, states(current_state)%l, states(current_state)%m, &
            "| ", states(current_state)%shell_index
        print *, "   ......................................:]"
        endif 
        ! Read weight data for this state (nkpt × nband entries)
        ! Read data line by line
        do i = 1, nband * nkpt
            read(10, '(a)', iostat=ios) buffer
            if (ios /= 0) exit
            ! Parse kpt, band, weight
            read(buffer, *, iostat=ios) ikpt, ibnd, weight
            if (ikpt < 1 .or. ikpt > nkpt) then
                print *, "Warning: k-point index out of range:", ikpt
                stop
            endif
            
            if (ibnd < 1 .or. ibnd > nband) then
                print *, "Warning: band index out of range:", ibnd
                stop
            endif
            ! Store the weight in the 5D array
            atom_idx = states(current_state)%atom_index
            shall_idx = states(current_state)%shell_index
            l_type = states(current_state)%l + 1  ! Convert l to array index (s=1, p=2, etc.)
            dataofwieghts(atom_idx,shall_idx,ikpt,ibnd) = weight
        enddo
        ! Update previous values for shell index tracking
        prev_atom_index = states(current_state)%atom_index
        prev_atom_symbol = states(current_state)%atom_symbol
    enddo 



    ! Create output directory
    dirname = 'out_pBands/' // trim(filename) // '/'
    call execute_command_line('mkdir -p ' // trim(dirname), exitstat=stat)
    if (stat /= 0) then
        print *, "Error creating directory:", trim(dirname)
        stop
    endif


    ! Process each atom
    do i = 1, natoms
        ! Extract orbital names for this atom
        allocate(orbital_names(nshall))
        do j = 1, nshall
            current_state = (i-1)*nshall + j
            if (current_state > nstates) exit
            call orbitalsName(states(current_state)%orbital, states(current_state)%m, orbital_names(j))
        enddo
        ! Extract weights for this atom (shape: norb x nkpt x nband)
        allocate(atom_weights(nshall, nkpt, nband))
        atom_weights = dataofwieghts(i, 1:nshall, 1:nkpt, 1:nband)
        ! Write to file
        call write_atom_weights(trim(dirname), &
                                trim(states((i-1)*nshall + 1)%atom_symbol), &
                                i, orbital_names, atom_weights, nkpt, nband, nshall)

        deallocate(orbital_names, atom_weights)
    enddo
    ! In the main program, after processing all atoms:
    !call write_global_totals(trim(dirname), states, natoms, nshall, dataofwieghts, nkpt, nband)
    
    !call write_element_totals(trim(dirname), states, natoms, nshall, dataofwieghts, nkpt, nband)
    call write_element_totals(trim(dirname), states, natoms, nshall, dataofwieghts, nkpt, nband)
    close(10)
    print *, "Processing completed successfully. Output in: ", trim(dirname)

contains
   ! Add this helper function in the main program
    character(len=20) function int2str(i)
        integer, intent(in) :: i
        write(int2str, *) i
        int2str = adjustl(int2str)
    end function
end program process_weights
