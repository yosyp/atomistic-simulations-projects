!     Leonid Zhigilei, University of Virginia
!     This program implements the Ising model in 2D
!     MSE 4270/6270: Introduction to Atomistic Simulations
!     2D square lattice used in the model

! The program starts form the line 137 after "program MMC"
module MC_procedures 
    implicit none
    integer, parameter :: Nx = 100, Ny = 100        ! Numbers of the lattice sites in X and Y directions
    integer, parameter :: Nsteps = 10000000         ! Number of Monte Carlo iterations
    real, parameter :: T = 100.0                    ! Temperature [K]
    integer, parameter :: n_interactions = 3        ! Number of possible interactions: AA,BB, AB -> 3
                                                    ! If you change it, change also "get_interaction" function
    real, parameter :: E(n_interactions) = (/-0.05, -0.05, -0.25/) ! Energy of AA, BB, AB interactions [eV]
    integer, parameter :: Nwrite_s = 50000          ! Frequency of writing the output file to "data" directory
    integer, parameter :: Nwrite_E = 500            ! Frequency of writing to the output file EnNB.out
        
    integer, parameter :: n_neighbors = 4           ! Number of nearest neighbors
    real, parameter :: kT = T*8.617385E-05          ! kT [eV]
    integer, parameter :: Natoms = Nx*Ny            ! Total number of the lattice sites                                           
    integer, parameter :: dim = 2
                                                    ! Lattice displacement to nearest neighbors
    integer, parameter :: nn_vec(dim,n_neighbors) = reshape((/0, 1, 1, 0, -1, 0, 0, -1/),shape(nn_vec))
    
    integer :: types(Natoms)                        ! Types of atoms
    integer :: neighbors(n_neighbors, Natoms)       ! List of neighbors for all the atoms
    
    contains
    !get an energy of atom with id "atom_id"
    real function get_E(atom_id, N)
        implicit none
        integer, intent(in) :: atom_id
        integer, intent(out) :: N(n_interactions)   ! Number of interactions of a particular kind
        integer :: i, interaction_id

        get_E = 0.0
        N(:) = 0
        if((atom_id < 1).or.(atom_id > Natoms)) then
            write(*,*) "get_E: invalid input"
            return
        endif
        do i = 1, n_neighbors
            interaction_id = get_interaction(types(atom_id), types(neighbors(i,atom_id)))
            get_E = get_E + E(interaction_id)
            N(interaction_id) = N(interaction_id) + 1
        enddo

        return
    endfunction

    !get the total energy of the system
    real function get_Etot(N)
        implicit none
        integer, intent(out) :: N(n_interactions)       ! Number of interactions of a particular kind
        integer :: i, Ncur(n_interactions)

        N(:) = 0 
        get_Etot = 0.0
        do i = 1, Natoms
            get_Etot = get_Etot + get_E(i, Ncur)
            N(:) = N(:) + Ncur(:)
        enddo
        get_Etot = 0.5*get_Etot                         ! Double counting
        N(:) = N(:)/2

        return
    endfunction
    
    subroutine build_nl()
        implicit none
        integer :: pos(dim), pos_n(dim)                 ! position of an atom and its neighbor
        integer :: i, j
            
        do i = 1, Natoms
            pos(:) = get_pos(i)
            do j = 1, n_neighbors
                pos_n(:) = pos(:) + nn_vec(:,j)
                !apply periodic boundary condition, make sure that index is within [1,N_]
                !since Fortran "mod" of a negative number is negative, add N_
                pos_n(:) = mod(pos_n(:) + (/Nx, Ny/) - 1, (/Nx, Ny/)) + 1
                neighbors(j,i) = get_id(pos_n)
            enddo
        enddo
    
        return
    endsubroutine

    !get id of a particular interaction (1 - AA, 2 - BB, 3 - AB)
    integer pure function get_interaction(id1, id2)
        implicit none
        integer, intent(in) :: id1, id2

        if((id1 == 1).and.(id2 == 1)) then
            get_interaction = 1
        else if((id1 == 2).and.(id2 == 2)) then
            get_interaction = 2
        else
            get_interaction = 3
        endif

        return
    endfunction
    
    !return 2D index of an atom with id "atom_id"
    pure function get_pos(atom_id)
        integer, intent(in) :: atom_id
        integer :: get_pos(dim)

        get_pos(1) = mod(atom_id - 1, Nx) + 1
        get_pos(2) = int((atom_id - 1)/Nx) + 1

        return
    endfunction

    !return id of an atom based on the 2D index
    integer pure function get_id(pos)
        integer, intent(in) :: pos(dim)

        get_id = pos(1) + (pos(2) - 1)*Nx

        return
    endfunction
    
    subroutine swap(id1, id2)
        integer, intent(inout) :: id1, id2
        integer :: tmp
        
        tmp = types(id1)
        types(id1) = types(id2)
        types(id2) = tmp
        
        return
    endsubroutine
endmodule


program MMC
    use MC_procedures
    implicit none
    integer :: Nbonds(n_interactions)               ! Number of AA, BB, AB bonds
    integer :: NA, NB                               ! Number of A and B atoms
    integer :: Nbonds_i(n_interactions), Nbonds_f(n_interactions) ! Number of local bonds before and after the move  
    integer :: N1, N2                               ! Two atoms to chosen to swap
    real:: Etot                                     ! Total energy
    real:: E_i, E_f                                 ! Local energy before and after the move
    integer:: Nreject                               ! Number of rejected moves
    character*9 :: str                              ! String for making names for the output files
    integer :: file_E, file_s                       ! IDs of file for writing the output
    integer :: step                                 ! the number of a step
    real :: rnd, W                                  ! a random number and MC exponent
    logical :: accepted
    integer :: i, Nbonds_tmp(3)
    
!   Distributing initial atom types at random
    do i = 1, Natoms
        call random_number(rnd)
        if(rnd > 0.5) then
            types(i) = 1
        else 
            types(i) = 2
        endif
    enddo  
!   Calculating numbers of particles
    NA = 0; NB=0
    do i = 1, Natoms
        if(types(i) == 1) NA = NA + 1
        if(types(i) == 2) NB = NB + 1
    enddo
    write(*,*) "Number of atoms:", Natoms, " Type1:", NA, " Type2:", NB
!   Making list of neighbors for all atoms
    call build_nl()
!   Calculating the initial energy and number of bonds
    Etot = get_Etot(Nbonds)
    open(newunit = file_E, file = "EnNB.out")
    Nreject = 0
    
!   Metropolis method starts here *************************
    do step = 0, Nsteps
!       Choosing two atoms to swap (N1 and N2)
        do
            call random_number(rnd)
            N1 = int(rnd*Natoms) + 1
            call random_number(rnd)
            N2 = int(rnd*Natoms) + 1
            if(types(N1).ne.types(N2)) exit
        enddo
!       Calculating energy and number of affected bonds before the swap
        E_i = get_E(N1, Nbonds_tmp)
        Nbonds_i(:) = Nbonds_tmp(:)
        E_i = E_i + get_E(N2, Nbonds_tmp)
        Nbonds_i(:) = Nbonds_i(:) + Nbonds_tmp(:)
!       Swapping the atoms
        call swap(N1, N2)
!       Calculating energy and number of affected bonds after the swap
        E_f = get_E(N1, Nbonds_tmp)
        Nbonds_f(:) = Nbonds_tmp(:)
        E_f = E_f + get_E(N2, Nbonds_tmp)
        Nbonds_f(:) = Nbonds_f(:) + Nbonds_tmp(:)
!       Deciding if we accept the move/swap
        accepted = .false.
        if(E_f > E_i) then
            W = exp(-(E_f - E_i)/kT)
            call random_number(rnd)
            if(rnd < W) accepted = .true.
        else
            accepted = .true.
        endif
        if(accepted) then
            Etot = Etot + E_f - E_i
            Nbonds(:) = Nbonds(:) + Nbonds_f(:) - Nbonds_i(:)
        else
!           Swapping the atoms back to the original state
            call swap(N1, N2)
            Nreject = Nreject + 1
        endif
        if(mod(step, Nwrite_E) == 0) write(*,*) "Step=", step, " f_reject=", real(Nreject)/max(step,1)
!       Writing output        
        if(mod(step, Nwrite_E) == 0) then
            write(file_E, "(I10,1x,F12.2,3(1x,I8))") step, Etot, Nbonds(1:3)
            flush(file_E)
        endif
        if(mod(step, Nwrite_s) == 0) then
            write(str, FMT='(I9.9)') step
            open(newunit = file_s, file = './data/mc'//trim(str)//'.d')
            write(file_s, "(2(I5,1x),I1,1x,I7)") (get_pos(i), types(i), i, i = 1, Natoms)
            close(file_s)
        endif
    enddo
    close(file_E)

    write(*,*) T 

endprogram
    
           