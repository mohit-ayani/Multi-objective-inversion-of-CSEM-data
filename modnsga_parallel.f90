! This module contains subroutines to carry out the
! parallel NSGA II inversion of CSEM data
! Author: Mohit Ayani(mohitayani@gmail.com)
! Dated: 06/19/18
!***************************************

module nsga
	use dipole1d
	use runfile
	use ieee_arithmetic
	use mpi
	
	implicit none
	!-------------------------------
	! defining the data types
	!-------------------------------
	type individual
		real(8), dimension(:), allocatable :: model_parameter, obj
		integer :: rank, dominates_i_counter, i_dominates_counter
		integer, allocatable, dimension(:) :: i_dominates_individual
		real(8) :: crowd_dist
		integer :: index                    ! used in tourn_select routine to make sure
										    ! the population member is not repeated
	end type individual
	
	type front_t
		integer, dimension(:), allocatable :: member
	end type front_t
	
	!---------------------------------------
	!       define constants 
	!---------------------------------------
	real(8), parameter, public :: EPSILON = 1e-6
	real(8), parameter, public :: noise_floor = 1e-15
	integer, parameter, public :: master_sending_data = 1
	integer, parameter, public :: worker_sending_data = 2
	
	!---------------------------------------
	! define public variables, subroutines
	!---------------------------------------
	
	!
	! define GA variables
	!
	integer, public :: idum, popsize, igen, maxgen, mating_pool, nobj,nparam
	type(individual), dimension(:), allocatable, public :: parent_pop, mating_pop,&
	tourn_pop, child_pop, intermediate_pop, parent_pop_global
	type(individual), dimension(:),allocatable, public :: acceptable_models ! used to store models with misfit 1
	type(individual), public :: par1, par2, child1, child2
	type(front_t), dimension(:), allocatable, public :: front_parent, front_intermediate
	character(1), public :: c
	
	!
	! mpi variables
	!
	integer, public :: procnum, nproc, ierr, points_per_proc
	integer, public :: istart, iend
	integer, dimension(mpi_status_size) :: status
	real(8) :: start_time, end_time
	! variables for sending rank, crowd_dist to workers
	integer, public :: beg, fin, mycounts
	integer, public :: acceptable_member_count      ! used to store models with misfit 1
	integer, dimension(:), allocatable :: counts_masterproc, disp_masterproc
	real(8), dimension(:), allocatable :: tempcrowdist_global, tempcrowdist_local
	integer, dimension(:), allocatable :: temprank_global, temprank_local
	real(8), dimension(:,:), allocatable :: temparray_obj_parent, temparray_obj_child,&
	temparray_mp_parent, temparray_mp_child
	
	!
	! EM variables
	!
	complex(8), dimension(:,:), allocatable, public :: data_true,data_syn
	real(8), dimension(:,:),allocatable, public ::logamp_true,logamp_syn,&
	phase_true, phase_syn, stdamp, stdphase, stdreal, stdimag
	real(8), dimension(30) :: residual_real, residual_imag, residual_amp,residual_phase 
	
	!
	! public subroutines
	!
	public :: initialize_ga, init_defaults
	public :: norm_rand, initialize_em
	public :: calc_objective, send_data_to_master
	public :: recv_data_from_worker, non_dominated_sort
	public :: calc_crowding_dist, tourn_select
	public :: cross_over,mutation
	
	!---------------------------------------
	! defining private variables, subroutines
	!---------------------------------------
	!
	! private variables
	!
	real(8), public :: pcross, n_distribution_c, pmutate, n_distribution_m
	integer, private :: tournsize
	real(8), dimension(:), private, allocatable :: parmin, parmax
	logical, private :: RIGID
	
	!
	!private subroutines
	!
	private :: create_children,get_beta, get_delta
	
contains
    !=======================================
    !       subroutine init_defaults
    !To intialize dipole related parameters
    !=======================================
    subroutine init_defaults 

    use dipole1d     ! Model parameters are passed in here, fields passed out here
    !
    ! Specify some parameters required by Dipole1D:
    !
    HTmethod1D      = 'kk_ht_201' 
    outputdomain1D  = 'spatial'
    lbcomp          = .true.
    sdm1D           = 1.0   ! (Am), dipole moment 
    lUseSpline1D    = .false.
    linversion      = .false. ! Compute Derivatives with respect to sigma(layers)
    end subroutine init_defaults 
	!=======================================
	!   subroutine initialize_ga
	!=======================================
	subroutine initialize_ga
		implicit none
		!local variables
		integer :: i,j
		integer :: nseed
		integer, allocatable, dimension(:) :: seed
		real(8) :: rtemp
		! initialize the variables
		idum = 2891
		popsize = 1024
		maxgen = 2000
		nparam = 53
		nobj = 2 
		pcross = 0.9
		pmutate = 0.1
		n_distribution_c = 20
		n_distribution_m = 20
		acceptable_member_count = 0
		RIGID = .true.
		
		! define no. of points to be handled by each proc
		points_per_proc = (popsize + nproc -1 )/nproc
		istart = procnum*points_per_proc + 1
		iend = min((procnum + 1)*points_per_proc, popsize)
		print *,'istart = ',istart, 'iend = ',iend, 'procnum',procnum
		
		! allocate arrays for each proc. Array range for each
		! proc is 1 -->iend - istart + 1
		allocate(parent_pop(iend - istart + 1))
		allocate(child_pop(iend - istart + 1))
		allocate(parmin(nparam), parmax(nparam))
		do i = 1, iend - istart + 1
			allocate(parent_pop(i)%model_parameter(nparam))
			allocate(child_pop(i)%model_parameter(nparam))
			allocate(parent_pop(i)%obj(nobj))
			allocate(child_pop(i)%obj(nobj))
			parent_pop(i)%crowd_dist = 0.d0
			child_pop(i)%crowd_dist = 0.d0
		end do
		mating_pool = (iend - istart + 1)/2
		tournsize = 2
		allocate(tourn_pop(tournsize))
		allocate(mating_pop(mating_pool))
		allocate(par1%model_parameter(nparam),par2%model_parameter(nparam))
		allocate(child1%model_parameter(nparam),child2%model_parameter(nparam))
		allocate(par1%obj(nobj),par2%obj(nobj),child1%obj(nobj),child2%obj(nobj))
		par1%crowd_dist = 0.d0
		par2%crowd_dist = 0.d0
		child1%crowd_dist = 0.d0
		child2%crowd_dist = 0.d0
		parmin = -2.3010
		parmax = 1.d0
		
		!***this section is only for the master node ***
		if(procnum == 0)then
			allocate(parent_pop_global(popsize))
			do i = 1, popsize
				allocate(parent_pop_global(i)%obj(nobj))
				allocate(parent_pop_global(i)%model_parameter(nparam))
				parent_pop_global(i)%crowd_dist = 0.d0
			end do
			allocate(intermediate_pop(2*popsize))
			do i = 1, 2*popsize
				allocate(intermediate_pop(i)%model_parameter(nparam))
				allocate(intermediate_pop(i)%obj(nobj))
				intermediate_pop(i)%crowd_dist = 0.d0
			end do
		end if
		!*********** end section ****************
		
		!*** following arrays are used for scatter*
		!gather of rank,crowd_dist b/w slaves & master
		allocate(counts_masterproc(nproc),disp_masterproc(nproc))
		allocate(temprank_global(popsize),tempcrowdist_global(popsize))
		allocate(temprank_local(iend - istart + 1),tempcrowdist_local(iend -istart +1))
		if(.not. allocated(temparray_obj_parent)) allocate(temparray_obj_parent(points_per_proc,nobj))
		if(.not. allocated(temparray_obj_child)) allocate(temparray_obj_child(points_per_proc,nobj))
		if(.not. allocated(temparray_mp_parent)) allocate(temparray_mp_parent(points_per_proc,nparam))
		if(.not. allocated(temparray_mp_child)) allocate(temparray_mp_child(points_per_proc,nparam))
		!****** end section *****
		
		! using different seed for random no. generation
		call random_seed(size = nseed)
		allocate(seed(nseed))
		do i = 1, nseed
			seed(i) = idum + 37*(i-1) + 97*procnum
		end do
		call random_seed(put = seed)
		deallocate(seed)
		
		do i = 1, iend - istart + 1
			do j = 1, nparam
				call random_number(rtemp)
				parent_pop(i)%model_parameter(j) = parmin(j)*(1.d0 - rtemp) + parmax(j)*rtemp
				parent_pop(i)%index = istart + i - 1
			end do
		end do
	end subroutine initialize_ga
	!=======================================
	!       subroutine initialize_em
	!=======================================
	subroutine initialize_em
		implicit none
		!local variables
		integer :: i,j
		real(8) :: realpart, imagpart, noise1,noise2
		
		!allocate the arrays
		allocate(data_true(nTx,nFreq),data_syn(nTx, nFreq))
		allocate(stdreal(nTx,nFreq), stdimag(nTx, nFreq))
		allocate(stdamp(nTx, nFreq), stdphase(nTx, nFreq))
		allocate(logamp_true(nTx,nFreq), logamp_syn(nTx, nFreq))
		allocate(phase_true(nTx, nFreq), phase_syn(nTx,nFreq))
		data_true = modelresponse
	
		! adding normally distributed random noise to the data
		do i = 1, nTx
			do j = 1, nFreq
				realpart = real(data_true(i,j))
				imagpart = aimag(data_true(i,j))
				stdreal(i,j) = max(noise_floor/realpart ,0.05)*realpart
				stdimag(i,j) = max(noise_floor/imagpart, 0.05)*imagpart
				noise1 = stdreal(i,j)*norm_rand(0.0,1.0)
				noise2 = stdimag(i,j)*norm_rand(0.0,1.0)
				realpart = realpart + noise1
				imagpart = imagpart + noise2
				data_true(i,j) = cmplx(realpart, imagpart)
				logamp_true(i,j) = log10(abs(data_true(i,j)))
				phase_true(i,j) = DATAN2(aimag(data_true(i,j)),real(data_true(i,j)))*(180/3.1416)
				stdamp(i,j) = 0.05*(1/log(10.00))
				stdphase(i,j) = 0.05*(180/3.1416)
			end do
			
		end do
	end subroutine initialize_em
	!=======================================
	!       function norm_rand
	!=======================================
	function norm_rand(mean, std_dev)
	    real :: norm_rand
	    real, intent(in) :: mean, std_dev
	    real :: x, y, r
	    real, save :: spare
	    logical, save :: has_spare
	    ! use a spare saved from a previous run if one exists
	    if (has_spare) then
	        has_spare = .FALSE.
	        norm_rand = mean + (std_dev * spare)
	        return
	    else
	        r = 1.0
	        do while ( r >= 1.0 )
	            ! generate random number pair between 0 and 1
	            call random_number(x)
	            call random_number(y)
	            ! normalise random numbers to be in square of side-length = R
	            x = (x * 2.0) - 1.0
	            y = (y * 2.0) - 1.0
	            r = x*x + y*y
	        end do

	        ! calculate the co-efficient to multiply random numbers x and y
	        ! by to achieve normal distribution
	        r = sqrt((-2.0 * log(r)) / r)

	        norm_rand = mean + (std_dev * x * r)
	        spare = y * r
	        has_spare = .TRUE.
	        return
	    end if
	end function norm_rand
    !===================================
    ! subroutine looptx: to perform fwd
    ! modeling for different models
    !===================================
    subroutine looptx
      implicit none
      ! local variables
      integer :: i, iFreq, iTx
	  character :: c


     ! Loop over each transmitter:! 
     do iTx = 1,nTx
         !write(*,'(a24,i6)') 'Transmitter #:',iTx    
         !
         ! Assign Tx parameters:
         !
         xTx1D       = xTxIn(iTx)
         yTx1D       = yTxIn(iTx)
         zTx1D       = zTxIn(iTx)
         azimuthTx1D = azimuthIn(iTx)
         dipTx1D     = dipIn(iTx)
         !
         ! Inner loop of requested frequencies
         !
         do iFreq= 1,nFreq

             ! Get the requested frequency:
             !
             ftx1D = fTxIn(iFreq) 
             !write(*,*) 'Frequency #:',fTxIn(iFreq)            
             !
             ! Compute CSEM fields:
             !           
             call comp_dipole1D
             !   
             ! Output the response for the current transmitter, note that I convert jz to ez here
             ! The assumption is that receivers on a boundary are actually in the top layer.
             ! So a seafloor site would use the sea layer conductivity to convert to vertical
             ! current to electric field.
             ! sigsite is defined for each site in subroutine readrunfile

     
             do i=1,n1D
                 jz1D(i) = jz1D(i)/ (sigsite(i) + II*eps*2*pi*ftx1D)
				 !write(*,*)'sigsite(i)',sigsite(i)
                 !write(*,100) &
                  ! &   real(ex1D(i)),aimag(ex1D(i)),real(ey1D(i)),aimag(ey1D(i)),real(jz1D(i)),aimag(jz1D(i)), &
                  ! &   real(bx1D(i)),aimag(bx1D(i)),real(by1D(i)),aimag(by1D(i)),real(bz1D(i)),aimag(bz1D(i))
                 !write(*,*)'realex1D', real(ex1D(i)), aimag(ex1D(i))
				 !write(*,*)'ey1D',real(ey1D(i)),aimag(ey1D(i))
                 modelresponse(iTx, iFreq) = ey1D(i)    ! *** important: use ey1D because in Kerry's
                 !call fget(c)                          ! code the axis are differnet as compared to 
				 										! emgs data file ****** 
				 
             enddo
			 !call fget(c)
        
            100    format(12(1x,e15.8)) 

            !
            ! Go on to the next frequency
            ! 
         enddo
            !
            ! Go on to the next Tx
            ! 
      enddo 
    end subroutine looptx
	!=======================================
	! subroutine calc_objective(array_indiv)
	!=======================================
	subroutine calc_objective(array_indiv)
		implicit none
		type(individual), dimension(:), allocatable, intent(inout) :: array_indiv
		!local variables
		integer :: i,j,k
		real(8) :: misfit2, roughness, tempphase
		do i = 1, iend - istart + 1
			sig1D(1) = 1.d0/1d12
			sig1D(2) = 3.33
			do j = 1, nparam
				sig1D(j+2) = 10.d0**(array_indiv(i)%model_parameter(j))
			end do
		
	        ! Now lets store the conductivity at the site location for later use.  Sites
	        ! The assumption is that receivers on a boundary are actually in the top layer.
	        ! So a seafloor site would use the sea layer conductivity to convert to vertical
	        ! current to electric field.
	        do j = 1, n1D
	            sigsite(j) = sig1D(1)
	            do k = 2, nlay1D
	                if ( zlay1D(k) .lt. z1D(j) ) then
	                    sigsite(j) = sig1D(k)   
	                end if
	            end do
	        end do
		
			call looptx
			do j = 1, nTx
				do k = 1, nFreq
					data_syn(j,k) = modelresponse(j,k)
				end do
			end do
		
			! defining misfit in terms of log10(amp) & phase(degrees)
			do j = 1, nTx        ! when using freq = 1.0 hz, only first 24 points above 10-15
				do k = 1, nFreq
					logamp_syn(j,k) = log10(abs(data_syn(j,k)))
					tempphase = DATAN2(aimag(data_syn(j,k)),real(data_syn(j,k)))*(180/3.1416)
					!check for phase with unwrapping
					phase_syn(j,k) = tempphase
					if( abs(tempphase + 360d0 - phase_true(j,k)) .le. abs(tempphase - phase_true(j,k)) )phase_syn(j,k) = (tempphase + 360d0)
					if( abs(tempphase - 360d0 - phase_true(j,k)) .le. abs(tempphase - phase_true(j,k)) )phase_syn(j,k) = (tempphase - 360d0)
				end do
			end do
			misfit2  = 0.d0
			do k = 1, nFreq
				if(k == 1)then
					do j = 3, 30
						misfit2 = misfit2 + ((logamp_syn(j,k) - logamp_true(j,k))/stdamp(j,k))**2+ &
						((phase_syn(j,k) - phase_true(j,k))/stdphase(j,k))**2
					end do
				else if(k == 2)then
					do j = 3, 24    ! using this for 1.0 Hz frequency
						misfit2 = misfit2 + ((logamp_syn(j,k) - logamp_true(j,k))/stdamp(j,k))**2+ &
						((phase_syn(j,k) - phase_true(j,k))/stdphase(j,k))**2
					end do
				else
					do j = 3, 10
						misfit2 = misfit2 + ((logamp_syn(j,k) - logamp_true(j,k))/stdamp(j,k))**2+ &
						((phase_syn(j,k) - phase_true(j,k))/stdphase(j,k))**2
					end do
				end if
			end do
			misfit2 = sqrt(misfit2/(2.d0*real(28)))
		
			! copy the misfit2 to obj(1)
			array_indiv(i)%obj(1) = misfit2
		
			!calculate the roughness
			roughness = 0.d0
			do j = 1, nparam - 1
				roughness = roughness + (log10(sig1D(j+3)) - log10(sig1D(j+2)))**2
			end do
			roughness = sqrt(roughness)
			array_indiv(i)%obj(2) = roughness
			
		end do          ! do i = 1, iend - istart + 1
	end subroutine calc_objective
	!=======================================
	!   subroutine send_data_to_master
	!=======================================
	subroutine send_data_to_master(flag)
		implicit none
		integer, intent(in) :: flag
		!local variables
		integer :: i,j, pts
		
		!copy the obj matrix to the temparray_obj_parent
		pts = iend - istart + 1
		select case(flag)
		   case(1)
		    ! send obj to parent_pop_global
			do i = 1, pts
				do j = 1, nobj
					temparray_obj_parent(i,j) = parent_pop(i)%obj(j)
				end do
			end do
			call mpi_send(temparray_obj_parent,pts*nobj,mpi_double_precision,&
			0,worker_sending_data,mpi_comm_world, ierr)
		   case(2)
		    ! send obj to intermediate_pop
			do i = 1, pts
				do j = 1, nobj
					temparray_obj_parent(i,j) = parent_pop(i)%obj(j)
					temparray_obj_child(i,j) = child_pop(i)%obj(j)
				end do
			end do
			call mpi_send(temparray_obj_parent, pts*nobj,mpi_double_precision,&
			0,worker_sending_data,mpi_comm_world,ierr)
			call mpi_send(temparray_obj_child, pts*nobj,mpi_double_precision,&
			0,worker_sending_data,mpi_comm_world, ierr)
			case(3)
			! send model_parameter to intermediate_pop
			do i = 1, pts
				do j = 1, nparam
					temparray_mp_parent(i,j) = parent_pop(i)%model_parameter(j)
					temparray_mp_child(i,j) = child_pop(i)%model_parameter(j)
				end do
			end do
			call mpi_send(temparray_mp_parent,pts*nparam,mpi_double_precision,&
			0,worker_sending_data,mpi_comm_world,ierr)
			call mpi_send(temparray_mp_child,pts*nparam,mpi_double_precision,&
			0,worker_sending_data,mpi_comm_world,ierr)
	    end select
		
	end subroutine send_data_to_master
	!=======================================
	!   subroutine recv_data_from_worker
	!=======================================
	subroutine recv_data_from_worker(i,flag)
		implicit none
		integer, intent(in) :: i
		integer, intent(in) :: flag       
		!local variable
		integer :: beg, fin,j, k
		
		beg = 1 + i*points_per_proc
		fin = (i+1)*points_per_proc
	
		select case(flag)
			case(1)
			 call mpi_recv(temparray_obj_parent, points_per_proc*nobj,mpi_double_precision,&
			 i,worker_sending_data, mpi_comm_world, status, ierr)
			 do k = beg, fin
				 do j = 1, nobj
					 parent_pop_global(k)%obj(j) = temparray_obj_parent(k - i*points_per_proc,j)
				 end do
			 end do
			 
			case(2)
 			 !recive obj in to intermediate_pop
  			 call mpi_recv(temparray_obj_parent,points_per_proc*nobj,mpi_double_precision,i,&
  			 worker_sending_data, mpi_comm_world,status, ierr)
  			 ! receive the child pop also
  			 call mpi_recv(temparray_obj_child,points_per_proc*nobj,mpi_double_precision,i,&
  			 worker_sending_data, mpi_comm_world,status, ierr)
  			 !store in the intermediate_pop array
  			 do k = beg,fin
  			 	do j = 1, nobj
  			 		intermediate_pop(k)%obj(j) = temparray_obj_parent(k - i*points_per_proc,j)
  			 		intermediate_pop(k+popsize)%obj(j) = temparray_obj_child(k - i*points_per_proc,j)
  			 	end do
  			 end do
			 
			case(3)
			 ! recive the model_parameter and save it in intermediate_pop
			 call mpi_recv(temparray_mp_parent, points_per_proc*nparam, mpi_double_precision,i,&
			 worker_sending_data,mpi_comm_world,status, ierr)
			 call mpi_recv(temparray_mp_child, points_per_proc*nparam, mpi_double_precision,i,&
			 worker_sending_data,mpi_comm_world,status, ierr)
			 do k = beg, fin
				 do j = 1, nparam
					 intermediate_pop(k)%model_parameter(j) = temparray_mp_parent(k - i*points_per_proc,j)
					 intermediate_pop(k+popsize)%model_parameter(j) = temparray_mp_child(k - i*points_per_proc,j)
				 end do
			 end do 
	    end select
	end subroutine recv_data_from_worker
	!=======================================
	!   subroutine non_dominated_sort      !
	!=======================================
	subroutine non_dominated_sort(array_indiv,size_array,front)
		implicit none
		type(individual), dimension(:), allocatable, intent(inout) :: array_indiv
		type(front_t), dimension(:), allocatable, intent(inout) :: front
		integer, intent(in) :: size_array
		!local variables
		integer :: i,j,k,m
		integer :: i_dominates_j, j_dominates_i, dom_equal
		integer :: loop_flag
		integer :: counter, front_count
		type(front_t) :: q
			
		allocate(q%member(0))
		allocate(front(0))
		!calculate the first front
		loop_flag = 0
		do while(loop_flag <= 1)
			do i = 1, size_array
				if(loop_flag == 1)then
					allocate(array_indiv(i)%i_dominates_individual(array_indiv(i)%i_dominates_counter))
				end if
				array_indiv(i)%i_dominates_counter = 0         ! counter which sees how many individuals ith individual dominates
				array_indiv(i)%dominates_i_counter = 0         ! counter which tracks how many individuals dominate i
				do j = 1, size_array
					i_dominates_j = 0                          ! used for comparing 'ith' and 'jth' individual with their   
					j_dominates_i = 0                          ! corresponding objective fun values            
					dom_equal = 0
					do k = 1, nobj
						if(array_indiv(i)%obj(k) .lt. array_indiv(j)%obj(k))then
							! idominates j
							i_dominates_j = i_dominates_j + 1
						else if(array_indiv(i)%obj(k) == array_indiv(j)%obj(k))then
							dom_equal = dom_equal + 1
						else
							! j dominates i
							j_dominates_i = j_dominates_i + 1
						end if
					end do  ! do k = 1, nobj

					if((i_dominates_j == 0) .and. (dom_equal /= nobj))then
						! j dominates i, so increment the dominates_i_counter
						array_indiv(i)%dominates_i_counter = array_indiv(i)%dominates_i_counter + 1
					else if((j_dominates_i == 0) .and. (dom_equal /= nobj))then
						! i dominates j, so increment the i_dominates_counter
						array_indiv(i)%i_dominates_counter = array_indiv(i)%i_dominates_counter + 1
						if(loop_flag == 1)then
							array_indiv(i)%i_dominates_individual(array_indiv(i)%i_dominates_counter) = j
						end if  ! if(loop_flag == 1)				
					end if
				end do   ! do j = 1, size_array
			
				! check after comparing with all the models in the population, which have the 
				! dominates_i_counter = 0, assign those individuals the rank 1  
				if(array_indiv(i)%dominates_i_counter == 0 .and. loop_flag == 1)then
					array_indiv(i)%rank = 1
					!debug: want to check how many individuals does this i dominated
					!print *,'i = ',i,'dominates',array_indiv(i)%i_dominates_counter, 'individuals'
					!print *,('k',k,'individual index',array_indiv(i)%i_dominates_individual(k),& 
					!new_line(c),k = 1, array_indiv(i)%i_dominates_counter)
				else
					if(loop_flag == 1)then
						!printing the dominates_i_counter of all the other members in the population
						!print *,'i',i ,'dominated by',array_indiv(i)%dominates_i_counter, 'individuals'
						!print *,'i',i, 'dominates',array_indiv(i)%i_dominates_counter, 'individuals'
						!print *,('k',k,'individual index',array_indiv(i)%i_dominates_individual(k),& 
						!new_line(c),k = 1, array_indiv(i)%i_dominates_counter)
					end if 
				end if   ! if(parent_pop(i)...)
			end do    ! do i = 1, popsize
			loop_flag = loop_flag + 1
		end do    ! do while(loop_flag <= 1)
		
		! add the index of the individuals belonging to the first front
		do i = 1, size_array
			if(array_indiv(i)%rank == 1)then
				q%member = [q%member, i]
			end if
		end do
		front = [front, q]
		
		!generating the subsequent fronts
		front_count = 1
		!print *,'front(',front_count,')%member',front(front_count)%member
		do while(allocated(q%member))
			deallocate(q%member)
			do k = 1, size(front(front_count)%member)
				do j = 1, size(array_indiv(front(front_count)%member(k))%i_dominates_individual)
					! making life easier by assigning of big varible to 'm' variable	
					m = array_indiv(front(front_count)%member(k))%i_dominates_individual(j)
					! decrease this mth individual's counter by 1
					array_indiv(m)%dominates_i_counter = array_indiv(m)%dominates_i_counter -1 
					if(array_indiv(m)%dominates_i_counter == 0)then
						if(.not. allocated(q%member))allocate(q%member(0))
						q%member = [q%member, m]
						!assign the rank also
						array_indiv(m)%rank = front_count + 1
					end if   !if (parent_pop(m)%dominates...
				end do  ! do j = 1, size(.....)
			end do  ! do k = 1, size(front(front_count)%member)
			if(allocated(q%member))then
				front = [front, q]
				front_count = front_count + 1
				!print *,'front(',front_count,')%member',front(front_count)%member
			end if
		end do  ! do while (allocated(q%individual))
		
		print *,'total no of fronts',size(front)
	end subroutine non_dominated_sort
	!=======================================
	!      subroutine non_dominated_sort2
	! this subroutine is intel compiler specific
	! as intel compiler doesn't allows to use the
	! unallocated arrays copying so used the
	! subroutine move_alloc in this
	!=======================================
	subroutine non_dominated_sort2(array_indiv,size_array,front)
		implicit none
		type(individual), dimension(:), allocatable, intent(inout) :: array_indiv
		type(front_t), dimension(:), allocatable, intent(inout) :: front
		integer, intent(in) :: size_array
		!local variables
		integer :: i,j,k,m
		integer :: i_dominates_j, j_dominates_i, dom_equal
		integer :: loop_flag
		integer :: counter, front_count,member_count
		type(front_t) :: q, temp
		type(front_t), dimension(:), allocatable :: temp_front
			
		
		!calculate the first front
		loop_flag = 0
		do while(loop_flag <= 1)
			do i = 1, size_array
				if(loop_flag == 1)then
					allocate(array_indiv(i)%i_dominates_individual(array_indiv(i)%i_dominates_counter))
				end if
				array_indiv(i)%i_dominates_counter = 0         ! counter which sees how many individuals ith individual dominates
				array_indiv(i)%dominates_i_counter = 0         ! counter which tracks how many individuals dominate i
				do j = 1, size_array
					i_dominates_j = 0                          ! used for comparing 'ith' and 'jth' individual with their   
					j_dominates_i = 0                          ! corresponding objective fun values            
					dom_equal = 0
					go to 205
					
					!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					!11/29/18: trying the different sorting method
					! this method is not working. So going down(goto 205)
					! also, in this method models are more rougher
					if(abs(array_indiv(i)%obj(1) - 1d0) .lt. (abs(array_indiv(j)%obj(1) - 1d0)))then
						! i dominates j
						i_dominates_j = i_dominates_j + 1
					else if(abs(array_indiv(i)%obj(1) - 1d0) .gt. (abs(array_indiv(j)%obj(1) - 1d0)))then
						! j dominates i
						j_dominates_i = j_dominates_i + 1
					else
						! both equal
						dom_equal = dom_equal + 1
					end if
				
					if(array_indiv(i)%obj(2) .lt. array_indiv(j)%obj(2))then
						! i dominates j
						i_dominates_j = i_dominates_j + 1
					else if(array_indiv(i)%obj(2) == array_indiv(j)%obj(2))then
						! both equal
						dom_equal = dom_equal + 1
					else
						! j dominates i
						j_dominates_i = j_dominates_i +1
					end if
					
				
					goto 211
					!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	205				do k = 1, nobj	
						if(k == 1)then		
							if((array_indiv(i)%obj(k) .gt. 1d0) .and. (array_indiv(j)%obj(k) .gt. 1d0))then
								! select the one with the least misfit
								if(array_indiv(i)%obj(k) .lt. array_indiv(j)%obj(k))then
									! i dominates j
									i_dominates_j = i_dominates_j + 1
								else if(array_indiv(i)%obj(k) == array_indiv(j)%obj(k))then
									! equal
									dom_equal = dom_equal + 1
								else
									! j dominates i
									j_dominates_i = j_dominates_i + 1
								end if
							else if ((array_indiv(i)%obj(k) .le. 1d0) .and. (array_indiv(j)%obj(k) .gt. 1d0))then
								! i dominates j
								i_dominates_j = i_dominates_j + 1
								
							else if ((array_indiv(i)%obj(k) .gt. 1d0) .and. (array_indiv(j)%obj(k) .le. 1d0))then
								! j dominates i
								j_dominates_i = j_dominates_i + 1
							else
								! both less than one so equal
								dom_equal = dom_equal + 1
							end if
						end if      ! if (k== 1)
						
						if(k == 2)then
							!12/02/18: trying the roughness dominates method when both RMS
							! less than 1 
							! output: doesn't works, made the more rougher models propagate, 
							! so discarded that method
							if((array_indiv(i)%obj(1) .le. 1d0) .and. (array_indiv(j)%obj(1) .le. 1d0))then
								dom_equal = 0
								j_dominates_i = 0
								i_dominates_j = 0
								if(array_indiv(i)%obj(k) .lt. array_indiv(j)%obj(k))then
									! i dominates j
									i_dominates_j = i_dominates_j + 2
								else if(array_indiv(i)%obj(k) == array_indiv(j)%obj(k))then
									! both equal
									dom_equal = dom_equal + 2
								else
									! j dominates i
									j_dominates_i = j_dominates_i + 2
								end if	
							else
								if(array_indiv(i)%obj(k) .lt. array_indiv(j)%obj(k))then
									! i dominates j
									i_dominates_j = i_dominates_j + 1
								else if(array_indiv(i)%obj(k) == array_indiv(j)%obj(k))then
									! both equal
									dom_equal = dom_equal + 1
								else
									! j dominates i
									j_dominates_i = j_dominates_i +1
								end if	
							end if
						end if   ! if(k== 2)then
					end do  ! do k = 1, nobj

	211				if((i_dominates_j == 0) .and. (dom_equal /= nobj))then
						! j dominates i, so increment the dominates_i_counter
						array_indiv(i)%dominates_i_counter = array_indiv(i)%dominates_i_counter + 1
					else if((j_dominates_i == 0) .and. (dom_equal /= nobj))then
						! i dominates j, so increment the i_dominates_counter
						array_indiv(i)%i_dominates_counter = array_indiv(i)%i_dominates_counter + 1
						if(loop_flag == 1)then
							array_indiv(i)%i_dominates_individual(array_indiv(i)%i_dominates_counter) = j
						end if  ! if(loop_flag == 1)				
					end if
				end do   ! do j = 1, size_array
			
				! check after comparing with all the models in the population, which have the 
				! dominates_i_counter = 0, assign those individuals the rank 1  
				if(array_indiv(i)%dominates_i_counter == 0 .and. loop_flag == 1)then
					array_indiv(i)%rank = 1
					!debug: want to check how many individuals does this i dominated
					!print *,'i = ',i,'dominates',array_indiv(i)%i_dominates_counter, 'individuals'
					!print *,('k',k,'individual index',array_indiv(i)%i_dominates_individual(k),& 
					!new_line(c),k = 1, array_indiv(i)%i_dominates_counter)
				else
					if(loop_flag == 1)then
						!printing the dominates_i_counter of all the other members in the population
						!print *,'i',i ,'dominated by',array_indiv(i)%dominates_i_counter, 'individuals'
						!print *,'i',i, 'dominates',array_indiv(i)%i_dominates_counter, 'individuals'
						!print *,('k',k,'individual index',array_indiv(i)%i_dominates_individual(k),& 
						!new_line(c),k = 1, array_indiv(i)%i_dominates_counter)
					end if 
				end if   ! if(parent_pop(i)...)
			end do    ! do i = 1, popsize
			loop_flag = loop_flag + 1
		end do    ! do while(loop_flag <= 1)
		
		! add the index of the individuals belonging to the first front
		member_count = 0
		do i = 1, size_array
			if(array_indiv(i)%rank == 1)then
				member_count = member_count + 1
				allocate(temp%member(1:member_count))
				if(allocated(q%member))temp%member(1:member_count - 1) = q%member
				call move_alloc(temp%member,q%member)
				q%member(member_count) = i
			end if
		end do
		!assigning the temp_front 
		allocate(temp_front(1))
		allocate(temp_front(1)%member(size(q%member)))
		temp_front(1)%member = q%member
		call move_alloc(temp_front, front)
		
		!generating the subsequent fronts
		front_count = 1
		!print *,'front(',front_count,')%member',front(front_count)%member
		do while(allocated(q%member))
			deallocate(q%member)
			member_count = 0
			do k = 1, size(front(front_count)%member)
				do j = 1, size(array_indiv(front(front_count)%member(k))%i_dominates_individual)
					! making life easier by assigning of big varible to 'm' variable	
					m = array_indiv(front(front_count)%member(k))%i_dominates_individual(j)
					! decrease this mth individual's counter by 1
					array_indiv(m)%dominates_i_counter = array_indiv(m)%dominates_i_counter -1 
					
					if(array_indiv(m)%dominates_i_counter == 0)then
						member_count = member_count + 1
						allocate(temp%member(1:member_count))
						if(allocated(q%member))temp%member(1:member_count -1) = q%member
						call move_alloc(temp%member,q%member)
						q%member(member_count) = m
						!assign the rank also
						array_indiv(m)%rank = front_count + 1
					end if   !if (parent_pop(m)%dominates...
				end do  ! do j = 1, size(.....)
			end do  ! do k = 1, size(front(front_count)%member)
			
			if(allocated(q%member))then
				front_count = front_count + 1
				allocate(temp_front(1:front_count))
				if(allocated(front))then
					temp_front(1:front_count -1) = front
				end if
				call move_alloc(temp_front,front)
				front(front_count) = q
			end if
		end do  ! do while (allocated(q%individual))
		
		print *,'total no of fronts',size(front)
	end subroutine non_dominated_sort2
	!=======================================
	! subroutine calc_crowding_dist
	!=======================================
	subroutine calc_crowding_dist(array_indiv, front)
		implicit none
		type(individual),dimension(:),allocatable, intent(inout) :: array_indiv
		type(front_t), dimension(:),allocatable, intent(inout) :: front
		!local variables
		integer :: i,j,k,irow,krow
		real(8),dimension(:,:),allocatable :: temparray
		real(8),dimension(2) :: buf
		real(8) :: obj_max,obj_min, next_obj,prev_obj
		real(8) :: inf
	    IF (ieee_support_inf(inf)) THEN
	      inf = -1.d0*ieee_value(inf,  ieee_negative_inf)
	    END IF
	    
		do i = 1, size(front)
			allocate(temparray(size(front(i)%member),2))
			do j = 1,nobj
				do k = 1, size(front(i)%member)
					temparray(k,1) = array_indiv(front(i)%member(k))%obj(j)
					temparray(k,2) = front(i)%member(k)
				end do     ! do k = 1, size(...)
				!sort the temparray
				do irow = 1, size(front(i)%member)
					krow = minloc(temparray(irow:size(front(i)%member),1),dim = 1) + irow - 1
					buf(:) = temparray(irow,:)
					temparray(irow,:) = temparray(krow,:)
					temparray(krow,:) = buf(:)
				end do     ! do irow = 1,...
				!******debug:*****
				!print *, '*******sorted pop on obj(',j,')************'
				!print *,(temparray(k,2),k = 1, size(front(i)%member))
				!do k = 1, size(front(i)%member)
				!	front(i)%member(k) = temparray(k,2)
				!	print *,'front(i)%member(k)',front(i)%member(k)
				!end do     ! do k = 1, size..
				!****** end debug***
				!calculate the crowding distance 
				obj_min = array_indiv(front(i)%member(1))%obj(j)
				obj_max = array_indiv(front(i)%member(size(front(i)%member)))%obj(j)
				!print *,'obj_max',obj_max,'obj_min',obj_min,'diff',obj_max -obj_min
				array_indiv(front(i)%member(1))%crowd_dist = inf
				array_indiv(front(i)%member(size(front(i)%member)))%crowd_dist = inf
				
				do k = 2, size(front(i)%member) - 1
					next_obj = array_indiv(front(i)%member(k+1))%obj(j)
					prev_obj = array_indiv(front(i)%member(k-1))%obj(j)
					if(obj_max - obj_min == 0)then
						array_indiv(front(i)%member(k))%crowd_dist = inf
					else
						array_indiv(front(i)%member(k))%crowd_dist = array_indiv(front(i)%member(k))%crowd_dist + &
						(next_obj - prev_obj)/(obj_max - obj_min)
					end if  ! if(obj_max - obj_min == 0)
				end do   ! do k = 2, size(front(i)%member) - 1
				!write(*,*)(array_indiv(front(i)%member(k))%obj(j),array_indiv(front(i)%member(k))%crowd_dist,&
				!new_line(c),k = 1, size(front(i)%member))
			end do         ! do j = 1, nobj
			deallocate(temparray)
			!call fget(c)
		end do             ! do i = 1, size(front)
	end subroutine calc_crowding_dist
	!=======================================
	!     subroutine calc_crowding_dist2
	! this one calculates the crowding dist
	! in the model space
	!=======================================
	subroutine calc_crowding_dist2(array_indiv, front)
		implicit none
		type(individual),dimension(:),allocatable, intent(inout) :: array_indiv
		type(front_t), dimension(:),allocatable, intent(inout) :: front
		!local variables
		integer :: i,j,k,irow,krow
		real(8),dimension(:,:),allocatable :: temparray
		real(8),dimension(2) :: buf
		real(8) :: param_max,param_min, next_param,prev_param
		real(8) :: inf
	    IF (ieee_support_inf(inf)) THEN
	      inf = -1.d0*ieee_value(inf,  ieee_negative_inf)
	    END IF
	    
		do i = 1, size(front)
			allocate(temparray(size(front(i)%member),2))
			do j = 1,nparam
				do k = 1, size(front(i)%member)
					temparray(k,1) = array_indiv(front(i)%member(k))%model_parameter(j)
					temparray(k,2) = front(i)%member(k)
				end do     ! do k = 1, size(...)
				!sort the temparray
				do irow = 1, size(front(i)%member)
					krow = minloc(temparray(irow:size(front(i)%member),1),dim = 1) + irow - 1
					buf(:) = temparray(irow,:)
					temparray(irow,:) = temparray(krow,:)
					temparray(krow,:) = buf(:)
				end do     ! do irow = 1,...
				!******debug:*****
				!print *, '*******sorted pop on model_parameter(',j,')************'
				!print *,(temparray(k,2),k = 1, size(front(i)%member))
				!do k = 1, size(front(i)%member)
				!	front(i)%member(k) = temparray(k,2)
				!	print *,'front(i)%member(k)',front(i)%member(k)
				!end do     ! do k = 1, size..
				!****** end debug***
				!calculate the crowding distance 
				param_min = array_indiv(front(i)%member(1))%model_parameter(j)
				param_max = array_indiv(front(i)%member(size(front(i)%member)))%model_parameter(j)
				!print *,'obj_max',obj_max,'obj_min',obj_min,'diff',obj_max -obj_min
				array_indiv(front(i)%member(1))%crowd_dist = inf
				array_indiv(front(i)%member(size(front(i)%member)))%crowd_dist = inf
				
				do k = 2, size(front(i)%member) - 1
					next_param = array_indiv(front(i)%member(k+1))%model_parameter(j)
					prev_param = array_indiv(front(i)%member(k-1))%model_parameter(j)
					if(param_max - param_min == 0)then
						array_indiv(front(i)%member(k))%crowd_dist = inf
					else
						array_indiv(front(i)%member(k))%crowd_dist = array_indiv(front(i)%member(k))%crowd_dist + &
						(next_param - prev_param)/(param_max - param_min)
					end if  ! if(param_max - param_min == 0)
				end do   ! do k = 2, size(front(i)%member) - 1
				!write(*,*)(array_indiv(front(i)%member(k))%model_parameter(j),array_indiv(front(i)%member(k))%crowd_dist,&
				!new_line(c),k = 1, size(front(i)%member))
			end do         ! do j = 1, nparam
			deallocate(temparray)
			!call fget(c)
		end do             ! do i = 1, size(front)
	end subroutine calc_crowding_dist2
	!=======================================
	!   subroutine tourn_select
	!=======================================
	subroutine tourn_select
		implicit none
		integer :: i,j,k
		real(8) :: rtemp
		type(individual) :: winner, pick
		
		!we have to repeat the same process untill
		! the mating_pool is full
		do k = 1, mating_pool
			do i = 1, tournsize
				call random_number(rtemp)
				j = 1 + (iend-istart+1)*rtemp
				tourn_pop(i) = parent_pop(j)
			end do          ! do i = 1, tournsize
			
			!making sure that the same candidates are not chosen in tourn_pop
			do i = 1, tournsize
				do while(any(tourn_pop(1:i-1)%index == tourn_pop(i)%index))
					call random_number(rtemp)
					j = 1 + (iend-istart+1)*rtemp
					tourn_pop(i) = parent_pop(j)
				end do      ! do while
			end do          ! do i = 1, tournsize
			
			!compare the rank of the tourn_pop members
			winner = tourn_pop(1)
			do i = 2, tournsize
				pick = tourn_pop(i)
				!compare on the basis of rank firstly
				if(pick%rank .lt. winner%rank)then
					winner = pick
				elseif(pick%rank == winner%rank)then
					!compare the crowding distance
					if(pick%crowd_dist .gt. winner%crowd_dist)then
						winner = pick
					else if(pick%crowd_dist == winner%crowd_dist)then
						!both have same rank and same distance then 
						!randomly select either one of them
						call random_number(rtemp)
						if(rtemp .lt. 0.5)winner = pick
					end if  ! if(pick%crowd_dist == winner%crowd_dist)
				end if      ! if(pick%rank .lt. winner%rank)
			end do          ! do i = 2, tournsize
			
			! assign the winner to mating pool
			mating_pop(k) = winner
		end do              ! do k = 1, mating_pool
		
		!making sure that the mating pop also doesn't has same candidates
		do k = 1, mating_pool
			do while(any(mating_pop(1:k-1)%index == mating_pop(k)%index))
				do i = 1, tournsize
					call random_number(rtemp)
					j = 1 + (iend-istart+1)*rtemp
					tourn_pop(i) = parent_pop(j)
				end do          ! do i = 1, tournsize
				!making sure that the same candidates are not chosen
				do i = 1, tournsize
					do while(any(tourn_pop(1:i-1)%index == tourn_pop(i)%index))
						call random_number(rtemp)
						j = 1 + (iend-istart+1)*rtemp
						tourn_pop(i) = parent_pop(j)
					end do      ! do while
				end do          ! do i = 1, tournsize
			
				!compare the rank of the tourn_pop members
				winner = tourn_pop(1)
				do i = 2, tournsize
					pick = tourn_pop(i)
					!compare on the basis of rank firstly
					if(pick%rank .lt. winner%rank)then
						winner = pick
					elseif(pick%rank == winner%rank)then
						!compare the crowding distance
						if(pick%crowd_dist .gt. winner%crowd_dist)then
							winner = pick
						else if(pick%crowd_dist == winner%crowd_dist)then
							!both have same rank and same distance then 
							!randomly select either one of them
							call random_number(rtemp)
							if(rtemp .lt. 0.5)winner = pick
						end if  ! if(pick%crowd_dist == winner%crowd_dist)
					end if      ! if(pick%rank .lt. winner%rank)
				end do          ! do i = 2, tournsize
			
				! assign the winner to mating pool
				mating_pop(k) = winner
			end do          ! do while
		end do              ! do k = 1, mating_pool
	end subroutine tourn_select
	!=======================================
	!    subroutine tourn_select2
	!=======================================
	subroutine tourn_select2
		implicit none
		integer :: i,j,k
		real(8) :: rtemp
		type(individual) :: winner, pick
		
		!we have to repeat the same process untill
		! the mating_pool is full
		do k = 1, mating_pool
			do i = 1, tournsize
				call random_number(rtemp)
				j = 1 + (iend-istart+1)*rtemp
				tourn_pop(i) = parent_pop(j)
			end do          ! do i = 1, tournsize
			
			!making sure that the same candidates are not chosen in tourn_pop
			do i = 1, tournsize
				do while(any(tourn_pop(1:i-1)%index == tourn_pop(i)%index))
					call random_number(rtemp)
					j = 1 + (iend-istart+1)*rtemp
					tourn_pop(i) = parent_pop(j)
				end do      ! do while
			end do          ! do i = 1, tournsize
			
			! compare the misfit of the tourn_pop members
			winner = tourn_pop(1)
			do i = 2, tournsize
				pick = tourn_pop(i)
				!compare on the basis of misfit firstly
				if((pick%obj(1) .gt. 1d0) .and. (winner%obj(1) .gt. 1d0))then
					if(pick%obj(1) .lt. winner%obj(1))then
						winner = pick
					end if
				else if((pick%obj(1) .le. 1d0) .and. (winner%obj(1) .gt. 1d0))then
					winner = pick
				else if((pick%obj(1) .le. 1d0) .and. (winner%obj(1) .le. 1d0))then
					! compare roughness
					if(pick%obj(2) .lt. winner%obj(2))then
						winner = pick
					end if
				end if
			end do          ! do i = 2, tournsize
			
			! assign the winner to mating pool
			mating_pop(k) = winner
		end do              ! do k = 1, mating_pool
		
		!making sure that the mating pop also doesn't has same candidates
		do k = 1, mating_pool
			do while(any(mating_pop(1:k-1)%index == mating_pop(k)%index))
				do i = 1, tournsize
					call random_number(rtemp)
					j = 1 + (iend-istart+1)*rtemp
					tourn_pop(i) = parent_pop(j)
				end do          ! do i = 1, tournsize
				!making sure that the same candidates are not chosen
				do i = 1, tournsize
					do while(any(tourn_pop(1:i-1)%index == tourn_pop(i)%index))
						call random_number(rtemp)
						j = 1 + (iend-istart+1)*rtemp
						tourn_pop(i) = parent_pop(j)
					end do      ! do while
				end do          ! do i = 1, tournsize
			
				! compare the misfit of the tourn_pop members
				winner = tourn_pop(1)
				do i = 2, tournsize
					pick = tourn_pop(i)
					!compare on the basis of misfit firstly
					if((pick%obj(1) .gt. 1d0) .and. (winner%obj(1) .gt. 1d0))then
						if(pick%obj(1) .lt. winner%obj(1))then
							winner = pick
						end if
					else if((pick%obj(1) .le. 1d0) .and. (winner%obj(1) .gt. 1d0))then
						winner = pick
					else if((pick%obj(1) .le. 1d0) .and. (winner%obj(1) .le. 1d0))then
						! compare roughness
						if(pick%obj(2) .lt. winner%obj(2))then
							winner = pick
						end if
					end if
				end do          ! do i = 2, tournsize
			
				! assign the winner to mating pool
				mating_pop(k) = winner
			end do          ! do while
		end do              ! do k = 1, mating_pool
		
	end subroutine tourn_select2
	!=======================================
	!   subroutine cross_over
	!=======================================
	subroutine cross_over(par1,par2,child1,child2)			
		implicit none
		type(individual),intent(inout) :: par1, par2, child1, child2
		!local variables                                    
		real(8) :: difference, x_mean, beta, u,rtemp
		integer :: site, k, x_s
		x_s = 0
		call random_number(rtemp)
		if(rtemp .lt. pcross)then
			do site = 1, nparam
				call random_number(rtemp)
				if(rtemp .le. 0.5)then
					call create_children(par1%model_parameter(site), &
					par2%model_parameter(site),child1%model_parameter(site),&
					child2%model_parameter(site),parmin(site),parmax(site))
				else
					! copy the corresponding parameters as it is
					child1%model_parameter(site) = par1%model_parameter(site)
					child2%model_parameter(site) = par2%model_parameter(site)		
				end if          ! if(rtemp .le. 0.5)
			end do              ! do site = 1, nparam
		else
			! copy the whole parameter as it is 
			child1 = par1
			child2 = par2
		end if                  ! if(rtemp .lt. pcross)
	end subroutine cross_over
	
	!=======================================
	!     subroutine create_children
	!=======================================
	subroutine create_children(p1,p2,c1,c2,low, high)
		implicit none
		real(8),intent(inout) :: p1, p2
		real(8),intent(out) :: c1, c2
		real(8),intent(in) :: low,high
		!local variables
		real(8) ::  difference, x_mean, beta,alpha,distance
		real(8) :: randvar,temp,rtemp, umax
		integer :: flag
		real(8) :: inf
	    IF (ieee_support_inf(inf)) THEN
	      inf = -1.d0*ieee_value(inf,  ieee_negative_inf)
	    END IF
		flag = 0
		if(p1 .gt. p2)then
			temp = p1
			p1 = p2
			p2 = temp
			flag = 1
		end if
		x_mean = (p1 + p2)*0.5
		difference = p2 - p1
		if((p1 - low) .lt. (high - p2))then
			distance = p1 - low
		else 
			distance = high - p2
		end if
		
		if(distance .lt. 0)then
			distance = 0.0
		end if
		
		if(RIGID .and. (difference .gt. EPSILON))then
			alpha = 1.0 + (2.0 *distance/difference)
			umax = 1.0 - (0.5/(alpha**(n_distribution_c + 1.0)))
			call random_number(rtemp)
			randvar = umax*rtemp
		else
			call random_number(randvar)
		end if    				! if(RIGID)
		call get_beta(randvar,beta)
		if(abs(difference*beta) .gt. inf)then
			beta = inf/difference
		end if
		c2 = x_mean + beta*0.5*difference
		c1 = x_mean - beta*0.5*difference
		
		if(c2 < low)then
			c2 = low
		end if
		if(c2 > high)then
			c2 = high
		end if
		
		if(c1 > high)then
			c1 = high
		end if
		if(c2 < low)then
			c2 = low
		end if
		if(flag ==1)then
			temp = c1
			c1 = c2
			c2 = temp
		end if	
	end subroutine create_children
	!=======================================
	!      subroutine get_beta
	!=======================================
	subroutine get_beta(u,beta)
		implicit none
		real(8), intent(out) :: beta
		real(8),intent(inout) :: u
		if(1.0 - u .lt. EPSILON)then
			u = 1.0 - EPSILON
		end if
		if(u .lt. 0)then
			u = 0.0
		end if
		if(u .lt. 0.5)then
			beta = (2.0*u)**(1.0/(n_distribution_c + 1.0))
		else 
			beta = (0.5/(1.0 - u))**(1.0/(n_distribution_c + 1.0))
		end if
	end subroutine get_beta
	!=======================================
	!      subroutine mutation
	!=======================================
	subroutine mutation(par)
		implicit none
		type(individual),intent(inout) :: par
		!local variables
		real(8) :: distance, rtemp,distance1
		real(8) :: delta_l, delta_u, u,x, delta
		integer :: i
		type(individual) :: temp_indiv
		!allocate(par)
		temp_indiv  = par
		do i = 1, nparam
			call random_number(rtemp)
			if(rtemp .lt. pmutate)then
				!x = par%model_parameter(i)
				x = temp_indiv%model_parameter(i)
				distance1 = parmin(i) - x
				delta_l = distance1/(parmax(i) - parmin(i))
				if(delta_l .lt. -1.0)then
					delta_l = -1.0
				end if
				distance1 = parmax(i) - x
				delta_u = distance1/(parmax(i) - parmin(i))
				if(delta_u .gt. 1.0)then
					delta_u = 1.0
				endif
				
				if(-1.0*delta_l .lt. delta_u)then
					delta_u = -1.0*delta_l
				else
					delta_l = -1.0*delta_u
				endif
				call random_number(rtemp)
				u = rtemp
				! calculation of the actual delta value
				call get_delta(u,delta_l, delta_u, delta)
				x = x + delta*(parmax(i) - parmin(i))
				par%model_parameter(i) = x
			endif           ! if(rtemp .lt. pmutate)
		end do 				! do i = 1, nparam
		
	end subroutine mutation
		
	!=======================================
	! subroutine get_delta
	!=======================================
	subroutine get_delta(u,delta_l, delta_u,delta)
		implicit none
		real(8),intent(inout) ::u
		real(8),intent(in) :: delta_l, delta_u
		real(8),intent(out) :: delta
		!local variables
		real(8) :: aa
		if(u .ge. 1.0 - 1.0e-9)then
			delta = delta_u
		else if( u .le. 0.0+ 1.0e-9)then
			delta = delta_l
		else
			if(u .le. 0.5)then
				aa = 2.0*u + (1.0-2.0*u)*(1+delta_l)**(n_distribution_m + 1.0)
				delta = aa**(1.0/(n_distribution_m + 1.0)) - 1.0
			else
				aa = 2.0*(1-u) + 2.0*(u-0.5)*(1-delta_u)**(n_distribution_m + 1.0)
				delta = 1.0 - aa**(1.0/(n_distribution_m + 1.0))
			end if
		end if 
		if(delta <-1.0 .or. delta > 1.0)then
			print *,'error in mutation'
			!stop
		end if
	end subroutine get_delta
	!=======================================
	!       subroutine replace_pop
	!=======================================
	subroutine replace_pop(front,intermediate_pop,parent_pop_global)
		implicit none
		type(front_t), dimension(:),allocatable, intent(in) :: front
		type(individual),allocatable, dimension(:), intent(in) :: intermediate_pop
		type(individual),allocatable, dimension(:),intent(inout) :: parent_pop_global
		! local variables
		integer :: pop_count,front_counter,member_count,i,irow,krow,members_need
		real(8),dimension(:,:),allocatable :: temparray
		real(8), dimension(2) :: buf
		pop_count = 0
		front_counter = 1
		do while(pop_count .lt. popsize)
			!get the members of the ith front
			member_count = size(front(front_counter)%member)
			if(member_count + pop_count .lt. popsize)then
				! add all the members of this front 
				do i = 1, member_count
					pop_count = pop_count + 1
					parent_pop_global(pop_count) = intermediate_pop(front(front_counter)%member(i))
				end do        
			else
				! calculate how many members you want from this front
				members_need = popsize - pop_count
				! sort the members on the basis of crowding distance
				allocate(temparray(member_count,2))
				do i = 1, member_count
					temparray(i,1) = intermediate_pop(front(front_counter)%member(i))%crowd_dist
					temparray(i,2) = front(front_counter)%member(i)
				end do          ! do i = 1, size(....)
				! sort the temparray with maximum crowd_distance on top
				do irow = 1, size(front(front_counter)%member)
					krow = maxloc(temparray(irow:size(front(front_counter)%member),1),dim =1) + irow -1
					buf(:) = temparray(irow,:)
					temparray(irow,:) = temparray(krow,:)
					temparray(krow,:) = buf(:)
				end do      ! do irow
				! now take the members from the temparray and enter them to the parent_pop
				do i = 1, members_need
					pop_count = pop_count + 1
					parent_pop_global(pop_count) = intermediate_pop(int(temparray(i,2)))
				end do          ! do i = 1, members_need
				deallocate(temparray)
			end if              ! if(member_count + pop_count .lt. popsize)
			front_counter = front_counter + 1
		end do                  ! do while(pop_count .lt. popsize)
	end subroutine replace_pop
	!=======================================
	!       subroutine replace_pop2
	!=======================================
	subroutine replace_pop2(front,intermediate_pop,parent_pop_global)
		implicit none
		type(front_t), dimension(:),allocatable, intent(in) :: front
		type(individual),allocatable, dimension(:), intent(in) :: intermediate_pop
		type(individual),allocatable, dimension(:),intent(inout) :: parent_pop_global
		! local variables
		integer :: pop_count,front_counter,member_count,i,irow,krow,members_need
		real(8),dimension(:,:),allocatable :: temparray
		real(8), dimension(2) :: buf
		pop_count = 0
		front_counter = 1
		do while(pop_count .lt. popsize)
			!get the members of the ith front
			member_count = size(front(front_counter)%member)
			if(member_count + pop_count .lt. popsize)then
				! add all the members of this front 
				do i = 1, member_count
					pop_count = pop_count + 1
					parent_pop_global(pop_count) = intermediate_pop(front(front_counter)%member(i))
				end do        
			else
				! calculate how many members you want from this front
				members_need = popsize - pop_count
				! sort the members on the basis of crowding distance
				!********
				
				!********
				allocate(temparray(member_count,2))
				do i = 1, member_count
					temparray(i,1) = intermediate_pop(front(front_counter)%member(i))%crowd_dist
					temparray(i,2) = front(front_counter)%member(i)
				end do          ! do i = 1, size(....)
				! sort the temparray with maximum crowd_distance on top
				do irow = 1, size(front(front_counter)%member)
					krow = maxloc(temparray(irow:size(front(front_counter)%member),1),dim =1) + irow -1
					buf(:) = temparray(irow,:)
					temparray(irow,:) = temparray(krow,:)
					temparray(krow,:) = buf(:)
				end do      ! do irow
				! now take the members from the temparray and enter them to the parent_pop
				do i = 1, members_need
					pop_count = pop_count + 1
					parent_pop_global(pop_count) = intermediate_pop(int(temparray(i,2)))
				end do          ! do i = 1, members_need
				deallocate(temparray)
			end if              ! if(member_count + pop_count .lt. popsize)
			front_counter = front_counter + 1
		end do                  ! do while(pop_count .lt. popsize)
		
	end subroutine replace_pop2
	!=======================================
	!       subroutine send_data_to_worker
	!=======================================
	subroutine send_data_to_worker(i)
		implicit none
		integer, intent(in) :: i        !this is procnum
		!local variables
		integer :: beg, fin, k, j
		
		beg = 1 + i*points_per_proc
		fin = (i+1)*points_per_proc
		
		do k = beg, fin
			do j = 1, nparam
				temparray_mp_parent(k-i*points_per_proc,j) = parent_pop_global(k)%model_parameter(j)
			end do
			do j = 1, nobj
				temparray_obj_parent(k-i*points_per_proc,j) = parent_pop_global(k)%obj(j)
			end do
		end do
		call mpi_send(temparray_mp_parent,points_per_proc*nparam,mpi_double_precision,&
		i,master_sending_data, mpi_comm_world,ierr)
		call mpi_send(temparray_obj_parent,points_per_proc*nobj,mpi_double_precision,&
		i,master_sending_data, mpi_comm_world, ierr)
	end subroutine send_data_to_worker
	!=======================================
	!   subroutine recv_data_from_master
	!=======================================
	subroutine recv_data_from_master
		implicit none
		!local variables
		integer :: i, j
		
		call mpi_recv(temparray_mp_parent,points_per_proc*nparam,mpi_double_precision,&
		0, master_sending_data, mpi_comm_world,status,ierr)
		call mpi_recv(temparray_obj_parent, points_per_proc*nobj,mpi_double_precision,&
		0,master_sending_data, mpi_comm_world, status, ierr)
		
		!copy it to the parent_pop%model_parameter
		do i = 1, iend -istart + 1
			do j = 1, nparam
				parent_pop(i)%model_parameter(j) = temparray_mp_parent(i,j)
			end do
			do j = 1, nobj
				parent_pop(i)%obj(j) = temparray_obj_parent(i,j)
			end do
		end do
	end subroutine recv_data_from_master
	!=======================================
	!   subroutine acceptable_models
	! this subroutine is used to store all the 
	! models which have the misfit value of 1.
	! It's very messy for now.
	!=======================================
	subroutine calc_acceptable_models
		implicit none
		!local variables
		integer :: i,j,k, member_count
		type(individual), dimension(:), allocatable :: temp_acceptable,second_temp_acceptable
		
		! first count the accepatable models
		member_count = 0
		do i = 1, popsize
			if((parent_pop_global(i)%obj(1) .le. 1.05d0) .and. (parent_pop_global(i)%obj(1) .ge. 0.95d0))then
				member_count = member_count + 1
			end if
		end do
		!print *,'member_count',member_count
		
		! again loop over the popsize to fill the temp_acceptable, reinitialize the counter
		! but only if member_count .gt. 0
		if(member_count .gt. 0)then
			! now allocate a temp_model array which stores all these models
			allocate(temp_acceptable(member_count))
			member_count = 0
			do i = 1, popsize
				if((parent_pop_global(i)%obj(1) .le. 1.05d0) .and. (parent_pop_global(i)%obj(1) .ge. 0.95d0))then
					member_count = member_count + 1
					temp_acceptable(member_count) = parent_pop_global(i)
				end if
			end do
		end if
		
		! now copy this temp_acceptable to the acceptable_models 
		if(allocated(acceptable_models))then
			! now again reallocate the acceptables model array with the increased size
			acceptable_member_count = acceptable_member_count + member_count
			!print *,'acceptable_member_count',acceptable_member_count
			if(acceptable_member_count > 0)then
				allocate(second_temp_acceptable(acceptable_member_count))
				if(allocated(acceptable_models))second_temp_acceptable(1:acceptable_member_count - member_count) = acceptable_models
				call move_alloc(second_temp_acceptable,acceptable_models)
			end if
			
			do j = 1, member_count   
				acceptable_models(acceptable_member_count - member_count + j) = temp_acceptable(j)
			end do
		else
			if(member_count .gt. 0)then
				acceptable_member_count = member_count
				allocate(acceptable_models(acceptable_member_count))
				acceptable_models = temp_acceptable
			end if
		end if
		if(allocated(temp_acceptable))deallocate(temp_acceptable)
		
	end subroutine calc_acceptable_models
end module nsga
