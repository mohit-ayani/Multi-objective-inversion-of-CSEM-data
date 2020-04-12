! This is the main program which implements the 
! parallel NSGA II on marine CSEM data
! Author: Mohit Ayani(mohitayani@gmail.com)
! Dated: 06/19/18
!***************************************
program main
	use nsga
	
	implicit none
	!local variables
	integer :: i, j, flag, child_counter
	real(8) :: rtemp
	logical :: exist
	character*12 :: filename
	! initialize mpi routines
	call mpi_init(ierr)
	call mpi_comm_size(mpi_comm_world, nproc, ierr)
	call mpi_comm_rank(mpi_comm_world, procnum, ierr)
	call mpi_barrier(mpi_comm_world, ierr)
	
	! start the timer
	start_time = mpi_Wtime()
	
	! initialize the population
	call init_defaults_Dipole1D                    ! defined in dipole1D.f90
	call init_defaults                             ! defined in modnsga_parallel.f90
	
	! take turns to read RUNFILE and calling CallDipole1D
	if(procnum ==0)then
		call initialize_ga
		call CallDipole1D

		! send the message to next proc to read the file
		if(nproc > 1)then
			call mpi_send(mpi_bottom,0,mpi_integer,1,4,&
			mpi_comm_world, ierr)
		endif
	else
		! receive the message to go ahead & call initialize
		call mpi_recv(mpi_bottom,0,mpi_integer,procnum-1,4,&
		mpi_comm_world,status,ierr)
		call initialize_ga
		call CallDipole1D
		if(procnum < nproc-1)then
			call mpi_send(mpi_bottom, 0, mpi_integer,procnum + 1,4,&
			mpi_comm_world, ierr)
		end if
	end if                  ! if (procnum == 0)
	
	! initalize em data and add random_noise to it
	call initialize_em
	if(procnum == 0)then
		open(unit = 34, file = 'data_noise.txt',status = 'unknown', action = 'write')
		do j = 1, nFreq
			do i = 1, nTx
				print *,data_true(i,j),'i',i,'j',j
				write(34,*)logamp_true(i,j), phase_true(i,j)
			end do
		end do
		close(34)
	end if
	! when we want to invert the anisotorpic data using anisotropic inversion
	! then use the file data_for_inversion to read it first and store it in
	!if(procnum == 0)then
	!	open(unit = 53, file = 'data_for_inversion.txt',status = 'unknown',action = 'read')
	!	do i = 1, 40
	!		read(53,*)data_true(i,1), stdamp(i,1),stdphase(i,1)
	!		print *, stdamp(i,1),stdphase(i,1)
	!	end do
	!	close(53)
	!	! send the message for other procs to read the file
	!	if(nproc > 1)then
	!		call mpi_send(mpi_bottom,0,mpi_integer,1,4,&
	!		mpi_comm_world, ierr)
	!	endif
	!else
	!	! receive the message to go ahead & call initialize
	!	call mpi_recv(mpi_bottom,0,mpi_integer,procnum-1,4,&
	!	mpi_comm_world,status,ierr)
	!	
	!	open(unit = 53, file = 'data_for_inversion.txt',status = 'unknown',action = 'read')
	!	do i = 1, 40
	!		read(53,*)data_true(i,1), stdamp(i,1),stdphase(i,1)
	!	end do
	!	close(53)
	!	if(procnum < nproc-1)then
	!		call mpi_send(mpi_bottom, 0, mpi_integer,procnum + 1,4,&
	!		mpi_comm_world, ierr)
	!	end if
	!end if 
	
	! evaluate the objective funcn
	call calc_objective(parent_pop)
	call mpi_barrier(mpi_comm_world, ierr)
	
	! send obj matrix to master node
	if(procnum /= 0)then
		flag = 1            ! 1 is for sending to parent_pop_global
		call send_data_to_master(flag)
	end if
	if(procnum == 0)then
		flag = 1
		do i = 1, iend - istart + 1
			do j = 1, nobj
				parent_pop_global(i)%obj(j) = parent_pop(i)%obj(j)
			end do
		end do
		
		do i = 1, nproc - 1
			call recv_data_from_worker(i, flag)
		end do
		
		! perform non-dominated sorting
		call non_dominated_sort2(parent_pop_global,popsize,front_parent)
		
		! calculate the crowding distance
		call calc_crowding_dist(parent_pop_global, front_parent)
		
		!call calc_crowding_dist2(parent_pop_global, front_parent)
		!store the rank and crowd-distance in the temparrays
		do i = 1, popsize
			tempcrowdist_global(i) = parent_pop_global(i)%crowd_dist
			temprank_global(i) = parent_pop_global(i)%rank
		end do
		
		! generate count and disp arrays for the slaves
		do i = 1, nproc
			beg = 1 + (i-1)*points_per_proc
			fin = min(points_per_proc + (i-1)*points_per_proc,popsize)
			counts_masterproc(i) = fin - beg + 1
			disp_masterproc(i) = (i-1)*points_per_proc
		end do
	end if                  ! if(procnum == 0)
	
	call mpi_bcast(counts_masterproc, nproc,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_bcast(disp_masterproc, nproc, mpi_integer,0,mpi_comm_world,ierr)
	mycounts = counts_masterproc(procnum + 1)
	
	!use scatterv to distribute ranks and crowd dist
	call mpi_scatterv(tempcrowdist_global,counts_masterproc,disp_masterproc,&
	mpi_double_precision,tempcrowdist_local,mycounts,mpi_double_precision,0,&
	mpi_comm_world,ierr)
	call mpi_scatterv(temprank_global,counts_masterproc,disp_masterproc,&
	mpi_integer,temprank_local, mycounts, mpi_integer,0,mpi_comm_world,ierr)
	
	!copy these localrank and localdist to the parent_pop%dist and parent_pop
	!rank
	parent_pop(1:iend-istart+1)%rank = temprank_local
	parent_pop(1:iend-istart+1)%crowd_dist = tempcrowdist_local
	
	!*** gen loop starts here
	do igen = 1, maxgen
		call tourn_select
		
		!perform cross-over and mutation
		child_counter = 1                       ! keeps track of child_pop
		do i = 1, mating_pool
			call random_number(rtemp)
			j = 1 + rtemp*mating_pool
			par1 = mating_pop(j)
			call random_number(rtemp)
			j = 1 + rtemp*mating_pool
			par2 = mating_pop(j)
			!make sure that the same parent is not selected
			do while(par1%index == par2%index)
				call random_number(rtemp)
				j = 1 + rtemp*mating_pool
				par2 = mating_pop(j)
			end do
			!** using iteration adaptive crossover rate
			pcross = 0.7 - 0.1*real(igen/maxgen)
		    n_distribution_c = 1d0 + 19d0*real(igen/maxgen)
			call cross_over(par1,par2,child1,child2)
			!pmutate = 1.d0/real(nparam) + real(igen/maxgen)*(1.d0 - 1.d0/real(nparam)) 
			pmutate = 0.01
			n_distribution_m = 100d0 + real(igen)
			call mutation(child1)
			call mutation(child2)
			child_pop(child_counter) = child1
			child_pop(child_counter + 1) = child2
			child_counter = child_counter + 2
		end do                  ! do i = 1, mating_pool
		
		! calculate the objective value of child_pop
		call calc_objective(child_pop)
		call mpi_barrier(mpi_comm_world,ierr)
		
		!perform elitism
		! workers send their obj matrix to master node
		if(procnum /= 0)then
			flag = 2    ! sending to the intermediate pop
			call send_data_to_master(flag)
		end if
		if(procnum == 0)then
			flag = 2
			do i = 1, iend - istart + 1
				do j = 1, nobj
					intermediate_pop(i)%obj(j) = parent_pop(i)%obj(j)
					intermediate_pop(i+popsize)%obj(j) = child_pop(i)%obj(j)
				end do
			end do
			do i = 1, nproc - 1
				call recv_data_from_worker(i,flag)
			end do
		end if
		call mpi_barrier(mpi_comm_world,ierr)
		
		if(procnum == 0)then
			!since the parent_pop(i)%i_dominates_individual are already allocated which 
			! in turn will also allocate the child_pop(i)%i_dominates_individual if there
			! is no cross-over and the parents are copied as it is to the child_pop.
			!thus, we have to deallocate the intermediate_pop(i)%i_dominates_individual
			! for its population before calling the subroutine non_dominated_sort. 
			do i = 1, 2*popsize
				if(allocated(intermediate_pop(i)%i_dominates_individual))then
					deallocate(intermediate_pop(i)%i_dominates_individual)
				end if
			end do 	
	
			! do non-dominated sorting of the intermediate pop
			call non_dominated_sort2(intermediate_pop, 2*popsize,front_intermediate)
	
			!calculate the crowding distance
			call calc_crowding_dist(intermediate_pop,front_intermediate)
			!call calc_crowding_dist2(intermediate_pop, front_intermediate)
		end if	            ! if(procnum == 0)
		if(procnum /= 0)then
			flag = 3        ! 3 is for sending model_parameter to master
			call send_data_to_master(flag)
		end if

		if(procnum == 0)then
			flag = 3
			do i = 1, iend - istart + 1
				do j = 1, nparam
					intermediate_pop(i)%model_parameter(j) = parent_pop(i)%model_parameter(j)
					intermediate_pop(i+popsize)%model_parameter(j) = child_pop(i)%model_parameter(j)
				end do
			end do
			do i = 1, nproc - 1
				call recv_data_from_worker(i, flag)
			end do
		end if
		call mpi_barrier(mpi_comm_world,ierr)
	
		if(procnum == 0)then
	
			! generate the parent_pop for next gen
			call replace_pop(front_intermediate, intermediate_pop,parent_pop_global)
	
			! distribute the parent_pop to the workers
			!a: copy on its own node
			do i = 1, iend - istart + 1
				parent_pop(i) = parent_pop_global(i)
			end do
	
			!b: send to the workers
			do i = 1, nproc - 1
				call send_data_to_worker(i)
			end do
			do i = 1, popsize
				tempcrowdist_global(i) = parent_pop_global(i)%crowd_dist
				temprank_global(i) = parent_pop_global(i)%rank
			end do
		end if                      ! if(procnum == 0)
		
		! slaves send their model parameter array to master proc
		if(procnum /= 0)then
			call recv_data_from_master
		end if

		call mpi_scatterv(tempcrowdist_global,counts_masterproc,disp_masterproc,&
		mpi_double_precision,tempcrowdist_local,mycounts,mpi_double_precision,0,&
		mpi_comm_world,ierr)

		call mpi_scatterv(temprank_global,counts_masterproc,disp_masterproc,&
		mpi_integer,temprank_local, mycounts, mpi_integer,0,mpi_comm_world,ierr)
		call mpi_barrier(mpi_comm_world,ierr)
		
		do i = 1, nproc 
			if(i - 1 == procnum)then
				!copy these localrank and localdist to the parent_pop%dist
				! and parent_pop%rank
				parent_pop(1:iend-istart + 1)%rank = temprank_local
				parent_pop(1:iend-istart + 1)%crowd_dist = tempcrowdist_local
			end if
		end do
	
		!rescale the index of the parent_pop 
		do i = 1, iend - istart + 1
			parent_pop(i)%index = istart + i - 1
		end do
		call mpi_barrier(mpi_comm_world,ierr)
		if(procnum == 0)then
			deallocate(front_intermediate)
			print *,'igen',igen
			!inquire(file="subset_pop_canonical.txt", exist=exist)
			 ! if (exist) then
			  !  open(30, file="subset_pop_canonical.txt", status="old", position="append", action="write")
			  !else
			  !  open(30, file="subset_pop_canonical.txt", status="new", action="write")
			  !end if
			!inquire(file="subset_obj_canonical.txt",exist=exist)
			!if(exist)then
			!	open(53,file='subset_obj_canonical.txt',status='old',position='append',action='write')
			!else
			!	open(53,file='subset_obj_canonical.txt',status='new',action ='write')
			!end if
			call calc_acceptable_models
			goto 789
			if(igen == 1)then
				open(unit = 33, file = 'init_front.txt',status = 'unknown', action = 'write')
				do i = 1, popsize
					write(33,*)parent_pop_global(i)%obj
				end do
				close(33)
			end if
			if(mod(igen, 100) == 0)then
				write(filename,'("front",I3,".txt")')igen
				open(unit = igen,file = filename,status='unknown', action = 'write')
				do i = 1, popsize
					write(igen,*)parent_pop_global(i)%obj
				end do
				close(igen)
			end if
	789		if(mod(igen,10) == 0)then
				open(unit = 30, file = 'subset_pop_twenty_t1.txt',status ='unknown',action = 'write')
				open(unit = 53, file = 'subset_obj_twenty_t1.txt',status ='unknown',action = 'write')
				do i = 1, acceptable_member_count
					write(30,100)acceptable_models(i)%model_parameter
					write(53,*)acceptable_models(i)%obj
				end do
				close(30)
				close(53)
			end if
		
			if(mod(igen, 10) == 0)then
				open(unit = 4, file = 'temp_pop_twenty_t1.txt',status = 'unknown', action = 'write')
				open(unit = 5, file = 'temp_obj_twenty_t1.txt',status = 'unknown',action ='write')
				do i = 1, popsize
					!if(parent_pop_global(i)%rank == 1)then
					    write(5,*)parent_pop_global(i)%obj
					    write(4,100)parent_pop_global(i)%model_parameter
						!end if
				end do
				close(4)
				close(5)
			end if
		end if
	end do              ! do igen = 1, maxgen
	
	!write the final pop and final objective to files
	!if(igen == maxgen + 1)then
	!	if(procnum == 0)then
	!		open(unit = 4, file = 'final_pop_canonical.txt',status = 'unknown', action = 'write')
	!		open(unit = 5, file = 'final_obj_canonical.txt',status = 'unknown',action ='write')
	!		do i = 1, popsize
	!			if(parent_pop_global(i)%rank == 1)then
	!				! storing only those models which have misfit -1 <0
	!				if(parent_pop_global(i)%obj(1) .le. 1.d0)then
	!				    write(5,*)parent_pop_global(i)%obj
	!				    write(4,100)parent_pop_global(i)%model_parameter
	!				end if
	!			end if
	!		end do
	!		close(4)
	!		close(5)
	!	end if
	!end if
	
	100 format(53(1x,f18.14))
	
	if(procnum == 0)then
		end_time = mpi_Wtime()
		print *,'execution time(s)',end_time - start_time
	end if

	call mpi_finalize(ierr)
end program main