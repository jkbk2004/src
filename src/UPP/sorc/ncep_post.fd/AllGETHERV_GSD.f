      SUBROUTINE AllGETHERV(GRID1)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    AllGETHERV VERT INTRP OF MODEL LVLS TO PRESSURE
!   PRGRMMR: MING HU           ORG: GSD     DATE: 2012-01-01
!
! ABSTRACT:
!   .
!
! PROGRAM HISTORY LOG:
!
!    21-09-02  Bo Cui - Decompose UPP in X direction               
      
     use ctlblk_mod, only : im,jm,num_procs,me,jsta,jend,ista,iend,mpi_comm_comp

     implicit none

     include "mpif.h"

!
    integer i,j,ij
    integer ierr

     REAL GRID1(IM,JM)
     REAL ibufrecv(IM*JM)
     REAL ibufsend((iend-ista+1)*(jend-jsta+1))
     integer SENDCOUNT,RECVCOUNTS(num_procs),DISPLS(num_procs)
!
!     write(*,*) 'check mpi', im,jm,num_procs,me,jsta,jend
     SENDCOUNT=(iend-ista+1)*(jend-jsta+1)
     call MPI_ALLGATHER(SENDCOUNT, 1, MPI_INTEGER, RECVCOUNTS,1 , &
                MPI_INTEGER, mpi_comm_comp, ierr)
     DISPLS(1)=0
     do i=2,num_procs
         DISPLS(i)=DISPLS(i-1)+RECVCOUNTS(i-1)
     enddo
!
!     write(*,*) me,'RECVCOUNTS=',RECVCOUNTS
!     write(*,*) me,'DISPLS=',DISPLS
! 
     ij=0
     ibufsend=0.0
     do j=jsta,jend
        do i=ista,iend
           ij=ij+1
           ibufsend(ij)=GRID1(i,j)
        enddo
     enddo
     if(ij /= RECVCOUNTS(me+1)) then
        write(*,*) 'Error: send account is not equal to receive account',me,ij,RECVCOUNTS(me+1)
     endif
  
     call MPI_ALLGATHERV(ibufsend, ij, MPI_REAL, ibufrecv, RECVCOUNTS,DISPLS, &
                MPI_REAL, mpi_comm_comp, ierr)

     ij=0
     do j=1,JM
        do i=1,IM
           ij=ij+1
           GRID1(i,j)=ibufrecv(ij)
        enddo
     enddo
!
!     END OF ROUTINE.
!
      RETURN
      END
