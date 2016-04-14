module CPU_timer
USE IFPORT
implicit none
!#########################################
type,public:: secnds_timer
    private
    real(4):: t1
contains
    procedure,public:: tic
    procedure,public:: toc
endtype
!########################################
!type,public:: omp_timer
!    private
!    real(4):: t1
!contains
!    procedure,public:: tic
!    procedure,public:: toc
!endtype
!#########################################
private:: tic,toc
contains
subroutine tic(this)
implicit none
    class(secnds_timer):: this
    this.t1 = secnds(0.0)
endsubroutine
!**************************************
subroutine toc(this,label)
implicit none
    class(secnds_timer):: this
    character(*):: label
    this.t1 = secnds(this.t1)
    write(*,'(A,A,F18.6,A)')label, ' costs ', this.t1, ' s'
endsubroutine
!################################################################
!subroutine tic(this)
!implicit none
!    class(omp_timer):: this
!    this.t1 = omp_get_wtime()
!endsubroutine
!***********************************
!subroutine toc(this,label)
!implicit none
!    class(omp_timer):: this
!    character(*):: label
!    this.t1 = omp_get_wtime(this.t1)
!    write(*,'(A,A,F10.6,A)')label, ' costs ', this.t1, ' s'
!endsubroutine



















!!#####################################
subroutine Timer_Start()
implicit none
    real(8):: t1
    t1 = timef()
endsubroutine
subroutine Timer_End(label)
implicit none
    real(8):: t1
    character(*):: label
    write(*,'(A,A,F10.6,A)') label, ' costs ', timef(), ' s'
endsubroutine
!!#######################################


endmodule