MODULE OPERATORS_MOD
!######################
!##################
implicit none
    INTERFACE OPERATOR(.w.)
        module procedure wave
        module procedure wave1
    ENDINTERFACE
    INTERFACE operator(.mt.)
        module procedure mtm_m
        module procedure mtm_m3
        module procedure mtm_v
        module procedure mtv_m
        module procedure mtv_v
    ENDINTERFACE
    interface exportf
        module procedure exportf_m
        module procedure exportf_mm
        module procedure exportf_v
        module procedure exportf_vv
    endinterface
    interface importf
        module procedure importf_m
        module procedure importf_v
        module procedure importf_mI
        module procedure importf_vI
    endinterface
    interface zeros
        module procedure zeros_m
        module procedure zeros_v
    endinterface
contains
!#################################################################################    
    FUNCTION  mtm_m(mat1,mat2)
        real(8),intent(in):: mat1(:,:)
        real(8),intent(in):: mat2(:,:)
        real(8):: mtm_m(size(mat1,1),size(mat2,2)) 
            mtm_m = matmul(mat1,mat2)
    ENDFUNCTION mtm_m
    
    FUNCTION  mtm_m3(mat1,mat2)
        real(8),intent(in):: mat1(:,:)
        real(8),intent(in):: mat2(:,:,:)
        real(8):: mtm_m3(size(mat1,1),size(mat2,2),size(mat2,3)) 
        integer:: ii,jj
        mtm_m3=0d0
        do ii=1,size(mat1,1)
            do jj=1,size(mat1,2)
               mtm_m3(ii,:,:) =mtm_m3(ii,:,:)+mat1(ii,jj)*mat2(jj,:,:)
            enddo
        enddo         
    ENDFUNCTION mtm_m3
    
    FUNCTION mtm_v(mat,vec)
        real(8),intent(in):: mat(:,:),vec(:)
        real(8):: mtm_v(size(mat,1))      
            mtm_v = matmul(mat,vec)
    ENDFUNCTION mtm_v
    
    FUNCTION mtv_m(vec,mat)
        real(8),intent(in):: vec(:),mat(:,:)
        real(8):: mtv_m(size(mat,2))       
            mtv_m = matmul(vec,mat)
    ENDFUNCTION mtv_m
    
    FUNCTION mtv_v(vec1,vec2)
        real(8),intent(in):: vec1(:),vec2(:)
        real(8):: mtv_v(size(vec1,1),size(vec2,1))  
        real(8):: v1(size(vec1,1),1),v2(1,size(vec2,1))
            v1(:,1)=vec1;
            V2(1,:)=vec2;
            mtv_v = matmul(v1,v2)
    ENDFUNCTION mtv_v
!################################################################################

    function wave(x)
    implicit none
        real(8),intent(in):: x(3)
        real(8):: wave(3,3)
        wave(1,1) = 0.d0
        wave(1,2) = -x(3)
        wave(1,3) = x(2)
        wave(2,1) = x(3)
        wave(2,2) = 0.d0
        wave(2,3) = -x(1)
        wave(3,1) = -x(2)
        wave(3,2) = x(1)
        wave(3,3) = 0.d0
    endfunction
    
    function wave1(x)
    implicit none
        real(8),intent(in):: x(3,1)
        real(8):: wave1(3,3)
        wave1(1,1) = 0.d0
        wave1(1,2) = -x(3,1)
        wave1(1,3) = x(2,1)
        wave1(2,1) = x(3,1)
        wave1(2,2) = 0.d0
        wave1(2,3) = -x(1,1)
        wave1(3,1) = -x(2,1)
        wave1(3,2) = x(1,1)
        wave1(3,3) = 0.d0
    endfunction
 !###########################################################################   
!    SUBROUTINE dimrank(mat,lett)
!    implicit none
!        real(8):: mat(:,:)
!        character(*),optional:: lett
!        character(31):: lett1
!        if (.not.present(lett)) then
!            lett1 = 'matrix_01'
!        else
!            lett1 = lett
!        endif
!        print*,'Dimension size of ',trim(lett1), ' : (',size(mat,1),',',size(mat,2),')'
!        print*,'Rank of ',trim(lett1),' : ',rank(mat)
!    ENDSUBROUTINE dimrank

    SUBROUTINE exportf_m(mat,file)
    implicit none
        real(8):: mat(:,:)
        character(*),optional:: file
        integer:: ma,mb,i,j
        if (present(file)) then
            open(11,file=file)
        else
            open(11,file='001.txt')
        endif

        ma = size(mat,1)
        mb = size(mat,2)
        do i=1,ma
            do j=1,mb
                write(11,'(E20.10,\)')mat(i,j)
            enddo
            write(11,*)
        enddo   
        close(11)
    ENDSUBROUTINE exportf_m

    SUBROUTINE exportf_mm(mat,file,a)
    implicit none
        real(8):: mat(:,:)
        character(*),optional:: file
        integer:: ma,mb,i,j,a
        if (present(file)) then
            open(11,file=file,ACCESS='append')
        else
            open(11,file='001.txt')
        endif

        ma = size(mat,1)
        mb = size(mat,2)
        do i=1,ma
            do j=1,mb
                write(11,'(E20.10,\)')mat(i,j)
            enddo
            write(11,*)
        enddo   
        close(11)
    ENDSUBROUTINE exportf_mm
    SUBROUTINE exportf_v(mat,file)
    implicit none
        real(8):: mat(:)
        character(*),optional:: file
        integer:: ma,i
        if (present(file)) then
            open(11,file=file)
        else
            open(11,file='001.txt')
        endif

        ma = size(mat)
        do i=1,ma
                write(11,'(E40.30)')mat(i)
        enddo
        close(11);
    ENDSUBROUTINE exportf_v
SUBROUTINE exportf_vv(mat,file,a)
    implicit none
        real(8):: mat(:)
        character(*),optional:: file
        integer:: ma,i,a
        if (present(file)) then
            open(11,file=file,ACCESS='append')
        else
            open(11,file='001.txt')
        endif

        ma = size(mat)
        do i=1,ma
                write(11,'(E40.30)')mat(i)
        enddo
        close(11);
ENDSUBROUTINE exportf_vv
    SUBROUTINE importf_mI(mat,file)
    implicit none
        integer:: mat(:,:)
        character(*):: file
        integer:: ma,mb,i,j
        ma = size(mat,1)
        mb = size(mat,2)
        open(11,file=file,STATUS='OLD')
        do i = 1,ma
            read(11,*) (mat(i,j),j=1,mb)
        enddo
        close(11)
    ENDSUBROUTINE importf_mI
    SUBROUTINE importf_vI(mat,file)
        implicit none
        integer:: mat(:)
        character(*):: file
        integer:: ma,i
        ma = size(mat,1)
        open(11,file=file,STATUS='OLD')
        do i = 1,ma
            read(11,*) mat(i)
        enddo
        close(11)
    ENDSUBROUTINE importf_vI
    SUBROUTINE importf_m(mat,file)
    implicit none
        real(8):: mat(:,:)
        character(*):: file
        integer:: ma,mb,i,j
        ma = size(mat,1)
        mb = size(mat,2)
        open(11,file=file,STATUS='OLD')
        do i = 1,ma
            read(11,*) (mat(i,j),j=1,mb)
        enddo
        close(11)
    ENDSUBROUTINE importf_m
    SUBROUTINE importf_v(mat,file)
        implicit none
        real(8):: mat(:)
        character(*):: file
        integer:: ma,i
        ma = size(mat,1)
        open(11,file=file,STATUS='OLD')
        do i = 1,ma
            read(11,*) mat(i)
        enddo
        close(11)
    ENDSUBROUTINE importf_v
!##########################################################################    
function zeros_m(m,n)
implicit none
    integer:: m,n
    real(8):: zeros_m(m,n)
    zeros_m = 0.d0
endfunction zeros_m

function zeros_v(m)
implicit none
    integer:: m
    real(8):: zeros_v(m)
    zeros_v = 0.d0
endfunction zeros_v

         
ENDMODULE operators_mod
!!#################################################################
!!#################################################################



