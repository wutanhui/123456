SUBROUTINE Acc0(y,dy,ddy,r,dr,omega,Le,ele_num_all,RouA,Iz,Gra,a0,P,Mp,Cm,Rou,Fn,F_num)
!    use body_module
!    use Element_module
!    use Model_module
!   use integration_perimeter
    use OPERATORS_MOD   
    use BLAS95
    use LAPACK95
    implicit none
    real(8),allocatable::r(:),dr(:),omega(:),P(:),Mp(:);
    real(8),allocatable::Z(:,:),z1(:),y(:),dy(:),ddy(:),Z_U(:,:),Z1_U(:);
    real(8)::Le,RouA,Iz(3,3),Gra(3),a0(3),Cm(3,3),Rou;
    real(8)::A1(3,3)
    integer::ele_num_all,F_num
    integer,allocatable::FN(:),IPIV(:)
    INTERFACE
        SUBROUTINE recur(Z,z1,r,dr,omega,y,Le,ele_num_all,dy,A1,RouA,Iz,Gra,a0,P,Mp,Cm,Rou)  
            real(8),allocatable:: r(:),dr(:),omega(:),P(:),Mp(:),y(:),dy(:),Z(:,:),Z1(:)
            real(8)::Le,RouA,Iz(3,3),Gra(3),a0(3),A1(3,3),Cm(3,3),Rou
            integer::ele_num_all
        END SUBROUTINE   
    END INTERFACE
    !F_num=3*ele_num_all;
    allocate(Z(3*ele_num_all+6,3*ele_num_all+6),Z1(3*ele_num_all+6))
    allocate(Z_U(F_num,F_num),Z1_u(F_num),ipiv(F_num))



    !FN=[7:3*ele_num_all+6];!!¸Ä

    A1=reshape([  cos(y(5))*cos(y(6)), sin(y(4))*sin(y(5))*cos(y(6))+cos(y(4))*sin(y(6)),-cos(y(4))*sin(y(5))*cos(y(6))+sin(y(4))*sin(y(6)),&
                &-cos(y(5))*sin(y(6)),-sin(y(4))*sin(y(5))*sin(y(6))+cos(y(4))*cos(y(6)), cos(y(4))*sin(y(5))*sin(y(6))+sin(y(4))*cos(y(6)),&
                & sin(y(5)),          -sin(y(4))*cos(y(5)),                               cos(y(4))*cos(y(5))],[3,3]);
    !A1=reshape([0d0,0d0,-1d0,0d0,1d0,0d0,1d0,0d0,0d0],[3,3])
    !A1=reshape([0d0,1d0,0d0,-1d0,0d0,0d0,0d0,0d0,1d0],[3,3])

    !call recur(Z,z1,r,dr,omega,y,Le,ele_num_all,dy,A1,RouA,Iz,Gra,a0,P,Mp,Cm,Rou)
    
    !dx(1:ele_num_all+6)=x(ele_num_all+7:2*ele_num_all+12);
    !
    !
    Z_U=Z(FN,FN)
    !call exportf(Z,'Z.txt')
    z1_U=z1(FN)
    !call potrf(Z_U)
    !call potrs(Z_U,z1_U)
    call getrf(Z_U,ipiv)
    call getrs(Z_U,ipiv,Z1_U)
    ddy(Fn)=z1_U;
    !call exportf(Z,'Z.txt')
    !call exportf(Z1,'Z1.txt')
    call exportf(Z1_U,'Z1_u.txt')
    
    
endsubroutine Acc0