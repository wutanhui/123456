
    
    
    
    !subroutine check(Time,Time_I,y,dy,ddy,r,dr,omega,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ)
!    use body_module
!    use functions
!    implicit none
!    real(8),allocatable::x(:),dx(:),r(:),dr(:),omega(:),con_value(:);
!    real(8),allocatable::y(:),dy(:),ddy(:),Z_U(:,:),Z1_u(:),Phi_C(:,:),gama_w(:),lumda_temp(:),Phi_C2(:,:),gama_w2(:),lumda_temp2(:),Equ(:),lumda_Quater_temp(:);
!    real(8)::Time,Gra(3),a0(3);
!    !real(8)::A1(3,3)
!    integer::ele_num_all,F_num_all,con_num_all,con_num_act,body_num,II,JJ,KK,LL,MM,OO,HH,NN,I_begin,I_end,Time_I;
!    integer,allocatable::Con_num2
!    integer,allocatable::FN(:),ipiv(:),Con_DOF(:,:),Con_list(:,:),CON_dof2(:),Act_Cn(:)
!    !type(Body_type),allocatable:: Body(:)
!    INTERFACE
!        SUBROUTINE recur_g(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,Gra,a0,Mp,Contact_node_num,Contact_node,Body_G,GG,Phi_quater,Array,Array_P)  
!                  !recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Gra,a0,P,Mp)    
!            real(8),allocatable:: r(:),dr(:),omega(:),Mp(:),y(:),dy(:),Z(:,:),Z1(:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),gama_w2(:),Phi_C3(:,:),Body_G(:,:),GG(:,:),Phi_quater(:),Array(:,:),Array_P(:,:,:)
!            real(8)::Time,Le,RouA,Iz(3,3),Gra(3),a0(3),A1(3,3),Cm(3,3),Rou
!            integer::ele_num,Con_num,CON_dof(:,:),Con_num2,CON_dof2(:),Contact_node_num,Contact_node(:)
!            !integer,allocatable::
!        END SUBROUTINE   
!    END INTERFACE
!    
!    
!    I_begin=1
!    do II=1,body_num
!        I_end=I_begin+3*body(II).ele_num+6
!        body(II).y_temp=y(I_begin:I_end)
!        body(II).dy_temp=dy(I_begin:I_end)
!        body(II).ddy_temp=ddy(I_begin:I_end)
!        I_begin=I_end+1;
!    enddo
!
!    
!    !con_num_all=sum(body(:).con_num)
!    allocate(Z_U(F_num_all+body_num+con_num_act,F_num_all+body_num+con_num_act),Z1_u(F_num_all+body_num+con_num_act),ipiv(F_num_all+body_num+con_num_act))
!    Z_U=0d0;
!    z1_U=0d0;
!
!    do II=1,Body_num
!        body(II).A1=dir_cos_q(Body(II).y_temp(4:7))
!
!        call recur_g(Time,body(II).Z,body(II).z1,body(II).r_temp,body(II).dr_temp,body(II).omega_temp,body(II).y_temp,body(II).dy_temp,body(II).Le,body(II).A1,body(II).ele_num,body(II).RouA,body(II).Iz,body(II).Cm,body(II).Rou,body(II).Con_num,body(II).CON_dof,body(II).Phi_C,body(II).gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,body(II).Phi_C3,Gra,a0,body(II).Mp,body(II).contact_node_num,body(II).Contact_node,Body(II).G,Body(II).GG,Body(II).Phi_quater,body(II).Array,body(II).Array_P)
!            
!    enddo
!   ! call  contact_detect(body_num)
!    
!
!    
!    I_begin=1
!    do II=1,Body_num   
!        I_end=I_begin+body(II).f_num-1
!        !Z_U(I_begin:I_end,I_begin:I_end)=body(II).Z(body(II).Fn,body(II).Fn)
!        Z1_U(I_begin:I_end)=body(II).Z1(body(II).Fn)
!        I_begin=I_end+1
!    enddo
! 
!    JJ=0;    
!    do OO=1,Con_num_all
!        if (Con_list(OO,10)==1) then
!            JJ=JJ+1
!        endif
!        if (Con_list(OO,1)==0.and.Con_list(OO,10)==1) then
!            II=Con_list(OO,2)
!            LL=Con_list(OO,8)
!            !if (body(II).con_num>0) then
!            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
!            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
!            z1_U(f_num_all+body_num+JJ)=-con_value(OO)+body(II).gama_w(LL)
!            
!            z1_d(f_num_all+body_num+JJ)=body(II).gama_d(LL)
!            z1_v(f_num_all+body_num+JJ)=body(II).gama_v(LL)
!        endif
!        
!        if (Con_list(OO,1)==1.and.Con_list(OO,10)==1) then
!            II=Con_list(OO,2)
!            KK=Con_list(OO,5)
!            LL=Con_list(OO,8)
!            MM=Con_list(OO,9)
!                    
!            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
!            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-body(KK).Phi_C(MM,body(KK).Fn)
!            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
!            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)=-body(KK).Phi_C(MM,body(KK).Fn)
!            !
!            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)-body(KK).gama_w(MM)
!
!        endif
!        if (Con_list(OO,1)==2.and.Con_list(OO,10)==1) then
!            II=Con_list(OO,2)
!            LL=Con_list(OO,8)               
!            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
!            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
!            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)
!            z1_d(f_num_all+body_num+JJ)=body(II).gama_d(LL)
!            z1_v(f_num_all+body_num+JJ)=body(II).gama_v(LL)
!        endif
!        
!        if (Con_list(OO,1)==102.and.Con_list(OO,10)==1) then
!            II=Con_list(OO,2)
!            LL=Con_list(OO,8)
!            !if (body(II).con_num>0) then
!            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
!            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
!            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)
!            z1_d(f_num_all+body_num+JJ)=body(II).gama_d(LL)
!            z1_v(f_num_all+body_num+JJ)=body(II).gama_v(LL)
!        endif
!        
!                
!        !if (Con_list(OO,1)==103.and.Con_list(OO,10)==1) then
!        !    II=Con_list(JJ,2)
!        !    HH=Con_list(JJ,3)
!        !    LL=Con_list(JJ,8)
!        !    
!        !    KK=Con_list(JJ,5)
!        !    NN=Con_list(JJ,6)
!        !    MM=Con_list(JJ,9)            
!        !    
!        !    
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -matmul(transpose(body(II).Array_P(:,body(II).Fn,LL)),body(KK).Array(:,MM))
!        !    
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
!        !    
!        !    Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -matmul(transpose(body(II).Array_P(:,body(ii).Fn,LL)),body(KK).Array(:,MM))
!        !    
!        !    Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
!        !    
!        !    !z1_U(f_num_all+body_num+JJ)=2d0*dot_product(cross(body(KK).omega_temp(3*NN-2:3*NN),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),body(KK).Array(:,MM))-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),body(KK).Array(:,MM))-dot_product(cross(body(KK).omega_temp(3*NN-2:3*NN),cross(body(kk).omega_temp(3*NN-2:3*NN),body(kk).Array(:,MM))),body(II).Array(:,LL))+dot_product(body(II).Array_gama(:,LL),body(KK).Array(:,MM))+dot_product(body(II).Array(:,LL),body(KK).Array_gama(:,MM))
!        !    
!        !endif
!        !if (Con_list(OO,1)==104.and.Con_list(OO,10)==1) then
!        !    II=Con_list(OO,2)
!        !    LL=Con_list(OO,8)
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
!        !                                                                                                                                                                  !
!        !    Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
!        !    
!        !    !z1_U(f_num_all+body_num+JJ)=-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),[body(II).r_temp(1),body(II).r_temp(2),0d0])+dot_product(body(II).Array_gama(:,LL),[body(II).r_temp(1),body(II).r_temp(2),0d0])-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL)),[body(II).dr_temp(1),body(II).dr_temp(2),0d0])![1d0,0d0,0d0])!
!        !    
!        !endif
!        !
!        !if (Con_list(OO,1)==105.and.Con_list(OO,10)==1) then
!        !    II=Con_list(OO,2)
!        !    LL=Con_list(OO,8)
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
!        !                                                                                                                                                                  !
!        !    Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
!        !    
!        !    z1_U(f_num_all+body_num+JJ)=dot_product(body(II).Array(:,LL),[0d0,0d0,1d0])![1d0,0d0,0d0])!
!        !    
!        !endif        
!        !if (Con_list(OO,1)==3.and.Con_list(OO,10)==1) then
!        !    II=Con_list(OO,2)
!        !    KK=Con_list(OO,5)
!        !    LL=Con_list(OO,8)
!        !    MM=Con_list(OO,9)
!        !    !if (body(II).con_num>0) then
!        !            
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= ((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(II).Phi_C(LL,body(II).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(II).Phi_C3(LL,body(II).Fn))
!        !    Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(KK).Phi_C(MM,body(KK).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(KK).Phi_C3(MM,body(KK).Fn))
!        !    !Z_U(f_num_all+sum(body(1:-1).con_num)+1:F_num_all+sum(body(1:II-1).con_num)+body(II).con_num,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=body(II).Phi_C(:,body(II).Fn)
!        !    Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))
!        !    Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))
!        !    
!        !    z1_U(f_num_all+body_num+JJ)=-(body(II).dr_temp(3*Con_list(OO,3)-2)-body(kk).dr_temp(3*Con_list(OO,6)-2))**2-(body(II).dr_temp(3*Con_list(OO,3)-1)-body(kk).dr_temp(3*Con_list(OO,6)-1))**2
!        !    
!        !    !ddy(f_num_all+JJ)=body(II).lumda(LL)body(II).gama_w(LL)+body(KK).gama_w(MM)!
!        !endif
!        
!        
!    enddo 
!    
!    
!
!        
!    
!    
!endsubroutine
!
!subroutine recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,Gra,a0,Mp,Contact_node_num,Contact_node,Body_G,GG,Phi_quater,Array,Array_P)
!!    use body_module
!!    use Element_module
!!    use Model_module
!    use functions
!    use OPERATORS_MOD
!    implicit none
!    
!
!    real(8)::Time,Le,RouA,Iz(3,3),Gra(3),a0(3),Cm(3,3),Rou;
!    real(8)::A1(3,3)
!    integer::ele_num,Con_num,Con_num2,contact_node_num;
!    integer,allocatable::CON_dof(:,:),Con_dof2(:),contact_node(:);
!    
!    integer:: ii,kk,jj,ii_contact,LL
!    real(8):: Eye(3,3),A2(3,3)
!    real(8):: U(9,3),Roup_T(3,3),omega_T(3,3),Me(9,9),Fe(9),Fec(6)!,g_temp(9)
!    real(8):: phi(3,3),Pi(3),beta1(3),beta2(3),Ki(3,3),omega_r(3),Kqq2(3,3),ggg(6);
!    
!    real(8)::R_Q(3,4)
!    real(8),allocatable::Phi_quater(:)
!    
!    real(8),allocatable::r(:),dr(:),omega(:),y(:),dy(:),Mp(:);
!    real(8),allocatable::Z(:,:),z1(:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),Phi_C3(:,:),gama_w2(:),Body_G(:,:);
!    
!    real(8),allocatable:: G(:,:),g0(:,:),ones(:),M(:,:),F(:)
!    real(8),allocatable:: GG(:,:),gg0(:,:)!为末端点约束而建
!    real(8),allocatable:: a(:),da(:),F1(:),T(:,:),Roup(:,:),dRoup(:,:),beta(:,:)
!    
!    real(8):: A2_P1(3,3),A2_P2(3,3),A2_P3(3,3)
!    real(8):: A0_P1(3,3),A0_P2(3,3),A0_P3(3,3),A0_P0(3,3)
!
!    real(8),allocatable:: A1_P1(:,:),A1_P2(:,:),A1_P3(:,:),A1_P0(:,:)
!    real(8),allocatable:: Array(:,:),Array_P(:,:,:),AI(:,:)
!    
!    real(8)::temp_T(3,3)
!    
!    allocate( G(9*ele_num,3*ele_num+7),g0(9*ele_num,ele_num),M(9*ele_num,9*ele_num),F(9*ele_num),F1(3*ele_num+7),beta(9,ele_num+1))
!    allocate( gg0(6,ele_num))
!    allocate( T(9*ele_num+9,1:9),a(3*ele_num),da(3*ele_num),Roup(3,ele_num),droup(3,ele_num),ones(ele_num))
!    
!    allocate( A1_P1(3*ele_num+3,3*ele_num+3),A1_P2(3*ele_num+3,3*ele_num+3),A1_P3(3*ele_num+3,3*ele_num+3),A1_P0(3*ele_num+3,3))
!    allocate( AI(3,3*ele_num+3))
!    
!    !allocate(con_dof(con_num),Phi_C(con_num,3*ele_num_all+6),gama_w(con_num))
!    ones=1.d0;
!    temp_T=0d0;
!    !!
! !   write(*,*) ones
!    U=0d0; U(7,1)=1d0;U(8,2)=1d0;U(9,3)=1d0;
!    G=0d0;
!    GG=0d0;
!    G(1,1)=1d0;G(2,2)=1d0;G(3,3)=1d0;
!    !Ki=reshape([1d0,0d0,0d0,0d0,cos(y(4)),sin(y(4)),sin(y(5)),-cos(y(5))*sin(y(4)),cos(y(5))*cos(y(4))],[3,3])
!    R_Q=2d0*R_Quater(y(4:7))
!    !Ki=reshape([cos(y(5))*cos(y(6)),-cos(y(5))*sin(y(6)),sin(y(5)),sin(y(6)),cos(y(6)),0d0,0d0,0d0,1d0],[3,3])
!    G(4:6,4:7)=R_Q;
!    !G(4,4)=1d0;G(5,5)=1d0;G(6,6)=1d0;
!    G(7,8)=1d0;G(8,9)=1d0;G(9,10)=1d0;
!    g0=0d0;
! !   write(*,*) ones
!    r(1:3)=y(1:3);
!    dr(1:3)=dy(1:3);
!    omega(1:3)=R_Q.mt.dy(4:7);
!         
!    a = y(8:3*ele_num+7);
!    da=dy(8:3*ele_num+7);
!    Phi=0d0;
!    Z=0.d0;
!    Z1=0.d0;
!    F1=0d0;
!    T=0d0;
!    Eye=0d0;Eye(1,1)=1d0;Eye(2,2)=1d0;Eye(3,3)=1d0;
!    M=0d0;
!    F=0d0;
!    beta=0d0
!   ! p=0d0;
!    fec=0d0;
!    ii_contact=1;
!    Phi_c=0d0
!    !--------------根节点方向余弦，对欧拉四元数的导数-----------------!
!    
!    AI=0d0;
!    
!    Ai(1:3,1:3)=A1;
!    Array=0d0;Array_P=0d0
!    !--------------根节点方向余弦，对欧拉四元数的导数-----------------!
!    
!    
!    do ii=1,ele_num
!        call Phi_Pi(droup(:,ii),Roup(:,ii),PHI,Pi,a(3*ii-2:3*ii),da(3*ii-2:3*ii),Le)
!        Roup(:,ii) =A1.mt.Roup(:,ii)
!        dRoup(:,ii)=A1.mt.dRoup(:,ii)
!        PHI=A1.mt.Phi
!        Pi =A1.mt.Pi
!        r(3*ii+1:3*ii+3)=r(3*ii-2:3*ii)+Roup(:,ii);
!        !!!!!!!!!!!!!!!!!!!!!!!!!!约束
!       
!        !!!!!!!!!!!!!!!!!!!!!!!!!!
!        Ki=K_R(a(3*ii-2:3*ii))
!        Kqq2(1,3)= da(3*ii-1)*cos(a(3*ii-1))
!        Kqq2(2,2)=-da(3*ii-2)*sin(a(3*ii-2))
!        Kqq2(2,3)=-da(3*ii-2)*cos(a(3*ii-2))*cos(a(3*ii-1))+da(3*ii-1)*sin(a(3*ii-2))*sin(a(3*ii-1))
!        Kqq2(3,2)= da(3*ii-2)*cos(a(3*ii-2))
!        Kqq2(3,3)=-da(3*ii-2)*sin(a(3*ii-2))*cos(a(3*ii-1))-da(3*ii-1)*cos(a(3*ii-2))*sin(a(3*ii-1))  
!        !Ki=reshape([cos(a(3*ii-1))*cos(a(3*ii)),-cos(a(3*ii-1))*sin(a(3*ii)),sin(a(3*ii-1)),sin(a(3*ii)),cos(a(3*ii)),0d0,0d0,0d0,1d0],[3,3])
!        omega_r=A1.mt.Ki.mt.da(3*ii-2:3*ii)
!        call Tensor(Roup_T,Roup(:,ii));
!        T(9*ii+1:9*ii+3,1:9)=reshape([eye,-Roup_T,phi],[3,9]);
!        T(9*ii+4:9*ii+6,4:6)=Eye;
!        T(9*ii+4:9*ii+6,7:9)=A1.mt.Ki;
!        !Eye=A1*Ki
!        !Eye=matmul(A1,Ki)
!        dr(3*ii+1:3*ii+3)=dRoup(:,ii)-(Roup_T.mt.omega(3*ii-2:3*ii))+dr(3*ii-2:3*ii);!T(9*i)       
!        
!        omega(3*ii+1:3*ii+3)=omega(3*ii-2:3*ii)+omega_r;!
!        
!        call Tensor(omega_T,omega(3*ii-2:3*ii));
!        beta1=(omega_T.mt.omega_T.mt.Roup(:,ii))+2d0*(omega_T.mt.droup(:,ii))+PI
!        
!        beta2=(A1.mt.Kqq2.mt.da(3*ii-2:3*ii))+(omega_T.mt.omega_r)!!!改
!        
!        A2=Dir_COS(a(3*ii-2:3*ii))
!        beta(:,ii+1)=[beta1,beta2,0d0,0d0,0d0];
!        A1=A1.mt.A2!
!        Ai(1:3,3*ii+1:3*ii+3)=A1;        
!        do kk=1,ii+1
!            if (kk<ii+1.and.ii<ele_num) then
!                if (kk==1) then
!                    G(9*ii+1:9*ii+9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
!                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
!                else
!                    G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
!                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
!                endif
!            elseif (kk==ii+1.and.ii<ele_num) then
!                G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=U;
!                g0(9*ii+1:9*ii+9,kk)=beta(:,ii+1);
!            endif;
!            
!            if (ii==ele_num) then
!                if (kk==1) then
!                    GG(1:9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
!                    gG0(1:6,kk)=T(9*ii+1:9*ii+6,1:9).mt.g0(9*ii-8:9*ii,kk);
!                elseif(1<kk.and.kk<ii+1) then
!                    GG(1:9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
!                    gg0(1:6,kk)=T(9*ii+1:9*ii+6,1:9).mt.g0(9*ii-8:9*ii,kk);
!                else
!                    ggg=beta(1:6,ii+1);
!                endif
!            endif
!
!        enddo
!    
!    enddo
!
!    Phi_Quater=0d0;
!    Phi_Quater(4)=2d0*y(4)
!    Phi_Quater(5)=2d0*y(5)
!    Phi_Quater(6)=2d0*y(6)
!    Phi_Quater(7)=2d0*y(7)
!
!    Body_G=G;
!        
!    do ii=1,con_num
!        jj=Con_dof(ii,1);
!       
!        if (Con_dof(ii,2)==0.or.Con_dof(ii,2)==1) then
!             kk=con_dof(ii,3);
!                if (jj<=9*ele_num) then 
!                    Phi_C(ii,:)=G(jj,:);
!                    gama_d(ii)=-r(3*kk-3+con_dof(ii,4))
!                    gama_v(ii)=-dr(3*kk-3+con_dof(ii,4))
!                else
!                    Phi_C(ii,:)=GG(jj-9*ele_num,:);                    
!                    gama_d(ii)=-r(3*kk-3+con_dof(ii,4))
!                    gama_v(ii)=-dr(3*kk-3+con_dof(ii,4))
!                endif
!        endif
!        if (Con_dof(ii,2)==2) then
!            kk=con_dof(ii,3);
!            if (jj<=9*ele_num) then 
!                Phi_C(ii,:)=2d0*r(3*kk-2)*G(jj,:)+2d0*r(3*kk-1)*G(jj+1,:);
!                gama_d(ii)=r(3*kk-2)**2+r(3*kk-1)**2-0.3d0**2
!                gama_v(ii)=2d0*dr(3*kk-2)**2+2d0*dr(3*kk-1)**2
!            else
!                Phi_C(ii,:)=2d0*r(3*kk-2)*GG(jj-9*ele_num,:)+2d0*r(3*kk-1)*GG(jj+1-9*ele_num,:); 
!                gama_d(ii)=r(3*kk-2)**2+r(3*kk-1)**2-0.3d0**2
!                gama_v(ii)=2d0*dr(3*kk-2)**2+2d0*dr(3*kk-1)**2
!            endif
!        endif
!
!        ! if (Con_dof(ii,2)==103.or.Con_dof(ii,2)==104.or.Con_dof(ii,2)==105) then
!        !    LL=Con_dof(ii,3)
!        !    if(Con_dof(ii,4)==1) then
!        !        Array(:,ii)=AI(:,3*LL-2)
!        !    elseif(Con_dof(ii,4)==2) then
!        !        Array(:,ii)=AI(:,3*LL-1)
!        !    elseif(Con_dof(ii,4)==3) then
!        !        Array(:,ii)=AI(:,3*LL)
!        !    endif
!        !     
!        !    call tensor(temp_T,Array(:,ii))           
!        !    if (jj<=9*ele_num) then 
!        !        Array_P(:,:,ii)=temp_T.mt.G(jj:jj+2,:)
!        !    else
!        !        Array_P(:,:,ii)=temp_T.mt.GG(4:6,:)
!        !    endif
!        ! endif
!        !
!        !if (Con_dof(ii,2)==3) then
!        !    kk=con_dof(ii,3);
!        !    if (jj<=9*ele_num) then 
!        !        Phi_C(ii,:)=2d0*G(jj,:)
!        !        Phi_C3(ii,:)=2d0*G(jj+1,:);
!        !        !gama_w(ii)=r(3*kk-1)**2+r(3*kk-2)**2
!        !    else
!        !        Phi_C(ii,:)=2d0*GG(jj-9*ele_num,:)
!        !        Phi_C3(ii,:)=2d0*GG(jj+1-9*ele_num,:); 
!        !        !gama_w(ii)=r(3*kk-1)**2+r(3*kk-2)**2
!        !    endif
!        !endif
!    enddo
!endsubroutine recur_g