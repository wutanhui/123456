SUBROUTINE Equation(Time,Time_I,y,dy,ddy,r,dr,omega,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ)
    !ele_num_body,body_num,RouA_body,Iz_body,Cm_body,Rou_body,Fn_body,F_num_body,Con_num_body,CON_dof_Body,
!    use body_module
!    use Element_module
!    use Model_module
    use OPERATORS_MOD   
    use BLAS95
    use LAPACK95
    use body_module
    use functions
    !use lin_sol_gen_int
    implicit none
    real(8),allocatable::x(:),dx(:),r(:),dr(:),omega(:),con_value(:);
    real(8),allocatable::y(:),dy(:),ddy(:),Z_U(:,:),Z1_u(:),Phi_C(:,:),gama_w(:),lumda_temp(:),Phi_C2(:,:),gama_w2(:),lumda_temp2(:),Equ(:),lumda_Quater_temp(:);
    real(8)::Time,Gra(3),a0(3);
    !real(8)::A1(3,3)
    integer::ele_num_all,F_num_all,con_num_all,con_num_act,body_num,II,JJ,KK,LL,MM,OO,HH,NN,I_begin,I_end,Time_I;
    integer,allocatable::Con_num2
    integer,allocatable::FN(:),ipiv(:),Con_DOF(:,:),Con_list(:,:),CON_dof2(:),Act_Cn(:)
    !type(Body_type),allocatable:: Body(:)
    INTERFACE
        SUBROUTINE recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,Gra,a0,Mp,Contact_node_num,Contact_node,Body_G,GG,Phi_quater,Array,Array_P)  
                  !recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Gra,a0,P,Mp)    
            real(8),allocatable:: r(:),dr(:),omega(:),Mp(:),y(:),dy(:),Z(:,:),Z1(:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),gama_w2(:),Phi_C3(:,:),Body_G(:,:),GG(:,:),Phi_quater(:),Array(:,:),Array_P(:,:,:)
            real(8)::Time,Le,RouA,Iz(3,3),Gra(3),a0(3),A1(3,3),Cm(3,3),Rou
            integer::ele_num,Con_num,CON_dof(:,:),Con_num2,CON_dof2(:),Contact_node_num,Contact_node(:)
            !integer,allocatable::
        END SUBROUTINE   
    END INTERFACE
    INTERFACE
        SUBROUTINE Abs_Coor(y,r,Roup,AI,Ele_num,Le)
            real(8),allocatable:: r(:),y(:),dy(:),Roup(:),AI(:,:)
            real(8)::Le
            integer::ele_num
        END SUBROUTINE   
    END INTERFACE
    
    INTERFACE
        SUBROUTINE Abs_Vel(y,dy,dr,dRoup,Phi,AI,Ele_num,Le)
            real(8),allocatable:: y(:),dr(:),dy(:),dRoup(:),Phi(:,:),AI(:,:)
            real(8)::Le
            integer::ele_num
        END SUBROUTINE   
    END INTERFACE
    INTERFACE
        SUBROUTINE Abs_omega(y,dy,omega,omega_r,Ki,AI,Ele_num)
            real(8),allocatable:: y(:),dy(:),omega(:),omega_r(:),AI(:,:),Ki(:,:)
            integer::ele_num
        END SUBROUTINE   
    END INTERFACE
    INTERFACE
        SUBROUTINE recur_matrix(G,g0,y,dy,r,dr,AI,Phi,Roup,droup,Ki,omega,omega_r,Ele_num,Le)
            real(8),allocatable::G(:,:),g0(:,:),y(:),dy(:),r(:),dr(:),AI(:,:),Phi(:,:),Roup(:),dRoup(:),Ki(:,:),omega(:),omega_r(:)
            real(8)::Le
            integer::ele_num
        END SUBROUTINE   
    END INTERFACE
    
    !!由状态量X到各物体的分量！
    I_begin=1
    do II=1,body_num
        I_end=I_begin+3*body(II).ele_num+6
        body(II).y_temp=y(I_begin:I_end)
        body(II).dy_temp=dy(I_begin:I_end)
        body(II).ddy_temp=ddy(I_begin:I_end)
        I_begin=I_end+1;
    enddo

    
    !con_num_all=sum(body(:).con_num)
    allocate(Z_U(F_num_all+body_num+con_num_act,F_num_all+body_num+con_num_act),Z1_u(F_num_all+body_num+con_num_act),ipiv(F_num_all+body_num+con_num_act))
    Z_U=0d0;
    z1_U=0d0;
    
     do II=1,Body_num    
         call Abs_Coor(body(II).y_temp,body(II).r_temp,body(II).Roup_temp,body(II).AI_temp,body(II).Ele_num,body(II).Le)
         call Abs_Vel(body(II).y_temp,body(II).dy_temp,body(II).dr_temp,body(II).dRoup_temp,body(II).Phi_temp,body(II).AI_temp,body(II).Ele_num,body(II).Le)
         call Abs_Omega(body(II).y_temp,body(II).dy_temp,body(II).omega_temp,body(II).omega_r_temp,body(II).Ki_temp,body(II).AI_temp,body(II).Ele_num)
         call recur_matrix(body(II).G,body(II).g0,body(II).y_temp,body(II).dy_temp,body(II).r_temp,body(II).dr_temp,body(II).AI_temp,body(II).Phi_temp,body(II).Roup_temp,body(II).droup_temp,body(II).Ki_temp,body(II).omega_temp,body(II).omega_r_temp,body(II).Ele_num,body(II).Le)
     enddo
    
    
    

    do II=1,Body_num
        body(II).A1=dir_cos_q(Body(II).y_temp(4:7))

        call recur(Time,body(II).Z,body(II).z1,body(II).r_temp,body(II).dr_temp,body(II).omega_temp,body(II).y_temp,body(II).dy_temp,body(II).Le,body(II).A1,body(II).ele_num,body(II).RouA,body(II).Iz,body(II).Cm,body(II).Rou,body(II).Con_num,body(II).CON_dof,body(II).Phi_C,body(II).gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,body(II).Phi_C3,Gra,a0,body(II).Mp,body(II).contact_node_num,body(II).Contact_node,Body(II).G,Body(II).GG,Body(II).Phi_quater,body(II).Array,body(II).Array_P)
            !recur(Time,Z,z1,r,dr,omega,         y,              dy,              Le,         A1,         ele_num,         RouA,         Iz,         Cm,         Rou,         Con_num,         CON_dof,         Phi_C,         gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Gra,a0,P,Mp)
    enddo
    !!!%%%%%%%%%%%%%% 加入接触检测段对相应的z1时行修改%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   ! call  contact_detect(body_num)
    

    
    I_begin=1
    do II=1,Body_num   
        I_end=I_begin+body(II).f_num-1
        Z_U(I_begin:I_end,I_begin:I_end)=body(II).Z(body(II).Fn,body(II).Fn)
        Z1_U(I_begin:I_end)=body(II).Z1(body(II).Fn)
        I_begin=I_end+1
    enddo
    
    JJ=0;    
    do OO=1,Con_num_all
        if (Con_list(OO,10)==1) then
            JJ=JJ+1
        endif
        if (Con_list(OO,1)==0.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            !if (body(II).con_num>0) then
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            z1_U(f_num_all+body_num+JJ)=-con_value(OO)+body(II).gama_w(LL)
            
          
        endif
        
        if (Con_list(OO,1)==1.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            !if (body(II).con_num>0) then
                    
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-body(KK).Phi_C(MM,body(KK).Fn)
            !Z_U(f_num_all+sum(body(1:-1).con_num)+1:F_num_all+sum(body(1:II-1).con_num)+body(II).con_num,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=body(II).Phi_C(:,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)=-body(KK).Phi_C(MM,body(KK).Fn)
            !
            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)-body(KK).gama_w(MM)
            
            !ddy(f_num_all+JJ)=body(II).lumda(LL)
            !endif
        endif
        if (Con_list(OO,1)==2.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)               
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)
        endif
        
        if (Con_list(OO,1)==102.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            !if (body(II).con_num>0) then
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)
        endif
        
                
        if (Con_list(OO,1)==103.and.Con_list(OO,10)==1) then
            II=Con_list(JJ,2)
            HH=Con_list(JJ,3)
            LL=Con_list(JJ,8)
            
            KK=Con_list(JJ,5)
            NN=Con_list(JJ,6)
            MM=Con_list(JJ,9)            
            
            
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -matmul(transpose(body(II).Array_P(:,body(II).Fn,LL)),body(KK).Array(:,MM))
            
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
            
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -matmul(transpose(body(II).Array_P(:,body(ii).Fn,LL)),body(KK).Array(:,MM))
            
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
            
            !z1_U(f_num_all+body_num+JJ)=2d0*dot_product(cross(body(KK).omega_temp(3*NN-2:3*NN),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),body(KK).Array(:,MM))-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),body(KK).Array(:,MM))-dot_product(cross(body(KK).omega_temp(3*NN-2:3*NN),cross(body(kk).omega_temp(3*NN-2:3*NN),body(kk).Array(:,MM))),body(II).Array(:,LL))+dot_product(body(II).Array_gama(:,LL),body(KK).Array(:,MM))+dot_product(body(II).Array(:,LL),body(KK).Array_gama(:,MM))
            
        endif
        if (Con_list(OO,1)==104.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
                                                                                                                                                                          !
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
            
            !z1_U(f_num_all+body_num+JJ)=-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL))),[body(II).r_temp(1),body(II).r_temp(2),0d0])+dot_product(body(II).Array_gama(:,LL),[body(II).r_temp(1),body(II).r_temp(2),0d0])-dot_product(cross(body(II).omega_temp(3*HH-2:3*HH),body(II).Array(:,LL)),[body(II).dr_temp(1),body(II).dr_temp(2),0d0])![1d0,0d0,0d0])!
            
        endif

        if (Con_list(OO,1)==105.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
                                                                                                                                                                          !
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
            
            z1_U(f_num_all+body_num+JJ)=dot_product(body(II).Array(:,LL),[0d0,0d0,1d0])![1d0,0d0,0d0])!
            
        endif        
        if (Con_list(OO,1)==3.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            !if (body(II).con_num>0) then
                    
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= ((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(II).Phi_C(LL,body(II).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(II).Phi_C3(LL,body(II).Fn))
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(KK).Phi_C(MM,body(KK).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(KK).Phi_C3(MM,body(KK).Fn))
            !Z_U(f_num_all+sum(body(1:-1).con_num)+1:F_num_all+sum(body(1:II-1).con_num)+body(II).con_num,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=body(II).Phi_C(:,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))
            
            z1_U(f_num_all+body_num+JJ)=-(body(II).dr_temp(3*Con_list(OO,3)-2)-body(kk).dr_temp(3*Con_list(OO,6)-2))**2-(body(II).dr_temp(3*Con_list(OO,3)-1)-body(kk).dr_temp(3*Con_list(OO,6)-1))**2
            
            !ddy(f_num_all+JJ)=body(II).lumda(LL)body(II).gama_w(LL)+body(KK).gama_w(MM)!
        endif
        
        
    enddo    
    
    
    do ii=1,body_num
       Z_U(ii+F_num_all,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num)) = body(II).Phi_Quater(body(II).Fn)
       Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),ii+F_num_all) = body(II).Phi_Quater(body(II).Fn)
       Z1_U(ii+F_num_all)=-2*dot_product(body(ii).dy_temp(4:7),body(ii).dy_temp(4:7));
    enddo
        
    Equ=Z1_U;
    call gesv(Z_U,Equ) 
    
    ddy(Fn)=Equ(1:f_num_all)
    lumda_quater_temp=Equ(F_num_all+1:F_num_all+body_num)
    lumda_temp(act_cn)=Equ(F_num_all+body_num+1:F_num_all+body_num+con_num_act)
    !dx(3*ele_num_all+6*body_num+Fn)=z1_U(1:F_num_all);
    !if (con_num>0) then
    !lumda_temp=z1_U(f_num_all+1:F_num_all+con_num_all)
    !    if (con_num2>0) then
    !        lumda_temp2(con_dof2)=z1_U(f_num+con_num2+1:F_num+con_num+con_num2)
    !    endif
    !endif
    !if (con_num==0.and.con_num2>0) then
    !    lumda_temp2(con_dof2)=z1_U(f_num+1:F_num+con_num2)
    !endif
    !deallocate(gw)
    !call exportf(Equ,'Equ.txt') 
    !call exportf(Z1_u,'Z1_u3.txt') 
    !call potrf(Z_U)
    !call potrs(Z_U,z1_U)
    !call getrf(Z_U,ipiv)
    !call getrs(Z_U,ipiv,Z1_U)
    !dx(3*ele_num_all+13:6*ele_num_all+12)=z1_U;
    !call exportf(Z_U,'Z03.txt')
    !call exportf(Z_U,'Z13.txt')
    !call exportf(Z1_U,'Z1.txt')
    !call exportf(A1,'A13.txt')
    !call exportf(r,'r.txt')
endsubroutine equation

subroutine recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,Gra,a0,Mp,Contact_node_num,Contact_node,Body_G,GG,Phi_quater,Array,Array_P)
!    use body_module
!    use Element_module
!    use Model_module
    use functions
    use OPERATORS_MOD
    implicit none
    INTERFACE
        SUBROUTINE kinetic(y,dy,r,dr,Rou,dRou,Omega,Omega_r,AI,Phi,Pi,Beta,Ele_num,Le)
            real(8)::y(:),dy(:),r(:),dr(:),Rou(:),dRou(:),Omega(:),Omega_r(:),AI(:,:),PHi(:,:),Pi(:),Beta(:,:),Le
            integer::Ele_num
        END SUBROUTINE   
    END INTERFACE
    INTERFACE
        SUBROUTINE abs_coor(y,r,Rou,AI,Ele_num,Le)
            real(8)::y(:),r(:),Rou(:),AI(:,:),Le
            integer::Ele_num
        END SUBROUTINE   
    END INTERFACE

    real(8)::Time,Le,RouA,Iz(3,3),Gra(3),a0(3),Cm(3,3),Rou;
    real(8)::A1(3,3)
    integer::ele_num,Con_num,Con_num2,contact_node_num;
    integer,allocatable::CON_dof(:,:),Con_dof2(:),contact_node(:);
    
    integer:: ii,kk,jj,ii_contact,LL
    real(8):: Eye(3,3),A2(3,3)
    real(8):: U(9,3),Roup_T(3,3),omega_T(3,3),Me(9,9),Fe(9),Fec(6)!,g_temp(9)
    real(8):: phi(3,3),Pi(3),beta1(3),beta2(3),Ki(3,3),omega_r(3),Kqq2(3,3),ggg(6);
    
    real(8)::R_Q(3,4)
    real(8),allocatable::Phi_quater(:)
    
    real(8),allocatable::r(:),dr(:),omega(:),y(:),dy(:),Mp(:);
    real(8),allocatable::Z(:,:),z1(:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),Phi_C3(:,:),gama_w2(:),Body_G(:,:);
    
    real(8),allocatable:: G(:,:),g0(:,:),ones(:),M(:,:),F(:)
    real(8),allocatable:: GG(:,:),gg0(:,:)!为末端点约束而建
    real(8),allocatable:: a(:),da(:),F1(:),T(:,:),Roup(:,:),dRoup(:,:),beta(:,:)
    
    real(8):: A2_P1(3,3),A2_P2(3,3),A2_P3(3,3)
    real(8):: A0_P1(3,3),A0_P2(3,3),A0_P3(3,3),A0_P0(3,3)

    real(8),allocatable:: A1_P1(:,:),A1_P2(:,:),A1_P3(:,:),A1_P0(:,:)
    real(8),allocatable:: Array(:,:),Array_P(:,:,:),AI(:,:)
    
    real(8),allocatable:: Phi0(:,:),Pi0(:),Roup0(:),dRoup0(:),Omega_r0(:)
    real(8)::temp_T(3,3)
    
    allocate( G(9*ele_num,3*ele_num+7),g0(9*ele_num,ele_num),M(9*ele_num,9*ele_num),F(9*ele_num),F1(3*ele_num+7),beta(9,ele_num+1))
    allocate( gg0(6,ele_num))
    allocate( T(9*ele_num+9,1:9),a(3*ele_num),da(3*ele_num),Roup(3,ele_num),droup(3,ele_num),ones(ele_num))
    
    allocate( A1_P1(3*ele_num+3,3*ele_num+3),A1_P2(3*ele_num+3,3*ele_num+3),A1_P3(3*ele_num+3,3*ele_num+3),A1_P0(3*ele_num+3,3))
    allocate( AI(3,3*ele_num+3))
    
    
    allocate(Phi0(3,3*ele_num),Pi0(3*ele_num+3),Roup0(3*ele_num),dRoup0(3*ele_num),Omega_r0(3*ele_num))
    
    
    
    !allocate(con_dof(con_num),Phi_C(con_num,3*ele_num_all+6),gama_w(con_num))
    ones=1.d0;
    temp_T=0d0;
    !!
 !   write(*,*) ones
    U=0d0; U(7,1)=1d0;U(8,2)=1d0;U(9,3)=1d0;
    G=0d0;
    GG=0d0;
    G(1,1)=1d0;G(2,2)=1d0;G(3,3)=1d0;
    !Ki=reshape([1d0,0d0,0d0,0d0,cos(y(4)),sin(y(4)),sin(y(5)),-cos(y(5))*sin(y(4)),cos(y(5))*cos(y(4))],[3,3])
    R_Q=2d0*R_Quater(y(4:7))
    !Ki=reshape([cos(y(5))*cos(y(6)),-cos(y(5))*sin(y(6)),sin(y(5)),sin(y(6)),cos(y(6)),0d0,0d0,0d0,1d0],[3,3])
    G(4:6,4:7)=R_Q;
    !G(4,4)=1d0;G(5,5)=1d0;G(6,6)=1d0;
    G(7,8)=1d0;G(8,9)=1d0;G(9,10)=1d0;
    g0=0d0;
 !   write(*,*) ones
    r(1:3)=y(1:3);
    dr(1:3)=dy(1:3);
    omega(1:3)=R_Q.mt.dy(4:7);
         
    a = y(8:3*ele_num+7);
    da=dy(8:3*ele_num+7);
    Phi=0d0;
    Z=0.d0;
    Z1=0.d0;
    F1=0d0;
    T=0d0;
    Eye=0d0;Eye(1,1)=1d0;Eye(2,2)=1d0;Eye(3,3)=1d0;
    M=0d0;
    F=0d0;
    beta=0d0
   ! p=0d0;
    fec=0d0;
    ii_contact=1;
    Phi_c=0d0
    !--------------根节点方向余弦，对欧拉四元数的导数-----------------!
    A0_P1=Dir_cos_Q1(y(4:7))
    A0_P2=Dir_cos_Q2(y(4:7))
    A0_P3=Dir_cos_Q3(y(4:7))
    A0_P0=Dir_cos_Q0(y(4:7))
    A1_P1=0d0;
    A1_P2=0d0;
    A1_P3=0d0;
    A1_P0=0d0;
    A1_P1(1:3,1:3)=A0_P1;
    A1_P2(1:3,1:3)=A0_P2;
    A1_P3(1:3,1:3)=A0_P3;
    A1_P0(1:3,1:3)=A0_P0;
    
    AI=0d0;
    
    Ai(1:3,1:3)=A1;
    Array=0d0;Array_P=0d0
    !--------------根节点方向余弦，对欧拉四元数的导数-----------------!
    
    call abs_coor(y,r,Roup0,AI,Ele_num,Le)
    call Kinetic(y,dy,r,dr,Roup0,dRoup0,Omega,Omega_r0,AI,Phi0,Pi0,Beta,Ele_num,Le)
    
    
    
    do ii=1,ele_num
        call Phi_Pi(droup(:,ii),Roup(:,ii),PHI,Pi,a(3*ii-2:3*ii),da(3*ii-2:3*ii),Le)
        Roup(:,ii) =A1.mt.Roup(:,ii)
        dRoup(:,ii)=A1.mt.dRoup(:,ii)
        PHI=A1.mt.Phi
        Pi =A1.mt.Pi
        r(3*ii+1:3*ii+3)=r(3*ii-2:3*ii)+Roup(:,ii);
        !!!!!!!!!!!!!!!!!!!!!!!!!!约束
        Fec=0d0;
        if (contact_node_num>0.and.ii<ele_num) then
            if (contact_node_num==1) then
                if (r(3*ii+3)>1d0.and.r(3*ii+3)<1.21d0) then!.and.mod(ii,10)==0
                    Fec(1:3)=-[1d5*(r(3*ii+1)-0.3d0),1d5*(r(3*ii+2)),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];
                endif
            endif
            if (contact_node_num==2) then
                if (r(3*ii+3)>1.0d0.and.r(3*ii+3)<1.21d0) then!.and.mod(ii,10)==0
                    Fec(1:3)=-[1d4*(r(3*ii+1)+0.15d0),1d2*(r(3*ii+2)-0.259807621135332d0),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];;
                endif
            endif
            if (contact_node_num==3) then
                if (r(3*ii+3)>1.0d0.and.r(3*ii+3)<1.21d0) then!.and.mod(ii,10)==0
                    Fec(1:3)=-[1d4*(r(3*ii+1)+0.15d0),1d2*(r(3*ii+2)+0.259807621135332d0),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];;
                endif
            endif
            !if (sqrt(r(3*ii+1)**2+r(3*ii+2)**2)<0.2990d0) then
            !    Fec(1:3)=Fec(1:3)+0.5d-4*[r(3*ii+1),r(3*ii+2),0d0];
            !endif
             if (r(3*ii+3)<ii*0.001d0) then
                Fec(1:3)=Fec(1:3)-1d4*[0d0,0d0,(r(3*ii+3)-ii*0.001d0)];
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        Ki=K_R(a(3*ii-2:3*ii))
        Kqq2(1,3)= da(3*ii-1)*cos(a(3*ii-1))
        Kqq2(2,2)=-da(3*ii-2)*sin(a(3*ii-2))
        Kqq2(2,3)=-da(3*ii-2)*cos(a(3*ii-2))*cos(a(3*ii-1))+da(3*ii-1)*sin(a(3*ii-2))*sin(a(3*ii-1))
        Kqq2(3,2)= da(3*ii-2)*cos(a(3*ii-2))
        Kqq2(3,3)=-da(3*ii-2)*sin(a(3*ii-2))*cos(a(3*ii-1))-da(3*ii-1)*cos(a(3*ii-2))*sin(a(3*ii-1))  
        !Ki=reshape([cos(a(3*ii-1))*cos(a(3*ii)),-cos(a(3*ii-1))*sin(a(3*ii)),sin(a(3*ii-1)),sin(a(3*ii)),cos(a(3*ii)),0d0,0d0,0d0,1d0],[3,3])
        omega_r=A1.mt.Ki.mt.da(3*ii-2:3*ii)
        call Tensor(Roup_T,Roup(:,ii));
        T(9*ii+1:9*ii+3,1:9)=reshape([eye,-Roup_T,phi],[3,9]);
        T(9*ii+4:9*ii+6,4:6)=Eye;
        T(9*ii+4:9*ii+6,7:9)=A1.mt.Ki;
        !Eye=A1*Ki
        !Eye=matmul(A1,Ki)
        dr(3*ii+1:3*ii+3)=dRoup(:,ii)-(Roup_T.mt.omega(3*ii-2:3*ii))+dr(3*ii-2:3*ii);!T(9*i)       
        
        omega(3*ii+1:3*ii+3)=omega(3*ii-2:3*ii)+omega_r;!
        
        call Tensor(omega_T,omega(3*ii-2:3*ii));
        beta1=(omega_T.mt.omega_T.mt.Roup(:,ii))+2d0*(omega_T.mt.droup(:,ii))+PI
        
        beta2=(A1.mt.Kqq2.mt.da(3*ii-2:3*ii))+(omega_T.mt.omega_r)!!!改
        
        A2=Dir_COS(a(3*ii-2:3*ii))
        !!--------------单元相对转动，产生的方向余弦对本单元相对转角的偏导数------20150405------------------------!
        A2_P1=Dir_cos_P1(a(3*ii-2:3*ii))
        A2_P2=Dir_cos_P2(a(3*ii-2:3*ii))
        A2_P3=Dir_cos_P3(a(3*ii-2:3*ii))
        !--------------单元相对转动，产生的方向余弦对本单元相对转角的偏导数-----20150405-------------------------!
        !---------------求方向余旋的偏导数-------------!
        do kk=0,ii               
            if(ii==kk) then
                A1_P1(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1,A2_P1)
                A1_P2(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1,A2_P2)
                A1_P3(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1,A2_P3)
            else
                A1_P1(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3),A2)
                A1_P2(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3),A2)
                A1_P3(3*ii+1:3*ii+3,3*kk+1:3*kk+3)=matmul(A1_P3(3*ii-2:3*ii,3*kk+1:3*kk+3),A2)
            endif
        enddo
        A1_P0(3*ii+1:3*ii+3,1:3)=matmul(A1_P0(3*ii-2:3*ii,1:3),A2)
        
        
        
        beta(:,ii+1)=[beta1,beta2,0d0,0d0,0d0];
        call Element_MF(Me,Fe,a(3*ii-2:3*ii),RouA,Le,A1,omega(3*ii-2:3*ii),da(3*ii-2:3*ii),Mp(6*ii-5:6*ii)+Fec,Iz,Gra,a0,Roup(:,ii),Phi,Cm,Ki,Rou)!
        M(9*ii-8:9*ii,9*ii-8:9*ii)=Me;
        F(9*ii-8:9*ii)=Fe;
        !A1=reshape([cos(x(4))*cos(x(6)),sin(x(4))*sin(x(5))*cos(x(6))+cos(x(5))*sin(x(6)),-cos(x(4))*sin(x(5))*cos(x(6))+sin(x(5))*sin(x(6)),-cos(x(5))*sin(x(6)),sin(x(4))*sin(x(5))*sin(x(6))+cos(x(5))*cos(x(6)),cos(x(4))*sin(x(5))*sin(x(6))+sin(x(5))*cos(x(6)),sin(x(5)),-sin(x(4))*cos(x(5)),cos(x(4))*cos(x(5))],[3,3]);
        A1=A1.mt.A2!
        Ai(1:3,3*ii+1:3*ii+3)=A1;        
        do kk=1,ii+1
            if (kk<ii+1.and.ii<ele_num) then
                if (kk==1) then
                    G(9*ii+1:9*ii+9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
                else
                    G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
                endif
            elseif (kk==ii+1.and.ii<ele_num) then
                G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=U;
                g0(9*ii+1:9*ii+9,kk)=beta(:,ii+1);
            endif;
            
            if (ii==ele_num) then
                if (kk==1) then
                    GG(1:9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
                    gG0(1:6,kk)=T(9*ii+1:9*ii+6,1:9).mt.g0(9*ii-8:9*ii,kk);
                elseif(1<kk.and.kk<ii+1) then
                    GG(1:9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
                    gg0(1:6,kk)=T(9*ii+1:9*ii+6,1:9).mt.g0(9*ii-8:9*ii,kk);
                else
                    ggg=beta(1:6,ii+1);
                endif
            endif
            
            !组集部分。
            !if(kk==1) then
            !    Z(1:9,1:9)=Z(1:9,1:9)+matmul(transpose(G(9*ii-8:9*ii,1:9)),matmul(Me,G(9*ii-8:9*ii,1:9)));
            !    
            !    !Z1(1:9)=Z1(1:9)+matmul(transpose(G(9*ii-8:9*ii,1:9)),matmul(Me,-g0(9*ii-8:9*ii,1))+Fe);
            !elseif(kk<=ele_num_all) then
            !    Z(1:9,3*kk+4:3*kk+6)=Z(1:9,3*kk+4:3*kk+6)+matmul(transpose(G(9*ii-8:9*ii,1:9)),matmul(Me,G(9*ii-8:9*ii,3*kk+4:3*kk+6)));
            !    Z(3*kk+4:3*kk+6,1:9)=transpose(Z(1:9,3*kk+4:3*kk+6))
            !
            !    !Z1(1:9)=Z1(1:9)+matmul(transpose(G(9*ii-8:9*ii,1:9)),matmul(Me,-g0(9*ii-8:9*ii,kk)))
            !    do jj=2,kk
            !            Z(3*jj+4:3*jj+6,3*kk+4:3*kk+6)=Z(3*jj+4:3*jj+6,3*kk+4:3*kk+6)+matmul(transpose(G(9*ii-8:9*ii,3*jj+4:3*jj+6)),matmul(Me,G(9*ii-8:9*ii,3*kk+4:3*kk+6)));
            !            !if (jj/=kk) then
            !            Z(3*kk+4:3*kk+6,3*jj+4:3*jj+6)=transpose(Z(3*jj+4:3*jj+6,3*kk+4:3*kk+6));
            !            !endif
            !            !
            !    enddo
            !    !g_temp=0d0
            !    !do jj=2,ii
            !    !    !Z1(3*kk+4:3*kk+6)=Z1(3*kk+4:3*kk+6)+dot_product((G(9*ii-8:9*ii,3*kk+4:3*kk+6)),matmul(Me,-g0(9*ii-8:9*ii,jj)));
            !    !    g_temp=g_temp+g0(9*ii-8:9*ii,jj)
            !    !enddo
            !    !!Z1(1:9)=Z1(1:9)+matmul(transpose(G(9*ii-8:9*ii,1:9)),matmul(Me,-g_temp))
            !    !Z1(3*kk+4:3*kk+6)=Z1(3*kk+4:3*kk+6)+matmul(transpose(G(9*ii-8:9*ii,3*kk+4:3*kk+6)),Fe-matmul(Me,g_temp));
            !endif
        enddo
    
    enddo
   
    !call exportf(Z,'Z31.txt')
    !print *,maxval(abs(Z-matmul(transpose(G),matmul(M,G))))
    Z=matmul(transpose(G),matmul(M,G));
    Z1=matmul(transpose(G),matmul(M,-matmul(g0,ones))+F);!+F1;
    Phi_Quater=0d0;
    Phi_Quater(4)=2d0*y(4)
    Phi_Quater(5)=2d0*y(5)
    Phi_Quater(6)=2d0*y(6)
    Phi_Quater(7)=2d0*y(7)
    !Z1=transpose(G).mt.F
    Body_G=G;
        
    do ii=1,con_num
        jj=Con_dof(ii,1);
       
        if (Con_dof(ii,2)==0.or.Con_dof(ii,2)==1) then
            !if(mod(jj,9)<4) then
             kk=con_dof(ii,3);
                if (jj<=9*ele_num) then 
                    Phi_C(ii,:)=G(jj,:);
                    gama_w(ii)=-dot_product(g0(jj,:),ones)
                    if(ii>2) then
                        gama_w(ii)=gama_w(ii)
                    endif
                else
                    Phi_C(ii,:)=GG(jj-9*ele_num,:);                    
                    gama_w(ii)=-dot_product(gg0(jj-9*ele_num,:),ones)-ggg(jj-9*ele_num)
                endif
        endif
        if (Con_dof(ii,2)==2) then
            kk=con_dof(ii,3);
            if (jj<=9*ele_num) then 
                Phi_C(ii,:)=2d0*r(3*kk-2)*G(jj,:)+2d0*r(3*kk-1)*G(jj+1,:);
                gama_w(ii)=-2d0*r(3*kk-2)*dot_product(g0(jj,:),ones)-2d0*r(3*kk-1)*dot_product(g0(jj+1,:),ones)-2d0*dr(3*kk-2)**2-2d0*dr(3*kk-1)**2;
            else
                Phi_C(ii,:)=2d0*r(3*kk-2)*GG(jj-9*ele_num,:)+2d0*r(3*kk-1)*GG(jj+1-9*ele_num,:); 
                gama_w(ii)=-2d0*r(3*kk-2)*(dot_product(gg0(jj-9*ele_num,:),ones)+ggg(jj-9*ele_num))-2d0*r(3*kk-1)*(dot_product(gg0(jj+1-9*ele_num,:),ones)+ggg(jj-9*ele_num+1))-2d0*dr(3*kk-2)**2-2d0*dr(3*kk-1)**2;
            endif
        endif

        ! if (Con_dof(ii,2)==103.or.Con_dof(ii,2)==104.or.Con_dof(ii,2)==105) then
        !    LL=Con_dof(ii,3)
        !    if(Con_dof(ii,4)==1) then
        !        Array(:,ii)=AI(:,3*LL-2)
        !    elseif(Con_dof(ii,4)==2) then
        !        Array(:,ii)=AI(:,3*LL-1)
        !    elseif(Con_dof(ii,4)==3) then
        !        Array(:,ii)=AI(:,3*LL)
        !    endif
        !     
        !    call tensor(temp_T,Array(:,ii))           
        !    if (jj<=9*ele_num) then 
        !        Array_P(:,:,ii)=temp_T.mt.G(jj:jj+2,:)
        !    else
        !        Array_P(:,:,ii)=temp_T.mt.GG(4:6,:)
        !    endif
        ! endif
        !
        !if (Con_dof(ii,2)==3) then
        !    kk=con_dof(ii,3);
        !    if (jj<=9*ele_num) then 
        !        Phi_C(ii,:)=2d0*G(jj,:)
        !        Phi_C3(ii,:)=2d0*G(jj+1,:);
        !        !gama_w(ii)=r(3*kk-1)**2+r(3*kk-2)**2
        !    else
        !        Phi_C(ii,:)=2d0*GG(jj-9*ele_num,:)
        !        Phi_C3(ii,:)=2d0*GG(jj+1-9*ele_num,:); 
        !        !gama_w(ii)=r(3*kk-1)**2+r(3*kk-2)**2
        !    endif
        !endif
    enddo
    !print *,'TEST'
    !gama_w(1)=gama_w(1)-32d0*cos(8d0*time)
    !gama_w(2)=gama_w(2)+32d0*sin(8d0*time)
    !!if (time<0.01d0) then
    !!    gama_w(3)=gama_w(3)+100
    !!endif
    !
    !!!!
    !do ii=1,con_num2
    !    jj=Con_dof2(ii);
    !    Phi_C2(ii,:)=G(jj,:);
    !    gama_w2(ii)=dot_product(g0(jj,:),ones);
    !enddo
    !do ii=1,con_num
    !    jj=Con_dof(ii);
    !    Phi_C(ii,:)=G(jj,:);
    !    gama_w(ii)=dot_product(g0(jj,:),ones);
    !enddo
    
    !call exportf(Z,'Z3.txt')
    !call exportf(M,'M3.txt')
    !call exportf(F,'F3.txt')
    !call exportf(Phi_C,'PHi_C.txt')
    !call exportf(G,'G3.txt')
    !call exportf(G0,'G03.txt')
    !call exportf(Z1,'Z13.txt')
    !call exportf(T,'T3.txt')
    !call exportf(A1,'A1.txt')
    !write(*,*) a
endsubroutine recur



subroutine Element_MF(Me,Fe,a,RouA,Le,A1,omega,da,Mp,Iz,Gra,a0,Roup,Phi,Cm,Ki,Rou)
    use OPERATORS_MOD
    use BLAS95
    use LAPACK95
    implicit none
    real(8):: Me(9,9),Me1(9,9),Me2(9,9),Fe(9),A1(3,3),omega(3),Mp(6),Gra(3),eye(3,3),Roup(3),Phi(3,3),Cm(3,3),Ki(3,3),Rou
    real(8):: a(3),Le,RouA,da(3),Iz(3,3),a0(3)
    real(8):: Fw1(9),Fw2(9),Fu(3),Fp(9)
    integer:: ipiv(9)
    call mass_force_E(Me1,Me2,Fw1,Fw2,Fu,a,Le,omega,A1,da,Gra,Cm,Iz,a0)
    
    call force_p(Fp,Mp,A1,Roup,Phi,Ki)
    
    !call exportf(RouA*Me,'Me30.txt')
    eye=reshape([1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0],[3,3]);
    Me(1:3,1:3)=RouA*Le*Eye;
    Me(1:3,4:6)=RouA*matmul(A1,matmul(Me1(1:3,4:6),transpose(A1)));
    Me(1:3,7:9)=RouA*matmul(A1,Me1(1:3,7:9));
    Me(4:6,1:3)=transpose(Me(1:3,4:6));
    Me(4:6,4:6)=matmul(matmul(A1,RouA*Me1(4:6,4:6)+Rou*Me2(4:6,4:6)),transpose(A1));
    Me(4:6,7:9)=matmul(A1,RouA*Me1(4:6,7:9)+Rou*Me2(4:6,7:9));
    Me(7:9,1:3)=transpose(Me(1:3,7:9));
    Me(7:9,4:6)=transpose(Me(4:6,7:9));
    Me(7:9,7:9)=RouA*Me1(7:9,7:9)+Rou*Me2(7:9,7:9);
    !call force_f(Fw,a,Le,omega,A1,da,Gra)
    
    
    !Fu=[0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,1d0];!!!
    !A1=matmul(A1,matmul(Me1(1:3,4:6),transpose(A1)))
    !call exportf(A1,'A13.txt')
        
    
    Fe=Fp-[0d0,0d0,0d0,0d0,0d0,0d0,Fu]-RouA*Fw1-Rou*Fw2;
    !call exportf(Fe,'Fe31.txt') 
    !call exportf(Me,'Me3.txt')
    !!!call potrf(Me)
    !!!call potrs(Me,Fe)
    !call getrf(Me,ipiv)
    !call getrs(Me,ipiv,Fe)
    !call exportf(Me1,'Me13.txt')
    !call exportf(Me2,'Me23.txt')
    !call exportf(Fe,'DDq.txt') 

endsubroutine Element_MF

subroutine Mass_force_E(M1,M2,Fw1,Fw2,Fu,a,Le,omega,A1,da,Gra,Cm,Iz,a0)
!   force_f(Fw,a,Le,omega,A1,da,Gra)
    implicit none
    real(8)::M1(9,9),M2(9,9),Fw1(9),Fw2(9),Fu(3),a(3),Le,omega(3),A1(3,3),da(3),Gra(3),Cm(3,3),Iz(3,3),a0(3);
    real(8)::Fue1(3),Fwe1(9),Fwe2(9),Me1(9,9),Me2(9,9),g(6),w(6),aa;
    integer::kk
    ! W=[0.236926885056189;0.236926885056189;0.478628670499366;0.478628670499366;0.568888888888889];
    ! g=[0.906179845938664;-0.906179845938664;0.538469310105683;-0.538469310105683;0];
     !g=[0.861136311594053,0.339981043584856,-0.861136311594053,-0.339981043584856];
     !W=[0.347854845137454,0.652145154862546,0.347854845137454,0.652145154862546];Me1=0;
     g=[-0.932469514203152d0,-0.661209386466265d0,-0.238619186083197d0,0.932469514203152d0,0.661209386466265d0,0.238619186083197d0];
    W=[0.17132449237917d0, 0.360761573048139d0,0.467913934572691d0, 0.17132449237917d0,  0.360761573048139d0, 0.467913934572691d0];
    !g=[-0.9602898564975363d0,-0.7966664774136268d0,-0.525532409916329d0,-0.1834346424956498d0,0.1834346424956498d0,0.525532409916329d0,0.7966664774136268d0,0.9602898564975363d0]
    !W=[0.1012285362903768d0,0.2223810344533745d0,0.3137066458778874d0,0.362683783378362d0,0.362683783378362d0,0.3137066458778874d0,0.2223810344533745d0,0.1012285362903768d0]
    !g=[-0.9894009349916499d0,-0.9445750230732326d0,-0.8656312023878318d0,-0.755404408355003d0,-0.6178762444026438d0,-0.4580167776572274d0,-0.2816035507792589d0,-0.09501250983763744d0,0.09501250983763744d0,0.2816035507792589d0,0.4580167776572274d0,0.6178762444026438d0,0.755404408355003d0,0.8656312023878318d0,0.9445750230732326d0,0.9894009349916499d0]
    !W=[0.02715245941175406d0,0.06225352393864778d0,0.0951585116824929d0,0.1246289712555339d0,0.1495959888165768d0,0.1691565193950026d0,0.1826034150449236d0,0.1894506104550685d0,0.1894506104550685d0,0.1826034150449236d0,0.1691565193950026d0,0.1495959888165768d0,0.1246289712555339d0,0.0951585116824929d0,0.06225352393864778d0,0.02715245941175406d0]
    
    
    M1=0.d0;
    M2=0d0;
    Fw1=0d0;
    Fw2=0d0;
    Fu=0d0;
    do kk=1,6
            aa=W(kk)*Le/2d0;
            call Mass_force_d(Me1,Me2,Fwe1,Fwe2,Fue1,g(kk)*Le/2d0+Le/2d0,a,Le,omega,A1,da,Gra,Cm,Iz,a0)
            M1 =M1+Me1*aa;
            M2 =M2+Me2*aa;
            !call force_d(Fe1,g(kk)*Le/2d0+Le/2d0,a,Le,omega,A1,da,Gra)
            Fw1=Fw1+Fwe1*aa;
            Fw2=Fw2+Fwe2*aa;
            Fu =Fu +Fue1*aa;
    enddo

endsubroutine Mass_force_E

subroutine Mass_force_d(Me1,Me2,Fwe1,Fwe2,Fue1,kesi,a,Le,omega,A1,da,Gra,Cm,Iz,a0)
    use OPERATORS_MOD
    use functions
    implicit none    
    real(8)::Fue1(3),Fwe1(9),Fwe2(9),Me1(9,9),Me2(9,9),kesi,a(3),Le,omega(3),A1(3,3),A2(3,3),da(3),Gra(3),Cm(3,3),Iz(3,3),a0(3)
    real(8)::droup(3),Roup(3),Phi(3,3),PI(3),Roup_T(3,3),omega_T(3,3),W1(3),W2(3),W3(3)
    real(8)::q(3),dq(3),Kq(3,3),Kqq(3,3),Kqq1(3,3),Kqq2(3,3),Ap(3,3),kapa(3)

    call Phi_Pi_Kesi(droup,Roup,Phi,Pi,a,da,Kesi,Le)
    call Tensor(Roup_T,Roup)
    !
    Me1=0d0;
    Me1(1:3,4:6)=-Roup_T;
    Me1(1:3,7:9)= Phi;
    Me1(4:6,4:6)=-Roup_T.mt.Roup_T;
    Me1(4:6,7:9)= Roup_T.mt.Phi;
    Me1(7:9,7:9)= transpose(Phi).mt.Phi;
    
    call tensor(omega_T,omega);
    W1=(omega_T.mt.Omega_T.mt.A1.mt.Roup)+2d0*(omega_T.mt.A1.mt.droup)+(A1.mt.Pi)+Gra;!!!
    !call exportf(Roup,'Roup.txt') 
    call Tensor(Roup_T,(A1.mt.Roup));
    !call exportf(Roup_T,'Roup_T.txt') 
    W2=Roup_T.mt.W1;
    !W2=matmul(matmul(matmul(matmul(matmul(omega_T,A1),Roup_T),Roup_T),transpose(A1)),omega)+matmul(Roup_T,2d0*matmul(matmul(omega_T,matmul(A1,Phi)),da)+matmul(A1,Pi)+Gra);
    W3=transpose(Phi).mt.transpose(A1).mt.W1;
    Fwe1=[W1,W2,W3];
    
    !call exportf(Fwe1,'Fwe1.txt')     
    q=kesi/Le*a;
    A2=Dir_COS(q(1:3));
    dq=kesi/Le*da;
    !Kq=reshape([1d0,0d0,0d0,0d0,cos(q(1)),sin(q(1)),sin(q(2)),-cos(q(2))*sin(q(1)),cos(q(2))*cos(q(1))],[3,3])
    !Kq=K_R(q(1:3))
    !Kqq1=0d0;
    !Kqq1(1,2)= q(3)*cos(q(2))
    !Kqq1(2,1)=-q(2)*sin(q(1))-q(3)*cos(q(1))*cos(q(2))
    !Kqq1(2,2)= q(3)*sin(q(1))*sin(q(2))
    !Kqq1(3,1)= q(2)*cos(q(1))-q(3)*sin(q(1))*cos(q(2))
    !Kqq1(3,2)=-q(3)*cos(q(1))*sin(q(2))
    !尝试，应变计算为截面的局部坐标系下，而不是单元的坐标系下。
    Kq=K_B(q(1:3))    
    Kqq1=0d0;
    Kqq1(1,2)=-q(1)*sin(q(2))*cos(q(3))
    Kqq1(1,3)= q(2)*cos(q(3))-q(1)*cos(q(2))*sin(q(3))
    Kqq1(2,2)= q(1)*sin(q(2))*sin(q(3))
    Kqq1(2,3)=-q(2)*sin(q(3))-q(1)*cos(q(2))*cos(q(3))
    Kqq1(3,2)= q(1)*cos(q(2))
    Kqq=transpose(Kqq1+Kq)/Le
    Fue1=(Kqq.mt.Cm.mt.Kq.mt.(a-a0))/Le
    !kapa=Kq.mt.(a-a0)/Le
    !Fue1=Kqq(:,1)*Cm(1,1)*kapa(1)+Kqq(:,2)*Cm(2,2)*kapa(2)+Kqq(:,3)*Cm(3,3)*kapa(3)
    
    !!if  (Fue1(2)/=0) then
    !!    print*,fue1(2)
    !!endif
    !!call exportf(Kqq1,'K1.txt')
    !!call exportf(Kq,'K.txt')
    !!call exportf(Cm,'Cm.txt')
    !!call exportf(fue1,'f.txt')
    !!!!考虑截面绕中线转动后广义惯性力的变化！ 
    Kq=K_R(q(1:3))
    Kqq2=0d0;
    Kqq2(1,3)= dq(2)*cos(q(2))
    Kqq2(2,2)=-dq(1)*sin(q(1))
    Kqq2(2,3)=-dq(1)*cos(q(1))*cos(q(2))+dq(2)*sin(q(1))*sin(q(2))
    Kqq2(3,2)= dq(1)*cos(q(1))
    Kqq2(3,3)=-dq(1)*sin(q(1))*cos(q(2))-dq(2)*cos(q(1))*sin(q(2))
    Ap=A1.mt.A2
    Fwe2=0d0
    Fwe2(4:6)=kesi/Le*matmul(matmul(matmul(matmul(matmul(matmul(Ap,Iz),transpose(Ap)),omega_T),A1),Kq),da)+kesi/Le*matmul(matmul(matmul(matmul(Ap,Iz),transpose(A2)),Kqq2),da);
    Fwe2(7:9)=(kesi/Le)**2*matmul(matmul(matmul(matmul(matmul(matmul(matmul(transpose(Kq),A2),Iz),transpose(Ap)),omega_T),A1),Kq),da)+(kesi/Le)**2*matmul(matmul(matmul(transpose(Kq),A2),matmul(Iz,matmul(transpose(A2),Kqq2))),da);   
    !!考虑截面绕中线转动后质量阵的变化！
    Me2=0d0
    Me2(4:6,4:6)=matmul(matmul(A2,Iz),transpose(A2));
    Me2(4:6,7:9)=kesi/Le*matmul(matmul(A2,Iz),matmul(transpose(A2),Kq));!Iz!
    !Me2(4:6,4:6)=Iz;
    Me2(7:9,7:9)=(kesi/Le)**2d0*matmul(matmul(transpose(Kq),A2),matmul(Iz,matmul(transpose(A2),Kq)))
    !!Me2=0d0
    endsubroutine Mass_force_d



subroutine force_p(Fp,Mp,A1,Roup,Phi,Ki)
    use OPERATORS_MOD
    implicit none
    real(8)::Fp(9),Mp(6),A1(3,3),Ki(3,3);
    real(8)::Roup(3),Phi(3,3),T(6,9),Roup_T(3,3);
    
    T=0d0;
    T(1,1)=1d0;T(2,2)=1d0;T(3,3)=1d0;
    call tensor(Roup_T,Roup);
    T(1:3,4:6)=-Roup_T;
    T(1:3,7:9)=Phi;
    T(4:6,4:6)=reshape([1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0],[3,3]);
    T(4:6,7:9)=A1.mt.Ki;
    Fp=transpose(T).mt.Mp;
endsubroutine force_p

subroutine Tensor(rdx_T,rdx)
    real(8)::rdx_T(3,3),rdx(3)
    rdx_T=reshape([0.d0,rdx(3),-rdx(2),-rdx(3),0.d0,rdx(1),rdx(2),-rdx(1),0.d0],[3,3])
endsubroutine Tensor



subroutine Phi_Pi(droup,Roup,Phi,Pi,a,da,Le)
    use OPERATORS_MOD
    implicit none
        real(8):: droup(3),Roup(3),phi(3,3),Pi(3),a(3),da(3),Le
        real(8):: bi1,bi2,dbi1,dbi2
        bi1=a(1)+a(2)
        bi2=a(1)-a(2)
        dbi1=da(1)+da(2)
        dbi2=da(1)-da(2)
        
        Roup(1)= Le*( (a(2))/2d0-(a(2)**3)/24d0+(a(2)**5)/720d0-(a(2)**7)/720d0/56d0+a(2)**9/720d0/56d0/90d0)
        Roup(2)=-Le/2d0*( (bi1+bi2)/2d0-(bi1**3+bi2**3)/24d0+(bi1**5+bi2**5)/720d0-(bi1**7+bi2**7)/720d0/56d0+(bi1**9+bi2**9)/720d0/56d0/90d0)
        Roup(3)= Le+Le/2d0*(-(bi1**2+bi2**2)/6d0 +(bi1**4+bi2**4)/120d0 -(bi1**6+bi2**6)/120d0/42d0+(bi1**8+bi2**8)/120d0/42d0/72d0)
        
        Phi=0d0;
        Phi(1,2)= Le/2d0+Le*(-(a(2)**2)/8d0+(a(2)**4)/144d0-(a(2)**6)/720d0/8d0+a(2)**8/720d0/56d0/10d0)
        Phi(2,1)=-Le/2d0-Le/2d0*(-(bi1**2+bi2**2)/8d0+(bi1**4+bi2**4)/144d0-(bi1**6+bi2**6)/720d0/8d0+(bi1**8+bi2**8)/720d0/56d0/10d0)
        Phi(2,2)=       -Le/2d0*(-(bi1**2-bi2**2)/8d0+(bi1**4-bi2**4)/144d0-(bi1**6-bi2**6)/720d0/8d0+(bi1**8-bi2**8)/720d0/56d0/10d0)
        Phi(3,1)= Le/2d0*(-(bi1+bi2)/3d0+(bi1**3+bi2**3)/30d0-(bi1**5+bi2**5)/840d0+(bi1**7+bi2**7)/120d0/42d0/9d0)
        Phi(3,2)= Le/2d0*(-(bi1-bi2)/3d0+(bi1**3-bi2**3)/30d0-(bi1**5-bi2**5)/840d0+(bi1**7-bi2**7)/120d0/42d0/9d0)
               
        droup=Phi.mt.da; 
        
        pi=0d0;
        pi(1)= Le*(-a(2)/4d0+(a(2)**3)/36d0-(a(2)**5)/120d0/8d0+a(2)**7/720d0/70d0)*da(2)**2
        Pi(2)=-Le/2d0*(-(bi1*dbi1**2+bi2*dbi2**2)/4d0+(bi1**3*dbi1**2+bi2**3*dbi2**2)/36d0-(bi1**5*dbi1**2+bi2**5*dbi2**2)/120d0/8d0+(bi1**7*dbi1**2+bi2**7*dbi2**2)/720d0/70d0)
        Pi(3)= Le/2d0*(-(dbi1**2+dbi2**2)/3d0+(bi1**2*dbi1**2+bi2**2*dbi2**2)/10d0-(bi1**4*dbi1**2+bi2**4*dbi2**2)/168d0+(bi1**6*dbi1**2+bi2**6*dbi2**2)/120d0/6d0/9d0)
endsubroutine

subroutine Phi_Pi_Kesi(droup,Roup,Phi,Pi,a,da,Kesi,Le)
    use OPERATORS_MOD
    implicit none
        real(8):: droup(3),roup(3),phi(3,3),Pi(3),a(3),da(3),kesi,Le
        real(8):: x,bi1,bi2,dbi1,dbi2
        bi1=a(1)+a(2)
        bi2=a(1)-a(2)
        dbi1=da(1)+da(2)
        dbi2=da(1)-da(2)
        x=kesi/Le
        Roup(1)= Le*( x**2*(a(2))/2d0-x**4*(a(2)**3)/24d0+x**6*(a(2)**5)/720d0-x**8*(a(2)**7)/720d0/56d0+x**10*a(2)**9/720d0/56d0/90d0)
        Roup(2)=-Le/2d0*( x**2*(bi1+bi2)/2d0-x**4*(bi1**3+bi2**3)/24d0+x**6*(bi1**5+bi2**5)/720d0-x**8*(bi1**7+bi2**7)/720d0/56d0+x**10*(bi1**9+bi2**9)/720d0/56d0/90d0)
        Roup(3)= kesi+Le/2d0*(-x**3*(bi1**2+bi2**2)/6d0 +x**5*(bi1**4+bi2**4)/120d0 -x**7*(bi1**6+bi2**6)/120d0/42d0+x**9*(bi1**8+bi2**8)/120d0/42d0/72d0)
        
        Phi=0d0;
        Phi(1,2)= Le*x**2/2d0+Le*(-x**4*(a(2)**2)/8d0+x**6*(a(2)**4)/144d0-x**8*(a(2)**6)/720d0/8d0+x**10*a(2)**8/720d0/56d0/10d0)
        Phi(2,1)=-Le*x**2/2d0-Le/2d0*(-x**4*(bi1**2+bi2**2)/8d0+x**6*(bi1**4+bi2**4)/144d0-x**8*(bi1**6+bi2**6)/720d0/8d0+x**10*(bi1**8+bi2**8)/720d0/56d0/10d0)
        Phi(2,2)=            -Le/2d0*(-x**4*(bi1**2-bi2**2)/8d0+x**6*(bi1**4-bi2**4)/144d0-x**8*(bi1**6-bi2**6)/720d0/8d0+x**10*(bi1**8-bi2**8)/720d0/56d0/10d0)
        Phi(3,1)= Le/2d0*(-x**3*(bi1+bi2)/3d0+x**5*(bi1**3+bi2**3)/30d0-x**7*(bi1**5+bi2**5)/840d0+x**9*(bi1**7+bi2**7)/120d0/42d0/9d0)
        Phi(3,2)= Le/2d0*(-x**3*(bi1-bi2)/3d0+x**5*(bi1**3-bi2**3)/30d0-x**7*(bi1**5-bi2**5)/840d0+x**9*(bi1**7-bi2**7)/120d0/42d0/9d0)

        droup=Phi.mt.da;
        pi=0d0;
        Pi(1)= Le*(-x**4*(a(2))/4d0+x**6*(a(2)**3)/36d0-x**8*(a(2)**5)/120d0/8d0+x**10*a(2)**7/720d0/70d0)*da(2)**2
        Pi(2)=-Le/2d0*(-x**4*(bi1*dbi1**2+bi2*dbi2**2)/4d0+x**6*(bi1**3*dbi1**2+bi2**3*dbi2**2)/36d0-x**8*(bi1**5*dbi1**2+bi2**5*dbi2**2)/120d0/8d0+x**10*(bi1**7*dbi1**2+bi2**7*dbi2**2)/720d0/70d0)
        Pi(3)= Le/2d0*(-x**3*(dbi1**2+dbi2**2)/3d0+x**5*(bi1**2*dbi1**2+bi2**2*dbi2**2)/10d0-x**7*(bi1**4*dbi1**2+bi2**4*dbi2**2)/168d0+x**9*(bi1**6*dbi1**2+bi2**6*dbi2**2)/120d0/6d0/9d0)

endsubroutine






