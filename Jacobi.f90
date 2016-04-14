SUBROUTINE Jacobi(Time,Time_I,y,dy,r,dr,omega,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,a0,lumda_temp,Act_Cn,lumda_temp2,lumda_Quater_temp,Jac)
    use OPERATORS_MOD   
    use BLAS95
    use LAPACK95
    use integration_parameter
    use body_module
    use functions
    !use lin_sol_gen_int
    implicit none
    real(8),allocatable::x(:),dx(:),r(:),dr(:),omega(:),P(:),Mp(:);
    real(8),allocatable::y(:),dy(:),Z_U(:,:),Z1_u(:),Phi_C(:,:),gama_w(:),lumda_temp(:),Phi_C2(:,:),gama_w2(:),lumda_temp2(:),lumda_Quater_temp(:);
    real(8)::Time,Gra(3),a0(3);
    real(8)::Temp_T1(3,3),Temp_T2(3,3),Temp3(3),Temp4(3);
    !integer::F_U(:),F_W(:)
    
    real(8),allocatable::Jac(:,:),Jac_u(:,:),Jac_w(:,:),H_T_P(:,:),Jac_u_temp(:,:),Z1_u_u_temp(:),Z1_Jac_U_temp(:,:),Z1_Jac_w_temp(:,:)
    
    real(8),allocatable::Phi_C_P_all(:,:),Phi_T_P(:,:)
    !real(8)::A1(3,3);
    integer::ele_num_all,F_num_all,con_num_all,Con_num_Act,body_num,II,JJ,KK,LL,MM,NN,OO,I_begin,I_end,Time_I;
    integer,allocatable::Con_num2
    integer,allocatable::FN(:),ipiv(:),Con_DOF(:,:),Con_list(:,:),CON_dof2(:),Act_Cn(:)
    INTERFACE
        SUBROUTINE recur_J(Time,Z,z1,Z1_J,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,a0,Contact_node_num,Contact_node,Body_G,GG,Phi_C_P,Phi_quater,Array,Array_P,Array_PP,Array_G,Mp)  
                  !recur(Time,Z,z1,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Gra,a0,P,Mp)    
            real(8),allocatable:: r(:),dr(:),omega(:),y(:),dy(:),Z(:,:),Z1(:),Z1_J(:,:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),gama_w2(:),Phi_C3(:,:),Body_G(:,:),GG(:,:),Phi_C_P(:,:,:),Phi_quater(:),Array(:,:),Array_P(:,:,:),Array_PP(:,:,:,:),Array_G(:,:,:),Mp(:)
            real(8)::Time,Le,RouA,Iz(3,3),a0(3),A1(3,3),Cm(3,3),Rou
            integer::ele_num,Con_num,CON_dof(:,:),Con_num2,CON_dof2(:),Contact_node_num,Contact_node(:)
        END SUBROUTINE   
    END INTERFACE
    !allocate(Jac_w(con_num_all,F_num_all-con_num_all),Jac_U(con_num_all,con_num_all))
    !allocate(Jac_u_temp(con_num_all,con_num_all),Z1_u_u_temp(con_num_all),Z1_Jac_U_temp(con_num_all,F_num_all),Z1_Jac_w_temp(con_num_all,F_num_all))
    allocate(Phi_C_P_all(F_num_all,F_num_all))
    !allocate(IpIv(con_num_all))
    !
    !
    !allocate(Jacu_T_Fu(Con_num_all))Jac_w_P1(Con_num_all,F_num_all-con_num_all),Jac_w_P2(Con_num_all,F_num_all-con_num_all),Jac_w_P3(Con_num_all,F_num_all-con_num_all),
    allocate(Phi_T_P(F_num_all,F_num_all))
    !allocate(Jac_u_P1(Con_num_all,con_num_all),Jac_u_P2(Con_num_all,con_num_all),Jac_u_P3(Con_num_all,con_num_all),Jacu_u_T_fu1(F_num_all-con_num_all),Jacu_u_T_fu2(F_num_all-con_num_all),Jacu_u_T_fu3(F_num_all-con_num_all))
    allocate(Z_U(F_num_all+body_num+con_num_act,F_num_all+body_num+con_num_act),Z1_u(F_num_all+body_num+con_num_act))    
    
    I_begin=1
    do II=1,body_num
        I_end=I_begin+3*body(II).ele_num+6
        body(II).y_temp=y(I_begin:I_end)
        body(II).dy_temp=dy(I_begin:I_end)
        I_begin=I_end+1        
    enddo
    
    ipiv=0
    Z_U=0d0;
    z1_U=0d0;
    do II=1,Body_num
        body(II).A1=Dir_COS_Q(body(II).y_temp(4:7))
        call recur_J(Time,body(II).Z,body(II).z1,body(II).z1_J,body(II).r_temp,body(II).dr_temp,body(II).omega_temp,body(II).y_temp,body(II).dy_temp,body(II).Le,body(II).A1,body(II).ele_num,body(II).RouA,body(II).Iz,body(II).Cm,body(II).Rou,body(II).Con_num,body(II).CON_dof,body(II).Phi_C,body(II).gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,body(II).Phi_C3,a0,body(II).contact_node_num,body(II).Contact_node,Body(II).G,Body(II).GG,body(II).Phi_C_P,body(II).Phi_quater,body(II).Array,body(II).Array_P,body(II).Array_PP,body(II).Array_G,body(II).MP)
            !recur(Time,Z,z1,r,dr,omega,         y,              dy,              Le,         A1,         ele_num,         RouA,         Iz,         Cm,         Rou,         Con_num,         CON_dof,         Phi_C,         gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Gra,a0,P,Mp)
    enddo
    
    I_begin=1
    do II=1,Body_num   
        I_end=I_begin+body(II).f_num-1
        Z_U(I_begin:I_end,I_begin:I_end)=body(II).Z1_J(body(II).Fn,body(II).Fn)
        Z1_U(I_begin:I_end)=body(II).Z1(body(II).Fn)
        I_begin=I_end+1
    enddo
    JJ=0
    do OO=1,Con_num_all
        if (Con_list(OO,10)==1) then
            JJ=JJ+1
        endif
        if (Con_list(OO,1)==0.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)   

        endif
        
        if (Con_list(OO,1)==1.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-body(KK).Phi_C(MM,body(KK).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)=-body(KK).Phi_C(MM,body(KK).Fn)
        endif 
        
        if (Con_list(OO,1)==2.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)               
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
            z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)-0.3**2
        endif
        
        
        if (Con_list(OO,1)==102.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= body(II).Phi_C(LL,body(II).Fn)
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= body(II).Phi_C(LL,body(II).Fn)
        endif
        
        
        
        if (Con_list(OO,1)==103.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -matmul(transpose(body(II).Array_P(:,body(II).Fn,LL)),body(KK).Array(:,MM))
            
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
            
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -matmul(transpose(body(II).Array_P(:,body(ii).Fn,LL)),body(KK).Array(:,MM))
            
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= -matmul(body(II).Array(:,LL),body(KK).Array_P(:,body(KK).Fn,MM))
            
            !z1_U(f_num_all+body_num+JJ)=dot_product(body(II).Array(:,LL),body(KK).Array(:,MM))
            
        endif
        if (Con_list(OO,1)==104.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
                        
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[body(II).r_temp(1),body(II).r_temp(2),0d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
            
            !z1_U(f_num_all+body_num+JJ)=dot_product(body(II).Array(:,LL),[body(II).r_temp(1),body(II).r_temp(2),0d0])
            
        endif
        if (Con_list(OO,1)==105.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
                        
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= -[0d0,0d0,1d0].mt.body(II).Array_P(:,body(II).Fn,LL)![1d0,0d0,0d0].mt.body(II).Array_P(:,body(II).Fn,LL)!
            
            !z1_U(f_num_all+body_num+JJ)=dot_product(body(II).Array(:,LL),[body(II).r_temp(1),body(II).r_temp(2),0d0])
            
        endif
        
        if (Con_list(OO,1)==3.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            !if (body(II).con_num>0) then
                    
            Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= ((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(II).Phi_C(LL,body(II).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(II).Phi_C3(LL,body(II).Fn))
            Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-((body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(KK).Phi_C(MM,body(KK).Fn)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(KK).Phi_C3(MM,body(KK).Fn))
            Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))
            Z_U(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),f_num_all+body_num+JJ)= Z_U(f_num_all+body_num+JJ,sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))
            
            !z1_U(f_num_all+body_num+JJ)=body(II).gama_w(LL)+body(KK).gama_w(MM)-0.03d0**2!-(body(II).dr_temp(3*Con_list(JJ,3)-2)-body(kk).dr_temp(3*Con_list(JJ,6)-2))**2-(body(II).dr_temp(3*Con_list(JJ,3)-1)-body(kk).dr_temp(3*Con_list(JJ,6)-1))**2
            
        endif
        
    enddo
    
    do ii=1,body_num
       Z_U(ii+F_num_all,sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num)) = body(II).Phi_Quater(body(II).Fn)
       Z_U(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),ii+F_num_all) = body(II).Phi_Quater(body(II).Fn)
       !Z1_U(ii+F_num_all)=dot_product(body(ii).y_temp(4:7),body(ii).y_temp(4:7))-1;
    enddo
    
    JAC=Z_U   
    !----------------------------------------求PHIq矩阵的偏导------------------------------------------!
       
    !Jac_u=jac(F_num_all+1:f_num_all+con_num_all,F_u)
    !Jac_w=jac(F_num_all+1:f_num_all+con_num_all,F_W)
    !
    !
    !Jac_u_temp=transpose(Jac_u)
    !
    !call getrf(Jac_U_temp,ipiv)
    !
    !Jacu_T_Fu=Z1_U(F_u)
    !
    !call getrs(Jac_U_temp,ipiv,Jacu_T_Fu)
    !
    Phi_C_P_all=0d0;
    Phi_T_P=0d0
    JJ=0
    do OO=1,con_num_all  
        if (Con_list(OO,10)==1) then
            JJ=JJ+1
        endif
        Phi_C_P_all=0d0;
        !------------------------------约束方程偏导的偏导-----------------------------------------------!
        if (Con_list(OO,1)==0.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)                    
            Phi_C_P_all(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=  body(II).Phi_C_P(LL,body(II).Fn,body(II).Fn)
        
        endif
        if (Con_list(OO,1)==1.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2) 
            LL=Con_list(OO,8)
            Phi_C_P_all(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=  body(II).Phi_C_P(LL,body(II).Fn,body(II).Fn)
   
            II=Con_list(OO,5)              
            LL=Con_list(OO,9)
            Phi_C_P_all(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= -body(II).Phi_C_P(LL,body(II).Fn,body(II).Fn)
        endif
        if (Con_list(OO,1)==2.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            LL=Con_list(OO,8)                    
            Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=  body(II).Phi_C_P(LL,body(II).Fn,body(II).Fn)
        endif
        if (Con_list(OO,1)==103.and.Con_list(OO,10)==1) then
                II=Con_list(OO,2)
                KK=Con_list(OO,5)
                LL=Con_list(OO,8)
                MM=Con_list(OO,9)
                call Tensor(Temp_T1,body(II).Array(:,LL))
                call Tensor(Temp_T2,body(KK).Array(:,MM))
                Temp3=body(II).Array(:,LL)
                Temp4=body(KK).Array(:,MM)
                Phi_C_P_All(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= transpose(body(KK).Array_P(:,body(KK).Fn,mm)).mt.body(II).Array_P(:,body(II).Fn,LL)
                Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))= transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.body(KK).Array_P(:,body(KK).Fn,mm)
                Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-transpose(transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.Temp_T2.mt.body(II).Array_G(:,body(II).Fn,LL))-(Temp4(1)*body(II).Array_PP(1,body(II).Fn,body(II).Fn,LL)+Temp4(2)*body(II).Array_PP(2,body(II).Fn,body(II).Fn,LL)+Temp4(3)*body(II).Array_PP(3,body(II).Fn,body(II).Fn,LL))
                Phi_C_P_All(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-transpose(transpose(body(KK).Array_P(:,body(KK).Fn,mm)).mt.Temp_T1.mt.body(KK).Array_G(:,body(KK).Fn,mm))-(Temp3(1)*body(KK).Array_PP(1,body(KK).Fn,body(KK).Fn,MM)+Temp3(2)*body(KK).Array_PP(2,body(KK).Fn,body(KK).Fn,MM)+Temp3(3)*body(KK).Array_PP(3,body(KK).Fn,body(KK).Fn,MM))
            
        endif        
        
        if (Con_list(OO,1)==104.and.Con_list(OO,10)==1) then
                II=Con_list(OO,2)
                LL=Con_list(OO,8)
                !call Tensor(Temp_T1,body(II).Array(:,LL))
                call Tensor(Temp_T2,[body(II).r_temp(1),body(II).r_temp(2),0d0])               
                Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-transpose(transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.Temp_T2.mt.body(II).Array_G(:,body(II).Fn,LL))-(body(II).r_temp(1)*body(II).Array_PP(1,body(II).Fn,body(II).Fn,LL)+body(II).r_temp(2)*body(II).Array_PP(2,body(II).Fn,body(II).Fn,LL))
                !call Tensor(Temp_T2,[1d0,0d0,0d0])
                !Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-transpose(transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.Temp_T2.mt.body(II).Array_G(:,body(II).Fn,LL))-(body(II).Array_PP(1,body(II).Fn,body(II).Fn,LL))
        endif        
         
        if (Con_list(OO,1)==105.and.Con_list(OO,10)==1) then
                II=Con_list(OO,2)
                LL=Con_list(OO,8)
                !call Tensor(Temp_T1,body(II).Array(:,LL))
                call Tensor(Temp_T2,[0d0,0d0,1d0])               
                Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-transpose(transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.Temp_T2.mt.body(II).Array_G(:,body(II).Fn,LL))-body(II).Array_PP(3,body(II).Fn,body(II).Fn,LL)
                !call Tensor(Temp_T2,[1d0,0d0,0d0])
                !Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-transpose(transpose(body(II).Array_P(:,body(II).Fn,LL)).mt.Temp_T2.mt.body(II).Array_G(:,body(II).Fn,LL))-(body(II).Array_PP(1,body(II).Fn,body(II).Fn,LL))
        endif
        
        if (Con_list(OO,1)==3.and.Con_list(OO,10)==1) then
            II=Con_list(OO,2)
            KK=Con_list(OO,5)
            LL=Con_list(OO,8)
            MM=Con_list(OO,9)
            Phi_C_P_All(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))=-2d0*( (body(KK).Array_G(1,body(KK).Fn,MM).mt.body(II).Array_G(1,body(II).Fn,LL))+(body(KK).Array_G(2,body(KK).Fn,MM).mt.body(II).Array_G(2,body(II).Fn,LL)))
            Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-2d0*( (body(II).Array_G(1,body(II).Fn,LL).mt.body(KK).Array_G(1,body(KK).Fn,MM))+(body(II).Array_G(2,body(II).Fn,LL).mt.body(KK).Array_G(2,body(KK).Fn,MM)))
            Phi_C_P_All(sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num),sum(body(1:II-1).F_num)+1:sum(body(1:II).F_num))= 2d0*( (body(II).Array_G(1,body(II).Fn,LL).mt.body(II).Array_G(1,body(II).Fn,LL))+(body(II).Array_G(2,body(II).Fn,LL).mt.body(II).Array_G(2,body(II).Fn,LL))+(body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(II).Array_PP(1,body(II).Fn,body(II).Fn,LL)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(II).Array_PP(2,body(II).Fn,body(II).Fn,LL))
            Phi_C_P_All(sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num),sum(body(1:KK-1).F_num)+1:sum(body(1:KK).F_num))=-2d0*(-(body(KK).Array_G(1,body(KK).Fn,MM).mt.body(KK).Array_G(1,body(KK).Fn,MM))-(body(KK).Array_G(2,body(KK).Fn,MM).mt.body(II).Array_G(2,body(KK).Fn,MM))+(body(II).r_temp(3*Con_list(OO,3)-2)-body(kk).r_temp(3*Con_list(OO,6)-2))*body(KK).Array_PP(1,body(KK).Fn,body(KK).Fn,MM)+(body(II).r_temp(3*Con_list(OO,3)-1)-body(kk).r_temp(3*Con_list(OO,6)-1))*body(KK).Array_PP(2,body(KK).Fn,body(KK).Fn,MM))
            
        endif
         
        if (Con_list(OO,10)==1) then
            Phi_T_P=Phi_T_P+Phi_C_P_all*lumda_temp(OO); 
        endif
    enddo            
            
    !------------------------------约束方程偏导的偏导-----------------------------------------------!
    ! 
    do ii=1,body_num
        Jac(sum(body(1:II-1).F_num)+4,sum(body(1:II-1).F_num)+4)= Jac(sum(body(1:II-1).F_num)+4,sum(body(1:II-1).F_num)+4)+2d0*lumda_quater_temp(II)
        Jac(sum(body(1:II-1).F_num)+5,sum(body(1:II-1).F_num)+5)= Jac(sum(body(1:II-1).F_num)+5,sum(body(1:II-1).F_num)+5)+2d0*lumda_quater_temp(II)
        Jac(sum(body(1:II-1).F_num)+6,sum(body(1:II-1).F_num)+6)= Jac(sum(body(1:II-1).F_num)+6,sum(body(1:II-1).F_num)+6)+2d0*lumda_quater_temp(II)
        Jac(sum(body(1:II-1).F_num)+7,sum(body(1:II-1).F_num)+7)= Jac(sum(body(1:II-1).F_num)+7,sum(body(1:II-1).F_num)+7)+2d0*lumda_quater_temp(II)
    enddo
    !
    
    
    Jac(1:f_num_all,1:f_num_all)= Jac(1:f_num_all,1:f_num_all)+Phi_T_P;

    
    
    !call exportf(H_T_P,'H_T_P.txt')
    !
    !call exportf(body(1).Phi_C_P1(:,:,1),'Phi_C_P11.txt')
    !call exportf(body(1).Phi_C_P1(:,:,2),'Phi_C_P12.txt')
    !call exportf(body(1).Phi_C_P1(:,:,3),'Phi_C_P13.txt')
    !call exportf(Z1_u,'Z1_u3.txt')
    !call potrf(Z_U)
    !call potrs(Z_U,z1_U)
    !call getrf(Z_U,ipiv)
    !call getrs(Z_U,ipiv,Z1_U)
    !dx(3*ele_num_all+13:6*ele_num_all+12)=z1_U;
    !call exportf(Z_U,'Z3.txt')
    !call exportf(Z1,'Z13.txt')
    !call exportf(Z1_U,'Z1u3.txt')
    !call exportf(Jac,'JAC.txt')
    !call exportf(r,'r.txt')
endsubroutine Jacobi

subroutine recur_J(Time,Z,z1,z1_J,r,dr,omega,y,dy,Le,A1,ele_num,RouA,Iz,Cm,Rou,Con_num,CON_dof,Phi_C,gama_w,Con_num2,CON_dof2,Phi_C2,gama_w2,Phi_C3,a0,Contact_node_num,Contact_node,Body_G,GG,Phi_C_P,Phi_quater,Array,Array_P,Array_PP,Array_G,Mp)
!    use body_module
!    use Element_module
!    use Model_module
    use functions
    use OPERATORS_MOD
    implicit none
    real(8),allocatable::r(:),dr(:),omega(:),y(:),dy(:)
    real(8),allocatable::Z(:,:),z1(:),z1_J(:,:),Phi_C(:,:),gama_w(:),Phi_C2(:,:),Phi_C3(:,:),gama_w2(:),Body_G(:,:),Mp(:);

    real(8)::Time,Le,RouA,Iz(3,3),a0(3),Cm(3,3),Rou;
    real(8)::A1(3,3)
    integer::ele_num,Con_num,Con_num2,contact_node_num;
    integer,allocatable::CON_dof(:,:),Con_dof2(:),contact_node(:);
    
    integer:: ii,kk,jj,ii_contact,LL
    real(8):: Eye(3,3),A2(3,3)
    real(8):: U(9,3),Roup_T(3,3),omega_T(3,3),Me(9,9),Fe(9),Fe_P(9,3),Fe_c(9,3),Fec(6)
    real(8):: phi(3,3),Pi(3),Ki(3,3)
    
    real(8)::R_Q(3,4)
    real(8),allocatable::Phi_quater(:)
    
    real(8),allocatable:: G(:,:),g0(:,:),ones(:),M(:,:),F(:),F_P(:,:)
    real(8),allocatable:: GG(:,:),gg0(:,:)!为末端点约束而建
    real(8),allocatable:: a(:),da(:),F1(:),T(:,:),Roup(:,:),dRoup(:,:),beta(:,:)
    
    !-------------------------rou，phi，Ki偏导数------------20150318------------------!
    real(8):: Roup_temp(3);
    real(8):: phi_P1(3,3),phi_P2(3,3),phi_P3(3,3)
    real(8):: Ki_P1(3,3),Ki_P2(3,3),Ki_P3(3,3)
    real(8):: Roup_P1_T(3,3),roup_P2_T(3,3),roup_P3_T(3,3)
    real(8):: Phi0(3,3)
    real(8):: A2_P1(3,3),A2_P2(3,3),A2_P3(3,3)
    real(8):: A0_P1(3,3),A0_P2(3,3),A0_P3(3,3),A0_P0(3,3)
    real(8):: K1_P1(3,3),K1_P2(3,3),K1_P3(3,3)
    real(8),allocatable:: G_P1(:,:,:), G_P2(:,:,:), G_P3(:,:,:),G_P0(:,:),T_P1(:,:),T_P2(:,:),T_P3(:,:),T_P0(:,:)
    real(8),allocatable:: A1_P1(:,:),A1_P2(:,:),A1_P3(:,:),A1_P0(:,:)
    real(8),allocatable:: Array(:,:),Array_P(:,:,:),AI(:,:),Array_G(:,:,:)
    
    real(8),allocatable:: A1_P0P0(:,:,:),A1_P0P1(:,:,:),A1_P0P2(:,:,:),A1_P0P3(:,:,:),A1_P1P1(:,:,:),A1_P1P2(:,:,:),A1_P1P3(:,:,:),A1_P2P2(:,:,:),A1_P2P3(:,:,:),A1_P3P3(:,:,:)
    real(8):: A2_P1P1(3,3),A2_P1P2(3,3),A2_P1P3(3,3),A2_P2P2(3,3),A2_P2P3(3,3),A2_P3P3(3,3)
    real(8),allocatable:: Array_PP(:,:,:,:)
    
    real(8),allocatable:: GG_P1(:,:,:), GG_P2(:,:,:), GG_P3(:,:,:), GG_P0(:,:)
    real(8),allocatable:: Phi_C_P(:,:,:)!0(:,:),Phi_C_P1(:,:,:), Phi_C_P2(:,:,:), Phi_C_P3(:,:,:)
    !-------------------------rou，phi，Ki偏导数------------20150318------------------!!20150322加约束方程的二阶雅可比
    
    !----------------------Analytical Jacobi----------------20150318------------------!!20150322加约束方程的二阶雅可比
    
    real(8)::Temp_T(3,3)
    !-----------------------------------------------------------------------------20150714
    real(8)::S1(3,3),S2(3,3),S3(3,3);
    real(8)::H1(3,6),H2(3,6),H3(3,6),H4(3,6),H5(3,6),H6(3,6),H7(3,6),H8(3,6),H9(3,6)
    !real(8)::T_PP1(9,9),T_PP2(9,9),T_PP3(9,9)
    real(8),allocatable::T_PP1(:,:),T_PP2(:,:),T_PP3(:,:)
    !
    
    allocate( T_PP1(9*ele_num+9,9*ele_num+9),T_PP2(9*ele_num+9,9*ele_num+9),T_PP3(9*ele_num+9,9*ele_num+9))
    allocate( G_P1(9*ele_num,3*ele_num+7,ele_num+1),G_P2(9*ele_num,3*ele_num+7,ele_num+1),G_P3(9*ele_num,3*ele_num+7,ele_num+1),G_P0(9*ele_num,3*ele_num+7))
    allocate( T_P1(9*ele_num+9,9*ele_num+9),T_P2(9*ele_num+9,9*ele_num+9),T_P3(9*ele_num+9,9*ele_num+9),T_P0(9*ele_num+9,9))
    allocate( A1_P1(3*ele_num+3,3*ele_num+3),A1_P2(3*ele_num+3,3*ele_num+3),A1_P3(3*ele_num+3,3*ele_num+3),A1_P0(3*ele_num+3,3))
    
    allocate( A1_P0P0(3,3,ele_num+1),A1_P0P1(3*ele_num+3,3,ele_num+1),A1_P0P2(3*ele_num+3,3,ele_num+1),A1_P0P3(3*ele_num+3,3,ele_num+1))
    allocate( A1_P1P1(3*ele_num+3,3*ele_num+3,ele_num+1),A1_P1P2(3*ele_num+3,3*ele_num+3,ele_num+1),A1_P1P3(3*ele_num+3,3*ele_num+3,ele_num+1),A1_P2P2(3*ele_num+3,3*ele_num+3,ele_num+1),A1_P2P3(3*ele_num+3,3*ele_num+3,ele_num+1),A1_P3P3(3*ele_num+3,3*ele_num+3,ele_num+1))
    allocate( AI(3,3*ele_num+3))
    
    allocate( GG_P1(9,3*ele_num+7,ele_num+1),GG_P2(9,3*ele_num+7,ele_num+1),GG_P3(9,3*ele_num+7,ele_num+1),GG_P0(9,3*ele_num+7))

    !---------------------Analytical Jacobi----------------20150318------------------!!20150322加约束方程的二阶雅可比
    
    allocate( G(9*ele_num,3*ele_num+7),g0(9*ele_num,3*ele_num),M(9*ele_num,9*ele_num),F(9*ele_num),F_P(9*ele_num,3*ele_num+7),F1(3*ele_num+6),beta(9,3*ele_num+1))
    allocate( gg0(6,ele_num))
    allocate( T(9*ele_num+9,1:9),a(3*ele_num),da(3*ele_num),Roup(3,ele_num),droup(3,ele_num),ones(3*ele_num))
    
    ones=1.d0;
    !!
 !   write(*,*) ones
    U=0d0; U(7,1)=1d0;U(8,2)=1d0;U(9,3)=1d0;
    
    R_Q=2d0*R_Quater(y(4:7))
    !------------------------------------------！
    S1=0d0;S2=0d0;S3=0d0
    S1(2,3)=-1d0;S1(3,2)= 1d0
    S2(1,3)= 1d0;S2(3,1)=-1d0
    S3(1,2)=-1d0;S3(2,1)= 1d0
    T_PP1=0d0;T_PP2=0d0;T_PP3=0d0;
    H1=0d0;H2=0d0;H3=0d0;H4=0d0;H5=0d0;h6=0d0;h7=0d0;H8=0d0;H9=0d0;
    !------------------------------------------!
    G_P1=0d0
    G_P2=0d0
    G_P3=0d0
    G_P0=0d0
    GG_P1=0d0
    GG_P2=0d0
    GG_P3=0d0
    GG_P0=0d0
    
    G_P1(4,4,1)=-2d0;G_P1(5,7,1)=-2d0;G_P1(6,6,1)= 2d0;
    
    G_P2(4,7,1)= 2d0;G_P2(5,4,1)=-2d0;G_P2(6,5,1)=-2d0;
    
    G_P3(4,6,1)=-2d0;G_P3(5,5,1)= 2d0;G_P3(6,4,1)=-2d0;
    
    G_P0(4,5)=2d0;G_P0(5,6)=2d0;G_P0(6,7)=2d0;
        
    T_P0=0d0;T_P1=0d0;T_P2=0d0;T_P3=0d0;
    F_P=0d0;
    
    !------------------------------------------!
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
    
    Ai(1:3,1:3)=A1;
    Array=0d0;Array_P=0d0
    !--------------根节点方向余弦，对欧拉四元数的导数-----------------!
    
    !--------------根节点方向余弦，对欧拉四元数的二阶导数-----------------!
    !A1_P0P0=0d0;A1_P0P1=0d0;A1_P0P2=0d0;A1_P0P3=0d0;
    !A1_P0P0(1,1,1)= 4d0; A1_P0P0(2,2,1)= 4d0; A1_P0P0(3,3,1)=4d0;
    !A1_P0P1(2,3,1)=-2d0; A1_P0P1(3,2,1)= 2d0;
    !A1_P0P2(1,3,1)= 2d0; A1_P0P2(3,1,1)=-2d0;
    !A1_P0P3(1,2,1)=-2d0; A1_P0P3(2,1,1)= 2d0;
    !
    !A1_P1P1=0d0;A1_P1P2=0d0;A1_P1P3=0d0
    !A1_P1P1(1,1,1)= 4d0;
    !A1_P1P2(1,2,1)= 2d0; A1_P1P2(2,1,1)= 2d0;
    !A1_P1P3(1,3,1)= 2d0; A1_P1P3(3,1,1)= 2d0;
    !          
    !A1_P2P2=0d0;A1_P2P3=0d0;A1_P3P3=0d0
    !A1_P2P2(2,2,1)= 4d0;
    !A1_P2P3(2,3,1)= 2d0;A1_P2P3(3,2,1)= 2d0;
    !A1_P3P3(3,3,1)= 4d0;
    
    Array_PP=0d0;
    !--------------根节点方向余弦，对欧拉四元数的二阶导数-----------------!
  
    
    G=0d0;
    G(1,1)=1d0;G(2,2)=1d0;G(3,3)=1d0;
    G(4:6,4:7)=R_Q;
    G(7,8)=1d0;G(8,9)=1d0;G(9,10)=1d0;
    g0=0d0;
 !  write(*,*) ones
    r(1:3)=y(1:3);
    dr(1:3)=dy(1:3);
    omega(1:3)=R_Q.mt.dy(4:7);
    a = y(8:3*ele_num+7);
    da=dy(8:3*ele_num+7);
    Phi=0d0;
    Z=0d0;
    Z1=0d0;
    Z1_J=0d0;
    F1=0d0;
    T=0d0;

    
    Eye=0d0;Eye(1,1)=1d0;Eye(2,2)=1d0;Eye(3,3)=1d0;
    M=0d0;
    F=0d0;
    F_P=0d0
    
    beta=0d0
    Roup_temp=0d0;
    do ii=1,ele_num
        call Phi_Pi(droup(:,ii),Roup_temp,PHI0,Pi,a(3*ii-2:3*ii),da(3*ii-2:3*ii),Le)
        call Phi_P(Phi_P1,Phi_P2,Phi_P3,a(3*ii-2:3*ii),Le)
        Roup(:,ii)=matmul(A1,Roup_temp)
        PHI=matmul(A1,Phi0)
        r(3*ii+1:3*ii+3)=r(3*ii-2:3*ii)+Roup(:,ii); 
         Fec=0d0;
        !if (contact_node_num>0.and.ii<ele_num) then
        !    !if (contact_node_num==1) then
        !    !    if (r(3*ii+3)>1.06d0.and.r(3*ii+3)<1.51d0) then!.and.mod(ii,10)==0
        !    !        Fec(1:3)=-[1d4*(r(3*ii+1)-0.3d0),1d4*(r(3*ii+2)),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];
        !    !    endif
        !    !endif
        !    !if (contact_node_num==2) then
        !    !    if (r(3*ii+3)>1.06d0.and.r(3*ii+3)<1.51d0) then!.and.mod(ii,10)==0
        !    !        Fec(1:3)=-[1d4*(r(3*ii+1)+0.15d0),1d4*(r(3*ii+2)-0.26d0),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];;
        !    !    endif
        !    !endif
        !    !if (contact_node_num==3) then
        !    !    if (r(3*ii+3)>1.06d0.and.r(3*ii+3)<1.51d0) then!.and.mod(ii,10)==0
        !    !        Fec(1:3)=-[1d4*(r(3*ii+1)+0.15d0),1d4*(r(3*ii+2)+0.26d0),0d0]!+[10d0*(dr(3*ii+1)),10d0*(dr(3*ii+2)),0d0];;
        !    !    endif
        !    !endif
        !    if (sqrt(r(3*ii+1)**2+r(3*ii+2)**2)<0.2990d0) then
        !        Fec(1:3)=Fec(1:3)+0.5d-4*[r(3*ii+1),r(3*ii+2),0d0];
        !    endif
        !    !if (r(3*ii+3)<0d0) then
        !    !    Fec(1:3)=Fec(1:3)-0.075d0*[0d0,0d0,r(3*ii+3)];
        !    !endif
        !endif
        
        A2=Dir_COS(a(3*ii-2:3*ii))

        !--------------单元相对转动，产生的方向余弦对本单元相对转角的偏导数------20150405------------------------!
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
        
        !---------------------------------方向余旋的两阶导数--------------------------------
        
        !A2_P1P1=Dir_cos_P1P1(a(3*ii-2:3*ii))
        !A2_P1P2=Dir_cos_P1P2(a(3*ii-2:3*ii))
        !A2_P1P3=Dir_cos_P1P3(a(3*ii-2:3*ii))
        !A2_P2P2=Dir_cos_P2P2(a(3*ii-2:3*ii))
        !A2_P2P3=Dir_cos_P2P3(a(3*ii-2:3*ii))
        !A2_P3P3=Dir_cos_P3P3(a(3*ii-2:3*ii))
        !       
        !A1_P0P0(1:3,1:3,ii+1)=A1_P0P0(1:3,1:3,ii).mt.A2          
        !do jj=0,ii  
        !    if(jj==ii) then
        !        A1_P0P1(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0(3*ii-2:3*ii,1:3).mt.A2_P1
        !        A1_P0P2(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0(3*ii-2:3*ii,1:3).mt.A2_P2
        !        A1_P0P3(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0(3*ii-2:3*ii,1:3).mt.A2_P3
        !    else                    
        !        A1_P0P1(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0P1(3*jj+1:3*jj+3,1:3,ii).mt.A2
        !        A1_P0P2(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0P2(3*jj+1:3*jj+3,1:3,ii).mt.A2
        !        A1_P0P3(3*jj+1:3*jj+3,1:3,ii+1)=A1_P0P3(3*jj+1:3*jj+3,1:3,ii).mt.A2 
        !    endif
        !     do kk=0,ii
        !        if(jj==ii.and.kk==ii) then
        !            A1_P1P1(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P1P1.mt.A1
        !            A1_P1P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P1P2.mt.A1
        !            A1_P1P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P1P3.mt.A1
        !            A1_P2P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P2P2.mt.A1
        !            A1_P2P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P2P3.mt.A1
        !            A1_P3P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A2_P3P3.mt.A1
        !        elseif(jj==ii.and.kk<ii) then
        !            A1_P1P1(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P1
        !            A1_P1P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P2
        !            A1_P1P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P3
        !            A1_P2P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P2
        !            A1_P2P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P3
        !            A1_P3P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P3(3*ii-2:3*ii,3*kk+1:3*kk+3).mt.A2_P3
        !        elseif(ii==kk.and.jj<ii) then
        !            A1_P1P1(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P1
        !            A1_P1P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P2
        !            A1_P1P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P3
        !            A1_P2P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P2
        !            A1_P2P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P3
        !            A1_P3P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P3(3*ii-2:3*ii,3*jj+1:3*jj+3).mt.A2_P3                
        !        elseif(ii<kk.and.jj<ii) then      
        !            A1_P1P1(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1P1(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !            A1_P1P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !            A1_P1P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P1P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !            A1_P2P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2P2(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !            A1_P2P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P2P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !            A1_P3P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii+1)=A1_P3P3(3*jj+1:3*jj+3,3*kk+1:3*kk+3,ii).mt.A2
        !        endif
        !    enddo
        !enddo                             
        !---------------------------------方向余旋的两阶导数--------------------------------
        
        Ki=K_R(a(3*ii-2:3*ii))
        Ki_P1=K_R_P1(a(3*ii-2:3*ii))
        Ki_P2=K_R_P2(a(3*ii-2:3*ii))
        Ki_P3=0d0
        !
        call Tensor(Roup_T,Roup(:,ii));
        
        T(9*ii+1:9*ii+3,1:9)=reshape([eye,-Roup_T,phi],[3,9]);
        T(9*ii+4:9*ii+6,4:6)=Eye;
        T(9*ii+4:9*ii+6,7:9)=A1.mt.Ki;
        !
        call Tensor(Roup_P1_T,PHI(:,1))
        call Tensor(Roup_P2_T,PHI(:,2))
        call Tensor(Roup_P3_T,PHI(:,3))
        !
        !
        !
        !!---------------求T的偏导数-------------!
        call tensor(T_P0(9*ii+1:9*ii+3,4:6),-matmul(A1_P0(3*ii-2:3*ii,1:3),Roup_temp))
        T_P0(9*ii+1:9*ii+3,7:9)=matmul(A1_P0(3*ii-2:3*ii,1:3),Phi0)
        T_P0(9*ii+4:9*ii+6,7:9)=matmul(A1_P0(3*ii-2:3*ii,1:3),Ki)
        
        do kk=0,ii
            if(kk<ii) then              
                call tensor(T_P1(9*ii+1:9*ii+3,9*kk+4:9*kk+6),-matmul(A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3),Roup_temp))
                            T_P1(9*ii+1:9*ii+3,9*kk+7:9*kk+9)= matmul(A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3),Phi0)
                            T_P1(9*ii+4:9*ii+6,9*kk+7:9*kk+9)= matmul(A1_P1(3*ii-2:3*ii,3*kk+1:3*kk+3),Ki)
                
                call tensor(T_P2(9*ii+1:9*ii+3,9*kk+4:9*kk+6),-matmul(A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3),Roup_temp))
                            T_P2(9*ii+1:9*ii+3,9*kk+7:9*kk+9)= matmul(A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3),Phi0)
                            T_P2(9*ii+4:9*ii+6,9*kk+7:9*kk+9)= matmul(A1_P2(3*ii-2:3*ii,3*kk+1:3*kk+3),Ki)
                
                call tensor(T_P3(9*ii+1:9*ii+3,9*kk+4:9*kk+6),-matmul(A1_P3(3*ii-2:3*ii,3*kk+1:3*kk+3),Roup_temp))
                            T_P3(9*ii+1:9*ii+3,9*kk+7:9*kk+9)= matmul(A1_P3(3*ii-2:3*ii,3*kk+1:3*kk+3),Phi0)
                            T_P3(9*ii+4:9*ii+6,9*kk+7:9*kk+9)= matmul(A1_P3(3*ii-2:3*ii,3*kk+1:3*kk+3),Ki)
            else                
                T_P1(9*ii+1:9*ii+3,9*kk+4:9*kk+6)=-(Roup_P1_T)
                T_P1(9*ii+1:9*ii+3,9*kk+7:9*kk+9)= matmul(A1,Phi_P1)
                T_P1(9*ii+4:9*ii+6,9*kk+7:9*kk+9)= matmul(A1,Ki_P1)
                
                T_P2(9*ii+1:9*ii+3,9*kk+4:9*kk+6)=-(Roup_P2_T)
                T_P2(9*ii+1:9*ii+3,9*kk+7:9*kk+9)= matmul(A1,Phi_P2)
                T_P2(9*ii+4:9*ii+6,9*kk+7:9*kk+9)= matmul(A1,Ki_P2)
                
                T_P3(9*ii+1:9*ii+3,9*kk+4:9*kk+6)=-(Roup_P3_T)
            endif
        enddo    
        !---------------求T的偏导数-------------!
        
        !call Tensor(omega_T,omega(3*ii-2:3*ii));
        H1(:,1:3)=-S1.mt.Roup_T
        H1(:,4:6)= S1.mt.Phi
        H2(:,1:3)=-S2.mt.Roup_T
        H2(:,4:6)= S2.mt.Phi
        H3(:,1:3)=-S3.mt.Roup_T
        H3(:,4:6)= S3.mt.Phi
        
        H4(:,1:3)=-Roup_P1_T
        H4(:,4)=A1.mt.PHi_P1(:,1);H4(:,5)=A1.mt.PHi_P2(:,1);H4(:,6)=A1.mt.PHi_P3(:,1);
        H5(:,1:3)=-Roup_P2_T
        H5(:,4)=A1.mt.PHi_P1(:,2);H5(:,5)=A1.mt.PHi_P2(:,2);H5(:,6)=A1.mt.PHi_P3(:,2);
        H6(:,1:3)=-Roup_P3_T
        H6(:,4)=A1.mt.PHi_P1(:,3);H6(:,5)=A1.mt.PHi_P2(:,3);H6(:,6)=A1.mt.PHi_P3(:,3);
        
        call tensor(H7(:,1:3),-A1.mt.Ki(:,1))
        H7(:,4)=A1.mt.Ki_P1(:,1);H7(:,5)=A1.mt.Ki_P2(:,1);H7(:,6)=A1.mt.Ki_P3(:,1);
        call tensor(H8(:,1:3),-A1.mt.Ki(:,2))
        H8(:,4)=A1.mt.Ki_P1(:,2);H8(:,5)=A1.mt.Ki_P2(:,2);H8(:,6)=A1.mt.Ki_P3(:,2);
        call tensor(H9(:,1:3),-A1.mt.Ki(:,3))
        H9(:,4)=A1.mt.Ki_P1(:,3);H9(:,5)=A1.mt.Ki_P2(:,3);H9(:,6)=A1.mt.Ki_P3(:,3);
        do kk=0,ii
            T_PP1(9*ii+1:9*ii+3,9*kk+4)=H1.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+1:9*ii+3,9*kk+5)=H2.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+1:9*ii+3,9*kk+6)=H3.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+1:9*ii+3,9*kk+7)=H4.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+1:9*ii+3,9*kk+8)=H5.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+1:9*ii+3,9*kk+9)=H6.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+4:9*ii+6,9*kk+7)=H7.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+4:9*ii+6,9*kk+8)=H8.mt.G(9*ii-5:9*ii,3*kk+5)
            T_PP1(9*ii+4:9*ii+6,9*kk+9)=H9.mt.G(9*ii-5:9*ii,3*kk+5)
            
            T_PP2(9*ii+1:9*ii+3,9*kk+4)=H1.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+1:9*ii+3,9*kk+5)=H2.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+1:9*ii+3,9*kk+6)=H3.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+1:9*ii+3,9*kk+7)=H4.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+1:9*ii+3,9*kk+8)=H5.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+1:9*ii+3,9*kk+9)=H6.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+4:9*ii+6,9*kk+7)=H7.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+4:9*ii+6,9*kk+8)=H8.mt.G(9*ii-5:9*ii,3*kk+6)
            T_PP2(9*ii+4:9*ii+6,9*kk+9)=H9.mt.G(9*ii-5:9*ii,3*kk+6)
            
            T_PP3(9*ii+1:9*ii+3,9*kk+4)=H1.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+1:9*ii+3,9*kk+5)=H2.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+1:9*ii+3,9*kk+6)=H3.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+1:9*ii+3,9*kk+7)=H4.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+1:9*ii+3,9*kk+8)=H5.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+1:9*ii+3,9*kk+9)=H6.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+4:9*ii+6,9*kk+7)=H7.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+4:9*ii+6,9*kk+8)=H8.mt.G(9*ii-5:9*ii,3*kk+7)
            T_PP3(9*ii+4:9*ii+6,9*kk+9)=H9.mt.G(9*ii-5:9*ii,3*kk+7)
            
        enddo 
        

        !T_PP1(1:3,4)=H1.mt.G(9*ii-5:9*ii,8)
        !T_PP1(1:3,5)=H2.mt.G(9*ii-5:9*ii,8)
        !T_PP1(1:3,6)=H3.mt.G(9*ii-5:9*ii,8)
        !T_PP1(1:3,7)=H4.mt.G(9*ii-5:9*ii,8)
        !T_PP1(1:3,8)=H5.mt.G(9*ii-5:9*ii,8)
        !T_PP1(1:3,9)=H6.mt.G(9*ii-5:9*ii,8)
        !T_PP1(4:6,7)=H7.mt.G(9*ii-5:9*ii,8)
        !T_PP1(4:6,8)=H8.mt.G(9*ii-5:9*ii,8)
        !T_PP1(4:6,9)=H9.mt.G(9*ii-5:9*ii,8)       
        !
        !T_PP2(1:3,4)=H1.mt.G(9*ii-5:9*ii,9)
        !T_PP2(1:3,5)=H2.mt.G(9*ii-5:9*ii,9)
        !T_PP2(1:3,6)=H3.mt.G(9*ii-5:9*ii,9)
        !T_PP2(1:3,7)=H4.mt.G(9*ii-5:9*ii,9)
        !T_PP2(1:3,8)=H5.mt.G(9*ii-5:9*ii,9)
        !T_PP2(1:3,9)=H6.mt.G(9*ii-5:9*ii,9)
        !T_PP2(4:6,7)=H7.mt.G(9*ii-5:9*ii,9)
        !T_PP2(4:6,8)=H8.mt.G(9*ii-5:9*ii,9)
        !T_PP2(4:6,9)=H9.mt.G(9*ii-5:9*ii,9)       
        !
        !T_PP3(1:3,4)=H1.mt.G(9*ii-5:9*ii,10)
        !T_PP3(1:3,5)=H2.mt.G(9*ii-5:9*ii,10)
        !T_PP3(1:3,6)=H3.mt.G(9*ii-5:9*ii,10)
        !T_PP3(1:3,7)=H4.mt.G(9*ii-5:9*ii,10)
        !T_PP3(1:3,8)=H5.mt.G(9*ii-5:9*ii,10)
        !T_PP3(1:3,9)=H6.mt.G(9*ii-5:9*ii,10)
        !T_PP3(4:6,7)=H7.mt.G(9*ii-5:9*ii,10)
        !T_PP3(4:6,8)=H8.mt.G(9*ii-5:9*ii,10)
        !T_PP3(4:6,9)=H9.mt.G(9*ii-5:9*ii,10)         
        

        call Element_MF_J(Me,Fe,Fe_P,Mp(6*ii-5:6*ii)+Fec,a(3*ii-2:3*ii),a0,RouA,Le,A1,omega(3*ii-2:3*ii),da(3*ii-2:3*ii),Iz,Cm,Rou,Phi,Ki,Roup(:,ii))
        
        M(9*ii-8:9*ii,9*ii-8:9*ii)=Me;
        F_P(9*ii-8:9*ii,3*ii+5:3*ii+7)=Fe_P!+Fe_c;
        F(9*ii-8:9*ii)=Fe
        A1=A1.mt.A2    !        
        Ai(1:3,3*ii+1:3*ii+3)=A1;
        do kk=1,ii+1
            if (kk<ii+1.and.ii<ele_num) then
                if (kk==1) then
                    G(9*ii+1:9*ii+9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
                else
                    G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
                endif
            elseif (kk==ii+1.and.ii<ele_num) then
                G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=U;
            endif;
            
            if (ii==ele_num) then
                if (kk==1) then
                    GG(1:9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
                elseif(1<kk.and.kk<ii+1) then
                    GG(1:9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
                endif
            endif
           
        enddo
    
    enddo
    !--------------------传递函数G的偏导数--------------------!
      !!G(ii,jj)对q(kk)的偏导数!    
    
      T_P1=T_PP1;
      T_P2=T_PP2;
      T_P3=T_PP3;
    do ii=1,ele_num           
        do jj=1,ele_num    
            if(ii<ele_num) then
                if (jj==1) then
                    G_P0(9*ii+1:9*ii+9,1:10)         =matmul(T_P0(9*ii+1:9*ii+9,1:9),G(9*ii-8:9*ii,1:10))         +matmul(T(9*ii+1:9*ii+9,1:9),G_P0(9*ii-8:9*ii,1:10));
                else
                    G_P0(9*ii+1:9*ii+9,3*jj+5:3*jj+7)=matmul(T_P0(9*ii+1:9*ii+9,1:9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P0(9*ii-8:9*ii,3*jj+5:3*jj+7));                  
                endif
            else
                if (jj==1) then
                    GG_P0(1:9,1:10)                  =matmul(T_P0(9*ii+1:9*ii+9,1:9),G(9*ii-8:9*ii,1:10))         +matmul(T(9*ii+1:9*ii+9,1:9),G_P0(9*ii-8:9*ii,1:10));                   
                else
                    GG_P0(1:9,3*jj+5:3*jj+7)         =matmul(T_P0(9*ii+1:9*ii+9,1:9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P0(9*ii-8:9*ii,3*jj+5:3*jj+7));                 
                endif
            endif             
            do kk=0,ele_num               
                if(ii<ele_num) then
                    if (jj==1) then
                        G_P1(9*ii+1:9*ii+9,1:10,kk+1)=         matmul(T_P1(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))+         matmul(T(9*ii+1:9*ii+9,1:9),G_P1(9*ii-8:9*ii,1:10,kk+1));
                        G_P2(9*ii+1:9*ii+9,1:10,kk+1)=         matmul(T_P2(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))+         matmul(T(9*ii+1:9*ii+9,1:9),G_P2(9*ii-8:9*ii,1:10,kk+1));
                        G_P3(9*ii+1:9*ii+9,1:10,kk+1)=         matmul(T_P3(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))+         matmul(T(9*ii+1:9*ii+9,1:9),G_P3(9*ii-8:9*ii,1:10,kk+1));                   
                    else
                        G_P1(9*ii+1:9*ii+9,3*jj+5:3*jj+7,kk+1)=matmul(T_P1(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P1(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));
                        G_P2(9*ii+1:9*ii+9,3*jj+5:3*jj+7,kk+1)=matmul(T_P2(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P2(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));
                        G_P3(9*ii+1:9*ii+9,3*jj+5:3*jj+7,kk+1)=matmul(T_P3(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P3(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));                    
                    endif
                else
                    if (jj==1) then
                        GG_P1(1:9,1:10,kk+1)=         matmul(T_P1(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))         +matmul(T(9*ii+1:9*ii+9,1:9),G_P1(9*ii-8:9*ii,1:10,kk+1));
                        GG_P2(1:9,1:10,kk+1)=         matmul(T_P2(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))         +matmul(T(9*ii+1:9*ii+9,1:9),G_P2(9*ii-8:9*ii,1:10,kk+1));
                        GG_P3(1:9,1:10,kk+1)=         matmul(T_P3(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,1:10))         +matmul(T(9*ii+1:9*ii+9,1:9),G_P3(9*ii-8:9*ii,1:10,kk+1));                    
                    else
                        GG_P1(1:9,3*jj+5:3*jj+7,kk+1)=matmul(T_P1(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P1(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));
                        GG_P2(1:9,3*jj+5:3*jj+7,kk+1)=matmul(T_P2(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P2(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));
                        GG_P3(1:9,3*jj+5:3*jj+7,kk+1)=matmul(T_P3(9*ii+1:9*ii+9,9*kk+1:9*kk+9),G(9*ii-8:9*ii,3*jj+5:3*jj+7))+matmul(T(9*ii+1:9*ii+9,1:9),G_P3(9*ii-8:9*ii,3*jj+5:3*jj+7,kk+1));                    
                    endif
                endif 
            enddo
        enddo
    enddo
    !!!--------------------传递函数G的偏导数--------------------!
    Z1_J(:,4)=matmul(transpose(G_P0(:,:)),F)
    do kk=0,ele_num
        Z1_J(:,3*kk+5)=matmul(transpose(G_P1(:,:,kk+1)),F)
        Z1_J(:,3*kk+6)=matmul(transpose(G_P2(:,:,kk+1)),F)
        Z1_J(:,3*kk+7)=matmul(transpose(G_P3(:,:,kk+1)),F)
    enddo
    !
    !call exportf(G,'G3.txt')
    !call exportf(G_P1(:,:,1),'G_P11.txt')
    !call exportf(G_P1(:,:,2),'G_P12.txt')
    !call exportf(G_P1(:,:,3),'G_P13.txt')
    !call exportf(G_P1(:,:,4),'G_P14.txt')
    !call exportf(G_P1(:,:,5),'G_P15.txt')
    !call exportf(G_P1(:,:,6),'G_P16.txt')
    !call exportf(GG_P1(:,:,6),'GG_P11.txt')
    !call exportf(GG_P1(:,:,1),'GG_P11.txt')
    !call exportf(GG_P1(:,:,2),'GG_P12.txt')
    !call exportf(GG_P1(:,:,3),'GG_P13.txt')
    !call exportf(GG_P1(:,:,4),'GG_P14.txt')
    !call exportf(GG_P1(:,:,5),'GG_P15.txt')
    !call exportf(GG_P1(:,:,6),'GG_P16.txt')
    !call exportf(T_P1,'T_P1.txt')
    !call exportf(T_P2,'T_P2.txt')
    !call exportf(T_P3,'T_P3.txt')
    !call exportf(A1_P1,'A1_P1.txt')
    !call exportf(G_P2,'G_P1.txt')
    !call exportf(G_P1,'G_P1.txt')
    !call exportf(Z1,'Z1.txt')
    !call exportf(T_PP1,'T_PP1.txt')
    !call exportf(T_PP2,'T_PP2.txt')
    !call exportf(T_PP3,'T_PP3.txt')
    !
    !call exportf(T_P0,'T_P0.txt')
    !call exportf(T_P1,'T_P1.txt')
    !call exportf(T_P2,'T_P2.txt')
    !call exportf(T_P3,'T_P3.txt')
    
    !Z=matmul(transpose(G),matmul(M,G));
    Z1_J=Z1_J+matmul(transpose(G),F_P);    
    Z1=matmul(transpose(G),F);
    
    Phi_Quater=0d0;
    Phi_Quater(4)=2d0*y(4)
    Phi_Quater(5)=2d0*y(5)
    Phi_Quater(6)=2d0*y(6)
    Phi_Quater(7)=2d0*y(7)
    !print *,'TEST'
    Phi_C_P=0d0;
    do ii=1,con_num
        jj=Con_dof(ii,1);
        !print *,II,JJ
        if (Con_dof(ii,2)==0.or.Con_dof(ii,2)==1) then
                if (jj<=9*ele_num) then 
                    Phi_C(ii,:)=G(jj,:);
                    Phi_C_P(ii,:,5:3*ele_num+7:3)=G_P1(jj,:,:);
                    Phi_C_P(ii,:,6:3*ele_num+7:3)=G_P2(jj,:,:);
                    Phi_C_P(ii,:,7:3*ele_num+7:3)=G_P3(jj,:,:);
                    Phi_C_P(ii,:,4)              =G_P0(jj,:)
                else!if (jj<=9*ele_num+9) then 
                    Phi_C(ii,:)=GG(jj-9*ele_num,:);
                    Phi_C_P(ii,:,5:3*ele_num+7:3)=GG_P1(jj-9*ele_num,:,:);
                    Phi_C_P(ii,:,6:3*ele_num+7:3)=GG_P2(jj-9*ele_num,:,:);
                    Phi_C_P(ii,:,7:3*ele_num+7:3)=GG_P3(jj-9*ele_num,:,:);
                    Phi_C_P(ii,:,4)              =GG_P0(jj-9*ele_num,:)
                endif
                !Phi_C_P(ii,:,:)=Phi_C_P(ii,:,:)
                !call exportf(Phi_C_P(ii,:,:),'PhiCP1')
                !call exportf(Phi_C(ii,:),'PhiC')
                !call exportf(G_P1(jj,:,:),'GG_P1')
                !call exportf(G_P2(jj,:,:),'GG_P2')
                !call exportf(G_P0(:,:),'GG_P0')
                !call exportf(T,'T')
                !call exportf(T_P0,'T_P0')
                !call exportf(T_P1,'T_P1')
                !call exportf(T_P2,'T_P2')
                !call exportf(T_P3,'T_P3')
                !call exportf(G_P1(:,:,1),'G_P11')
                !call exportf(G_P2(:,:,1),'G_P21')
                !call exportf(G_P3(:,:,1),'G_P31')
                !call exportf(G_P1(:,:,2),'G_P12')
                !call exportf(G_P2(:,:,2),'G_P22')
                !call exportf(G_P3(:,:,2),'G_P32')
                !call exportf(G,'G')
        endif
        if (Con_dof(ii,2)==2) then
            kk=con_dof(ii,3);
            if (jj<=9*ele_num) then 
                Phi_C(ii,:)=2d0*r(3*kk-2)*G(jj,:)+2d0*r(3*kk-1)*G(jj+1,:);
                
                Phi_C_P(ii,:,5:3*ele_num+7:3)=2d0*r(3*kk-2)*G_P1(jj,:,:)+2d0*r(3*kk-1)*G_P1(jj+1,:,:);
                Phi_C_P(ii,:,6:3*ele_num+7:3)=2d0*r(3*kk-2)*G_P2(jj,:,:)+2d0*r(3*kk-1)*G_P2(jj+1,:,:);
                Phi_C_P(ii,:,7:3*ele_num+7:3)=2d0*r(3*kk-2)*G_P3(jj,:,:)+2d0*r(3*kk-1)*G_P3(jj+1,:,:);
                Phi_C_P(ii,:,4)              =2d0*r(3*kk-2)*G_P0(jj,:)  +2d0*r(3*kk-1)*G_P0(jj+1,:);
                
                Phi_C_P(ii,:,:)=2d0*(G(jj,:).mt.G(jj,:))+2d0*(G(jj+1,:).mt.G(jj+1,:))+Phi_C_P(ii,:,:)
            
            else
                Phi_C(ii,:)=2d0*r(3*kk-2)*GG(jj-9*ele_num,:)+2d0*r(3*kk-1)*GG(jj+1-9*ele_num,:);
                
                Phi_C_P(ii,:,5:3*ele_num+7:3)=2d0*r(3*kk-2)*GG_P1(jj-9*ele_num,:,:)+2d0*r(3*kk-1)*GG_P1(jj+1-9*ele_num,:,:);
                Phi_C_P(ii,:,6:3*ele_num+7:3)=2d0*r(3*kk-2)*GG_P2(jj-9*ele_num,:,:)+2d0*r(3*kk-1)*GG_P2(jj+1-9*ele_num,:,:);
                Phi_C_P(ii,:,7:3*ele_num+7:3)=2d0*r(3*kk-2)*GG_P3(jj-9*ele_num,:,:)+2d0*r(3*kk-1)*GG_P3(jj+1-9*ele_num,:,:);
                Phi_C_P(ii,:,4)              =2d0*r(3*kk-2)*GG_P0(jj-9*ele_num,:)  +2d0*r(3*kk-1)*GG_P0(jj+1-9*ele_num,:);
                
                Phi_C_P(ii,:,:)=Phi_C_P(ii,:,:)+2d0*(GG(jj-9*ele_num,:).mt.GG(jj-9*ele_num,:))+2d0*(GG(jj+1-9*ele_num,:).mt.GG(jj+1-9*ele_num,:))
            endif
        endif

         if (Con_dof(ii,2)==103.or.Con_dof(ii,2)==104.or.Con_dof(ii,2)==105) then
            LL=Con_dof(ii,3)
            if(Con_dof(ii,4)==1) then
                Array(:,ii)=AI(:,3*LL-2)
            elseif(Con_dof(ii,4)==2) then
                Array(:,ii)=AI(:,3*LL-1)
            elseif(Con_dof(ii,4)==3) then
                Array(:,ii)=AI(:,3*LL)
            endif
            
            call tensor(temp_T,Array(:,ii))            
            if (jj<=9*ele_num) then 
                Array_P(:,:,ii)=temp_T.mt.G(jj:jj+2,:)
                Array_G(:,:,ii)=G(jj:jj+2,:)
                Array_PP(:,:,4,ii)=temp_T.mt.G_P0(jj:jj+2,:)
                Array_PP(:,:,5:3*ele_num+7:3,ii)=temp_T.mt.G_P1(jj:jj+2,:,:)
                Array_PP(:,:,6:3*ele_num+7:3,ii)=temp_T.mt.G_P2(jj:jj+2,:,:)
                Array_PP(:,:,7:3*ele_num+7:3,ii)=temp_T.mt.G_P3(jj:jj+2,:,:)
            else
                Array_P(:,:,ii)=temp_T.mt.GG(4:6,:)
                Array_G(:,:,ii)=GG(4:6,:)
                Array_PP(:,:,4,ii)=(temp_T.mt.GG_P0(4:6,:))
                Array_PP(:,:,5:3*ele_num+7:3,ii)=(temp_T.mt.GG_P1(4:6,:,:))
                Array_PP(:,:,6:3*ele_num+7:3,ii)=(temp_T.mt.GG_P2(4:6,:,:))
                Array_PP(:,:,7:3*ele_num+7:3,ii)=(temp_T.mt.GG_P3(4:6,:,:))
            endif
         endif
         
         if (Con_dof(ii,2)==3) then
              if (jj<=9*ele_num) then 
                 Phi_C(ii,:)=2d0*G(jj,:)
                 Phi_C3(ii,:)=2d0*G(jj+1,:); 
                 Array_P(1:2,:,ii)=G(jj:jj+2,:)
                 Array_PP(1:2,:,4,ii)=G_P0(jj:jj+1,:)
                 Array_PP(1:2,:,5:3*ele_num+7:3,ii)=G_P1(jj:jj+1,:,:)
                 Array_PP(1:2,:,6:3*ele_num+7:3,ii)=G_P2(jj:jj+1,:,:)
                 Array_PP(1:2,:,7:3*ele_num+7:3,ii)=G_P3(jj:jj+1,:,:)
              else
                 Phi_C(ii,:)=2d0*GG(1,:)
                 Phi_C3(ii,:)=2d0*GG(2,:);   
                 Array_P(1:2,:,ii)=GG(1:2,:)
                 Array_PP(1:2,:,4,ii)=GG_P0(1:2,:)
                 Array_PP(1:2,:,5:3*ele_num+7:3,ii)=GG_P1(1:2,:,:)
                 Array_PP(1:2,:,6:3*ele_num+7:3,ii)=GG_P2(1:2,:,:)
                 Array_PP(1:2,:,7:3*ele_num+7:3,ii)=GG_P3(1:2,:,:) 
            endif
         endif
    enddo
    !call exportf(Z1_J,'Z_J.txt')
    !!call exportf(M,'M3.txt')
    !!call exportf(F,'F3.txt')
    !!call exportf(Me,'Me3.txt')
    !call exportf(G,'G3.txt')
    !!call exportf(Z1_J,'Z1_J.txt')
    !!call exportf(G,'G.txt')
    !call exportf(T,'T.txt')
    !call exportf(Z1,'Z13.txt')
    !!call exportf(T,'T3.txt')
    !call exportf(A1,'A1.txt')
    !call exportf(Array_PP(:,:,4,1),'Array_PP1.txt')
    !call exportf(Array_PP(:,:,5,1),'Array_PP2.txt')
    !call exportf(Array_PP(:,:,6,1),'Array_PP3.txt')
    !call exportf(Array_PP(:,:,7,1),'Array_PP4.txt')
    !call exportf(Array_PP(:,:,8,1),'Array_PP5.txt')
    !call exportf(Array_PP(:,:,9,1),'Array_PP6.txt')
    !call exportf(PHi_C,'PHI_C.txt')
    !write(*,*) a
endsubroutine recur_j


subroutine Element_MF_J(Me,Fe,Fe_P,Mp,a,a0,RouA,Le,A1,omega,da,Iz,Cm,Rou,Phi,Ki,Roup)
    !use OPERATORS_MOD
    implicit none
    real(8):: Me(9,9),Me1(9,9),Me2(9,9),Fe(9),Fe_P(9,3),A1(3,3),omega(3),eye(3,3),Cm(3,3),Rou,Mp(6),Ki(3,3),Phi(3,3),Roup(3)
    real(8):: a(3),Le,RouA,da(3),Iz(3,3),a0(3)
    real(8):: Fu(3),Fu_P(3,3),FP(9)
    call force_p(Fp,Mp,A1,Roup,Phi,Ki)
    call mass_force_E_J(Me1,Me2,Fu,Fu_P,a,a0,Le,da,Cm,Iz)
    Fe=0d0
    Fe(7:9)=Fu;
    Fe=-Fe+Fp;
    Fe_P(7:9,:)=-Fu_P;
    
endsubroutine Element_MF_j

subroutine Mass_force_E_J(M1,M2,Fu,Fu_p,a,a0,Le,da,Cm,Iz)
!   force_f(Fw,a,Le,omega,A1,da,Gra)
    implicit none
    real(8)::M1(9,9),M2(9,9),Fu_P(3,3),Fu(3),a(3),a0(3),Le,da(3),Cm(3,3),Iz(3,3);
    real(8)::Fue1(3),Fue1_P(3,3),Me1(9,9),Me2(9,9),g(6),w(6),aa;
    integer::kk
    ! W=[0.236926885056189;0.236926885056189;0.478628670499366;0.478628670499366;0.568888888888889];
    ! g=[0.906179845938664;-0.906179845938664;0.538469310105683;-0.538469310105683;0];
    ! g=[0.861136311594053,0.339981043584856,-0.861136311594053,-0.339981043584856];
     !W=[0.347854845137454,0.652145154862546,0.347854845137454,0.652145154862546];Me1=0;
    g=[0.932469514203152d0,0.661209386466265d0,0.238619186083197d0,-0.932469514203152d0,-0.661209386466265d0,-0.238619186083197d0];
    W=[0.17132449237917d0,0.360761573048139d0,0.467913934572691d0,0.17132449237917d0,0.360761573048139d0,0.467913934572691d0];
    !g=[-0.9602898564975363d0,-0.7966664774136268d0,-0.525532409916329d0,-0.1834346424956498d0,0.1834346424956498d0,0.525532409916329d0,0.7966664774136268d0,0.9602898564975363d0]
    !W=[0.1012285362903768d0,0.2223810344533745d0,0.3137066458778874d0,0.362683783378362d0,0.362683783378362d0,0.3137066458778874d0,0.2223810344533745d0,0.1012285362903768d0]
    !g=[-0.9894009349916499d0,-0.9445750230732326d0,-0.8656312023878318d0,-0.755404408355003d0,-0.6178762444026438d0,-0.4580167776572274d0,-0.2816035507792589d0,-0.09501250983763744d0,0.09501250983763744d0,0.2816035507792589d0,0.4580167776572274d0,0.6178762444026438d0,0.755404408355003d0,0.8656312023878318d0,0.9445750230732326d0,0.9894009349916499d0]
    !W=[0.02715245941175406d0,0.06225352393864778d0,0.0951585116824929d0,0.1246289712555339d0,0.1495959888165768d0,0.1691565193950026d0,0.1826034150449236d0,0.1894506104550685d0,0.1894506104550685d0,0.1826034150449236d0,0.1691565193950026d0,0.1495959888165768d0,0.1246289712555339d0,0.0951585116824929d0,0.06225352393864778d0,0.02715245941175406d0]
    

    Fu=0d0;
    Fu_P=0d0;
    do kk=1,6
            aa=W(kk)*Le/2d0;
            call Mass_force_d_J(Me1,Me2,Fue1,Fue1_P,g(kk)*Le/2d0+Le/2d0,a,a0,Le,da,Cm,Iz)
            Fu=Fu+Fue1*aa;
            Fu_P=Fu_P+Fue1_P*aa;
    enddo

endsubroutine Mass_force_E_j

subroutine Mass_force_d_J(Me1,Me2,Fue1,Fue1_P,kesi,a,a0,Le,da,Cm,Iz)
    use functions
    use OPERATORS_MOD
    implicit none    
    real(8)::Fue1(3),Fue1_P(3,3),Me1(9,9),Me2(9,9),kesi,a(3),a0(3),Le,A2(3,3),da(3),Cm(3,3),Iz(3,3)
    real(8)::droup(3),Roup(3),Phi(3,3),PI(3),Roup_T(3,3)
    real(8)::q(3),dq(3),Kq(3,3),Kqq(3,3),Kqq1(3,3),Kqq2(3,3),Ap(3,3),kapa(3)
    
    !call Phi_Pi_Kesi(droup,Roup,Phi,Pi,a,da,Kesi,Le)

    !call Tensor(Roup_T,Roup)

    q=kesi/Le*a;
    !Kq=K_R(q(1:3))
    !Kqq1=0d0;
    !Kqq1(1,2)= q(3)*cos(q(2))
    !Kqq1(2,1)=-q(2)*sin(q(1))-q(3)*cos(q(1))*cos(q(2))
    !Kqq1(2,2)= q(3)*sin(q(1))*sin(q(2))
    !Kqq1(3,1)= q(2)*cos(q(1))-q(3)*sin(q(1))*cos(q(2))
    !Kqq1(3,2)=-q(3)*cos(q(1))*sin(q(2))
    
    Kq=K_B(q(1:3))
    Kqq1=0d0;
    Kqq1(1,2)=-q(1)*sin(q(2))*cos(q(3))
    Kqq1(1,3)= q(2)*cos(q(3))-q(1)*cos(q(2))*sin(q(3))
    Kqq1(2,2)= q(1)*sin(q(2))*sin(q(3))
    Kqq1(2,3)=-q(2)*sin(q(3))-q(1)*cos(q(2))*cos(q(3))
    Kqq1(3,2)= q(1)*cos(q(2))
    
    Kqq=transpose(Kqq1+Kq)/Le
    Fue1=matmul(matmul(matmul(Kqq,Cm),Kq),a-a0)/Le
    !kapa=Kq.mt.(a-a0)/Le
    Fue1_P=matmul(matmul(Kqq,Cm),Kq)/Le

    !Fue1=Kqq(:,1)*Cm(1,1)*kapa(1)+Kqq(:,2)*Cm(2,2)*kapa(2)+Kqq(:,3)*Cm(3,3)*kapa(3)
    
    !--------------------2015.03.14,对雅可比进行添加，原来下式有误--------------------！
    Fue1_P(:,1)=Fue1_P(:,1)+matmul(matmul(matmul(reshape([0d0,-sin(q(2))*cos(q(3)),-cos(q(2))*sin(q(3)),&
                                                         &0d0, sin(q(2))*sin(q(3)),-cos(q(2))*cos(q(3)),&
                                                         &0d0, cos(q(2))          ,0d0                   ],[3,3]),cm),Kq),a-a0)/Le*kesi/Le**2
    Fue1_P(:,2)=Fue1_P(:,2)+matmul(matmul(matmul(reshape([-sin(q(2))*cos(q(3)),-q(1)*cos(q(2))*cos(q(3)), cos(q(3))+q(1)*sin(q(2))*sin(q(3)),&
                                                         & sin(q(2))*sin(q(3)), q(1)*cos(q(2))*sin(q(3)),-sin(q(3))+q(1)*sin(q(2))*cos(q(3)),&
                                                         & cos(q(2))          ,-q(1)*sin(q(2))          ,0d0                   ],[3,3]),cm),Kq),a-a0)/Le*kesi/Le**2
    
    Fue1_P(:,3)=Fue1_P(:,3)+matmul(matmul(matmul(reshape([-cos(q(2))*sin(q(3)), q(1)*sin(q(2))*sin(q(3))+cos(q(3)),-q(2)*sin(q(3))-q(1)*cos(q(2))*cos(q(3)),&
                                                         &-cos(q(2))*cos(q(3)), q(1)*sin(q(2))*cos(q(3))-sin(q(3)),-q(2)*cos(q(3))+q(1)*cos(q(2))*sin(q(3)),&
                                                         &0d0                 , 0d0                               ,0d0                                      ],[3,3]),cm),Kq),a-a0)/Le*kesi/Le**2
    
    Fue1_P(:,2)=Fue1_P(:,2)+matmul(matmul(Kqq,Cm),[-(a(1)-a0(1))*sin(q(2))*cos(q(3)),(a(1)-a0(1))*sin(q(2))*sin(q(3)),(a(1)-a0(1))*cos(q(2))])/Le*kesi/Le;
                                                     
    Fue1_P(:,3)=Fue1_P(:,3)+matmul(matmul(Kqq,Cm),[-(a(1)-a0(1))*cos(q(2))*sin(q(3))+(a(2)-a0(2))*cos(q(3)),-(a(1)-a0(1))*cos(q(2))*cos(q(3))-(a(2)-a0(2))*sin(q(3)),0d0])/Le*kesi/Le;
    !
endsubroutine Mass_force_d_j



subroutine Phi_P(Phi_P1,Phi_P2,Phi_P3,a,Le)
!----------Phi的偏导数求解，对三个角分别进行求偏导,word雅可比阵20150314------------------------------------------!
        real(8):: phi_P1(3,3),phi_P2(3,3),phi_P3(3,3),a(3),Le
        real(8):: bi1,bi2
        bi1=a(1)+a(2)
        bi2=a(1)-a(2)

        Phi_P1=0d0;
        Phi_P2=0d0;
        Phi_P3=0d0;  
        !--------------------------------------------------------------------------------------------------------!
        Phi_P1(1,2)= 0d0
        
        Phi_P1(2,1)=-Le/2d0*(-(bi1+bi2)/4d0+(bi1**3+bi2**3)/36d0-(bi1**5+bi2**5)/120d0/8d0+(bi1**7+bi2**7)/120d0/56d0/10d0)
        Phi_P1(2,2)=-Le/2d0*(-(bi1-bi2)/4d0+(bi1**3-bi2**3)/36d0-(bi1**5-bi2**5)/120d0/8d0+(bi1**7-bi2**7)/120d0/56d0/10d0)
        
        Phi_P1(3,1)= Le/2d0*(-2d0/3d0+(bi1**2+bi2**2)/10d0-(bi1**4+bi2**4)/168d0+(bi1**6+bi2**6)/120d0/6d0/9d0)
        Phi_P1(3,2)= Le/2d0*(         (bi1**2-bi2**2)/10d0-(bi1**4-bi2**4)/168d0+(bi1**6-bi2**6)/120d0/6d0/9d0)
        !----------------------------对alph的偏导结束------------------------------------------------------------!
        Phi_P2(1,2)= Le*(-a(2)/4d0+a(2)**3/36d0-a(2)**5/120d0/8d0+a(2)**7/120d0/56d0/10d0)
        
        Phi_P2(2,1)=-Le/2d0*(-(bi1-bi2)/4d0+(bi1**3-bi2**3)/36d0-(bi1**5-bi2**5)/120d0/8d0+(bi1**7-bi2**7)/120d0/56d0/10d0)
        Phi_P2(2,2)=-Le/2d0*(-(bi1+bi2)/4d0+(bi1**3+bi2**3)/36d0-(bi1**5+bi2**5)/120d0/8d0+(bi1**7+bi2**7)/120d0/56d0/10d0)
        
        Phi_P2(3,1)= Le/2d0*(         (bi1**2-bi2**2)/10d0-(bi1**4-bi2**4)/168d0+(bi1**6-bi2**6)/120d0/6d0/9d0)
        Phi_P2(3,2)= Le/2d0*(-2d0/3d0+(bi1**2+bi2**2)/10d0-(bi1**4+bi2**4)/168d0+(bi1**6+bi2**6)/120d0/6d0/9d0)
        !----------------------------对beta的偏导结束------------------------------------------------------------!       
                 
              
endsubroutine