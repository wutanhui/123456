SUBROUTINE Initial_DATA
        use body_module
        use Element_module
        use Model_module
        use integration_parameter
        use OPERATORS_MOD
        implicit none
        interface 
            subroutine recur_post(r_temp,dr_temp,ddr_temp,y_temp,dy_temp,ddy_temp,omega_temp,d_omega_temp,Le,ele_num_all)
                real(8),allocatable::r_temp(:),dr_temp(:),ddr_temp(:),y_temp(:),dy_temp(:),ddy_temp(:),omega_temp(:),d_omega_temp(:)
                real(8)::Le
                integer::ele_num_all
            end subroutine 
        end interface
            integer i,K,II,I_begin,I_end


            call Allocate_space;
            Gra=[0d0,0d0,0d0];
            k=0;

            do I=1,body_num
                if (I<=3) then 
                    body(I).F_num=3*body(I).ele_num+7;
                    allocate(body(I).Fn(Body(I).F_num));
                    body(I).FN=[1:(3*body(I).ele_num+7)];
                    !body(I).con_dof=[1,2,3,4,5,6,91,92,93,181,182,183, 271,272,273];
                    body(I).Contact_node_num=I;
                    !allocate(body(I).contact_node(body(I).Contact_node_num));
                    !body(I).contact_node=[11,21,31]
                else
                    body(I).F_num=3*body(I).ele_num+7
                    allocate(body(I).Fn(Body(I).F_num))
                    body(I).FN=[1:(3*body(I).ele_num+7)]
                    !body(I).con_dof=[1,2,3,46,47,48]
                    body(I).Contact_node_num=0;
                endif
                !初始位置

                body(I).a=0.d0;
                !body(I).a(1:3*body(I).ele_num:3)=0.1d0;
                !body(I).a(2:3*body(I).ele_num:3)=0.2d0;
                !body(1).a(3*body(1).ele_num-5)=1d-2;
                !body(2).a(3*body(2).ele_num-5)=2d-2;
                !body(3).a(3*body(3).ele_num-5)=3d-2;
                !a(1:3*ele_num_all/2:3)=-1d0*3.1415926d0/ele_num_all;a(2:3*ele_num_all:3)=0.0d0;a(3:3*ele_num_all:3)=0.0d0;!a(3*ele_num_all-2)=0.01d0;
                !a(1:3*ele_num_all/2:3)=4d0*3.1415926d0/ele_num_all;
                !call importf(a(121:3*ele_num_all),'y0723.txt')
            
                body(I).da=0.d0;
                !body(I).da(2:3*body(I).ele_num:3)=0.d0
                body(I).dda=0.d0;
                
                !a0=[0.00d0,-0.000d0,0.0d0]

                
                !omega1=[8.660254083d0,0d0,-4.999999923d0];
                
                body(I).r=0d0;
                body(I).dr=0d0;
                body(I).ddr=0d0;
                body(I).omega=0d0;
                body(I).domega=0d0;
                body(I).y=0d0;
                body(I).dy=0.0d0;
                body(I).ddy=0d0;
                body(I).Mp=0.d0;
                body(I).y(1:3*body(I).ele_num+7,1)=[body(I).r0, body(I).Quater ,body(I).a]
                body(I).con_dof=-1;
            enddo
            !body(1).Mp(58)=1d0;
            !body(2).Mp(59)=1d0;
            !body(3).Mp(58)=-1d0;
            !!从文件中读取数据
            !open(11,file='x_initial.txt',STATUS='OLD')
            !do I=1,body_num
            !    read(11,*) (body(I).y(K,1),K=1,3*body(I).ele_num+6)
            !enddo
            !close(11)
            
            !body(1).dy(1:3,1)=[0d0,0.1d0,0d0]*1d1;
            !body(1).dy(8,1)=2d0/body(1).Le*1d0
            !body(1).dy(11:(3*body(1).ele_num+6):3,1)=-(-1)**[2:body(1).ele_num]*4d0/body(1).Le*1d0
            !body(1).dy(1:3,1)=[0d0,0.1d0,0d0]*1d1;
            !body(2).dy(1:3,1)=[-0.086602540378444d0,-0.05d0,0d0]*1d1;
            !body(3).dy(1:3,1)=[0.086602540378444d0,-0.05d0,0d0]*1d1;
            !
            !body(1).dy(8,1)=-2d0/0.375d0*1d1
            !body(1).dy(11:(3*body(1).ele_num+6):3,1)=(-1)**[2:body(1).ele_num]*4d0/0.375d0*1d1
            !!
            !body(2).dy(8,1)=-2d0/0.375d0*1d1
            !body(2).dy(11:(3*body(2).ele_num+6):3,1)=(-1)**[2:body(2).ele_num]*4d0/0.375d0*1d1
            !
            !body(3).dy(8,1)=-2d0/0.375d0*1d1
            !body(3).dy(11:(3*body(3).ele_num+6):3,1)=(-1)**[2:body(3).ele_num]*4d0/0.375d0*1d1
            a0=0d0;
            dr1=[0.d0,0.d0,0.d0];
            omega1=[0d0,00.0d0,0.d0];     
            !r1=[0.3d0,0.0d0,0.d0];
            !r2=[-0.15d0,0.26d0,0.d0];
            !r3=[-0.15d0,-0.26d0,0.d0];
            !body(1).y(1:3*body(1).ele_num+6,1)=[r1,0d0,0d0,3.1415926d0/2,body(1).a];
            !body(2).y(1:3*body(2).ele_num+6,1)=[r2,0d0,0d0,3.14159260d0/6*7,body(2).a];
            !body(3).y(1:3*body(3).ele_num+6,1)=[r3,0d0,0d0,3.14159260d0/6*11,body(3).a];
                !!!从文档里读
                !call importf(y(:,1),'y.txt')
                !!!
                !dy(:,1)=[dr1,omega1,da]
            y=0d0;
            dy=0d0;
            ddy=0d0;
            I_begin=1
            !Do i=1,body_num
            !    I_end=I_begin+3*body(I).ele_num+5;
            !    x(I_begin:I_end,1)=[body(I).y(:,1)];
            !    x(I_begin+3*ele_num_all+6*body_num:I_end+3*ele_num_all+6*body_num,1)=body(I).dy(:,1)
            !    I_begin=I_end+1;
            !enddo
            Do i=1,body_num
                I_end=I_begin+3*body(I).ele_num+6;
                y(I_begin:I_end,1)=[body(I).y(:,1)];
                dy(I_begin:I_end,1)=body(I).dy(:,1)
                I_begin=I_end+1;
            enddo
            ddy=0d0;
            !!!!!!!!!!!!!!!全局约束划去点
            F_num_all=sum(body(:).F_num);
            allocate(Fn(sum(body(:).F_num)))
            I_begin=1
            do I=1,body_num
                I_end=I_begin+body(I).f_num-1
                Fn(I_begin:I_end)=body(I).Fn+3*sum(body(1:(I-1)).ele_num)+(I-1)*7
                I_begin=I_end+1        
            enddo
            
            dx(:,1)=0.d0;
            call importf(y(:,1),'y0.txt')
            !dx(1:3*ele_num_all+6*body_num,1)=dy(:,1);
            !P=0.d0;
            !!!!!!!!!!!!!!!!加力
            !body(1).Mp(6*body(1).ele_num-5) =1.693d1/3;
            !body(2).Mp(6*body(2).ele_num-5)=-1.693d1/3;
            !body(3).Mp(6*body(3).ele_num-5)=1.693d1/3;
            !%%%%%%%%%%%%%%%%
            I_begin=1
            do II=1,body_num
                I_end=I_begin+3*body(II).ele_num+6
                body(II).y=y(I_begin:I_end,:)
                body(II).dy=dy(I_begin:I_end,:)
                body(II).ddy=ddy(I_begin:I_end,:)
                I_begin=I_end+1;
            enddo
    
            do ii=1,body_num
                body(II).r_temp=0d0;body(II).dr_temp=0d0;body(II).ddr_temp=0d0;body(II).y_temp=0d0;body(II).dy_temp=0d0;body(II).ddy_temp=0d0;body(II).omega_temp=0d0;body(II).domega_temp=0d0;
                body(II).y_temp=body(II).y(:,1)
                body(II).dy_temp=body(II).dy(:,1)
                body(II).ddy_temp=body(II).ddy(:,1)
                call recur_post(body(II).r_temp,body(II).dr_temp,body(II).ddr_temp,body(II).y_temp,body(II).dy_temp,body(II).ddy_temp,body(II).omega_temp,body(II).domega_temp,body(II).Le,body(II).ele_num)
                body(II).r(:,1)=body(II).r_temp
                body(II).dr(:,1)=body(II).dr_temp
                body(II).omega(:,1)=body(II).omega_temp
            enddo
            !%%%%%%%%%%%%%%%%%%%%%%%%%%
            con_value=0d0;
            II=1
            do I=1,con_num_all
                if (con_list(I,1)==0) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+con_list(I,4);!约束了物体的第几个自由度
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)!约束类型
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)!约束的第几个节点
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)!约束的第几个节点的第几自由度
                    !if(con_list(I,4)<3) then
                        con_value(I)=body(con_list(I,2)).r(3*con_list(I,3)-3+con_list(I,4),1)
                    !else
                        !con_value(I)=body(con_list(I,2)).y(3*con_list(I,3)-3+con_list(I,4),1)
                        !con_value(i)=3.75d0
                    !endif
                elseif (con_list(I,1)==1) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+con_list(I,4);
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                    body(con_list(I,5)).con_dof(con_list(I,9),1)=9*con_list(I,6)-9+con_list(I,7);
                    body(con_list(I,5)).con_dof(con_list(I,9),2)=con_list(I,1)
                    body(con_list(I,5)).con_dof(con_list(I,9),3)=con_list(I,6)
                    body(con_list(I,5)).con_dof(con_list(I,9),4)=con_list(I,7)
                elseif (con_list(I,1)==2) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+con_list(I,4);
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                elseif (con_list(I,1)==3) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+con_list(I,4);
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                    body(con_list(I,5)).con_dof(con_list(I,9),1)=9*con_list(I,6)-9+con_list(I,7);
                    body(con_list(I,5)).con_dof(con_list(I,9),2)=con_list(I,1)
                    body(con_list(I,5)).con_dof(con_list(I,9),3)=con_list(I,6)
                    body(con_list(I,5)).con_dof(con_list(I,9),4)=con_list(I,7)
                elseif (con_list(I,1)==103) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+4;
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                    body(con_list(I,5)).con_dof(con_list(I,9),1)=9*con_list(I,6)-9+4;
                    body(con_list(I,5)).con_dof(con_list(I,9),2)=con_list(I,1)
                    body(con_list(I,5)).con_dof(con_list(I,9),3)=con_list(I,6)
                    body(con_list(I,5)).con_dof(con_list(I,9),4)=con_list(I,7)
                elseif (con_list(I,1)==104) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+4;
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                elseif (con_list(I,1)==105) then
                    body(con_list(I,2)).con_dof(con_list(I,8),1)=9*con_list(I,3)-9+4;
                    body(con_list(I,2)).con_dof(con_list(I,8),2)=con_list(I,1)
                    body(con_list(I,2)).con_dof(con_list(I,8),3)=con_list(I,3)
                    body(con_list(I,2)).con_dof(con_list(I,8),4)=con_list(I,4)
                endif  
                if(con_list(I,10)==1) then
                   Act_Cn(II)=I;
                   II=II+1;
                endif
            enddo
            
       
            !!!!!!!!!!!!!!!!!
            !Mp=-0d0
            !Mp(2)=-0.2930d0;
            !Mp(6*ele_num_all-1)=0.2930d0;
            !integration_perimeter
            Rad=1d0;
            alpha_f=Rad/(Rad+1d0);!0d0!(1d0-Rad)/(Rad+1d0);!
            alpha_m=(2d0*Rad-1d0)/(Rad+1d0);!0d0!
            beta=(1d0-alpha_m+alpha_f)**2d0/4d0;
            gama=0.5d0-alpha_m+alpha_f;
            gama_1=gama/del_t/beta
            beta_1=(1d0-alpha_m)/beta/(1d0-alpha_f)/del_t**2
            !!!约速参数
            
            !if (con_num>0) then
            !    !con_dof_body(ii)=[(9*16+1):(9*19+1):9,(9*16+3):(9*19+3):9]
            !    con_dof_body(ii)(1)=1
            !    con_dof_body(ii)(2)=2
            !    !con_dof_body(ii)(3)=9*ele_num_all-6
            !endif
            
    endsubroutine Initial_Data
    
    subroutine Allocate_space
        use Element_module
        use body_module
        use Model_module
        implicit none
        integer::I

        call import_data

        
        ele_num_all=sum(body(:).ele_num);
        Node_num_all=sum(body(:).node_num);
        
        do I=1,Body_num
            allocate( Body(I).a(3*Body(I).ele_num),Body(I).da(3*Body(I).ele_num),Body(I).dda(3*Body(I).ele_num))
            allocate( Body(I).r(3*Body(I).Node_num,n_end+1),Body(I).dr(3*Body(I).Node_num,n_end+1),Body(I).ddr(3*Body(I).Node_num,n_end+1),Body(I).omega(3*Body(I).Node_num,n_end+1),Body(I).domega(3*Body(I).Node_num,n_end+1))
            allocate( Body(I).y(3*Body(I).Node_num+4,n_end+1),Body(I).dy(3*Body(I).Node_num+4,n_end+1),Body(I).ddy(3*Body(I).Node_num+4,n_end+1))
            allocate( Body(I).y_temp(3*Body(I).Node_num+4),Body(I).dy_temp(3*Body(I).Node_num+4))
            
            allocate( Body(I).Roup_temp(3*Body(I).Ele_num),Body(I).dRoup_temp(3*Body(I).Ele_num),Body(I).AI_temp(3,3*Body(I).Ele_num+3))
            
            allocate( Body(I).Phi_temp(3,9*Body(I).Ele_num),Body(I).Ki_temp(3,3*Body(I).Ele_num),Body(I).omega_r_temp(3*Body(I).Ele_num))
            
            allocate( Body(I).G(9*Body(I).ele_num+9,3*Body(I).ele_num+10),Body(I).g0(9*Body(I).ele_num+9,Body(I).ele_num+1))
            
            allocate( Body(I).Z(3*Body(I).Node_num+4,3*Body(I).Node_num+4),Body(I).Z1(3*Body(I).Node_num+4))
            
            allocate( Body(I).r_temp(3*Body(I).Node_num),Body(I).dr_temp(3*Body(I).Node_num),Body(I).ddr_temp(3*Body(I).Node_num),Body(I).omega_temp(3*Body(I).Node_num),Body(I).domega_temp(3*Body(I).Node_num))

            allocate( Body(I).r_temp1(3*Body(I).Node_num,1),Body(I).dr_temp1(3*Body(I).Node_num,1),Body(I).ddr_temp1(3*Body(I).Node_num,1),Body(I).omega_temp1(3*Body(I).Node_num,1),Body(I).domega_temp1(3*Body(I).Node_num,1),Body(I).y_temp1(3*Body(I).Node_num+4,1),Body(I).dy_temp1(3*Body(I).Node_num+4,1),Body(I).ddy_temp1(3*Body(I).Node_num+4,1))

            allocate( Body(I).Mp(6*Body(I).Node_num))
            allocate( Body(I).Con_DOF(Body(I).con_num,4),Body(I).gama_w(Body(I).con_num))
            allocate( Body(I).Phi_C(Body(I).Con_num,3*Body(I).Ele_num+7),Body(I).Phi_C3(Body(I).Con_num,3*Body(I).Ele_num+7))
            allocate( Body(I).GG(9,3*Body(I).Ele_num+7),Body(I).F_contact(9*Body(I).Ele_num+9) )
            
            allocate( Body(I).Z1_J(3*Body(I).Ele_num+7,3*Body(I).Ele_num+7))
            
            allocate( Body(I).Ai(3,3*Body(I).Ele_num+3))
            
            allocate( Body(I).Array(3,Body(I).con_num),Body(I).Array_P(3,3*Body(I).ele_num+7,Body(I).con_num),Body(I).Array_G(3,3*Body(I).ele_num+7,Body(I).con_num))

            allocate( Body(I).Array_PP(3,3*Body(I).ele_num+7,3*Body(I).ele_num+7,Body(I).con_num))
            
            !allocate( Body(I).Phi_C_P1(Body(I).Con_num,3*Body(I).Ele_num+7,Body(I).Ele_num+1),Body(I).Phi_C_P2(Body(I).Con_num,3*Body(I).Ele_num+7,Body(I).Ele_num+1),Body(I).Phi_C_P3(Body(I).Con_num,3*Body(I).Ele_num+7,Body(I).Ele_num+1))
            allocate( Body(I).Phi_C_P(Body(I).Con_num,3*Body(I).Ele_num+7,3*Body(I).Ele_num+7))    
            allocate( BodY(I).Phi_quater(3*Body(I).ele_num+7))
            !allocate( Body(I).Phi_C_P0(Body(I).Con_num,3*Body(I).Ele_num+7))

        enddo
        
        !ALLOcate( X(6*ele_num_all+12*Body_num,n_end+1),dx(6*ele_num_all+12*Body_num,n_end+1))
        ALLOcate( y(3*ele_num_all+7*Body_num,n_end+1),dy(3*ele_num_all+7*Body_num,n_end+1),ddy(3*ele_num_all+7*Body_num,n_end+1),lumda(con_num_all,n_end+1))
        allocate( P(9*ele_num_all),Mp(6*ele_num_all))
        allocate( lumda_quater(body_num,N_end+1))
    endsubroutine allocate_space         
    
    
    subroutine import_data
        use Element_module
        use body_module
        use Model_module
        use OPERATORS_MOD
        use integration_parameter
        implicit none
        integer I,J,K;
        
        open(11,file='import_data_1_5.dat',STATUS='OLD')
        read(11,*) Body_num,con_num_all,Del_t0,N_end
        
        allocate( body(body_num))
        do I=1,Body_num
            read(11,*) body(I).ele_num,body(I).Con_num,body(I).Node_num
            

            read(11,*) (Body(I).r0(K),k=1,3)
            read(11,*) (Body(I).Quater(K),k=1,4)
            
            read(11,*) body(I).length,body(I).E,body(I).Rou,body(I).S
            body(I).RouA=body(I).Rou*body(I).S
            do J=1,3
                read(11,*) (body(I).Iz(J,K),k=1,3)
            enddo
            body(I).Le=body(I).length/body(I).ele_num
            body(I).Cm=0d0;
            body(I).Cm(1,1)=body(I).Iz(1,1)*body(I).E
            body(I).Cm(2,2)=body(I).Iz(2,2)*body(I).E
            body(I).Cm(3,3)=body(I).Iz(3,3)*body(I).E/2.6d0
        enddo
        
        close(11)
        
        allocate(con_list(con_num_all,10),Con_Value(con_num_all))
        if(con_num_all>0) then
               call importf(con_list(:,:),'con_list_1124.txt')
        endif
        Con_num_Act=sum(con_list(:,10))
        allocate(Act_Cn(con_num_act))
        !Act_Cn=[1:107,345:442,680:777,1015:1038]
        
    endsubroutine import_data
    
    
    
