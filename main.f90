!  recur.f90 
!
!  FUNCTIONS:
!  recur - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: recur
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program main
        use body_module
        use Element_module
        use Model_module
        use CPU_timer
        use OPERATORS_MOD
        use BLAS95
        use LAPACK95
        use integration_parameter
        !use functions
        !use lin_sol_gen_int
    implicit none
    
    integer ii,kkk,jj,kk,IIB,I_begin,I_end,I,ijk,ijk1,ll,mm,nn,oo,A11,A22,A33
    real(8)::pi
    real(8),allocatable:: a_temp(:),y_temp(:),dy_temp(:),ddy_temp(:),y_temp1(:),dy_temp1(:),ddy_temp1(:),r_temp(:),dr_temp(:),omega_temp(:),Equ(:),Equ1(:),Equ2(:),Equ3(:),Jac(:,:),lumda_temp(:),lumda_temp2(:),lumda_quater_temp(:);
    real(8),allocatable:: x_temp(:),dx_temp(:),k1(:),k2(:),k3(:),k4(:)
    real(8),allocatable:: equ_temp(:),Jac_temp(:,:)
    integer,allocatable:: ipiv(:)
    character(10):: r_name,dr_name,omega_name
    type(secnds_timer):: tt
    INTERFACE
        SUBROUTINE Acc0(y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,Le,ele_num_all,RouA,Iz,Gra,a0,P,Mp,Cm,Rou,Fn,F_num)  
            real(8),allocatable:: y_temp(:),dy_temp(:),ddy_temp(:),r_temp(:),dr_temp(:),omega_temp(:),P(:),Mp(:)
            real(8)::Le,RouA,Iz(3,3),Gra(3),a0(3),Cm(3,3),Rou
            integer::ele_num_all,F_num,Fn(:)
        END SUBROUTINE   
    END INTERFACE
     INTERFACE
        SUBROUTINE Jacobi(Time,Time_I,y_temp,dy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,a0,lumda_temp,Act_Cn,lumda_temp2,lumda_Quater_temp,Jac)  
            real(8),allocatable:: y_temp(:),dy_temp(:),ddy_temp(:),r_temp(:),dr_temp(:),omega_temp(:),lumda_temp(:),lumda_temp2(:),lumda_Quater_temp(:),Jac(:,:)
            real(8)::time,Le,RouA,Iz(3,3),Cm(3,3),Rou,a0(3)
            integer::Time_I,body_num,Fn(:),F_num_all,ele_num_all,Con_num_Act,Con_num_all,con_list(:,:),Act_Cn(:)
        END SUBROUTINE   
     END INTERFACE
     INTERFACE
        SUBROUTINE equation(Time,Time_I,y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num,ele_num_all,Con_num_all,Con_num_Act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,equ)  
        !,Le_body,ele_num_all,ele_num_body,body_num,RouA_body,Iz_body,Cm_body,Rou_body,Fn,F_num_body,Con_num_body,con_dof_body,lumda_temp2(:)
            real(8),allocatable:: y_temp(:),dy_temp(:),ddy_temp(:),r_temp(:),dr_temp(:),omega_temp(:),con_value(:),P(:),Mp(:),lumda_temp(:),equ(:),lumda_Quater_temp(:)
            !real(8),allocatable:: Le_body(:),RouA_body(:),Iz_body(:,:,:),Cm_body(:,:,:),Rou_body(:)
            real(8)::Time,Gra(3),a0(3)
            !type(Body_type),allocatable:: Body(:)
            integer::Time_I,ele_num_all,F_num,body_num,Fn(:),Con_list(:,:),Con_num_all,Con_num_Act,Act_Cn(:)!,F_num,Con_num_Body(:),con_dof_body(:,:),ele_num_body(:)
        END SUBROUTINE   
     END INTERFACE
    interface 
        subroutine recur_post(r_temp,dr_temp,ddr_temp,y_temp,dy_temp,ddy_temp,omega_temp,d_omega_temp,Le,ele_num_all)
            real(8),allocatable::r_temp(:),dr_temp(:),ddr_temp(:),y_temp(:),dy_temp(:),ddy_temp(:),omega_temp(:),d_omega_temp(:)
            real(8)::Le
            integer::ele_num_all
        end subroutine 
    end interface
    print *,'star'
    call initial_DATA;
    call tt.tic;
   
    allocate(k1(2*(3*ele_num_all+7*body_num)),k2(2*(3*ele_num_all+7*body_num)),k3(2*(3*ele_num_all+7*body_num)),k4(2*(3*ele_num_all+7*body_num)),time1(N_end+1))
    !allocate(x_temp(6*ele_num_all+12*body_num),dx_temp(6*ele_num_all+12*body_num),r_temp(3*ele_num_all+3*body_num),dr_temp(3*ele_num_all+3*body_num),omega_temp(3*ele_num_all+3*body_num),lumda(con_num_all,n_end+1),lumda_temp2(9*ele_num_all+9))
    allocate(dy_temp(3*ele_num_all+7*body_num),ddy_temp(3*ele_num_all+7*body_num),y_temp(3*ele_num_all+7*body_num),a_temp(3*ele_num_all+6*body_num),r_temp(3*ele_num_all+3*body_num),dr_temp(3*ele_num_all+3*body_num),omega_temp(3*ele_num_all+3*body_num))
    ALLOCATE(ipiv(F_num_ALL+body_num+con_num_act),Jac(F_num_all+body_num+con_num_act,F_num_all+body_num+con_num_act))
    
    allocate(lumda_quater_temp(body_num))
    
    allocate(Equ(F_num_all+body_num+con_num_act),Equ1(F_num_all+body_num+con_num_act),Equ2(3*ele_num_all+6*body_num+con_num_act),Equ3(3*ele_num_all+6*body_num+con_num_act))
    
    allocate(Jac_temp(F_num_all+body_num+con_num_act,F_num_all+body_num+con_num_act),Equ_temp(F_num_all+con_num_act))
    if (con_num_all>0) then
       allocate(lumda_temp(con_num_all))
    endif
    !dx_temp=0d0;r_temp=0d0;dr_temp=0d0;omega_temp=0d0;x_temp=0d0;lumda_temp=0d0;lumda_temp2=0d0    
    y_temp=0d0;dy_temp=0d0;ddy_temp=0d0;r_temp=0d0;dr_temp=0d0;omega_temp=0d0;
    ipiv = 0
    !FN=[1:3*ele_num+6];
    Jac=0d0;
    a_temp=0d0;
    !y=0d0;
    y_temp=y(:,1);
    dy_temp=dy(:,1)
    ddy_temp=ddy(:,1)
    Equ=0d0
    equ1=0d0;equ2=0d0
    lumda_quater=0d0
    lumda=0d0;
    time1=0d0;
    ijk=0
    ijk1=0
    Pi=1d0*3.1415926d0
    nn=0
    oo=1
    !do jj=1,Mn
    !y(1,1)=0.3d0*cos(1d1/3d0*del_t0/100d0)
    !y(2,1)=0.3d0*sin(1d1/3d0*del_t0/100d0)     
    !y(1+307,1)=0.3d0*cos(1d1/3d0*del_t0/100d0+2.094395102393195d0)
    !y(2+307,1)=0.3d0*sin(1d1/3d0*del_t0/100d0+2.094395102393195d0)
    !        
    !y(1+614,1)=0.3d0*cos(1d1/3d0*del_t0/100d0+4.188790204786391d0)
    !y(2+614,1)=0.3d0*sin(1d1/3d0*del_t0/100d0+4.188790204786391d0)
    
    !y(9,1)=0.0001d0
    !y(316,1)=0.0001d0
    !y(623,1)=0.0001d0
    !time1(1)=1.92d0;
    !lumda_quater(1,1)=-5.018790710000000d-09
    !call importf(lumda(:,1),'lumda0.txt');
    !call importf(y(:,1),'y0.txt');

        
        do iib=104,con_num_all
            if(body(1).r(3*con_list(iib,3),1)<1.2d0.and.body(1).r(3*con_list(iib,3),1)>0.95d0) then                
                con_list(iib,10)=1
            endif
            if(body(1).r(3*con_list(iib,3),1)>1.2d0.or.body(1).r(3*con_list(iib,3),1)<0.95d0) then                
                con_list(iib,10)=0
            endif
        enddo
        con_num_act=sum(con_list(:,10))
        deallocate(Jac)                          
        deallocate(Equ,Equ1,Equ2,Equ3)             
        deallocate(Jac_temp,Equ_temp,Act_cn)                    
        ALLOCATE(Jac(F_num_all+body_num+con_num_Act,F_num_all+body_num+con_num_Act))                    
        allocate(Equ(F_num_all+body_num+con_num_Act),Equ1(F_num_all+body_num+con_num_Act),Equ2(F_num_all+body_num+con_num_Act),Equ3(F_num_all+body_num+con_num_Act))              
        allocate(Jac_temp(F_num_all+con_num_Act+body_num,F_num_all+con_num_Act+body_num),Equ_temp(F_num_all+con_num_Act+body_num),Act_cn(con_num_Act))
        
        mm=1
        do ll=1,con_num_all
            if(con_list(ll,10)==1) then
                Act_Cn(mm)=ll;
                mm=mm+1;
            endif
        enddo
        
        
        con_value=0d0
        Mn=4000
        Del_t=Del_t0
        body(1).Mp(6*body(1).ele_num-3)=1d0;
    do jj=1,Mn
        do ii=1,n_end
           
            time1(ii+1)=time1(ii)+del_t;
            !Mp=0d0!
            !if(time1(ii)<4)then
            !    body(1).Mp(6*body(1).ele_num-1)=0.293*(1d0-COS(3.1415926/0.2*time1))/2;
            !    body(1).Mp(6*body(1).ele_num-1)=293d0/0.2d0*time1(ii)/12d0;
            !else
                !
            !endif
            y_temp   =  y(:,ii)
            dy_temp  = dy(:,ii)
            lumda_temp(Act_cn)= 0d0;
            lumda_quater_temp=0d0
            !con_value(1)=0.3d0*1d2/9d0*cos(1d1/3d0*time1(ii+1))
            if(time1(ii)<0.3) then
                con_value(2)=-2d0
            else
                con_value(2)=0d0
            endif
            !con_value(4)=0.3d0*1d2/9d0*cos(1d1/3d0*time1(ii+1)+2.094395102393195d0)
            !con_value(5)=0.3d0*1d2/9d0*sin(1d1/3d0*time1(ii+1)+2.094395102393195d0)
            !con_value(7)=0.3d0*1d2/9d0*cos(1d1/3d0*time1(ii+1)+4.188790204786391d0)
            !con_value(8)=0.3d0*1d2/9d0*sin(1d1/3d0*time1(ii+1)+4.188790204786391d0)
            call equation(Time1(ii+1),II+1,y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ1)
            ddy(:,ii)=ddy_temp;
            lumda(Act_cn,ii)=lumda_temp(Act_cn);
            lumda_quater(:,ii)=lumda_quater_temp
            do kk=1,body_num
                body(kk).r(:,ii)= body(kk).r_temp;
                body(kk).dr(:,ii)=body(kk).dr_temp;
                body(kk).omega(:,ii)=body(kk).omega_temp;
            enddo
                        
            k1=[del_t*dy_temp,del_t*ddy_temp];
            y_temp=y(:,ii)+del_t*dy_temp/2d0;
            dy_temp=dy(:,ii)+del_t*ddy_temp/2d0;
            
            call equation(Time1(ii+1),II+1,y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ1)
            k2=[del_t*dy_temp,del_t*ddy_temp];
            y_temp=y(:,ii)+del_t*dy_temp/2d0;
            dy_temp=dy(:,ii)+del_t*ddy_temp/2d0;
            
            call equation(Time1(ii+1),II+1,y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ1)
            k3=[del_t*dy_temp,del_t*ddy_temp];
            y_temp=y(:,ii)+del_t*dy_temp;
            dy_temp=dy(:,ii)+del_t*ddy_temp;
            
            call equation(Time1(ii+1),II+1,y_temp,dy_temp,ddy_temp,r_temp,dr_temp,omega_temp,body_num,Fn,F_num_all,ele_num_all,Con_num_all,Con_num_Act,Con_list,con_value,Gra,a0,lumda_temp,Act_Cn,lumda_Quater_temp,Equ1)           
            k4=[del_t*dy_temp,del_t*ddy_temp];
            
            y(:,ii+1)=y(:,ii)+(k1(1:3*ele_num_all+7*body_num)+2.d0*k2(1:3*ele_num_all+7*body_num)+2.d0*k3(1:3*ele_num_all+7*body_num)+k4(1:3*ele_num_all+7*body_num))/6.d0
            
            dy(:,ii+1)=dy(:,ii)+(k1(3*ele_num_all+7*body_num+1:6*ele_num_all+14*body_num)+2.d0*k2(3*ele_num_all+7*body_num+1:6*ele_num_all+14*body_num)+2.d0*k3(3*ele_num_all+7*body_num+1:6*ele_num_all+14*body_num)+k4(3*ele_num_all+7*body_num+1:6*ele_num_all+14*body_num))/6.d0
            
            If(mod(ii,50)==1) then
                I_begin=1
                do IIB=1,body_num
                    I_end=I_begin+3*body(IIB).ele_num+6
                    body(IIB).y_temp=y_temp(I_begin:I_end)
                    body(IIB).dy_temp=dy_temp(I_begin:I_end)
                    body(IIB).ddy_temp=ddy_temp(I_begin:I_end)
                    I_begin=I_end+1;
                enddo
        
                do iiB=1,body_num
                    body(IIB).r_temp=0d0;body(IIB).dr_temp=0d0;body(IIB).ddr_temp=0d0;body(IIB).omega_temp=0d0;body(IIB).domega_temp=0d0;
                    call recur_post(body(IIB).r_temp,body(IIB).dr_temp,body(IIB).ddr_temp,body(IIB).y_temp,body(IIB).dy_temp,body(IIB).ddy_temp,body(IIB).omega_temp,body(IIB).domega_temp,body(IIB).Le,body(IIB).ele_num)
                    write(omega_name,'(I3)') IIB
                    body(IIB).r_temp1(:,1)=body(IIB).r_temp;
                    body(IIB).dr_temp1(:,1)=body(IIB).dr_temp;
                    body(IIB).ddr_temp1(:,1)=body(IIB).ddr_temp;
                    body(IIB).omega_temp1(:,1)=body(IIB).omega_temp;
                    body(IIB).domega_temp1(:,1)=body(IIB).domega_temp;
                    body(IIB).y_temp1(:,1)=body(IIB).y_temp;
                    body(IIB).dy_temp1(:,1)=body(IIB).dy_temp
                    body(IIB).ddy_temp1(:,1)=body(IIB).ddy_temp
                    if (ii==1.and.jj==1) then
                        call exportf(transpose(body(IIB).r_temp1),'r'//trim(adjustl(omega_name))//'.txt')
                        call exportf(transpose(body(IIB).dr_temp1),'dr'//trim(adjustl(omega_name))//'.txt')
                        call exportf(transpose(body(IIB).ddr_temp1),'ddr'//trim(adjustl(omega_name))//'.txt')
                        call exportf(transpose(body(IIB).omega_temp1),'omega'//trim(adjustl(omega_name))//'.txt')
                        call exportf(transpose(body(IIB).y_temp1),'y'//trim(adjustl(omega_name))//'.txt')    
                        call exportf(transpose(body(IIB).dy_temp1),'dy'//trim(adjustl(omega_name))//'.txt')
                        call exportf(transpose(body(IIB).ddy_temp1),'ddy'//trim(adjustl(omega_name))//'.txt')
                    else
                        call exportf(transpose(body(IIB).r_temp1),'r'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).dr_temp1),'dr'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).ddr_temp1),'ddr'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).omega_temp1),'omega'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).y_temp1),'y'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).dy_temp1),'dy'//trim(adjustl(omega_name))//'.txt',1)
                        call exportf(transpose(body(IIB).ddy_temp1),'ddy'//trim(adjustl(omega_name))//'.txt',1)
                    endif
                enddo
                 print *,act_cn
            endif
            if (mod(ii,1000)==0) then
                    print *,ii;
                   
            endif            
        !call exportf(time1,'time.txt')
            if (isnan(y(7,ii))) then
                    print *, ii
                    print *,'NAN'
                    goto 10
            endif
        ijk=0 
        !!!!!!!!!!!!!!!!!!!!!!!!!!
        do iib=104,con_num_all
            if(body(1).r(3*con_list(iib,3),ii)<1.2d0.and.body(1).r(3*con_list(iib,3),ii)>0.95d0) then                
                con_list(iib,10)=1
            endif
            if(body(1).r(3*con_list(iib,3),ii)>1.2d0.or.body(1).r(3*con_list(iib,3),ii)<0.95d0) then                
                con_list(iib,10)=0
            endif
        enddo
        
        do iib=3,103
            if(body(1).r(3*con_list(iib,3),ii)>1.2d0) then                
                con_list(iib,10)=0
            endif
        enddo
        con_num_act=sum(con_list(:,10))
        deallocate(Jac)                          
        deallocate(Equ,Equ1,Equ2,Equ3)             
        deallocate(Jac_temp,Equ_temp,Act_cn)                    
        ALLOCATE(Jac(F_num_all+body_num+con_num_Act,F_num_all+body_num+con_num_Act))                    
        allocate(Equ(F_num_all+body_num+con_num_Act),Equ1(F_num_all+body_num+con_num_Act),Equ2(F_num_all+body_num+con_num_Act),Equ3(F_num_all+body_num+con_num_Act))              
        allocate(Jac_temp(F_num_all+con_num_Act+body_num,F_num_all+con_num_Act+body_num),Equ_temp(F_num_all+con_num_Act+body_num),Act_cn(con_num_Act))
        
        mm=1
        do ll=1,con_num_all
            if(con_list(ll,10)==1) then
                Act_Cn(mm)=ll;
                mm=mm+1;
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
        if (Body(1).r(3*con_list(2,3),ii)>1.1d0) then
                Body(1).con_dof(2,3)=  Body(1).con_dof(2,3)-1
                body(1).con_dof(2,1)=9*Body(1).con_dof(2,3)-9+Body(1).con_dof(2,4);    
                print *, 'hello' ,Body(1).con_dof(2,3)
        endif
                
        
        
        ! call check
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        enddo!!!!!!!!!!一千步小循环结束
        print *, jj
        !
        ! do kk=1,body_num 
        !    write(omega_name,'(I3)') KK+2
        !    call exportf(transpose(body(KK).r(:,1:n_end:50)),'r'//trim(adjustl(omega_name))//'.txt',1)
        !    call exportf(transpose(body(KK).dr(:,1:n_end:50)),'dr'//trim(adjustl(omega_name))//'.txt',1)
        !    call exportf(transpose(body(KK).omega(:,1:n_end:50)),'omega'//trim(adjustl(omega_name))//'.txt',1)
        !enddo
        !call exportf(transpose(dy(:,1:n_end:50)),'dy.txt',1)
        !call exportf(transpose(y(:,1:n_end:50)),'y.txt',1)
        call exportf(time1(1:n_end:50),'time.txt',1)
        call exportf(transpose(lumda(:,1:n_end:50)),'lumda.txt',1)
        call exportf(transpose(lumda_quater(:,1:n_end:50)),'lumda_quater.txt',1)
        !call post
        y(:,1)=y(:,n_end+1);
        dy(:,1)=dy(:,n_end+1);
        time1(1)=time1(n_end+1)
    enddo!!!!!!!!大循环结束
10  call tt.toc('loop')
!100    call tt.toc('loop')
    !call exportf(r(:,1:n_end:10),'r.txt',1)
    !call exportf(dr(:,1:n_end:10),'dr.txt',1)
    !call exportf(x(:,1:n_end:1),'x.txt')
    !call exportf(dx(:,1:n_end:1),'dx.txt')
    !call exportf(y(:,1:ii:10),'y.txt')
    !call exportf(dy(:,1:ii:10),'dy.txt')
    !call exportf(ddy(:,1:ii:10),'ddy.txt')
    !call exportf(time1(1:ii:10),'time.txt')
    !call exportf(dx(:,1:ii:50),'dx.txt')
    ! Variables

    ! Body of recur
    !call post

    end program main

