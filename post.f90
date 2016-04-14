subroutine post
    use Model_module
    use body_module
    use OPERATORS_MOD
    implicit none
    interface 
        subroutine recur_post(r_temp,dr_temp,ddr_temp,y_temp,dy_temp,ddy_temp,omega_temp,d_omega_temp,Le,ele_num_all)
            real(8),allocatable::r_temp(:),dr_temp(:),ddr_temp(:),y_temp(:),dy_temp(:),ddy_temp(:),omega_temp(:),d_omega_temp(:)
            real(8)::Le
            integer::ele_num_all
        end subroutine 
    end interface
    integer::ii,jj,kk,I_begin,I_end
    character(10):: omega_name
    
    
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
        do jj=1,n_end,1
            body(II).y_temp=body(II).y(:,jj)
            body(II).dy_temp=body(II).dy(:,jj)
            body(II).ddy_temp=body(II).ddy(:,jj)
            call recur_post(body(II).r_temp,body(II).dr_temp,body(II).ddr_temp,body(II).y_temp,body(II).dy_temp,body(II).ddy_temp,body(II).omega_temp,body(II).domega_temp,body(II).Le,body(II).ele_num)
            body(II).r(:,jj)=body(II).r_temp
            body(II).dr(:,jj)=body(II).dr_temp
            body(II).omega(:,jj)=body(II).omega_temp
        enddo
        write(omega_name,'(I3)') II
        call exportf(transpose(body(ii).r(:,1:n_end:1)),'r'//trim(adjustl(omega_name))//'.txt',1)
        call exportf(transpose(body(ii).dr(:,1:n_end:1)),'dr'//trim(adjustl(omega_name))//'.txt',1)
        call exportf(transpose(body(ii).omega(:,1:n_end:1)),'omega'//trim(adjustl(omega_name))//'.txt',1)
        call exportf(transpose(body(ii).y(:,1:n_end:1)),'y'//trim(adjustl(omega_name))//'.txt',1)
    enddo
    
    close(11);
    close(12)
    close(13)
    close(14)
    close(15)
    !call exportf(r,'rr.txt')
    !call exportf(dr,'drr.txt')
    call exportf(time1,'time.txt')
endsubroutine post    
    
subroutine recur_post(r,dr,ddr,y,dy,ddy,omega,d_omega,Le,ele_num_all)
    use functions
    use OPERATORS_MOD
    implicit none
    integer::ii,ele_num_all
    real(8)::Le
    real(8)::r(:),dr(:),ddr(:),y(:),dy(:),ddy(:),omega(:),d_omega(:)
    
    !上为传入程序的，下为子程序中用到的
    real(8)::beta1(3),beta2(3),Ki(3,3),A1(3,3),A2(3,3),Phi(3,3),Pi(3),omega_r(3),Omega_T(3,3),Roup_T(3,3),Kqq2(3,3)
    real(8)::R_Q(3,4)
    real(8),allocatable::a(:),da(:),dda(:),roup(:,:),droup(:,:)
    allocate(a(3*ele_num_all),da(3*ele_num_all),dda(3*ele_num_all),roup(3,ele_num_all),droup(3,ele_num_all))
    
    r(1:3)=y(1:3);
    dr(1:3)=dy(1:3);
    !Ki=reshape([1d0,0d0,0d0,0d0,cos(y(4)),sin(y(4)),sin(y(5)),-cos(y(5))*sin(y(4)),cos(y(5))*cos(y(4))],[3,3])
    
    R_Q=2d0*R_Quater(y(4:7))
    
    Kqq2=0d0
    omega(1:3)=R_Q.mt.dy(4:7);   
    d_omega(1:3)=R_Q.mt.ddy(4:7);
    
    a  =   y(8:3*ele_num_all+7);
    da =  dy(8:3*ele_num_all+7);
    dda= ddy(8:3*ele_num_all+7);
    A1=dir_cos_q(y(4:7))
    
    do ii=1,ele_num_all
        call Phi_Pi(droup(:,ii),Roup(:,ii),PHI,Pi,a(3*ii-2:3*ii),da(3*ii-2:3*ii),Le)
        Roup(:,ii)=A1.mt.Roup(:,ii)
        dRoup(:,ii)=A1.mt.dRoup(:,ii)
        PHI=A1.mt.Phi
        Pi=A1.mt.Pi
        Ki=K_R(a(3*ii-2:3*ii))
        
        omega_r=A1.mt.Ki.mt.da(3*ii-2:3*ii)
        call Tensor(Roup_T,Roup(:,ii));
        call Tensor(omega_T,omega(3*ii-2:3*ii));

        Kqq2(1,3)= da(3*ii-1)*cos(a(3*ii-1))
        Kqq2(2,2)=-da(3*ii-2)*sin(a(3*ii-2))
        Kqq2(2,3)=-da(3*ii-2)*cos(a(3*ii-2))*cos(a(3*ii-1))+da(3*ii-1)*sin(a(3*ii-2))*sin(a(3*ii-1))
        Kqq2(3,2)= da(3*ii-2)*cos(a(3*ii-2))
        Kqq2(3,3)=-da(3*ii-2)*sin(a(3*ii-2))*cos(a(3*ii-1))-da(3*ii-1)*cos(a(3*ii-2))*sin(a(3*ii-1))
        
        A2=Dir_COS(a(3*ii-2:3*ii))

        beta1=(omega_T.mt.(omega_T.mt.Roup(:,ii)))+2d0*(omega_T.mt.droup(:,ii))+PI
        
        beta2=(A1.mt.Kqq2.mt.da(3*ii-2:3*ii))+(omega_T.mt.omega_r)!!!改
        
        
        
        r(3*ii+1:3*ii+3)=r(3*ii-2:3*ii)+Roup(:,ii);
        dr(3*ii+1:3*ii+3)=dr(3*ii-2:3*ii)-(Roup_T.mt.omega(3*ii-2:3*ii))+dRoup(:,ii);
        ddr(3*ii+1:3*ii+3)=ddr(3*ii-2:3*ii)-(Roup_T.mt.d_omega(3*ii-2:3*ii))+(PHI.mt.dda(3*ii-2:3*ii))+beta1;
        
        omega(3*ii+1:3*ii+3)=omega(3*ii-2:3*ii)+omega_r;
        d_omega(3*ii+1:3*ii+3)=d_omega(3*ii-2:3*ii)+(A1.mt.Ki.mt.dda(3*ii-2:3*ii))+beta2;
        
        A1=A1.mt.A2;
    enddo
 

endsubroutine recur_post