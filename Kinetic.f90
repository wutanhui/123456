subroutine Kinetic(y,dy,r,dr,Roup,dRoup,Omega,Omega_r,AI,Phi,Pi,Beta,Ele_num,Le)
    use functions
    use OPERATORS_MOD
    implicit none
    
    real(8),allocatable::y(:),dy(:),r(:),dr(:),Roup(:),dRoup(:),Omega(:),Omega_r(:),AI(:,:),PHi(:,:),Pi(:),Beta(:,:)
    
    real(8)::Le,Roup_T(3,3),Omega_T(3,3),Ki(3,3),Kqq2(3,3)
    integer::Ele_num,ii

    do ii=1,ele_num
        
        
        
    enddo
    
        
endsubroutine kinetic


subroutine recur_matrix(G,g0,y,dy,r,dr,AI,Phi,Roup,droup,Ki,omega,omega_r,Ele_num,Le)
    use functions
    use OPERATORS_MOD
    implicit none
    integer::Ele_num,ii,jj,kk
    real(8),allocatable::G(:,:),g0(:,:),y(:),dy(:),r(:),dr(:),AI(:,:),Phi(:,:),Roup(:),dRoup(:),Ki(:,:),omega(:),omega_r(:),T(:,:)
    real(8)::R_Q(3,4),A1(3,3),Phi_T(3,3),eye(3,3),Roup_T(3,3),Kqq2(3,3),beta1(3),beta2(3),Pi(3),Le
    
    
    allocate(  T(9*ele_num+9,1:9))
    Eye=0d0;Eye(1,1)=1d0;Eye(2,2)=1d0;Eye(3,3)=1d0;
    G=0d0;
    G0=0d0;
    G(1,1)=1d0;G(2,2)=1d0;G(3,3)=1d0;
    R_Q=2d0*R_Quater(y(4:7))
    G(4:6,4:7)=R_Q;
    G(7,8)=1d0;G(8,9)=1d0;G(9,10)=1d0;
    T=0d0;
    do ii=1,ele_num
        call Tensor(Roup_T,Roup(3*ii-2:3*ii))
        T(9*ii+1:9*ii+3,1:9)=reshape([eye,-Roup_T,phi(1:3,3*ii-2:3*ii)],[3,9]);
        T(9*ii+4:9*ii+6,4:6)=Eye;
        T(9*ii+4:9*ii+6,7:9)=AI(1:3,3*ii-2:3*ii).mt.Ki(1:3,3*ii-2:3*ii);   
        
        Kqq2(1,3)= dy(3*ii+6)*cos(y(3*ii+6))
        Kqq2(2,2)=-dy(3*ii+5)*sin(y(3*ii+5))
        Kqq2(2,3)=-dy(3*ii+5)*cos(y(3*ii+5))*cos(y(3*ii+6))+dy(3*ii+6)*sin(y(3*ii+5))*sin(y(3*ii+6))
        Kqq2(3,2)= dy(3*ii+5)*cos(y(3*ii+5))
        Kqq2(3,3)=-dy(3*ii+5)*sin(y(3*ii+5))*cos(y(3*ii+6))-dy(3*ii+6)*cos(y(3*ii+5))*sin(y(3*ii+6))  
        
        call Relative_Pi(Pi,y(3*ii+5:3*ii+7),dy(3*ii+5:3*ii+7),Le)
        
        beta1=cross(omega(3*ii-2:3*ii),cross(omega(3*ii-2:3*ii),Roup(3*ii-2:3*ii)))+2d0*cross(omega(3*ii-2:3*ii),droup(3*ii-2:3*ii))+PI
        
        beta2=(AI(1:3,3*ii-2:3*ii).mt.Kqq2.mt.dy(3*ii+5:3*ii+7))+cross(omega(3*ii-2:3*ii),omega_r(3*ii-2:3*ii))
        
        
        do kk=1,ii+1
            if (kk<ii+1) then
                if (kk==1) then
                    G(9*ii+1:9*ii+9,1:10)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,1:10);
                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
                else
                    G(9*ii+1:9*ii+9,3*kk+5:3*kk+7)=T(9*ii+1:9*ii+9,1:9).mt.G(9*ii-8:9*ii,3*kk+5:3*kk+7);
                    g0(9*ii+1:9*ii+9,kk)=T(9*ii+1:9*ii+9,1:9).mt.g0(9*ii-8:9*ii,kk);
                endif
            elseif (kk==ii+1) then
                G(9*ii+7:9*ii+9,3*kk+5:3*kk+7)=eye;
                g0(9*ii+1:9*ii+6,kk)=[beta1,beta2];
            endif;
        enddo
    enddo
    

endsubroutine recur_matrix


subroutine Abs_Omega(y,dy,omega,omega_r,Ki,AI,Ele_num)
    use functions
    use OPERATORS_MOD
    implicit none    
    real(8),allocatable::y(:),dy(:),omega(:),omega_r(:),AI(:,:),Ki(:,:)
    real(8):: R_Q(3,4)
    real(8):: Le
    integer::Ele_num,ii
    
    R_Q=2d0*R_Quater(y(4:7))
    omega(1:3)=R_Q.mt.dy(4:7);
    do ii=1,ele_num
        Ki(1:3,3*ii-2:3*ii)=K_R(y(3*ii+5:3*ii+7))
        omega_r(3*ii-2:3*ii)=AI(1:3,3*ii-2:3*ii).mt.Ki(1:3,3*ii-2:3*ii).mt.dy(3*ii+5:3*ii+7)
        omega(3*ii+1:3*ii+3)=omega(3*ii-2:3*ii)+omega_r(3*ii-2:3*ii);      
    enddo       
endsubroutine Abs_omega



subroutine Abs_Vel(y,dy,dr,dRoup,Phi,AI,Ele_num,Le)
    use OPERATORS_MOD
    implicit none    
    real(8),allocatable::y(:),dy(:),dr(:),dRoup(:),Phi(:,:),AI(:,:)
    real(8)::Le
    integer::Ele_num,ii
    dr(1:3)=dy(1:3);
    do ii=1,ele_num
        call Relative_Phi(dRoup(3*ii-2:3*ii),Phi(:,9*ii-9:9*ii),y(3*ii+5:3*ii+7),dy(3*ii+5:3*ii+7),Le)
        dRoup (3*ii-2:3*ii)=AI(1:3,3*ii-2:3*ii).mt.dRoup(3*ii-2:3*ii)               
        dr(3*ii+1:3*ii+3)=dr(3*ii-2:3*ii)+dRoup(3*ii-2:3*ii);      
    enddo       
endsubroutine Abs_Vel



subroutine Abs_Coor(y,r,Roup,AI,Ele_num,Le)
    use functions
    use OPERATORS_MOD
    implicit none    
    real(8),allocatable::y(:),r(:),Roup(:),AI(:,:)
    real(8)::Le
    integer::Ele_num,ii
    r(1:3)=y(1:3);
    AI(1:3,1:3)=dir_cos_Q(y(4:7))
    do ii=1,ele_num
        call Relative_coor(Roup(3*ii-2:3*ii),y(3*ii+5:3*ii+7),Le)
        Roup (3*ii-2:3*ii)=AI(1:3,3*ii-2:3*ii).mt. Roup(3*ii-2:3*ii)               
        AI(1:3,3*ii+1:3*ii+3)=AI(1:3,3*ii-2:3*ii).mt.dir_cos(y(3*ii+5:3*ii+7))
        r(3*ii+1:3*ii+3)=r(3*ii-2:3*ii)+Roup(3*ii-2:3*ii);      
    enddo       
endsubroutine Abs_Coor



subroutine Relative_coor(Roup,a,Le)
    implicit none        
    real(8):: Roup(3),a(3),Le
    real(8):: bi1,bi2
    bi1=a(1)+a(2)
    bi2=a(1)-a(2)
        
    Roup(1)= Le*( (a(2))/2d0-(a(2)**3)/24d0+(a(2)**5)/720d0-(a(2)**7)/40320d0+a(2)**9/3628800d0)
    Roup(2)=-Le/2d0*( (bi1+bi2)/2d0-(bi1**3+bi2**3)/24d0+(bi1**5+bi2**5)/720d0-(bi1**7+bi2**7)/40320d0+(bi1**9+bi2**9)/3628800d0)
    Roup(3)= Le+Le/2d0*(-(bi1**2+bi2**2)/6d0 +(bi1**4+bi2**4)/120d0 -(bi1**6+bi2**6)/5040d0+(bi1**8+bi2**8)/362880d0)
    
endsubroutine Relative_coor

subroutine Relative_Phi(dRoup,Phi,a,da,Le)
    use OPERATORS_MOD
    implicit none
        real(8):: dRoup(3),phi(3,3),a(3),da(3),Le
        real(8):: bi1,bi2,dbi1,dbi2
        bi1=a(1)+a(2)
        bi2=a(1)-a(2)
        dbi1=da(1)+da(2)
        dbi2=da(1)-da(2)
        
        Phi=0d0;
        Phi(1,2)= Le/2d0+Le*(-(a(2)**2)/8d0+(a(2)**4)/144d0-(a(2)**6)/5760d0+a(2)**8/403200d0)
        Phi(2,1)=-Le/2d0-Le/2d0*(-(bi1**2+bi2**2)/8d0+(bi1**4+bi2**4)/144d0-(bi1**6+bi2**6)/5760d0+(bi1**8+bi2**8)/403200d0)
        Phi(2,2)=       -Le/2d0*(-(bi1**2-bi2**2)/8d0+(bi1**4-bi2**4)/144d0-(bi1**6-bi2**6)/5760d0/8d0+(bi1**8-bi2**8)/403200d0)
        Phi(3,1)= Le/2d0*(-(bi1+bi2)/3d0+(bi1**3+bi2**3)/30d0-(bi1**5+bi2**5)/840d0+(bi1**7+bi2**7)/45360d0)
        Phi(3,2)= Le/2d0*(-(bi1-bi2)/3d0+(bi1**3-bi2**3)/30d0-(bi1**5-bi2**5)/840d0+(bi1**7-bi2**7)/45360d0)
        
        dRoup=phi.mt.da
endsubroutine Relative_Phi

subroutine Relative_Pi(Pi,a,da,Le)
    implicit none
        real(8):: Pi(3),a(3),da(3),Le
        real(8):: bi1,bi2,dbi1,dbi2
        bi1=a(1)+a(2)
        bi2=a(1)-a(2)
        dbi1=da(1)+da(2)
        dbi2=da(1)-da(2)

        pi(1)= Le*(-a(2)/4d0+(a(2)**3)/36d0-(a(2)**5)/960d0+a(2)**7/50400d0)*da(2)**2
        Pi(2)=-Le/2d0*(-(bi1*dbi1**2+bi2*dbi2**2)/4d0+(bi1**3*dbi1**2+bi2**3*dbi2**2)/36d0-(bi1**5*dbi1**2+bi2**5*dbi2**2)/960d0+(bi1**7*dbi1**2+bi2**7*dbi2**2)/50400d0)
        Pi(3)= Le/2d0*(-(dbi1**2+dbi2**2)/3d0+(bi1**2*dbi1**2+bi2**2*dbi2**2)/10d0-(bi1**4*dbi1**2+bi2**4*dbi2**2)/168d0+(bi1**6*dbi1**2+bi2**6*dbi2**2)/6480d0)
endsubroutine Relative_Pi



!subroutine Recur_Matrix(G,g,r,dr,beta,Roup,Phi,Pi,Ki,Ele_num)
!
!
!
!endsubroutine Recur_Matix