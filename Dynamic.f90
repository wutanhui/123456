subroutine body_MF(Z,z0,y,dy,A1,Phi,Ki,Cm,Rou,Iz,Gra,RouA,Le)
    use OPERATORS_MOD
    use functions
    implicit none
    
    integer:: ii,jj,kk
    
    
    
    do ii=1,ele_num
    
        call Element_MF(Me,Fe,y(3*ii+5:3*ii+7),RouA,Le,A1,omega(3*ii-2:3*ii),dy(3*ii+5:3*ii+7),0,Iz,Gra,a0,Roup(:,ii),Phi,Cm,Ki,Rou)
        M(9*ii-8:9*ii,9*ii-8:9*ii)=Me;
        F(9*ii-8:9*ii)=Fe;
    enddo
    
    Z=matmul(transpose(G),matmul(M,G));
    Z0=matmul(transpose(G),matmul(M,-matmul(g0,ones))+F)
 endsubroutine body_MF