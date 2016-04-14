SUBROUTINE contact_detect(body_num)!2014.11.17
    use body_module
    implicit none
    integer::ii,jj,kk,ll,mm,contact_note,body_num;
    real(8)::distance1,distance2,length,normal(3);!distance1点到直线的法向距离（只有正），distance2点到在直线上投影点位置到线段起始点距离（有正有负），normal点到直线的法线方向
     real(8)::kesi1,kesi2,delta
     !real(8)::F_contact(3)
    real(8),allocatable:: r1(:),r2(:),r3(:),F_contact1(:),F_contact2(:),F_contact3(:)
    allocate(r1(3*body(1).node_num),r2(3*body(2).node_num),r3(3*body(3).node_num),F_contact1(9*body(1).ele_num),F_contact2(9*body(2).ele_num),F_contact3(9*body(3).ele_num))
   
    
    
    r1=body(1).r_temp
    r2=body(2).r_temp    
    r3=body(3).r_temp  
    
    F_contact1=0d0;F_contact2=0d0;F_contact3=0d0;
    do ii=1,body_num
        body(ii).F_contact=0d0;
    enddo
    
    !do ii=11,body(1).Node_num,1       
    !    do jj=1,body(2).Ele_num
    !        contact_note=0
    !        if(r1(3*ii-2)>r3(3*jj-2)-0.03d0.and.r1(3*ii-2)<=r3(3*jj+1)+0.03d0.and.r1(3*ii-1)>r3(3*jj-1)-0.03d0.and.r1(3*ii-1)<r3(3*jj+2)+0.03d0.and.r1(3*ii)>r3(3*jj)-0.03d0.and.r1(3*ii)<r3(3*jj+3)+0.03d0) then
    !            call detect_node_line(contact_note,distance1,distance2,length,normal,r1(3*ii-2:3*ii),r3(3*ii-2:3*ii),r3(3*jj+1:3*jj+3))
    !            if(contact_note==1) then!检测到接触，记录下梁1的接触点A，梁2的接触单元B，并计算出接触力的大小，最终施加到A点，及B单元的1，2节点。
    !                F_contact1(9*ii-8:9*ii-6)= F_contact1(9*ii-8:9*ii-6)-2.5d3*normal*(distance1-0.04d0)
    !                F_contact3(9*jj-8:9*jj-6)= F_contact3(9*jj-8:9*jj-6)+2.5d3*normal*(distance1-0.04d0)*(1d0-distance2/length)
    !                F_contact3(9*jj+1:9*jj+3)= F_contact3(9*jj+1:9*jj+3)+2.5d3*normal*(distance1-0.04d0)*distance2/length
    !            endif
    !        endif
    !    enddo
    !enddo
    !do ii=11,body(2).Node_num,1       
    !    do jj=1,body(3).Ele_num
    !        contact_note=0
    !        if(r2(3*ii-2)>r1(3*jj-2)-0.03d0.and.r2(3*ii-2)<=r1(3*jj+1)+0.03d0.and.r2(3*ii-1)>r1(3*jj-1)-0.03d0.and.r2(3*ii-1)<r1(3*jj+2)+0.03d0.and.r2(3*ii)>r1(3*jj)-0.03d0.and.r2(3*ii)<r1(3*jj+3)+0.03d0) then
    !            call detect_node_line(contact_note,distance1,distance2,length,normal,r2(3*ii-2:3*ii),r1(3*ii-2:3*ii),r1(3*jj+1:3*jj+3))
    !            if(contact_note==1) then
    !                F_contact2(9*ii-8:9*ii-6)= F_contact2(9*ii-8:9*ii-6)-2.5d3*normal*(distance1-0.04d0)
    !                F_contact1(9*jj-8:9*jj-6)= F_contact1(9*jj-8:9*jj-6)+2.5d3*normal*(distance1-0.04d0)*(1d0-distance2/length)
    !                F_contact1(9*jj+1:9*jj+3)= F_contact1(9*jj+1:9*jj+3)+2.5d3*normal*(distance1-0.04d0)*distance2/length
    !            endif
    !        endif
    !    enddo
    !enddo    
    !do ii=11,body(3).Node_num,1       
    !    do jj=1,body(1).Ele_num
    !        contact_note=0
    !        if(r3(3*ii-2)>r2(3*jj-2)-0.03d0.and.r3(3*ii-2)<=r2(3*jj+1)+0.03d0.and.r3(3*ii-1)>r2(3*jj-1)-0.03d0.and.r3(3*ii-1)<r2(3*jj+2)+0.03d0.and.r3(3*ii)>r2(3*jj)-0.03d0.and.r3(3*ii)<r2(3*jj+3)+0.03d0) then
    !            call detect_node_line(contact_note,distance1,distance2,length,normal,r3(3*ii-2:3*ii),r2(3*ii-2:3*ii),r2(3*jj+1:3*jj+3))
    !            if(contact_note==1) then
    !                F_contact3(9*ii-8:9*ii-6)= F_contact3(9*ii-8:9*ii-6)-2.5d3*normal*(distance1-0.04d0)!5*(body(3).dr(3*ii-2:3*ii)-body(2).dr(3*jj-2:3*jj)*(1d0-distance2/length)-body(2).dr(3*jj+1:3*jj+3)*distance2/length)
    !                F_contact2(9*jj-8:9*jj-6)= F_contact2(9*jj-8:9*jj-6)+2.5d3*normal*(distance1-0.04d0)*(1d0-distance2/length)!-5*(body(3).dr(3*ii-2:3*ii)-body(2).dr(3*jj-2:3*jj)*(1d0-distance2/length)-body(2).dr(3*jj+1:3*jj+3)*distance2/length)*(1d0-distance2/length)
    !                F_contact2(9*jj+1:9*jj+3)= F_contact2(9*jj+1:9*jj+3)+2.5d3*normal*(distance1-0.04d0)*distance2/length!-5*(body(3).dr(3*ii-2:3*ii)-body(2).dr(3*jj-2:3*jj)*(1d0-distance2/length)-body(2).dr(3*jj+1:3*jj+3)*distance2/length)*distance2/length
    !            endif
    !        endif
    !    enddo
    !enddo    
    !body(1).z1=body(1).z1+matmul(transpose(body(1).G),F_contact1)
    !body(2).z1=body(2).z1+matmul(transpose(body(2).G),F_contact2)
    !body(3).z1=body(3).z1+matmul(transpose(body(3).G),F_contact3)
    
    do ii=4,body_num-9,3
        do jj=ii,ii+2
            do kk=1,body(jj).ele_num
                do ll=ii+3,ii+5
                    do mm=1,body(ll).ele_num
                        call detect_line_line(contact_note,kesi1,kesi2,normal,body(jj).r_temp(3*kk-2:3*kk),body(jj).r_temp(3*kk+1:3*kk+3),body(ll).r_temp(3*mm-2:3*mm),body(ll).r_temp(3*mm+1:3*mm+3))
                        if (contact_note==1) then
                            delta=norm2(normal)-0.04d0;
                            body(jj).F_contact(9*kk-8:9*kk-6)=body(jj).F_contact(9*kk-8:9*kk-6)-normal*(1-kesi1)*delta*5d4
                            body(jj).F_contact(9*kk+1:9*kk+3)=body(jj).F_contact(9*kk+1:9*kk+3)-normal*kesi1    *delta*5d4
                            body(ll).F_contact(9*mm-8:9*mm-6)=body(ll).F_contact(9*mm-8:9*mm-6)+normal*(1-kesi2)*delta*5d4
                            body(ll).F_contact(9*mm+1:9*mm+3)=body(ll).F_contact(9*mm+1:9*mm+3)+normal*kesi2    *delta*5d4
                             !print *, 'contact'
                        endif
                    enddo
                enddo
            enddo
        enddo
    enddo
    
    do ii=4,body_num-6
        body(ii).z1=body(ii).z1+matmul(transpose(body(ii).G),body(ii).F_contact(1:9*body(ii).ele_num))+matmul(transpose(body(ii).GG),body(ii).F_contact(9*body(ii).ele_num+1:9*body(ii).ele_num+9))        
    enddo
    
endsubroutine contact_detect

subroutine detect_node_line(contact_note,distance1,distance2,length,normal1,a,b1,b2)!!!2014.11.17
    implicit none
    interface 
       function cross(a,b)
            real(8)::a(3),b(3),cross(3);
       end function
    end interface
    
    integer::contact_note
   ! real(8)::cross
    real(8)::distance1,distance2,length;
    real(8)::tangent1(3),tangent2(3),normal1(3),normal2(3);
    real(8)::a(3),b1(3),b2(3);
    
    length=norm2(b2-b1)
    tangent1=(b2-b1)/norm2(b2-b1);
    tangent2=(a -b1)/norm2(a -b1);
    normal2=cross(tangent1,tangent2)/norm2(cross(tangent1,tangent2));
    normal1=cross(normal2,tangent1);
    distance1=dot_product(a-b1,normal1);!
    distance2=dot_product(a-b1,tangent1);
    
    !contact_note=0;
    if(distance1<0.04d0.and.distance2<(norm2(b2-b1)).and.distance2>=0) then
        contact_note=1;
    endif

endsubroutine 

 function cross(a,b)
    implicit none
    real(8)::A(3),B(3);
    REAL(8)::cross(3)
    Cross(1) = A(2) * B(3) - A(3) * B(2)
    Cross(2) = A(3) * B(1) - A(1) * B(3)
    Cross(3) = A(1) * B(2) - A(2) * B(1)
    !return
 endfunction cross
 
 subroutine detect_line_line(contact_note,s,t,normal,P1,P2,Q1,Q2)!!!2014.11.24
    implicit none
    !interface 
    !   function cross(a,b)
    !        real(8)::a(3),b(3),cross(3);
    !   end function
    !end interface
    
    integer::contact_note
   ! real(8)::cross
    real(8)::s,t;
    real(8)::P1(3),P2(3),Q1(3),Q2(3);
    real(8)::d1(3),d2(3),r(3)
    real(8)::normal(3)
    real(8)::a,b,c,d,e,f;
    
    d1=P2-P1;
    d2=Q2-Q1;
    r= P1-Q1;
    a=dot_product(d1,d1)
    b=dot_product(d1,d2)
    c=dot_product(d1,r)
    e=dot_product(d2,d2)
    f=dot_product(d2,r)
    d=a*e-b**2
    contact_note=0
    if (d/=0) then
        s=(b*f-c*e)/d;
        t=(a*f-b*c)/d;
    endif    
    if (d==0) then
       s=0;
       t=f/e;
    endif
    !print *, 'contact0'
    if (s>=0d0.and.s<1d0.and.t>=0d0.and.t<1d0) then
        normal=r+s*d1-t*d2;
        !print *, 'contact1'
        if(norm2(normal)<=0.04d0) then
            contact_note=1
            !print *, 'contact2'
        endif
    endif

endsubroutine 
 
 