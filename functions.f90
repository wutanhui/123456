 module functions
    implicit none
    
    interface cross
        module procedure cross
    endinterface
    interface Dir_Cos
        module procedure Dir_cos
    endinterface
    !#############################
    interface Dir_Cos_Q
        module procedure dir_cos_Q
    endinterface 
    interface R_Quater
        module procedure R_Quater
    endinterface
    interface L_Quater
        module procedure L_Quater
    endinterface
    interface Dir_Cos_Q1
        module procedure dir_cos_Q1
    endinterface
    interface Dir_Cos_Q2
        module procedure dir_cos_Q2
    endinterface
    interface Dir_Cos_Q3
        module procedure dir_cos_Q3
    endinterface
    interface Dir_Cos_Q0
        module procedure dir_cos_Q0
    endinterface
    !#############################
    interface K_B
        module procedure K_B
    endinterface
    interface K_R
        module procedure K_R
    endinterface
    !################################
    interface dir_Cos_P1
        module procedure dir_Cos_P1
    endinterface   dir_Cos_P1
    interface dir_Cos_P2
        module procedure dir_Cos_P2
    endinterface   dir_Cos_P2
    interface dir_Cos_P3
        module procedure dir_Cos_P3
    endinterface   dir_Cos_P3  
    !################################
    
    !################################
    interface dir_Cos_P1P1
        module procedure dir_Cos_P1P1
    endinterface   dir_Cos_P1P1
    interface dir_Cos_P1P2
        module procedure dir_Cos_P1P2
    endinterface   dir_Cos_P1P2
    interface dir_Cos_P1P3
        module procedure dir_Cos_P1P3
    endinterface   dir_Cos_P1P3
    interface dir_Cos_P2P2
        module procedure dir_Cos_P2P2
    endinterface   dir_Cos_P2P2
    interface dir_Cos_P2P3
        module procedure dir_Cos_P2P3
    endinterface   dir_Cos_P2P3
    interface dir_Cos_P3P3
        module procedure dir_Cos_P3P3
    endinterface   dir_Cos_P3P3
    !################################
    interface K_R_P1
        module procedure K_R_P1
    endinterface
    interface K_R_P2
        module procedure K_R_P2
    endinterface
    !################################
  contains
    
    function cross(a,b)
        implicit none
        real(8)::A(3),B(3);
        REAL(8)::cross(3)
        Cross(1) = A(2) * B(3) - A(3) * B(2)
        Cross(2) = A(3) * B(1) - A(1) * B(3)
        Cross(3) = A(1) * B(2) - A(2) * B(1)
    endfunction cross
    

 
    function Dir_COS(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS(3,3)
        Dir_COS=reshape([cos(a(2))*cos(a(3)), sin(a(1))*sin(a(2))*cos(a(3))+cos(a(1))*sin(a(3)), -cos(a(1))*sin(a(2))*cos(a(3))+sin(a(1))*sin(a(3)),&
                       &-cos(a(2))*sin(a(3)),-sin(a(1))*sin(a(2))*sin(a(3))+cos(a(1))*cos(a(3)),  cos(a(1))*sin(a(2))*sin(a(3))+sin(a(1))*cos(a(3)),&
                       & sin(a(2)),          -sin(a(1))*cos(a(2)),                                cos(a(1))*cos(a(2))],[3,3])
    endfunction Dir_COS
 
    function Dir_COS_Q(A)
        implicit none
        real(8)::A(4);
        REAL(8)::Dir_COS_Q(3,3)
        Dir_COS_Q=reshape([2d0*(a(1)**2+a(2)**2)-1d0, 2d0*(a(2)*a(3)+a(1)*a(4)),  2d0*(a(2)*a(4)-a(1)*a(3)),&
                         & 2d0*(a(2)*a(3)-a(1)*a(4)), 2d0*(a(1)**2+a(3)**2)-1d0,  2d0*(a(3)*a(4)+a(1)*a(2)),&
                         & 2d0*(a(2)*a(4)+a(1)*a(3)), 2d0*(a(3)*a(4)-a(1)*a(2)),  2d0*(a(1)**2+a(4)**2)-1d0],[3,3])
    endfunction Dir_COS_Q 

    function Dir_COS_Q0(A)
        implicit none
        real(8)::A(4);
        REAL(8)::Dir_COS_Q0(3,3)
        Dir_COS_Q0=reshape([4d0*a(1), 2d0*a(4),  -2d0*a(3),&
                         & -2d0*a(4), 4d0*a(1),   2d0*a(2),&
                         &  2d0*a(3),-2d0*a(2),   4d0*a(1)],[3,3])
    endfunction Dir_COS_Q0 
    function Dir_COS_Q1(A)
        implicit none
        real(8)::A(4);
        REAL(8)::Dir_COS_Q1(3,3)
        Dir_COS_Q1=reshape([4d0*a(2),  2d0*a(3),  2d0*a(4),&
                         &  2d0*a(3),  0d0,       2d0*a(1),&
                         &  2d0*a(4), -2d0*a(1),  0d0],[3,3])
    endfunction Dir_COS_Q1
    function Dir_COS_Q2(A)
        implicit none
        real(8)::A(4);
        REAL(8)::Dir_COS_Q2(3,3)
        Dir_COS_Q2=reshape([0d0,      2d0*a(2),  -2d0*a(1),&
                         &  2d0*a(2), 4d0*a(3),   2d0*a(4),&
                         &  2d0*a(1), 2d0*a(4),   0d0     ],[3,3])
    endfunction Dir_COS_Q2
    function Dir_COS_Q3(A)
        implicit none
        real(8)::A(4);
        REAL(8)::Dir_COS_Q3(3,3)
        Dir_COS_Q3=reshape([0d0,      2d0*a(1),   2d0*a(2),&
                         & -2d0*a(1), 0d0,        2d0*a(3),&
                         &  2d0*a(2), 2d0*a(3),   4d0*a(4) ],[3,3])
    endfunction Dir_COS_Q3
    function R_Quater(A)
        implicit none
        real(8)::A(4);
        REAL(8)::R_Quater(3,4)
        R_Quater=reshape([-A(2), -A(3), -A(4),&
                         & A(1),  A(4), -A(3),&
                         &-A(4),  A(1),  A(2),&
                         & A(3), -A(2),  A(1)],[3,4])
    endfunction R_Quater
    
    function L_Quater(A)
        implicit none
        real(8)::A(4);
        REAL(8)::L_Quater(3,4)
        L_Quater=reshape([-A(2), -A(3), -A(4),&
                         & A(1), -A(4),  A(3),&
                         & A(4),  A(1), -A(2),&
                         &-A(3),  A(2),  A(1)],[3,4])
    endfunction L_Quater
    
    function K_R(A)
        implicit none
        real(8)::A(3);
        REAL(8)::K_R(3,3)           
        K_R=reshape([1d0,0d0,0d0,0d0,cos(A(1)),sin(A(1)),sin(A(2)),-cos(A(2))*sin(A(1)),cos(A(2))*cos(A(1))],[3,3])
    endfunction K_R
    
    function K_B(A)
        implicit none
        real(8)::A(3);
        REAL(8)::K_B(3,3)
        
        K_B=reshape([cos(A(2))*cos(A(3)),-cos(A(2))*sin(A(3)),sin(A(2)),&
                    &sin(A(3)),           cos(A(3)),          0d0,&
                    &0d0,                 0d0,                1d0],[3,3])
    endfunction K_B
     
    function Dir_COS_P1(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P1(3,3)
        Dir_COS_P1=reshape([0d0, cos(A(1))*sin(A(2))*cos(A(3))-sin(A(1))*sin(A(3)), sin(A(1))*sin(A(2))*cos(A(3))+cos(A(1))*sin(A(3)),&
                           &0d0,-cos(A(1))*sin(A(2))*sin(A(3))-sin(A(1))*cos(A(3)),-sin(A(1))*sin(A(2))*sin(A(3))+cos(A(1))*cos(A(3)),&
                           &0d0,-cos(A(1))*cos(A(2)),                              -sin(A(1))*cos(A(2))],[3,3])
    endfunction Dir_COS_P1
    
    function Dir_COS_P2(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P2(3,3)
        Dir_COS_P2=reshape([-sin(A(2))*cos(A(3)), sin(A(1))*cos(A(2))*cos(A(3)),-cos(A(1))*cos(A(2))*cos(A(3)),&
                           & sin(A(2))*sin(A(3)),-sin(A(1))*cos(A(2))*sin(A(3)), cos(A(1))*cos(A(2))*sin(A(3)),&
                           & cos(A(2)),           sin(A(1))*sin(A(2)),          -cos(A(1))*sin(A(2))],[3,3])
    endfunction Dir_COS_P2

    function Dir_COS_P3(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P3(3,3)
        Dir_COS_P3=reshape([-cos(A(2))*sin(A(3)),-sin(A(1))*sin(A(2))*sin(A(3))+cos(A(1))*cos(A(3)), cos(A(1))*sin(A(2))*sin(A(3))+sin(A(1))*cos(A(3)),&
                           &-cos(A(2))*cos(A(3)),-sin(A(1))*sin(A(2))*cos(A(3))-cos(A(1))*sin(A(3)), cos(A(1))*sin(A(2))*cos(A(3))-sin(A(1))*sin(A(3)),&
                           & 0d0,                 0d0,                                               0d0],[3,3])
    endfunction Dir_COS_P3
 !________________方向余弦的二阶导_________________________!
    function Dir_COS_P1P1(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P1P1(3,3)
        Dir_COS_P1P1=reshape([0d0,-sin(A(1))*sin(A(2))*cos(A(3))-cos(A(1))*sin(A(3)), cos(A(1))*sin(A(2))*cos(A(3))-sin(A(1))*sin(A(3)),&
                             &0d0, sin(A(1))*sin(A(2))*sin(A(3))-cos(A(1))*cos(A(3)),-cos(A(1))*sin(A(2))*sin(A(3))-sin(A(1))*cos(A(3)),&
                             &0d0, sin(A(1))*cos(A(2)),                              -cos(A(1))*cos(A(2))],[3,3])
    endfunction Dir_COS_P1P1
    
    function Dir_COS_P1P2(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P1P2(3,3)
        Dir_COS_P1P2=reshape([0d0, cos(A(1))*cos(A(2))*cos(A(3)), sin(A(1))*cos(A(2))*cos(A(3)),&
                             &0d0,-cos(A(1))*cos(A(2))*sin(A(3)),-sin(A(1))*cos(A(2))*sin(A(3)),&
                             &0d0, cos(A(1))*sin(A(2)),           sin(A(1))*sin(A(2))],[3,3])
    endfunction Dir_COS_P1P2
    
    function Dir_COS_P1P3(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P1P3(3,3)
        Dir_COS_P1P3=reshape([0d0,-cos(A(1))*sin(A(2))*sin(A(3))-sin(A(1))*cos(A(3)),-sin(A(1))*sin(A(2))*sin(A(3))+cos(A(1))*cos(A(3)),&
                             &0d0,-cos(A(1))*sin(A(2))*cos(A(3))+sin(A(1))*sin(A(3)),-sin(A(1))*sin(A(2))*cos(A(3))-cos(A(1))*sin(A(3)),&
                             &0d0, 0d0,                                               0d0],[3,3])
    endfunction Dir_COS_P1P3   
    
     function Dir_COS_P2P2(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P2P2(3,3)
        Dir_COS_P2P2=reshape([-cos(A(2))*cos(A(3)),-sin(A(1))*sin(A(2))*cos(A(3)), cos(A(1))*sin(A(2))*cos(A(3)),&
                             & cos(A(2))*sin(A(3)), sin(A(1))*sin(A(2))*sin(A(3)),-cos(A(1))*sin(A(2))*sin(A(3)),&
                             &-sin(A(2)),           sin(A(1))*cos(A(2)),          -cos(A(1))*cos(A(2))],[3,3])
     endfunction Dir_COS_P2P2
    function Dir_COS_P2P3(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P2P3(3,3)
        Dir_COS_P2P3=reshape([ sin(A(2))*sin(A(3)),-sin(A(1))*cos(A(2))*sin(A(3)), cos(A(1))*cos(A(2))*sin(A(3)),&
                             & sin(A(2))*cos(A(3)),-sin(A(1))*cos(A(2))*cos(A(3)), cos(A(1))*cos(A(2))*cos(A(3)),&
                             & 0d0,                 0d0,                           0d0],[3,3])
    endfunction Dir_COS_P2P3
    function Dir_COS_P3P3(A)
        implicit none
        real(8)::A(3);
        REAL(8)::Dir_COS_P3P3(3,3)
        Dir_COS_P3P3=reshape([-cos(A(2))*cos(A(3)),-sin(A(1))*sin(A(2))*cos(A(3))-cos(A(1))*sin(A(3)), cos(A(1))*sin(A(2))*cos(A(3))-sin(A(1))*sin(A(3)),&
                             & cos(A(2))*sin(A(3)), sin(A(1))*sin(A(2))*sin(A(3))-cos(A(1))*cos(A(3)),-cos(A(1))*sin(A(2))*sin(A(3))-sin(A(1))*cos(A(3)),&
                             & 0d0,                 0d0,                                               0d0],[3,3])
    endfunction Dir_COS_P3P3    
  !________________方向余弦的二阶导_________________________!   
    function K_R_P1(A)
        implicit none
        real(8)::A(3);
        REAL(8)::K_R_P1(3,3)
        K_R_P1=reshape([0d0,  0d0,                  0d0,      &
                       &0d0, -sin(A(1)),            cos(A(1)),&
                       &0d0, -cos(A(2))*cos(A(1)), -cos(A(2))*sin(A(1))],[3,3])
    endfunction K_R_P1
    
    function K_R_P2(A)
        implicit none
        real(8)::A(3);
        REAL(8)::K_R_P2(3,3)
        K_R_P2=reshape([0d0,       0d0,                 0d0,&
                       &0d0,       0d0,                 0d0,&
                       &cos(A(2)), sin(A(2))*sin(A(1)),-sin(A(2))*cos(A(1))],[3,3])
    endfunction K_R_P2   
 endmodule functions