MODULE Element_module
        real(8):: W1(5),G1(5)
        real(8):: N(5,5),dN(5,5)!,Me2(30,30),Me1(15,15,5)!Bk(30,9,11),
        real(8):: Ele_length(11)
        


    end
    
Module body_module
        !real(8),ALLOCATable:: E_body(:),Rou_body(:),RouA_body(:)      
        !real(8),ALLOCATable:: Iz_body(:,:,:),length_body(:),Cm_body(:,:,:)
        Type Body_type
            integer:: Con_num,Ele_num,Node_num,F_num,Contact_node_num,F_U_num,F_W_num;
            real(8):: E,Rou,RouA,S;
            real(8):: Iz(3,3),length,Cm(3,3);
            real(8):: Le;
            real(8):: r0(3),Theta0(3),Quater(4);
            integer,allocatable:: Fn(:),Con_dof(:,:),Contact_node(:),F_u(:),F_W(:);
            real(8),ALLOCATable:: r(:,:),dr(:,:),ddr(:,:),omega(:,:),domega(:,:),gama_w(:),Phi_C(:,:),Phi_C3(:,:),Phi_C_P(:,:,:)!,Phi_C_P1(:,:,:),Phi_C_P2(:,:,:),Phi_C_P3(:,:,:)!
            real(8),ALLOCATable:: Phi_C_P0(:,:)
            real(8),ALLOCATable:: a(:),da(:),dda(:)!
            real(8),Allocatable:: y(:,:),dy(:,:),ddy(:,:)!
            real(8),Allocatable:: Mp(:)
            real(8),allocatable:: y_temp(:),dy_temp(:),ddy_temp(:),r_temp(:),dr_temp(:),ddr_temp(:),omega_temp(:),domega_temp(:),r_temp1(:,:),dr_temp1(:,:),ddr_temp1(:,:),omega_temp1(:,:),domega_temp1(:,:),y_temp1(:,:),dy_temp1(:,:),ddy_temp1(:,:);
            real(8),allocatable:: Roup_temp(:),dRoup_temp(:),Phi_temp(:,:),Ki_temp(:,:),AI_temp(:,:),omega_r_temp(:)
            real(8),allocatable:: Z(:,:),z1(:),G(:,:),GG(:,:),F_contact(:),z1_J(:,:),H_P(:,:),g0(:,:);
            real(8),allocatable:: Phi_quater(:)
            real(8),allocatable:: Ai(:,:),Array(:,:),Array_P(:,:,:),Array_PP(:,:,:,:),Array_G(:,:,:)
            real(8):: A1(3,3)
        endtype
        type(Body_type),allocatable:: Body(:)
    end
    
    MODULE Model_module
        INTEGER:: N_end,Mn,Body_num,Node_num_All,ele_num_all
        INTEGER:: Con_num_all,F_num_all,Con_num_Act
        !real(8),allocatable:: Le_body(:)
        real(8):: a0(3),r1(3),dr1(3),omega1(3),Gra(3);
        real(8):: A1(3,3);
        INTEGER,ALLOCATable:: FN(:),Con_list(:,:),Act_Cn(:)
        !real(8),ALLOCATable:: Gex(:),Le(:)!
        !integer,Allocatable:: Con_num_body(:),con_dof_body(:,:),ele_num_body(:),Node_num_body(:);
        !real(8),ALLOCATable:: r(:,:),dr(:,:),omega(:,:)!
        !real(8),ALLOCATable:: a(:),da(:),dda(:)!
        !real(8),Allocatable:: y(:,:),dy(:,:),ddy(:,:)!
        real(8),Allocatable:: time1(:),P(:),Mp(:),Con_Value(:);!
        !equation中的参数
        real(8),Allocatable:: X(:,:),dx(:,:),lumda(:,:),lumda2(:,:),Y(:,:),dY(:,:),DDY(:,:),lumda_quater(:,:)!,Z(:,:),z1(:)
        !recur0 中的参数
!        real(8),Allocatable:: Roup(:,:),dRoup(:,:),F1(:),beta(:,:)
!        real(8),Allocatable:: T(:,:)
    end
    module integration_parameter
        real(8):: Rad,alpha_f,alpha_m,beta,gama,del_t,del_t0,gama_1,beta_1;
            !ck=1-alpha_f;
            !c0=(1-alpha_m)/beta/delta_t^2;
            !c1=ck*gama/beta/delta_t;
            !c2=delta_t*c0;
            !c3=c2*delta_t/2-1;
            !c4=ck*gama/beta-1;
            !c5=ck*(gama/2/beta-1)*delta_t;        
    end