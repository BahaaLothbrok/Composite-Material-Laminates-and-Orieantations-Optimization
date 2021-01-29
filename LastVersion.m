% MATLAB CODE FOR OPTIMIZATION OF LAYERS

clc
close all
clear vars
clear all

%% Givens
k=0; %No. Layers
e_1=0;   %Equal to one when k=0;
sym=1;
rtheta=[-45 0 45 90];  %Input
% theta=[45 0 45];
theta1=rtheta;
theta2=rtheta;
theta3=rtheta;
theta4=rtheta;
theta5=rtheta;
theta6=rtheta;
theta7=rtheta;
theta8=rtheta;
v12=0.28;
E1=2e7;
E2=3e6;
G12=8e5;
% F=[30e3 20e3 15e3 40000 20000 10000];
F=[0 20e3 15e3 4000 2000 0];

%% For Loops [Proceeding]
if k==0
    k=1;
    e_1=1;
end

%% One Layer
if k==1
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    
    t=1;
    Av_RV=0;
    for x_1=1:length(theta1)
        for i = 1:k
            thetaa=theta1(x_1) ;
            [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
            A          = A + Q_bar*tlamina(i);
            B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
            D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
        end
        
        C = [A B; B D]  ;
        
        Eps = F/C;
        Eps0 = Eps(1:3)       ;    % midplane strain
        k0   = Eps(4:6)       ;   % midplane curvature
        Epsk_loc = zeros(3,N);
        S_loc = zeros(3,N);
        A = zeros(3);
        B = zeros(3);
        D = zeros(3);
        
        %         disp('For theta')
        %         disp(thetaa)
        %         disp(' *** (3) Tsai-Hill Stress Criterion ***')
        THV=zeros(1,k);
        for ic = 1:N
            % each lamina strain in laminate coordinates
            Epsk = Eps0 + Z(ic)*k0;
            % lmainate transformation matrix
            m = cosd(thetaa(ic));  n = sind(thetaa(ic));
            T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
            % each lamina strain in local coordinates
            Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
            % each lamina stress in local coordinates
            S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
            %     %______________________________________________________________________
            %    elseif CN==4
            Xt=15e4;%str2double(get(handles.Xt,'string'));
            Yt=5e4;%str2double(get(handles.Yt,'string'));
            Xc=-Xt;%str2double(get(handles.Xc,'string'));
            Yc=-Yt;%str2double(get(handles.Yc,'string'));
            S=7e4;%str2double(get(handles.S,'string'));
            
            % Tsai-Hill Stress Criterion
            
            
            if S_loc(1) >= 0
                X = Xt;
            elseif S_loc(1) <= 0
                X = Xc;
            end
            
            if S_loc(2) >= 0
                Y = Yt;
            elseif S_loc(2) <= 0
                Y = Yc;
            end
            
            TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
            THV(ic)=TH;
            if TH > 1
                %fprintf(2,' Failed Lamina ')
                %disp(['( ' num2str(ic) ' )'])
                test3(1,ic)={' failed'};
                vdd3(ic)=0;
            else
                %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                
                test3(1,ic)={' is safe'};
                vdd3(ic)=1;
                
            end
            %disp(TH)
        end
        if THV(1)<1
            R_total=sum(THV);
            
            Av_RV(t)=R_total/k; %Vector of Average R
            theta_AV(:,t)=[thetaa(1)];
            t=t+1;
            
            
        end
    end
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=2;
    end
end

%% Two Layers
if k==2
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    t=1;
    Av_RV=0;
    for x_2=1:length(theta2)
        for x_1=1:length(theta1)
            for i = 1:k
                thetaa=[theta2(x_2) theta1(x_1)] ;
                [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                A          = A + Q_bar*tlamina(i);
                B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
            end
            
            C = [A B; B D]  ;
            
            Eps = F/C;
            Eps0 = Eps(1:3)       ;    % midplane strain
            k0   = Eps(4:6)       ;   % midplane curvature
            Epsk_loc = zeros(3,N);
            S_loc = zeros(3,N);
            A = zeros(3);
            B = zeros(3);
            D = zeros(3);
            
            %disp('For theta')
            %disp(thetaa)
            %disp(' *** (3) Tsai-Hill Stress Criterion ***')
            THV=zeros(1,k);
            for ic = 1:N
                % each lamina strain in laminate coordinates
                Epsk = Eps0 + Z(ic)*k0;
                % lmainate transformation matrix
                m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                % each lamina strain in local coordinates
                Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                % each lamina stress in local coordinates
                S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                %     %______________________________________________________________________
                %    elseif CN==4
                Xt=15e4;%str2double(get(handles.Xt,'string'));
                Yt=5e4;%str2double(get(handles.Yt,'string'));
                Xc=-Xt;%str2double(get(handles.Xc,'string'));
                Yc=-Yt;%str2double(get(handles.Yc,'string'));
                S=7e4;%str2double(get(handles.S,'string'));
                
                % Tsai-Hill Stress Criterion
                
                
                if S_loc(1) > 0
                    X = Xt;
                elseif S_loc(1) < 0
                    X = Xc;
                end
                
                if S_loc(2) > 0
                    Y = Yt;
                elseif S_loc(2) < 0
                    Y = Yc;
                end
                
                TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                THV(ic)=TH;
                if TH > 1
                    %fprintf(2,' Failed Lamina ')
                    %disp(['( ' num2str(ic) ' )'])
                    test3(1,ic)={' failed'};
                    vdd3(ic)=0;
                else
                    %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                    
                    test3(1,ic)={' is safe'};
                    vdd3(ic)=1;
                    
                end
                %disp(TH)
            end
            if THV(1)<1 && THV(2)<1
                R_total=sum(THV);
                
                Av_RV(t)=R_total/k; %Vector of Average R
                theta_AV(:,t)=[thetaa(1) thetaa(2)];
                t=t+1;
                
                
            end
            
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=3;
    end
end

%% Three Layers
if k==3
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    
    t=1;
    Av_RV=0;
    for x_3=1:length(theta3)
        for x_2=1:length(theta2)
            for x_1=1:length(theta1)
                for i = 1:k
                    thetaa=[theta3(x_3) theta2(x_2) theta1(x_1)] ;
                    [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                    A          = A + Q_bar*tlamina(i);
                    B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                    D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                end
                
                C = [A B; B D]  ;
                
                Eps = F/C;
                Eps0 = Eps(1:3)       ;    % midplane strain
                k0   = Eps(4:6)       ;   % midplane curvature
                Epsk_loc = zeros(3,N);
                S_loc = zeros(3,N);
                A = zeros(3);
                B = zeros(3);
                D = zeros(3);
                
                
                %disp('For theta')
                %disp(thetaa)
                %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                THV=zeros(1,k);    %Vector of TH of each lamina
                for ic = 1:N
                    % each lamina strain in laminate coordinates
                    Epsk = Eps0 + Z(ic)*k0;
                    % lmainate transformation matrix
                    m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                    T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                    % each lamina strain in local coordinates
                    Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                    % each lamina stress in local coordinates
                    S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                    %     %______________________________________________________________________
                    %    elseif CN==4
                    Xt=15e4;%str2double(get(handles.Xt,'string'));
                    Yt=5e4;%str2double(get(handles.Yt,'string'));
                    Xc=-Xt;%str2double(get(handles.Xc,'string'));
                    Yc=-Yt;%str2double(get(handles.Yc,'string'));
                    S=7e4;%str2double(get(handles.S,'string'));
                    
                    % Tsai-Hill Stress Criterion
                    
                    
                    if S_loc(1) > 0
                        X = Xt;
                    elseif S_loc(1) < 0
                        X = Xc;
                    end
                    
                    if S_loc(2) > 0
                        Y = Yt;
                    elseif S_loc(2) < 0
                        Y = Yc;
                    end
                    
                    TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                    THV(ic)=TH;
                    if TH > 1
                        %fprintf(2,' Failed Lamina ')
                        %disp(['( ' num2str(ic) ' )'])
                        test3(1,ic)={' failed'};
                        vdd3(ic)=0;
                    else
                        %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                        
                        test3(1,ic)={' is safe'};
                        vdd3(ic)=1;
                        
                    end
                    %disp(TH)
                end
                
                if THV(1)<1 && THV(2)<1 && THV(3)<1
                    R_total=sum(THV);
                    
                    Av_RV(t)=R_total/k; %Vector of Average R
                    theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3)];
                    %                     disp('For theta')
                    %                     disp(theta_AV(:,t)')
                    %                     disp('Average R')
                    %                     disp(Av_RV(t))
                    %                     disp('Laminate is Safe')
                    %                     disp('----------------------------------------------------')
                    
                    %theta_2(t)=thetaa(2);
                    %theta_3(t)=thetaa(3);
                    t=t+1;
                    
                    
                end
                
                
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=4;
    end
    
end

%% Four Layers
if k==4
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    t=1;
    Av_RV=0;
    for x_4 = 1:length(theta4)
        for x_3=1:length(theta3)
            for x_2=1:length(theta2)
                for x_1=1:length(theta1)
                    for i = 1:k
                        thetaa=[theta4(x_4) theta3(x_3) theta2(x_2) theta1(x_1)] ;
                        [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                        A          = A + Q_bar*tlamina(i);
                        B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                        D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                    end
                    
                    C = [A B; B D]  ;
                    
                    Eps = F/C;
                    Eps0 = Eps(1:3)       ;    % midplane strain
                    k0   = Eps(4:6)       ;   % midplane curvature
                    Epsk_loc = zeros(3,N);
                    S_loc = zeros(3,N);
                    A = zeros(3);
                    B = zeros(3);
                    D = zeros(3);
                    
                    %disp('For theta')
                    %disp(thetaa)
                    %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                    THV=zeros(1,k);    %Vector of TH of each lamina
                    for ic = 1:N
                        % each lamina strain in laminate coordinates
                        Epsk = Eps0 + Z(ic)*k0;
                        % lmainate transformation matrix
                        m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                        T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                        % each lamina strain in local coordinates
                        Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                        % each lamina stress in local coordinates
                        S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                        %     %______________________________________________________________________
                        %    elseif CN==4
                        Xt=15e4;%str2double(get(handles.Xt,'string'));
                        Yt=5e4;%str2double(get(handles.Yt,'string'));
                        Xc=-Xt;%str2double(get(handles.Xc,'string'));
                        Yc=-Yt;%str2double(get(handles.Yc,'string'));
                        S=7e4;%str2double(get(handles.S,'string'));
                        
                        % Tsai-Hill Stress Criterion
                        
                        
                        if S_loc(1) > 0
                            X = Xt;
                        elseif S_loc(1) < 0
                            X = Xc;
                        end
                        
                        if S_loc(2) > 0
                            Y = Yt;
                        elseif S_loc(2) < 0
                            Y = Yc;
                        end
                        
                        TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                        THV(ic)=TH;
                        if TH > 1
                            %fprintf(2,' Failed Lamina ')
                            %disp(['( ' num2str(ic) ' )'])
                            test3(1,ic)={' failed'};
                            vdd3(ic)=0;
                        else
                            %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                            
                            test3(1,ic)={' is safe'};
                            vdd3(ic)=1;
                            
                        end
                        %disp(TH)
                    end
                    if THV(1)<1 && THV(2)<1 && THV(3)<1 && THV(4)<1
                        R_total=sum(THV);
                        
                        Av_RV(t)=R_total/k; %Vector of Average R
                        theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3) thetaa(4)];
                        t=t+1;
                        
                        
                    end
                    
                end
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=5;
    end
    
end

%% Five Layers
if k==5
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    t=1;
    Av_RV=0;
    for x_5 = 1:length(theta5)
        for x_4 = 1:length(theta4)
            for x_3=1:length(theta3)
                for x_2=1:length(theta2)
                    for x_1=1:length(theta1)
                        for i = 1:k
                            thetaa=[theta5(x_5) theta4(x_4) theta3(x_3) theta2(x_2) theta1(x_1)] ;
                            [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                            A          = A + Q_bar*tlamina(i);
                            B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                            D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                        end
                        
                        C = [A B; B D]  ;
                        
                        Eps = F/C;
                        Eps0 = Eps(1:3)       ;    % midplane strain
                        k0   = Eps(4:6)       ;   % midplane curvature
                        Epsk_loc = zeros(3,N);
                        S_loc = zeros(3,N);
                        A = zeros(3);
                        B = zeros(3);
                        D = zeros(3);
                        
                        %disp('For theta')
                        %disp(thetaa)
                        %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                        THV=zeros(1,k);    %Vector of TH of each lamina
                        for ic = 1:N
                            % each lamina strain in laminate coordinates
                            Epsk = Eps0 + Z(ic)*k0;
                            % lmainate transformation matrix
                            m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                            T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                            % each lamina strain in local coordinates
                            Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                            % each lamina stress in local coordinates
                            S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                            %     %______________________________________________________________________
                            %    elseif CN==4
                            Xt=15e4;%str2double(get(handles.Xt,'string'));
                            Yt=5e4;%str2double(get(handles.Yt,'string'));
                            Xc=-Xt;%str2double(get(handles.Xc,'string'));
                            Yc=-Yt;%str2double(get(handles.Yc,'string'));
                            S=7e4;%str2double(get(handles.S,'string'));
                            
                            % Tsai-Hill Stress Criterion
                            
                            
                            if S_loc(1) > 0
                                X = Xt;
                            elseif S_loc(1) < 0
                                X = Xc;
                            end
                            
                            if S_loc(2) > 0
                                Y = Yt;
                            elseif S_loc(2) < 0
                                Y = Yc;
                            end
                            
                            TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                            THV(ic)=TH;
                            if TH > 1
                                %fprintf(2,' Failed Lamina ')
                                %disp(['( ' num2str(ic) ' )'])
                                test3(1,ic)={' failed'};
                                vdd3(ic)=0;
                            else
                                %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                                
                                test3(1,ic)={' is safe'};
                                vdd3(ic)=1;
                                
                            end
                            %disp(TH)
                        end
                        if THV(1)<1 && THV(2)<1 && THV(3)<1 && THV(4)<1&& THV(5)<1
                            R_total=sum(THV);
                            
                            Av_RV(t)=R_total/k; %Vector of Average R
                            theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3) thetaa(4) thetaa(5)];
                            t=t+1;
                            
                            
                        end
                        
                    end
                end
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=6;
    end
    
end

%% Six Layes
if k==6
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    t=1;
    Av_RV=0;
    for x_6 = 1:length(theta6)
        for x_5 = 1:length(theta5)
            for x_4 = 1:length(theta4)
                for x_3=1:length(theta3)
                    for x_2=1:length(theta2)
                        for x_1=1:length(theta1)
                            for i = 1:k
                                thetaa=[theta6(x_6) theta5(x_5) theta4(x_4) theta3(x_3) theta2(x_2) theta1(x_1)] ;
                                [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                                A          = A + Q_bar*tlamina(i);
                                B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                                D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                            end
                            
                            C = [A B; B D]  ;
                            
                            Eps = F/C;
                            Eps0 = Eps(1:3)       ;    % midplane strain
                            k0   = Eps(4:6)       ;   % midplane curvature
                            Epsk_loc = zeros(3,N);
                            S_loc = zeros(3,N);
                            A = zeros(3);
                            B = zeros(3);
                            D = zeros(3);
                            
                            %disp('For theta')
                            %disp(thetaa)
                            %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                            THV=zeros(1,k);    %Vector of TH of each lamina
                            for ic = 1:N
                                % each lamina strain in laminate coordinates
                                Epsk = Eps0 + Z(ic)*k0;
                                % lmainate transformation matrix
                                m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                                T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                                % each lamina strain in local coordinates
                                Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                                % each lamina stress in local coordinates
                                S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                                %     %______________________________________________________________________
                                %    elseif CN==4
                                Xt=15e4;%str2double(get(handles.Xt,'string'));
                                Yt=5e4;%str2double(get(handles.Yt,'string'));
                                Xc=-Xt;%str2double(get(handles.Xc,'string'));
                                Yc=-Yt;%str2double(get(handles.Yc,'string'));
                                S=7e4;%str2double(get(handles.S,'string'));
                                
                                % Tsai-Hill Stress Criterion
                                
                                
                                if S_loc(1) > 0
                                    X = Xt;
                                elseif S_loc(1) < 0
                                    X = Xc;
                                end
                                
                                if S_loc(2) > 0
                                    Y = Yt;
                                elseif S_loc(2) < 0
                                    Y = Yc;
                                end
                                
                                TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                                THV(ic)=TH;
                                
                                if TH > 1
                                    %fprintf(2,' Failed Lamina ')
                                    %disp(['( ' num2str(ic) ' )'])
                                    test3(1,ic)={' failed'};
                                    vdd3(ic)=0;
                                else
                                    %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                                    
                                    test3(1,ic)={' is safe'};
                                    vdd3(ic)=1;
                                    
                                end
                                %disp(TH)
                            end
                            if THV(1)<1 && THV(2)<1 && THV(3)<1 ...
                                    && THV(4)<1&& THV(5)<1&& THV(6)<1
                                R_total=sum(THV);
                                
                                Av_RV(t)=R_total/k; %Vector of Average R
                                theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3)...
                                    thetaa(4) thetaa(5) thetaa(6)];
                                t=t+1;
                                
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=7;
    end
    
end

%% Seven Layers
if k==7
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    
    t=1;
    Av_RV=0;
    for x_7 = 1:length(theta7)
        for x_6 = 1:length(theta6)
            for x_5 = 1:length(theta5)
                for x_4 = 1:length(theta4)
                    for x_3=1:length(theta3)
                        for x_2=1:length(theta2)
                            for x_1=1:length(theta1)
                                for i = 1:k
                                    thetaa=[theta7(x_7) theta6(x_6) theta5(x_5) theta4(x_4) theta3(x_3) theta2(x_2) theta1(x_1)] ;
                                    [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                                    A          = A + Q_bar*tlamina(i);
                                    B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                                    D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                                end
                                
                                C = [A B; B D]  ;
                                
                                Eps = F/C;
                                Eps0 = Eps(1:3)       ;    % midplane strain
                                k0   = Eps(4:6)       ;   % midplane curvature
                                Epsk_loc = zeros(3,N);
                                S_loc = zeros(3,N);
                                A = zeros(3);
                                B = zeros(3);
                                D = zeros(3);
                                
                                %disp('For theta')
                                %disp(thetaa)
                                %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                                THV=zeros(1,k);
                                for ic = 1:N
                                    % each lamina strain in laminate coordinates
                                    Epsk = Eps0 + Z(ic)*k0;
                                    % lmainate transformation matrix
                                    m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                                    T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                                    % each lamina strain in local coordinates
                                    Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                                    % each lamina stress in local coordinates
                                    S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                                    %     %______________________________________________________________________
                                    %    elseif CN==4
                                    Xt=15e4;%str2double(get(handles.Xt,'string'));
                                    Yt=5e4;%str2double(get(handles.Yt,'string'));
                                    Xc=-Xt;%str2double(get(handles.Xc,'string'));
                                    Yc=-Yt;%str2double(get(handles.Yc,'string'));
                                    S=7e4;%str2double(get(handles.S,'string'));
                                    
                                    % Tsai-Hill Stress Criterion
                                    
                                    
                                    if S_loc(1) > 0
                                        X = Xt;
                                    elseif S_loc(1) < 0
                                        X = Xc;
                                    end
                                    
                                    if S_loc(2) > 0
                                        Y = Yt;
                                    elseif S_loc(2) < 0
                                        Y = Yc;
                                    end
                                    
                                    TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                                    THV(ic)=TH;
                                    if TH > 1
                                        %fprintf(2,' Failed Lamina ')
                                        %disp(['( ' num2str(ic) ' )'])
                                        test3(1,ic)={' failed'};
                                        vdd3(ic)=0;
                                    else
                                        %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                                        
                                        test3(1,ic)={' is safe'};
                                        vdd3(ic)=1;
                                        
                                    end
                                    %disp(TH)
                                end
                                
                                if THV(1)<1 && THV(2)<1 && THV(3)<1 ...
                                        && THV(4)<1&& THV(5)<1&& THV(6)<1 ...
                                        && THV(7)<1
                                    R_total=sum(THV);
                                    
                                    Av_RV(t)=R_total/k; %Vector of Average R
                                    theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3)...
                                        thetaa(4) thetaa(5) thetaa(6) thetaa(7)];
                                    t=t+1;
                                    
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates Failed')
        disp('-----------------------------------------------------------')
        k=8;
    end
    
    
end

%% Eight Layers
if k==8
    
    tlamina=0.1*ones(1,k);
    N      = k;
    thickness      = sum(tlamina);
    Nh   = N+1;
    h    = zeros(1, Nh);                 % lamina buttom coordinate from the bottom of the laminate
    Z    = zeros(1, N);                  % lamina coordinates from the midplane
    h(1) = -thickness/2;
    h(Nh)= thickness/2;
    Z(1) = -thickness/2 + tlamina(1)/2;
    
    A = zeros(3);
    B = zeros(3);
    D = zeros(3);
    
    for i = 2:N
        Z(i) = Z(i-1) + tlamina(i-1)/2 + tlamina(i)/2;
        h(i) = Z(i-1) + 0.5*tlamina(i-1);
    end
    
    t=1;
    Av_RV=0;
    for x_8 = 1:length(theta8)
        for x_7 = 1:length(theta7)
            for x_6 = 1:length(theta6)
                for x_5 = 1:length(theta5)
                    for x_4 = 1:length(theta4)
                        for x_3=1:length(theta3)
                            for x_2=1:length(theta2)
                                for x_1=1:length(theta1)
                                    for i = 1:k
                                        thetaa=[theta8(x_8) theta7(x_7) theta6(x_6) theta5(x_5) theta4(x_4) theta3(x_3) theta2(x_2) theta1(x_1)] ;
                                        [T_theta,Q_ij,Q_bar]=Lamina(thetaa(i),v12,E1,E2,G12);
                                        A          = A + Q_bar*tlamina(i);
                                        B          = B + (1/2)*Q_bar*(h(i+1)^2 - h(i)^2);
                                        D          = D + (1/3)*Q_bar*(h(i+1)^3 - h(i)^3);
                                    end
                                    
                                    C = [A B; B D]  ;
                                    
                                    Eps = F/C;
                                    Eps0 = Eps(1:3)       ;    % midplane strain
                                    k0   = Eps(4:6)       ;   % midplane curvature
                                    Epsk_loc = zeros(3,N);
                                    S_loc = zeros(3,N);
                                    A = zeros(3);
                                    B = zeros(3);
                                    D = zeros(3);
                                    
                                    %disp('For theta')
                                    %disp(thetaa)
                                    %disp(' *** (3) Tsai-Hill Stress Criterion ***')
                                    THV=zeros(1,k);
                                    for ic = 1:N
                                        % each lamina strain in laminate coordinates
                                        Epsk = Eps0 + Z(ic)*k0;
                                        % lmainate transformation matrix
                                        m = cosd(thetaa(ic));  n = sind(thetaa(ic));
                                        T = [m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
                                        % each lamina strain in local coordinates
                                        Epsk_loc = T*[Epsk(1); Epsk(2); 0.5*Epsk(3)];
                                        % each lamina stress in local coordinates
                                        S_loc    = Q_ij*[Epsk_loc(1); Epsk_loc(2); 2*Epsk_loc(3)];
                                        %     %______________________________________________________________________
                                        %    elseif CN==4
                                        Xt=15e4;%str2double(get(handles.Xt,'string'));
                                        Yt=5e4;%str2double(get(handles.Yt,'string'));
                                        Xc=-Xt;%str2double(get(handles.Xc,'string'));
                                        Yc=-Yt;%str2double(get(handles.Yc,'string'));
                                        S=7e4;%str2double(get(handles.S,'string'));
                                        
                                        % Tsai-Hill Stress Criterion
                                        
                                        
                                        if S_loc(1) > 0
                                            X = Xt;
                                        elseif S_loc(1) < 0
                                            X = Xc;
                                        end
                                        
                                        if S_loc(2) > 0
                                            Y = Yt;
                                        elseif S_loc(2) < 0
                                            Y = Yc;
                                        end
                                        
                                        TH = (S_loc(1)/X)^2 - S_loc(1)*S_loc(2)/X^2 ...
                                            + S_loc(2)^2/Y^2 + S_loc(3)^2/S^2;
                                        THV(ic)=TH;
                                        if TH > 1
                                            %fprintf(2,' Failed Lamina ')
                                            %disp(['( ' num2str(ic) ' )'])
                                            test3(1,ic)={' failed'};
                                            vdd3(ic)=0;
                                        else
                                            %disp(['lamina ' num2str(ic) ' is safe', 'r'])
                                            
                                            test3(1,ic)={' is safe'};
                                            vdd3(ic)=1;
                                            
                                        end
                                        %disp(TH)
                                    end
                                    
                                    if THV(1)<1 && THV(2)<1 && THV(3)<1 ...
                                            && THV(4)<1&& THV(5)<1&& THV(6)<1 ...
                                            && THV(7)<1 && THV(8)<1
                                        R_total=sum(THV);
                                        
                                        Av_RV(t)=R_total/k; %Vector of Average R
                                        theta_AV(:,t)=[thetaa(1) thetaa(2) thetaa(3)...
                                            thetaa(4) thetaa(5) thetaa(6) thetaa(7) ...
                                            thetaa(8)];
                                        t=t+1;
                                        
                                        
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    if e_1==1 & Av_RV==0
        disp(['For Number of layers ' num2str(k)])
        disp('All Laminates for all number of layers from 1 to 8 Failed')
        disp('-----------------------------------------------------------')
        k=3;
    end 
end

%% Displaying The Best Laminates
RTheta=zeros(k+1,t-1);
for i=1:t-1
    RTheta(1,i)=Av_RV(i);
    RTheta(2:k+1,i)=theta_AV(:,i);
end



e=0;
if Av_RV==0
    disp('All Laminated Failed')
else
    UAv_RV=unique(Av_RV);  %Unique values of Average R
    SAV_R=sort(UAv_RV);   %Sorted Average R
    disp('Sorted R Ascending')
    disp('-------------------------------------------------------------------')
    if sym==0
        for i=1:length(UAv_RV)
            if i==1
                disp('The best Laminate Orientation are')
            elseif i==2
                disp('The Second best Laminate Orientation are')
            elseif i==3
                disp('The Third best Laminate Orientation are')
            else
                disp(['The ' num2str(i) 'th best Laminate Orientation are '])
            end
            disp('Orientation')
            disp(RTheta(2:k+1,find(SAV_R(i)==RTheta(1,:)))')
            disp('Average R')
            disp(SAV_R(i))
            disp('----------------------------------------------------------------')
        end
    elseif sym==1
        
        for i=1:t-1
            
            if k==1
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==2 && RTheta(2,i)==RTheta(k+1,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==3 && RTheta(2,i)==RTheta(k+1,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
                
            elseif k==4 && RTheta(2,i)==RTheta(k+1,i) && ...
                    RTheta(3,i)==RTheta(k,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==5 && RTheta(2,i)==RTheta(k+1,i) && ...
                    RTheta(3,i)==RTheta(k,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==6 && RTheta(2,i)==RTheta(k+1,i) && ...
                    RTheta(3,i)==RTheta(k,i) && RTheta(4,i)==RTheta(k-1,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==7 && RTheta(2,i)==RTheta(k+1,i) && ...
                    RTheta(3,i)==RTheta(k,i) && RTheta(4,i)==RTheta(k-1,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            elseif k==8 && RTheta(2,i)==RTheta(k+1,i) && ...
                    RTheta(3,i)==RTheta(k,i)&& RTheta(4,i)==RTheta(k-1,i) ...
                    && RTheta(5,i)==RTheta(k-2,i)
                disp('Orientation')
                disp(RTheta(2:k+1,i)')
                disp('Average R')
                disp(RTheta(1,i))
                disp('----------------------------------------------------------------')
                e=e+1;
            end
        end
        
        
        
    end
end


if e==0  && sym==1
    disp('No Safe Symmetric Laminates')
end
