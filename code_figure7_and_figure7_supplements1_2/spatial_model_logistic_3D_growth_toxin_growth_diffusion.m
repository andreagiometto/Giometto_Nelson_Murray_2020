% Written by Andrea Giometto, MIT license
% This function runs the spatial model with growth-dependent diffusion
% coefficient for yeast cells, growth-dependent toxin production rate and 3D toxin
% diffusion. This code was used to generate the model simulations reported
% in Fig. 7 supplement 1E. The initial conditions are idealized and do not
% reflect experimental intial conditions

function [] = spatial_model_logistic_3D_growth_toxin_growth_diffusion()
    % radius of the inoculum (cm)
    r0=0.5;
    tic
    
    R=1.5; % cm (radial size of the domain)
    h=0.36; % cm (height of the agar)
    Dt=0.003; % cm^2/h (toxin diffusion coefficient)
    g=0.36; % 1/h (maximum growth rate)
    K=2.3e+08*h; % g/cm^2 (spatial carrying capacity, note the multiplication by h)
    a1=1.57e-09;  % mL/(cells h^2) (toxin production coefficient)
    a2=4.75e-10;
    b=2.5e-09; % mL/(cells h) (toxin attachment coefficient)
    Dy=3.1460e-07; % cells diffusion coefficient cm^2/h
    Dfactor=1; % scaling factor for yeast diffusion coefficient

    K0=100000; % initial reference cell density cells/mL    

    transfers=6; % number of transfers
    tfin=48*transfers; % end time of the experiments (48 h for each growth period)
    tsample=zeros(1+transfers,1); % times at which to record population profile
    for i=1:transfers
        tsample(i+1,1)=48*i;
    end
    dilution=1;
    dilutionF=10000; % dilution factor during replica plating

    dt=0.0001; % temporal integration step
    dr=0.0025; % radial integration step
    dz=0.005; % altitudinal integration step

    N=round(R/dr); % number of points in the radial direction
    M=round(h/dz); % number of points in the altitudinal direction
    T=ceil(tfin/dt); % number of time points

    n1=zeros(N,1); % population profile of strain K1
    n2=zeros(N,1); % population profile of strain K2
    c1=zeros(M,N); % concentration profile of toxin c1
    c2=zeros(M,N); % concentration profile of toxin c2

    r=linspace(0,R,N);
    z=linspace(0,h,M);

    % variable to store profile through time
    sample=1;
    c1sample=zeros(M,N,length(tsample));
    c2sample=zeros(M,N,length(tsample));
    n1sample=zeros(N,length(tsample));
    n2sample=zeros(N,length(tsample));

    % set initial condition
    for j=1:round(r0/dr)
        n1(j,1)=K0/2;
    end
    for j=round(r0/dr)+1:N
        n2(j,1)=K0/2;
    end
    n1T=n1;
    n2T=n2;
    c1T=c1;
    c2T=c2;

    rinvmatrix=zeros(M,N);
    for i=1:M
        for j=2:N
            rinvmatrix(i,j)=1/r(j);
        end
    end

    percent=1;
    for t=1:T
        % display time advancement
        if t/T>=percent/10
            percent=percent+1;
            t/T*10
        end

        % concentrations update
        c1(2:M-1,2:N-1)=c1T(2:M-1,2:N-1)+dt*Dt*(rinvmatrix(2:M-1,2:N-1).*(c1T(2:M-1,3:N)-c1T(2:M-1,1:N-2))/(2*dr)+(c1T(2:M-1,3:N)+c1T(2:M-1,1:N-2)-2*c1T(2:M-1,2:N-1))/(dr^2)+(c1T(1:M-2,2:N-1)+c1T(3:M,2:N-1)-2*c1T(2:M-1,2:N-1))/(dz^2));
        c2(2:M-1,2:N-1)=c2T(2:M-1,2:N-1)+dt*Dt*(rinvmatrix(2:M-1,2:N-1).*(c2T(2:M-1,3:N)-c2T(2:M-1,1:N-2))/(2*dr)+(c2T(2:M-1,3:N)+c2T(2:M-1,1:N-2)-2*c2T(2:M-1,2:N-1))/(dr^2)+(c2T(1:M-2,2:N-1)+c2T(3:M,2:N-1)-2*c2T(2:M-1,2:N-1))/(dz^2));

        for j=2:(N-1) % top - bottom boundary

            c1(1,j)=c1(2,j); % bottom
            % growth dependent toxin production rate
            c1(M,j)=(c1(M-1,j)+a1*dz/Dt*n1(j,1))*(1-(n1T(j,1)+n2T(j,1))/K)/(1+b*dz/Dt*(n1(j,1)+n2(j,1)));
            c2(1,j)=c2(2,j); % bottom
            % growth dependent toxin production rate
            c2(M,j)=(c2(M-1,j)+a2*dz/Dt*n2(j,1))*(1-(n1T(j,1)+n2T(j,1))/K)/(1+b*dz/Dt*(n1(j,1)+n2(j,1)));
        end
        for i=1:M % left - right boundary: no flux 
            c1(i,1)=c1(i,2);
            c1(i,N)=c1(i,N-1);
            c2(i,1)=c2(i,2);
            c2(i,N)=c2(i,N-1);
        end
        for j=2:N-1 % n1, n2 update
            n1(j,1)=n1T(j,1)+dt*g*n1T(j,1)*(1-(n1T(j,1)+n2T(j,1))/K)+dt*Dfactor*Dy*(...
            1/r(j)*(n1T(j+1,1)-n1T(j-1,1))/(2*dr)*(1-(n1T(j,1)+n2T(j,1))/K)...
            +2*((1-(n1T(j+1,1)+n2T(j+1,1))/K)-(1-(n1T(j-1,1)+n2T(j-1,1))/K))/(2*dr)*(n1T(j+1,1)-n1T(j-1,1))/(2*dr)...
            +(1-(n1T(j,1)+n2T(j,1))/K)*(n1T(j+1,1)-2*n1T(j,1)+n1T(j-1,1))/(dr^2)...
            +1/r(j)*n1T(j,1)*((1-(n1T(j+1,1)+n2T(j+1,1))/K)-(1-(n1T(j-1,1)+n2T(j-1,1))/K))/(2*dr)...
            +n1T(j,1)*(1-(n1T(j+1,1)+n2T(j+1,1))/K+1-(n1T(j-1,1)+n2T(j-1,1))/K-2*(1-(n1T(j,1)+n2T(j,1))/K))/(dr^2)...
            )-dt*n1T(j,1)*c2T(M,j);
            n2(j,1)=n2T(j,1)+dt*g*n2T(j,1)*(1-(n1T(j,1)+n2T(j,1))/K)+dt*Dfactor*Dy*(...
            1/r(j)*(n2T(j+1,1)-n2T(j-1,1))/(2*dr)*(1-(n1T(j,1)+n2T(j,1))/K)...
            +2*((1-(n1T(j+1,1)+n2T(j+1,1))/K)-(1-(n1T(j-1,1)+n2T(j-1,1))/K))/(2*dr)*(n2T(j+1,1)-n2T(j-1,1))/(2*dr)...
            +(1-(n1T(j,1)+n2T(j,1))/K)*(n2T(j+1,1)-2*n2T(j,1)+n2T(j-1,1))/(dr^2)...
            +1/r(j)*n2T(j,1)*((1-(n1T(j+1,1)+n2T(j+1,1))/K)-(1-(n1T(j-1,1)+n2T(j-1,1))/K))/(2*dr)...
            +n2T(j,1)*(1-(n1T(j+1,1)+n2T(j+1,1))/K+1-(n1T(j-1,1)+n2T(j-1,1))/K-2*(1-(n1T(j,1)+n2T(j,1))/K))/(dr^2)...
            )-dt*n2T(j,1)*c1T(M,j);
        end
        n1(1,1)=n1(2,1);
        n1(N,1)=n1(N-1,1);
        n1T=n1;
        n2(1,1)=n2(2,1);
        n2(N,1)=n2(N-1,1);
        n2T=n2;
        c1T=c1;
        c2T=c2;

        % store profiles at sampling time points
        if t*dt>=tsample(sample)
            c1sample(:,:,sample)=c1;
            c2sample(:,:,sample)=c2;
            n1sample(:,sample)=n1;
            n2sample(:,sample)=n2;
            sample=sample+1;
        end

        % perform dilution
        if t*dt>dilution*48
            c1=zeros(M,N);
            c2=zeros(M,N);
            c1T=c1;
            c2T=c2;
            n1=n1/dilutionF;
            n2=n2/dilutionF;
            n1T=n1;
            n2T=n2;
            dilution=dilution+1;
        end
    end
    toc

    % save simulations data
    save(['data_models/model_output_growth_logistic_propto_grate_Dgrowth_ideal_inoculum_diffusion_Dfactor_',num2str(Dfactor),'_Dt_',num2str(Dt),'_K_',num2str(K),'_a1_',num2str(a1),'_a2_',num2str(a2),'_transfers_',num2str(transfers),'_r0_',num2str(r0),'_dilutionF_',num2str(dilutionF),'_R_',num2str(R),'_dr_',num2str(dr),'_dz_',num2str(dz),'_dt_',num2str(dt),'_n0IC-n0IM_',num2str(K0/2),'.mat'])

end