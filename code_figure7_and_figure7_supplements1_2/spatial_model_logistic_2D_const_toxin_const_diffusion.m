% Written by Andrea Giometto, MIT license
% This function runs the spatial model with growth-dependent diffusion
% coefficient for yeast cells, constant toxin production rate and 2D toxin
% diffusion. This code was used to generate the model simulations reported
% in Fig. 7 supplement 1B. The initial conditions are idealized and do not
% reflect experimental intial conditions

function [] = spatial_model_logistic_2D_const_toxin_const_diffusion()
    % radius of the inoculum (cm)
    r0=0.5;
    tic
    
    R=1.5; % cm (radial size of the domain)
    h=0.36; % cm (height of the agar --- only used here to compute K)
    Dt=0.003; % cm^2/h (toxin diffusion coefficient)
    g=0.36; % 1/h (maximum growth rate)
    K=2.3e+08*h; % g/cm^2 (spatial carrying capacity, note the multiplication by h)
    a1=2.36e-09;  % mL/(cells h^2) (toxin production coefficient)
    a2=9.60e-10;
    b=1.2e-08; % mL/(cells h) (toxin attachment coefficient)
    Dy=3.1460e-07; % cells diffusion coefficient cm^2/h
    Dfactor=1; % scaling factor for yeast diffusion coefficient

    K0=100000; % initial reference cell density cells/mL    

    transfers=6; % number of transfers
%     transfers=1+13; % number of transfers
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
    T=ceil(tfin/dt); % number of time points

    n1=zeros(N,1); % population profile of strain K1
    n2=zeros(N,1); % population profile of strain K2
    c1=zeros(N,1); % concentration profile of toxin c1
    c2=zeros(N,1); % concentration profile of toxin c2

    r=linspace(0,R,N);

    % variable to store profile through time
    sample=1;
    c1sample=zeros(N,length(tsample));
    c2sample=zeros(N,length(tsample));
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

    percent=1;
    for t=1:T
        % display time advancement
        if t/T>=percent/10
            percent=percent+1;
            t/T*10
        end

        % left - right boundary: no flux 
        c1(1,1)=c1(2,1);
        c1(N,1)=c1(N-1,1);
        c2(1,1)=c2(2,1);
        c2(N,1)=c2(N-1,1);
        
        for j=2:N-1
            % n1, n2 update
            n1(j,1)=n1T(j,1)+dt*g*n1T(j,1)*(1-(n1T(j,1)+n2T(j,1))/K)+dt*Dfactor*Dy*(1/r(j)*(n1T(j+1,1)-n1T(j-1,1))/(2*dr)+(n1T(j+1,1)-2*n1T(j,1)+n1T(j-1,1))/(dr^2))-dt*n1T(j,1)*c2T(j,1);
            n2(j,1)=n2T(j,1)+dt*g*n2T(j,1)*(1-(n1T(j,1)+n2T(j,1))/K)+dt*Dfactor*Dy*(1/r(j)*(n2T(j+1,1)-n2T(j-1,1))/(2*dr)+(n2T(j+1,1)-2*n2T(j,1)+n2T(j-1,1))/(dr^2))-dt*n2T(j,1)*c1T(j,1);
            % concentrations update
            c1(j,1)=c1T(j,1)+dt*a1*n1T(j,1)-b*(n1T(j,1)+n2T(j,1)).*c1T(j,1)+dt*Dt*(1/r(j)*(c1T(j+1,1)-c1T(j,1))/(2*dr)+(c1T(j+1,1)+c1T(j-1,1)-2*c1T(j,1))/(dr^2));
            c2(j,1)=c2T(j,1)+dt*a2*n2T(j,1)-b*(n1T(j,1)+n2T(j,1)).*c2T(j,1)+dt*Dt*(1/r(j)*(c2T(j+1,1)-c2T(j,1))/(2*dr)+(c2T(j+1,1)+c2T(j-1,1)-2*c2T(j,1))/(dr^2));
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
            c1sample(:,sample)=c1;
            c2sample(:,sample)=c2;
            n1sample(:,sample)=n1;
            n2sample(:,sample)=n2;
            sample=sample+1;
        end

        % perform dilution
        if t*dt>dilution*48
            c1=zeros(N,1);
            c2=zeros(N,1);
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
    save(['data_models/model_output_growth_logistic2D_ideal_inoculum_diffusion_Dfactor_',num2str(Dfactor),'_Dt_',num2str(Dt),'_K_',num2str(K),'_a1_',num2str(a1),'_a2_',num2str(a2),'_transfers_',num2str(transfers),'_r0_',num2str(r0),'_dilutionF_',num2str(dilutionF),'_R_',num2str(R),'_dr_',num2str(dr),'_dt_',num2str(dt),'_n0IC-n0IM_',num2str(K0/2),'.mat'])

end