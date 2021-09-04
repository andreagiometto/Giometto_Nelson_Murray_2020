% Written by Andrea Giometto, MIT license
% This function runs the spatial model with growth-dependent diffusion
% coefficient for yeast cells. This code was used to generate the model
% simulations reported in Fig. 7 The initial conditions reflect the initial
% conditions in the experiments

function [] = spatial_model_experimental_inoculi(r0index)
    % initial density in the interior (M: K1, C: K2b)
    rhoIM0s=[0.214203,0.251825,0.222762,0.303108,0.299436,0.234341,0.349156,0.879377,0.625939,0.447904,0.579982,0.643444,0.631881,0.733991,0.649776,0.600602,0.652104,0.6219,0.602539,0.673123,0.68176,0.7224];
    rhoIC0s=[0.417637,0.321751,0.446792,0.532461,0.434654,0.496755,0.54775,0.572893,0.606651,0.566116,0.640378,0.481546,0.434229,0.509429,0.505914,0.481838,0.639716,0.54229,0.467361,0.64752,0.66321,0.59774];
    % initial density in the ring (M: K1, C: K2b)
    rhoRM0s=[0.23,0.59,0.44,0.45,0.53,0.45,0.53,0.605,0.61,0.58,0.576,0.705,0.755,0.792,0.673,0.656,0.625,0.672,0.719,0.781,0.812,0.698];
    rhoRC0s=[0.53,0.38,0.72,0.69,0.62,0.66,0.59,0.716,0.757,0.715,0.664,0.623,0.554,0.419,0.655,0.63,0.483,0.711,0.594,0.552,0.712,0.644];
    % width of the ring (cm)
    widths=[163.914,175.831,256.477,272.129,238.649,191.046,248.639,325.119,323.818,353.226,371.471,327.69,315.565,327.814,373.138,392.242,354.093,420.912,441.792,433.572,487.91,475.23]/10000.0;
    % radius of the inoculum (cm)
    radii0=[1531.25,1483.5,1485.25,1563.75,1525.25,1185.25,1425.25,2038.25,2090,2134.5,2521,2300.5,2363.25,2669,2845.25,2954.75,2592.5,2785.5,2542.75,3128.25,3052.75,3034.00]/10000.0;

    r0=radii0(r0index);
    tic
    
    R=1.5; % cm (radial size of the domain)
    h=0.36; % cm (height of the agar)
    Dg=0.024; %cm^2/h (glucose diffusion coefficient)
    Dt=0.003; % cm^2/h (toxin diffusion coefficient)
    gmax=0.3146; % 1/h (maximum growth rate)
    KS=2.0e-05; % g/mL (Monod growth constant)
    Y=1.2634e+10; % cells/g glucose (yield)
    cg0=0.02; % g glucose/cm^3 (initial concentration of glucose)
    a1=1.76e-09;  % mL/(cells h^2) (toxin production coefficient)
    a2=0.85e-09; % mL/(cells h^2) (toxin production coefficient)
    b=3.2e-09; % mL/(cells h) (toxin attachment coefficient)
    cellD=10^-3; % cells diameter
    Dfactor=1; % scaling factor for yeast diffusion coefficient

    K0=100000; % initial reference cell density cells/mL    
    ring=widths(r0index); % width of the ring
    n0IM=rhoIM0s(r0index)*K0; % density of invading strain in the interior of the droplet
    n0RM=rhoRM0s(r0index)*K0; % density of invading strain in the ring
    n0IC=rhoIC0s(r0index)*K0; % density of resident strain in the interior of the droplet
    n0RC=rhoRC0s(r0index)*K0; % density of resident strain in the ring

    transfers=1+13; % number of transfers
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
    cg=cg0*ones(M,N); % concentration profile of glucose

    r=linspace(0,R,N);
    z=linspace(0,h,M);

    % variable to store profile through time
    sample=1;
    c1sample=zeros(M,N,length(tsample));
    c2sample=zeros(M,N,length(tsample));
    cgsample=cg0*ones(M,N,length(tsample));
    n1sample=zeros(N,length(tsample));
    n2sample=zeros(N,length(tsample));

    % set initial condition
    for j=1:round((r0-ring)/dr)
        n1(j,1)=n0IM;
        n2(j,1)=n0IC;
    end
    for j=round((r0-ring)/dr+1):round(r0/dr)
        n1(j,1)=n0RM;
        n2(j,1)=n0RC;
    end
    for j=round(r0/dr)+1:N
        n2(j,1)=n0IC;
    end
    n1T=n1;
    n2T=n2;
    c1T=c1;
    c2T=c2;
    cgT=cg;

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
        cg(2:M-1,2:N-1)=cgT(2:M-1,2:N-1)+dt*Dg*(rinvmatrix(2:M-1,2:N-1).*(cgT(2:M-1,3:N)-cgT(2:M-1,1:N-2))/(2*dr)+(cgT(2:M-1,3:N)+cgT(2:M-1,1:N-2)-2*cgT(2:M-1,2:N-1))/(dr^2)+(cgT(1:M-2,2:N-1)+cgT(3:M,2:N-1)-2*cgT(2:M-1,2:N-1))/(dz^2));

        for j=2:(N-1) % top - bottom boundary
            % Note: cg must be updated before c1 and c2
            cg(1,j)=cg(2,j); % bottom
            cg(M,j)=-(dz*gmax*(n1(j,1)+n2(j,1)))/(2*Dg*Y)+(cg(M-1,j)-KS)/2+sqrt(cg(M-1,j)*KS+(dz*gmax/(2*Y*Dg)*(n1(j,1)+n2(j,1))+1/2*(KS-cg(M-1,j)))^2); % top
            c1(1,j)=c1(2,j); % bottom
            c1(M,j)=(c1(M-1,j)+a1*dz/Dt*n1(j,1)*cg(M,j)/(cg(M,j)+KS))/(1+b*dz/Dt*(n1(j,1)+n2(j,1)));
            c2(1,j)=c2(2,j); % bottom
            c2(M,j)=(c2(M-1,j)+a2*dz/Dt*n2(j,1)*cg(M,j)/(cg(M,j)+KS))/(1+b*dz/Dt*(n1(j,1)+n2(j,1)));
        end
        for i=1:M % left - right boundary: no flux 
            c1(i,1)=c1(i,2);
            c1(i,N)=c1(i,N-1);
            c2(i,1)=c2(i,2);
            c2(i,N)=c2(i,N-1);
            cg(i,1)=cg(i,2);
            cg(i,N)=cg(i,N-1);
        end
        for j=2:N-1 % n1, n2 update
            n1(j,1)=n1T(j,1)+dt*gmax*n1T(j,1)*cg(M,j)/(cg(M,j)+KS)+dt*Dfactor*cellD^2*gmax*cg(M,j)/(cg(M,j)+KS)/r(j)*(n1T(j+1,1)-n1T(j-1,1))/(2*dr)+2*dt*Dfactor*cellD^2*gmax*(cg(M,j+1)/(cg(M,j+1)+KS)-cg(M,j-1)/(cg(M,j-1)+KS))/(2*dr)*(n1T(j+1,1)-n1T(j-1,1))/(2*dr)+dt*Dfactor*cellD^2*gmax*cg(M,j)/(cg(M,j)+KS)*(n1T(j+1,1)-2*n1T(j,1)+n1T(j-1,1))/(dr^2)+dt*Dfactor*cellD^2*n1T(j,1)*(cg(M,j+1)/(cg(M,j+1)+KS)+cg(M,j-1)/(cg(M,j-1)+KS)-2*cg(M,j)/(cg(M,j)+KS))/(dr^2)+dt*Dfactor*cellD^2/r(j)*n1T(j,1)*(cg(M,j+1)/(cg(M,j+1)+KS)-cg(M,j-1)/(cg(M,j-1)+KS))/(2*dr)-dt*n1T(j,1)*c2T(M,j);
            n2(j,1)=n2T(j,1)+dt*gmax*n2T(j,1)*cg(M,j)/(cg(M,j)+KS)+dt*Dfactor*cellD^2*gmax*cg(M,j)/(cg(M,j)+KS)/r(j)*(n2T(j+1,1)-n2T(j-1,1))/(2*dr)+2*dt*Dfactor*cellD^2*gmax*(cg(M,j+1)/(cg(M,j+1)+KS)-cg(M,j-1)/(cg(M,j-1)+KS))/(2*dr)*(n2T(j+1,1)-n2T(j-1,1))/(2*dr)+dt*Dfactor*cellD^2*gmax*cg(M,j)/(cg(M,j)+KS)*(n2T(j+1,1)-2*n2T(j,1)+n2T(j-1,1))/(dr^2)+dt*Dfactor*cellD^2*n2T(j,1)*(cg(M,j+1)/(cg(M,j+1)+KS)+cg(M,j-1)/(cg(M,j-1)+KS)-2*cg(M,j)/(cg(M,j)+KS))/(dr^2)+dt*Dfactor*cellD^2/r(j)*n2T(j,1)*(cg(M,j+1)/(cg(M,j+1)+KS)-cg(M,j-1)/(cg(M,j-1)+KS))/(2*dr)-dt*n2T(j,1)*c1T(M,j);
        end
        n1(1,1)=n1(2,1);
        n1(N,1)=n1(N-1,1);
        n1T=n1;
        n2(1,1)=n2(2,1);
        n2(N,1)=n2(N-1,1);
        n2T=n2;
        c1T=c1;
        c2T=c2;
        cgT=cg;

        % store profiles at sampling time points
        if t*dt>=tsample(sample)
            c1sample(:,:,sample)=c1;
            c2sample(:,:,sample)=c2;
            cgsample(:,:,sample)=cg;
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
            cg=cg0*ones(M,N);
            cgT=cg;
            n1=n1/dilutionF;
            n2=n2/dilutionF;
            n1T=n1;
            n2T=n2;
            dilution=dilution+1;
        end
    end
    toc

    % save simulations data
    save(['model_output_growth_diffusion_cellD_',num2str(cellD),'_Dfactor_',num2str(Dfactor),'_Dg_',num2str(Dg),'_Dt_',num2str(Dt),'_KS_',num2str(KS),'_a1_',num2str(a1),'_a2_',num2str(a2),'_transfers_',num2str(transfers),'_r0_',num2str(r0),'_ring_',num2str(ring),'_dilutionF_',num2str(dilutionF),'_R_',num2str(R),'_dr_',num2str(dr),'_dz_',num2str(dz),'_dt_',num2str(dt),'_n0IC_',num2str(n0IC),'_n0IM_',num2str(n0IM),'_n0RC_',num2str(n0RC),'_n0RM_',num2str(n0RM),'.mat'])

end