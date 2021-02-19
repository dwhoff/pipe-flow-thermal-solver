clear,clc
close all

global fl_obj
%%%%%%%%%%%%%%%%%%%% - Pipe Flow Thermal Solver - %%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This code solves heat transfer problems in pipe flows. A description of
% the main settings is given below.
%
% fluid: defines the fluid object for look-up in Cantera
%        options: 'air', 'water', ..., 'user defined'
%        note: if 'user defined', fluid properties must be input manually.
%              This option must be selected if Cantera is not installed.
% shape: defines the shape of the pipe wall
%        options: 'round', 'rectangle'
% BC: defines the type of boundary condition used for the internal fluid
%        options: 'prescribed temp', 'conjugate'
%        note: if 'prescribed temp', the wall temperature must be input
%              manually
% wall_conduction: defines the discretization used in the solid domain
%        options: 'thin', 'thick'
%        note: BC must be set to 'conjugate'
% transient: toggles steady or unsteady analysis
%        options: true, false
%        note: BC must be set to 'conjugate'
% radiation: toggles radiation effects on the outer pipe wall
%        options: true, false
%        note: BC must be set to 'conjugate'
% write: toggles option to write data to text file
%        options: true, false
%
% Based on Stanford's course, ME352C: Convective Heat Transfer
% Written by Davis Hoffman, 2021

% Settings
fluid = 'water';
shape = 'round';
BC = 'conjugate';
wall_conduction = 'thin';
transient = false;
radiation = true;
write = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS - Parameters
D = 5e-3; % pipe inner diameter (if shape='round'), m
b1 = D; % pipe width (if shape='rectangle'), m
b2 = D; % pipe heigth (if shape='rectangle'), m
t = 2e-3; % pipe thickness, m
P = 101e3; % fluid pressure, Pa
T_in = 80; % inlet temperature, deg C
L = 2; % pipe length, m
u_bar = 0.5; % bulk velocity, m/s
rho_p = 5000; % pipe wall density, kg/m^3
kp = 3; % pipe wall thermal conductivity, W/m.K
cpp = 4000; % pipe wall specific heat capacity, J/kg.K
ho = 50; % external convection coefficient, W/m^2.K
Tinf = 0; % external temperature for convection boundary condition, deg C
eps = 0.85; % outer pipe wall emissivity
sigma = 5.67e-8; % Stefan-Boltzmann constant, W/m^2.K^4
Trad = 0; % external temperature for radiation, deg C
N = 4; % # of elements in radial direction (if wall-conduction='thick')
M = 36; % # of elements in axial direction
To_i = 0; % initial pipe wall temperature (if transient=true), deg C
dt = 1; % time step (if transient=true), s
time_f = 10; % final time (if transient=true), s
tol = 1e-3; % wall temperature tolerance, K
writename = 'output.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert temperatures to Kelvin scale
To_i = To_i+273.15;
T_in = T_in+273.15;
Tinf = Tinf+273.15;
Trad = Trad+273.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize fluid properties at inlet temperature
if strcmp(fluid,'water')
    fl_obj = Water;
elseif strcmp(fluid,'air')
    fl_obj = Air;
else
    fl_obj = GRI30('Multi'); % open Cantera ideal gas object
end
if strcmp(fluid,'user defined')
    disp('Enter fluid properties: ')
    rho = input(sprintf('\tDensity (kg/m^3): '));
    nu = input(sprintf('\tKinematic viscosity (m^2/s): '));
    cp = input(sprintf('\tSpecific heat (J/kg.K): '));
    k = input(sprintf('\tThermal conductivity (W/m.K): '));
else
    [rho,nu,cp,k] = get_fluid_props(fluid,P,T_in);
end
alpha = k/rho/cp; % thermal diffusivity, m^2/s
if strcmp(shape,'rectangle')
    D = 2*b1*b2/(b1+b2); % compute hydraulic diameter (if shape~='round')
end
Re = u_bar*D/nu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
dx = L/M; % pipe discretization length, m
if strcmp(wall_conduction,'thin')
    N = 1;
end
dr = t/N;
xc = dx/2:dx:L-dx/2; % cell center coordinates
rc = D/2+(dr/2:dr:t-dr/2);
xe = 0:dx:L; % cell face coordinates
re = D/2+(0:dr:t);
[X,R] = meshgrid(xc,rc);
if strcmp(shape,'round')
    Acond = pi*(re(2:end).^2-re(1:end-1).^2); % conduction area, m^2
    Aconvi = pi*D*dx; % convection area on inner pipe surface, m^2
    Aconvo = pi*(D+2*t)*dx; % convection area on outer pipe surface, m^2
    Across = pi/4*D^2; % internal cross section area, m^2
    m = rho_p*reshape(repmat(Acond,M,1),M*N,1)*dx; % element mass, m^2
elseif strcmp(shape,'rectangle')
    Acond = (b1+t)*(b2+t)-b1*b2;
    Aconvi = 2*(b1+b2)*dx; % convection area on inner pipe surface, m^2
    Aconvo = 2*((b1+2*t)+(b2+2*t))*dx; % convection area on outer pipe surface, m^2
    Across = b1*b2; % internalcross section area, m^2
    m = rho_p*reshape(repmat(Acond,M,1),M,1)*dx; % element mass, m^2
end
Bx = k*Acond/dx; % axial conductance parameter, W/K
Br = 2*pi*k*dx./log(rc(2:end)./rc(1:end-1)); % radial conductance parameter, W/K
if transient
    Nt = ceil(time_f/dt); % number of time steps
else
    Nt = 1;
end

% form BC vector (external heat transfer on elements, in Watts)
q_dot = zeros(M,1); % INPUT (when BC='conjugate')

% Initialize vector of eigenvalues and G coefficients
G = []; lambda_sq = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create block diagonal conduction matrix
A_cond = zeros(N*M);
for i = 1:N
    tmpx = Bx(i).*(diag(2*ones(M,1))-diag(ones(M-1,1),1)-diag(ones(M-1,1),-1));
    tmpx(1,1) = tmpx(1,1)/2;
    tmpx(end,end) = tmpx(end,end)/2;
    if i == 1 && strcmp(wall_conduction,'thick')
        tmpr = Br(i).*diag(ones(M,1));
    elseif i == N && strcmp(wall_conduction,'thick')
        tmpr = Br(i-1).*diag(ones(M,1));
    elseif i>1 && i<N && strcmp(wall_conduction,'thick')
        tmpr = (Br(i)+Br(i-1)).*diag(ones(M,1));
    else
        tmpr = zeros(M);
    end
    A_cond((i-1)*M+1:i*M,(i-1)*M+1:i*M) = tmpx+tmpr;
end
tmp = repmat(Br,M,1);
tmp = tmp(:);
tmpr = -diag(tmp,M)-diag(tmp,-M);
A_cond = A_cond+tmpr;

% create wall temperature vector, in K
if strcmp(BC,'prescribed temp')
    T = 10+10*sin(2*pi*xc'/L)+T_in; % INPUT (when BC='prescribed temp')
elseif strcmp(BC,'conjugate')
    T = To_i*ones(N*M,1);
end

clear tmpx tmpr Bx Br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate through time steps
time = 0;
figure(1); set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.125, 0.25, 0.75, 0.5]);
for j = 1:Nt
    fprintf('Timestep %d of %d, t = %.2f s\n',j,Nt,time)
    
    T_prev = T;
    T_old = zeros(N*M,1);
    resid = Inf;
    converged = 0;
    iter = 0;
    % iterate solver until film temperature converges
    while ~converged    
        dx_plus = alpha*dx/(2*u_bar*(D/2)^2);

        % form qb vector by summing infinite series until converged
        a = zeros(M,1);
        resid = Inf;
        a_tol = 1e-3;
        i = 1;
        while resid>a_tol
            if i>numel(lambda_sq) 
                [G(i),lambda_sq(i)] = get_eigenvalues(i,shape,Re);
            end
            a = a-G(i)/lambda_sq(i)*(exp(-lambda_sq(i)*(1:M)'*dx_plus)- ...
                exp(-lambda_sq(i)*(0:M-1)'*dx_plus));
            if i>1
                resid = abs(a(1)-a1_old)/a1_old;
            end
            a1_old = a(1);
            i = i+1;
        end
        qb = 4*pi*k*dx/dx_plus*a;

        %%
        % form linear system for the wall temperature
        % Steady: A*T=b
        % Transient: m*cp*dT/dt=-A*T+b

        % create lower triangular convection matrix
        tmp1 = zeros(M);
        tmp2 = zeros(M);
        for i = 1:M
            tmp1(i:end,i) = qb(1:end-i+1);
        end
        tmp2(2:end,:) = tmp1(1:end-1,:);
        A_conv1 = zeros(N*M);
        A_conv2 = zeros(N*M);
        A_conv1(1:M,1:M) = tmp1;
        A_conv2(1:M,1:M) = tmp2;

        % create external heat transfer matrix/vector
        tmp = diag(ho*Aconvo*ones(M,1));
        A_convo = zeros(N*M);
        A_convo(end-M+1:end,end-M+1:end) = tmp;
        b_convo = zeros(N*M,1);
        b_convo(end-M+1:end) = ho*Aconvo*Tinf;

        % create radiation heat transfer matrix/vector
        if radiation
            Rrad = 1./(eps*sigma*Aconvo*(T(end-M+1:end)+Trad).*(T(end-M+1:end).^2+Trad^2));
            tmp = diag(1./Rrad);
            A_rad = zeros(N*M);
            A_rad(end-M+1:end,end-M+1:end) = tmp;
            b_rad = zeros(N*M,1);
            b_rad(end-M+1:end) = Trad./Rrad;
        else
            A_rad = zeros(N*M);
            b_rad = zeros(N*M,1);
        end

        % create vector of external heat transfers
        b_ext = zeros(N*M,1);
        b_ext(1:M) = q_dot;

        % create vector that absorbs the inlet temperature terms
        b_inlet = zeros(N*M,1);
        b_inlet(1:M) = T_in*qb;

        if transient
            f = @(T) 1./(m*cpp).*(-(A_cond+A_conv1-A_conv2+A_convo+A_rad)*T+ ...
                (b_ext+b_convo+b_rad+b_inlet));
            T = advance_RK4(f,T_prev,dt);
        else
            if strcmp(BC,'conjugate')
                T = (A_cond+A_conv1-A_conv2+A_convo+A_rad)\ ...
                    (b_ext+b_convo+b_rad+b_inlet);
            end
        end
        
        resid = sqrt(mean((T-T_old).^2));
        if resid<tol
            converged = 1;
        else
            T_old = T;
        end
        
        % get convective heat flux into fluid
        qconv = (A_conv1(1:M,1:M)-A_conv2(1:M,1:M))*(T(1:M)-T_in);
        Tm = [T_in; T_in+cumsum(qconv)/(rho*cp*u_bar*Across)];

        % get fluid properties at the average mixed-mean temperature
        Tfilm = mean(Tm);
        if ~strcmp(fluid,'user defined')
            [rho,nu,cp,k] = get_fluid_props(fluid,P,Tfilm);
            alpha = k/rho/cp; % thermal diffusivity, m^2/s
        end
        
        iter = iter+1;

        if iter==1
            fprintf('\tIteration: %d\n',iter)
        else
            fprintf('\tIteration: %d, Residual: %.2e\n',iter,resid)
        end
    end
    time = time+dt;
    
    figure(1)
    if  strcmp(wall_conduction,'thin')
        plot(xc*1e3,T,'k-')
        xlabel('$x$ (mm)','interpreter','latex')
        ylabel('$T_o(x)$ ($^o$C)','interpreter','latex')
    elseif  strcmp(wall_conduction,'thick')
        pcolor(X*1e3,R*1e3,reshape(T,M,N)'); shading interp; colormap hot; colorbar
        xlabel('$x$ (mm)','interpreter','latex')
        ylabel('$r$ (mm)','interpreter','latex')
    end
    title('Wall temperature distribution (K)')
    plotfixer
    drawnow
end
To = T(end-M+1:end);

clear T_old T_prev tmp tmp1 tmp2 a1_old A_cond A_conv1 A_conv2 ...
    A_convo A_rad b_convo b_ext b_inlet b_rad q_dot qb Rrad m ...
    fl_obj G lambda_sq j

%%
% check errors
Re = u_bar*D/nu; % Reynolds number
if Re>2300
    disp('Warning: Re exceeds 2300')
end
h_max = max([3.66*k/D ho]);
Bi = h_max*t/kp; % Biot number
if Bi>0.1 && strcmp(wall_conduction,'thin') && strcmp(BC,'conjugate')
    disp('Warning: Bi exceeds 0.1. Consider a thick wall analysis for improved results.')
end
Br = nu*rho*u_bar^2/(k*abs(mean(T(1:M))-T_in)); % Brinkman number
if Br>11/4800
    disp('Warning: Br exceeds ?. Viscous heating may be important.')
end
tau_solid = min(1/4*rho_p*cpp*((D+2*t)^2-D^2)*[1/3.66/k 1/ho/(D+2*t)]); % solid time response, s
tau_adv = L/u_bar; % advection (flow-through) time scale
if tau_solid<tau_adv && transient
    disp('Warning: Quasi-steady assumption may be invalid.')
end

%%
% plot
if strcmp(wall_conduction,'thin')
    figure(1)
    plot(xc*1e3,To-273.15,'k-'); hold on
    plot((0:dx:L)*1e3,Tm-273.15,'k--')
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$T$ ($^o$C)','interpreter','latex')
    scale = axis;
    axis([scale(1) 1e3*L scale(3) scale(4)])
    legend({'$T_o(x)$','$T_m(x)$'},'interpreter','latex')
    legend boxoff
    plotfixer
elseif strcmp(wall_conduction,'thick')
    figure(1)
    TT = reshape(T,M,N)';
    pcolor(X*1e3,R*1e3,TT); shading interp; colormap hot; colorbar
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$r$ (mm)','interpreter','latex')
    title('Wall temperature distribution (K)')

    figure(2)
    plot(xc*1e3,To-273.15,'k-'); hold on
    plot((0:dx:L)*1e3,Tm-273.15,'k--')
    xlabel('$x$ (mm)','interpreter','latex')
    ylabel('$T$ ($^o$C)','interpreter','latex')
    scale = axis;
    axis([scale(1) 1e3*L scale(3) scale(4)])
    legend({'$T_o(x)$','$T_m(x)$'},'interpreter','latex')
    legend boxoff
end

figure(3)
plot(xc*1e3,qconv/Aconvi,'k-'); hold on
xlabel('$x$ (mm)','interpreter','latex')
ylabel('$\dot{q}''''_o$ (W/m$^2$)','interpreter','latex')
scale = axis;
axis([scale(1) 1e3*L scale(3) scale(4)])

figure(4)
plot(xc*1e3,qconv/Aconvi./(T(1:M)-(Tm(2:end)+Tm(1:end-1))/2),'k-'); hold on
xlabel('$x$ (mm)','interpreter','latex')
ylabel('$h$ (W/m$^2$K)','interpreter','latex')
scale = axis;
axis([scale(1) 1e3*L scale(3) scale(4)])

plotfixer

%%
% write to file
if write
    To_tmp = To;
    To = zeros(M+1,1);
    To(2:end-1) = (To_tmp(1:end-1)+To_tmp(2:end))/2;
    To(1) = 1.5*To_tmp(1)-0.5*To_tmp(2);
    To(end) = 1.5*To_tmp(end)-0.5*To_tmp(end-1);
    qconv_tmp = qconv;
    qconv = zeros(M+1,1);
    qconv(2:end-1) = (qconv_tmp(1:end-1)+qconv_tmp(2:end))/2;
    qconv(1) = 1.5*qconv_tmp(1)-0.5*qconv_tmp(2);
    qconv(end) = 1.5*qconv_tmp(end)-0.5*qconv_tmp(end-1);
    fileID = fopen(writename,'w');
    fprintf(fileID,'x (m)\t\tTo (K)\t\tTm (K)\t\tq (W/m^2)\n');
    for i = 1:M+1
        fprintf(fileID,'%.4e\t%.4e\t%.4e\t%.4e\n',xe(i),To(i),Tm(i),qconv(i)/Aconvi);
    end
    fclose(fileID);
end
clear To_tmp qconv_tmp fileID