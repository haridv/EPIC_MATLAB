%%%% Vijay Harid 08/17/21 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Attempt at developing at a Narrowband PIC simulation for triggered
%%%% emission modeling. This version attempts to simulate Chorus

%%%% _stat_v1: This iteration uses first order PIC (linear interpolation) %%%
%%%% This forces the NB aspect by keeping track of amplitude, phase, frequency, and wavenumber %%%
%%%% This version changes right BC to randomize phase for inflow
%%%% This version also extends the simulation space for particles and limits
%%%% the wave to an interior region...this stops weird particle BC issues from impacting
%%%% the waves

%%%% This version loops thourgh particle instead of grid points...expected
%%%% to be considerably faster.

%%%% This version caluclates B_input as a fraction of the trapping
%%%% threshold amplitude

%%%% This version uses linear interpolation to particle positions (instead
%%%% of spline)

%%%% This version updates particle interpolation correctly by adding 1 to index (incorrect in
%%%% old versions)

%%% This version updates does real annd imaginary part of wave instead of
%%% amplitiude and phase

%%%% This version removes the v_tr limitation for calculating resonant
%%%% currents.

%%%% This version switches to a new gyrophase angle that is only relative
%%%% to the incident wave. This allows us to not calculate the frequency
%%%% and wave number of the new wave but rather put it into the lorentz
%%%% force as a phase directly.

%%% This version feeds back phase to the particles

%%%%%%%%%%%%%%% v1_2_Chorus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This version increase gradient of background field to be aritifically
%%% sharp and correspondingly decreases size simulation space

%%% This version also increases bandwidth to 1kHz around carrier

%%% This version changes distribution function to be properly normalized in
%%% velocity space (all parameters are subsequently changed to be closer to
%%% Tao, Katoh, Hikishima

%%% This version makes density such that wp = 5*wc

%%% This version also saves figures when plotting

%%% This version tries to simulation Chorus by switching off injected
%%% waves and switches to inflow BC

%%% This version also increases growth rate dramatically to force a chorus
%%% Both the aniostropy and density are increased from previous sims

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc
close all;
%%%%% TEST NUMBER %%%%

addpath('C:\Users\haridv\Dropbox\Particle Pusher\NB_WPI')
T_num = 1; %%%test number
f_save_dir = ['EPIC_stat_TEST_v1_2_Chorus_',num2str(T_num)]; %%%% name of folder to save data
mkdir(f_save_dir)

%%%%%%% Physical constants %%%%%%%%%%%%%

ep=8.854e-12; %%% permittivity of free space
mu0 = 4*pi*1e-7; %%% permeability of free space
q=1.619e-19; %%% electronic charge
m_e=9.109e-31; %%% electron mass
c=3e8; %%% free space speed of light
i=sqrt(-1);

%%%%%%%

%%%%%%% Cold Plasma Density Properties
% Nc=5e8; %%% cold plasma density (1/m^3)
% Nc=20e8; %%% cold plasma density (1/m^3)

%%%%%% Background B-field Paramter and z-grid
B00=5.3027e-07; %%% Equatorial B-field
% sigma=5.1886e-15; %%% Parabolic B-field curvature'
sigma= (5.1886e-15)*1024; %%% Parabolic B-field curvature (sharper version)
zeq=0; %%% equatorp position
% zmin = -2000e3; %%%% minimum z-location
% zmax = 2000e3; %%%% maximum z-location
zmin = -300e3; %%%% minimum z-location
zmax = 300e3; %%%% maximum z-location
% zmin = -100e3; %%%% minimum z-location
% zmax = 100e3; %%%% maximum z-location
% z_min_wave  = 0.5*zmin; %%% min particle location
% z_max_wave  = 0.5*zmax;%%% max particle location
z_min_wave  = 0.8*zmin; %%% min particle location
z_max_wave  = 0.8*zmax;%%% max particle location
% z_min_wave  = 1*zmin; %%% min particle location
% z_max_wave  = 1*zmax;%%% max particle location



%%% Initial distribution paramters %%%
beta = 0.3;
%%%


% sigma= 2*5.1886e-13; %%% Parabolic B-field curvature —> artificailly high
% L=4.9; %%% L_shell
% Re = 6370000; %%% Earth radius
% hm=3.5*Re;
% eps_m = (1/L)*(Re+hm)/Re;
% alpha_eq = asin(sqrt( eps_m^3/sqrt(1+3*(1-eps_m)) )); %%% Equaotrial pitch angle
% alpha_lc =  asin(sqrt(abs(B0/B00)).*sin(alpha_eq));
% alpha_p = 40*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Input Magnetic field Signal Paramters %%%

f_0=3e3; %%% input frequency
% f_0=4.9e3; %%% input frequency
w=2*pi*f_0;
% B_input = (1e-2)*10^-12;  %%% Input wave magnetic field
% B_input = (1)*10^-12;  %%% Input wave magnetic field
% B_low = B_input/1e5; %%% lowest amplitude
% dur = 0.5;
% dur = 6000/min(wc);
w_0=2*pi*f_0;
f_ramp = 0; %%% rano frequency for rising or falling tones
% w_min=w_0-2*pi*20;
% w_max=w_0+2*pi*20;

w_min=w_0-2*pi*1000;
w_max=w_0+2*pi*1000;
% dt=2*zmax/(0.1*c);
% dt=8.3058e-06; %%% time step

%%%%%%%%%%%%%


%%%% phase space grid dimensions
Num_vll = 200;
% Num_vll = 400;
Num_vL = 5;
Num_phi = 50;
% Num_phi = 12;
Num_z = 500; %%% number of spatial grid points for wave/currents
% Num_z = 1000; %%% number of spatial grid points for wave/currents
Num_z_part = Num_z; %%% number of spatial grid points for particles
% Num_alpha = 4;
%%%%
Num_p = Num_vll*Num_vL*Num_phi*Num_z_part; %%% Total number of particles
[num2str(Num_p/1e6), ' million particles']
% Num_alpha = 4;
%%%%

%%% Create spatial Grid and gyrophase grid
z_val=linspace(zmin,zmax,Num_z); %%% z-grid
dz_val = z_val(2)-z_val(1);
%%%% Get indices for limiting wave propoagtion space
[~,indz_min] = min(abs(z_val-z_min_wave));
[~,indz_max] = min(abs(z_val-z_max_wave));

%%%% Create phase vector (particles)
phi  = linspace(0*pi,2*pi,Num_phi+1);
phi = phi(1:end-1); %%% remove double count at 2*pi
%%%


%%% Create physics paramter vectors that depened on spatial grid %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0=B00*(1+sigma*(z_val-zeq).^2); %%% Background magnetic field
d_B0=2*B00*sigma*(z_val-zeq); %%% field gradient
wc=abs(q*B0/m_e); %%% gyrofrequency
d_wc=q*d_B0/m_e;
del = (max(wc) - min(wc))/(min(wc)*abs(z_val(end)*z_val(end/2))); %%%% map to magnetic bottle configuration
% wp=sqrt(Nc*q*q/m_e/ep); %%% cold plasma frequency
wp= 5*min(wc); %%% cold plasma frequency
Nc = wp*wp*ep*m_e/q/q; %%% cold plasma density
%%% Initialize wave/spatial paramters
w = w_0 + 0*z_val;
% k = sqrt((w/c)^2 - w*wp*wp./(w-wc)/c/c); %%% cold plasma whistler wavenumber
v_g = 2*c*sqrt(w).*(wc-w).^(3/2)./wc./wp;
v_p = c*sqrt(w).*sqrt(wc-w)./wp;
k = w./v_p;
v_res = (wc - w)./k; %%% resonance velocity
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % %%%%% Create parallel velocity grid based on limits of gyroresonance
% and adiabatic motion
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Resonance velocity due to lowest possible frequency
vg_check1 = 2*c*sqrt(w_min).*(wc-w_min).^(3/2)./wc/wp; %%% lowest frequency group velocity
vp_check1 = c*sqrt(w_min)*sqrt(wc-w_min)/wp; %%% lowest frequency phase velocity
k_check1 = w_min./vp_check1; %%% wave number correspondonding lowest grequnecy
vresz_check1 = ((wc - w_min)./k_check1); %%%% non-relativistic resonance due to lowest frequency
vll_max = max(vresz_check1); %%% top edge of grid in vll

%%%%% Resonance velocity due to higest possible frequency %%%%%%
vg_check2 = 2*c*sqrt(w_max)*(wc-w_max).^(3/2)./wc/wp;
vp_check2 = c*sqrt(w_max)*sqrt(wc-w_max)/wp;
k_check2 = w_max./vp_check2;
vresz_check2 = ((wc - w_max)./k_check2); %%%% non-relativistic resonance due to highest frequency
% veq = vresz_check2;
% wc_end = wc(end);
% wc_eq = wc(round(length(wc)/2));
% vll_min = veq*sqrt(1-(wc_end/wc_eq)*sin(alpha_eq*pi/180)^2); %%% bottom edge of grid in vll based on adiabatic motion
vll_min = min(vresz_check2);
% vll_min = min(v_res)-(vll_max-min(v_res));

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_ll = linspace(vll_min,vll_max,Num_vll);

%%%%% Create vL grid and physical paramters that depend on vL %%%%

if Num_vL == 1 %%%% if only one particle
    v_L = min(v_res)*tan(70*pi/180);
else
    v_L = linspace(min(v_res)*tan(30*pi/180),min(v_res)*tan(70*pi/180),Num_vL);
end


%%% Calculate threshold amplitude
B_th = abs((1./(abs(q)*min(k)*max(v_L/m_e))).*(0.5*k.*max(v_L.^2)./2./wc +1.5*v_res).*d_wc); %%% S-paramters for df/dt = 0
%%%%

%%% Determine time step
% dt = 2*pi/5/max(k*(vll_max-vll_min)); %%% time_step: fifth of the maximum gyrofrequency
% dt = 2*pi/5/max(wc); %%% time_step: fifth of the maximum gyrofrequency
dt = min([2*pi/5/max(k*(vll_max-vll_min)), dz_val/min(v_g),dz_val/min(v_res)]); %%% min time step of all relavent physics
%%%
dur = 6000/min(wc);
M = round(dur/(dt));
t = linspace(0,dur,M);
f_in = f_0 + f_ramp*t;
% % B_input = 0.3*B_th(indz_max); %%% input is some fraction of threshold amplitude at input region
% B_input = 0.9*B_th(indz_max); %%% input is some fraction of threshold amplitude at input region
% % B_low = B_input/1e5; %%% lowest amplitude
B_low = 0; %%% lowest amplitude
% B_in = B_input+0*t; %%% input wave amplitude
% % %%% Make Short Pulse
% B_in = B_in.*tukeywin(length(t),.1)';
% % B_in = B_in.*tukeywin(length(t),0.3)';
% % T_extra = .2; %%% extra simulation time
% % T_extra = (z_max_wave-z_min_wave)./max(v_g); %%% extra simulation time to propagate wave through space (lower bound)
% % T_extra = (z_max_wave-z_min_wave)./min(v_g); %%% extra simulation time to propagate wave through space (lower bound)
% T_extra = 2*(z_max_wave-z_min_wave)./min(v_g); %%% extra simulation time to propagate wave through space (upper bound)
%
% % % % B_low = B_input/1000; %%% lowest amplitude
% B_in(B_in<B_low)=B_low; %%% set lowest amplitude
% B_in = [B_in, B_low*ones(1,round(T_extra/dt))]; %%%% add T_extra zeros to allow wave to pass through
% f_in = [f_in, f_in(end)*ones(1,round(T_extra/dt))];
% t = 0:dt:length(B_in)*dt - dt;  %%% fix time vector
% M = length(t); %%% fix length of time vector
% phi_in = (pi/4)*ones(1,M); %%% input phase
% phi_in = (0)*ones(1,M); %%% input phase
% B_inr = B_in.*cos(phi_in); %%% input real component
% B_ini = B_in.*sin(phi_in); %%% input imaginary component
% B_inr = 0.*cos(phi_in); %%% input real component
% B_ini = 0.*sin(phi_in); %%% input imaginary component


%%% Initialize waves and currents for first time step
% B_w = 0*z_val + B_low; %%% small initial wave
% B_w = 0*z_val + B_input; %%% large initial wave (FOR TESTING)
B_w = 0*z_val; %%% small initial wave
E_w = B_w.*w./k; %%% Plane wave relationship (good approx for coherent waves)

Jr = zeros(1,length(z_val));
Ji = zeros(1,length(z_val));
% B_wr = zeros(1,length(z_val))+B_inr(1);
% B_wi = zeros(1,length(z_val))+B_ini(1);
B_wr = zeros(1,length(z_val));
B_wi = zeros(1,length(z_val));
B_w = sqrt(B_wr.^2 + B_wi.^2); %%% conpute wave amplitude
phi_w = phase(B_wr+1i*B_wi); %%% compute wave phase
E_w = B_w.*w./k; %%% Plane wave relationship (good approx for coherent waves)
% w_tr_max = sqrt(abs(q)*max(k)*max(B_in)*max(v_L/m_e)); %%% non-relativsitc
% v_tr = 2*w_tr_max/min(k); %%%% non-relatvistic trapping width


%%% IMPORTANT: This if the Fraction of trapping width away from resonance
%%% to search for particles (ensures resonant particles are dominating).
% v_tr_fr = 1.5; %%%
% v_tr_fr = 100; %%%



%%
%%%% Create Particle initial conditions based on phase-space grid points %%%
% z_ind = 1:10:length(z_val);
% % [~,z_ind] = min(abs(z_val-100e3)); %%% single particle
% z_part = z_val(z_ind); %%%% particle positions
% z_part = linspace(z_min_part,z_max_part,Num_z_part); %%% particle positions
z_part = linspace(zmin,zmax,Num_z_part); %%% particle positions
dz_p = z_part(2)-z_part(1); %%% grid spacing
dphi_p = phi(2)-phi(1); %%% grid spacing
dvll_p = v_ll(2)-v_ll(1); %%% grid spacing
dvL_p = v_L(2)-v_L(1); %%% grid spacing


[V_ll_all,V_L_all,PHI_all,Z_all] = ndgrid(v_ll,v_L,phi,z_part); %%%% generate 4D grid
% [V_ll_w,V_L_w,PHI_w,Z_w] = ndgrid(v_ll,v_L,phi,z_val); %%%% generate 4D grid
% v_ll = linspace(0.5*c,0.75*c,Num_vll);
% Z_G=[z_val(end/2) z_val(3*end/4)];
% V_L_all = V_ll_all.*tan(ALPHA_all);

%%%% particle coordinates in phase-space
s(1,:) = V_ll_all(:)';
s(2,:) = V_L_all(:)';
s(3,:) = PHI_all(:)';
s(4,:) = Z_all(:)';

%%%% randomized particle coordinates in phase-space
% s(1,:) = vll_min+(vll_max-vll_min)*rand(1,Num_p);
% s(2,:) = min(v_L)+(max(v_L)-min(v_L))*rand(1,Num_p);
% s(3,:) = 2*pi*rand(1,Num_p);
% s(4,:) = min(z_val)+(max(z_val)-min(z_val))*rand(1,Num_p);

%%%% relabel according to row number...is this redundant?
% vll = f(1,:);
% vL = f(2,:);
% phii = f(3,:);
% z = f(4,:);

%%%% Store initial particle coordinates
% vll_init = s(1,:);
% vL_init = s(2,:);
% phii_init = s(3,:);
% z_init = s(4,:);

%%
%

%%%%% Track particles in time (make this a loop for each tine step)
%%%% Note*: Throws away particles that mirror (vll<0) or leave the space
%%%% abs(z)>max(z_val)
% dt = -abs(dt);  %%% make deltat negative (scattering backwards)
% t_start = round(2*M/3); %%% time index to start running backwards
% t_step = t_start;
% Z_w_f = Z(min(find(B_W(t_start,:) > 4e-11)));
% Z_w_b = Z_w_f + dur*v_g(round(mean(find(B_W(t_start,:) > 4e-11))));

%%%% Re-initialize particle coordinates
% s(1,:) = V_ll_all(:)';
% s(2,:) = V_L_all(:)';
% s(3,:) = PHI_all(:)';
% s(4,:) = Z_all(:)';

%%%% initial distribution function as seen by particles
v_ll_th = .2*c;
% v_L_th  = .35*c;
v_L_th  = .45*c;
Nh = (1e-2)*min(Nc);

F_part = bi_max_vel(s(1,:),s(2,:),s(4,:),B00,del,beta,v_ll_th,v_L_th, Nh); %%% get PSD
N_part = F_part.*s(2,:)*dphi_p*dvll_p*dvL_p; %%%% density of each super particle
% w_part = dz_val/dz_p; %%% particle weight (making each super-particle wider but shorter to have the same total charge)
w_part = 1; %%% no weight (accounted for with different densities)

s_temp = s;

%%%% Initialize Wave matrices %%%
%%%% Create Wave fields/parameters over space-time grid
% B_W = zeros(M,length(z_val)); %%% B-field amplitude matrix
B_Wr = zeros(M,length(z_val)); %%% B-field real amplitude matrix
B_Wi = zeros(M,length(z_val)); %%% B-field imaginary amplitude matrix
% B_W = zeros(M,length(z_val)); %%% B-field amplitude
% PHI_W = zeros(M,length(z_val)); %%% B-field amplitude
% B_W(1,:) = B_w; %%% first time step
% B_Wr(1,:) = B_w; %%% first time step
% B_Wi(1,:) = 0; %%% first time step
% B_wr = B_w; %%% first time step
% B_wi = 0*B_w; %%% first time step

% E_W = zeros(M,length(z_val)); %%% E-field amplitude matrix
% W_W = zeros(M,length(z_val)); %%% Frequency matrix
% K_W = zeros(M,length(z_val)); %%%
JR = zeros(M,length(z_val)); %%%
JI = zeros(M,length(z_val)); %%%

Jr = zeros(1,length(z_val)); %%% first time step
Ji = zeros(1,length(z_val)); %%% first time step

% t_st_plot = 10; %%% time step to plot
t_st_plot = floor((M/1000)); %%% time step to plot




% ind_keep = 1:size(f,2); %%% initial array of indices to keep tracking
counter = 1;
figure(2);
set(gcf, 'Position', get(0, 'Screensize')); %%% Make fullscreen

for t_step = 1:M %%% time loop for tracking particles
    tic
    %%%%% Load wave data from space-time matrices
    %         f_temp = f(:,ind_keep); %%% temporary
    
    
    %     B_w = B_W(t_step,:);
    %     E_w = E_W(t_step,:);
    %     w = W_W(t_step,:);
    %     k = K_W(t_step,:);
    %     w= w(:);
    %     k = k(:);
    %%%% scatter particles forwards with RK4
    %     k1  = dt*push_RK4(s,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
    %     k2  = dt*push_RK4(s+k1/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
    %     k3  = dt*push_RK4(s+k2/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
    %     k4  = dt*push_RK4(s+k3,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
    %     s = s + k1/6 + (k2+k3)/3 + k4/6;
    
    %         k1  = dt*push_RK4(f(:,ind_keep),wc,d_wc,B_w,E_w,w,k,z_val,q,m,c);
    %         k2  = dt*push_RK4(f(:,ind_keep)+k1/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m,c);
    %         k3  = dt*push_RK4(f(:,ind_keep)+k2/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m,c);
    %         k4  = dt*push_RK4(f(:,ind_keep)+k3,wc,d_wc,B_w,E_w,w,k,z_val,q,m,c);
    %         f(:,ind_keep) = f(:,ind_keep) + k1/6 + (k2+k3)/3 + k4/6;
    % %         t_step = t_step-1;
    %         if (t_step == 0)
    %             break
    %         end
    
    %%% Push particles with RK4
    k1  = dt*L_force_stat(s,wc,d_wc,B_w,E_w,phi_w,w,k,z_val,q,m_e,c);
    k2  = dt*L_force_stat(s+k1/2,wc,d_wc,B_w,E_w,phi_w,w,k,z_val,q,m_e,c);
    k3  = dt*L_force_stat(s+k2/2,wc,d_wc,B_w,E_w,phi_w,w,k,z_val,q,m_e,c);
    k4  = dt*L_force_stat(s+k3,wc,d_wc,B_w,E_w,phi_w,w,k,z_val,q,m_e,c);
    s = s + k1/6 + (k2+k3)/3 + k4/6;
    %%%
    if sum(isnan(s))>0
        keyboard;
    end
    
    
    %%%% find particles that are above or below simulation vll-space or
    %%%% exiting spatial simulation space (pseudo-periodic)
    %%% Need to fix this to make faster and more accurate
    %     tic
    ind_leave = find((s(1,:)<=min(v_ll))|(s(1,:)>=max(v_ll))|(s(4,:)>=max((z_val))));
    s(4,ind_leave) = -s_temp(4,ind_leave); %%% reintroduce particles at their negative position (adiabatic motion)
    s(1,ind_leave) = s_temp(1,ind_leave);
    s(2,ind_leave) = s_temp(2,ind_leave);
    
    s(3,ind_leave) = 2*pi*rand(1,length(ind_leave)); %%% randomize incoming phase when re-introducing particles on other side of space
    s_temp = s; %%% re-initialize particles in temporary array
    % %     %     toc
    
    %%%% loop through particles and calculate Je and Jb
    Jr = zeros(1,length(z_val));
    Ji = zeros(1,length(z_val));
    
    %     tic
    for ss = 1:length(s)
        
        %%% iterpolate spatial quantities to particle position (linear)
        %         B_w_part = interp1(z_val,B_w,s(4,ss),'linear');
        %         E_w_part = interp1(z_val,E_w,s(4,ss),'linear');
        %         w_part = interp1(z_val,w,s(4,ss),'linear');
        %         k_part = interp1(z_val,k,s(4,ss),'linear');
        %         wc_part = interp1(z_val,wc,s(4,ss),'linear');
        %         d_wc_part = interp1(z_val,wc,s(4,ss),'linear');
        %%%
        
        %       k1  = dt*push_RK4(s,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
        %     k2  = dt*push_RK4(s+k1/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
        %     k3  = dt*push_RK4(s+k2/2,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
        %     k4  = dt*push_RK4(s+k3,wc,d_wc,B_w,E_w,w,k,z_val,q,m_e,c);
        %     s = s + k1/6 + (k2+k3)/3 + k4/6;
        
        %%% Step particles forward with RK4 (NOTE*: UPDATE THIS to make fields inteprolated at each intermediate step!!!!)
        %         k1 = dt*L_force(s,wc_part,d_wc_part,Bw_part,Ew_part,w_part,k_part,q,m);
        %         k2 = dt*L_force(s+k1/2,wc_part,d_wc_part,Bw_part,Ew_part,w_part,k_part,q,m);
        %         k3 = dt*L_force(s+k2/2,wc_part,d_wc_part,Bw_part,Ew_part,w_part,k_part,q,m);
        %         k4 = dt*L_force(s+k3,wc_part,d_wc_part,Bw_part,Ew_part,w_part,k_part,q,m);
        %         s = s + k1/6 + (k2+k3)/3 + k4/6;
        %         %%%
        
        
        z_pos = s(4,ss); %%% get particle position
        vll_pos = s(1,ss); %%% get particle parallel velocity
        
        %%% get up and lower indices
        %         z_up = ceil((z_pos - z_val(1))/dz_val);
        z_up = ceil((z_pos - z_val(1))/dz_val)+1;
        
        
        
        %         z_down = floor((z_pos - z_val(1))/dz_val);
        %%%% get weights to neighboring points
        %         w_z_up = abs(z_val(z_up)-z_pos)/dz_val; %%%
        %         w_z_down = abs(dz_val-del_z_up*dz_val)/dz_val; %%%
        
        
        if (z_up >= indz_min)&&(z_up<=indz_max) %%%if particle inside wave-allowed space
            w_tr_max = sqrt(abs(q)*(k(z_up))*B_w(z_up)*max(v_L/m_e)); %%% get trapping frequency (rad/sec)
            v_tr = 2*w_tr_max/(k(z_up)); %%%% non-relatvistic trapping velocity
            vres_loc = v_res(z_up); %%% get local resonance velocity
            
            %%% note this assumes same resonance veclocity for both left (Down) and
            %%% right(up) grid points...its uses the one to the right of
            %%% particle. CHANGE THIS FOR NEXT ITERATION.
            
            w_z_down = abs(z_val(z_up)-z_pos)/dz_val; %%%
            w_z_up = abs(dz_val-w_z_down*dz_val)/dz_val; %%%
            %         %%% Find local v_tr here
            
            %         B_th = abs((1./(abs(q)*(k(zz))*max(v_L/m_e))).*(0.5*k(zz).*max(v_L.^2)./2./wc +1.5*vres_loc).*d_wc(zz)); %%% S-paramters for df/dt = 0
            
            
            %%% Only calculate currents if particle is within some number of trapping widths from resonance
            %             if (abs(s(1,ss)-vres_loc)<v_tr_fr*v_tr) %%% only calculate currents if particle is within some number of trapping widths from resonance
            %%% Accumulate currents (Je)
            Jr(z_up) = Jr(z_up) -q*w_z_up*w_part.*sum(N_part(ss).*s(2,ss).*sin(s(3,ss))); %%% update Je for upper index
            Jr(z_up-1) = Jr(z_up-1) -q*w_z_down*w_part.*sum(N_part(ss).*s(2,ss).*sin(s(3,ss))); %%% update Je for lower index
            
            %%% Accumulate currents (Jb)
            Ji(z_up) = Ji(z_up) -q*w_z_up*w_part.*sum(N_part(ss).*s(2,ss).*cos(s(3,ss))); %%% update Je for upper index
            Ji(z_up-1) = Ji(z_up-1) -q*w_z_down*w_part.*sum(N_part(ss).*s(2,ss).*cos(s(3,ss))); %%% update Je for lower index
            %             end
            
            
            
        end
    end
    %     toc
    
    
    [t_step, M]
    
    
    D = v_g*dt; %%% distance moved by wavefronts in one time step
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wave equation only inside simulation space.
    
    %%% Propagate wave forward with currents (FSL method)
    %     Je = zeros(1,length(z_val));
    %     Jb = zeros(1,length(z_val));
    
    
    %%% LINEAR INTERPOLATION OPTION
    
    %         B_w(indz_min:indz_max-1)  = interp1(z_val(indz_min:indz_max),B_w(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*Je(indz_min:indz_max-1)*dt;
    %
    %         B_w(1:indz_min-1) = B_low; %%% explicitly remove field outside sim space (left)
    %         B_w(indz_max+1:end) = B_low; %%% explicitly remove field outside sim space (right)
    %         B_w(B_w<0) = B_low;
    %
    %         %%% Calculate phase
    %         phi_w(indz_min:indz_max-1)  = interp1(z_val(indz_min:indz_max),phi_w(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*(Jb(indz_min:indz_max-1)./B_w(indz_min:indz_max-1))*dt;
    %
    %%% SPLINE INTERPOLATION OPTION
    
    %%%% Amplitude and Phae Method
    %     B_w(indz_min:indz_max-1)  = spline(z_val(indz_min:indz_max),B_w(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*Je(indz_min:indz_max-1)*dt;
    
    %     B_w(1:indz_min-1) = B_low; %%% explicitly remove field outside sim space (left)
    %     B_w(indz_max+1:end) = B_low; %%% explicitly remove field outside sim space (right)
    %     B_w(B_w<0) = B_low;
    
    %%%% Calculate phase
    %     phi_w(indz_min:indz_max-1)  = spline(z_val(indz_min:indz_max),phi_w(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*(Jb(indz_min:indz_max-1)./B_w(indz_min:indz_max-1))*dt;
    
    %%%% Real and imaginary method
    B_wr(indz_min:indz_max-1)  = spline(z_val(indz_min:indz_max),B_wr(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*Jr(indz_min:indz_max-1)*dt;
    B_wi(indz_min:indz_max-1)  = spline(z_val(indz_min:indz_max),B_wi(indz_min:indz_max),z_val(indz_min:indz_max-1)+D(indz_min:indz_max-1))-0.5*mu0*v_g(indz_min:indz_max-1).*Ji(indz_min:indz_max-1)*dt;
    
    
    
    B_wr(1:indz_min-1) = 0; %%% explicitly remove field outside sim space (left)
    B_wr(indz_max+1:end) = 0; %%% explicitly remove field outside sim space (right)
    B_wi(1:indz_min-1) = 0; %%% explicitly remove field outside sim space (left)
    B_wi(indz_max+1:end) = 0; %%% explicitly remove field outside sim space (right)
    %%%%
    
    
    %%% first order inflow boundary condition (take on the value downstream
    %%% of boundary)
    B_wr(indz_max) = B_wr(indz_max-1);
    B_wi(indz_max) = B_wi(indz_max-1);
    
    %%% Inject wave from boundary
    %
    %     B_wr(indz_max) = B_inr(t_step); %%%  real part injected
    %     B_wi(indz_max) = B_ini(t_step); %%% imaginary part injected
    
    %     B_w(indz_max) = B_in(t_step); %%% wave amplitude injected
    %     phi_w(indz_max) = phi_in(t_step);
    
    B_w = sqrt(B_wr.^2 + B_wi.^2); %%% conpute wave amplitude
    phi_w = phase(B_wr+1i*B_wi); %%% compute wave phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% When to plot
    if ((rem(t_step,t_st_plot)==0)||(t_step == 1))
        %      plot(z_val/1000,B_w)
        %
        %         subplot(221)
        %         plot(z_val/1000,v_res/c,'r',z_val/1000,vresz_check1/c,'r--',z_val/1000,vresz_check2/c,'r--') %%% plot wave amplitude
        %         %         plot(z_val/1000,v_res/c, z_val/1000,(v_res+v_tr)/c,'r--',z_val/1000,(v_res-v_tr)/c,'r--') %%% plot wave amplitude
        %         %         plot(z_val/1000,v_res/c,'r--','linewidth',2) %%% plot wave amplitude
        %         hold on;
        %         scatter(s(4,:)/1000,s(1,:)/c,'k') %%% plot parallel velopcity
        %         %         scatter(s(4,indz)/1000,s(1,:)/c,'k') %%% plot parallel velocity close to resonance
        %         xlabel('z [km]')
        %         ylabel('v_{ll}/c')
        %         title(['\delta z_{min} = ',num2str(min(v_ll)*t_step*dt/1e3),' km'])
        %         set(gca,'fontsize',16)
        %         drawnow
        %         hold off;
        
        figure(2);
        
        subplot(221)
        %         plot(z_val/1000,phi_w,'r') %%% plot phase of wave
        plot(z_val*min(wc)/c,phi_w,'r')
        %         xlabel('z [km]')
        xlabel('z\omega_c/c')
        ylabel('\phi_{w} [rad]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
        title(['t*\omega_c = ' num2str(t_step*dt*min(wc))])
        set(gca,'fontsize',16)
        %         ylim([0 1.2*B_input]/1e-12)
        drawnow
        %
        %         plot(z_val/1000,k,'r') %%% plot wave amplitude
        %         % plot(z_val/1000,dk+k,'r') %%% plot wave amplitude
        %         xlabel('z [km]')
        %         ylabel('Wavenumber [rad/m]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
        %         set(gca,'fontsize',16)
        %         %         ylim([0 1.2*B_input]/1e-12)
        %         drawnow
        
        
        
        subplot(222) %%% Real and imaginary parts
        %         plot(z_val/1000,phi_w,'r') %%% plot phase of wave
        %         plot(z_val/1000,B_th/1e-12,z_val/1000,B_wr/1e-12,'r',z_val/1000,B_wi/1e-12,'k',z_val/1000,B_w/1e-12,'b--') %%% plot wave amplitude
        plot(z_val*min(wc)/c,B_th/1e-12,z_val*min(wc)/c,B_wr/1e-12,'r',z_val*min(wc)/c,B_wi/1e-12,'k',z_val*min(wc)/c,B_w/1e-12,'b--') %%% plot wave amplitude
        % xlabel('z [km]')
        xlabel('z\omega_c/c')
        ylabel('B_{w} [pT]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
        title(['t*\omega_c = ' num2str(t_step*dt*min(wc))])
        set(gca,'fontsize',16)
        legend('B_th','B_r','B_i','|B_w|')
        %         ylim([0 1.2*B_input]/1e-12)
        drawnow
        
        
        %         plot(z_val/1000,phi_w,'r') %%% plot phase of wave
        %         xlabel('z [km]')
        %         ylabel('\phi_{w} [rad]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
        %         set(gca,'fontsize',16)
        %         %         ylim([0 1.2*B_input]/1e-12)
        %         drawnow
        %
        
        subplot(223)
        plot(z_val*min(wc)/c,phase(Jr+1i*Ji),'k') %%% plot currents
        %         xlabel('z [km]')
        xlabel('z\omega_c/c')
        ylabel('\Phi_J [rad]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
       title(['t*\omega_c = ' num2str(t_step*dt*min(wc))])
        set(gca,'fontsize',16)
        drawnow
        
        subplot(224)
        %         plot(z_val/1000,Jr,'k',z_val/1000,Ji,'r',z_val/1000,abs(Jr+1i*Ji)) %%% plot currents
        plot(z_val*min(wc)/c,Jr,'k',z_val*min(wc)/c,Ji,'r',z_val*min(wc)/c,abs(Jr+1i*Ji)) %%% plot currents
        xlabel('z\omega_c/c')
        %         xlabel('z [km]')
        ylabel('J [A/m^2]')
        %         title(['t = ' num2str(t_step*dt), ' sec'])
        title(['t*\omega_c = ' num2str(t_step*dt*min(wc))])
        legend('J_r','J_i','|J|')
        set(gca,'fontsize',16)
        drawnow
        
        
        
        %         length(ind_leave); %%% for TESTING
        %         if (isempty(ind_keep)) %%% if all particles have left.
        %             break;
        %         end
        
        fsave_name =  [f_save_dir,'/',num2str(counter),'.tif']; %%% name to save figure as
        saveas(gcf,fsave_name) %%% save figure as tif
        
        counter = counter + 1; %%% increment counter
    end
    
    
    
    
    
    
    %%% Update wavenumber using derivative of phase
    %     k =
    
    %%% Use cold plasma relationship (HACK! may not be correct!)
    %%% This is possibly reduntant
    v_g = 2*c*sqrt(w).*(wc-w).^(3/2)./wc./wp;
    v_p = c*sqrt(w).*sqrt(wc-w)./wp;
    k = w./v_p;
    %%% Update E-field use plane wave relationship
    E_w = B_w.*w./k;
    
    %%% Save fields/currents in matrix
    %     B_W(t_step,:) = B_w;
    %     PHI_W(t_step,:) = phi_w;
    B_Wr(t_step,:) = B_wr;
    B_Wi(t_step,:) = B_wi;
    JR(t_step,:) = Jr;
    JI(t_step,:) = Ji;
    %     B_W(t_step,:) = B_w;
    %     E_W(t_step,:) = E_w;
    %     W_W(t_step,:) = w;
    %     K_W(t_step,:) = k;
    %      B_w_real(1:end-1)  = spline(z_val_0(1:end),B_w_real(1:end),z_val_0(1:end-1)+D(1:end-1))-0.5*mu*v_g(1:end-1).*(JE_g(t_step,1:end-1)).*cos(Phi(t_step-1,1:end-1))*dt;
    %      B_w_imag(1:end-1)  = spline(z_val_0(1:end),B_w_imag(1:end),z_val_0(1:end-1)+D(1:end-1))-0.5*mu*v_g(1:end-1).*(JB_g(t_step,1:end-1)).*sin(Phi(t_step-1,1:end-1))*dt; %%%% Spline interpolate onto grid
    %
    toc
    [num2str(100*t_step/M),' % Complete']
end
save([f_save_dir,'.mat'])
% save(['EPIC_stat_TEST_v1_2_',num2str(T_num),'.mat'])
%%%% end particle tracking %%%
