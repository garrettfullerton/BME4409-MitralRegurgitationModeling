% BME 4409 Final Project
% Modeling Mitral Regurgitation in the Cardiac System
% Authors: Garrett Fullerton, Duc Tran, Andrea Lopez, Jennifer Noa
% Date: 4/1/2022
%
% This paper is an implementation of the Lumped-Parameter Cardiovascular
% System Model described by Paeme et al. in "Mathematical multi-scale model 
% of the cardiovascular system including mitral valve dynamics. Application 
% to ischemic mitral insufficiency."
%% Initialize parameters
clc; clear; close all;
healthy = 1; % healthy = 0 indicates mitral insufficiency
t_max = 20; % seconds
t_step = 0.01; 

% Table 1 values
% TODO: look at table 6 values instead ?
lvf = struct('Ees', 2.8798, 'Vd', 0, 'V0', 0, 'lambda', 0.033, 'P0', 0.1203);
rvf = struct('Ees', 0.5850, 'Vd', 0, 'V0', 0, 'lambda', 0.023, 'P0', 0.2157);
spt = struct('Ees', 48.7540, 'Vd', 2.00, 'V0', 2.00, 'lambda', 0.435, 'P0', 1.1101);
pcd = struct('Ees', nan, 'Vd', nan, 'V0', 200.00, 'lambda', 0.030, 'P0', 0.5003);
vc = struct('Ees', 0.0059, 'Vd', 0);
pa = struct('Ees', 0.3690, 'Vd', 0);
pu = struct('Ees', 0.0073, 'Vd', 0);
ao = struct('Ees', 0.0, 'Vd', 0);

% Table 2 values
mt = struct('R', 0.0158, 'L', 7.6968e-5);
tc = struct('R', 0.0237, 'L', 8.0093e-5);
av = struct('R', 0.0180, 'L', 1.2189e-4);
pv = struct('R', 0.0055, 'L', 1.4868e-4);

% Table 3 values
HR = 60; % bpm
Vtot = 5.5; % L
Pth = -4; % mmHg

switch healthy
    case 1
        Ks = 0.05; % 1/mmHg
        Amax = 1.1; % cm^2
        w = 30; % rad/s
        D = 10; % 1/rad
    case 0
        Ks = 0.05; % 1/mmHg
        Amax = 1.3; % cm^2
        w = 30; % rad/s
        D = 10; % 1/rad
        dc = 0.2;
end


%% Define the hemodynamic model
% This model had too many parameters for us to properly manage using
% MATLAB's differential equation solvers so we instead just solved for the
% values at each time point.

% The following initial conditions need to be figured out...
% Then I think this will work?
% Qpul = 0.001;
% Qmt = 0.001;
% Vpu = pcd.V0;
% Vlv = lvf.V0;
% Ppu = pcd.P0;
% A = Amax;
% dAdt = 0.001;
% Qav = 0.001;
% Pao = spt.P0;
% Qsys = 0.001;
% Vao = spt.V0;
% Qtc = 0.001;
% Vvc = 0.001;
% Vrv = rvf.V0;
% Pvc = 0.001;
% Qpv = 0.001;
% Ppa = 0.001;
% Vpa = 0.001;
% figure; hold on;

for t = 1:t_step:t_max
dVpudt = heaviside(Qpul) * Qpul - Qmt; % eqn 26
Vpu = Vpu + dVpudt * t_step;

[~, ~, ~, Plv] = driver(t, lvf.Ees, Vlv, lvf.Vd, lvf.V0, lvf.P0, lvf.lambda, HR);
dQmtdt = ((1/mt.L)*((Ppu-Plv)-Qmt*mt.R)) + Qmt*A/dAdt; % eqn 27
Qmt = Qmt + dQmtdt * t_step;

dVlvdt = Qmt - heaviside(Qav) * Qav; % eqn 28
Vlv = Vlv + dVlvdt * t_step;

dQavdt = heaviside(heaviside(Plv-Pao)+heaviside(Qav)-0.5) * ((1/av.L)*((Plv-Pao)-Qav*av.R)); %eqn 29
Qav = Qav + dQavdt * t_step;

dVaodt = heaviside(Qav) * Qav - heaviside(Qsys) * Qsys; % eqn 30
Vao = Vao + dVaodt * t_step;

dVvcdt = heaviside(Qsys) * Qsys - heaviside(Qtc) * Qtc; % eqn 31
Vvc = Vvc + dVvcdt * t_step;

[~, ~, ~, Prv] = driver(t, rvf.Ees, Vrv, rvf.Vd, rvf.V0, rvf.P0, rvf.lambda, HR);
dQtcdt = heaviside(heaviside(Pvc-Prv)+heaviside(Qtc)-0.5) * ((1/tc.L)*((Pvc-Prv)-Qtc*tc.R)); % eqn 32
Qtc = Qtc + dQtcdt * t_step;

dVrvdt = heaviside(Qtc) * Qtc - heaviside(Qpv) * Qpv; % eqn 33
Vrv = Vrv + dVrvdt * t_step;

dQpvdt = heaviside(heaviside(Prv-Ppa)+heaviside(Qpv)-0.5) * ((1/pv.L)*((Prv-Ppa)-Qpv*pv.R)); % eqn 34
Qpv = Qpv + dQpvdt * t_step;

dVpadt = heaviside(Qpv) * Qpv - heaviside(Qpul) * Qpul; % eqn 35
Vpa = Vpa + dVpadt * t_step;

switch healthy
    case 1
        % Healthy mitral valve area
        dAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A) - 0.5) * dAdt; % eqn 36
        A = A + dAdt * t_step;
        
        ddAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A) - 0.5) * ((Amax-A)*power(w,2)*Ks*(Ppu-Plv) - 2*D*dAdt*w - A*power(w,2)); % eqn 37
        dAdt = dAdt + ddAdt * t_step;
    case 0
        % Mitral insufficiency
        dAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A-dc) - 0.5) * dAdt; % eqn 38
        A = A + dAdt * t_step;
        
        ddAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A-dc) - 0.5) * ((Amax-A)*power(w,2)*Ks*(Ppu-Plv) - 2*D*dAdt*w - A*power(w,2)); % eqn 39
        dAdt = dAdt + ddAdt * t_step;
end
if isnan(Plv), break; end
scatter(Plv, Vlv)
end






%% Key terms from paper
% lv = left ventricle
% rv = right ventricle
% vc = vena cava
% ao = aorta
% pa = pulmonary artery
% pu = pulmonary veins
% Rsys = systemic circulation network resistance
% Rpul = pulmonary circulation network resistance
% E(t) = elastance = P(t) / V(t)
% mt = mitral valve
% tc = tricuspid valve
% av = aortic valve
% pv = pulmonary valve
% lvf = left ventricle free wall
% rvf = right ventricle ffree wall
% spt = septum free wall
% pcd = pericardium
% HR = heart rate
% Vtot = total blood volume
% Pth = thoracic cavity pressure
% Ks = static gain factor
% Amax = maximum mitral valve area
% w = eigen frequency
% D = damping factor
% dc = defect of closure
% V0 = volume at zero pressure
% Vd = unstressed volume chamber

