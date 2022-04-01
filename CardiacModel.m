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
dVpudt = heaviside(Qpul) * Qpul - Qmt; % eqn 26
dQmtdt = ((1/Lmt)*((Ppu-Plv)-Qmt*Rmt)) + Qmt*A/dAdt; % eqn 27
dVlvdt = Qmt - heaviside(Qav) * Qav; % eqn 28
dQavdt = heaviside(heaviside(Plv-Pao)+heaviside(Qav)-0.5) * ((1/Lav)*((Plv-Pao)-Qav*Rav)); %eqn 29
dVaodt = heaviside(Qav) * Qav - heaviside(Qsys) * Qsys; % eqn 30
dVvcdt = heaviside(Qsys) * Qsys - heaviside(Qtc) * Qtc; % eqn 31
dQtcdt = heaviside(heaviside(Pvc-Prv)+heaviside(Qtc)-0.5) * ((1/Ltc)*((Pvc-Prv)-Qtc*Rtc)); % eqn 32
dVrvdt = heaviside(Qtc) * Qtc - heaviside(Qpv) * Qpv; % eqn 33
dQpvdt = heaviside(heaviside(Prv-Ppa)+heaviside(Qpv)-0.5) * ((1/Lpv)*((Prv-Ppa)-Qpv*Rpv)); % eqn 34
dVpadt = heaviside(Qpv) * Qpv - heaviside(Qpul) * Qpul; % eqn 35

switch healthy
    case 1
        % Healthy mitral valve area
        dAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A) - 0.5) * dAdt; % eqn 36
        ddAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A) - 0.5) * ((Amax-A)*power(w,2)*Ks*(Ppu-Plv) - 2*d*dAdt*w - A*power(w,2)); % eqn 37
    case 0
        % Mitral insufficiency
        dAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A-dc) - 0.5) * dAdt; % eqn 38
        ddAdt = heaviside(heaviside(Ppu-Plv) + heaviside(A-dc) - 0.5) * ((Amax-A)*power(w,2)*Ks*(Ppu-Plv) - 2*d*dAdt*w - A*power(w,2)); % eqn 39
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

