% Defining driver function
% Used to define the pressure-volume relationship during filling and
% ejection

function [Pes, Ped, e, Pt] = driver(t, Ees, V, Vd, V0, P0, lambda, HR)
% eqn 2
% end systolic pressure-volume relationship
Pes = Ees * (V - Vd);

% eqn 3
% end diastolic pressure-volume relationship
Ped = P0 * (exp(lambda * (V-V0)) - 1);

% eqn 4
% driver function, below parameters defined by Smith et al.
% Smith BW, Chase JG, Shaw GM, Nokes RI: Experimentally verified minimal 
% cardiovascular system model for rapid diagnostic assistance. Control 
% Engineering Practice 2005, 13:1183-1193.
% ----------
% number of gaussians N = 1
% magnitude A_1 = 1
% width B_1 = 80 s^-2
% delay C_1 = 0.27 s
% cardiac cycle duration D_1 = 1 s (D_1 = 60/HR)
% ----------
A_1 = 1;
B_1 = 80;
C_1 = 0.27;
D_1 = 60 / HR;
e = A_1 * exp(-B_1 * power(((mod(t, D_1))-C_1),2));

% eqn 5
% pressure in the ventricle at any time t of a cardiac cycle
Pt = (e * Ees * (V - Vd)) + ((1 - e) * P0 * (exp(lambda * (V - V0)) - 1));
end
