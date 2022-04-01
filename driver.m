% Defining driver function
% Used to define the pressure-volume relationship during filling and
% ejection

function [Pes, Ped, e] = driver(Ees, V, Vd, V0, P0, lambda, HR, C)
% eqn 2
% end systolic pressure-volume relationship
Pes = Ees * (V - Vd);

% eqn 3
% end diastolic pressure-volume relationship
Ped = P0 * (exp(lambda * (V-V0)) - 1);

% eqn 4
% driver function

end
