% Find the instantaneous flow rate through the mitral valve
% eqn 20
function dQmtdt = instFlow(Ppu, Plv, Lmt, Qmt, Rmt, dAdt, A)
t1 = (Ppu - Plv) / Lmt;
t2 = Qmt * Rmt / Lmt;
t3 = Qmt * dAdt / A;
dQmtdt = t1 - t2 + t3;
end