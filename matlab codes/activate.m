function [t, phi] = activate(Tc, dt, cycle)

phi = [];
Ts = 2*Tc/5;
tsys = 0:dt:Ts;
phi1 = [sin(pi*mod(tsys, Tc)/Ts) zeros(1, (Tc-Ts)/dt-1)];
t = 0:dt:cycle*Tc-dt;
for i = 1:cycle;
    phi = [phi phi1];
end
axis auto
grid on
plot(t, phi);
