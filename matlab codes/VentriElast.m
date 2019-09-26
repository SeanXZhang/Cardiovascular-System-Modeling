function [t, En] = VentriElast(Tc, dt, cycle)

t1 = 0:dt:Tc-dt;
t = 0:dt:cycle*Tc-dt;
tn = t1/(0.2+0.1555*Tc);
E1 = 1.553174 * (tn/0.7).^1.9 ./ (1+(tn/0.7).^1.9) ./ (1+(tn/1.173474).^21.9);
En = [];
for i = 1:cycle;
    En = [En E1];
end
plot(t, En)
