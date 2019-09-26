function [new_Rp] = generate_Rp(c,reflection_coefficient)
% Authors: Xihao Li, Catherine Kocia, Michael Chen
c = c*1000;

new_Rp = (-c-reflection_coefficient*c)/(reflection_coefficient-1)*0.5;

end