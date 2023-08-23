function [df] = calc_delta(func_str, vars_struct)
%{ 
usage example: 
to calculate the RSS (root-sum-square) of all error factors 
of the potential energy U(x, y) = (x^2)*y + x/y, enter:
calc_delta('(x^2)*y + x/y',  
    struct('names', {'x', 'y'}, 'values', {2, 3}, 'deltas', {0.1, 0.01}))
%}
func = str2sym(func_str);
sum = 0;
for i = 1:length(vars_struct)
    dx = vars_struct(i).deltas;
    if dx <= 0
        continue;
    end
    x = str2sym(vars_struct(i).names);
    dfdx = diff(func, x);
    diff_val = vpa(subs(dfdx, {vars_struct.names}, {vars_struct.values}));
    sum = sum + (diff_val * dx) ^ 2;
end
df = sqrt(sum);
end

