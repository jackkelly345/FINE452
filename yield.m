function output = yield(param, m)
%outputs nelson siegel fit of yield for given maturity
output = param(1) + (param(2)*param(3)).*((1-exp(-m./param(4)))./(m./param(4))) - param(3)*exp(-m./param(4));

output;
end