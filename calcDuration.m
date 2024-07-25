function output = calcDuration(prePayTree,shortTree)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

p0 = prePayTree(1,1);
pU = prePayTree(1,2);
pD = prePayTree(2,2);
rU = shortTree(1,2);
rD = shortTree(2,2);

Duration = (-1/p0)*((pU-pD)/(rU-rD));

output = Duration;
end