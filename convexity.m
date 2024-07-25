function output = convexity(prePayTree,shortTree)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

p0 = prePayTree(1,1);
pU = prePayTree(1,2);
pD = prePayTree(2,2);
rU = shortTree(1,2);
rD = shortTree(2,2);

Duration = (-1/p0)*((pU-pD)/(rU-rD));

pUU = prePayTree(1,3);
pUD = prePayTree(2,3);
rUU = shortTree(1,3);
rUD = shortTree(2,3);

DU = (-1/pU)*((pUU-pUD)/(rUU-rUD));

pDD = prePayTree(1,3);
rDD = shortTree(1,3);

DD = (-1/pD)*((pUD-pDD)/(rUD-rDD));

Convexity = Duration^2 - (DU - DD)/(rU-rD);

output = Convexity;

end