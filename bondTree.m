function output = bondTree(param,i,shortTree)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%now need annualized rate eg divide by 12
%volatility adjusted by sqrt 12 
%tree step needs to be adjusted as well 

%make trees w guesses 
N=i;

guessTree = shortTree;

thisBottomR=param(1);
thisSigma=param(2);

treeStep=exp(2*thisSigma*sqrt(1/12));

%populate next step in shortTree
guessTree(N,N)=thisBottomR;
for index=N-1:-1:1
   guessTree(index,N) = guessTree(index+1,N)*treeStep;
end


%
%need to loop to make whole time slice for larger i
%guessTree(N, N) = param(1);
%for index = 1:N
 %   guessTree(N-index, N) = guessTree(N-index-1, N)*exp(2*param(2));
%end

%Building i-month Bond Tree

NYearBondTree = NaN(N+1,N+1);
NYearBondTree(:,N+1) = ones(N+1,1);

for col = N:-1:1
    for row = col:-1:1
        NYearBondTree(row, col) = (0.5*NYearBondTree(row+1, col+1) + 0.5*NYearBondTree(row, col+1))*exp(-guessTree(row,col)*1/12);
    end
end

modelYield1 = -1/(N-1)*log(NYearBondTree(1,2));
modelYield2 = -1/(N-1)*log(NYearBondTree(2,2));

modelPrice = NYearBondTree(1,1);
modelVol = 0.5*log(modelYield1/modelYield2)*1/sqrt(1/12);
output = [modelPrice modelVol];
end


