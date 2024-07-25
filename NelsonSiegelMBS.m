%import data for the treasury yield curve

yieldData=dataset('File','yieldAndVolatility2000.csv','Delimiter',',')

%par vector contains, respectively, beta0, beta1, beta2, and tau

startPar=[1,1,1,1];

m = yieldData.m;

%create error functions between nelsonsiegel yield and vol and observed
%data, 

%Alias the error function with observed data and parameters.

errorNelsonSiegelYield=@(par) sum((yield(par,m)-yieldData.y).^2);
errorNelsonSiegelVol=@(par) sum((volatility(par,m)-yieldData.vol).^2);


%Run optimization
[parMinYield,errorMin] = fminsearch(errorNelsonSiegelYield, startPar);
[parMinVol,errorMin]   = fminsearch(errorNelsonSiegelVol, startPar);


%Alias new function fittedYield that takes as input some maturity m and 
%returns fitted yields at the estimated parameters
fittedFunctionYield= @(m) yield(parMinYield,m);
fittedFunctionVol  = @(m) volatility(parMinVol,m);

%added fitted values to table
yieldData.nelsonSiegelYield= fittedFunctionYield(m);
yieldData.nelsonSiegelVol= fittedFunctionVol(m);

subplot(2,1,1)
plot(yieldData.m,yieldData.y,'*',yieldData.m,yieldData.nelsonSiegelYield);
legend('Observed Yields','Nelson Siegel Fitted Yields')
subplot(2,1,2)
plot(yieldData.m,yieldData.vol,'*',yieldData.m,yieldData.nelsonSiegelVol);
legend('Observed Vols','Nelson Siegel Fitted Vols')


