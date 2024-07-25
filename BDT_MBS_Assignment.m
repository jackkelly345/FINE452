
% JE & EZ - fits BDT tree in Excel Example file

clear all
clc

%% Calibrate Nelson Siegel Function

yieldData=dataset('File','yieldAndVolatility2000.csv','Delimiter',',');

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

%% Set parameters

%End date, in years
T=30;

%Number of steps in the binomial tree
N=360;

%set time step size
dt=T/N;


%m for vector of maturities
m=(dt:dt:T)';

%
observedData=dataset(m);
%% EXCEL REPLICATION
%%Commented out sections below to replicating excel results here.
observedData.price = exp(-fittedFunctionYield(m).*m);
%first cell empty below
observedData.volatility=fittedFunctionVol(m);

%Create empty tree. the (i,j) node represents time i and node j of the
%short rate tree

%create 360x360 matrix, fill lower diagonal of matrix
shortTree=NaN(N,N);

%Starting Guesses for fminsearch
%min diff between first bond price and model bond price
startR=.01;
startParam=[.01,.01];



%% Create first node in tree
i=1;
thisPrice= observedData.price(1);
thisMaturity= 1;

%No volatility matching for first node
%
thisError=@(r) (thisPrice-exp(-r*1/12))^2;

%plug in reason first R
thisR=fminsearch(thisError,startR);

shortTree(i,1)=thisR;


%% Moving on in the tree


for i = 2:N
        i/N;
        %inputs given from observed data
        thisPrice= observedData.price(i);
        thisMaturity=i/12;
        thisVolatility=observedData.volatility(i);

        thisPriceAndVol=[thisPrice,thisVolatility];
        
        %param is starting guess for the tree

        %first build short rate tree, add new slice for time period from
        %building outwards, 

        %then using this data, construct bond tree to find price, plug into
        %fmin search to minize error w market data 

        %build seperate function bondtree, outputs model price and model vol, give it params
        %ie short rate vol and lowest short rate guesses, give it current
        %tree, i to keep track of where it is 

        thisError=@(param) sum((thisPriceAndVol - bondTree(param, i, shortTree)).^2);
        thisParam=fminsearch(thisError,startParam);

        %Build next step in tree with the right parameters
        thisBottomR=thisParam(1);
        thisSigma=thisParam(2);

        %make start param (start guesses) based on previous tree slice 
        startParam = [thisParam(1), thisParam(2)];

        treeStep=exp(2*thisSigma*sqrt(1/12));

        %populate next step in shortTree
        shortTree(i,i)=thisBottomR;
        for j=i-1:-1:1
            shortTree(j,i) = shortTree(j+1,i)*treeStep;
        end

end


%% Value 30 Year MBS with face value $100,000,000, mortgage rate 7%

mtge_rt = 0.07/12;
int_rt = 0.05;
mtge_value = 100000000;

%Create payment schedule

time = 0:360;
principal_outst = zeros(N+1,1);
interest_pmt = zeros(N+1, 1);
principal_pmt = zeros(N+1, 1);
sch_pmt = zeros(N+1, 1);

paymentSchedule = table(time', principal_outst, interest_pmt, principal_pmt, sch_pmt);

paymentSchedule.principal_outst(1) = mtge_value;
paymentSchedule.interest_pmt(1) = 0;
paymentSchedule.principal_pmt(1) = 0;
paymentSchedule.sch_pmt(1) = 0;

fixed_pmt = mtge_value/((1/mtge_rt)*(1-1/(1+mtge_rt)^360));

for i = 2:N+1
    paymentSchedule.interest_pmt(i) =  mtge_rt*paymentSchedule.principal_outst(i-1);
    paymentSchedule.sch_pmt(i) = fixed_pmt;
    paymentSchedule.principal_pmt(i) = paymentSchedule.sch_pmt(i)-paymentSchedule.interest_pmt(i);
    paymentSchedule.principal_outst(i) = paymentSchedule.principal_outst(i-1) - paymentSchedule.principal_pmt(i);

end 

%Create tree for non-pre payable MBS

noPrePayTree = NaN(N+1,N+1);
noPrePayTree(:,N+1) = fixed_pmt*ones(N+1,1);

for col = N:-1:1
    for row = col:-1:1
       noPrePayTree(row, col) = (0.5*noPrePayTree(row+1, col+1) + 0.5*noPrePayTree(row, col+1))*exp(-shortTree(row,col)*1/12);
       if (row ~= 1 | col ~= 1)
       noPrePayTree(row, col) = noPrePayTree(row, col) + fixed_pmt;    
       end 
    end
end

%Create tree for pre payable MBS

prePayTree = NaN(N+1,N+1);
prePayTree(:,N+1) = fixed_pmt*ones(N+1,1);

for col = N:-1:1
    for row = col:-1:1
        contValue = (0.5*prePayTree(row+1, col+1) + 0.5*prePayTree(row, col+1))*exp(-shortTree(row,col)*1/12);
        outstPrincpal = paymentSchedule.principal_outst(col);
        prePayTree(row, col) = min(contValue,outstPrincpal);
       if (row ~= 1 | col ~= 1)
       prePayTree(row, col) = prePayTree(row, col) + fixed_pmt;    
       end 
    end
end



 %% Simulate Different Interest Rates

 % multiply short rate tree so current rates range from 2-10%

 rates = 0.02:0.01:0.1;
 treeMultis = rates./shortTree(1,1);

 %Create table to store MBS valuations, convexities and durations
 noPrePayValues = NaN(9,1);
 prePayValues = NaN(9,1);
 convexitiesNoPre = NaN(9,1);
 convexitiesPrePay = NaN(9,1);
 durationsNoPre = NaN(9,1);
 durationsPrePay = NaN(9,1);

 plotData = table(rates',noPrePayValues,prePayValues, convexitiesNoPre,convexitiesPrePay,durationsNoPre,durationsPrePay);


 for i=1:length(treeMultis)
    %adjust interest rates accordingly
    shortTreeAdj = shortTree.*treeMultis(i);
    
    %create non prepay tree using new rates 
    noPrePayTreeAdj = NaN(N+1,N+1);
    noPrePayTreeAdj(:,N+1) = fixed_pmt*ones(N+1,1);

    for col = N:-1:1
        for row = col:-1:1
           noPrePayTreeAdj(row, col) = (0.5*noPrePayTreeAdj(row+1, col+1) + 0.5*noPrePayTreeAdj(row, col+1))*exp(-shortTreeAdj(row,col)*1/12);
           if (row ~= 1 | col ~= 1)
           noPrePayTreeAdj(row, col) = noPrePayTreeAdj(row, col) + fixed_pmt;    
           end 
        end
    end

    %save valuation for plotting

    plotData.noPrePayValues(i) = noPrePayTreeAdj(1,1);

    %create prepay tree using new rates
    
    prePayTreeAdj = NaN(N+1,N+1);
    prePayTreeAdj(:,N+1) = fixed_pmt*ones(N+1,1);
    
    for col = N:-1:1
        for row = col:-1:1
            contValue = (0.5*prePayTreeAdj(row+1, col+1) + 0.5*prePayTreeAdj(row, col+1))*exp(-shortTreeAdj(row,col)*1/12);
            outstPrincpal = paymentSchedule.principal_outst(col);
            prePayTreeAdj(row, col) = min(contValue,outstPrincpal);
           if (row ~= 1 | col ~= 1)
           prePayTreeAdj(row, col) = prePayTreeAdj(row, col) + fixed_pmt;    
           end 
        end
    end
    
    %save valuation for plotting
    plotData.prePayValues(i) = prePayTreeAdj(1,1);

    plotData.convexitiesNoPre(i) = convexity(noPrePayTreeAdj, shortTreeAdj);
    plotData.durationsNoPre(i) = calcDuration(noPrePayTreeAdj, shortTreeAdj);

    plotData.convexitiesPrePay(i) = convexity(prePayTreeAdj, shortTreeAdj);
    plotData.durationsPrePay(i) = calcDuration(prePayTreeAdj, shortTreeAdj);

 end





subplot(3,2,1)
plot(yieldData.m, yieldData.y,'*',yieldData.m,yieldData.nelsonSiegelYield);
ylabel("Yields")
xlabel("Maturity (Years)")
legend('Observed Yields','Nelson Siegel Fitted Yields')

subplot(3,2,2)
plot(yieldData.m,yieldData.vol,'*',yieldData.m,yieldData.nelsonSiegelVol);
ylabel("Volatilities")
xlabel("Maturity (Years)")
legend('Observed Vols','Nelson Siegel Fitted Vols')

subplot(3,2,3)
plot(plotData.Var1, plotData.prePayValues, plotData.Var1, plotData.noPrePayValues)
ylabel("MBS Valuation ($)")
xlabel("Interest Rates")
legend('With prepayment option','Without prepayment option')

subplot(3,2,4)
plot(plotData.Var1, plotData.convexitiesPrePay, plotData.Var1, plotData.convexitiesNoPre)
ylabel("Convexity")
xlabel("Interest Rates")
legend('With prepayment option','Without prepayment option')

subplot(3,2,5)
plot(plotData.Var1, plotData.durationsPrePay, plotData.Var1, plotData.durationsNoPre)
ylabel("Duration")
xlabel("Interest Rates")
legend('With prepayment option','Without prepayment option')






