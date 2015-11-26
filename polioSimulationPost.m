function [newCases, activeCases, cumulativeCases] = polioSimulationPost( N, CaseIntroPopulationType, caseIntroSites, nbCases, nbDays)
% This code tests the post vaccination campaign new cases, active cases and cumulative cases

% N represents number of iterations. All other parameters are the same as
% in polioCase.m

% each output is a 3-D matrix where the first dimension indicates
%the vaccination level (it is the variable k in this example where k represents one of the vaccination levels "j" in the example below
% the second dimension represents the trial number from 1 to N 
%and the third one the number of cases in the given day (from
%day 0 to nbDays)

k=0;
Newpost=[];
Activepost=[];
Detectedpost=[];
Adultspost=[];
Childrenpost=[];
SiteCasespost=[];
Populationpost=[];

for j=0.80:0.025:1 %vaccination level starting from 80% to 100%
    k=k+1;


for i=1:N
    
 Indv=[];
[Indv nHhold] = poliomodel_Lebpost(Indv,j*ones(6,1));
clc
[Indv] = poliomodel_Syrpost(Indv,j*ones(6,1));
clc
[Indv] = poliomodel_Palpost(Indv,j*ones(6,1));
clc

[newCasesTemp, activeCasesTemp, cumulativeCasesTemp] = polioCaseLight(Indv,nHhold, CaseIntroPopulationType, caseIntroSites, nbCases, nbDays);
clc
newCases(k,i,:)=newCasesTemp;
activeCases(k,i,:)=activeCasesTemp;
cumulativeCases(k,i,:)= cumulativeCasesTemp;
end

end
end