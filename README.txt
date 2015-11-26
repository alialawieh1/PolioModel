This directory contains code and work related to the polio model. 

in order to generate the properties of the individual in Lebanon you need to use these commands in MATLAB:

Indv=[];
[ Indv nHhold] = poliomodel_Leb(Indv);
clc
[ Indv ] = poliomodel_Syr(Indv);
clc
[ Indv  ] = poliomodel_Pal(Indv);
clc

after generating individual properties of people living in Lebanon, here is an example of introduction of two infected syrian cases in North and Bekaa:
use these MATLAB commands to try the code:

CaseIntroPopulationType=[2 2];
caseIntroSites=[3 4];
nbCases=2;
nbDays=20;
[newCases, activeCases, detectedCases, Sites, populations,nbAdults,nbChildren,trackIndv, trackDays] = polioCase(Indv,nHhold, CaseIntroPopulationType, caseIntroSites, nbCases, nbDays);

polioSimulationPost.m contains an example on how to use the post vaccination campaign data