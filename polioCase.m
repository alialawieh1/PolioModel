function [newCases, activeCases, detectedCases, Sites, populations,nbAdults,nbChildren,trackIndv, trackDays ] = polioCase(Indv,nHhold, CaseIntroPopulationType, caseIntroSites, nbCases, nbDays)

% newCases is a vector of length nbDays+1. newCases(i) indicates how many new infected cases we have on day i
% activeCases vector of length nbDays+1, where activeCases(i) indicates how many active infecting cases we have on day i
% detectedCases vector of length nbDays+1, where detectedCases(i) indicates how many detected infecting cases we have till day i
% Sites 2-D matrix of size (nbDays+1)x4 where Site(i,j) represents how many
% new cases we have on day i ,site j
% populations  2-D matrix of size (nbDays+1)x3 where Site(i,j) represents how many
% new cases we have on day i ,population type j
% nbAdults  is a vector of length nbDays+1. nbAdults(i) indicates how many new adults got infected on day i
% nbChildren  is a vector of length nbDays+1. nbChildren(i) indicates how many new children got infected on day i

% trackIndv 2-D matrix of 7 columns indicating for each new infected case the following criteria:
% popType, Site, householdType, which household , which person, isAdult, has access

% trackDays 2-D matrix of 8 columns indicating for each new infected case the following
% criteria: global day since day 0,  days of infection of this indv, latency period, incubation , disease time, symptomsFlag, paralysis flag, detected

% Indv represents studied country people. it is generated using the
% following files: poliomodel_Leb.m, poliomodel_Syr.m, and poliomodel_Pal.m
% -- see README.txt for the generation of Indv and  nHhold parameter.
% CaseIntroPopulationType, caseIntroSites, nbCases, nbDays

%CaseIntroPopulationType is a vector representing the population type of the introduced infected case. the values in the vector are the integers 1,2,or 3 indicating respectively
%lebanese, syrian, palestinian
%caseIntroSites is a vector of integers 1,2,3,4 indicating respectively
%Beirut, Bekaa, North, and South, same length as number of cases. this
%parameter indicate where each infected case entered.

% nbCases is an integer indicating how many infected cases are introduced
% nbDays indicate on how many days we are studying the propagation of
% infection after the introduction of the infected cases.

symptomDist=[zeros(50,1);ones(5,1);zeros(45,1);];      % Symtpoms probability      0.05
paralysisDist=[zeros(595,1); ones(5,1);zeros(400,1)];  % paralysis probability     0.005

mu=[2 10 35 7 0.1 2 2 1 1];
RangL=[0.1 0.1 20 3 0 1 1 0 0];
RangH=[7 20 50 20 0.5 5 4 3 3];
F{1}=[];

for i=1:length(mu)
    Fa= mu(i)+(mu(i)-RangL(i))/2 *randn(1000,1);
    Fb=mu(i)+(RangH(i)-mu(i))/2 *randn(1000,1);
    inda=find(Fa<mu(i) & Fa>=0);
    indb=find(Fb>=mu(i));
    
    F{i}=[Fa(inda);Fb(indb)]; % generate skew normal distributions of mode mu, and lower boundary
    % RangL, and higher boundary RangH for each of the parameters indicated at
    % the end of the code.
end

tempTrackIndv1=[];% 2-D matrix of 7 columns indicating for each new infected case the following criteria:
% popType, Site, householdType, which household , which person, isAdult, has access

tempTrackIndv2=[];
% 2-D matrix of 7 columns indicating for each new infected case the following
% criteria: days of infection, latency period, incubation , disease time, symptomsFlag, paralysis flag, detected

for i=1:nbCases
    sympFlag=randsample(symptomDist,1) ;
    ParalyzedFlag=randsample(paralysisDist,1);
    if ParalyzedFlag==1
        sympFlag=1;
    end
    tempTrackIndv1=[tempTrackIndv1; CaseIntroPopulationType(i) caseIntroSites(i) -1 0 1 0 0]; % the introduced cases are usually children
    tempTrackIndv2=[tempTrackIndv2; 0 round(randsample(F{1},1)) round(randsample(F{2},1)) round(randsample(F{3},1)) sympFlag  ParalyzedFlag 0];
end

communicationPool{1}=[ones(1,4) 2 3]; %each population is more likely to communicate within itself than with others by a ratio of 2/3 internally, and 1/3 with other populations
communicationPool{2}=[ones(1,4)*2 1 3];
communicationPool{3}=[ones(1,4)*3 2 1];
sitesPool{1}=[2 3 4];%each infected individual is more likely to infect people within the same site, but if he infects from other sites the infection probability is 1/3 for each of the three other sites
sitesPool{2}=[1 3 4];
sitesPool{3}=[1 2 4];
sitesPool{4}=[1 2 3];

trackIndv=[];
trackDays=[];
tempIm=0;
tempAdult=0;
tempHasAccess=0;
rateInfection=0;
randInf=0;
trackIndv=[trackIndv;tempTrackIndv1];
trackDays=[trackDays;tempTrackIndv2(1,1)*ones(length(tempTrackIndv2(:,1)),1) tempTrackIndv2];
fprintf('Day %d ***\n',tempTrackIndv2(1,1));
fprintf('Cumulative number of cases %d ***\n',length(tempTrackIndv2(:,1)));

Sites=[];
populations=[];
activeCases(1)=0;

detectedCases(1)=0;
newCases(1)=nbCases;
tempPreviousCases=sum(newCases);
Sites(1,1:4)=zeros(1,4);
populations(1,1:3)=zeros(1,3);
for nb1=1:nbCases
    Sites(1,caseIntroSites(nb1))=Sites(1,caseIntroSites(nb1))+1;
    populations(1,CaseIntroPopulationType(nb1))= populations(1,CaseIntroPopulationType(nb1))+1;
end
nbAdults(1)=0;
nbChildren(1)=nbCases;

while (nbDays>tempTrackIndv2(1,1))
    
    tempTrackIndv2(:,1)=tempTrackIndv2(:,1)+1;
    N1=length(tempTrackIndv1(:,1));
    j1=tempTrackIndv2(1,1);
    Sites(j1+1,:)=zeros(1,4);
    populations(j1+1,:)=zeros(1,3);
    nbAdults(j1+1)=0;
    nbChildren(j1+1)=0;
    fprintf('Day %d ***\n',tempTrackIndv2(1,1));
    for i=1:N1
        %if the patient is within the incubation time, before the end of the
        %disease , and if not paralysed he continues to infect other people.
        flagInfect=0;
        if (tempTrackIndv2(i,1)<tempTrackIndv2(i,2) )
            continue; %if number of infected days< latency period he doesnt infect anyone
        elseif ((tempTrackIndv2(i,1)==tempTrackIndv2(i,2)) &&(i>nbCases)) % after latency period he infects his household if they are not immuned
            for i1=1:tempTrackIndv1(i,3)
                if i1~=tempTrackIndv1(i,5)%if not same individuel
                    
                    tempIm= Indv(tempTrackIndv1(i,2),tempTrackIndv1(i,1),tempTrackIndv1(i,3),tempTrackIndv1(i,4),i1,2);
                    
                    if tempIm<=0
                        if (tempIm==0)%if partially immuned
                            rateInfection=randsample(F{5},1);
                            partialPool=[ones(floor(rateInfection*100),1); zeros(floor((1-rateInfection)*100),1)];
                            if (~isempty(partialPool))
                                
                                randInf=randsample(partialPool,1);
                            end
                        end
                    end
                    if (randInf ||tempIm<0)% if not immuned or if partially immuned and was picked up to be infected
                        
                        tempAdult=Indv(tempTrackIndv1(i,2),tempTrackIndv1(i,1),tempTrackIndv1(i,3),tempTrackIndv1(i,4),i1,1);
                        tempHasAccess=Indv(tempTrackIndv1(i,2),tempTrackIndv1(i,1),tempTrackIndv1(i,3),tempTrackIndv1(i,4),i1,3);
                        if isempty(tempHasAccess)
                            tempHasAccess=0;
                        end
                        tempPop=tempTrackIndv1(i,1);
                        tempSite=tempTrackIndv1(i,2);
                        tempHholdType=tempTrackIndv1(i,3);
                        tempHhold=tempTrackIndv1(i,4);
                        tempInd=tempTrackIndv1(i,5);
                        tempRow=[tempPop tempSite tempHholdType tempHhold tempInd tempAdult tempHasAccess];
                        
                        [ ind1,loc1]=ismember(tempRow,tempTrackIndv1,'rows');%if he was not already infected add him to the infected list
                        if (ind1==0)
                            tempTrackIndv1=[tempTrackIndv1;tempRow];
                            Sites(j1+1, tempSite)=Sites(j1+1, tempSite)+1;
                            populations(j1+1, tempPop)=populations(j1+1, tempPop)+1;
                            
                            nbAdults(j1+1)=nbAdults(j1+1)+tempAdult;
                            nbChildren(j1+1)=nbChildren(j1+1)+ (~tempAdult);
                            sympFlag=randsample(symptomDist,1) ;
                            ParalyzedFlag=randsample(paralysisDist,1);
                            if ParalyzedFlag==1
                                sympFlag=1;
                            end
                            tempTrackIndv2=[tempTrackIndv2; 0 round(randsample(F{1},1)) round(randsample(F{2},1)) round(randsample(F{3},1)) sympFlag  ParalyzedFlag 0];
                            
                        end
                        
                    end
                    
                else
                    continue;
                end
            end
            
        elseif (tempTrackIndv2(i,7)==1)%detected
            continue;%if detected he cannot infect other people
        elseif (tempTrackIndv2(i,1)>tempTrackIndv2(i,3) && tempTrackIndv2(i,1)<tempTrackIndv2(i,4))
            if (tempTrackIndv2(i,6)==1)%paralizedFlag
                continue;%if paralyzed he cannot infect other people
                
            elseif (tempTrackIndv2(i,5)==1 && tempTrackIndv1(i,7)==1)%symptomatic and has access to medical care
                continue;%if symptomatic and has access to medcare he cannot infect other people
            else
                flagInfect=1;
            end
            
            %1- if symptomatic and has access to medical care detected and stops infecting
            %2- if symptomatic and paralysis stops infecting
            %3- else he contiues to infect
            
        elseif (flagInfect==1 || (tempTrackIndv2(i,1)>tempTrackIndv2(i,2) && tempTrackIndv2(i,1)<tempTrackIndv2(i,3)))
            %here the person has infectivity=1 and symptoms=0 and detection=0
            
            
            flagInfect=0;
            
            if tempTrackIndv1(i,6)==1 % check if adult
                nbOnSite= round(randsample(F{7},1));%same province number of people to meet if u r adult
                nbOutSite =round(randsample(F{9},1));%outside province number of people to meet if u r adult
            else
                nbOnSite= round(randsample(F{6},1));%same province number of people to meet if u r child
                nbOutSite =round(randsample(F{8},1));%outside province number of people to meet if u r child
            end         
            tempSite=tempTrackIndv1(i,2);
            for k=1:nbOnSite+nbOutSite
                tempPop= randsample(communicationPool{tempTrackIndv1(i,1)},1);
                if (k>nbOnSite)
                    tempSite=  randsample(sitesPool{tempTrackIndv1(i,2)},1);
                end
                
                tempHholdType = randi([1 8],1,1);
                N=0;
                
                N= nHhold(tempSite,tempPop,tempHholdType);
                
                tempHhold = randi([1 N],1,1);
                tempInd = randi([1 tempHholdType],1,1);
                flagImmune=0;
                tempAdult=0;
                
                flagImmune=Indv(tempSite,tempPop,tempHholdType,tempHhold,tempInd,2);
                tempAdult=Indv(tempSite,tempPop,tempHholdType,tempHhold,tempInd,1);
                tempHasAccess=Indv(tempSite,tempPop,tempHholdType,tempHhold,tempInd,3);
                
                if flagImmune>0
                    continue;
                end
                
                randInf=0;
                if flagImmune==0
                    
                    rateInfection=randsample(F{5},1);
                    partialPool=[ones(floor(rateInfection*100),1); zeros(floor((1-rateInfection)*100),1)];
                    if (~isempty(partialPool))
                        
                        randInf=randsample(partialPool,1);
                    end
                end
                if randInf==1 || flagImmune<0
                    if isempty(tempHasAccess)
                        tempHasAccess=0;
                    end
                    tempRow=[tempPop tempSite tempHholdType tempHhold tempInd tempAdult tempHasAccess];
                    [ ind1,loc1]=ismember(tempRow,tempTrackIndv1,'rows');
                    
                    if (ind1==0)
                        tempTrackIndv1=[tempTrackIndv1;tempRow];
                        Sites(j1+1,tempSite)=Sites(j1+1,tempSite)+1;
                        populations(j1+1, tempPop)=populations(j1+1, tempPop)+1;
                        nbAdults(j1+1)=nbAdults(j1+1)+tempAdult;
                        nbChildren(j1+1)=nbChildren(j1+1)+ (~tempAdult);
                        
                        sympFlag=randsample(symptomDist,1) ;
                        ParalyzedFlag=randsample(paralysisDist,1);
                        if ParalyzedFlag==1
                            sympFlag=1;
                        end
                        tempTrackIndv2=[tempTrackIndv2; 0 round(randsample(F{1},1)) round(randsample(F{2},1)) round(randsample(F{3},1)) sympFlag  ParalyzedFlag 0];
                        
                    end
                    
                end
            end
           
        elseif ((tempTrackIndv2(i,1)>=tempTrackIndv2(i,4))&& (i>nbCases))
            %the person becomes immuned after finishing the disease time
            Indv(tempSite,tempPop,tempHholdType,tempHhold,tempInd,2)=1;
            continue;
        end
        
    end
    
    trackIndv=tempTrackIndv1;
    trackDays=[tempTrackIndv2(1,1)*ones(length(tempTrackIndv2(:,1)),1) tempTrackIndv2];
    fprintf('Cumulative number of cases %d ***\n',length(tempTrackIndv2(:,1)));
    index1=find(trackDays(:,1)==tempTrackIndv2(1,1),1,'last');
    
    detectedCases(j1+1)=sum(trackDays(1:index1,8))+sum(detectedCases(1:j1));
    tempActive=[];
    for j=1:length(tempTrackIndv2(:,1))
        if ((tempTrackIndv2(j,1)>=tempTrackIndv2(j,2)) &&(tempTrackIndv2(j,1)<=tempTrackIndv2(j,4)))
            tempActive=[tempActive;j];         
        end
    end
    
    activeCases(j1+1)=length(tempActive);
    [index2]=find(trackDays(1:index1,2)==0); 
    newCases(j1+1)=numel(index2);
    tempPreviousCases=sum(newCases);
end

end
%%

%%

% Feature                   Mean	RangL	RangH
% 1-Latency period            2       0.1     7					# infected with time
% 2-Incubation                10      0.1     20					# detected with time	Symp=1	Access=1
% 3-Disease (Not Immune)      35      20      50					# paralyzed with time
% 4-Disease (Partial Immune)	7       3       20					Day of alarm
% Household infection rate
% 	Immune                  0       0       0
% 	5-Partial                 0.1     0       0.5
% 	Non                     1       1       1
% Contacts in province
% 	6-child                   2       1       10  67% same country, 33% other country..
% 	7-Adult                   2       1       4
% Contacts outside province
% 	8-child                   1       0       3
% 	9-Adult                   1       0       3
% Symtpoms probability      0.05
% paralysis                 0.005
% Access to Medcare
% 	Syrian                  75%
% 	Pal                     75%
% 	Leb                     98%
%%
%if paralyzed he stop infecting people. and paralysis may happen after the
%incubation period, not before it
%%