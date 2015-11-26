function [ Indv nHhold] = poliomodel_Leb(Indv)

% Indv as input should be empty at the first place and after running the
% code the lebanese population with all its characteristics will be present
% in the output Indv parameter.

% the output parameter Indv is 6-D matrix, where for Indv(i,j,k,z,m,n) i
% represents the site, j the population type, k the household type, z which
% household, m which person in this household, and n property of this
% person, if 1 it refers to the propertie isAdult, 2 refers to the
% immunization property, and 3 refers to the access to medcare property.
% for any of these values of n if the value in Indv(i,j,k,z,m,n) is one it
% means that the individual has this property.

%nHhold: returns the characteristics of the households for each population
%type, at each site.
% it is a 3-D matrix nHhold(i,j,k) where the value represents how many people are in a given site i, given population type j, and given household type k.
% So i represents the site (from 1 to 4 respectively representing Beirut, Bekaa, North, and South.
% j represents the population type (1 to 3 representing respectively
% Lebanese, Syrian, and Palastinian populations)
% and k represent the type of household, from 1 to 8 indicating how many
% people are in this household.

popLeb = 5882562; % popLeb : total number lebanese; Same for popSyr and popPal
popSyr = 1435840;
popPal = 455000;

percentLeb=[0.499 0.129 0.216 0.157]; %percentage of lebanese people in each site
percentSyr=[0.295 0.351 0.247 0.118];
percentPal=[0.164 0.077 0.21 0.549];
for k=1:4
    
    nbLeb(k) = percentLeb(k)*popLeb; %number of lebanese people in each site
    nbSyr(k) = percentSyr(k)*popSyr;
    nbPal(k) = percentPal(k)*popPal;
end


averageHhold=[3.99 4.58 4.74 4.384];%average of household in each site
percentHhold=[7.3 14.5 15.2 19.7 18.5 12.2 6.4 6.2]/100; %percentage of households having 1 to 8 people at home..
for k=1:4
    countLeb=0;
    countSyr=0;
    countPal=0;
    for j=1:8
        nHhold(k,1,j)= round(percentHhold(j)*nbLeb(k)/averageHhold(k));
        countLeb = countLeb+nHhold(k,1,j)*j;
        nHhold(k,2,j)= round(percentHhold(j)*nbSyr(k)/averageHhold(k));
        countSyr = countSyr+nHhold(k,2,j)*j;
        nHhold(k,3,j)= round(percentHhold(j)*nbPal(k)/averageHhold(k));
        countPal = countPal + nHhold(k,3,j)*j;
    end
    TotalInd(k,1)= countLeb;%number of individual in each site
    TotalInd(k,2)= countSyr;
    TotalInd(k,3)= countPal;
end

ImmSyrAdult= 0.91; % percentage of immunized syrian adults
ImmLebAdult= 0.95;
ImmPalAdult= 0.93;

ImmSyrChild= 0.325;
ImmLebChild= 0.95;
ImmPalChild= 0.63;

PImmSyrAdult= 0;% PI stands for partially immune
PImmLebAdult= 0;
PImmPalAdult= 0;

PImmSyrChild= 0.475;
PImmLebChild= 0;
PImmPalChild= 0.125;


percentadultleb = 0.627;
percentadultSyr = 0.468;
percentadultPal = 0.672;
N=[];

defaultImmAdult=0;

if ImmLebAdult>0.5% set the default based on the percntage of immunized. so that we optimize how to set the properties for each person in the population
    defaultImmAdult=1;
else
    defaultImmAdult=-1;
end

defaultImmChild=0;

if ImmLebChild>0.5
    defaultImmChild=1;
else
    defaultImmChild=-1;
end

defaultAccess=1;
%%
for i=1:4%number of sites
    for k1=1:8 %number of household types
        
        N(k1)=nHhold(i,1,k1);
    end
    
    TotalAdult(i,1) = percentadultleb*TotalInd(i,1);
    
    TotalChild(i,1) = (1 - percentadultleb)*TotalInd(i,1);
    
    Nadult=TotalAdult(i,1);
    NadultRemPerHholdType=floor((Nadult-(sum(N)+(N(2)*(0.85))))/6);
    
    for k=1:8
        
        fprintf('Site %d household type %d. ***\n',i,k);
        switch k
            case 1
                % ALL individuals in single person hholds are adults as per UNDP:
                % Demographic characteristics of residents page 4
                Indv(i,1,k,1:N(k),1,1)= 1;%isAdult
                Indv(i,1,k,1:N(k),1,2)= defaultImmAdult;%isImmune
                Indv(i,1,k,1:N(k),1,3) = defaultAccess;%access to health care
            case 2
                % first person in a household of two people is adult, the
                % second in in 85% of cases adult
                Indv(i,1,k,1:N(k),1,1)= 1;
                Indv(i,1,k,1:N(k),1,2) = defaultImmAdult;
                Indv(i,1,k,1:N(k),1,3)= defaultAccess;
                
                Indv(i,1,k,1:floor(N(k)*(0.85)),2,1)=1;
                Indv(i,1,k,1:floor(N(k)*(0.85)),2,2)=defaultImmAdult;
                
                Indv(i,1,k,floor(N(k)*(0.85)):N(k),2,1)=0;
                Indv(i,1,k,floor(N(k)*(0.85)):N(k),2,2)=defaultImmChild;
                Indv(i,1,k,floor(N(k)*(0.85)):N(k),2,3)=defaultAccess;
                
            otherwise
                TempMatrix=[];
                % first person is adult always
                Indv(i,1,k,1:N(k),1,1)=1;
                Indv(i,1,k,1:N(k),1,2)= defaultImmAdult;
                Indv(i,1,k,1:N(k),1,3)= defaultAccess;
                % we put the default to be a child, then we pick randomly
                % those who are adult based on how many adults are left to be assigned.
                Indv(i,1,k,1:N(k),2:k,1)=0;
                Indv(i,1,k,1:N(k),2:k,2)=defaultImmChild;
                Indv(i,1,k,1:N(k),2:k,3)=defaultAccess;
                
                RandInd = randi([2 k],3*NadultRemPerHholdType,1);
                RandHhold = randi([1 N(k)],3*NadultRemPerHholdType,1);
                
                TempMatrix=[TempMatrix;RandHhold RandInd];
                TempMatrix= unique(TempMatrix,'rows');
                TempMatrix=TempMatrix(1:NadultRemPerHholdType,:);
                
                for k2=1:NadultRemPerHholdType
                    Indv(i,1,k,TempMatrix(k2,1),TempMatrix(k2,2),1)=1;
                    Indv(i,1,k,TempMatrix(k2,1),TempMatrix(k2,2),2)= defaultImmAdult;
                    Indv(i,1,k,TempMatrix(k2,1),TempMatrix(k2,2),3)= defaultAccess;
                end
        end
        
    end
end

fprintf('finalizing 1/5.... Please wait ***\n');
N11=sum(TotalInd(:,1));
NbImmLebAdults=ImmLebAdult*percentadultleb*N11;
NbPImmLebAdults=PImmLebAdult*percentadultleb*N11;
NbNoImmLebAdults=(1-ImmLebAdult-PImmLebAdult)*percentadultleb*N11;

NbImmLebChild=ImmLebChild*(1-percentadultleb)*N11;
NbPImmLebChild=PImmLebChild*(1-percentadultleb)*N11;
NbNoImmLebChild=(1-ImmLebChild-PImmLebChild)*(1-percentadultleb)*N11;


chooseSite=[1 1 1 2 2 3 3 4 4];
%%assign immunization of each individual

Nb=0;
if defaultImmAdult==1
    Nb=NbNoImmLebAdults;
else
    Nb=NbImmLebAdults;
end

for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),1,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            
            if (Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,1) == 1 && Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2) == defaultImmAdult)
                Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2)=-defaultImmAdult;
                flag1=1;
                break;
            end
        end
    end
end

fprintf('finalizing 2/5.... Please wait ***\n');

Nb=NbPImmLebAdults;
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),1,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,1) == 1 && Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2) == defaultImmAdult)
                Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2)=0;
                flag1=1;
                break;
            end
        end
    end
end

Nb=0;
if defaultImmChild==1
    Nb=NbNoImmLebChild;
else
    Nb=NbImmLebChild;
end
fprintf('finalizing 3/5.... Please wait ***\n');

for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([2 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),1,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,1) == 0 && Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2) == defaultImmChild)
                Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2)=-defaultImmChild;
                flag1=1;
                break;
            end
        end
    end
end
Nb=NbPImmLebChild;
fprintf('finalizing 4/5.... Please wait ***\n');
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([2 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),1,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,1) == 0 && Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2) == defaultImmChild)
                Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,2)=0;
                flag1=1;
                break;
            end
        end
    end
end
%% assign medcare access to each individual
NbLebNoAccess=0.02*N11;
Nb=NbLebNoAccess;
fprintf('finalizing 5/5.... Please wait ***\n');
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),1,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if ( Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,3)== defaultAccess)
                Indv(chooseSite(RandIndSite),1,RandHhold,N2,nbIndv,3)=0;
                flag1=1;
                break;
            end
        end
    end
end
fprintf('Done......***\n');
end