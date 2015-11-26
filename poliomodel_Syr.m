function [ Indv ] = poliomodel_Syr(Indv)
%same as we did for poliomodel_Leb but now for the syrians. the input
%Indv is now filled for Leb population

popLeb = 5882562;
popSyr = 1435840;
popPal = 455000;

percentLeb=[0.499 0.129 0.216 0.157];
percentSyr=[0.295 0.351 0.247 0.118];
percentPal=[0.164 0.077 0.21 0.549];
for k=1:4
    
    nbLeb(k) = percentLeb(k)*popLeb;
    nbSyr(k) = percentSyr(k)*popSyr;
    nbPal(k) = percentPal(k)*popPal;
end

averageHhold=[3.99 4.58 4.74 4.384];
percentHhold=[7.3 14.5 15.2 19.7 18.5 12.2 6.4 6.2]/100;
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
    TotalInd(k,1)= countLeb;
    TotalInd(k,2)= countSyr;
    TotalInd(k,3)= countPal;
end

ImmSyrAdult= 0.91;
ImmLebAdult= 0.95;
ImmPalAdult= 0.93;

ImmSyrChild= 0.325;
ImmLebChild= 0.95;
ImmPalChild= 0.63;

PImmSyrAdult= 0;
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
if ImmSyrAdult>0.5
    defaultImmAdult=1;
else
    defaultImmAdult=-1;
end
defaultImmChild=0;
if ImmSyrChild>0.5
    defaultImmChild=1;
else
    defaultImmChild=-1;
end
defaultAccess=1;



for i=1:4
    for k1=1:8
        
        N(k1)=nHhold(i,2,k1);
    end
    
    TotalAdult(i,2) = percentadultSyr*TotalInd(i,2);
    
    TotalChild(i,2) = (1 - percentadultSyr)*TotalInd(i,2);
    
    Nadult=TotalAdult(i,2);
    NadultRemPerHholdType=floor((Nadult-(sum(N)+(N(2)*(0.85))))/6);
    
    for k=1:8
        
        fprintf('Site %d household type %d. ***\n',i,k);
        switch k
            case 1
                Indv(i,2,k,1:N(k),1,1)= 1;
                Indv(i,2,k,1:N(k),1,2)= defaultImmAdult;
                Indv(i,2,k,1:N(k),1,3) = defaultAccess;
                
            case 2
                
                Indv(i,2,k,1:N(k),1,1)= 1;
                Indv(i,2,k,1:N(k),1,2) = defaultImmAdult;
                Indv(i,2,k,1:N(k),1,3)= defaultAccess;
                
                Indv(i,2,k,1:floor(N(k)*(0.85)),2,1)=1;
                Indv(i,2,k,1:floor(N(k)*(0.85)),2,2)=defaultImmAdult;
                
                
                Indv(i,2,k,floor(N(k)*(0.85)):N(k),2,1)=0;
                Indv(i,2,k,floor(N(k)*(0.85)):N(k),2,2)=defaultImmChild;
                Indv(i,2,k,floor(N(k)*(0.85)):N(k),2,3)=defaultAccess;
                L1= length(floor(N(k)*(0.85)):N(k));
                
            otherwise
                TempMatrix=[];
                
                
                Indv(i,2,k,1:N(k),1,1)=1;
                Indv(i,2,k,1:N(k),1,2)= defaultImmAdult;
                Indv(i,2,k,1:N(k),1,3)= defaultAccess;
                
                Indv(i,2,k,1:N(k),2:k,1)=0;
                Indv(i,2,k,1:N(k),2:k,2)=defaultImmChild;
                Indv(i,2,k,1:N(k),2:k,3)=defaultAccess;
                
                RandInd = randi([2 k],3*NadultRemPerHholdType,1);
                RandHhold = randi([1 N(k)],3*NadultRemPerHholdType,1);
                
                
                TempMatrix=[TempMatrix;RandHhold RandInd];
                TempMatrix= unique(TempMatrix,'rows');
                TempMatrix=TempMatrix(1:NadultRemPerHholdType,:);
                for k2=1:NadultRemPerHholdType
                    Indv(i,2,k,TempMatrix(k2,1),TempMatrix(k2,2),1)=1;
                    Indv(i,2,k,TempMatrix(k2,1),TempMatrix(k2,2),2)= defaultImmAdult;
                    Indv(i,2,k,TempMatrix(k2,1),TempMatrix(k2,2),3)= defaultAccess;
                end
                
                
        end
        
    end
end

N11=sum(TotalInd(:,2));
NbImmSyrAdults=ImmSyrAdult*percentadultSyr*N11;
NbPImmSyrAdults=PImmSyrAdult*percentadultSyr*N11;
NbNoImmSyrAdults=(1-ImmSyrAdult-PImmSyrAdult)*percentadultSyr*N11;

NbImmSyrChild=ImmSyrChild*(1-percentadultSyr)*N11;
NbPImmSyrChild=PImmSyrChild*(1-percentadultSyr)*N11;
NbNoImmSyrChild=(1-ImmSyrChild-PImmSyrChild)*(1-percentadultSyr)*N11;


chooseSite=[1 1 1 2 2 3 3 4 4];
Nb=0;
if defaultImmAdult==1
    Nb=NbNoImmSyrAdults;
else
    Nb=NbImmSyrAdults;
end
fprintf('finalizing 1/5.... Please wait ***\n');
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),2,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            
            if (Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,1) == 1 && Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2) == defaultImmAdult)
                Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2)=-defaultImmAdult;
                flag1=1;
                break;
            end
        end
    end
end

fprintf('Finalizing 2/5.... Please wait ***\n');
Nb=NbPImmSyrAdults;
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),2,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,1) == 1 && Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2) == defaultImmAdult)
                Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2)=0;
                flag1=1;
                break;
            end
        end
    end
end

Nb=0;
if defaultImmChild==1
    Nb=NbNoImmSyrChild;
else
    Nb=NbImmSyrChild;
end
fprintf('Finalizing 3/5.... Please wait ***\n');

for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([2 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),2,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,1) == 0 && Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2) == defaultImmChild)
                Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2)=-defaultImmChild;
                flag1=1;
                break;
            end
        end
    end
end


fprintf('Finalizing 4/5.... Please wait ***\n');
Nb=NbPImmSyrChild;
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([2 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),2,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if (Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,1) == 0 && Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2) == defaultImmChild)
                Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,2)=0;
                flag1=1;
                break;
            end
        end
    end
end

NbSyrNoAccess=0.02*N11;

Nb=NbSyrNoAccess;
fprintf('Finalizing 5/5.... Please wait ***\n');
for i=1:Nb
    flag1=0;
    while (flag1==0)
        RandIndSite = randi(9,1);
        RandHhold = randi([1 8],1,1);
        
        N1=nHhold(chooseSite(RandIndSite),2,RandHhold);
        N2=randi([1 N1],1,1);
        
        for nbIndv=1:RandHhold
            if ( Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,3)== defaultAccess)
                Indv(chooseSite(RandIndSite),2,RandHhold,N2,nbIndv,3)=0;
                flag1=1;
                break;
            end
        end
    end
end
fprintf('Done......***\n');
end