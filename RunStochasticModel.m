function [ybar,zbar,Ntot] = RunStochasticModel(Nmax,beta,delta,tmax)

Vtype = [10,10]; % variant type (y,z)
Ntype = 1e3;     % number of Vtype host populations

cUp = 1e-3; % probability z = z + 1
cDn = 1e-5; % probability z = z - 1

eUp = 1e-3; % probability y = y + 1
eDn = 1e-5; % probability y = y - 1

ybar = zeros(tmax,1); ybar(1) = dot(Vtype(:,1),Ntype)/sum(Ntype);
zbar = zeros(tmax,1); zbar(1) = dot(Vtype(:,2),Ntype)/sum(Ntype);
Ntot = zeros(tmax,1); Ntot(1) = sum(Ntype);
mutZ = [];

for t = 2:tmax
    
    % births by HGT
    
    birthVtype = [];
    birthNtype = [];
    
    for n = 1:length(Ntype)
        
        betaSTAR = beta*(1 - sum(Ntype)/Nmax);
        births = binornd(sum(poissrnd(betaSTAR,Ntype(n),1)),exp(-Vtype(n,2)/5));
        
        if births > 0
            
            if isempty(birthVtype)
                birthVtype = [birthVtype;[0,Vtype(n,2)]]; %#ok<AGROW>
                birthNtype = [birthNtype;births];         %#ok<AGROW>
            else
                
                idx = find(and(birthVtype(:,1) == 0, birthVtype(:,2) == Vtype(n,2)));
                
                if isempty(idx)
                    birthVtype = [birthVtype;[0,Vtype(n,2)]]; %#ok<AGROW>
                    birthNtype = [birthNtype;births];         %#ok<AGROW>
                else
                    birthNtype(idx) = birthNtype(idx) + births; %#ok<AGROW>
                end
                
            end
            
        end
        
    end
    
    % deaths by gene loss
    
    for n = 1:length(Ntype)
        
        deaths = binornd(Ntype(n),delta*exp(-Vtype(n,1)/5));
        Ntype(n) = Ntype(n) - deaths; %#ok<AGROW>
       
    end
    
    % remove empties lost to death
    
    idx = find(Ntype == 0);
    Vtype(idx,:) = [];
    Ntype(idx) = []; %#ok<AGROW>
    
    % mutations (persistors only)
    
    cneVtype = [];
    cneNtype = [];
    
    for n = 1:length(Ntype)
        
        % indispensability
        vec = mnrnd(Ntype(n),[eDn,1-eDn-eUp,eUp]);
        
        if vec(1) > 0 % y = y - 1 >= 0
            Ntype(n) = Ntype(n) - vec(1);
            cneNtype(length(cneNtype)+1,1) = vec(1);                               %#ok<AGROW>
            cneVtype(size(cneVtype,1)+1,:) = [max([0,Vtype(n,1) - 1]),Vtype(n,2)]; %#ok<AGROW>
        end
        
        if vec(3) > 0 % y = y + 1 <= 20
            Ntype(n) = Ntype(n) - vec(3);
            cneNtype(length(cneNtype)+1,1) = vec(3);                                %#ok<AGROW>
            cneVtype(size(cneVtype,1)+1,:) = [min([50,Vtype(n,1) + 1]),Vtype(n,2)]; %#ok<AGROW>
        end
        
        % connectivity
        vec = mnrnd(Ntype(n),[cDn,1-cDn-cUp,cUp]);
        
        if vec(1) > 0 % z = z - 1 >= 0
            Ntype(n) = Ntype(n) - vec(1);
            cneNtype(length(cneNtype)+1,1) = vec(1);                               %#ok<AGROW>
            cneVtype(size(cneVtype,1)+1,:) = [Vtype(n,1),max([0,Vtype(n,2) - 1])]; %#ok<AGROW>
            mutZ = [mutZ;t,-1]; %#ok<AGROW>
        end
        
        if vec(3) > 0 % z = z + 1 <= 20
            Ntype(n) = Ntype(n) - vec(3);
            cneNtype(length(cneNtype)+1,1) = vec(3);                                %#ok<AGROW>
            cneVtype(size(cneVtype,1)+1,:) = [Vtype(n,1),min([50,Vtype(n,2) + 1])]; %#ok<AGROW>
        end
        
    end
    
    % remove empties lost to mutation
    idx = find(Ntype == 0);
    Vtype(idx,:) = [];
    Ntype(idx) = []; %#ok<AGROW>
   
    % combine arrays
    
    for n = 1:length(birthNtype) % births
        
        [isIN,inIDX] = ismember(birthVtype(n,:),Vtype,'rows');
        
        if isIN == 1
            Ntype(inIDX) = Ntype(inIDX) + birthNtype(n);
        else
            Vtype = [Vtype;birthVtype(n,:)];    %#ok<AGROW>
            Ntype = [Ntype;birthNtype(n)];      %#ok<AGROW>
        end
        
        
    end
    
    for n = 1:length(cneNtype) % mutants
        
        [isIN,inIDX] = ismember(cneVtype(n,:),Vtype,'rows');
        
        if isIN == 1
            Ntype(inIDX) = Ntype(inIDX) + cneNtype(n);
        else
            Vtype = [Vtype;cneVtype(n,:)];    %#ok<AGROW>
            Ntype = [Ntype;cneNtype(n)];      %#ok<AGROW>
        end
        
        
    end
    
    % update means
    
    ybar(t) = dot(Vtype(:,1),Ntype)/sum(Ntype);
    zbar(t) = dot(Vtype(:,2),Ntype)/sum(Ntype);
    Ntot(t) = sum(Ntype);
    
end

%% END