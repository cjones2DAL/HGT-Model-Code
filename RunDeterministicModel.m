function [n,exT] = RunDeterministicModel(Nmax,beta,delta,tmax)

% **************************************************
% wild type
% **************************************************
pD(1) = exp(-10/5); pB(1) = exp(-10/5); % parent
pD(2) = 1;          pB(2) = exp(-10/5); % offspring
% **************************************************

%%  run system

n = zeros(tmax,2);
n(1,1) = 1000;
n(1,2) = 0;

q = zeros(tmax,2);
q(1,1) = n(1,1)/sum(n(1,:));
q(1,2) = n(1,2)/sum(n(1,:));

exT = nan*ones(1,2);

for t = 2:tmax
    
    % fitnesses
    
    Ntot = sum(n(t-1,:));
    
    wP = zeros(1,2); % persistence
    wM = zeros(1,2); % multiplication
    for i = 1:2
        wP(i) = (1 - delta*pD(i));
        wM(i) = beta*(1-Ntot/Nmax)*pB(i);
    end
    
    n(t,1) = n(t-1,1)*wP(1);
    n(t,2) = n(t-1,1)*wM(1) + n(t-1,2)*(wP(2) + wM(2));

    if isnan(exT(1))
        if n(t,1) < 1
            exT(1) = t;
        end
    end
    
     if isnan(exT(2))
        if n(t,1) + n(t,2) < 1
            exT(2) = t;
        end
    end
    
    wbar = dot(q(t-1,:),wP + wM);
    
    q(t,1) = q(t-1,1)*wP(1)/wbar;
    q(t,2) = q(t-1,1)*wM(1)/wbar + q(t-1,2)*(wP(2) + wM(2))/wbar;
    
end


%% END