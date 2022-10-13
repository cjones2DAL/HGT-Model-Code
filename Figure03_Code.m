% Stochastic HGT Model
% Used to genererate Fig. 3
% 25 April 2022
% by C T Jones

% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:

clc
clear
close all

delta = 0.01;
tmax = 1e4;

Nmax = 1e4;
beta = 0.02:0.02:0.12;

nTrials = 50;

Y = nan*ones(nTrials,length(beta));
Z = nan*ones(nTrials,length(beta));
N = nan*ones(nTrials,length(beta));

Ybar = nan*ones(length(Nmax),length(beta));
Zbar = nan*ones(length(Nmax),length(beta));
Nbar = nan*ones(length(Nmax),length(beta));
Pext = nan*ones(length(Nmax),length(beta));

for a = 1:length(Nmax)
    
    for b = 1:length(beta)
        
        disp(['Nmax = ' num2str(Nmax(a)) ', beta = ' num2str(beta(b))])
        
        tic
        for n = 1:nTrials
            
            if mod(n,nTrials/10) == 0
                disp([num2str(n) '/' num2str(nTrials)])
            end
            
            [ybar,zbar,Ntot] = RunStochasticModel(Nmax(a),beta(b),delta,tmax);
            
            Y(n,b) = ybar(end);
            Z(n,b) = zbar(end);
            N(n,b) = Ntot(end);
            
        end
        toc
        
    end
    
    for b = 1:length(beta)
        
        idx = find(N(:,b) > 0);
        
        Pext(a,b) = (nTrials - length(idx))/nTrials;
        
        Ybar(a,b) = mean(Y(idx,b),'omitnan');
        Zbar(a,b) = mean(Z(idx,b),'omitnan');
        Nbar(a,b) = mean(N(idx,b),'omitnan');
         
    end
    
end

%% generate figure

figure(1),clf

hold on
plot(beta,Ybar,'o-','LineWidth',2,'Color',0.60*ones(1,3),'LineWidth',2,'MarkerSize',12)
plot(beta,Zbar,'o-','LineWidth',2,'Color',0.00*ones(1,3),'LineWidth',2,'MarkerSize',12)
hold off

set(gca,'XTick',beta,'XTickLabel',beta,'YTick',[0,10,20])
xlabel('$\beta$','Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex','FontSize',20,'LineWidth',2)
grid on, box on

LG = legend('mean $\bar{y}$','mean $\bar{z}$');
LG.Location = 'northeast';
LG.Interpreter = 'Latex';
LG.FontSize = 20;
grid on, box on

set(gcf,'Pos',[100,200,700,300])

%% END