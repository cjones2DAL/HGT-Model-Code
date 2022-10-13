% Stochastic HGT Model
% Used to genererate Fig. 2
% 25 April 2022
% by C T Jones

% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:

%% initialize

clc
clear
close all

% ***********************************
% variable parameters
% ***********************************

delta = 0.01;     % gene-loss parameter
beta1 = 0.06;     % birth-rate parameter, population 1
beta2 = 0.08;     % birth-rate parameter, population 2
Nmax = 1e4;       % max size of the metapopulation
tmax = 2e4;       % the number of mappings

%% run system

tic
disp(['beta = ' num2str(beta1)])
[ybar1,zbar1,Ntot1] = RunStochasticModel(Nmax,beta1,delta,tmax);
toc

tic
disp(['beta = ' num2str(beta2)])
[ybar2,zbar2,Ntot2] = RunStochasticModel(Nmax,beta2,delta,tmax);
toc

%% generate two-panel figure

tStr1 = ['$\beta$ = ' num2str(beta1) ', ' num2str(sum(Ntot1(end))) ' populations'];

t = length(Ntot1);

figure(1),clf
subplot(211)

set(gca,'TicklabelInterpreter','Latex','FontSize',20,'LineWidth',1)

hold on
plot(1:t,ybar1(1:t),'Color',0.50*ones(1,3),'LineWidth',2)
plot(1:t,zbar1(1:t),'Color',0.00*ones(1,3),'LineWidth',2)
hold off

set(gca,'YLim',[0,35]);
text(400,19,'a','FontSize',30,'Interpreter','Latex')
title(tStr1,'Interpreter','Latex','FontSize',20)
LG = legend('$\bar{y}$','$\bar{z}$');
LG.Location = 'southeast';
LG.Interpreter = 'Latex';
LG.FontSize = 15;
grid on, box on

tStr2 = ['$\beta$ = ' num2str(beta2) ', ' num2str(sum(Ntot2(end))) ' populations'];

subplot(212)

set(gca,'TicklabelInterpreter','Latex','FontSize',20,'LineWidth',1)

hold on
plot(1:t,ybar2(1:t),'Color',0.50*ones(1,3),'LineWidth',2)
plot(1:t,zbar2(1:t),'Color',0.00*ones(1,3),'LineWidth',2)
hold off

set(gca,'YLim',[0,15]);
xlabel('mappings','Interpreter','Latex')
text(400,12.5,'b','FontSize',30,'Interpreter','Latex')
title(tStr2,'Interpreter','Latex','FontSize',20)
LG = legend('$\bar{y}$','$\bar{z}$');
LG.Location = 'northeast';
LG.Interpreter = 'Latex';
LG.FontSize = 15;
grid on, box on


set(gcf,'Pos',[200,80,650,500])

%% END


