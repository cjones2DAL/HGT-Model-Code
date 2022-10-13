% Deterministic HGT Model
% Used to genererate Fig. 1
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

delta = 0.01;     % gene-loss parameter
beta1 = 0.06;     % birth-rate parameter, population 1
beta2 = 0.08;     % birth-rate parameter, population 2
Nmax = 1e4;       % max size of the metapopulation
tmax = 2e4;       % the number of mappings

tic
disp(['beta = ' num2str(beta1)])
[n1,exT1] = RunDeterministicModel(Nmax,beta1,delta,tmax);
toc

tic
disp(['beta = ' num2str(beta2)])
[n2,exT2] = RunDeterministicModel(Nmax,beta2,delta,tmax);
toc

%% make figures

figure(1),clf

subplot(211)

hold on
plot(n1(:,1),'Color',0.50*ones(1,3),'LineWidth',2)
plot(n1(:,2),'Color',0.00*ones(1,3),'LineWidth',2)
hold off

T1 = text(9000,200,'a','FontSize',35,'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex','FontSize',20,'LineWidth',1)
ylabel('populations','Interpreter','Latex')
grid on, box on, axis tight

LG = legend('$(y,z)=(10,10)$','$(y,z)=(0,10)$');
LG.Interpreter = 'Latex';
LG.FontSize = 15;

subplot(212)

hold on
plot(n2(:,1),'Color',0.50*ones(1,3),'LineWidth',2)
plot(n2(:,2),'Color',0.00*ones(1,3),'LineWidth',2)
hold off

text(9000,400,'b','FontSize',35,'Interpreter','Latex')
set(gca,'TickLabelInterpreter','Latex','FontSize',20,'LineWidth',1)
xlabel('mappings','Interpreter','Latex')
ylabel('populations','Interpreter','Latex')
grid on, box on, axis tight

LG = legend('$(y,z)=(10,10)$','$(y,z)=(0,10)$');
LG.Interpreter = 'Latex';
LG.FontSize = 15;


set(gcf,'Pos',[200,80,650,500])

%% END
