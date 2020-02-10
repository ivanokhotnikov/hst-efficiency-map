function printFit
load('Data.mat','cat','fit');
y=@(c,x) c(1).*x+c(2);
y1=@(c,x) c(1).*exp(x.*(-c(2)))+c(3);
Dp=[cat.brevini.Dp;cat.rexroth.Dp;cat.eaton.Dp;cat.parker.Dp;...
    cat.danfoss.Dp;cat.kawasaki.Dp;cat.hydac.Dp;cat.linde.Dp;cat.poclain.Dp];
Dm=[cat.rexroth.Dm;cat.eaton.Dm;cat.danfoss.Dm;cat.hydac.Dm;...
    cat.parker.Dm;cat.kawasaki.Dm;cat.linde.Dm;cat.poclain.Dm];
xp=linspace(min(Dp),max(Dp),length(Dp));
xm=linspace(min(Dm),max(Dm),length(Dm));
%% PUMP MASS
figure('Name','Pump mass objective','NumberTitle','off',...
    'Position',[50 50 650 450])
scatter(cat.eaton.Dp,cat.eaton.mp,20,'Linewidth',1);
set(groot,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
scatter(cat.rexroth.Dp,cat.rexroth.mp,20,'+','Linewidth',1);
scatter(cat.danfoss.Dp,cat.danfoss.mp,20,'*','Linewidth',1);
scatter(cat.parker.Dp,cat.parker.mp,20,'s','Linewidth',1);
scatter(cat.kawasaki.Dp,cat.kawasaki.mp,20,'x','Linewidth',1);
scatter(cat.hydac.Dp,cat.hydac.mp,20,'p','Linewidth',1);
scatter(cat.linde.Dp,cat.linde.mp,20,'h','Linewidth',1);
scatter(cat.poclain.Dp,cat.poclain.mp,20,'^','Linewidth',1);
scatter(cat.brevini.Dp,cat.brevini.mp,20,'d','Linewidth',1);
plot(xp, y(fit.C(:,1),xp),'-','Linewidth',0.75,'Color',[0 0.4470 0.7410]);
plot(xp, y(fit.C(:,1),xp)+fit.rmse(1),'--','Linewidth',0.75,'Color',...
    [0.8500 0.3250 0.0980]);
plot(xp, y(fit.C(:,1),xp)-fit.rmse(1),'--','Linewidth',0.75,'Color',...
    [0.9290 0.6940 0.1250]);
legend({'Eaton','Rexroth','Danfoss','Parker','Kawasaki','Hydac',...
    'Linde','Poclain','Brevini','Fit','Fit + \(\sigma_{mp}\)',...
    'Fit -- \(\sigma_{mp}\)'},'Location','northwest',...
    'Interpreter','latex');
ylim([0 600]);
yticks(0:100:600);
xlabel('Pump displacement, cc/rev');
ylabel('Pump mass, kg');
hold off;
grid on;
%% MOTOR MASS
figure('Name','Motor mass objective','NumberTitle','off',...
    'Position',[50 50 650 450])
scatter(cat.eaton.Dm,cat.eaton.mm,20,'Linewidth',1);
set(groot,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
scatter(cat.rexroth.Dm,cat.rexroth.mm,20,'+','Linewidth',1);
scatter(cat.danfoss.Dm,cat.danfoss.mm,20,'*','Linewidth',1);
scatter(cat.parker.Dm,cat.parker.mm,20,'s','Linewidth',1);
scatter(cat.kawasaki.Dm,cat.kawasaki.mm,20,'x','Linewidth',1);
scatter(cat.hydac.Dm,cat.hydac.mm,20,'p','Linewidth',1);
scatter(cat.linde.Dm,cat.linde.mm,20,'h','Linewidth',1);
scatter(cat.poclain.Dm,cat.poclain.mm,20,'^','LineWidth',1);
plot(xm, y(fit.C(:,2),xm),'-','Linewidth',0.75,'Color',[0 0.4470 0.7410]);
plot(xm, y(fit.C(:,2),xm)+fit.rmse(2),'--','Linewidth',0.75,'Color',...
    [0.8500 0.3250 0.0980]);
plot(xm, y(fit.C(:,2),xm)-fit.rmse(2),'--','Linewidth',0.75,'Color',...
    [0.9290 0.6940 0.1250]);
legend({'Eaton','Rexroth','Danfoss','Parker','Kawasaki','Hydac',...
    'Linde','Poclain','Fit','Fit + \(\sigma_{mm}\)','Fit -- \(\sigma_{mm}\)'},...
    'Location','northwest','Interpreter','latex');
hold off;
xlabel('Motor displacement, cc/rev');
ylabel('Motor mass, kg');
ylim([0 500]);
yticks(0:100:500);
xlim([0 800]);
xticks(0:100:800);
grid on;
%% MASS MODELS ERRORS
figure('Name','Mass models errors','NumberTitle','off',...
    'Position',[50 50 650 450])
set(groot,'defaultTextInterpreter','latex');
scatter(Dp,fit.resp(:,1),20,'Linewidth',1);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
xlabel('Displacement, cc/rev');
ylabel('Error, kg');
hold on
scatter(Dm,fit.resm(:,1),20,'+','Linewidth',1);
hold off;
legend({'Pump','Motor'},'Location','northeast','Interpreter','latex');
grid on;
%% HISTOGRAM OF THE ERRORS
figure('Name','Histogram of the mass models errors','NumberTitle','off',...
    'Position',[50 50 650 450]);
set(groot,'defaultTextInterpreter','latex');
histogram(fit.resp(:,1),'Normalization','count','BinWidth',5);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
histogram(fit.resm(:,1),'Normalization','count','BinWidth',5);
hold off
ylabel('Frequency');
xlabel('Error, kg');
grid on;
legend({'Pump','Motor'},'Location','northeast','Interpreter','latex');
box off
%% PUMP SPEED
figure('Name','Pump speed constraint','NumberTitle','off',...
    'Position',[50 50 650 450])
scatter(cat.eaton.Dp,cat.eaton.np,20,'Linewidth',1);
set(groot,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
scatter(cat.rexroth.Dp,cat.rexroth.np,20,'+','Linewidth',1);
scatter(cat.danfoss.Dp,cat.danfoss.np,20,'*','Linewidth',1);
scatter(cat.parker.Dp,cat.parker.np,20,'s','Linewidth',1);
scatter(cat.kawasaki.Dp,cat.kawasaki.np,20,'x','Linewidth',1);
scatter(cat.hydac.Dp,cat.hydac.np,20,'p','Linewidth',1);
scatter(cat.linde.Dp,cat.linde.np,20,'h','Linewidth',1);
scatter(cat.poclain.Dp,cat.poclain.np,20,'^','Linewidth',1);
scatter(cat.brevini.Dp,cat.brevini.np,20,'d','Linewidth',1);
plot(xm, y1(fit.Cn(:,1),xm),'-','Linewidth',1,'Color',[0 0.4470 0.7410]);
plot(xm, y1(fit.Cn(:,1),xm)+fit.rmse(3),'--','Linewidth',1,'Color',...
    [0.8500 0.3250 0.0980]);
plot(xm, y1(fit.Cn(:,1),xm)-fit.rmse(3),'--','Linewidth',1,'Color',...
    [0.9290 0.6940 0.1250]);
legend({'Eaton','Rexroth','Danfoss','Parker','Kawasaki','Hydac','Linde',...
    'Poclain','Brevini','Fit','Fit + \(\sigma_{np}\)',...
    'Fit -- \(\sigma_{np}\)'},'Location','northeast','Interpreter','latex');
hold off;
xlabel('Pump displacement, cc/rev');
ylabel('Pump speed, rpm');
grid on;
%% MOTOR SPEED
figure('Name','Motor speed constraint','NumberTitle','off',...
    'Position',[50 50 650 450])
set(groot,'defaultTextInterpreter','latex');
scatter(cat.eaton.Dm,cat.eaton.nm,20,'Linewidth',1);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
scatter(cat.rexroth.Dm,cat.rexroth.nm,20,'+','Linewidth',1);
scatter(cat.danfoss.Dm,cat.danfoss.nm,20,'*','Linewidth',1);
scatter(cat.parker.Dm,cat.parker.nm,20,'s','Linewidth',1);
scatter(cat.kawasaki.Dm,cat.kawasaki.nm,20,'x','Linewidth',1);
scatter(cat.hydac.Dm,cat.hydac.nm,20,'p','Linewidth',1);
scatter(cat.linde.Dm,cat.linde.nm,20,'h','Linewidth',1);
scatter(cat.poclain.Dm,cat.poclain.nm,20,'^','Linewidth',1);
plot(xm, y1(fit.Cn(:,2),xm),'-','Linewidth',1,'Color',[0 0.4470 0.7410]);
plot(xm, y1(fit.Cn(:,2),xm)+fit.rmse(4),'--','Linewidth',1,'Color',...
    [0.8500 0.3250 0.0980]);
plot(xm, y1(fit.Cn(:,2),xm)-fit.rmse(4),'--','Linewidth',1,'Color',...
    [0.9290 0.6940 0.1250]);
legend({'Eaton','Rexroth','Danfoss','Parker','Kawasaki','Hydac','Linde',...
    'Poclain','Fit','Fit + \(\sigma_{nm}\)','Fit -- \(\sigma_{nm}\)'},...
    'Location','northeast','Interpreter','latex');
hold off;
xlabel('Motor displacement, cc/rev');
ylabel('Motor speed, rpm');
grid on;
%% SPEED MODELS ERRORS
figure('Name','Speed models errors','NumberTitle','off',...
    'Position',[50 50 650 450])
set(groot,'defaultTextInterpreter','latex');
scatter(Dp,fit.resp(:,2),20,'Linewidth',1);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
xlabel('Displacement, cc/rev');
ylabel('Error, rpm');
hold on;
scatter(Dm,fit.resm(:,2),20,'+','Linewidth',1);
hold off;
grid on;
legend({'Pump','Motor'},'Location','northeast','Interpreter','latex');
%% HISTOGRAM OF THE ERRORS
figure('Name','Histogram of the speed models errors','NumberTitle','off',...
    'Position',[50 50 650 450]);
set(groot,'defaultTextInterpreter','latex');
histogram(fit.resp(:,2),'Normalization','count','BinWidth',50);
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',12);
hold on;
histogram(fit.resm(:,2),'Normalization','count','BinWidth',50);
hold off
ylabel('Frequency');
xlabel('Error, rpm');
ylim([0 14])
yticks(0:2:14);
grid on;
legend({'Pump','Motor'},'Location','northeast','Interpreter','latex');
box off;
end
