width = 1; % plot linewidth
x = linspace(0,50,51);

RK3full(51,:) = RK3full(1,:);
LFfull(51,:) = LFfull(1,:);
FTBSfull(51,:) = FTBSfull(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 1A: EXPERIMENT #1 (solution propogation: FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(3,1,1);
for i = 1:3
    plot(x,FTBSfull(:,i),'k','LineWidth',width);
    hold on
end
for i = 199:201
    plot(x,FTBSfull(:,i),':k','LineWidth',width);
end
plot(x,FTBSfull(:,1),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(3,1,2);
for i = 1:3
    plot(x,LFfull(:,i),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 199:201
    plot(x,LFfull(:,i),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
plot(x,LFfull(:,1),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(3,1,3);
for i = 1:3
    plot(x,RK3full(:,i),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 199:201
    plot(x,RK3full(:,i),':','Color',[0 0.447 0.741],'LineWidth',width);
    plot(x,exactsolution(:,i),':k','LineWidth',width);
    hold on
end
plot(x,RK3full(:,1),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 1B: EXPERIMENT #1 (log(A) and log(error): FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
loglog(exp1.summary.FTBS(:,1),exp1.summary.FTBS(:,3),'LineWidth',width);
hold on
loglog(exp1.summary.LF(:,1),exp1.summary.LF(:,3),'LineWidth',width);
loglog(exp1.summary.RK3(:,1),exp1.summary.RK3(:,3),'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3'},'Location','Southeast');


subplot(1,2,2)
timestep = linspace(0,1000,1001);
semilogy(timestep,exp1.summary.FTBS(:,2),'LineWidth',width);
hold on
semilogy(timestep,exp1.summary.LF(:,2),'LineWidth',width);
semilogy(timestep,exp1.summary.RK3(:,2),'LineWidth',width);
xlabel('t_n');
ylabel('A_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3'},'Location','Southwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2A: EXPERIMENT #2 (solution propogation: FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1);
for i = 1:3
    plot(x,exp2.dt2.full.FTBS(:,i),'k','LineWidth',width);
    hold on
end
for i = 15:17
    plot(x,exp2.dt2.full.FTBS(:,i),':k','LineWidth',width);
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,2);
for i = 1:3
    plot(x,exp2.dt2.full.LF(:,i),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 9:12
    plot(x,exp2.dt2.full.LF(:,i),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,3);
for i = 1:3
    plot(x,exp2.dt2.full.RK3(:,i),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 154:156
    plot(x,exp2.dt2.full.RK3(:,i),':','Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,4);
for i = 1:3
    plot(x,exp2.dt2.full.AB3(:,i),'Color',[0.929 0.694 0.125],'LineWidth',width);
    hold on
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2B: EXPERIMENT #2 (log(A) and log(error): FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
loglog(FTBSsummary(:,1),FTBSsummary(:,3),'LineWidth',width);
hold on
loglog(LFsummary(:,1),LFsummary(:,3),'LineWidth',width);
loglog(RK3summary(:,1),RK3summary(:,3),'LineWidth',width);
loglog(AB3summary(:,1),AB3summary(:,3),'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northwest');

subplot(1,2,2)
timestep = linspace(0,800,801);
semilogy(timestep,FTBSsummary(:,2),'LineWidth',width);
hold on
semilogy(timestep,LFsummary(:,2),'LineWidth',width);
semilogy(timestep,RK3summary(:,2),'LineWidth',width);
semilogy(timestep,AB3summary(:,2),'LineWidth',width);
xlabel('t_n');
ylabel('A_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 3A: EXPERIMENT #3 (solution propogation: FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1);
for i = 1
    plot(x,exp3.full.FTBS(:,i),'k','LineWidth',width);
    hold on
end
for i = 200
    plot(x,FTBSfull(:,i),':k','LineWidth',width);
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,2);
for i = 1
    plot(x,LFfull(:,i),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 200
    plot(x,LFfull(:,i),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,3);
for i = 1
    plot(x,RK3full(:,i),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 50
    plot(x,RK3full(:,i),':','Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,4);
for i = 1
    plot(x,AB3full(:,i),'Color',[0.929 0.694 0.125],'LineWidth',width);
    hold on
end
for i = 3
    plot(x,AB3full(:,i),':','Color',[0.929 0.694 0.125],'LineWidth',width);
    hold on
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 3B: EXPERIMENT #3 (log(A) and log(error): FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
loglog(exp3.summary.FTBS(:,1),exp3.summary.FTBS(:,3),'LineWidth',width);
hold on
loglog(exp3.summary.LF(:,1),exp3.summary.LF(:,3),'LineWidth',width);
loglog(exp3.summary.RK3(:,1),exp3.summary.RK3(:,3),'LineWidth',width);
loglog(exp3.summary.AB3(:,1),exp3.summary.AB3(:,3),'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northwest');

subplot(1,2,2)
timestep = linspace(0,1000,1001);
semilogy(timestep,exp3.summary.FTBS(:,2),'LineWidth',width);
hold on
semilogy(timestep,exp3.summary.LF(:,2),'LineWidth',width);
semilogy(timestep,exp3.summary.RK3(:,2),'LineWidth',width);
semilogy(timestep,exp3.summary.AB3(:,2),'LineWidth',width);
xlabel('t_n');
ylabel('A_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 4A: EXPERIMENT #4 (solution propogation: FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,1);
for i = 1:3
    plot(x,FTBSfull(:,i),'k','LineWidth',width);
    hold on
end
for i = 15:17
    plot(x,FTBSfull(:,i),':k','LineWidth',width);
end
plot(x, FTBSfull(:,1),'k','LineWidth',width+0.5);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,2);
for i = 1:3
    plot(x,exp4.full.LF(:,i),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 399:401
    plot(x,exp4.full.LF(:,i),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
plot(x, exp4.full.LF(:,1),'k','LineWidth',width+0.5);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,3);
for i = 1:3
    plot(x,RK3full(:,i),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 399:401
    plot(x,RK3full(:,i),':','Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
plot(x, RK3full(:,1),'k','LineWidth',width+0.5);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,4);
for i = 1:3
    plot(x,AB3full(:,i),'Color',[0.929 0.694 0.125],'LineWidth',width);
    hold on
end
for i = 399:401
    plot(x,AB3full(:,i),':','Color',[0.929 0.694 0.125],'LineWidth',width);
    hold on
end
plot(x, AB3full(:,1),'k','LineWidth',width+0.5);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 4B: EXPERIMENT #4 (log(A) and log(error): FTBS, Leap-Frog, RK3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
loglog(FTBSsummary(:,1),FTBSsummary(:,3),'LineWidth',width);
hold on
loglog(LFsummary(:,1),LFsummary(:,3),'LineWidth',width);
loglog(RK3summary(:,1),RK3summary(:,3),'LineWidth',width);
loglog(AB3summary(:,1),AB3summary(:,3),'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northwest');

subplot(1,2,2)
timestep = linspace(0,2000,2001);
semilogy(timestep,FTBSsummary(:,2),'LineWidth',width);
hold on
semilogy(timestep,LFsummary(:,2),'LineWidth',width);
semilogy(timestep,RK3summary(:,2),'LineWidth',width);
semilogy(timestep,AB3summary(:,2),'LineWidth',width);
xlabel('t_n');
ylabel('A_u(t)');
ax = gca;
ax.FontSize = 13;
legend({'FTBS','Leap-Frog','RK3','AB3'},'Location','Northeast');