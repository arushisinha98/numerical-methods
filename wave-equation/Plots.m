width = 1; % plot linewidth
x = linspace(0,50,51);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 1A: EXPERIMENT #1 (exact solution propogation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(x,exp1.exact(1,:),x,exp1.exact(2,:),x,exp1.exact(3,:),'LineWidth',width);
legend({'time_1 = 10.0 s','time_2 = 20.0 s','time_3 = 30.0 s'},'Location','northeast');
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;
title('Figure 1A: Exact & Numerical Solutions for \Deltat = 2.00 s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 1B: EXPERIMENT #1 (log(A) and log(error))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = exp1.error;
y2 = exp1.A;
time = linspace(0,200,101);
timestep = linspace(0,100,101);

figure
title('Figure 1B: Error and amplitude for numerical solutions');
subplot(1,2,1)
loglog(time,y1,'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca;
ax.FontSize = 13;

subplot(1,2,2)
semilogy(timestep,y2,'LineWidth',width);
xlabel('t_n');
ylabel('A_u(t)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2A: EXPERIMENT #2 (exact solution propogation)
% dt = 0.25, 0.5, 0.75, 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,2,4);
for i = 1:5
    plot(x,exp2.dt1.approx(i,:),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 37:41
    plot(x,exp2.dt1.approx(i,:),':','Color',[0 0.447 0.741],'LineWidth',width);
end
plot(x,exp2.dt1.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,2);
for i = 1:5
    plot(x,exp2.dt2.approx(i,:),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 77:81
    plot(x,exp2.dt2.approx(i,:),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
plot(x,exp2.dt2.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,1);
for i = 1:5
    plot(x,exp2.dt3.approx(i,:),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 157:161
    plot(x,exp2.dt3.approx(i,:),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
plot(x,exp2.dt2.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,2,3);
for i = 1:5
    plot(x,exp2.dt4.approx(i,:),'Color',[0 0.447 0.741],'LineWidth',width);
    hold on
end
for i = 50:54
    plot(x,exp2.dt4.approx(i,:),':','Color',[0 0.447 0.741],'LineWidth',width);
end
plot(x,exp2.dt4.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2B: EXPERIMENT #2 (log(A) and log(error))
% dt = 0.25, 0.5, 0.75, 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = exp2.dt3.error;
y2 = exp2.dt2.error;
y3 = exp2.dt4.error;
y4 = exp2.dt1.error;

figure
subplot(1,2,1)
loglog(linspace(0,200,(200/0.25)+1),y1,'LineWidth',width);
hold on
loglog(linspace(0,200,(200/0.50)+1),y2,'LineWidth',width);
loglog(linspace(0,200,(200/0.75)+1),y3,'LineWidth',width);
loglog(linspace(0,200,201),y4,'LineWidth',width);
legend({'\Deltat = 0.25 s','\Deltat = 0.50 s','\Deltat = 0.75 s','\Deltat = 1.00 s'},'Location','northwest');
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca
ax.FontSize = 13;

y1 = exp2.dt3.A;
y2 = exp2.dt2.A;
y3 = exp2.dt4.A;
y4 = exp2.dt1.A;

subplot(1,2,2)
semilogy(linspace(0,800,801),y1,'LineWidth',width);
hold on
semilogy(linspace(0,400,401),y2,'LineWidth',width);
semilogy(linspace(0,266,267),y3,'LineWidth',width);
semilogy(linspace(0,200,201),y4,'LineWidth',width);
legend({'\Deltat = 0.25 s','\Deltat = 0.50 s','\Deltat = 0.75 s','\Deltat = 1.00 s'},'Location','northeast');
xlabel('Time-step (t_n)');
ylabel('A_u(t)');
ax = gca
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 3A: EXPERIMENT #3 (exact solution propogation)
% dt = 4, 5, 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for i = 1:5
    plot(x,exp3.dt1.approx(i,:),'Color',[0.85 0.325 0.098],'LineWidth',width);
    hold on
end
for i = 11:15
    plot(x,exp3.dt1.exact(i,:),'k','LineWidth',width);
    plot(x,exp3.dt1.approx(i,:),':','Color',[0.85 0.325 0.098],'LineWidth',width);
end
for i = 1:5
    plot(x,exp3.dt1.exact(i,:),'k','LineWidth',width);
end
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;
title('Figure 3A: Exact & Numerical Solutions for \Deltat = 4.00 s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 3B: EXPERIMENT #3 (log(A) and log(error))
% dt = 4, 5, 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = exp3.dt1.error;
y2 = exp3.dt3.error;
y3 = exp3.dt2.error;

figure
subplot(1,2,1)
loglog(linspace(0,200,(200/4)+1),y1,'LineWidth',width);
hold on
loglog(linspace(0,200,(200/5)+1),y2,'LineWidth',width);
loglog(linspace(0,200,(200/8)+1),y3,'LineWidth',width);
legend({'\Deltat = 4.00 s','\Deltat = 5.00 s','\Deltat = 8.00 s'},'Location','northwest');
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca
ax.FontSize = 13;

y1 = exp3.dt1.A;
y2 = exp3.dt3.A;
y3 = exp3.dt2.A;

subplot(1,2,2)
semilogy(linspace(0,50,51),y1,'LineWidth',width);
hold on
semilogy(linspace(0,40,41),y2,'LineWidth',width);
semilogy(linspace(0,25,26),y3,'LineWidth',width);
legend({'\Deltat = 4.00 s','\Deltat = 5.00 s','\Deltat = 8.00 s'},'Location','northeast');
xlabel('Time-step (t_n)');
ylabel('A_u(t)');
ax = gca
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 4A: EXPERIMENT #4 (exact solution propogation)
% dt = 1, 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(exp4.dt1.exact(:,1))
    for j = 1:length(exp4.dt1.exact(1,:))
        if isnan(exp4.dt1.exact(i,j))
            exp4.dt1.exact(i,j) = 0;
        end
    end
end

for i = 1:length(exp4.dt1.approx(:,1))
    for j = 1:length(exp4.dt1.approx(1,:))
        if isnan(exp4.dt1.approx(i,j))
            exp4.dt1.approx(i,j) = 0;
        end
    end
end

for i = 1:length(exp4.dt2.exact(:,1))
    for j = 1:length(exp4.dt2.exact(1,:))
        if isnan(exp4.dt2.exact(i,j))
            exp4.dt2.exact(i,j) = 0;
        end
    end
end

for i = 1:length(exp4.dt2.approx(:,1))
    for j = 1:length(exp4.dt2.approx(1,:))
        if isnan(exp4.dt2.approx(i,j))
            exp4.dt2.approx(i,j) = 0;
        end
    end
end

figure
subplot(2,1,1);
for i = 1:3
    hold on
    plot(x,exp4.dt1.approx(i,:),'Color',[0.85 0.325 0.098],'LineWidth',width);
end
plot(x,exp4.dt1.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

subplot(2,1,2);
for i = 1:3
    hold on
    plot(x,exp4.dt2.approx(i,:),'Color',[0 0.447 0.741],'LineWidth',width);
end
plot(x,exp4.dt2.exact(1,:),'k','LineWidth',width);

xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 4B: EXPERIMENT #4 (log(A) and log(error))
% dt = 1, 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = exp4.dt1.error;
y2 = exp4.dt2.error;

figure
subplot(1,2,1)
loglog(linspace(0,200,(200/2)+1),y2,'LineWidth',width);
hold on
loglog(linspace(0,200,(200/1)+1),y1,'LineWidth',width);
legend({'\Deltat = 2.00 s','\Deltat = 1.00 s'},'Location','northeast');
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca
ax.FontSize = 13;

y1 = exp4.dt1.A;
y2 = exp4.dt2.A;

subplot(1,2,2)
semilogy(linspace(0,100,101),y2,'LineWidth',width);
hold on
semilogy(linspace(0,200,201),y1,'LineWidth',width);
legend({'\Deltat = 2.00 s','\Deltat = 1.00 s'},'Location','northeast');
xlabel('Time-step (t_n)');
ylabel('A_u(t)');
ax = gca
ax.FontSize = 13;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 5A: EXPERIMENT #5 (exact solution propogation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(x,exp5.exact(1,:),'k','LineWidth',width);
hold on
for i = 1:5
    plot(x,exp5.approx(i,:),'Color',[0 0.447 0.741],'LineWidth',width);
end
for i = 31:35
    plot(x,exp5.approx(i,:),':','Color',[0 0.447 0.741],'LineWidth',width);
end
plot(x,exp5.exact(1,:),'k','LineWidth',width);
xlabel('x');
ylabel('u(x,t_n)');
ax = gca;
ax.FontSize = 13;
title('Figure 5A: Numerical Solutions for Square Wave with \Deltat = 1.00 s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 5B: EXPERIMENT #5 (log(A) and log(error))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1 = exp5.error;
figure
subplot(1,2,1)
loglog(linspace(0,200,201),y1,'LineWidth',width);
xlabel('Time / s');
ylabel('\epsilon_u(t)');
ax = gca
ax.FontSize = 13;

y1 = exp5.A;

subplot(1,2,2)
semilogy(linspace(0,200,201),y1,'LineWidth',width);
xlabel('Time-step (t_n)');
ylabel('A_u(t)');
ax = gca
ax.FontSize = 13;
