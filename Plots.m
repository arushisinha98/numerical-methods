width = 0.5; % plot linewidth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 1: EULER-FORWARD z(t) VS t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot solution for dt = 10.0s
x1 = EF10s(:,1);
y1 = EF10s(:,2);

% plot solution for dt = 1.0s
x2 = EF1s(:,1);
y2 = EF1s(:,2);

% plot analytical solution
x3 = linspace(0,2000,2001);
y3 = 100*cos(sqrt((9.81/290)*(0.005))*x3);

figure
plot(x1,y1,x2,y2,x3,y3,'k','LineWidth',width);
legend({'\deltat = 10.0 s','\deltat =   1.0 s','z_p(t) = z_0cos(Nt)'},'Location','northwest');
xlabel('Time / s');
ylabel('z_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 1: z_p(t) for exact solution and numerical solutions using EF Scheme');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 2: EULER-FORWARD w(t) VS z(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [];
y1 = [];
x2 = [];
y2 = [];

% plot solution for dt = 10.0s
x1 = EF10s(:,2);
y1 = EF10s(:,3);

% plot solution for dt = 1.0s
x2 = EF1s(:,2);
y2 = EF1s(:,3);

figure
plot(x1,y1,x2,y2,'LineWidth',width);
legend({'dt = 10.0 s','dt = 1.0 s'},'Location','northwest');
xlabel('z_p(t)');
ylabel('w_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 2: w_p(t) for numerical solutions using EF Scheme');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 3: EULER-FORWARD error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [];
y1 = [];
x2 = [];
y2 = [];

% plot solution for dt = 10.0s
x1 = EF10s(:,1);
y1 = EF10s(:,4);

% plot solution for dt = 1.0s
x2 = EF1s(:,1);
y2 = EF1s(:,4);

figure
loglog(x1,y1,x2,y2,'LineWidth',width);
legend({'dt = 10.0 s','dt = 1.0 s'},'Location','northwest');
xlabel('Time / s');
ylabel('\epsilon_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 3: \epsilon_p(t) for numerical solutions using EF Scheme');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 4: LEAP-FROG z(t) VS t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [];
y1 = [];
x2 = [];
y2 = [];

% plot solution for dt = 10.0s
x1 = LF10s(:,1);
y1 = LF10s(:,2);

% plot solution for dt = 1.0s
x2 = LF1s(:,1);
y2 = LF1s(:,2);

% plot analytical solution
x3 = linspace(0,2000,2001);
y3 = 100*cos(sqrt((9.81/290)*(0.005))*x3);

figure
plot(x1,y1,x2,y2,x3,y3,'k','LineWidth',width);
legend({'dt = 10.0 s','dt =   1.0 s','z_p(t) = z_0cos(Nt)'},'Location','northwest');
xlabel('Time / s');
ylabel('z_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 4: z_p(t) for exact solution and numerical solutions using LF Scheme');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 5: LEAP-FROG w(t) VS z(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [];
y1 = [];
x2 = [];
y2 = [];

% plot solution for dt = 10.0s
x1 = LF10s(:,2);
y1 = LF10s(:,3);

% plot solution for dt = 1.0s
x2 = LF1s(:,2);
y2 = LF1s(:,3);

figure
plot(x1,y1,x2,y2,'LineWidth',width);
legend({'dt = 10.0 s','dt =   1.0 s'},'Location','northwest');
xlabel('z_p(t)');
ylabel('w_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 5: w_p(t) for numerical solutions using LF Scheme');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURE 6: LEAP-FROG error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [];
y1 = [];
x2 = [];
y2 = [];

% plot solution for dt = 10.0s
x1 = LF10s(:,1);
y1 = LF10s(:,4);

% plot solution for dt = 1.0s
x2 = LF1s(:,1);
y2 = LF1s(:,4);

figure
plot(x1,y1,x2,y2,'LineWidth',width);
legend({'dt = 10.0 s','dt = 1.0 s'},'Location','northwest');
xlabel('Time / s');
ylabel('\epsilon_p(t)');
ax = gca;
ax.FontSize = 13;
title('Figure 6: \epsilon_p(t) for numerical solutions using LF Scheme');

