%% Dependencies:
%%             Octave 4.0.0
%%
%%             Package Name  | Version
%%             --------------+---------
%%                 geometry  |   2.1.1
%%                 symbolic  |   2.7.1

clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load  packages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load symbolic;
pkg load geometry;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constant variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global saving_foler visualize
visualize = false;
saving_foler = '/tmp/output'; %% don't change this
mkdir("/tmp", "output");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_new = dynamics(x, u, dt_s)
  A = [1 dt_s;
       0 1];
  B = [0.5 * dt_s^2;
       dt_s];
  x_new = A * x + B * u;
end

function u = PD_controller(x, K)
  u = - K * [x(1) x(2)]';
end

function traj = sim_tvc(x0, controller, duration, dt_s)
  %% time varying controller
  time = 0;
  x = x0;
  traj = [];
  time_step = 0;
  while time < duration
    u = controller(x, time_step);
    traj = [traj; time x' u'];
    x = dynamics(x, u, dt_s);

    time = time + dt_s;
    time_step += 1;
  end
end

function traj = sim(x0, controller, duration, dt_s)
  time = 0;
  x = x0;
  traj = [];
  while time < duration
    u = controller(x);
    traj = [traj; time x' u'];
    x = dynamics(x, u, dt_s);

    time = time + dt_s;
  end
end

function show_traj(traj)
  plot(traj(:,1), traj(:,2:end), "LineWidth", 3);
  legend('position', 'speed', 'u')
  xlim([0 10]);
  xlabel("Time [s]")
  ylabel("States");
  title("Time response");
end

function [LL, T] = get_system_eigenvalues()
  %% general second-order discrete-time state system
  syms a1 a2 a3 a4 b1 b2 k1 k2 dt
  A = [a1 a2;
       a3 a4]
  B = [b1;
       b2]
  K = [k1 k2]
  T = A - B * K;
  LL =  eig(A - B*K);

  %% Specialize
  LL =  subs(LL, a1, 1);
  LL =  subs(LL, a2, dt);
  LL =  subs(LL, a3, 0);
  LL =  subs(LL, a4, 1);
  LL =  subs(LL, b1, dt^2/2);
  LL =  subs(LL, b2, dt);
  simplify(LL)

  %% Specialize
  T =  subs(T, a1, 1);
  T =  subs(T, a2, dt);
  T =  subs(T, a3, 0);
  T =  subs(T, a4, 1);
  T =  subs(T, b1, dt^2/2);
  T =  subs(T, b2, dt);
  simplify(T)
end

function plot_poles_gains(lambda1, lambda2, dt_s, k1, k2)
  %% Function to plot system pole for a given dt_x, k1, and k2
  %% @params : lambda1 and lambda2 are functions of dt_s, k1, k2

  global K1 K2 STABLE_REGION
  clf
  p1 = lambda1(dt_s, k1, k2);
  p2 = lambda2(dt_s, k1, k2);

  subplot(2,2,1);
  if size(STABLE_REGION, 1) > 0
    h = pcolor(K1, K2, STABLE_REGION); set(h, 'EdgeColor', 'none'); colormap('gray');
    hold on;
  end
  plot(k1, k2, "*", "LineWidth", 3);
  xlabel("k1"); ylabel("k2"); xlim([0 120]); ylim([0 25]); title("K selection");

  subplot(2,2,2);
  drawCircle(0,0,1.0);
  hold on
  plot(real(p1), imag(p1), "*", "LineWidth", 3,
       real(p2), imag(p2), "*", "LineWidth", 3);
  grid on; axis equal
  xlim([-2 2]); ylim([-2 2]); title(sprintf("Poles at k1 = %.1f  k2= %.1f", k1, k2))

  %% Time response
  x0 = [1 0]';
  controller = @(x) PD_controller(x, [k1 k2]);
  traj = sim(x0, controller, 10.0, dt_s);

  subplot(2,2,3:4); show_traj(traj);
end

function plot_multiple_poles_gains(lambda1, lambda2, dt_s, k1_list, k2_list, visualize = true, prefix='')
  %% lambda1 and lambda2 are functions of k1, k2, where K = [k1 k2]
  global saving_foler

  idx = 0;
  if !visualize
    h = gcf();
    set(h, 'visible','off')
  end

  for k1 = k1_list
    for k2 = k2_list
      plot_poles_gains(lambda1, lambda2, dt_s, k1, k2);
      if visualize
        refresh() %% drawnow()
        pause(0.2)
      else
        filename = sprintf("%s/%s%05d.png", saving_foler, prefix, idx);
        print(filename)
      end

      idx = idx + 1;
    end
  end
end

function plot_stable_region(lambda1, lambda2, dt_s, visualize)
  nk1 = 0:0.5:120;
  nk2 = 0:0.1:25;
  [K1, K2] = meshgrid(nk1,nk2);
  L1 = lambda1(dt_s, K1, K2);
  L2 = lambda2(dt_s, K1, K2);

  %% unstable region mask
  abs_mask = abs(L1) > 1 | abs(L2) > 1.0;
  %% underdamping region mask
  arg_mask = abs(arg(L1)) > 0 | abs(arg(L2)) > 0;

  STABLE_REGION = zeros(size(K1));
  STABLE_REGION(!abs_mask) = 0.5;
  STABLE_REGION(!arg_mask) = 1.0;

  if !visualize
    h = gcf();
    set(h, 'visible','off')
  end
  clf
  h = pcolor(K1, K2, STABLE_REGION); set(h, 'EdgeColor', 'none'); colormap('gray');
  xlabel("k1"); ylabel("k2"); title(sprintf("Stable and underdamped region when dt = %0.3f", dt_s));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LL, T] = get_system_eigenvalues();

%% Convert to function for plotting
lambda1 = function_handle(LL(1))
lambda2 = function_handle(LL(2))

%% variables for plotting
dt_s = 0.2;
nk1 = 0:0.5:120;
nk2 = 0:0.1:25;

global K1
global K2
[K1, K2] = meshgrid(nk1,nk2);
L1 = lambda1(dt_s, K1, K2);
L2 = lambda2(dt_s, K1, K2);

%% unstable region mask
abs_mask = abs(L1) > 1 | abs(L2) > 1.0;
%% underdamping region mask
arg_mask = abs(arg(L1)) > 0 | abs(arg(L2)) > 0;

global STABLE_REGION;
STABLE_REGION = zeros(size(K1));
STABLE_REGION(!abs_mask) = 0.5;
STABLE_REGION(!arg_mask) = 1.0;

################################################################################
# critical damping
################################################################################
syms dt k1 k2 
f1 = (LL(1) - LL(2) ) ^ 2
f2 = (LL(1) + LL(2) ) / 2;

## k1 in term of k2 when critical damping occurs
## the solution is not unique, but only the first one give stable system
k1_k2 = solve(f1 == 0, k1)

## the eigen values in term of k2 when critical damping occurs
lambda_k2 = subs(f2, k1, k1_k2)

dt_s = 0.2;
fun_lambda = @(k2) function_handle(lambda_k2)(dt_s, k2)
fun_k1 = @(k2) function_handle(k1_k2)(dt_s, k2)

k2 = 0:0.1:8;
close all
figure(1)
subplot(2,1,1);
plot(k2, fun_k1(k2)(1,:), "-r", "LineWidth", 3); hold on
plot(k2, fun_k1(k2)(2,:), "-b", "LineWidth", 3); grid on
title("Critical damping for dt = 0.2")
xlabel("k_2")
ylabel("k_1")
subplot(2,1,2);
plot(k2, fun_lambda(k2)(1,:), "-r", "LineWidth", 3); hold on
plot(k2, fun_lambda(k2)(2,:), "-b", "LineWidth", 3); grid on
xlabel("k_2")
ylabel("\lambda")
saveas(1, "critical_damping.png")

################################################################################
## plots
################################################################################

h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
[c h] = contourf(K1, K2, abs(L1), 0:0.2:2); colorbar();
clabel (c, h, 0:0.2:2, "fontsize", 12);
xlabel("k1"); ylabel("k2"); title("The first pose's magnitude");
grid on
if !visualize
  print(sprintf("%s/first_pole_abs.png", saving_foler));
end

h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
[c h] = contourf(K1, K2, abs(L2), 0:0.2:2); colorbar();
clabel (c, h, 0:0.2:2, "fontsize", 12);
xlabel("k1"); ylabel("k2"); title("The second pose's magnitude");
grid on
if !visualize
  print(sprintf("%s/second_pole_abs.png", saving_foler));
end

h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
[c, h] =  contourf(K1, K2, arg(L2));
colorbar(); xlabel("k1"); ylabel("k2");
title("The second pose's angle");
grid on;
if !visualize
  print(sprintf("%s/second_pole_arg.png", saving_foler));
end

%% Stable region
h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
A1 = abs(L1);
A1(abs_mask) = 3;
h = pcolor(K1, K2, A1); set(h, 'EdgeColor', 'none'); colorbar();
xlabel("k1"); ylabel("k2"); title("Stable region abs(Lambda_1)");
if !visualize
  print(sprintf("%s/stable_region_first_pole_abs.png", saving_foler));
end

h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
A1 = abs(L2);
A1(abs_mask) = 3;
h = pcolor(K1, K2, A1); set(h, 'EdgeColor', 'none'); colorbar();
xlabel("k1"); ylabel("k2"); title("Stable region abs(Lambda_2)");
grid on
if !visualize
  print(sprintf("%s/stable_region_second_pole_abs.png", saving_foler));
end

h = figure();
if !visualize
  set(h, 'visible','off')
end
clf
A1 = arg(L2);
A1(abs_mask) = pi;
h = pcolor(K1, K2, A1); set(h, 'EdgeColor', 'none'); colorbar();
xlabel("k1"); ylabel("k2"); title("Stable region arg(Lambda_2)");
if !visualize
  print(sprintf("%s/stable_region_second_pole_arg.png", saving_foler));
end

h = figure();
plot_stable_region(lambda1, lambda2, dt_s, visualize);
if !visualize
  print(sprintf("%s/stable_region.png", saving_foler));
end

h = figure();
dt_s = 0.2;
k1s = 100;
k2s = 0:0.5:21;
plot_multiple_poles_gains(lambda1, lambda2, dt_s, k1s, k2s, visualize, 'dt02_k100_');

k1s = 20;
k2s = 0:0.5:21;
plot_multiple_poles_gains(lambda1, lambda2, dt_s, k1s, k2s, visualize, 'dt02_k20_');

k1s = 10;
k2s = 0:0.5:21;
plot_multiple_poles_gains(lambda1, lambda2, dt_s, k1s, k2s, visualize, 'dt02_k10_');


%% Plot stable region changed with time
disp("Plot stable region changed with dt");
index = 0;
for dt_s = 0.01:0.01:0.2
  plot_stable_region(lambda1, lambda2, dt_s, false);
  print(sprintf("%s/stable_region_%05d.png", saving_foler, index));
  index = index + 1;
end


