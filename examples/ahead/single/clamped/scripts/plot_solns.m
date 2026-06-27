close all; clear all;

coords = csvread('clamped_single_gaussian-refe.csv');
x = coords(:, 1);
y = coords(:, 2);
z = coords(:, 3);
ind = find(x == min(x) & y == min(y));
pe = []; ke = [];
%pe2 = []; ke2 = [];
%K = csvread("stiffness.csv");
%M = csvread("mass.csv");
dispz1 = []; veloz1 = []; accez1 = [];
%fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig = figure();
save_figs = 0;
ctr = 1;
for i=0:10:1000
  if (i < 10)
    %pe_file_name = strcat('clamped_single_gaussian-potential-000', num2str(i), '.csv');
    %ke_file_name = strcat('clamped_single_gaussian-kinetic-000', num2str(i), '.csv');
    disp1_file_name = strcat('clamped_single_gaussian-disp-000', num2str(i), '.csv');
    velo1_file_name = strcat('clamped_single_gaussian-velo-000', num2str(i), '.csv');
    acce1_file_name = strcat('clamped_single_gaussian-acce-000', num2str(i), '.csv');
    time1_file_name = strcat('clamped_single_gaussian-time-000', num2str(i), '.csv');
  elseif (i < 100)
    %pe_file_name = strcat('clamped_single_gaussian-potential-00', num2str(i), '.csv');
    %ke_file_name = strcat('clamped_single_gaussian-kinetic-00', num2str(i), '.csv');
    disp1_file_name = strcat('clamped_single_gaussian-disp-00', num2str(i), '.csv');
    velo1_file_name = strcat('clamped_single_gaussian-velo-00', num2str(i), '.csv');
    acce1_file_name = strcat('clamped_single_gaussian-acce-00', num2str(i), '.csv');
    time1_file_name = strcat('clamped_single_gaussian-time-00', num2str(i), '.csv');
  elseif (i < 1000)
    %pe_file_name = strcat('clamped_single_gaussian-potential-0', num2str(i), '.csv');
    %ke_file_name = strcat('clamped_single_gaussian-kinetic-0', num2str(i), '.csv');
    disp1_file_name = strcat('clamped_single_gaussian-disp-0', num2str(i), '.csv');
    velo1_file_name = strcat('clamped_single_gaussian-velo-0', num2str(i), '.csv');
    acce1_file_name = strcat('clamped_single_gaussian-acce-0', num2str(i), '.csv');
    time1_file_name = strcat('clamped_single_gaussian-time-0', num2str(i), '.csv');
  else
    %pe_file_name = strcat('clamped_single_gaussian-potential-', num2str(i), '.csv');
    %ke_file_name = strcat('clamped_single_gaussian-kinetic-', num2str(i), '.csv');
    disp1_file_name = strcat('clamped_single_gaussian-disp-', num2str(i), '.csv');
    velo1_file_name = strcat('clamped_single_gaussian-velo-', num2str(i), '.csv');
    acce1_file_name = strcat('clamped_single_gaussian-acce-', num2str(i), '.csv');
    time1_file_name = strcat('clamped_single_gaussian-time-', num2str(i), '.csv');
  end
  %p = dlmread(pe_file_name);
  %k = dlmread(ke_file_name);
  d1 = dlmread(disp1_file_name);
  v1 = dlmread(velo1_file_name);
  ac1 = dlmread(acce1_file_name);
  t1 = dlmread(time1_file_name);
  dispz1 = [dispz1, d1(:, 3)];
  veloz1 = [veloz1, v1(:, 3)];
  accez1 = [accez1, ac1(:, 3)];

  subplot(3,1,1);
  ax = gca;
  plot(z(ind), dispz1(ind,ctr), '-b', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-disp');
  hold on;
  title(['displacement snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -0.001 0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,2);
  ax = gca;
  plot(z(ind), veloz1(ind,ctr), '-b', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-velo');
  hold on;
  title(['velocity snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -3e4*0.001 3e4*0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,3);
  ax = gca;
  plot(z(ind), accez1(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  xlabel('z');
  ylabel('z-acce');
  title(['acceleration snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -2.5e9*0.001 2.5e9*0.001]);
  hold on;
  ax.NextPlot = 'replaceChildren';
  %pause()
  pause(0.5)
  %pe = [pe; p];
  %ke = [ke; k];
 % disp_vec = zeros(length(d)*3,1);
 % disp_vec(1:3:end) = d(:,1);
 % disp_vec(2:3:end) = d(:,2);
 % disp_vec(3:3:end) = d(:,3);
 % velo_vec = zeros(length(v)*3,1);
 % velo_vec(1:3:end) = v(:,1);
 % velo_vec(2:3:end) = v(:,2);
 % velo_vec(3:3:end) = v(:,3);
 % ke2 = [ke2; 0.5*velo_vec'*M*velo_vec];
 % pe2 = [pe2; 0.5*disp_vec'*K*disp_vec];
  if (save_figs == 1)
    if (ctr < 10)
      filename = strcat('soln_000', num2str(ctr), '.png');
      filename2 = strcat('soln_000', num2str(ctr), '.fig');
    elseif (ctr < 100)
      filename = strcat('soln_00', num2str(ctr), '.png');
      filename2 = strcat('soln_00', num2str(ctr), '.fig');
    elseif (ctr < 1000)
      filename = strcat('soln_0', num2str(ctr), '.png');
      filename2 = strcat('soln_0', num2str(ctr), '.fig');
    else
      filename = strcat('soln_', num2str(ctr), '.png');
      filename2 = strcat('soln_', num2str(ctr), '.fig');
    end
    saveas(fig,filename)
    saveas(fig,filename2)
  end
  ctr = ctr + 1;
end
%te = pe + ke;
%te2 = pe2 + ke2;
%figure();
%plot(te);
%hold on;
%plot(te2)
%xlabel('snapshot #');
%ylabel('total energy');

