close all; clear all;

coords = csvread('clamped-1-refe.csv');
x = coords(:, 1);
y = coords(:, 2);
z = coords(:, 3);
ind = find(x == min(x) & y == min(y));
pe = []; ke = [];
%pe2 = []; ke2 = [];
%K = csvread("stiffness.csv");
%M = csvread("mass.csv");
dispz1 = []; veloz1 = []; accez1 = [];
dispz2 = []; veloz2 = []; accez2 = [];
disp1_computed = []; velo1_computed = []; acce1_computed = [];
disp2_computed = []; velo2_computed = []; acce2_computed = [];
%fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig = figure();
save_figs = 0;
ctr = 1;
for i=0:100:10000
  if (i < 10)
    %pe_file_name = strcat('clamped-1-potential-000', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-000', num2str(i), '.csv');
    disp1_file_name = strcat('clamped-1-disp-000', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-000', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-000', num2str(i), '.csv');
    time1_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-time-000', num2str(i), '.csv');
    disp2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-000', num2str(i), '.csv');
    velo2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-000', num2str(i), '.csv');
    acce2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-000', num2str(i), '.csv');
  elseif (i < 100)
    %pe_file_name = strcat('clamped-1-potential-00', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-00', num2str(i), '.csv');
    disp1_file_name = strcat('clamped-1-disp-00', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-00', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-00', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-00', num2str(i), '.csv');
    disp2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-00', num2str(i), '.csv');
    velo2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-00', num2str(i), '.csv');
    acce2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-00', num2str(i), '.csv');
  elseif (i < 1000)
    %pe_file_name = strcat('clamped-1-potential-0', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-0', num2str(i), '.csv');
    disp1_file_name = strcat('clamped-1-disp-0', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-0', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-0', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-0', num2str(i), '.csv');
    disp2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-0', num2str(i), '.csv');
    velo2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-0', num2str(i), '.csv');
    acce2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-0', num2str(i), '.csv');
  else
    %pe_file_name = strcat('clamped-1-potential-', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-', num2str(i), '.csv');
    disp1_file_name = strcat('clamped-1-disp-', num2str(i), '.csv');
    velo1_file_name = strcat('clamped-1-velo-', num2str(i), '.csv');
    acce1_file_name = strcat('clamped-1-acce-', num2str(i), '.csv');
    time1_file_name = strcat('clamped-1-time-', num2str(i), '.csv');
    disp2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-disp-', num2str(i), '.csv');
    velo2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-velo-', num2str(i), '.csv');
    acce2_file_name = strcat('../../dynamic-opinf-fom-neohookean/clamped-1-acce-', num2str(i), '.csv');
  end
  %p = dlmread(pe_file_name);
  %k = dlmread(ke_file_name);
  d1 = dlmread(disp1_file_name);
  v1 = dlmread(velo1_file_name);
  ac1 = dlmread(acce1_file_name);
  t1 = dlmread(time1_file_name);
  d2 = dlmread(disp2_file_name);
  v2 = dlmread(velo2_file_name);
  ac2 = dlmread(acce2_file_name);
  dispz1 = [dispz1, d1(:, 3)];
  veloz1 = [veloz1, v1(:, 3)];
  accez1 = [accez1, ac1(:, 3)];
  dispz2 = [dispz2, d2(:, 3)];
  veloz2 = [veloz2, v2(:, 3)];
  accez2 = [accez2, ac2(:, 3)];
  disp1_computed = [disp1_computed, d1(:,3)];
  velo1_computed = [velo1_computed, v1(:,3)];
  acce1_computed = [acce1_computed, ac1(:,3)];
  disp2_computed = [disp2_computed, d2(:,3)];
  velo2_computed = [velo2_computed, v2(:,3)];
  acce2_computed = [acce2_computed, ac2(:,3)];

  subplot(3,1,1);
  ax = gca;
  plot(z(ind), dispz1(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  plot(z(ind), dispz2(ind,ctr), '--g', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-disp');
  hold on;
  title(['displacement snapshot ', num2str(i+1), ' at time = ', num2str(t1)]);
  axis([min(z) max(z) -0.001 0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,2);
  ax = gca;
  plot(z(ind), veloz1(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  plot(z(ind), veloz2(ind,ctr), '--g', 'LineWidth', 2);
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
  plot(z(ind), accez2(ind,ctr), '--g', 'LineWidth', 2);
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

sz = size(disp1_computed);
numerator_disp = 0;
denomenator_disp = 0;
numerator_velo = 0;
denomenator_velo = 0;
numerator_acce = 0;
denomenator_acce = 0;
for i=1:sz(2)
  numerator_disp = numerator_disp + norm(disp1_computed(:,i) - disp2_computed(:,i))^2;
  denomenator_disp = denomenator_disp + norm(disp2_computed(:,i))^2;
  numerator_velo = numerator_velo + norm(velo1_computed(:,i) - velo2_computed(:,i))^2;
  denomenator_velo = denomenator_velo + norm(velo2_computed(:,i))^2;
  if (i > 1)
    numerator_acce = numerator_acce + norm(acce1_computed(:,i) - acce2_computed(:,i))^2;
    denomenator_acce = denomenator_acce + norm(acce2_computed(:,i))^2;
  end
end
dispz_relerr = sqrt(numerator_disp/denomenator_disp);
veloz_relerr = sqrt(numerator_velo/denomenator_velo);
accez_relerr = sqrt(numerator_acce/denomenator_acce);

fprintf('z-disp avg rel error = %e\n', dispz_relerr);
fprintf('z-velo avg rel error = %e\n', veloz_relerr);
fprintf('z-acce avg rel error = %e\n', accez_relerr);
