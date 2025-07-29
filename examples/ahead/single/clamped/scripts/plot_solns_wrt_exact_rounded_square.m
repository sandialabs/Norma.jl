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
dispz = []; veloz = []; accez = [];
disp_computed = []; velo_computed = []; acce_computed = [];
disp_exact = []; velo_exact = []; acce_exact = [];
%fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig = figure();
save_figs = 1;
ctr = 1;
for i=0:100:10000
  if (i < 10)
    %pe_file_name = strcat('clamped-1-potential-000', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-000', num2str(i), '.csv');
    disp_file_name = strcat('clamped-1-disp-000', num2str(i), '.csv');
    velo_file_name = strcat('clamped-1-velo-000', num2str(i), '.csv');
    acce_file_name = strcat('clamped-1-acce-000', num2str(i), '.csv');
    time_file_name = strcat('clamped-1-time-000', num2str(i), '.csv');
  elseif (i < 100)
    %pe_file_name = strcat('clamped-1-potential-00', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-00', num2str(i), '.csv');
    disp_file_name = strcat('clamped-1-disp-00', num2str(i), '.csv');
    velo_file_name = strcat('clamped-1-velo-00', num2str(i), '.csv');
    acce_file_name = strcat('clamped-1-acce-00', num2str(i), '.csv');
    time_file_name = strcat('clamped-1-time-00', num2str(i), '.csv');
  elseif (i < 1000)
    %pe_file_name = strcat('clamped-1-potential-0', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-0', num2str(i), '.csv');
    disp_file_name = strcat('clamped-1-disp-0', num2str(i), '.csv');
    velo_file_name = strcat('clamped-1-velo-0', num2str(i), '.csv');
    acce_file_name = strcat('clamped-1-acce-0', num2str(i), '.csv');
    time_file_name = strcat('clamped-1-time-0', num2str(i), '.csv');
  else
    %pe_file_name = strcat('clamped-1-potential-', num2str(i), '.csv');
    %ke_file_name = strcat('clamped-1-kinetic-', num2str(i), '.csv');
    disp_file_name = strcat('clamped-1-disp-', num2str(i), '.csv');
    velo_file_name = strcat('clamped-1-velo-', num2str(i), '.csv');
    acce_file_name = strcat('clamped-1-acce-', num2str(i), '.csv');
    time_file_name = strcat('clamped-1-time-', num2str(i), '.csv');
  end
  %p = dlmread(pe_file_name);
  %k = dlmread(ke_file_name);
  d = dlmread(disp_file_name);
  v = dlmread(velo_file_name);
  ac = dlmread(acce_file_name);
  t = dlmread(time_file_name);
  dispz = [dispz, d(:, 3)];
  veloz = [veloz, v(:, 3)];
  accez = [accez, ac(:, 3)];
  c = sqrt(1e9/1e3);
  T = 1e-3;
  a=0.0005; b=100; s=0.6;
  d_ex = a/2 * (tanh(-b*(z+0.5-s-c*t)) + tanh(b*(z-1+0.5+s-c*t))) + ...
         a/2 * (tanh(-b*(z+0.5-s+c*t)) + tanh(b*(z-1+0.5+s+c*t))) - ...
         a/2 * (tanh(-b*(z+0.5-s-c*(T-t))) + tanh(b*(z-1+0.5+s-c*(T-t)))) - ...
         a/2 * (tanh(-b*(z+0.5-s+c*(T-t))) + tanh(b*(z-1+0.5+s+c*(T-t))));
  v_ex = a/2*((1-tanh(b*(z+.5-s-c*t)).^2)*b*c-(1-tanh(b*(z-.5+s-c*t)).^2)*b*c)+a/2*(-(1-tanh(b*(z+.5-s+c*t)).^2)*b*c+(1-tanh(b*(z-.5+s+c*t)).^2)*b*c)-a/2*(-(1-tanh(b*(z+.5-s-c*(T-t))).^2)*b*c+(1-tanh(b*(z-.5+s-c*(T-t))).^2)*b*c)-a/2*((1-tanh(b*(z+.5-s+c*(T-t))).^2)*b*c-(1-tanh(b*(z-.5+s+c*(T-t))).^2)*b*c);
  a_ex = 1/2*a*(2*tanh(b*(z+.5-s-c*t)).*(1-tanh(b*(z+.5-s-c*t)).^2)*b^2*c^2-2*tanh(b*(z-.5+s-c*t)).*(1-tanh(b*(z-.5+s-c*t)).^2)*b^2*c^2)+1/2*a*(2*tanh(b*(z+.5-s+c*t)).*(1-tanh(b*(z+.5-s+c*t)).^2)*b^2*c^2-2*tanh(b*(z-.5+s+c*t)).*(1-tanh(b*(z-.5+s+c*t)).^2)*b^2*c^2)-1/2*a*(2*tanh(b*(z+.5-s-c*(T-t))).*(1-tanh(b*(z+.5-s-c*(T-t))).^2)*b^2*c^2-2*tanh(b*(z-.5+s-c*(T-t))).*(1-tanh(b*(z-.5+s-c*(T-t))).^2)*b^2*c^2)-1/2*a*(2*tanh(b*(z+.5-s+c*(T-t))).*(1-tanh(b*(z+.5-s+c*(T-t))).^2)*b^2*c^2-2*tanh(b*(z-.5+s+c*(T-t))).*(1-tanh(b*(z-.5+s+c*(T-t))).^2)*b^2*c^2);
  disp_computed = [disp_computed, d(:,3)];
  velo_computed = [velo_computed, v(:,3)];
  acce_computed = [acce_computed, ac(:,3)];
  disp_exact = [disp_exact, d_ex];
  velo_exact = [velo_exact, v_ex];
  acce_exact = [acce_exact, a_ex];

  subplot(3,1,1);
  ax = gca;
  plot(z(ind), dispz(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  plot(z(ind), d_ex(ind), '--g', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-disp');
  hold on;
  title(['displacement snapshot ', num2str(i+1), ' at time = ', num2str(t)]);
  axis([min(z) max(z) -0.001 0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,2);
  ax = gca;
  plot(z(ind), veloz(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  plot(z(ind), v_ex(ind), '--g', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-velo');
  hold on;
  title(['velocity snapshot ', num2str(i+1), ' at time = ', num2str(t)]);
  axis([min(z) max(z) -3e4*0.001 3e4*0.001]);
  ax.NextPlot = 'replaceChildren';
  subplot(3,1,3);
  ax = gca;
  plot(z(ind), accez(ind,ctr), '-b', 'LineWidth', 2);
  hold on;
  plot(z(ind), a_ex(ind), '--g', 'LineWidth', 2);
  xlabel('z');
  ylabel('z-acce');
  title(['acceleration snapshot ', num2str(i+1), ' at time = ', num2str(t)]);
  axis([min(z) max(z) -4e9*0.001 4e9*0.001]);
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

sz = size(disp_exact);
numerator_disp = 0;
denomenator_disp = 0;
numerator_velo = 0;
denomenator_velo = 0;
numerator_acce = 0;
denomenator_acce = 0;
for i=1:sz(2)
  numerator_disp = numerator_disp + norm(disp_computed(:,i) - disp_exact(:,i))^2;
  denomenator_disp = denomenator_disp + norm(disp_exact(:,i))^2;
  numerator_velo = numerator_velo + norm(velo_computed(:,i) - velo_exact(:,i))^2;
  denomenator_velo = denomenator_velo + norm(velo_exact(:,i))^2;
  if (i > 1)
    numerator_acce = numerator_acce + norm(acce_computed(:,i) - acce_exact(:,i))^2;
    denomenator_acce = denomenator_acce + norm(acce_exact(:,i))^2;
  end
end
dispz_relerr = sqrt(numerator_disp/denomenator_disp);
veloz_relerr = sqrt(numerator_velo/denomenator_velo);
accez_relerr = sqrt(numerator_acce/denomenator_acce);

fprintf('z-disp avg rel error = %e\n', dispz_relerr);
fprintf('z-velo avg rel error = %e\n', veloz_relerr);
fprintf('z-acce avg rel error = %e\n', accez_relerr);
