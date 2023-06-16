clearvars -except K; clc;  

global Kp v_dc_ref v_ac_ref 
global Power1 Power2 Power3 

omega  = 2*pi*50;           % reference omega, f = 50 Hz
jw =    [0      -omega;     % cross coupiling dq matrix
         omega   0 ];       % cross coupiling dq matrix

Power1 =  1e6;  
Power2 =  1e6;     
Power3 =  1e6;     


v_g = [563 0]';
v_dc_ref = 1100;
v_ac_ref = 563;
fs = 2000;
Ts = 1/fs;


%% Simulation

x_vsc1 =  [  0  0    1100    0    0    0    563   0   ]';   % initial state
x_vsc2 =  [  0  0    1100    0    0    0    563   0   ]';   % initial state
x_vsc3 =  [  0  0    1100    0    0    0    563   0   ]';   % initial state
x_PN   =  [  0  0   0  0   0  0  ]';   % initial state

x0 = [  x_vsc1  ;  x_vsc2  ;  x_vsc3  ;  x_PN ]';   % initial state
tspan = [0 1.7];    % simulation time


x0 = 1.0e+03 *   [ 1.6787
   -0.6657
    1.1000
    0.0023
    0.0010
   -0.0000
    0.4964
    0.0000
    1.5104
   -0.1128
    1.1000
    0.0020
    0.0011
   -0.0000
    0.5517
   -0.0000
    1.4409
    0.1535
    1.1000
    0.0019
    0.0011
   -0.0000
    0.5783
   -0.0000
    2.8237
    3.7043
    1.3960
    2.5985
    0.5559
    1.3381 ]; % initial state with simulation flat start (to be computed by static analysis)


disp('Simulation running...');
disp(' ')
[t,x] = ode15s(@odefile, tspan, x0);    % run the transient simulation
x(length(t),:);

disp('Simulation completed!');
disp(' ')
disp('Loading plots...');


%% State definition and prelevation - Converter 1
i_c1    = x(:,1:2);
v_dc1   = x(:,3);
phi_dc1 = x(:,4);
delta1  = x(:,5);
phi_q1  = x(:,6);
v_cc1 = x(:,7:8);

%% State definition and prelevation - Converter 2
i_c2    = x(:,9:10);
v_dc2   = x(:,11);
phi_dc2 = x(:,12);
delta2  = x(:,13);
phi_q2  = x(:,14);
v_cc2 = x(:,15:16);

%% State definition and prelevation - Converter 3
i_c3    = x(:,17:18);
v_dc3   = x(:,19);
phi_dc3 = x(:,20);
delta3  = x(:,21);
phi_q3  = x(:,22);
v_cc3 = x(:,23:24);

%% State definition and prelevation - Power Network
i_pccDQ  = x(:,25:26);
il_12DQ  = x(:,27:28);
il_23DQ  = x(:,29:30);


i_g1DQ = i_pccDQ - il_12DQ;
i_g2DQ = il_12DQ - il_23DQ;
i_g3DQ = il_23DQ;

for k = 1 : length(delta1)
    i_g1(:,k) = [cos(delta1(k)) sin(delta1(k)); -sin(delta1(k)) cos(delta1(k))] * i_g1DQ(k,:)';  % From global to local
    i_g2(:,k) = [cos(delta2(k)) sin(delta2(k)); -sin(delta2(k)) cos(delta2(k))] * i_g2DQ(k,:)';  % From global to local
    i_g3(:,k) = [cos(delta3(k)) sin(delta3(k)); -sin(delta3(k)) cos(delta3(k))] * i_g3DQ(k,:)';  % From global to local
    
end

i_g1 = i_g1';
i_g2 = i_g2';
i_g3 = i_g3';


v_g1 = Kp*(i_c1-i_g1) + v_cc1;
v_g2 = Kp*(i_c2-i_g2) + v_cc2; 
v_g3 = Kp*(i_c3-i_g3) + v_cc3; 

for k = 1 : length(delta1)
    i_pcc(:,k) = [cos(delta1(k)) sin(delta1(k)); -sin(delta1(k)) cos(delta1(k))] * i_pccDQ(k,:)';  % From global to local
end

v_pcc = v_g1;
i_pcc = i_pcc';


%% Display time-domain simulation results

% % 
% % 
% figure(1); 
% plot(t,v_dc); grid on; hold on;
% %
% figure(2); 
% subplot(2,1,1);
% plot(t,v_pcc(:,1)); grid on; hold on;
% subplot(2,1,2);
% plot(t,v_pcc(:,2)); grid on; hold on;
% %
% figure(3); 
% subplot(2,1,1);
% plot(t,i_pcc(:,1)); grid on; hold on;
% subplot(2,1,2);
% plot(t,i_pcc(:,2)); grid on; hold on;


figure(14); set(gca,'FontSize',24)
subplot(3,1,2); set(gca,'FontSize',24)
plot(t,v_pcc(:,1),'Linewidth',1.5); grid on; hold on; ylabel('v_{pcc_d} (V)');  axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
subplot(3,1,3); set(gca,'FontSize',24)
plot(t,v_pcc(:,2),'Linewidth',1.5); grid on; hold on; ylabel('v_{pcc_q} (V)'); xlabel('Time (s)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
set(gca,'FontSize',24)
figure(25); set(gca,'FontSize',24)
subplot(3,1,1); set(gca,'FontSize',24)
plot(t,v_dc1,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc_1} (V)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
subplot(3,1,2); set(gca,'FontSize',24)
plot(t,v_dc2,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc_2} (V)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
subplot(3,1,3); set(gca,'FontSize',24)
plot(t,v_dc3,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc_3} (V)'); xlabel('Time (s)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
set(gca,'FontSize',24)

figure(300); set(gca,'FontSize',24)
subplot(2,1,1); set(gca,'FontSize',24)
plot(t,i_g1(:,1), 'Linewidth',1.5); grid on; hold on; ylabel('i_{pcc_d} (A)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
subplot(2,1,2); set(gca,'FontSize',24)
plot(t,i_g1(:,2), 'Linewidth',1.5); grid on; hold on; ylabel('i_{pcc_q} (A)'); xlabel('Time (s)'); axis([0.9 1.4 -inf inf]); yticklabels('auto'); xticklabels('auto');
set(gca,'FontSize',24)

% figure(15); 
% subplot(3,1,1);
% plot(t,v_dc1,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc1} (V)');
% subplot(3,1,2);
% plot(t,i_g1(:,1), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g1_d} (A)');
% subplot(3,1,3); 
% plot(t,i_g1(:,2), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g1_q} (A)'); xlabel('Time (s)');
% 
% figure(24); 
% subplot(3,1,2);
% plot(t,v_g2(:,1),'Linewidth',1.5); grid on; hold on; ylabel('v_{g2_d} (V)');
% subplot(3,1,3);
% plot(t,v_g2(:,2),'Linewidth',1.5); grid on; hold on; ylabel('v_{g2_q} (V)'); xlabel('Time (s)');
% 
% figure(25); 
% subplot(3,1,1);
% plot(t,v_dc2,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc} (V)');
% subplot(3,1,2);
% plot(t,i_g2(:,1), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g2_d} (A)');
% subplot(3,1,3); 
% plot(t,i_g2(:,2), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g2_q} (A)'); xlabel('Time (s)');
% 
% 
% figure(34); 
% subplot(3,1,2);
% plot(t,v_g3(:,1),'Linewidth',1.5); grid on; hold on; ylabel('v_{g3_d} (V)');
% subplot(3,1,3);
% plot(t,v_g3(:,2),'Linewidth',1.5); grid on; hold on; ylabel('v_{g3_q} (V)'); xlabel('Time (s)');
% 
% figure(35); 
% subplot(3,1,1);
% plot(t,v_dc3,'Linewidth',1.5); grid on; hold on; ylabel('v_{dc3} (V)');
% subplot(3,1,2);
% plot(t,i_g3(:,1), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g3_d} (A)');
% subplot(3,1,3); 
% plot(t,i_g3(:,2), 'Linewidth',1.5); grid on; hold on; ylabel('i_{g3_q} (A)'); xlabel('Time (s)');
