function [dxdt] = odefile(t,x)

global Kp_PLL Ki_PLL Kp Ki Kp_DC Ki_DC Kp_AC L_g R_g v_dc_ref v_ac_ref T_c C_dc e
global Power1 Power2 Power3 

Power  = Power1 + Power2 + Power3;

omega  = 2*pi*50;           % reference omega, f = 50 Hz
jw =    [0      -omega;     % cross coupiling dq matrix
         omega   0 ];       % cross coupiling dq matrix

e = [563 0]';

%% State definition and prelevation - Converter 1
i_c1    = x(1:2);
v_dc1   = x(3);
phi_dc1 = x(4);
delta1  = x(5);
phi_q1  = x(6);
v_cc1   = x(7:8);

%% State definition and prelevation - Converter 2
i_c2    = x(9:10);
v_dc2   = x(11);
phi_dc2 = x(12);
delta2  = x(13);
phi_q2  = x(14);
v_cc2 = x(15:16);

%% State definition and prelevation - Converter 3
i_c3    = x(17:18);
v_dc3   = x(19);
phi_dc3 = x(20);
delta3  = x(21);
phi_q3  = x(22);
v_cc3 = x(23:24);

%% State definition and prelevation - Power Network
i_pccDQ  = x(25:26);
il_12DQ  = x(27:28);
il_23DQ  = x(29:30);


%% Parameters declaration

omega  = 2*pi*50;           % reference omega, f = 50 Hz
jw =    [0      -omega;     % cross coupiling dq matrix
         omega   0 ];       % cross coupiling dq matrix

fs = 2000;
Ts = 1/fs;

Power1 = 1e6;  % 1 MW
Power2 = 1e6;  % 1 MW
Power3 = 1e6;  % 1 MW
Power  = Power1 + Power2 + Power3;

v_dc_ref = 1100;
v_ac_ref = 563;
C_dc = 22e-3;

% Distribution lines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ll_12 = 0.01 / (2*pi*50);    
Rl_12 = 2*pi*50 * Ll_12 * 2.5;

Ll_23 = 0.02 / (2*pi*50);   
Rl_23 = 2*pi*50 * Ll_23 * 1.5;

% PLL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLL_D = 1/sqrt(2);
PLL_Tsettle = 0.2;  % original 0.2
Vamp = 563;

Kp_PLL = 2*4.6 /PLL_Tsettle;
Ti_PLL = PLL_Tsettle*(PLL_D^2) / 2.3;
Ki_PLL = Kp_PLL/Ti_PLL;

Kp_PLL = Kp_PLL / Vamp;
Ki_PLL = Ki_PLL / Vamp;

% Current Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lf = 0.1e-3;    % L Filter
Rf = 0.33*2*pi*50*Lf;

Kp  = Lf*1.8/(3*Ts);
Ti  = Lf/Rf;
Ki  = Kp/Ti;
T_c = Kp/Lf;

% Dc and AC voltage Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kp_DC = 0.12 * C_dc/Ts * 1.8;   % DC Voltage loop
Ti_DC = 17*Ts;                  % DC Voltage loop
Ki_DC = Kp_DC/Ti_DC / 4;        % DC Voltage loop

Kp_AC = -25; % AC Voltage droop


% Grid Impedance

SCR_des = 3.4; 
RXratio = 0.3;

Z_g = 690^2 / (Power * SCR_des);
X_g = sqrt(   Z_g^2 / ( RXratio^2 + 1 )   );

L_g = X_g / omega;
R_g = RXratio * X_g;


%% Disturbance Simulation - Voltage sag

phase_jump = -0*pi/9;
sag = 0.8;
v_after = [cos(phase_jump) sin(phase_jump); -sin(phase_jump) cos(phase_jump)] * [563*sag ; 0] ; % SAG + PHASE JUMP
% 
% % %%%%%%%%%% 
if(t > 1) % at t=1 the disturbance occurs
      e = v_after; % voltage sag
%       Power3 = 1e6;
end

%% Power network nodal equations, according to the Kirkhoff laws

i_g1DQ = i_pccDQ - il_12DQ;
i_g2DQ = il_12DQ - il_23DQ;
i_g3DQ = il_23DQ;

%% Equations Converter 1

iq_ref1 = 0;
Tdelta1 = [cos(delta1) sin(delta1); -sin(delta1) cos(delta1)]; % from DQ to dq

v_g1 = Kp*(i_c1-Tdelta1*i_g1DQ) + v_cc1;

i_ref1 =  [  Kp_DC*(v_dc1 - v_dc_ref)  + Ki_DC*phi_dc1 ;
             ( iq_ref1 + Kp_AC*(v_ac_ref - sqrt(v_g1'*v_g1)) ) ];   %%%%%%%%  Voltage controller
        
wcc = Kp / Lf;
Tdelta1 = [cos(delta1) sin(delta1); -sin(delta1) cos(delta1)]; % from DQ to dq
        
di_c1    = -wcc *i_c1 + wcc *i_ref1;
dv_dc1   = -3/2 * 1/C_dc * v_g1'*Tdelta1*i_g1DQ/v_dc1   +   1/C_dc*Power1 / v_dc1;
dphi_dc1 = v_dc1 - v_dc_ref;
ddelta1  = [0 Kp_PLL]*v_g1   +   Ki_PLL*phi_q1;
dphi_q1  = [0 1] * v_g1;
dv_cc1   = -Ki*Tdelta1*i_g1DQ + Ki*i_c1 ;

v_g1DQ   = Tdelta1'*v_g1;

dx1 = [ di_c1 ; dv_dc1 ; dphi_dc1 ; ddelta1 ; dphi_q1 ; dv_cc1 ];

   
%% Equations Converter 2

iq_ref2 = 0;
Tdelta2 = [cos(delta2) sin(delta2); -sin(delta2) cos(delta2)]; % from DQ to dq

v_g2 = Kp*(i_c2-Tdelta2*i_g2DQ) + v_cc2;

i_ref2 =  [    Kp_DC*(v_dc2 - v_dc_ref)  + Ki_DC*phi_dc2 ;
               ( iq_ref2 + Kp_AC*(v_ac_ref - sqrt(v_g2'*v_g2)) ) ];   %%%%%%%% % Voltage controller


        
di_c2    = -wcc *i_c2 + wcc *i_ref2;
dv_dc2   = -3/2 * 1/C_dc * v_g2'*Tdelta2*i_g2DQ/v_dc2   +   1/C_dc*Power2 / v_dc2;
dphi_dc2 = v_dc2 - v_dc_ref;
ddelta2  = [0 Kp_PLL]*v_g2   +   Ki_PLL*phi_q2;
dphi_q2  = [0 1] * v_g2;
dv_cc2   = -Ki*Tdelta2*i_g2DQ + Ki*i_c2;

v_g2DQ = Tdelta2'*v_g2;

dx2 = [ di_c2 ; dv_dc2 ; dphi_dc2 ; ddelta2 ; dphi_q2 ; dv_cc2 ];

%% Equations Converter 3

iq_ref3 = -500;
Tdelta3 = [cos(delta3) sin(delta3); -sin(delta3) cos(delta3)]; % from DQ to dq

v_g3 = Kp*(i_c3-Tdelta3*i_g3DQ) + v_cc3;

i_ref3 =  [    Kp_DC*(v_dc3 - v_dc_ref)  + Ki_DC*phi_dc3 ;
               ( iq_ref3 + Kp_AC*(v_ac_ref - sqrt(v_g3'*v_g3)) )  ];   %%%%%%%% % Voltage controller
        
        
di_c3    = -wcc *i_c3 + wcc *i_ref3;
dv_dc3   = -3/2 * 1/C_dc * v_g3'*Tdelta3*i_g3DQ/v_dc3   +   1/C_dc*Power3 / v_dc3;
dphi_dc3 = v_dc3 - v_dc_ref;
ddelta3  = [0 Kp_PLL]*v_g3   +   Ki_PLL*phi_q3;
dphi_q3  = [0 1] * v_g3;
dv_cc3 = -Ki*Tdelta3*i_g3DQ + Ki*i_c3 ;

v_g3DQ = Tdelta3'*v_g3;


dx3 = [ di_c3 ; dv_dc3 ; dphi_dc3 ; ddelta3 ; dphi_q3 ; dv_cc3 ];
   
%% Power network and grid equations, according to the Kirkhoff laws

di_pccDQ  = -(R_g)/(L_g)*i_pccDQ  - jw*i_pccDQ  + 1/(L_g) * (v_g1DQ - e); % Grid 
dil_12DQ  = -(Rl_12)/(Ll_12)*il_12DQ - jw*il_12DQ + 1/(Ll_12) * (v_g2DQ - v_g1DQ); % Transmission line 1-2
dil_23DQ  = -(Rl_23)/(Ll_23)*il_23DQ - jw*il_23DQ + 1/(Ll_23) * (v_g3DQ - v_g2DQ); % Transmission line 2-3

dxPN = [ di_pccDQ  ; dil_12DQ ; dil_23DQ ];  % states concatenation

dxdt = [ dx1 ; dx2 ; dx3 ; dxPN ];



