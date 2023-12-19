function [p_amb,T_amb]=atm_standard(z) 
% la quota z va data in metri
p_0 = 101325;
T_0 = 288.15;
gamma_a=1.4;
cp_a=1004.5;
R_a=cp_a-cp_a/gamma_a;
p_02 = 22632;
T_02 = 216.65;
if z <= 11000
    p_amb = p_0 * ((T_0 - 0.0065 * z)/T_0)^(9.80665/(R_a*0.0065));
    T_amb = T_0 - 0.0065*z;
else 
    p_amb = p_02 * exp(-9.80665*(z-11000)/(R_a * T_02));
    T_amb = T_02;
end