function a = attack_rate(b_a,T)
global T_min T_max e_a
a = b_a*(T-T_min).*(T_max-T).^e_a;
end