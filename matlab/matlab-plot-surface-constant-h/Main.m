clear
clc
global T_min T_max e_a e_h
T_min = 3.8;
T_max = 5.3;
T = T_min+0.1:0.05:T_max-0.1;
B = linspace(0,2468,length(T));

age = 6;

if age == 3
        b_a = 0.00934499611463646 ;
        b_h = 6.47171864976904 ;
        e_a = 0.999999958391371;
elseif age==4
        b_a = 0.00334371401193976 ;
        b_h = 2.84774641530102 ;
        e_a = 0.999999997632107 ;
elseif age==5
        b_a = 0.00283826519189607 ;
        b_h = 1.69775076937433 ;
        e_a = 0.999999999492532 ;
else
        b_a = 0.00386284507963447 ;
        b_h = 1.61430245667652;
        e_a = 0.998788422292012;
end

[X,Y] = meshgrid(T,B);
C = (attack_rate(b_a,X).*Y)./(1+attack_rate(b_a,X).*handleing_time(b_h,X).*Y);
surf(X,Y,C)
xlabel('Temperature (^{\circ}C)', 'Interpreter','tex')
ylabel('Capelin biomass (10^3 t)', 'Interpreter','tex')
zlabel('Consumption (kg/month)', 'Interpreter','tex')

xh = get(gca,'XLabel'); % Handle of the x label
set(xh, 'Units', 'Normalized')
pos = get(xh, 'Position');
set(xh, 'Position',pos.*[1,0,0],'Rotation',15)

yh = get(gca,'YLabel'); % Handle of the y label
set(yh, 'Units', 'Normalized')
pos = get(yh, 'Position');
set(yh, 'Position',pos.*[0.5,-0.5,0],'Rotation',-25)

str = sprintf('Cod age: %d yrs', age);
title(str)