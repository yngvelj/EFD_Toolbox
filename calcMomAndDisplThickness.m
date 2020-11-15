%% calcMomAndDisplThickness -   Calculates momentum and displacement
%%                              thickness based on boundary layer profile
%==========================================================================
% y =   distance to wall, can be flipped both ways as long as it is 
%       starting from y = 0 or ending at y=0
% u =   Corresponding velocity
% U =   Freestream velocity
%==========================================================================
function [delta_star, theta] = calcMomAndDisplThickness(y, u, U, fig)
% Change direction so that we start from y = 0 if vectors are wrong
[startval,start] = min(y);
stop = max(y);
if start == length(y) %Flip vectors
    y = flip(y);
    u = flip(u);
end
% Evaluate integrands
integrand_delta_star = 1-(u./U);
integrand_theta = (u./U).*(1-(u./U));
% Calculate integral
delta_star = trapz(y,integrand_delta_star);
theta = trapz(y,integrand_theta);
% Plot integrands
figure(fig);
hold on;
area(y,integrand_delta_star,"FaceColor", "b","FaceAlpha",0.1);
area(y,integrand_theta,"FaceColor", "r","FaceAlpha",0.5);
hold off;
grid on;
axis([0, stop, 0,1]);
str = sprintf("Momentum and displacement thicknesses for U_{\\infty} = %1.2f m/s",U);
title(str);
xlabel("y [mm]");
lgd1 = "$1-\frac{u}{U}$";
lgd2 = "$\frac{u}{U}(1-\frac{u}{U})$";
legend(lgd1,lgd2,"interpreter","latex");
% Show values in plot
str_theta = strcat("\theta = ",num2str(theta)," mm");
[max_theta, j] = max(integrand_theta);
x_theta = [0.6,0.25];
y_theta = [0.3,0.15];

str_delta_star = strcat("\delta^{*} = ",num2str(delta_star)," mm");
x_delta_star = [0.3,0.15];
y_delta_star = [0.7,0.5];

annotation("textarrow",x_delta_star,y_delta_star,"String",str_delta_star);
annotation("textarrow",x_theta,y_theta,"String",str_theta);


end



