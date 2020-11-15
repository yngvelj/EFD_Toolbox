%% ****************BOUNDARY LAYER ANALYSIS*********************************
% ------------------------------------------------------------------------
%   Description:
%               Collect PIV-data in order to plot the boundary layer and (yplus uplus) in 
%               front of the BUMP to check for stability.
%Author: Yngve Lilleeng Jenssen
%Date: october 2020
clear all
rho = 1000;
water = water_properties;


%Get PIV-data 7200 series (7201-7203) Boundary layer analysis
piv_bl = readimx('B00001_AvgV.vc7');

velnposin = showimx(piv_bl.Frames{1}); %Struct with position and mean velocity
close all;
% [jmin,ymin] = find_limits(velnpos); %Calculate where the measurement starts("wall")
velnpos = sortArray(velnposin); %Get velocity and position structure
velnpos.U(velnpos.U==0)=nan;    %If, velocity = 0 set value to nan
velnpos.U(velnpos.V==0)=nan;    %If, velocity = 0 set value to nan
U_inf = mean(velnpos.U(100,100:end),'omitnan');

%Calculate total velocity
U_tot = sqrt(velnpos.U.^2+velnpos.V.^2); %calculates the total 2D velocity

%Plot boundary layer profile: Raw-data

l{1} = sprintf('X = 0');
c = 0;
s = 5; %88%Start from iterregation window nr five from left
h = 5;%5
e = 50; %150
figure(2)
title('Boundary layer velocity')
ylabel('Y-position');
xlabel('Velocity u_{bar}/U_0');
hold on;
for i = s:h:e %
    c = c+1;
    x(c) =velnpos.X(i);
    l{c} = sprintf('X = %.1f',x(c)); %legend cell array
    plot(velnpos.U(i,:)/U_inf,abs(velnpos.Y(i,:)+abs(velnpos.Y(i,1))));
    

end
legend(l);
hold off

theta = calculateTheta(velnpos,U_inf); %Calculate momentum Thickness THETA

% theta(theta <2) = nan;
% for i = 1:length(theta)
%    if(theta(i)<2)
%        theta(i)  = [];
%        Xtheta(i) = [];
%    end
% end

theta = theta(s:end-3)*1e-3; %Extract wanted elements from theta vector
Xtheta = velnpos.X(s:end-3,1)*1e-3; %and corresponding X-position values
% theta(theta<1.5) = nan;
figure()
hold on;
plot(Xtheta,theta)
z = fit(Xtheta,(theta)','poly2'); %curve fit theta in order to have a smooth function for differentiation
plot(z);
xlabel('X-position');
ylabel('\Theta Momentum thickness [mm]')
legend('Momentum thickness \Theta','Curve fit');


Cf = 2*(2*z.p1*Xtheta+z.p2); %Calculate Cf from curve-fit
tau_wall = Cf*1/2*rho*U_inf^2; %Calculate tau wall from Cf

figure()
hold on;
plot(Xtheta,tau_wall);
xlabel('X-position');
ylabel('\tau_{wall}');

 %Calculate V_star Friction velocity
figure()
hold on;
c = 0;

for i = s:h:e
    c = c+1;
    u_star = sqrt(tau_wall(i)/rho);
    y_plus = abs(velnpos.Y(i,:)+abs(velnpos.Y(i,1)))*u_star./water.nu(20);
    u_plus = velnpos.U(i,:)./u_star;
    x(c) =Xtheta(i);
    l{c} = sprintf('X = %.1f',x(c));
    plot(y_plus,u_plus);
    

end
legend(l);
xlabel('y^+');
ylabel('u^+');
set(gca, 'XScale', 'log')
hold off;


% Plot the velocity field in a contour plot
figure()
contourf(velnpos.X,velnpos.Y,velnpos.U,30);
c = colorbar;
c.Label.String = 'Velocity u';
for i = 1:length(x)
   xline(x(i)); 
end


figure()
contourf(velnpos.X,velnpos.Y,velnpos.V,30);
c = colorbar;
c.Label.String = 'Velocity v';
for i = 1:length(x)
   xline(x(i)); 
end


%%

function velnpos = sortArray(velnposin) %Get rid of zero-elements
[r,c] = size(velnposin.U);

for i = 1:r
counter = 0;
    for j = 1:c
        if(velnposin.U(i,j) > 0) %Check only x-direction
            counter = counter +1;
            velnpos.U(i,counter) = velnposin.U(i,j);
            velnpos.V(i,counter) = velnposin.V(i,j);
            velnpos.X(i,counter) = velnposin.X(i,j);
            velnpos.Y(i,counter) = velnposin.Y(i,j);
        end
    end
end

end


function theta = calculateTheta(velnpos,Ue)
[r,c] = size(velnpos.X);

    for i = 1:r
        
        integrand = velnpos.U(i,:)/Ue.*(1-velnpos.U(i,:)/Ue);
        theta(i) = trapz(abs(velnpos.Y(i,:)+abs(velnpos.Y(i,1))),max(integrand,0));
    end
end


%% UNUSED
% function [jmin,ymin] = find_limits(velnpos)
% [r,c] = size(velnpos.X);
% foundmin = 0;
% foundmax = 0;
% 
% ymax = 0;
%     for i = 1:r
%         foundmin =0;
%         for j = 1:c
%             if(velnpos.U(i,j) > 0 && foundmin == 0)
%                 ymin(i) = velnpos.Y(i-1,j-1);
%                 jmin(i) = j;
%                 foundmin = 1;
%                 break;
%             end
%         end
%     end
% end