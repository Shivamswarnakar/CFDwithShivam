% Shivam_Swarnakar_184106011_ME704
% MATLAB code for an Explicit Method-based Solution of 2D unsteady state heat convection in a channel
clc;
clear;

L = 6;
H = 1;
rho = 1;
gamma = 0.02;
u = 1;
v = 0;

imax = 42; 
jmax = 22;

delx = L/(imax-1);
dely = H/(jmax-1);

deltd = 0.5*gamma*(1/(delx^2)+1/(dely)^2);
deltc = 1/max((abs(u)/delx),(abs(v)/dely));
delt = min(deltd,deltc);

eps = 0.00001;

%% FOU Scheme
err = 0;
max_err = 1 ;

phi_FU(1:imax,1:jmax) = zeros;
phi_FU(:,:) = 0;
phi_FU(1,:) = 1;

phi_FU_old = phi_FU;

aE = ((gamma*delt)/(delx^2))-(min(u,0)*delt)/delx;
aW = ((gamma*delt)/(delx^2))+(max(u,0)*delt)/delx;

aN = ((gamma*delt)/(dely^2))-(min(v,0)*delt)/dely;
aS = ((gamma*delt)/(dely^2))+(max(v,0)*delt)/dely;
aP = 1; 

omega_FU = 0.3;

n =0;
while (abs(max_err) >= eps)
    for i = 2:imax-1
            phi_FU(imax,:) = phi_FU(imax-1,:);
        for j = 2:jmax-1
            b = (1-(aE+aW+aN+aS))*phi_FU(i,j);
            phi_FU(i,j) = aE*phi_FU(i+1,j) + aW*phi_FU(i-1,j) + aN*phi_FU(i,j+1) + aS*phi_FU(i,j-1) + b;
            phi_FU(i,j) = phi_FU_old(i,j) + omega_FU*(phi_FU(i,j) - phi_FU_old(i,j));
        end
    end
    err = (phi_FU - phi_FU_old);
    max_err = max(max(err));
    phi_FU_old = phi_FU;
    n = n+1;
end

% % Plotting phi_FU

x = linspace(0,6,imax);
y = linspace(0,1,jmax);
contourf(x,y,phi_FU','LineStyle','none');
colormap('jet');
colorbar;
xlabel('x','FontSize',15);
ylabel('y','FontSize',15);
axis equal;
saveas(gcf,'FOU_contour_plot.png');

plot(phi_FU(8,:),y,phi_FU(15,:),y,phi_FU(22,:),y,phi_FU(19,:),y,'LineWidth',2);
xlabel('\theta','FontSize',15);
ylabel('y','FontSize',15);
legend('x/H = 1','x/H = 2','x/H = 3','x/H = 5','Position',[0.35,0.55,0,0]);
title('First order upwind scheme');
set(gca,'FontSize',15);
saveas(gcf,'FOU_Line_plot.png');

%% CD Scheme

err = 0;
max_err = 1 ;

phi_CD(1:imax,1:jmax) = zeros;
phi_CD(:,:) = 0;
phi_CD(1,:) = 1;

phi_CD_old = phi_CD;

aE = ((gamma*delt)/(delx^2))-((u*delt)/(2*delx));
aW = ((gamma*delt)/(delx^2))+((u*delt)/(2*delx));
aN = ((gamma*delt)/(dely^2))-((v*delt)/(2*dely));
aS = ((gamma*delt)/(dely^2))+((v*delt)/(2*dely));
omega_CD = 0.3;

n = 0;
while (abs(max_err) >= eps)
    
    for i = 2:imax-1
        phi_CD(imax,:) = phi_CD(imax-1,:);
        for j = 2:jmax-1
            b = (1-(aE+aW+aN+aS))*phi_CD(i,j);
            phi_CD(i,j) = (aE*phi_CD(i+1,j) + aW*phi_CD(i-1,j) + aN*phi_CD(i,j+1) + aS*phi_CD(i,j-1) + b);
            phi_CD(i,j) = phi_CD_old(i,j) + omega_CD*(phi_CD(i,j) - phi_CD_old(i,j));
        end
    end
    err = (phi_CD - phi_CD_old);
    max_err = max(max(err));
    phi_CD_old = phi_CD;
    n = n+1;
end

% %Plotting phi_CD

x = linspace(0,6,imax);
y = linspace(0,1,jmax);
contourf(x,y,phi_CD','LineStyle','none');
colormap('jet');
colorbar;
xlabel('x','FontSize',15);
ylabel('y','FontSize',15);
axis equal;
saveas(gcf,'CD_contour_plot.png');

plot(phi_CD(8,:),y,phi_CD(15,:),y,phi_CD(22,:),y,phi_CD(19,:),y,'LineWidth',2);
xlabel('\theta','FontSize',15);
ylabel('y','FontSize',15);
legend('x/H = 1','x/H = 2','x/H = 3','x/H = 5','Position',[0.35,0.55,0,0]);
title('Central difference scheme');
set(gca,'FontSize',15);
saveas(gcf,'CD_Line_plot.png');









