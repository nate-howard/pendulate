%
% pendulate.m
%
% created on: 21.01.2016
%     author: rungger
%
% see readme file for more information on the pendulate example
%
% you need to run ./pendulate binary first 
%
% so that the files: pendulate_ss.bdd 
%                    pendulate_obst.bdd
%                    pendulate_target.bdd
%                    pendulate_controller.bdd 
% are created
%

function pendulate
clear set
close all



%% simulation
% System parameters
A = 0.1; % Amplitude of vertical driving motion
L = 1; % Length of pendulum
g = 9.81; % Accleration of gravity
gamma = 0.01; % Linear drag term

% Precomputing constant terms
beta = A / L;
alpha = g / L;

% ODE
function dxdt = pendulate_ode(t,x,u)
  dxdt = zeros(2,1);

  dxdt(1)=x(2);
  dxdt(2)=(beta*(u^2)*cos(u*t) - alpha)*sin(x(1)) - gamma*x(2);
end


% target set
% L=[2 0 0; 0 1 0; 0 0 .1];
% c=[9.5; 0.6; 0];

% initial state
x0=[pi, -0.1];

% controller=SymbolicSet('bdd/pendulate_controller.bdd','projection',[1 2 3]);

y=x0;
v=[];
n_iter = 1000;
dt = 0.01;
u = 8 + (1:n_iter)*dt*0.5;
for i = 1:n_iter

  
  % if ( (y(end,:)-c')*L'*L*(y(end,:)'-c)<=1 )
  %   break;
  % end 

  % u=controller.getInputs(y(end,:));
  v=[v; u];
  [t, x]=ode45(@pendulate_ode, [0, dt], y(end,:),[],u(i));

  y=[y; x(end,:)];
end

time = (1:n_iter+1) * dt;

plot_trajectory(y, time, L)

end

function plot_trajectory(x, t, L)
    figure('Position', [100,100,1800,530])
    tiledlayout(1,3, 'TileSpacing','tight')
    
    ax1 = nexttile;
    scatter(x(:,1)/pi, x(:,2), 25, '.', 'CData',t)
    title('Pendulum Phase Portrait')
    xlabel('\theta')
    ylabel('$\dot{\theta}$', 'Interpreter','latex')
    colormap(ax1, winter)
    cb1 = colorbar;
    ylabel(cb1, 'Time (s)')

    ax2 = nexttile;
    scatter(t, x(:,1)/pi, 25, '.', 'CData',abs(x(:,2)))
    title('Pendulum \theta vs Time')
    xlabel('Time')
    ylabel('\theta')
    colormap(ax2, parula)
    cb2 = colorbar;
    ylabel(cb2, 'Angular Speed ($|\dot{\theta}|$)', 'Interpreter','latex')
    
    ax3 = nexttile;
    scatter(0, 0, 50, 'red')
    hold on
    pos = [sin(x(:,1)), -cos(x(:,1))] * L;
    % Color by w
    % scatter(pos(:,1), pos(:,2), 25, '.', 'CData',abs(x(:,2)))
    % colormap(ax3, parula)
    % cb3 = colorbar;
    % ylabel(cb3, 'Angular Speed ($|\dot{\theta}|$)', 'Interpreter','latex')
    % Color by t
    scatter(pos(:,1), pos(:,2), 25, '.', 'CData',t)
    colormap(ax3, winter)
    cb3 = colorbar;
    ylabel(cb3, 'Time(s)')

    xlim([-L, L])
    ylim([-L, L])
    title('Pendulum Position')
    hold off
end
