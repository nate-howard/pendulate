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
gamma = 0; % Linear drag term

% Precomputing constant terms
alpha = g / L;
beta = A / L;

% Averaging constants
avg_beta = 0.25 * beta^2;

% ODE
function dxdt = pendulate_ode(t,x,u)
  dxdt = zeros(3,1);

  dxdt(1)=x(2);
  dxdt(2)=(beta*(u^2)*cos(x(3)) - alpha)*sin(x(1)) - gamma*x(2);
  dxdt(3)=u;
end

function dxdt = pendulate_averaged_ode(t,x,u)
  dxdt = zeros(2,1);

  dxdt(1)=x(2);
  dxdt(2)=-alpha*sin(x(1)) - avg_beta*(u^2)*sin(2*x(1));
end

function new_r = pendulate_growth(r,x,u,dt)
    new_r = zeros(size(r));

    new_r(2) = r(2) + r(1)*abs(alpha*cos(x(1)) + 2*u*u*avg_beta*cos(2*x(1)))*dt;
    % new_r(2) = r(2) + r(1)*abs(alpha + 2*u*u*avg_beta)*dt; % Conservative bound
    new_r(1) = r(1) + new_r(2)*dt;
end

% Visualize growth bound
% x0 = [pi,0.8];
% bounds = [0.05, 0.05]/2;
% dt = 0.1;
% u = 50;
% compute_growth(@pendulate_averaged_ode, @pendulate_growth, x0, bounds, dt, u)


% Plot state, target, controller sets
target_set = SymbolicSet('bdd/stay_pendulate_target.bdd');
state_set = SymbolicSet('bdd/stay_pendulate_ss.bdd');
controller = SymbolicSet('bdd/stay_pendulate_controller.bdd', 'projection',[1, 2]);
% t_points = target_set.points;
% s_points = state_set.points;
% c_points = controller.points;
% figure(1)
% hold on;
% scatter(s_points(:,1),s_points(:,2),'.','blue')
% scatter(t_points(:,1)+.01,t_points(:,2)+.01,'.','red')
% scatter(c_points(:,1)-.01,c_points(:,2)-.01,'.','green')
% hold off;


% Averaged ode
x0=[2.8, -0.2];
filename = 'images/driving_showcase.gif';

n_iter = 150;
dt = 0.1;
u = 50;
y=x0;
v=[1];
time = [0];
for i = 1:n_iter

    % Get input from controller
    % try
    %     u=min(controller.getInputs(y(end,:)));
    % catch e
    %     warning('Point not in controller!!!');
    %     y(end,:)
    %     break
    % end
    v=[v; u];

    % Integrate DE's
    [t, x]=ode45(@pendulate_averaged_ode, [time(end), time(end) + dt], y(end,:),[],u);

    % Record result
    last_x = x(end,:);
    last_x(1) = wrapTo2Pi(last_x(1));

    time = [time; t(end)];
    y=[y; last_x];

    if i > n_iter/2
        u = 1;
    end
end

% y(end-10:end,:)

plot_trajectory(y, time, v, target_set, false)

if size(y,1) > n_iter
    animate_trajectory(y, v, L, target_set, dt/2, filename)
end

end

function compute_growth(ode, growth, x0, bounds, dt, u)
    offsets = [1,1; 1,0; 1,-1; 0,1; 0,-1; -1,1; -1,0; -1,-1];
    offset_x = x0 + (offsets.*bounds);
    orig_area = prod(2*bounds)


    [t, x]=ode45(ode, [0, dt], x0, [], u);
    new_x = x(end,:);

    % Integrate offset points
    new_offset_x = zeros(size(offset_x));
    for i = 1:8
        [t, x]=ode45(ode, [0, dt], offset_x(i,:), [], u);
        new_offset_x(i,:) = x(end,:);
    end

    % Convex hull
    [hull_idx, hull_area] = convhull(new_offset_x);
    hull_area

    % Rectangular bounds from offsets
    ub = max(new_offset_x);
    lb = min(new_offset_x);
    rect_area = (ub(1)-lb(1))*(ub(2)-lb(2))

    % Formal growth bounds
    corner_offsets = [1,1; 1,-1; -1,-1; -1,1];
    new_bounds = growth(bounds, new_x, u, dt)
    formal_area = prod(2*new_bounds)
    formal_offset_x = new_x + (corner_offsets.*new_bounds);

    hold on;
    patch(formal_offset_x(:,1), formal_offset_x(:,2), [1,.6,.6])

    % patch([lb(1),ub(1),ub(1),lb(1)], [lb(2),lb(2),ub(2),ub(2)], [.6,1,.6])

    plot(new_offset_x(hull_idx,1), new_offset_x(hull_idx,2))

    scatter(x0(1), x0(2), 25, 'blue')
    scatter(offset_x(:,1), offset_x(:,2), 25, 'magenta')
    scatter(new_x(1), new_x(2), 25, 'blue')
    scatter(new_offset_x(:,1), new_offset_x(:,2), 25, 'red')
    
    xlabel('\theta')
    ylabel('$\dot{\theta}$', 'Interpreter','latex')
    title('Growth Bounds')
    hold off;
end

function plot_trajectory(x, t, v, target_set, center_0)
    if center_0
        x(:,1) = wrapToPi(x(:,1));
    end

    figure('Position', [100,100,530,1130])
    tiledlayout(2,1, 'TileSpacing','tight')
    
    ax1 = nexttile;
    hold on;
    lb1 = min(target_set.points(:,1:2),[],1)' - (target_set.eta(1:2) / 2);
    ub1 = max(target_set.points(:,1:2),[],1)' + (target_set.eta(1:2) / 2);
    patch([lb1(1),ub1(1),ub1(1),lb1(1)]/pi, [lb1(2),lb1(2),ub1(2),ub1(2)], [.6,1,.6])

    % c_points = controller.points;
    % scatter(c_points(:,1)/pi,c_points(:,2),'.','green')

    scatter(x(:,1)/pi, x(:,2), 25, '.', 'CData',t)
    title('Pendulum Phase Portrait')
    xlabel('\theta')
    ylabel('$\dot{\theta}$', 'Interpreter','latex')
    colormap(ax1, winter)
    cb1 = colorbar;
    ylabel(cb1, 'Time (s)')
    hold off;

    ax2 = nexttile;
    hold on;
    lb2 = [0, lb1(1)];
    ub2 = [t(end), ub1(1)];
    patch([lb2(1),ub2(1),ub2(1),lb2(1)], [lb2(2),lb2(2),ub2(2),ub2(2)]/pi, [.6,1,.6])
    scatter(t, x(:,1)/pi, 25, '.', 'CData',v)
    title('Pendulum \theta vs Time')
    xlabel('Time')
    ylabel('\theta')
    colormap(ax2, copper)
    cb2 = colorbar;
    ylabel(cb2, 'Driving Angular Frequency')
    hold off;
    
    % ax3 = nexttile;
    % hold on;
    % patch([0, 1.5*L*sin(lb1(1)), 1.5*L*sin(ub1(1))], ...
    %       [0, -1.5*L*cos(lb1(1)), -1.5*L*cos(ub1(1))], [.6,1,.6])
    % scatter(0, 0, 50, 'red')
    % pos = [sin(x(:,1)), -cos(x(:,1))] * L;
    % % Color by w
    % % scatter(pos(:,1), pos(:,2), 25, '.', 'CData',abs(x(:,2)))
    % % colormap(ax3, parula)
    % % cb3 = colorbar;
    % % ylabel(cb3, 'Angular Speed ($|\dot{\theta}|$)', 'Interpreter','latex')
    % % Color by t
    % scatter(pos(:,1), pos(:,2), 25, '.', 'CData',v)
    % colormap(ax3, copper)
    % cb3 = colorbar;
    % ylabel(cb3, 'Driving Angular Frequency')
    % xlim([-L, L])
    % ylim([-L, L])
    % title('Pendulum Position')
    % hold off;
end

function animate_trajectory(x, v, L, target_set, delay, filename)
    fig = figure;

    % Target patch
    lb = min(target_set.points(:,1:2),[],1)' - (target_set.eta(1:2) / 2);
    ub = max(target_set.points(:,1:2),[],1)' + (target_set.eta(1:2) / 2);
    p = zeros(3,2);
    p(:,1) = [0, 1.5*L*sin(lb(1)), 1.5*L*sin(ub(1))];
    p(:,2) = [0, -1.5*L*cos(lb(1)), -1.5*L*cos(ub(1))];

    % Animated trajectory
    traj = zeros([2, size(x)]);
    traj(2,:,1) = sin(x(:,1)) * L;
    traj(2,:,2) = -cos(x(:,1)) * L;

    % Driving freq visualization
    cmap = copper(50);
    colors = cmap(int8(v),:);

    % Plot
    hold on;
    % patch(p(:,1), p(:,2), [.6,1,.6], 'EdgeColor',[.6,1,.6])
    pivot = plot(0, 0, 'o', 'MarkerSize',20, 'MarkerFaceColor',colors(1,:), 'MarkerEdgeColor',colors(1,:));
    pendulum = plot(traj(:,1,1), traj(:,1,2), 'LineWidth',5, 'Color','black');
    xlim([-L, L])
    ylim([-L, L])
    colormap(cmap)
    cb3 = colorbar;
    clim([1,50])
    ylabel(cb3, 'Driving Angular Frequency')
    axis off;
    hold off;

    if filename ~= ""
        im = frame2im(getframe(fig));
        [A,map] = rgb2ind(im,256);
        imwrite(A, map, filename, 'gif', 'LoopCount',Inf, 'DelayTime',delay)
    end

    % Animate
    for i = 1:size(x,1)
        pause(delay)
        pivot.MarkerFaceColor = colors(i,:);
        pivot.MarkerEdgeColor = colors(i,:);
        pendulum.XData = traj(:,i,1);
        pendulum.YData = traj(:,i,2);
        drawnow

        if filename ~= ""
            im = frame2im(getframe(fig));
            [A,map] = rgb2ind(im,256);
            imwrite(A, map, filename, 'gif', 'WriteMode','append', 'DelayTime',delay)
        end
    end

    if filename ~= ""
        fprintf("Animation saved to %s\n", filename)
    end
end
