%% Parameters
L=3;            % wheelbase
b=1.5;          % distance from rear wheel to center of mass

dt=0.01;        % time discretization
T=0:dt:6;       % time span

%% load reference trajectory
load('part1_traj_05_timestep.mat')

U_ref=interp1(0:0.05:6,[U,U(:,end)]',T)';

Y_ref=interp1(0:0.05:6,Y,T)';

% the reference trajectory Y_ref is given as a 3x601 double

%% 4.1 Discrete-time A and B matrices
%these are the system linearized in discrete time about the reference 
%trajectory i.e. x(i+1)=A_i*x_i+B_i*u_i

A=@(i) [0,0,-U_ref(1,i)*sin(Y_ref(3,i))-((b*U_ref(1,i)*tan(U_ref(2,i))*cos(Y_ref(3,i)))/L);
        0,0,U_ref(1,i)*cos(Y_ref(3,i))-((b*U_ref(1,i)*tan(U_ref(2,i))*sin(Y_ref(3,i)))/L);
        0,0,0]*0.01 + eye(3);

B=@(i) [cos(Y_ref(3,i))-((b*tan(U_ref(2,i))*sin(Y_ref(3,i)))/L), -(b*U_ref(1,i)*sin(Y_ref(3,i))*sec(U_ref(2,i))^2)/L;
        sin(Y_ref(3,i))+((b*tan(U_ref(2,i))*cos(Y_ref(3,i)))/L), (b*U_ref(1,i)*cos(Y_ref(3,i))*sec(U_ref(2,i))^2)/L;
        tan(U_ref(2,i))/L, (U_ref(1,i)*sec(U_ref(2,i))^2)/L]*0.01;

%% 4.2 Number of decision variables for colocation method
npred=10;

Ndec= 53;

%% 4.3 Write and test function to construct Aeq beq (equality constraints
%enforce x(i+1)=A_i*x_i+B_i*u_i for the prediction horizon) 
%check the values of Aeq and beq timestep 1 (t=0)
%with the given initial condition of Y(0)=[0.25;-0.25;-0.1];

init = [0.25;-0.25;-0.1];
[Aeq_test1, beq_test1]=eq_cons(1,A,B,init);

% 
%% 4.4 Write and test function to generate limits on inputs 
%check the value at timestep 1 (t=0)

[Lb_test1, Ub_test1] = bound_cons(1,U_ref);

%% 4.5 Simulate controller working from initial condition [0.25;-0.25;-0.1]
%use ode45 to between inputs

%Calculate max distance error between the actual and nominal trajectories,
%when the x value of the actual trajectory is between or equal to 3 and 4 m

Q = [1 0 0; 0 1 0; 0 0 0.5];
R = [0.1 0; 0 0.01];

H = blkdiag(Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,Q,R,R,R,R,R,R,R,R,R,R);
f = zeros(53,1);

Zc = init;
Z_store = zeros(length(T),3);
Z_store(1,:) = Zc;

% for i=1:591
%     % calculate inputs using MPC
%     % create matrices
%     [Aeqi, beqi]=eq_cons(i,A,B, Y(:, i)-Y_ref(:, i));
%     [Lbi, Ubi]=bound_cons(i,U_ref,input_range);
% 
%     z=quadprog(H,zeros(53, 1),[],[],Aeqi,beqi,Lbi,Ubi);  %solve the quadratic program
%     u_del=z(34:35, 1)+U_ref(:, i);
%     u=u_del(1);
%     delta_f=u_del(2);
%     [t,y_t]=ode45(@(t,y) q45(t,y,u, delta_f), 0:0.001:0.01, Y(:, i));  % Ask why this is returning size(t)x3 and not 1x3???
%     Y(:, i+1)=y_t(end, :);
% end

for i = 1:length(T)-10
       
    [Aeq,Beq] = eq_cons(i,A,B,Zc-Y_ref(:,i));
    [Lb,Ub] = bound_cons(i,U_ref);
    Zo = quadprog(H,f,[],[],Aeq,Beq,Lb,Ub);
    U = Zo(34:35,1)+U_ref(:,i);
    [t,Zn] = ode45(@(t,Zc) kinematic_bicycle(t,Zc,U,b,L),[0,dt],Zc);
    
    Zc = Zn(end,:)';
    
    Z_store(i+1,:) = Zc';

end

ii = Z_store(:, 1) >= 3 & Z_store(:, 1) <= 4;
Z_ = Z_store(ii, :);
Z_ref_ = Y_ref(:,ii)';
max_dist_error=max(sqrt((Z_(:,1) - Z_ref_(:,1)).^2 + (Z_(:,2) - Z_ref_(:,2)).^2 ));

% max_dist_error=

figure
plot(T,Z_store)

% I highly recommend you write functions for constraint generation and dynamics down here and 
%call them when needed above, for example,

function [Aeq,beq]=eq_cons(idx,A,B,init)
    % build matrix for A_i*x_i+B_i*u_i-x_{i+1}=0
    % in the form Aeq*z=beq
    % initial_idx specifies the time index of initial condition from the reference trajectory 
    % A and B are function handles above

    Aeq = [eye(3), zeros(50,3)';
           A(idx), -eye(3), zeros(27,3)', B(idx), zeros(18,3)';
           zeros(3,3)', A(idx+1), -eye(3), zeros(24,3)', zeros(2,3)', B(idx+1), zeros(16,3)';
           zeros(6,3)', A(idx+2), -eye(3), zeros(21,3)', zeros(4,3)', B(idx+2), zeros(14,3)';
           zeros(9,3)', A(idx+3), -eye(3), zeros(18,3)', zeros(6,3)', B(idx+3), zeros(12,3)';
           zeros(12,3)', A(idx+4), -eye(3), zeros(15,3)', zeros(8,3)', B(idx+4), zeros(10,3)';
           zeros(15,3)', A(idx+5), -eye(3), zeros(12,3)', zeros(10,3)', B(idx+5), zeros(8,3)';
           zeros(18,3)', A(idx+6), -eye(3), zeros(9,3)', zeros(12,3)', B(idx+6), zeros(6,3)';
           zeros(21,3)', A(idx+7), -eye(3), zeros(6,3)', zeros(14,3)', B(idx+7), zeros(4,3)';
           zeros(24,3)', A(idx+8), -eye(3), zeros(3,3)', zeros(16,3)', B(idx+8), zeros(2,3)';
           zeros(27,3)', A(idx+9), -eye(3), zeros(18,3)', B(idx+9)];

    beq = [init;zeros(30,1)];

end

function [Lb,Ub]=bound_cons(idx,U_ref)
% initial_idx is the index along uref the initial condition is at
    
    U = U_ref(:,idx:idx+9);
    
    U = reshape(U,[20,1]);

    Lb = [-Inf*ones(33,1);
          0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5;0;-0.5]-[zeros(33,1);U];

    Ub = [Inf*ones(33,1);
          1;0.5;1;0.5;1;0.5;1;0.5;1;0.5;1;0.5;1;0.5;1;0.5;1;0.5;1;0.5]-[zeros(33,1);U];

end

function ydot = kinematic_bicycle(t,y,U,b,L)

    psi = y(3);
    u = U(1);
    delta = U(2);

    ydot = [u*cos(psi) - ((b*u*tan(delta)*sin(psi))/L);
            u*sin(psi) + ((b*u*tan(delta)*cos(psi))/L);
            u*tan(delta)/L];

end


