clear variables
clc

dt=0.005;  % Sampling time in seconds
t_total=5;  % Simultion time in seconds
t=0:dt:t_total;

% Initializing simulation parameters
u_0 = [0, 0 ,0]';  % Initial Guess for rod's curvature at its base
P=eye(3,3);
T=[];
U=[];
r_end=[];   % Cartesian coordinates of the rod backbone
u_end=[];   % Instantaneous curvature of the rod's tip
for i=1:length(t)

V=0.01;    
q=[V*t(i) 2*pi*t(i)/10];  % rods insertion and rotation
f=[1*sin(2*pi*t(i)/10) 1*cos(2*pi*t(i)/10) 1*sin(2*pi*t(i)/10)]'; % External Force 

% solving rod's model (equations 8 and 9), inputs are: rods insertion and rotation,
% initial curvature, external forces, respectively.
[r,u, B] =moving_rod(q, u_0, f);  

if rem(i,50)==1
    % plotting rod's backbone
    p1=plot3(r(:,1), r(:,2), r(:,3), 'k', 'LineWidth', 1);
    hold on
end


A=V*eye(3,3);  % gamma (equation 3a)
Q=30*eye(3,3);
R=120*eye(3,3);  % inverse of R

% Observer (equations 11 and 12)
dP=-A.' *P -P*A- P*B*R*B.'*P+Q;
du0=-R*B.'*P*u;
P= P + dt.* dP;
u_0= u_0 + dt.* du0;

% recording rod's tip position
r_end=[r_end;r(end,:)];  % rod's tip position

end


% plot rod's tip position
p2 = plot3(r_end(:,1),r_end(:,2),r_end(:,3),'r', 'LineWidth', 3);
legend([p1 p2],'Rod Backbone','Rod Tip Trajectory','Location','NorthWest')
axis equal
grid on


%% main ode solver
function [r, u, dU0] = moving_rod(q, u_0, f)


% length of the rod
l=1e-3.*400;
% physical parameters of the rod (Table 1)
E=70e9;
r1 =(1.5)*1e-3;
r2=(1)*1e-3;
J= pi/2*(r1^4-r2^4);
I=  0.5 * J;
G=10e9; 
% pre-curvature of the rod
Ux=14;
Uy=5;
Uz=0;

q_0=[-0.2 0]';
B=q(1)+q_0(1);  % length of the rod before template
%initial angles
alpha=q(2)+q_0(2);



%% Solving ode for shape of the rod
r=[];
r0=[0 0 0]'; R0=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
R0=reshape(R0',[9,1]);
du0=reshape(eye(3,3)',[9,1]);    
s_span =[0 l+B]; 

% initial conditions
y0=zeros(15+9+9,1);
y0(1:3,1)=u_0;
y0(4:15,1)=[r0 ; R0];
y0(16:24,1)=du0;

[~,y] = ode23(@(s,y) ode(s,y,Ux,Uy,Uz,E.*I',G.*J, f,u_0), s_span, y0);
shape=[y(:,4),y(:,5),y(:,6)];
r=[r; shape];
u=y(end,1:3)';
dU0=reshape(y(end,16:24),[3,3])';
end


%% ODE (equations 8 and 9)
function dydt = ode(~,y,Ux,Uy,Uz,EI,GJ, f,u_0)

dydt=zeros(15+9+9,1);
% first 3 elements of y are curvatures along x, y, and z
% next 12 elements are r (position) and R (orientations), respectively
% next 9 elements are partial derivative of u with respect to u0
% next 9 are are partial derivative of R^T*f with respect to u0

% calculating 1st tube's curvatures in x and y direction
ux = y(1);
uy = y(2);
uz = y(3); 
U=[ux uy uz].';

GJ(EI==0)=0;
K =diag([EI EI GJ]);
u_hat = [0 -uz uy;uz 0 -ux;-uy ux 0];
e3=[0 0 1]';
e3_hat = [0 -1 0;1 0 0;0 0 0];
Us=[Ux Uy Uz].';

R=[y(7) y(8) y(9);y(10) y(11) y(12);y(13) y(14) y(15)];
B=[y(16) y(17) y(18);y(19) y(20) y(21);y(22) y(23) y(24)]; % partial derivative of u with respect to u0
C=[y(25) y(26) y(27);y(28) y(29) y(30);y(31) y(32) y(33)]; % partial derivative of R'f with respect to u0

% odes
% partial derivative of R'f with respect to u0
D=R.'*f;
C_hat=[0 -D(3) D(2);D(3) 0 -D(1);-D(2) D(1) 0];
dC= C_hat*B-u_hat*C;

% partial derivative of u with respect to u0
A=K*(U-Us);
A_hat=[0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
dB=-inv(K)*(-A_hat*B+u_hat*K*B+ e3_hat*C );

% curvature
du= B*u_0;

du = -inv(K)* ((u_hat) * K * (U-Us)+ e3_hat*R.'*f) ;
% R and r
dr = R*e3;
dR=R*(u_hat);


dydt(1:3,1)=du;
dydt(4:6)=dr(1:3);
dydt(7:15)=reshape(dR.',[9,1]);
dydt(16:24)=reshape(dB.',[9,1]);
dydt(25:33)=reshape(dC.',[9,1]);

end
