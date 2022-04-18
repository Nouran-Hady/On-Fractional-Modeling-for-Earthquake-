 format compact
 k_el = 1.7; 
 kon = 100; 
 koff = 2.5; 

clc;clear;close all

polys=[1 1 1];
orders=[2 1.5 0];
yx=@(x)3-x;

ydash_orders=[0 1];
ydash_values=[1 1];

L=length(orders);

m=max(orders);
D_1=D1(m+1);

z=find(ceil(orders)~=orders);
alpha=orders(z);
D2=Dalpha(m,alpha);

C=sym('c',[1 m+1]);

for n=2:m
    Dn=D_1;
    for j=2:n
        Dn=Dn*D_1;
    end
    D{n}=Dn;
end
D{1}=D_1;

G=G_fcn(yx,m);
x0=0;
P=pi_fcn(m);
Rm=sym(-G'*P);
syms x


for i=1:L
    o=orders(i);
    if o>0
        if ceil(o)~=o
            Rm=Rm+C*D2*P;
        else
            Rm=Rm+C*D{o}*P;
        end
    else
        Rm=Rm+C*P;
    end
end
eq1=int(Rm,x,0,1);

eqs=eq1;
p=subs(P,x,0);
for j=1:length(ydash_orders)
    ydo=ydash_orders(j);
    ydv=ydash_values(j);
    if ydo==0
        eq=sum(C*p)-ydv;
    else
        eq=C*D{ydo}*p-ydv;
    end
    eqs=[eqs;eq];
end

Cs=solve(eqs);
fn=fieldnames(Cs);
for n=1:m+1
    C_sol(1,n)=Cs.(fn{n});
    xlabel('time, s')
         ylabel('y(t)')
         title('graph of fractinal function')
        plot(C_sol); legend;
end
Yx=C_sol*P;

%% Showing Results
fx=sym('0');
syms D y(x)
for j=1:length(orders)
    fx=fx-polys(j)*(D^orders(j))*y(x);
end
fx=fx==yx(x);

disp('For equation :')
disp('')
pretty(fx)
disp('')
for j=1:length(ydash_orders)
    ydo=ydash_orders(j);
    if ydo>0
        zz=D^ydo*y(x0)==ydash_values(j);
         xlabel('time, s')
         ylabel('y(t)')
         title('graph of fractinal function')
         plot(zz); legend;
    else
        zz=y(x0)==ydash_values(j);
          xlabel('time, s')
         ylabel('y(t)')
         title('graph of fractinal function')
        plot(zz); legend;
    end
    pretty(zz)
end
disp('Exact Solution :')
disp('')
pretty(Yx)
disp('Where D is the differential Opeational')


% x = b+h;
x = [0:5:150];
p= [150:5:450];
k = [300:5:450];
b =[200:5:650];
c =[300:5:1000];

y = 0.260563308888 *x.^2 + x + 1;
f = 0.260563308888 *p.^2 + p + 1;
z= 0.260563308888 *k.^2 + k + 1;
g = 0.260563308888 *b.^2 + b + 1;
s = 0.260563308888 *c.^2 + c + 1;
disp(y);
disp(f);
disp(z);
disp(g);
disp(s);


% c = x;

% plot(x, y)
% plot(p, f)
% plot(k, z)
% plot(b, g)
% plot(c, s)

plot(x, y, p, f, k, z, b, g, c, s)

%  l = 2251799813685248;
%  o = x^2;
%  m = l*o;
%  n = 4321022448120783;
%  b = m/n;
%  h = x+1;
%  x = b+h;
%  x = [0:5:100];
% y = 0.2605633088835905028 * x.^2+x+1;;
% % c = x;
% plot( x, y)
% plot(x, y)

xlabel('time, s')
ylabel('c(t)')
title('graph of fractinal function')


% format compact
% k_el = 1.7; 
% kon = 100; 
% koff = 2.5; 
% % function definition
% myfun = @(t,y) [-kon*y(1)*y(2) + koff*y(3);
%     -kon*y(1)*y(2) + koff*y(3) - k_el*y(2);
%     kon*y(1)*y(2) - koff*y(3)];
% % using ODE
% tspan = [0, 50];
% tic; sol = ode23t(myfun,tspan,[4;2.5;0]); toc;
% figure(1); plot(sol.x,sol.y); legend; grid;
% % using fde12
% t0=0; tfinal=50; 
% h=2^(-6); alpha=1.0;
% tic; [t, qsol]=fde12(alpha,myfun,t0,tfinal,[4;2.5;0],h); toc;
% figure(2); plot(t,qsol(1,:)); legend; grid;
 
 
 







