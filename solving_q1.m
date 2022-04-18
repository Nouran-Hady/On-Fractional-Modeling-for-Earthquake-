clc;clear;close all

polys=[1 1 1];
orders=[2 1.5 0];
yx=@(x)1+x;

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
end
Yx=C_sol*P;

%% Showing Results
fx=sym('0');
syms D y(x)
for j=1:length(orders)
    fx=fx+polys(j)*(D^orders(j))*y(x);
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
    else
        zz=y(x0)==ydash_values(j);
    end
    pretty(zz)
end
disp('Exact Solution :')
disp('')
pretty(Yx)
disp('Where D is the differential Opeational')