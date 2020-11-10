clc
clear;
clc;
clear;

tic
h=1/400;  x0=-2; x1=4;  xx=x0:h:x1;        m=round((x1-x0)/h);
dt=1/200; t0=0; t1=10;  tt=t0:dt:t1; n=round((t1-t0-2*dt)/dt);

ci=0.2; %noise intensity

x00=0.2402; %initial point of stochstic THC model
xT=1.0687;  %final point of stochstic THC model



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Forward Fokker-Planck equation 
P1=zeros(n+1,m+1);
F_ba=1.1;
nu2=6.2;

f=@(x) (F_ba-x*(1+nu2*(1-x).^2)); %drift term

P1(1,:)=(pi*2*dt*ci^2)^(-1/2)*exp(-(xx-x00-f(x00)*dt).^2/(2*dt*ci^2)); % initial condition delta function

r1=dt/(4*h); r2=(dt*ci^2)/(4*h^2); 

F=zeros(1,m+1);
  for i=1:m+1
      F(1,i)=f(xx(i)); 
  end
  
 a11=(1+2*r2)*ones(1,m-1);
 a12=-r2*ones(1,m-2);
 A1=diag(a11,0)+diag(a12,-1)+diag(a12,1); 
 
 b11=r1*F(3:end-1);
 b12=-r1*F(2:end-2);
 B1=diag(b11,1)+diag(b12,-1);% B is also the coefficient matrix of P(k,2:end-1)

 c11=(1-2*r2)*ones(1,m-1);
 c12=r2*ones(1,m-2);
 C1=diag(c11,0)+diag(c12,-1)+diag(c12,1); 
 

 E1=zeros(n+1,m-1);
  for k=1:n
    E1(k,:)=[r1*F(1,1)*P1(k+1,1)+r2*P1(k+1,1)+r2*P1(k,1),zeros(1,m-3),-r1*F(1,end)*P1(k+1,m+1)+r2*P1(k+1,m+1)+r2*P1(k,m+1)];
    P1(k+1,2:end-1)=(inv(A1+B1)*((C1-B1)*P1(k,2:end-1)'+E1(k,:)'))';   
   end
  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %backward fokker-Planck equation
  
  P2=zeros(n+1,m+1);
  P2(n+1,:)=(2*pi*dt*ci^2)^(-1/2)*exp(-(xx-xT+f(xT)*dt).^2/(2*dt*ci^2));
  
  a21=(-1-2*r2)*ones(1,m-1);
  a22=r2*ones(1,m-2);
  A2=diag(a21,0)+diag(a22,-1)+diag(a22,1);
 
  b21=r1*F(1,2:end-2);
  b22=-r1*F(1,3:end-1);
  B2=diag(b21,1)+diag(b22,-1);

  c21=(-1+2*r2)*ones(1,m-1);
  c22=-r2*ones(1,m-2);
  C2=diag(c21,0)+diag(c22,-1)+diag(c22,1);

  D2=zeros(n+1,m-1);
  E2=zeros(n+1,m-1);

 for k=n:-1:1
     E2(k+1,:)=[-r2*(P2(k,1)+P2(k+1,1))+r1*F(1,2)*P2(k,1),zeros(1,m-3),-r1*F(1,m)*P2(k,m+1)-r2*(P2(k,m+1)+P2(k+1,m+1))];
     P2(k,2:end-1)=(inv(A2+B2)*((C2+B2)*P2(k+1,2:end-1)'+E2(k+1,:)'))';
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  PP=zeros(n+1,m+1);

 for i=1:n+1
    PP(i,:)=P1(i,:).*P2(i,:);
    [m,index]=max(PP(i,:));
    x_max(i)=xx(index);
 end
 

tt=tt(2:end-1);
plot(tt,x_max,'.-')
hold on 

 toc
% % %  
 
 
