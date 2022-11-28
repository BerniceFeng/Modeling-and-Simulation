clear
clc
close all

g=9.8; % gravity
num=2; % number of unit block
dt=0.001/num; % time step as a function of of number of unit block
a=6*num; % number of the deck block
b=3*num; % number of the tower block
m=0.8/num; % mass of each node

% stiffness and dampness of the deck
S_hori1=3000*num;
S_vert1=4000*num;
S_diag1=3000*num;
D_hori1=10*num;
D_vert1=10*num;
D_diag1=10*num;

% stiffness and dampness of the tower
S_hori2=3000*num;
S_vert2=7000*num;
S_diag2=3000*num;
D_hori2=10*num;
D_vert2=10*num;
D_diag2=10*num;

% stiffness and dampness of the main cable
S_tilt=4000*num;
S_curv=4000*num;
D_tilt=10*num;
D_curv=10*num;

% stiffness and dampness of the suspension cable
S_down1=4000*num;
S_down2=3000*num;
D_down1=10*num;
D_down2=10*num;

k_max=(2*a+2)+(4*b+4)+(6*num-5); % index of the node
l_max=(5*a+1)+(10*b+2)+(6*num-2)+(6*num-5); % index of the link

%rest length of the springs
R0hori1=1; 
R0vert1=1.2;
R0diag1=1.4;
R0hori2=1; 
R0vert2=1;
R0diag2=1.4;
R0tilt=1.3;
R0curv=1;
R0down1M=zeros(2*num-2,1);
R0down11M=zeros(2*num-2,1);
R0down2M=zeros(2*num-1,1);

X=zeros(k_max,2); %position setup
U=zeros(k_max,2); %velocity setup
jj=zeros(l_max,1); %node on one side of the spring
kk=zeros(l_max,1); %node on the other side of the spring

m_car = 50; %car mass
v_car = 0.1; %car velocity
F_car = -m_car*g; %force given by the car
t_car_start = 5; %starting time of the car

t_eq_start = 10; %starting time of the earthquake
t_eq_end = 15; %ending time of the earthquake

%% construction of the node system X

%construction of the deck
for k = 1:2*a+2 
    if mod(k,2) == 0
        X(k,1) = (k/2-1)*R0hori1;
        X(k,2) = 0;
    else
        X(k,1) = (k-1)*R0hori1/2;
        X(k,2) = R0vert1;
    end
end

%construction of the tower1
for k = 1:2*b+2
    if mod(k,2) == 0
        X(k+2*a+2,1) = (2*num-R0hori1)*R0hori2;
        X(k+2*a+2,2) = -num-1+k/2;
    else
        X(k+2*a+2,1) = (2*num)*R0hori2;
        X(k+2*a+2,2) = -num-1+(k+1)/2;
    end
end         

%construction of the tower2
for k = 1:2*b+2
    if mod(k,2) == 0
        X(k+2*a+2+2*b+2,1) = (4*num+R0hori1)*R0hori2;
        X(k+2*a+2+2*b+2,2) = -num-1+k/2;
    else
        X(k+2*a+2+2*b+2,1) = (4*num)*R0hori2;
        X(k+2*a+2+2*b+2,2) = -num-1+(k+1)/2;
    end
end

%construction of the cables
for k = 1:2*num-2
    X(k+2*a+2+4*b+4,1) = k;
    X(k+2*a+2+4*b+4,2) = k+1;
end  

for k = 1:2*num-2
    X(k+2*a+2+4*b+4+2*num-2,1) = a-k;
    X(k+2*a+2+4*b+4+2*num-2,2) = k+1;
end 
    
for k = 1:2*num-1
    X(k+2*a+2+4*b+4+4*num-4,1) = k+2*num;
    X(k+2*a+2+4*b+4+4*num-4,2) = num*2;
end 

M = ones(k_max,1).*m; % mass matrix of the nodes


%% construction of the spring system 

% deck

%horizontal
for blk = 1:a
    % This constructs the uppper part
    jj(blk) = 2*blk-1;
    kk(blk) = 2*blk+1;
    %This constructs the lower part
    jj(a+blk) = 2*blk;
    kk(a+blk) = 2*(blk+1);
end

%verticle
for blk = 1:a+1
    jj(2*a+blk)  = 2*blk-1;
    kk(2*a+blk)  = 2*blk;
end


%diagonal
for blk = 1:a
    jj(3*a+1+blk) =  2*blk-1;
    kk(3*a+1+blk) =  2*(blk+1);
end

for blk = 1:a
    jj(4*a+1+blk) =  2*blk;
    kk(4*a+1+blk) =  2*blk+1;
end


%construct the towers:

% tower NO1

%vertical left
for blk = 1:b
    jj(5*a+1+blk) = 2*a+2+2*blk-1;
    kk(5*a+1+blk)  = 2*a+2+2*blk+1;
end

%vertical right
for blk = 1:b
    jj(5*a+1+b+blk) = 2*a+2+2*blk;
    kk(5*a+1+b+blk)  = 2*a+2+2*(blk+1);
end

%horizontal

for blk = 1:b+1
    jj(5*a+1+2*b+blk) = 2*a+2+2*blk-1;
    kk(5*a+1+2*b+blk) = 2*a+2+2*blk;
end
    
%diagonal 1
for blk = 1:b
    jj(5*a+1+3*b+1+blk) =  2*a+2+2*blk-1;
    kk(5*a+1+3*b+1+blk) =  2*a+2+2*(blk+1);
end

%diagonal 2
for blk = 1:b
    jj(5*a+1+4*b+1+blk) =  2*a+2+2*blk;
    kk(5*a+1+4*b+1+blk) =  2*a+2+2*blk+1;
end


% tower NO2

%vertical left
for blk = 1:b
    jj(5*(a+b)+2+blk) = 2*(a+b)+4+2*blk-1;
    kk(5*(a+b)+2+blk)  = 2*(a+b)+4+2*blk+1;
end

%vertical right
for blk = 1:b
    jj(5*(a+b)+2+b+blk) = 2*(a+b)+4+2*blk;
    kk(5*(a+b)+2+b+blk)  = 2*(a+b)+4+2*(blk+1);
end

%horizontal
for blk = 1:b+1
    jj(5*(a+b)+2+2*b+blk) = 2*(a+b)+4+2*blk-1;
    kk(5*(a+b)+2+2*b+blk) = 2*(a+b)+4+2*blk;
end
    
%diagonal 1
for blk = 1:b
    jj(5*(a+b)+2+3*b+1+blk) =  2*(a+b)+4+2*blk-1;
    kk(5*(a+b)+2+3*b+1+blk) =  2*(a+b)+4+2*(blk+1);
end

%diagonal 2
for blk = 1:b
    jj(5*(a+b)+2+4*b+1+blk) =  2*(a+b)+4+2*blk;
    kk(5*(a+b)+2+4*b+1+blk) =  2*(a+b)+4+2*blk+1;
end


% main cable

%starting link
jj(5*(a+b)+2+5*b+2) = 1;
kk(5*(a+b)+2+5*b+2) = 1+2*a+2+4*b+4;

%middle part
for blk = 1:2*num-3
    jj(5*(a+b)+2+5*b+2+blk) = blk+2*a+2+4*b+4;
    kk(5*(a+b)+2+5*b+2+blk) = blk+1+2*a+2+4*b+4;
    jj(5*(a+b)+2+5*b+2+4*num-2+2*num-1+blk) = blk+2*a+2+4*b+4;
    kk(5*(a+b)+2+5*b+2+4*num-2+2*num-1+blk) = 2*blk+1;
end

%ending link
jj(5*(a+b)+2+5*b+2+2*num-2) = 2*num-3+1+2*a+2+4*b+4;
kk(5*(a+b)+2+5*b+2+2*num-2) = 2*a+2*b+4;

jj(5*(a+b)+2+5*b+2+4*num-2+2*num+2*num-3) = 2*num-3+1+2*a+2+4*b+4;
kk(5*(a+b)+2+5*b+2+4*num-2+2*num+2*num-3) = 2*(2*num-3)+3;

%starting link
jj(5*(a+b)+2+5*b+2+2*num-1) = 2*a+1;
kk(5*(a+b)+2+5*b+2+2*num-1) = 2*num-3+1+2*a+2+4*b+4+1;

%middle part
for blk = 1:2*num-3
    jj(5*(a+b)+2+5*b+2+2*num-1+blk) = blk+2*num-3+1+2*a+2+4*b+4;
    kk(5*(a+b)+2+5*b+2+2*num-1+blk) = blk+1+2*num-3+1+2*a+2+4*b+4;
    jj(5*(a+b)+2+5*b+2+4*num-2+2*num+2*num-3+blk) = ...
        blk+2*num-3+1+2*a+2+4*b+4;
    kk(5*(a+b)+2+5*b+2+4*num-2+2*num+2*num-3+blk) = 2*a-2*blk+1;
end

%ending link
jj(5*(a+b)+2+5*b+2+4*num-3) = 2*num-3+1+2*num-3+1+2*a+2+4*b+4;
kk(5*(a+b)+2+5*b+2+4*num-3) = 2*a+4*b+6;

jj(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+1) = ...
    (2*num-3)+1+2*num-3+1+2*a+2+4*b+4;
kk(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+1) = 2*a-2*(2*num-3)-1;

% suspension cable

%starting link
jj(5*(a+b)+2+5*b+2+4*num-2) = 2*a+2*b+3;
kk(5*(a+b)+2+5*b+2+4*num-2) = 2*num-3+1+2*num-3+1+2*a+2+4*b+4+1;

%middle part
for blk = 1:2*num-2
    jj(5*(a+b)+2+5*b+2+blk+4*num-2) = ...
        blk-1+2*num-3+1+2*num-3+1+2*a+2+4*b+4+1;
    kk(5*(a+b)+2+5*b+2+blk+4*num-2) = ...
        blk+2*num-3+1+2*num-3+1+2*a+2+4*b+4+1;
    jj(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+1+blk) = ...
        blk-1+2*num-3+1+2*num-3+1+2*a+2+4*b+4+1;
    kk(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+1+blk) = ...
        2*(2*num)+2*blk+1;
end

%ending link
jj(5*(a+b)+2+5*b+2+4*num-2+2*num-1) = ...
    2*num-2+2*num-3+1+2*num-3+1+2*a+2+4*b+4+1;
kk(5*(a+b)+2+5*b+2+4*num-2+2*num-1) = 2*a+4*b+5;

jj(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+2*num) = ...
    2*num-2-1+2*num-3+1+2*num-3+1+2*a+2+4*b+4+2;
kk(5*(a+b)+2+5*b+2+4*num-2+2*num+(2*num-3)*2+2*num) = ...
    2*(2*num)+2*(2*num-2)+3;

%rest lengths of the suspension cables
for i=1:2*num-2
    R0down1M(i)=0.9*i;
end

for i=1:2*num-2
    R0down11M(i)=0.9*i;
end

for i=1:2*num-1
    R0down2M(i)=0.25*(abs(a/2-(i+2*num)))*num/2+num;
end 

%% construction of stiffness, dampness and rest length
S=    [S_hori1*ones(2*a,1);S_vert1*ones(a+1,1);S_diag1*ones(2*a,1);...
    S_hori2*ones(2*b,1);S_vert2*ones(b+1,1);S_diag2*ones(2*b,1);...
    S_hori2*ones(2*b,1);S_vert2*ones(b+1,1);S_diag2*ones(2*b,1);...
    S_tilt*ones(4*num-2,1);S_curv*ones(2*num,1);...
    S_down1*ones(4*num-4,1);S_down2*ones(2*num-1,1)];
D=    [D_hori1*ones(2*a,1);D_vert1*ones(a+1,1);D_diag1*ones(2*a,1);...
    D_hori2*ones(2*b,1);D_vert2*ones(b+1,1);D_diag2*ones(2*b,1);...
    D_hori2*ones(2*b,1);D_vert2*ones(b+1,1);D_diag2*ones(2*b,1);...
    D_tilt*ones(4*num-2,1);D_curv*ones(2*num,1);...
    D_down1*ones(4*num-4,1);D_down2*ones(2*num-1,1)];
Rzero=[R0hori1*ones(2*a,1);R0vert1*ones(a+1,1);R0diag1*ones(2*a,1);...
    R0hori2*ones(2*b,1);R0vert2*ones(b+1,1);R0diag2*ones(2*b,1);...
    R0hori2*ones(2*b,1);R0vert2*ones(b+1,1);R0diag2*ones(2*b,1);...
    R0tilt*ones(4*num-2,1);R0curv*ones(2*num,1);...
    R0down1M;R0down11M;R0down2M]; 


%% draw the graph
myVideo = VideoWriter('myVideoFile'); %open video file
myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
open(myVideo)

pbaspect([a 1 1])
p=plot([X(jj(:),1),X(kk(:),1)]',[X(jj(:),2),X(kk(:),2)]');
hold on
axis equal
axis manual
axis([-1,a+1,-num-2,2*num+2])
drawnow
pause(0.01) %Pause and grab frame
frame = getframe(gcf); %get frame
writeVideo(myVideo, frame);


%% loop for simulation

for clock=1:1000
  t = clock*dt;
  DX = X(jj,:) - X(kk,:);
  DU = U(jj,:) - U(kk,:); 
  R = sqrt(sum(DX.^2,2));
  sum(DX.*DU,2)
  T = S.*(R-Rzero) + (D./R).*sum(DX.*DU,2);
  TR=T./R; 
  FF=[TR,TR].*DX; 

  F=zeros(k_max,2); 
  F(:,2) = - M*g;
  for link=1:l_max
    F(kk(link),:)=F(kk(link),:)+FF(link,:);
    F(jj(link),:)=F(jj(link),:)-FF(link,:);
  end

  
% when the car passes

  if(t_car_start == t)
      cp=plot(X(1,1), X(1,2),'o');
      tt=0;
  end

  if((t_car_start < t) && ((t < (2*a+2-1)/2*dt*100*num + t_car_start)))
      if mod(t,0.1) == 0
        delete(cp);
        tt=round((t-t_car_start)/dt)/100/(num);
        F((2*tt)+1,:) = F((2*tt)+1,:) + F_car;
        cp=plot(X(2*tt+1,1), X(2*tt+1,2),'o');
      else
        F((2*tt)+1,:) = F((2*tt)+1,:) + F_car;
      end
  end
  
  if t == (2*a+2-1)/2*dt*100*num + t_car_start
      delete(cp);
  end

%   if(t_car_start == t)
%       cp=plot(X(1,1), X(1,2),'o');
%       tt=0;
%   end
% 
%   if((t_car_start < t) && ((t < a/v_car + t_car_start)))
%       if mod(t,0.1) == 0
%         delete(cp);
%         x_car = v_car*(t-t_car_start);
%         for i=1:a
%             if abs(x_car-X(i,1)) < 0.5
%                 F((i)+1,:) = F((i)+1,:) + F_car;
%             end
%         end
%         cp=plot(X(i,1), X(i,2),'o');
%       else
%         F((2*tt)+1,:) = F((2*tt)+1,:) + F_car;
%       end
%   end
%   
%   if t == (2*a+2-1)/2*dt*100*num + t_car_start
%       delete(cp);
%   end
      
  U = U + dt*F./[M,M];  
  U(1,:)=0; 
  U(2,:)=0; 
  U(2*a+1,:)=0; 
  U(2*a+2,:)=0; 
  U(2*a+3,:)=0; 
  U(2*a+4,:)=0;
  U(2*(a+b)+5,:)=0;
  U(2*(a+b)+6,:)=0;
  X = X + dt*U; 

  if(t_eq_start < t) && (t < t_eq_end)
      X(2*a+3,1)=X(2*a+3,1)+sin(2*pi*(t-t_eq_start)*10)/400;
      X(2*a+4,1)=X(2*a+4,1)+sin(2*pi*(t-t_eq_start)*10)/400;
      X(2*a+2*b+5,1)=X(2*a+2*b+5,1)+sin(2*pi*(t-t_eq_start)*4)/400;
      X(2*a+2*b+6,1)=X(2*a+2*b+6,1)+sin(2*pi*(t-t_eq_start)*4)/400;
      U(2*a+3,1)=2*pi*cos(2*pi*(t-t_eq_start)*10)/400;
      U(2*a+4,1)=2*pi*cos(2*pi*(t-t_eq_start)*10)/400;
      U(2*a+2*b+5,1)=2*pi*cos(2*pi*(t-t_eq_start)*4)/400;
      U(2*a+2*b+6,1)=2*pi*cos(2*pi*(t-t_eq_start)*4)/400;
  end

  if(mod(clock,100)==0)
    c=0;
    for l=1:l_max
        c=c+1;
        p(c).XData=[X(jj(l),1),X(kk(l),1)];
        p(c).YData=[X(jj(l),2),X(kk(l),2)];
    end
    drawnow
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    
  end

end

close(myVideo)
