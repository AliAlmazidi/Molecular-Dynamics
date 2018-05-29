%% MD 2D; DT021A FYP; Ali Al-Mazidi
%  Version 4.1
%  23/04/2018

clc; clear; close all;

%% Computational Parameters Setting

N = 25;
d = 3;

T = 0.02;
R = 0.1;
V = 1;

drcut = 4;
nSteps = 500;

plotDelay = 0.05;
plotStep = 5;

%% MD Initialising 

N = round(N^(1/2))^2;
L = ceil(d*N^(1/2));
drcut2 = drcut*drcut;

xprev = zeros(1,N);     yprev = zeros(1,N);  
xcurr = zeros(1,N);     ycurr = zeros(1,N);
xnew = zeros(1,N);      ynew = zeros(1,N);

vx0 = zeros(1,N);       vy0 = zeros(1,N);
vx = zeros(1,N);        vy = zeros(1,N);

g = d/2:d:L-(d/2);
[gx,gy] = meshgrid(g,g);
gx = gx(:);
gy = gy(:);
gx = gx(1:N);
gy = gy(1:N);

for i = 1:N
    xcurr(i) = gx(i) + 2*(rand-0.5)*R;
    ycurr(i) = gy(i) + 2*(rand-0.5)*R;
    vx0(i) = 2*(rand-0.5)*V;
    vy0(i) = 2*(rand-0.5)*V;
    xprev(i) = xcurr(i) - vx0(i)*T;
    yprev(i) = ycurr(i) - vy0(i)*T;
end

figure;
p = scatter(xcurr,ycurr,'k','filled'); hold on;
ax = gca;   
ax.XLim = [0 L];   ax.XLabel.String = 'x';
ax.YLim = [0 L];   ax.YLabel.String = 'y';
ax.GridAlpha = 0.15;    ax.Box = 'on';
ax.Title.String = ({'2D Molecular Dynamics Simulation';'step = 0'});
disp('Press any key to start the simulation');
pause;

%% MD Simulation

for n = 1:nSteps
    
    dx = repmat(xcurr',[1,N])-repmat(xcurr,[N,1]);
    dy = repmat(ycurr',[1,N])-repmat(ycurr,[N,1]);

    dx = dx - L*round(dx/L);
    dy = dy - L*round(dy/L);

    dr2 = dx.^2 + dy.^2;

    [row,col,DR2] = find(triu(dr2,1));
    linearIndex = (col-1)*N+row;
    dx = dx(linearIndex);
    dy = dy(linearIndex);

    invDR2 = zeros(size(DR2));
    cfIndex = (DR2<drcut2);
    invDR2(cfIndex) = 1./DR2(cfIndex);
    f = 48*invDR2.^4.*(invDR2.^3 - 0.5);

    fx = f.*dx;
    fy = f.*dy;

    fx = full(sparse(row,col,fx,N,N));
    fy = full(sparse(row,col,fy,N,N));

    fx = sum(-fx'+fx,2);
    fy = sum(-fy'+fy,2);

    for i = 1:N
        xnew(i) = 2*xcurr(i) - xprev(i) + fx(i)*T^2;
        ynew(i) = 2*ycurr(i) - yprev(i) + fy(i)*T^2;
        
        if xnew(i) < 0
            xnew(i) = xnew(i) + L;
        end
        if xnew(i) > L
            xnew(i) = xnew(i) - L;
        end
        if ynew(i) < 0
            ynew(i) = ynew(i) + L;
        end
        if ynew(i) > L
            ynew(i) = ynew(i) - L;
        end
        
        vx(i) = (xnew(i) - xprev(i))/(2*T);
        vy(i) = (ynew(i) - yprev(i))/(2*T);
    end
    
    for i = 1:N
           xprev(i) = xcurr(i);
           yprev(i) = ycurr(i);
           xcurr(i) = xnew(i);
           ycurr(i) = ynew(i); 
            
    end
    
    if mod(n,plotStep) == 0       
        p.XData = xcurr;
        p.YData = ycurr;
        scatter(xcurr,ycurr,'k.');
        ax.Title.String = ({'2D Molecular Dynamics Simulation';['step = ',num2str(n)]});
        pause(plotDelay);
    end
    
end
 
%video = VideoWriter('Ali.avi','Uncompressed AVI');
%open (video)
%WriteVideo (video,f);
%close (video)
