clc; clear; close all;

%% Computational Parameters Setting

N = 20;
d = 2.5;

T = 0.02;
R = 0.25;
V = 0.25;

drcut = 4;
nSteps = 1000;

plotDelay = 0.02;
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
    vx0(i) = 1;
    vy0(i) = 1;
    xprev(i) = xcurr(i) - vx0(i)*T;
    yprev(i) = ycurr(i) - vy0(i)*T;
end

% figure;
% p = scatter(xcurr,ycurr,'k','filled'); 
% hold on;
% ax = gca;
% ax.XLim = [0 L];   ax.XLabel.String = 'x';
% ax.YLim = [0 L];   ax.YLabel.String = 'y';
% ax.GridAlpha = 0.15;    ax.Box = 'on';
% ax.Title.String = ({'2D Molecular Dynamics Simulation';'step = 0'});


%% MD Simulation
run = 1;
count = 0;
m = 0;

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
        count = count + 1;
        
%% TAKING THE VALUES OF SPEED AFTER EVERY 10 ITERATIONS
        
        if count == 10
        m = m + 1;
        v(:,m) = vx;
        count = 0;
        end
    
%%        
      for i = 1:N
           xprev(i) = xcurr(i);
           yprev(i) = ycurr(i);
           xcurr(i) = xnew(i);
           ycurr(i) = ynew(i); 
      end
    

%     if mod(n,plotStep) == 0       
%         p.XData = xcurr;
%         p.YData = ycurr;
%         scatter(xcurr(run),ycurr(run),'ko');
%         ax.Title.String = ({'2D Molecular Dynamics Simulation';['step = ',num2str(n)]});
%         pause(plotDelay);
%     end

end

%% SEPERATING INTERVALs 0-20, 20-40, 40-60

for n=1:20
    v_20(:,n) = v(:,n);
    v_40(:,n) = v(:,n+19);
    v_60(:,n) = v(:,n+39);
end

%% MAKING BINS FOR PROBABILITY CALCULATION

s = size(v_20);
n_bins = 16;

v_20 = sum(v_20,2)/20;

interval_20 = (max(v_20)-min(v_20))/n_bins;
probs_20 = zeros(n_bins,1);
j = 1;

for n = 1:n_bins
for i = 1:s(1)
    
   if  v_20(i) <= n*interval_20 && v_20(i)> (n-1)*interval_20
       probs_20(n) = probs_20(n) + 1;
   end
end    
end
probs_20 = probs_20/s(1);



s = size(v_40);

v_40 = sum(v_40,2)/20;

interval_40 = (max(v_40)-min(v_40))/n_bins;
probs_40 = zeros(n_bins,1);
j = 1;

for n = 1:n_bins
for i = 1:s(1)
    
   if  v_40(i) <= n*interval_40 && v_40(i)> (n-1)*interval_40
       probs_40(n) = probs_40(n) + 1;
   end
end    
end
probs_40 = probs_40/s(1);



s = size(v_60);

v_60 = sum(v_60,2)/20;

interval_60 = (max(v_60)-min(v_60))/n_bins;
probs_60 = zeros(n_bins,1);
j = 1;

for n = 1:n_bins
for i = 1:s(1)
    
   if  v_60(i) <= n*interval_60 && v_60(i)> (n-1)*interval_60
       probs_60(n) = probs_60(n) + 1;
   end
end    
end
probs_60 = probs_60/s(1);

%% PLOTTING VELOCITY DISTRIBUTION
up = (max(v_20)-interval_20);
bins = (min(v_20):interval_20:up)';

plot(bins,probs_20,'ro-')
xlabel('vx')
ylabel('P(vx)')
grid on
hold on
plot(bins,probs_40,'bo-')
hold on
plot(bins,probs_60,'mo-')
hold on
legend('T = 0-20','T = 20-40','T = 40-60')

