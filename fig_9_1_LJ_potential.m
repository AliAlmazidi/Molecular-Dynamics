clc; clear; close all;

epsilon = 1;
sigma = 3.4;
r = 0:0.1:4*sigma;
V = 4*epsilon * ((sigma./r).^12 - (sigma./r).^6);
plot(r/sigma,V/epsilon,'k','linewidth',3);
ylim([-1.5 2]);
xlabel('r / \sigma'); ylabel('V(r) / \epsilon'); title('Lennard-Jones potential');