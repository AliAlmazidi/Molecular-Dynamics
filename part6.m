clc
clear all
close all

%generating random nos
N(1) = 100;
x = rand(1,N);
k(1) = 1;

check(1) = sqrt(N)*abs(1/N * sum(x.^k) - 1/(k+1));
N(2) = 100;
x = rand(1,N(2));
k(2) = 3;

check(2) = sqrt(N(2))*abs(1/N(2) * sum(x.^k(2)) - 1/(k(2)+1));
N(3) = 100;
x = rand(1,N(3));
k(3) = 7;

check(3) = sqrt(N(3))*abs(1/N(3) * sum(x.^k(3)) - 1/(k(3)+1));
N(4) = 10000;
x = rand(1,N(4));
k(4) = 1;

check(4) = sqrt(N(4))*abs(1/N(4) * sum(x.^k(4)) - 1/(k(4)+1));
N(5) = 10000;
x = rand(1,N(5));
k(5) = 3;
check(5) = sqrt(N(5))*abs(1/N(5) * sum(x.^k(5)) - 1/(k(5)+1));

N(6) = 10000;
x = rand(1,N(6));
k(6) = 7;
check(6) = sqrt(N(6))*abs(1/N(6) * sum(x.^k(6)) - 1/(k(6)+1));

N(7) = 100000;
x = rand(1,N(7));
k(7) = 1;
check(7) = sqrt(N(7))*abs(1/N(7) * sum(x.^k(7)) - 1/(k(7)+1));

N(8) = 100000;
x = rand(1,N(8));
k(8) = 3;
check(8) = sqrt(N(8))*abs(1/N(8) * sum(x.^k(8)) - 1/(k(8)+1));

N(9) = 100000;
x = rand(1,N(9));
k(9) = 7;
check(9) = sqrt(N(9))*abs(1/N(9) * sum(x.^k(9)) - 1/(k(9)+1));


N = N';
k = k';
check = check';
disp(table(N,k,check));
