close all;
clear;

seq = 4.1212;
par = [seq 1.47 0.93 0.52 0.61 0.62];
dis = [seq 3.13 2.63];
hyb2 = [dis(2) 2.38 2.01 1.87 1.88 1.93];
hyb4 = [dis(3) 2.74 2.57 2.21 2.1 2.02];

speedup_par = zeros(1, length(par));
speedup_dis = zeros(1, length(dis));
speedup_hyb2 = zeros(1, length(hyb2));
speedup_hyb4 = zeros(1, length(hyb4));

for i = 1: length(par)
    speedup_par(i) = par(1)/ par(i);
end

for i = 1: length(dis)
    speedup_dis(i) = dis(1) / dis(i) ;
end

for i = 1: length(hyb2)
    speedup_hyb2(i) = hyb2(1) / hyb2(i);
end

for i = 1: length(hyb4)
    speedup_hyb4(i) = hyb4(1) / hyb4(i);
end

figure( 1 );
plot( [1 2 4 8 12 16], speedup_par, 'Linewidth', 0.8 );
title( 'Speedup of parallel version  ($t_{seq}/t_{par}$)', 'interpreter', 'latex', 'FontWeight', 'bold' );
ylabel( 'Speedup', 'interpreter', 'latex' );
xlabel( 'Number of threads', 'interpreter', 'latex' );
xticks([1 2 4 8 12 16]);

figure( 2 );
plot( [1 2 4], speedup_dis, 'Linewidth', 0.8 );
title( 'Speedup of distributed version  ($t_{seq}/t_{dis}$)', 'interpreter', 'latex', 'FontWeight', 'bold' );
ylabel( 'Speedup', 'interpreter', 'latex' );
xlabel( 'Number of nodes', 'interpreter', 'latex' );
xticks([1 2 4]);

figure( 3 );
plot( [1 2 4 8 12 16], speedup_hyb2, 'Linewidth', 0.8 );
title( 'Speedup of hybrid version compared to distributed version (2 nodes) ($t_{dis2}/t_{hyb2}$)', 'interpreter', 'latex', 'FontWeight', 'bold' );
ylabel( 'Speedup', 'interpreter', 'latex' );
xlabel( 'Number of threads per node', 'interpreter', 'latex' );
xticks([1 2 4 8 12 16]);

figure( 4 );
plot( [1 2 4 8 12 16], speedup_hyb4, 'Linewidth', 0.8 );
title( 'Speedup of hybrid version compared to distributed version (4 nodes) ($t_{dis4}/t_{hyb4}$)', 'interpreter', 'latex', 'FontWeight', 'bold' );
ylabel( 'Speedup', 'interpreter', 'latex' );
xlabel( 'Number of threads per node', 'interpreter', 'latex' );
xticks([1 2 4 8 12 16]);

