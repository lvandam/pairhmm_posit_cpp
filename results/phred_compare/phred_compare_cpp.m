clf;

data = csvread('phred_float_posit.txt', 1);

Q = data(:,1);
dE_f = data(:,2);
dE_f_log = data(:,3);
dE_p = data(:,4);
dE_p_log = data(:,5);

scatter(Q, dE_f_log, 'r');
hold on;
scatter(Q, dE_p_log, 'b');
xlabel('Phred Quality Score');
ylabel('log_{10}(relative error)');
grid on; grid minor;
legend({'float32','posit<32,2>'}, 'Location', 'northwest');
