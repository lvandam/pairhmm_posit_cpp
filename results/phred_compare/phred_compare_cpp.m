clf;

data = csvread('phred_float_posit.txt', 1);
pdfname = 'phred_float_posit';

Q = data(:,1);
da_f = data(:,2);
da_p2 = data(:,3);
da_p3 = data(:,4);

clf; h = figure;
set(h, 'units', 'pixels', 'position', [200, 200, 600, 200]);

hold on;
scatter(Q, da_f, 'r');
scatter(Q, da_p2, 'b');
scatter(Q, da_p3, 'g');
xlabel('Phred Quality Score');
ylabel('decimal accuracy');
grid on; grid minor;
legend({'float','posit<32,2>','posit<32,3>'}, 'Location', 'northeast');

% Save as PDF
% Tight boundaries
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'points ', 'PaperSize', [pos(3), pos(4)])
print(h, pdfname, '-dpdf', '-r0')