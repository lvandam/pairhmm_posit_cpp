initial_constants = [100, 50, 20, 10, 5, 2, 1];

file_prefix = '5cba05c_nofused_float_posit32-2/';

results_decimal_accuracy = zeros(1,3);

i = 1;
for initial_constant = initial_constants
    data = readtable([file_prefix num2str(initial_constant) '.txt']);
    
    da_float = data{strcmp(data{:,1}, 'result[32]'), 9};
    da_posit = data{strcmp(data{:,1}, 'result[32]'), 10};
    
    results_decimal_accuracy(i,:) = [initial_constant da_float da_posit];
    
    i = i + 1;
end

clf; h = figure;
set(h, 'units', 'pixels', 'position', [200, 200, 600, 200]);
hold on
scatter(results_decimal_accuracy(:,1), results_decimal_accuracy(:,2), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .7) % float
scatter(results_decimal_accuracy(:,1), results_decimal_accuracy(:,3), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .7) % posit
xticks(sort(initial_constants));
legend({'float', 'posit<32,2>'});
grid on; grid minor;
set(gca, 'TickLabelInterpreter', 'none');
ylabel('decimal accuracy');
xlabel('initial constant');

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
print(h, '32_data_nofused_initial_constants', '-dpdf', '-r0')