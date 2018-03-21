file_prefix_es2 = 'd195ce7_fused_float_posit32-2/';
file_prefix_es3 = 'd195ce7_fused_float_posit32-3/';
pdfname = '32_data_fused_d195ce7_initial_constants_float_posit_32-3_posit_32-3';

results_decimal_accuracy = zeros(1,4);

i = 1;
for initial_constant = 5:5:100
    data_es2 = readtable([file_prefix_es2 num2str(initial_constant) '.txt']);
    data_es3 = readtable([file_prefix_es3 num2str(initial_constant) '.txt']);
    
    da_float = data_es2{strcmp(data_es2{:,1}, 'result[32]'), 9};
    da_posit_es2 = data_es2{strcmp(data_es2{:,1}, 'result[32]'), 10};
    da_posit_es3 = data_es3{strcmp(data_es3{:,1}, 'result[32]'), 10};
    
    results_decimal_accuracy(i,:) = [initial_constant da_float da_posit_es2 da_posit_es3];
    
    i = i + 1;
end



clf; h = figure;
set(h, 'units', 'pixels', 'position', [200, 200, 600, 200]);
hold on
scatter(results_decimal_accuracy(:,1), results_decimal_accuracy(:,2), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .7) % float
scatter(results_decimal_accuracy(:,1), results_decimal_accuracy(:,3), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .7) % posit<32,2>
scatter(results_decimal_accuracy(:,1), results_decimal_accuracy(:,4), 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g', 'MarkerFaceAlpha', .7) % posit<32,3>
xticks(5:5:100);
legend({'float', 'posit<32,2>', 'posit<32,3>'});
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
print(h, pdfname, '-dpdf', '-r0')