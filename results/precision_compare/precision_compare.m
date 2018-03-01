clf;

data = readtable('pairhmm_values_32_data_2-1.txt');

datalength = length(data{:,1});

% Config (define different regions that should be plotted individually)
main_title = 'Dataset = 32\_data.txt - Initial constant = 2^{1}';
regions = [1 datalength; 4543 4638; 7039 7134; 7135 datalength]; % start and end indices of each region
titles = ["All intermediate values", "M[1], X[1], Y[1]", "M[27], X[27], Y[27]", "Result Accumulation"]; % region labels
showlabels = [false, true, true, true]; % toggle X-axis variable names for each region
wide = [true, true, true, true]; % horizontally wide plot for region

% Process
wide_count = nnz(wide == true);

total_subplots = 0;
for el = wide
    if el == true
        if mod(total_subplots+1,2) == 0
            if total_subplots == 0
                total_subplots = 2;
            else
                total_subplots = total_subplots + 3;
            end
        else
            total_subplots = total_subplots + 2;
        end
    else
        total_subplots = total_subplots + 1;
    end
end

subplot_cols = 2;
subplot_rows = ceil(total_subplots / 2);

subplotnum = 1;
for row = 1:size(regions,1)
    if wide(row)
        if mod(subplotnum,2) == 0
            subplotnum = subplotnum + 1;
        end
        subplot(subplot_rows, subplot_cols, [subplotnum subplotnum+1]);
        subplotnum = subplotnum + 2;
    else
        subplot(subplot_rows, subplot_cols, subplotnum);
        subplotnum = subplotnum + 1;
    end

    start_idx = regions(row,1);
    end_idx = regions(row,2);
    
    % Plot
    hold on
    scatter(start_idx:end_idx, data{start_idx:end_idx,4}, 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .7) % float
    scatter(start_idx:end_idx, data{start_idx:end_idx,5}, 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .7) % posit<32,2>
    legend({'float', 'posit<32,2>'});
    grid on; grid minor;
    set(gca, 'TickLabelInterpreter', 'none');
    if showlabels(row)
        xlabel('name');
        set(gca, 'xtick', start_idx:end_idx, 'xticklabel', data{start_idx:end_idx,1})
    end
    xtickangle(45);
    title(titles(row)); ylabel('log_{10}(|relative error|)');
end

p = mtit(main_title, 'FontSize', 14, 'Color', 'r', 'XOff', 0, 'YOff', .025);
