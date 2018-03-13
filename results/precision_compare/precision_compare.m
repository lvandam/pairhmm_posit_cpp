initial_constants = [100, 50, 20, 10, 5, 2, 1];

for initial_constant = initial_constants
    clf; set(gcf, 'units', 'pixels', 'position', [200, 200, 1500, 1000])
    
    data = readtable(['/Users/ldam/pairhmm_posit/cmake-build-debug/' num2str(initial_constant) '.txt']);
    datalength = length(data{:,1});

    % Config (define different regions that should be plotted individually)

    % 32_data
    main_title = ['Dataset = 32\_data.txt - Initial constant = 2^{' num2str(initial_constant) '}'];
    regions = [1 datalength; 4543 4638; 4543 4638; 7039 7134; 7135 datalength]; % start and end indices of each region
    titles = ["All intermediate values", "M[1][c]", "Y[1][c]", "M[27], X[27], Y[27]", "Result Accumulation"]; % region labels
    showlabels = [false, true, true, true, true]; % toggle X-axis variable names for each region
    wide = [true, false, false, true, true]; % horizontally wide plot for region
    filter_text = ["", "M", "Y", "", ""]; % filtered based on the first (name) column

    % 1_data
    % main_title = 'Dataset = 1\_data.txt - Initial constant = 2^{1}';
    % regions = [1 datalength]; % start and end indices of each region
    % titles = ["All intermediate values"]; % region labels
    % showlabels = [true]; % toggle X-axis variable names for each region
    % wide = [true]; % horizontally wide plot for region
    % filter_text = [""]; % filtered based on the first (name) column

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

        % Match with any given filter (start string)
        if filter_text(row) == ""
            data_x_filtered = find(ones(datalength,1) >= 0);
            data_matching = data_x_filtered >= start_idx & data_x_filtered <= end_idx;
            data_x_filtered = find(data_matching);
        else
            data_matching1 = strncmp(data{:,1}, filter_text(row), length(filter_text(row)));
            data_x_filtered = find(data_matching1 >= 0);
            data_matching2 = data_x_filtered >= start_idx & data_x_filtered <= end_idx;
            data_matching = data_matching1 .* data_matching2;
            data_x_filtered = find(data_matching);
        end

        data_float_filtered = data{data_x_filtered, 4};
        data_posit_filtered = data{data_x_filtered, 5};

        scatter(data_x_filtered, data_float_filtered, 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', .7) % float %data{start_idx:end_idx,4},
        scatter(data_x_filtered, data_posit_filtered, 30, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', .7) % posit<32,2> %data{start_idx:end_idx,5},
        legend({'float', 'posit<32,2>'});
        grid on; grid minor;
        set(gca, 'TickLabelInterpreter', 'none');
        if showlabels(row)
            set(gca, 'xtick', data_x_filtered, 'xticklabel', data{data_x_filtered,1})
        else
            set(gca, 'xtick', [])
        end
        xtickangle(45);
        title(titles(row)); ylabel('log_{10}(|relative error|)');
    end

    p = mtit(main_title, 'FontSize', 14, 'Color', 'r', 'XOff', 0, 'YOff', .025);
    
    % Save as PNG
    print(num2str(initial_constant), '-dpng');
end



