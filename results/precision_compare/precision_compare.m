initial_constants = [100, 50, 20, 10, 5, 2, 1];

[s, git_hash_string] = system('git rev-parse HEAD');
git_hash_string = git_hash_string(1:7);
[s, git_branch_name] = system('git rev-parse --abbrev-ref HEAD');
git_branch_name = strtrim(git_branch_name);

for initial_constant = initial_constants
    clf; set(gcf, 'units', 'pixels', 'position', [200, 200, 1500, 1000])
    
    data = readtable(['/Users/ldam/pairhmm_posit/cmake-build-debug/' num2str(initial_constant) '.txt']);
    datalength = length(data{:,1});

    % Config (define different regions that should be plotted individually)

    % 32_data
    set_name = '32_data';
    main_title = ['Dataset = ' strrep(set_name, '_', '\_') ' - Initial constant = 2^{' num2str(initial_constant) '}' newline 'git commit ' git_hash_string ' / ' git_branch_name];
    pngname = set_name+"_"+git_hash_string+"-"+git_branch_name+"_const"+num2str(initial_constant);
    regions = [1 datalength; 1 datalength; 1 datalength; 1 datalength; 1 datalength; 1 datalength]; % start and end indices of each region
    titles = ["All intermediate values", "Phred scores", "M[1][c]", "Y[1][c]", "M[27], X[27], Y[27]", "Result Accumulation"]; % region labels
    showlabels = [false, true, true, true, true, true]; % toggle X-axis variable names for each region
    wide = [true, true, false, false, true, true]; % horizontally wide plot for region
    
    filter_text1 = ["",     "p",    "M[1]", "Y[1]", "M[27]", "result"]; % filtered based on the first (name) column
    filter_text2 = ["NONE", "NONE", "NONE", "NONE", "X[27]", "NONE"];
    filter_text3 = ["NONE", "NONE", "NONE", "NONE", "Y[27]", "NONE"];

    % 1_data
%     set_name = '1_data';
%     main_title = ['Dataset = ' strrep(set_name, '_', '\_') ' - Initial constant = 2^{' num2str(initial_constant) '}' newline 'git commit ' git_hash_string ' / ' git_branch_name];
%     pngname = set_name+"_"+git_hash_string+"-"+git_branch_name+"_const"+num2str(initial_constant);
%     regions = [1 datalength]; % start and end indices of each region
%     titles = ["All intermediate values"]; % region labels
%     showlabels = [true]; % toggle X-axis variable names for each region
%     wide = [true]; % horizontally wide plot for region
%     filter_text1 = [""]; % filtered based on the first (name) column
%     filter_text2 = ["NONE"];
%     filter_text3 = ["NONE"];

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
        if filter_text1(row) == ""
            data_x_filtered = find(ones(datalength,1) >= 0);
            data_matching = data_x_filtered >= start_idx & data_x_filtered <= end_idx;
            data_x_filtered = find(data_matching);
        else
            data_matching1_1 = strncmp(data{:,1}, filter_text1(row), strlength(filter_text1(row)));
            data_matching1_2 = strncmp(data{:,1}, filter_text2(row), strlength(filter_text2(row)));
            data_matching1_3 = strncmp(data{:,1}, filter_text3(row), strlength(filter_text3(row)));
            
            data_matching1 = data_matching1_1 | data_matching1_2 | data_matching1_3;
                        
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
    print(pngname, '-dpng');
end

close;

