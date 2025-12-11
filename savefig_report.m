function savefig_report(filename, folder, varargin)
% savefig_report - Save figure with standard report settings (A4)
% Usage: savefig_report('myplot.pdf', './Figures')
%        savefig_report('myplot.pdf', './Figures', 'Width', 7, 'Height', 9)

% Default settings (A4: 210mm x 297mm = 8.27" x 11.69")
paper_width = 8.27;
paper_height = 11.69;

% Create folder if it doesn't exist
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Parse optional arguments
for i = 1:2:length(varargin)
    if strcmp(varargin{i}, 'Width')
        paper_width = varargin{i+1};
    elseif strcmp(varargin{i}, 'Height')
        paper_height = varargin{i+1};
    end
end


% Apply settings
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [paper_width, paper_height]);
set(gcf, 'PaperPositionMode', 'auto');
set(gca, 'FontSize', 14, 'LineWidth', 1.5);

% Make all text bold
set(get(gca, 'XLabel'), 'FontWeight', 'bold');
set(get(gca, 'YLabel'), 'FontWeight', 'bold');
set(get(gca, 'ZLabel'), 'FontWeight', 'bold');
set(get(gca, 'Title'), 'FontWeight', 'bold', 'FontSize', 18);


% Thicken plot lines
h = findobj(gca, 'Type', 'line');
set(h, 'LineWidth', 2);

% Save with high quality using exportgraphics
filepath = fullfile(folder, filename);
exportgraphics(gcf, filepath, 'Resolution', 300, 'BackgroundColor', 'white');
fprintf('Saved: %s\n', filepath);
end
