
%%%%%This package written by robert wilson, freely available at 
% https://github.com/d-r-b-o-b/UsefulMatlabFunctions%%

fontname = 'Helvetica';
fontsize = 16;
ABCfontsize = 24;
fontweight = 'normal';
linewidth = 3;

matlabGrey = [1 1 1]*215/255;

global arminn arminn_pink


arminn = [119, 171, 189]/256;
arminn_pink = [230, 57, 70]/256;
% lighten blue


global orange
orange = [0.906 0.463 0.247];
% colormap gray
% CC = colormap;
% CM = (CC).*repmat((1-[0.906 0.463 0.247]/0.906), size(CC,1),1);
% colormap(1-CM)

set(0, 'defaultfigurecolor', 'w', ...
    'defaultaxesfontsize', fontsize, ...
    'defaultaxesfontweight', fontweight, ...
    'defaultaxesfontname', fontname, ...
    'defaultaxestickdir', 'out', ...
    'defaultaxesbox', 'off', ...
    'defaultaxesydir', 'normal', ...
    'defaultlinelinewidth', linewidth, ...
    'defaultlinemarkersize', 30, ...
    'defaultfigureposition', [811   486   618   500])

