function [] = savePlot(fh,filename,w,h,fspec,res)
%--------------------------------------------------------------------------
% Purpose: REU 2016 figure printer,
%          This function saves a figure given a figure handle 'fh' 
%
% Inputs:  fh  =  figure handle
%          w   =  figure width (in)
%          h   =  figuer height (in)
%
% Author:  Dalton Anderson
% Date:    July 7, 2016
%--------------------------------------------------------------------------
w = w - 0.1;
h = h - 0.1;
% if file name and destination folder not specified
datadir = pwd; % current directory
if nargin < 5
%     fspec = '-dpng'; 
%     res = '-r1000';
    fspec = '-dpdf';
    res = '-r0';
end
 
set(fh,'Units','centimeters')
pos = get(fh,'position');
set(fh,'Units','centimeters',...
       'Position',[pos(1),pos(2),w,h],...
       'outerposition',[pos(1)-0.5,pos(2)-0.5,w+1,h+1])
pos = get(fh,'Position');
set(fh,'PaperUnits','centimeters',...
       'PaperPosition',[0,0,w,h],...
       'PaperSize',[w,h])
 
% Print figure
print(fh,[datadir,'/',filename],fspec,res)
% save2pdf([datadir,'/',filename],fh,600)
end