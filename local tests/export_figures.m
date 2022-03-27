figurepath = '../../figures/forces'; 
filepaths = dir(fullfile(cd, figurepath)); 

f =  findobj('type','figure');
set(0, 'defaultLegendInterpreter', 'latex'); 
set(0, 'defaultTextInterpreter', 'latex'); 
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(0, 'defaultColorbarTickLabelInterpreter', 'latex'); 
set(0, 'defaultLineLineWidth', 1.5); 

%%
for k = 1:length(f)
    figure(f(k)); 
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    if isempty(get(f(k), 'name')) 
        filename = fullfile(figurepath, sprintf('Figure %i', k)); 
    else
        filename = fullfile(figurepath, get(f(k), 'name')); 
    end
%     exportgraphics(f(k), [filename, '.pdf'], 'Resolution', 600); 
    exportgraphics(f(k), [filename, '.pdf'], 'ContentType', 'vector'); 
    exportgraphics(f(k), [filename, '.png'], 'Resolution', 600); 
    saveas(f(k), filename); 
end