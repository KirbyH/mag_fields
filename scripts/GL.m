function globalvars = GL(name, val)
% Global variable function equivalent. Stores default values in a structure
% G, and values can be updated one-by-one by use of the input parameters. 
% Influenced by: https://www.mathworks.com/matlabcentral/answers/369955-once-more-avoid-global-variable
% 
% INPUTS : 
%     name : name of new variable in structure
%     val : value of corresponding variable
% 
% OUTPUTS : 
%     G : Structure of all global variables
%     
% Kirby Heck
% 02/21/2021 

persistent G_
if nargin == 2
    G_.(name) = val; 
end

globalvars = G_; 
end

