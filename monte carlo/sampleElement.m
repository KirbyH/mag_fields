function Z = sampleElement(N)
% Samples one of the four common GCR elements: H, He, C, Fe from
% distributions given by Simpson, 1983
% 
% INPUTS : 
%     N (optional) : length of sample to return
% 
% OUTPUTS : 
%     Z : atomic number (AMU)

% ========== DISTRIBUTION OF ELEMENTS ============
el(1) = 540;    % H
el(2) = 26;     % He
el(6) = 2.2;    % C
el(26) = 0.12;  % Fe
% ================================================

el_pdf = el/sum(el); 
el_cdf = zeros(size(el_pdf)); 
el_cdf(1) = el_pdf(1); 
for ii = 2:length(el_pdf)
    el_cdf(ii) = el_cdf(ii-1) + el_pdf(ii); 
end

if ~exist('N', 'var')
    N = 1; 
end
ksi = rand(N,1); 
Z = zeros(N,1); 
for ii = 1:N
    Z(ii) = find(el_cdf > ksi(ii), 1); 
end

end
