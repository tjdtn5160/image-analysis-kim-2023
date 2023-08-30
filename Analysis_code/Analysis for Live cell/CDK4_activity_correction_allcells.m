function [corrected_cdk4] = CDK4_activity_correction_allcells(cdk4_trace,cdk2_trace)
%Corrects Cdk4 activity based on CDK2 activity from the same cell
robust_fit_varible1=0;%0.2786+0.2537;
robust_fit_variable2=(0.3279+0.1715)/2;

% length_of_trace=size(cdk4_trace,2);

corrected_cdk4=NaN(size(cdk4_trace));
corrected_cdk4=cdk4_trace - (robust_fit_variable2.*cdk2_trace + robust_fit_varible1);
end

