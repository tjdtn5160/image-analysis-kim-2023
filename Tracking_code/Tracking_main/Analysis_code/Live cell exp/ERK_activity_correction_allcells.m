function [corrected_erk] = ERK_activity_correction_allcells(erk_trace,cdk2_trace)
%Corrects Cdk4 activity based on CDK2 activity from the same cell
robust_fit_variable1=0;%(0.2548+0.4284)/2;
robust_fit_variable2=(0.3244+0.2840)/2;

% length_of_trace=size(cdk4_trace,2);

corrected_erk=NaN(size(erk_trace));
corrected_erk=erk_trace - (robust_fit_variable2.*cdk2_trace + robust_fit_variable1);
end

