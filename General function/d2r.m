function radians = d2r(varargin)
% Examples
% d2r(180,30,45,60,90)
% d2r(20:5:180), d2r(30)
degang = cell2mat(varargin);

if nargin < 1
    error('Invalid number of inputs');
end

if ischar(degang)
    error('Input must be a angle or a vector of angles')
end


radians = degang*pi/180;

