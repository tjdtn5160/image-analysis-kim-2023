function degrees = r2d(varargin)
% Examples: 
% r2d(pi,pi/6,pi/4,pi/2,2*pi/3);
% r2d(1.48), r2d(0:pi/6:pi)

radang = cell2mat(varargin);

if nargin < 1
    error('Invalid number of inputs');
end

if ischar(radang)
    error('Input must be a angle or a vector of angles')
end


degrees = radang*180/pi;

