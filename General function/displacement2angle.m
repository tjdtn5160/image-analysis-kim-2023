function ang=displacement2angle(dx,dy)
%displacement2angle returns the anlge (relative to East) for a given
%displacement (dx,dy)

ang=atan(dy/dx);
ang=ang/pi*180;
if dx<0
    ang=180+ang;
end

if ang<0
    ang=ang+360;
end