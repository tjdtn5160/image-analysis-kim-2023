function [rowstring,colstring,sitestring] = wellnum2strRCS_3(row,col,site)
rowstringtable = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
colstringtable = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12'};
sitestringtable = {'0000','0001' '0002' '0003' '0004' '0005' '0006' '0007' '0008' '0009' '0010' '0011' '0012' '0013' '0014' '0015' '0016' '0017' '0018' '0019' '0020' '0021' '0022' '0023' '0024' '0025' '0026' '0027' '0028' '0029' '0030' '0031' '0032' '0033' '0034' '0035' '0036'};
rowstring = rowstringtable{row};
colstring = colstringtable{col};
sitestring = sitestringtable{site};
end