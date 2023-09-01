function wellname=row_column_to_wellname(row,col)
%converts row and column numbers to the typical well names (eg. A06)

rowletters={'A','B','C','D','E','F','G','H'};
row_letter=rowletters{row};

if col<10;
   column_number=['0',num2str(col)]; 
else column_number=[num2str(col)];
end

wellname=[row_letter,column_number];

end

