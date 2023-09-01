function str=paddedNum2Str(num,len)
%This function converts the given integer to a string and pads it with
%zeros to make it reach the specified length

str=num2str(num);
if length(str)<len
    padlen=len-length(str);
    str1(1:padlen)='0';
    str=[str1 str];
end
