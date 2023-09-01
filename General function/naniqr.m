function [iqr_of_trace] = naniqr(vector_of_data)
new_vector=vector_of_data(~isnan(vector_of_data));
iqr_of_trace=iqr(new_vector);
end