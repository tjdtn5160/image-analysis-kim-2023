function [mad_of_trace] = nanmad(vector_of_data)
new_vector=vector_of_data(~isnan(vector_of_data));
mad_of_trace=mad(new_vector);
end

