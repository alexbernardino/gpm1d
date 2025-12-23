function table = create_table(Nw, Nm, Nv, discretize_std)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%% Generate the hypotheses table
 %m1 hypotheses are positive.

% partition the ]0,1[ range into Nw equidistant values for w1
delta_w = 1/(1+Nw);
w1_array = delta_w*linspace(1,Nw,Nw)';
w2_array = 1-w1_array; %w2 is uniquely determined from w1
%The upper bound for m1 is uniquely determined by w1 and w2
m1_ub_array = sqrt(w2_array)./sqrt(w1_array);
%resolutions of the m1 partitions of [0, m1_ub[ range into Nm equidistant values
delta_m_array = m1_ub_array/Nm;
%partition the [0, m1_ub_array(i)[ range into Nm equidistant values
m1_matrix = delta_m_array.*linspace(0,Nm-1,Nm);
%m2 uniquely determined from m1, w1, w2
m2_matrix = -m1_matrix.*w1_array./w2_array;
%The upper bound for v1 is uniquely determined by m1, m2, w1, w2
v1_ub_matrix = ones(Nw,Nm)./w1_array - (m1_matrix.^2)./w2_array;
s1_ub_matrix = sqrt(v1_ub_matrix);
if discretize_std 
    %partition the ] 0, s1_ub_delta(i,j) [ range into Nv equidistant values
    delta_s_matrix = s1_ub_matrix/(Nv+1);
    s1_tensor = delta_s_matrix.*reshape(logspace(log10(1),log10(Nv),Nv),1,1,[]);
    v1_tensor = s1_tensor.^2;
else %discretize the variance range
   %partition the ] 0, v1_ub_delta(i,j) [ range into Nv equidistant values
   delta_v_matrix = v1_ub_matrix/(Nv+1);
   v1_tensor = delta_v_matrix.*reshape(logspace(log10(1),log10(Nv),Nv),1,1,[]);
   s1_tensor = sqrt(v1_tensor);
end
%v2 uniquely determined by all other vals
v2_tensor = (ones(Nw,Nm,Nv) - w1_array.*v1_tensor - w1_array.*(m1_matrix.^2) - w2_array.*(m2_matrix.^2))./w2_array;
s2_tensor = sqrt(v2_tensor);

table.w1 = w1_array;
table.w2 = w2_array;
table.m1 = m1_matrix;
table.m2 = m2_matrix;
table.v1 = v1_tensor;
table.v2 = v2_tensor;
table.s1 = s1_tensor;
table.s2 = s2_tensor;
