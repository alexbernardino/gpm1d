% analisis of the parameter space of normalized 1D two-gaussian mixture
% model.

Nw = 30;
Nm = 30;
Nv = 30;

tab = create_table(Nw,Nm, Nv, 0);

%arrays to store the samples points
%weights
w1_array = [];
w2_array = [];
%means
m1_array = [];
m2_array = [];
%variances
v1_array = [];
v2_array = [];
%standard deviations
s1_array = [];
s2_array = [];
idx = 0;

% partition the ]0,1[ range into Nw equidistant values for w1
delta_w1 = 1/(1+Nw);
w1 = delta_w1:delta_w1:1-delta_w1;
w2 = 1-w1; %w2 is uniquely determined from w1
%The upper bound for m1 is uniquely determined by w1 and w2
m1_ub= sqrt(w2)./sqrt(w1);
%resolutions of the m1 partitions of [0, m1_ub[ range into Nm equidistant values
m1_delta = m1_ub/Nm;

for i = 1:length(w1)
    %partition the [0, m1_ub(i)[ range into Nm equidistant values
    m1 = 0:m1_delta(i):m1_ub(i)-m1_delta(i);
    m2 = -m1*w1(i)/w2(i); %m2 uniquely determined from m1, w1, w2
    %The upper bound for v1 is uniquely determined by m1, m2, w1, w2
    v1_ub = (1-w1(i)*m1.^2-w2(i)*m2.^2)/w1(i); %is it better to discretize var or std ?
    %s1_ub = sqrt(v1_ub);
    for j = 1:length(m1)
        %partition the ] 0, v1_ub(i,j) [ range into Nv equidistant values
        v1_delta = v1_ub(j)/(Nv+1);
        v1 = v1_delta: v1_delta : v1_ub(j)-v1_delta;
        v2 = (1 - w1(i)*v1 - w1(i)*m1(j)^2 - w2(i)*m2(j)^2)/w2(i); %v1 uniquely determined by all other vals
        %partition the ] 0, s1_ub(i,j) [ range into Nv equidistant values
        %s1_delta = s1_ub(j)/(Nv+1);
        %s1 = s1_delta: s1_delta : s1_ub(j)-s1_delta;
        %s2 = sqrt((1 - w1(i)*s1.^2 - w1(i)*m1(j)^2 - w2(i)*m2(j)^2)/w2(i)); %s1 uniquely determined by all other vals
        for k = 1:length(v1)
            idx = idx+1;
            w1_array(idx) = w1(i);
            w2_array(idx) = w2(i);
            m1_array(idx) = m1(j);
            m2_array(idx) = m2(j);
            v1_array(idx) = v1(k);
            v2_array(idx) = v2(k);
            %s1_array(idx) = s1(k);
            %s2_array(idx) = s2(k);
        end
    end
end

idx

%var1 upper bound - computed from the constraint that total variance must
%be unitary

%for i = 1:length(w1)
%    for j = 1:length(m1)
%        var2_ub(i,j) = 1/w1(i)+ (w1(i)-w2(i))/w2(i)*m1(j)^2;
%    end;
%end;
%figure
%mesh(w1,m1,var2_ub);
%figure
%contour(w1,m1,var2_ub, 100);

% figure;
% plot(w1_array(:), w2_array(:), '.');
% xlabel('w_1');
% ylabel('w_2');
% figure;
% plot(m1_array(:), m2_array(:), '.');
% xlabel('m_1');
% ylabel('m_2');
% figure;
% plot(sqrt(v1_array(:)), sqrt(v2_array(:)), '.');
% xlabel('\sigma_1');
% ylabel('\sigma_2');
% figure;
% plot(w1_array(:), m1_array(:), 'b.', w1_array(:), m2_array(:), 'g.');
% xlabel('w_1');
% ylabel('m_1, m_2');
% figure;
% plot(w1_array(:), sqrt(v1_array(:)), 'b.', w1_array(:), sqrt(v2_array(:)), 'g.');
% xlabel('w_1');
% ylabel('\sigma_1, \sigma_2');
% figure;
% plot(m1_array(:), sqrt(v1_array(:)), 'b.', m1_array(:), sqrt(v2_array(:)), 'g.');
% xlabel('m_1');
% ylabel('\sigma_1, \sigma_2');

%figure; hist(w1_array(:)); title('w_1');
%figure; hist(w2_array(:)); title('w_2');
%figure; hist(m1_array(:)); title('m_1');
%figure; hist(m2_array(:)); title('m_2');
%figure; hist(sqrt(v1_array(:))); title('\sigma_1');
%figure; hist(sqrt(v2_array(:))); title('\sigma_2');

figure; hist(tab.w1(:)); title('w_1');
figure; hist(tab.w2(:)); title('w_2');
figure; hist(tab.m1(:)); title('m_1');
figure; hist(tab.m2(:)); title('m_2');
figure; hist(tab.v1(:)); title('\sigma_1');
figure; hist(tab.v2(:)); title('\sigma_2');




