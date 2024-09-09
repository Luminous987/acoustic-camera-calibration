function outlier = finderror(measurement)
%%
addpath('func');
% load("./9.3measurement.mat");
outlier = [];
n=1;
for i = 1:7   
    measured_values =  measurement(:,i);
    % c = 10*mean(measured_values);
    d = 0.005;
    for k = 1:6:length(measured_values)
        s_51 = measured_values(k+4)-measured_values(k);
        s_53 = measured_values(k+4)-measured_values(k+2);
        s_62 = measured_values(k+5)-measured_values(k+1);
        s_64 = measured_values(k+5)-measured_values(k+3);
        if abs(s_51) > d || abs(s_53) > d || abs(s_62) > d || abs(s_64) > d
            if ismember(ceil(k/6),outlier) == 0 
                outlier(n) = ceil(k/6);
                n = n+1;               
            end
        end
    end
    % plot(measured_values, 'b-o', 'DisplayName', '测量值');
    hold on
end
outlier = sort(outlier, 'descend');
sort(outlier)
end