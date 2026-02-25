function [probability] = func_recPb (doaSource, doaEstimation, interval)

N_Trial = size(doaEstimation,1);
N_Source = size(doaEstimation,2);
hit = zeros(1,N_Trial);
doaSource = sort(doaSource,'ascend');
for ii = 1:N_Trial
    doaEst = doaEstimation(ii,:);
    doaEst = sort(doaEst,'ascend');
    temp = zeros(1,N_Source);
    for k = 1:N_Source
        if abs(doaEst(k)-doaSource(k))<interval/2
            temp(k) = temp(k) + 1;
        end
    end
    if sum(temp) == N_Source
        hit(ii) = 1;
    end
end
probability = mean(hit);
end