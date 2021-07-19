%test entropy function
result=[];
% i=0:0.01:1
% i=0:0.001:0.1
for i=0:0.01:1
    e=entropy(i);
    result=[result; i e];
end
function [entropy] = entropy(input)

    % expects probability distribution as a [] x 1 array 
    % outputs sum of every entropy of every probability
    temp=input;
    % temp is the probability
    entropy= -temp .* log2(temp);
%     entropy= sum(entropy);%output sum of all entropy
    
end
% function [entropy] = entropy(input)
% 
%     % expects probability distribution as a [] x 1 array 
%     % outputs sum of every entropy of every probability
%     temp=input;
%     % temp is the probability
%     entropy= -log2(temp);
% %     entropy= sum(entropy);%output sum of all entropy
%     
% end