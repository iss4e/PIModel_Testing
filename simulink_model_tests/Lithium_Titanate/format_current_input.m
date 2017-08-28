load('test_current_input.mat');

n = length(current);

newarray = zeros(2,n);

for i = 1:n
    newarray(:,i) = [i;current(i)];
end

save('test_current_formatted.mat', 'newarray', '-v7.3');