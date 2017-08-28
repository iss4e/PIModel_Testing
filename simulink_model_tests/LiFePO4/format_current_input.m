load('test_current_input.mat');

n = length(current);

newarray = zeros(2,n);

index = 1;
for i = 4002:n
    newarray(:,i) = [index*10;current(i)];
    index = index+1;
end

newarray = newarray(:,4002:n);

save('test_current_formatted.mat', 'newarray', '-v7.3');