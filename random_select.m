function [training_x,training_y,test_x,test_y]=random_select(x,y,size)
    mylength = length(y);
    total_index = [1:mylength];
	res = [];
	while length(res) ~= size
        tmp_res = unique(randi(mylength,1,size));
        tmp_res = sort(tmp_res);
        res = [res;tmp_res'];
        res = unique(res);
        if length(res) >= size
            res = res(1:size);
            res = sort(res);
        end
    end
    total_index(res) = [];
    training_x = x(res,:);
    test_x = x(total_index,:);
    training_y = y(res);
    test_y = y(total_index);
end