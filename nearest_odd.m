function nearestOddArray = nearest_odd(arr)
    % Initialize output array
    nearest_int = nearest(arr);
    for i = find(mod(nearest_int,2) == 0)
        cond = arr(i)<nearest_int(i);
        nearest_int(i) = (nearest_int(i)-1)*cond + (nearest_int(i)+1)*(~cond);
        assert(mod(nearest_int(i),2) ~= 0, "int: "+ arr(i))
    end
    nearestOddArray = nearest_int;
end