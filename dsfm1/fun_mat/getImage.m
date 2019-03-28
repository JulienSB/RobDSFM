function Image = getImage(Y, ni)
    % Create 4 dimensional representation of the data matrix.
    % Args:
    %   Y:      matrix representation of the data.
    %   ni:     vector of the number of voxels in each dimension.
    %           ni should have length 3.
    % Output:
    %   Image:  4-dimensional representation of the data.
    
    T = size(Y, 1);
    Image = zeros(ni(1), ni(2), ni(3), T);
    for t = 1:T
        Image(:, :, :, t) = reshape(Y(t,:), [ni(1), ni(2), ni(3)]); 
    end
end