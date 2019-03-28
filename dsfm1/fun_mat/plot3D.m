function plot3D(Image, time, plot_dim)
    % Create a plot of an Image (4-dimensional object).
    % Args:
    %   Image:      4-dimensional array, last dimension is time.
    % 	time:       time-point at which the image is displayed.
    % 	plot_dim:   dimensions of the matrix plot.

    dim = size(Image);
    if nargin < 2
       time = 1;
    end

    if nargin < 3
       plot_dim = [1 1];
    end

    for i = 1:min(plot_dim(1) * plot_dim(2), dim(3))
        subplot(plot_dim(1), plot_dim(2), i), imagesc(Image(:, :, i, time))
    end
end