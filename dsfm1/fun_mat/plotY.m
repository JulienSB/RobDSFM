function plotY(Y, ni, time, plot_dim)
    % Create a plot of the data (matrix representation).
    % Args:
    %   Y:          data matrix.
    % 	time:       time-point at which the image is displayed.
    % 	plot_dim:   dimensions of the matrix plot.

    T = size(Y, 1);
    if nargin < 3
       time = 1;
    end

    if nargin < 4
        if ni > 1
            plot_dim = [4 4];
        else
            plot_dim = [1 1];
        end
    end
    
    Image = getImage(Y, ni);
    plot3D(Image, time, plot_dim)
end