function [Y, X, Image,Data]=getDATA(xaxis, yaxis, zaxis, info, T, varargin)
    % Create data matrix from an fMRI object.
    % Args:
    %   xaxis:       
    %   yaxis:
    %   zaxis:  slices
    %   info:   an object produced by load_untouch_nii.
    %   T:      number of time points taken.
    % Output:
    %   Y:      a matrix representation of the data.
    %   X:      a matrix of the voxel index.
    %   Image:  a 4 dimensional matrix representation of the data.
    %   Data:   
    
    p = inputParser;
    defaultStart = 10;

    addRequired(p,'xaxis',@isnumeric);
    addRequired(p,'yaxis',@isnumeric);
    addRequired(p,'zaxis',@isnumeric);
    addRequired(p,'info');
    addRequired(p,'T',@isnumeric);
    addOptional(p,'start',defaultStart,@isnumeric);


    parse(p,xaxis,yaxis,zaxis,info,T,varargin{:});
    ni1=length(p.Results.xaxis);
    ni2=length(p.Results.yaxis);
    ni3=length(p.Results.zaxis);
    J=ni1*ni2*ni3;

    Y=zeros(p.Results.T,J);
    if(ni3==1) %2D
        Image=zeros(ni1,ni2,p.Results.T); 
        X=[kron(ones(1,ni2),1:ni1)',kron(1:ni2,ones(ni1,1)')'];
    else  % 3D
        Image=zeros(ni1,ni2,ni3,p.Results.T); 
        X=[kron(ones(1,ni2*ni3),1:ni1)',...
        kron(ones(ni3,1),kron(1:ni2,ones(ni1,1)')'),...
        kron(1:ni3,ones(ni1*ni2,1)')'];
    end



    %% nii.gz files:
    %%%%%%%%%%%%%%%%    

    fullData=p.Results.info.img;

    if ni3==1
        Data=fullData(p.Results.xaxis,p.Results.yaxis,p.Results.start:(p.Results.T+p.Results.start-1));
        Image=Data;
         for t=1:p.Results.T
             Data2=Data(:,:,t);
             Y(t,:) = Data2(:);
         end
    else
        Data=fullData(p.Results.xaxis,p.Results.yaxis,p.Results.zaxis,p.Results.start:(p.Results.T+p.Results.start-1));
        Image=Data;    
         for t=1:p.Results.T
             Data2=Data(:,:,:,t);
             Y(t,:) = Data2(:);
         end
    end
end





