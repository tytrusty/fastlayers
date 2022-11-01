function [P,W_RGB,W_RGBXY] = precompute(filename)
    
    if nargin < 1
        filename = '16x16.png';
    end

    A = imread(filename);

    % Flatten image into (W * H) x 3 matrix
    RGB = reshape(permute(A, [3 1 2]), 3, [])';
    
    [width,height] = size(A,[1 2]);

    % XY coordinates of image pixels
    [X,Y] = meshgrid(1:width, 1:height);

    RGBXY = double([RGB X(:) Y(:)]);
    RGBXY(:,1:3) = RGBXY(:,1:3) ./ 255;
    RGBXY(:,4) = RGBXY(:,4) ./ width;
    RGBXY(:,5) = RGBXY(:,5) ./ height;

    % Compute 5D convex hull
    F_5D = convhulln(RGBXY);
    convhull_verts = unique(F_5D(:));
    
    % Construct simplified convex hull
    V_RGB = RGBXY(convhull_verts, 1:3);
    F = convhull(V_RGB);
    target_faces = 30;
    ratio = target_faces / size(F,1);
    [V_3D,F_3D] = decimate_libigl(V_RGB,F,ratio,'Method','progressive-hulls');

    % Get colors by clamping to RGB cube
    P = min(max(V_3D,0),1.0); % Clamp to RGB cube
    
    % Delaunay triangulation of 3D RGB convex hull
    T_3D = delaunay(V_3D);

    % Barycentric coordinates of 5D convex hull vertices to 3D
    % triangulation
    [~,I,J,L] = in_element(V_3D, T_3D, V_RGB);
    I = repmat(I,[1,4]);
    J = T_3D(J,:);
    W_RGB = sparse(I(:),J(:),L(:));

    % Delaunay triangulation of convex hull
    V_5D = RGBXY(convhull_verts, :);
    T_5D = delaunayn(V_5D, {'Qt','Qbb','Qc','Qx'});
    
    % Barycentric coordinates weighting pixels to 5D convex hull
%     [J,L] = tsearchn(V_5D,T_5D,RGBXY);
%     I = [1:size(RGBXY,1)]';

    [~,I,J,L] = in_element(V_5D,T_5D, RGBXY, 'Epsilon', 1e-3);
    I = repmat(I,[1,6]);
    sum(isnan(J))
    J(isnan(J)) = 1;
    J = T_5D(J,:);
    W_RGBXY = sparse(I(:),J(:),L(:), size(RGBXY,1), size(V_5D,1));
    [~,name,~] = fileparts(filename);
    fn = ['cache' filesep name '.mat'];
    save(fn,'A','P','W_RGB','W_RGBXY','RGBXY','V_3D','F_3D')
end