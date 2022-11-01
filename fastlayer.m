function fastlayer(matfile)
    
    if nargin < 1
%         filename = 'ngc6543a.jpg';
        %filename = 'dog.jpg';
        matfile = 'cache/16x16.mat';
    end
    mat = load(matfile);
    A = mat.A;
    P = mat.P;
    V = mat.V_3D;
    F = mat.F_3D;
    RGBXY = mat.RGBXY;
    W_RGB = mat.W_RGB;
    W_RGBXY = mat.W_RGBXY;

    [width,height] = size(A,[1 2]);

    % Plot all pixel RGB values
    figure(1); clf; hold on; axis equal;
    scatter3(RGBXY(:,1),RGBXY(:,2),RGBXY(:,3),0.3, RGBXY(:,1:3));
    tsurf(F,V,'FaceColor','none')
    scatter3(V(:,1), V(:,2), V(:,3), 55, V,'filled');

    % Plot just the cvx hull
    figure(2); clf; hold on; axis equal;
    tsurf(F,V,'FaceColor','none')
    scatter3(V(:,1), V(:,2), V(:,3), 55, V,'filled',...
        'ButtonDownFcn',@oncontrolsdown);

    figure(3);
    plt = imshow(A);

    figure(4);
    I = W_RGBXY * W_RGB * P;
    I = reshape(I, [width, height, 3]);
    A_normalized = double(A) ./ 255;
    error = norm(I(:) - A_normalized(:))
    plt_rec = imshow(I);

    function oncontrolsdown(src,ev)
        disp('asdf')
        [~,i] = min(vecnorm(V - ev.IntersectionPoint,2,2));
        
        % Modify palette color
        c = uisetcolor(P(i,:), 'Select color');
        src.CData(i,:) = c;

        % If color changed by reasonable amount, reconstruct image
        if (norm(c - P(i,:)) > 1e-8)
            P(i,:) = c;
            I = W_RGBXY * W_RGB * P;
            I = reshape(I, [width, height, 3]);
            plt_rec.CData = I;
        end

    end

end