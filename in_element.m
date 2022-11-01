function [found,I,J,L] = in_element(V,F,P,varargin)
  % IN_ELEMENT test for each p in P whether it lies inside each f in F defined
  % over V.
  % 
  % [found,I,J,L] = in_element(V,F,P)
  % [found,I,J,L] = in_element(V,F,P,'ParameterName',ParameterValue,...)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim+1 list of element indices
  %   P  #P by dim list of query positions
  %   Optional:
  %     'Method' followed by one of the following {'knn'}:
  %       'brute-force' no acceleration O(#P * #F)
  %       'edge-walk' walk along edges ~O(#P * sqrt(#F)) Starting with a random
  %         barycenter step along the edges that intersect with the ray toward
  %         the query point. If a boundary is reached then search along all
  %         boundary edges and jump to farthest hit. **dim=2 only**
  %       'knn' use knnsearch to find closest element barycenters
  %       'spatial-hash' spatial hash on regular grid ~O(#P * sqrt(#F)) **dim=2
  %         only**
  %       'ray-stabbing'
  %     'First'  only keep first match {false}
  %     'Quiet' suppress warnings {false}
  %     'Epsilon' epsilon used for determining inclusion {eps}
  % Outputs:
  %   found  #P by 1 list of flags whether point was found
  %   I  #I list of indices into P of identified finds
  %   J  #I list of indices into F of identified finds
  %   L  #I list of barycentric corodinates of P(I,:) in F(J,:)
  %

  function [I] = in_element_brute_force(V,F,P)
    dim = size(V,2);
    assert(dim+1 == size(F,2));
  
    % number of elements 
    m = size(F,1);
    % number of query points 
    np = size(P,1);
    
    switch dim 
    case 3
      T = F;
      % tet face ares
      vol = abs(volume(V,T));
      allF = [ ...
        T(:,2) T(:,4) T(:,3); ...
        T(:,1) T(:,3) T(:,4); ...
        T(:,1) T(:,4) T(:,2); ...
        T(:,1) T(:,2) T(:,3); ...
        ];
      % List of tets for each face f of each tet t for each point p
      TP = cat(2, ...
        repmat(allF,[1 1 np]), ...
        permute(repmat(size(V,1)+(1:np),m*4,1),[1 3 2]));
      TP = reshape(permute(TP,[1 3 2]),m*4*np,dim+1);
      Pvol = abs(volume([V;P],TP));
      % Pvol(t,f,p) --> volume of tet t, face f with point p
      Pvol = reshape(Pvol,[size(T,1) dim+1 np]);
      % sumvol(p,t) --> sum of volumes of tets made with faces of t and p
      sumvol = permute(sum(Pvol,2),[3 1 2]);
      I = sparse(abs(bsxfun(@minus,sumvol,vol')) < sqrt(epsilon));
    case 2
      % triangle side lengths
      l = [ ...
        sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
        sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
        sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
        ];
  
      B = zeros([np m dim+1]);
      for ii = 1:(dim+1)
        jj = mod(ii+1,dim+1)+1;
        kk = mod(ii,dim+1)+1;
        ljj = pdist2(P,V(F(:,jj),:));
        lkk = pdist2(P,V(F(:,kk),:));
        
        % semiperimeters
        s = bsxfun(@plus,l(:,ii)',ljj + lkk)*0.5;
        % Heron's formula for area
        B(:,:,ii) = 2*sqrt(s.*(bsxfun(@minus,s,l(:,ii)').*(s-ljj).*(s-lkk)));
      end
      % sum of barycentric coordinates
      sumA = sum(B,3);
      % area of element
      dblA = doublearea(V,F);
      %% check whether sum is more than true are
      %I = ~bsxfun(@gt,sumA,dblA');
      I = sparse((bsxfun(@minus,sumA,dblA')) < sqrt(epsilon));
    case 1
      % To-do: this is actually being computed in a dense way
      I = sparse( ...
        ((P <= V(F(:,1))') & (P >= V(F(:,2))')) | ...
        ((P <= V(F(:,2))') & (P >= V(F(:,1))')));
    end
    %B1 = sparse(B(:,:,1));
    %B2 = sparse(B(:,:,2));
    %B3 = sparse(B(:,:,3));
  end

  function I = in_element_hash_helper(V,F,P)
    assert(size(F,2) == 3, ...
      'F must contain triangles for Method=''spatial-hash''');
    num_bins = ceil(sqrt(size(F,1)));
    bin_x = ceil(sqrt(num_bins));
    bin_y = ceil(num_bins/bin_x); 
    num_bins = bin_x*bin_y;
    % spatial hash
    function VH = hash(V,MN,MX,bin_x,bin_y)
      [~,X] = histc(V(:,1),linspace(MN(1),MX(1),bin_x));
      [~,Y] = histc(V(:,2),linspace(MN(2),MX(2),bin_y));
      VH = sub2ind([bin_x bin_y],X,Y);
    end
    %% http://stackoverflow.com/a/5929567/148668
    %primes = [ 40960001, 59969537 45212177];
    %hash = @(X) mod( ...
    %  bitxor( int32(V(:,1)*primes(1)), int32(V(:,2)*primes(2))), ...
    %  num_bins);
    MN = min([V;P]);
    MX = max([V;P]);
    VH = hash(V,MN,MX,bin_x,bin_y);
    PH = hash(P,MN,MX,bin_x,bin_y);
    % This is wrong for triangles that span more hash cells than their vertices:
    % any time a triangle lands on the corner....
    FH = sparse(repmat(1:size(F,1),1,size(F,2))',VH(F(:)),1,size(F,1),num_bins)~=0;
    [~,FHI,FHJ] = find(FH);
    [FHJX,FHJY] = ind2sub([bin_x,bin_y],FHJ);
    % This assumes that #P >> #F
    I = sparse(size(P,1),size(F,1));
    for h = 1:num_bins
      %Vh = V(VH==h,:);
      IFh = FH(:,h);
      if any(IFh)
        IPh = PH==h;
        Ph = P(IPh,:);
        if any(IPh)
          Fh = F(IFh,:);
          Ih = in_element_brute_force(V,Fh,Ph);
          I(IPh,IFh) = Ih;
        end
      end
    end
  end

  function [I,J,L] = ray_stab_recursive(U,G,P)
    [J,~,L] = ray_mesh_intersect( ...
      [P ones(size(P,1),1)],repmat([0 0 -1],size(P,1),1),U,G);
    hit = J~=0;
    I = find(hit);
    J = J(hit);
    L = L(hit,:);
    if ~isempty(I) 
      % only keep triangles without a hit
      KG = 0==accumarray(J,1,[size(G,1) 1]);
      U = U(repmat(KG,3,1),:);
      G = reshape(1:size(U,1),[],3);
      P = P(I,:);
      [Ir,Jr,Lr] = ray_stab_recursive(U,G,P);
      K = find(KG);
      I = [I;I(Ir)];
      J = [J;K(Jr)];
      L = [L;Lr];
    end
  end

  % default values
  method = [];
  L = [];
  first = false;
  quiet = false;
  epsilon = eps;
  polish = eps;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
  {'Method','First','Quiet','Epsilon','Polish'}, ...
    {'method','first','quiet','epsilon','polish'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  % assumes samples are all inside exactly one element
  BC = barycenter(V,F);

  % indices of points we haven't found yet
  IP = 1:size(P,1);
  I = [];
  J = [];

  prev_k = 0;
  k = 2;
  while true
    K = knnsearch(BC,P(IP,:),'K',k);
    K = K(:,prev_k+1:end);
    for ki = 1:size(K,2)
      switch size(F,2)
      case 3
        B = abs(barycentric_coordinates( ...
          P(IP,:), ...
          V(F(K(:,ki),1),:), ...
          V(F(K(:,ki),2),:), ...
          V(F(K(:,ki),3),:)));
      case 4
        B = abs(barycentric_coordinates( ...
          P(IP,:), ...
          V(F(K(:,ki),1),:), ...
          V(F(K(:,ki),2),:), ...
        V(F(K(:,ki),3),:), ...
        V(F(K(:,ki),4),:)));
      case 6
        B = abs(barycentric_coordinates( ...
          P(IP,:), ...
          V(F(K(:,ki),1),:), ...
          V(F(K(:,ki),2),:), ...
          V(F(K(:,ki),3),:), ...
          V(F(K(:,ki),4),:), ...
          V(F(K(:,ki),5),:), ...
          V(F(K(:,ki),6),:)));
      end
       % asdf= abs(sum(B,2)-1);
      found = abs(sum(B,2)-1)<sqrt(epsilon);
        %I(sub2ind(size(I),IP,K(:,ki)')) = found;
      I = [I;IP(found)'];
      J = [J;K(found,ki)];

        % DOESN'T REALLY PAY OFF TO COMPUTE THESE HERE
        %B1(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,1);
        %B2(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,2);
        %B3(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,3);
        %if size(F,2) == 4
        %  B4(sub2ind(size(I),IP,K(:,ki)')) = found.*B(:,4);
        %end
        % Peel off found
      IP = IP(~found);
      if isempty(IP)
        break;
      end
      K = K(~found,:);
    end

      %for p = 1:numel(IP)
      %  Kp = K(p,:);
      %  U = V(F(Kp,:),:);
      %  FKp = reshape((1:size(U,1))',size(U,1)/size(F,2),size(F,2));
      %  I(IP(p),Kp) = in_element_brute_force(U,FKp,P(IP(p),:));
      %end
      %IP = find(~any(I,2));

    if isempty(IP)
      break;
    end
    prev_k = k;
    k = min(prev_k*2,size(BC,1));
    if k == size(BC,1)
      if ~quiet
        warning('Some points not found');
      end
      break;
    end
  end

    %return;

  assert(first == false,'not (re)implemented')
  %if first
  %  % Only keep first
  %  [mI,J] = max(I,[],2);
  %  I = sparse(1:size(I,1),J,mI,size(I,1),size(I,2));
  %end

  vec = @(X) X(:);
  I = vec(I);
  J = vec(J);
  found = accumarray(I,1,[size(P,1),1])>0;

  % Compute barycentric coordinates
  if nargout > 3 && (isempty(L) || polish)
    corners = arrayfun(@(i) V(F(J,i),:),1:size(F,2),'UniformOutput',false);
    L = barycentric_coordinates(P(I,:),corners{:});
  end

end


