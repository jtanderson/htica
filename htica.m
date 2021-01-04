function [ orthogonalizer, Aest, orthSig ] = htica(X, varargin)
  %% HTICA The heavy-tailed ICA procedure
  %  Inputs:
  %    X: dim-by-samples where "dim" is the number of sensors
  %
  damp = true;
  algoirthm = 'pow3';
  orthmethod = 'covariance';
  only_orthogonalize = false;

  %% ----------------------- Process Arguments -------------------------

  if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
  else
    for i=1:2:(length(varargin)-1)
      if ~ischar (varargin{i}),
        error (['Unknown type of optional parameter name (parameter' ...
          ' names must be strings).']);
      end
      % change the value of parameter
      switch lower (varargin{i})
        case 'damp'
          damp = strcmpi(varargin{i+1}, 'true');
        case 'verbose'
          verbose = strcmpi(varargin{i+1}, 'true');
        case 'algorithm'
          algorithm = lower(varargin{i+1});
        case 'orthmethod'
          orthmethod = lower(varargin{i+1});
        case 'orthmatrix'
          orthmethod = 'oracle'
          orthogonalizer = varargin{i+1}
        case 'only'
          switch lower(varargin{i+1})
            case 'orthogonalize'
              only_orthogonalize = true;
          end
        otherwise
          error(['Unknown argument: ''' varargin{i} '''']);
      end;
    end;
  end

  %% --------------------------------------------------------------------

  [n,m] = size(X);

  switch orthmethod
    case 'covariance'
      C = (1/m) * (X * X');
      orthogonalizer = inv(sqrtm(C));
    case 'centroid'
      if verbose
        tic;
      end
      orthogonalizer = centroidOrthogonalizer(X);
      if verbose
        disp(['Time to orthogonalize via centroid body: ' num2str(toc)]);
      end
    case 'centroidscaling'
      if verbose
        tic;
      end
      orthogonalizer = centroidOrthogonalizer(X, 'scale');
      if verbose
        disp(['Time to orthogonalize via centroid body: ' num2str(toc)]);
      end
    case 'oracle'
      % The user gave it to us... hopefully
      if verbose
        disp(['Using provided orthogonalization matrix']);
      end
      if ~exist(orthogonalizer)
        error(['To use oracle orthogonalization, pass the orthogonalization matrix with the "orthmatrix" parameter'])
      end
    case 'none'
      orthogonalizer = eye(n);
    otherwise
      error(['Invalid orthogonalizer choice ''' algorithm '''']);
  end

  X = orthogonalizer * X;
  orthSig = X;
  Aest = eye(n); % just in case we skip the next step b/c of only_orthogonalize

  if ~only_orthogonalize
    % Perform damping, if necessary.
    if damp
      % Choose some constants to constrain the binary search.
      % Can be tweaked for different cases, see the paper for details.
      C2 = 3;
      R = 1;
      Kest = 0;
      cumest = 0;
      % Currently a bad idea to estimate K_{X_R} from the same samples
      % that we're going to use later, but can be fixed easily
      while Kest <= 0.5 || cumest <= 1/(n^C2)
        % Use a different Z every time
        Z = unifrnd(0,1,1,size(X,2));
        samplecount = size(S,2);
        % At termination, we have R large enough and already know the
        % values Exp[-Norm[x]^2/R^2 for each sample point x
        R = R*2;
        if R == Inf && verbose
          disp('Could not find large enough R...');
          disp(['Current Kest: ' num2str(Kest)]);
          disp(['Current cumest: ' num2str(cumest)]);
          error('Failed: R too large!');
        end
        Xthreshold = exp(-sum(X.^2,1)/R^2);
        Kest = mean(Xthreshold);
        Sthreshold = exp(-sum(X.^2,1)/R^2);
        tmp =  S(Z <= Sthreshold, :);
        cumest = min(abs(sum(tmp.^4, 1)/size(tmp,1) - 3*(sum(tmp.^2,1)/size(tmp,1))));
      end

      if verbose
        disp(['Chosen R: ' int2str(R)]);
      end

      % Reject samples based on the uniform samples z vs the damping
      firstcount = size(X,2);
      Z = unifrnd(0,1,1,size(X,2));
      X = X(:,Z <= Xthreshold);
      if verbose
        disp(['Samples remaining after rejection: ' ... 
          int2str(size(X,2)) ' out of ' int2str(firstcount) ...
          ' (' num2str(100*size(X,2)/firstcount) '%)']);
      end
    end

    % Finally, execute the desired ICA algorithm
    switch algorithm
      case 'pow3'
        % Run ICA with g(u) = u^3
        [~, Aest, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n);
        [whitesig, B, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n, 'only', 'white');

        % Sometimes the FastICA result has less columns than we
        % actually want. Try again, please...
        while size(Aest,2) < dim
          error('Degenerate!')
          size(Aest)
          Aest
          if verbose
            disp('Having to correct degeneracy...')
          end
          [~, Aest, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n);
          [~, B, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n, 'only', 'white');
        end
        disp(['Whitened center norm:' num2str(norm(mean(whitesig,2)))]);
      case 'tanh'
        % Run ICA with g(u) = tanh(u)
        [~, Aest, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n, 'g', 'tanh');

        % Sometimes the FastICA result has less columns that we
        % actually want. Try again, please...
        while size(Aest,2) < dim
          if verbose
            disp('Having to correct degeneracy...')
          end
          [~, Aest, ~] = fastica(X, 'verbose', 'off', 'numOfIC', n, 'g', 'tanh');
        end
      case 'fpca'
        % Run FPCA (naive for now). Note we need to transpose X
        V = phaseCorrect(naiveFPCA(X', 1.6));
        B = (X*X')/m;
        Aest = sqrtm(B)*V;
      case 'sobi'
        [Aest, ~] = sobi(X);
      case 'jade'
        Aest = inv(jadeR(X));
      case 'simple'
        cov = (X*X')/m;
        [V, ~] = eig(cov);
        Aest = V;
      case 'yeredor'
        [~,m1]=size(X);
        X1 = X(:,1:floor(m1/2));                
        X2 = X(:,floor(m1/2)+1:m1);

        [~,col1] = size(X1);
        [~,col2] = size(X2);
        M1 = (X1*X1')/col1;
        M2 = (X2*X2')/col2;

        M = M1*inv(M2);
        [V, ~] = eig(M);
        Aest = V;
      otherwise
        error(['Invalid algorithm choice ''' algorithm '''']);
    end

    Aest = inv(orthogonalizer) * Aest;

    % Normalize the results
    Aest = Aest*(inv(diag(rownorm(Aest'))));
  end
end
