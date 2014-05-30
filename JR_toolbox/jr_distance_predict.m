function d = jr_distance_predict(X, f)
% f = jr_distance_predict(X, f)

%% Prepare for calculations & error-check
[nr, nc] = size(X);                  % get dimentions of the data set
np       = size(f,3);

%% Predict glm values
d = NaN(nr,nc,np);
for c=1:nc                           % for each column representing a feature
    for i=1:np                       % go through all permutations of columns in d
        d(:,c,i) = glmval(f(:,c,i), X(:,c), 'probit');
    end
end
