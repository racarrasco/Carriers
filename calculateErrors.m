function[sigmap, correlation] = calculateErrors(tempSize, jacobian,rmse)
%%% Calculate the error on the fit parameters
%%%https://www.mathworks.com/help/pdf_doc/curvefit/curvefit.pdf
%%%page 7-52 through 7-53




%Unity weighting
weightones = ones(tempSize,1);
weights = diag(weightones);
%{
%weight the data by itself
weights = diag(t0_vTemp(1,:,end));
%}
hessian = jacobian.' *weights * jacobian;
covariance = inv(hessian).*rmse.^2;
sigmap = sqrt(diag(covariance));
correlation = zeros(size(covariance));
for i = 1:size(correlation,1)
    for j = 1:size(correlation,2)
        correlation(i, j) = covariance(i, j)./(sqrt(covariance(i, i))*...
            sqrt(covariance(j, j)));
    end
end

end