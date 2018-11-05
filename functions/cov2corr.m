function [ExpSigma, ExpCorrC] = cov2corr(ExpCovariance)
%COV2CORR Converts covariance to standard deviation and correlation coefficient.
%   Computes the volatilities of N random processes and the degree of
%   correlation between the processes.  
%
%   [ExpSigma, ExpCorrC] = cov2corr(ExpCovariance)
%
%   Input:
%     ExpCovariance: N by N covariance matrix, e.g. from COV or EWSTATS
%
%   Outputs:
%     ExpSigma     : 1 by N vector with the standard deviations of each process
%
%     ExpCorrC     : N by N matrix of correlation coefficients.  The
%     entries of ExpCorrC range from 1 (completely correlated) to -1
%     (completely anti-correlated).  A value of 0 in the (i,j) entry
%     indicates that the i'th and j'th processes are uncorrelated.
% 
%   ExpSigma(i) = sqrt( ExpCovariance(i,i) );
%   ExpCorrC(i,j) = ExpCovariance(i,j)/( ExpSigma(i)*ExpSigma(j) );
% 
%   See also EWSTATS, COV, CORRCOEF, STD, CORR2COV.

%    Copyright 1995-2006 The MathWorks, Inc.


%-----------------------------------------------------------------
% Argument checking
% ExpCovariance   [N by N]  with diag(ExpCovariance)>=0
% N      [scalar]
%-----------------------------------------------------------------
if nargin<1,
  error(message('finance:cov2corr:missingInput'))
end

if size(ExpCovariance,1)~=size(ExpCovariance,2)
  error(message('finance:cov2corr:invalidCovMatrixSize'))
else
  N = size(ExpCovariance,1);
end

if any( diag(ExpCovariance)<0 )
  error(message('finance:cov2corr:invalidCovMatrixSymmetry'))
end

%-----------------------------------------------------------------
% Simple correlation is ExpCovariance./( ExpSigma'*ExpSigma )
% ExpSigma [1 by N]
% ExpCorrC [N by N]
%-----------------------------------------------------------------
ExpSigma = sqrt(diag(ExpCovariance))'; 

% start with default correlation of identity for degenerate processes
ExpCorrC = eye(N);

% find processes which are not degenerate
IndPos = ExpSigma>0;

% Compute correlation only for non-degenerate processes
ExpCorrC(IndPos,IndPos) = ExpCovariance(IndPos,IndPos) ./ ... 
    (ExpSigma(IndPos)'*ExpSigma(IndPos));

% Force exact ones along the main diagonal.
ExpCorrC(sub2ind([N N], 1:N, 1:N)) = 1;

%-----------------------------------------------------------------
% end of function COV2CORR
%-----------------------------------------------------------------
