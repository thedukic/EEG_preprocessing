function check_duplicatefunc(functionName)
% CHECK_DUPLICATEFUNC Checks for multiple instances of a function on the MATLAB path.
%   Throws an error if multiple instances are found.

  paths = which('-all', functionName);
  
  if numel(paths) > 1
      errorMessage = ['Multiple instances of function "' functionName '" found on the MATLAB path:' newline];
      for i = 1:numel(paths)
          errorMessage = [errorMessage, paths{i}, newline]; 
      end

      % Add more info
      errorMessage = [errorMessage, 'Please, first, reset your MATLAB paths to default.', newline];
      % Open the path window
      pathtool;
      % Throw the error
      error(errorMessage); 
  end
end