function assertErrorId(fcn, expectedId)
%ASSERTERRORID Assert that a function handle throws with a given identifier.

try
    fcn();
catch ME
    if ~strcmp(ME.identifier, expectedId)
        error('assertErrorId:wrongId', ...
            'Expected id ''%s'', got ''%s'': %s', ...
            expectedId, ME.identifier, ME.message);
    end
    return;
end

error('assertErrorId:noError', 'Expected error id ''%s'' but call succeeded.', expectedId);
end
