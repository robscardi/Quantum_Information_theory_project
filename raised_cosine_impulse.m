function h = raised_cosine_impulse(t, beta, T)
    % raised_cosine: Computes the raised cosine impulse response.
    %
    % Inputs:
    %   t    - Time vector
    %   beta - Roll-off factor (0 <= beta <= 1)
    %   T    - Reciprocal of the symbol rate (symbol period)
    %
    % Output:
    %   h - Raised cosine impulse response

    % Initialize the impulse response vector
    h = zeros(size(t));

    % Calculate the raised cosine pulse
    for i = 1:length(t)
        if abs(1 - (2 * beta * t(i) / T)) < 1e-10 || abs(1 + 2*beta*t(i)/T) < 1e-10  % handle the singularity when t = T/(2*beta)
            h(i) = (pi/4) * sinc(1/(2*beta));
        elseif t(i) == 0
            h(i) = 1;
        else
            numerator = sin(pi * t(i) / T) .* cos(pi * beta * t(i) / T);
            denominator = (pi * t(i) / T) .* (1 - (2 * beta * t(i) / T).^2);
            h(i) = numerator / denominator;
        end
    end
end