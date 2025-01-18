function h = root_raised_cosine_impulse(t, beta, T)
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

    c = (1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/4*beta); 
    % Calculate the raised cosine pulse
    for i = 1:length(t)
        if abs(1 - (4 * beta * t(i) / T)) < 1e-10 || abs(1 + 4*beta*t(i)/T) < 1e-10  % handle the singularity when t = T/(2*beta)
            h(i) = (beta/(T*sqrt(2)))*c;
        elseif t(i) == 0
            h(i) = (1+beta*((4/pi)-1))/T;
        else
            numerator = sin(pi*(t(i)/T)*(1-beta)) + (4*beta*(t(i)/T)*cos(pi*(t(i)/T)*(1+beta)));
            denominator = pi*(t(i)/T)*(1-(4*beta*(t(i)/T)^2));
            h(i) =(1/T)* numerator / denominator;
        end
    end
end