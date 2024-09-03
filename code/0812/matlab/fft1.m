function frequencies = fftfreq(n, d)
    % Calculate the frequencies for the FFT output.
    % n: Number of data points
    % d: Sample spacing (inverse of sampling frequency)

    % Ensure that n is even (if it's not, round up to the nearest even number).
    if mod(n, 2) == 1
        n = n + 1;
    end

    % Calculate the spacing between frequencies.
    df = 1 / (n * d);

    % Generate the frequency array.
    frequencies = linspace(0, 1, n+1); % MATLAB's linspace includes both endpoints
    frequencies = frequencies(1:n) - 0.5; % Center frequencies around 0.0
    frequencies = frequencies / df;
end