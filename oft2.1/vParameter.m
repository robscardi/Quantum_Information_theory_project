function v = vParameter(diameter, wavelength, material1, material2)
% Returns the V-parameter for two-layer structures 
% 
% Diameter is in microns, wavelength is in
% nanometers. Materials are as understood by refrIndex.

v = diameter * pi ./ (wavelength / 1000) .* ...
    sqrt(refrIndex(material1, wavelength).^2 - ...
    refrIndex(material2, wavelength).^2);
