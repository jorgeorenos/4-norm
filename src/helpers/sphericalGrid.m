function [azGrid, polGrid, Xgrid, amplitudes] = sphericalGrid(Q)
%SPHERICALGRID Build a spherical grid, its 4-norm normalization, and amplitudes.
    azimuthCount = 180;
    polarCount = 90;
    azimuth = (0:azimuthCount) * (2*pi/azimuthCount);
    polar = (0:polarCount) * (pi/polarCount);
    [azGrid, polGrid] = meshgrid(azimuth, polar);
    sinPol = sin(polGrid);
    Xgrid = [reshape(cos(azGrid).*sinPol, 1, []); ...
        reshape(sin(azGrid).*sinPol, 1, []); ...
        reshape(cos(polGrid), 1, [])];
    Xgrid = Xgrid ./ sum(abs(Xgrid).^4, 1).^(1/4);
    amplitudes = reshape(sum(abs(Q * Xgrid).^4, 1).^(1/4), size(azGrid));
end
