function grid = snapgrid(spacing, geodata, param, first,grid)

% input:
%   spacing - lat/lon grid spacing, e.g. 0.25, 0.5, 1.0, etc.
%   geodata.lat, geodata.lon, geodata.satzen, geodata.sunzen
%   param - the parameter you want to snap to grid, 'rad'for radiance or
%           'cld' for cloud
%   first - 1=first; 0=not
%   grid  - regrid out put from a previous run. 
%   type  - 1=day; 0=night
% output:
%   regrid - data structure with regridded fields
%
% Nadia Smith, April 2010
% Variation on Nick Bearson's code
% =======================================================

% if this is the first granule/file, then create new grids
if first==1
    % day time
    grid.d_dist   = ones(180/spacing, 360/spacing)*NaN;
    grid.d_satzen = ones(180/spacing, 360/spacing)*NaN;
    grid.d_solzen = ones(180/spacing, 360/spacing)*NaN;
    grid.d_param  = ones(180/spacing, 360/spacing)*NaN;
    grid.d_type   = ones(180/spacing, 360/spacing)*NaN;
    % night time
    grid.n_dist   = ones(180/spacing, 360/spacing)*NaN;
    grid.n_satzen = ones(180/spacing, 360/spacing)*NaN;
    grid.n_solzen = ones(180/spacing, 360/spacing)*NaN;
    grid.n_param  = ones(180/spacing, 360/spacing)*NaN;
    grid.n_type   = ones(180/spacing, 360/spacing)*NaN;
end

for i=1:length(geodata.lat(:,1))
    for j=1:length(geodata.lat(1,:))

      % get the lat/lon of the data point
      lat = geodata.lat(i, j);
      lon = geodata.lon(i, j);

      if isnan(lat) == 0 && isnan(lon) == 0

        % adjust the lon/lat to be 0-360 and 0-180, and then
        % adjust them to our spacing grid for easier calculation
        adjlat = ((lat + 90) / spacing);
        adjlon = ((lon + 180) / spacing);

        % now find the closest grid crossing
        latbin = round(adjlat);
        lonbin = round(adjlon);

        % We don't want to reference the 0 spot, so we'll use the 180/360 instead (lat/lon)
        if latbin == 0
          latbin = 180 / spacing;
        end
        if lonbin == 0
          lonbin = 360 / spacing;
        end

        % determine our point's distance from the grid crossing
        % note that we're finding the hypotenuse distance, as that differs from just adding the two
        dist = sqrt(((latbin - adjlat)^2) + ((lonbin - adjlon)^2));
        if geodata.sunzen(i, j) < 84 % day time grid
            % now we're ready to get the current distance in that bin and see if ours is closer
            if (dist < (1/spacing)/4) && (isnan(grid.d_satzen(latbin, lonbin)) == 1 ||...
                    abs(geodata.satzen(i, j)) < abs(grid.d_satzen(latbin, lonbin)))
                grid.d_dist(latbin, lonbin) = dist;
                grid.d_satzen(latbin, lonbin) = geodata.satzen(i, j);
                grid.d_sunzen(latbin, lonbin) = geodata.sunzen(i, j);
                grid.d_param(latbin, lonbin) = param(i, j);  
            end
        else % night time grid
            % now we're ready to get the current distance in that bin and see if ours is closer
            if (dist < (1/spacing)/4) && (isnan(grid.n_satzen(latbin, lonbin)) == 1 ||...
                    abs(geodata.satzen(i, j)) < abs(grid.n_satzen(latbin, lonbin)))
                grid.n_dist(latbin, lonbin) = dist;
                grid.n_satzen(latbin, lonbin) = geodata.satzen(i, j);
                grid.n_sunzen(latbin, lonbin) = geodata.sunzen(i, j);
                grid.n_param(latbin, lonbin) = param(i, j);  
            end
        end
     end
  end
end



