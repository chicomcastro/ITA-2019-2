function [isVisible, R_sat_e, V_sat_e_inercial, azimute, elevation] = isSVVisibleFrom(geodesicCoord, SV, t, Sat_OrbitData)
    [R_sat_e, V_sat_e_inercial] = coordSV(SV,t,Sat_OrbitData);
    [azimute, elevation] = azimuthElevation(geodesicCoord(1), ...
        geodesicCoord(2), ...
        geodesicCoord(3), ...
        R_sat_e);
    
    if (elevation > 0)
        isVisible = true;
    else
        isVisible = false;
    end
end

