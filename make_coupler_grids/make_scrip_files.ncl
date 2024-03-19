begin

a = addfile("grids.nc","r")
b = addfile("masks.nc","r")

um1t_lat = a->$"um1t.lat"$
um1t_lon = a->$"um1t.lon"$
um1t_lon_corners = a->$"um1t.clo"$
um1t_lat_corners = a->$"um1t.cla"$
um1t_mask = b->$"um1t.msk"$
um1t_mask = 1 - um1t_mask

um1t_lat_c_reorder = um1t_lat_corners(ny1|:, nx1|:, nc1|:)
um1t_lon_c_reorder = um1t_lon_corners(ny1|:, nx1|:, nc1|:)

opt = True
opt@GridCornerLat = um1t_lat_c_reorder
opt@GridCornerLon = um1t_lon_c_reorder
opt@GridMask = um1t_mask
opt@ForceOverwrite = True

um1t_file = "um1t_scrip.nc"
curvilinear_to_SCRIP(um1t_file, um1t_lat, um1t_lon, opt)
delete(opt)

um1u_lat = a->$"um1u.lat"$
um1u_lon = a->$"um1u.lon"$
um1u_lon_corners = a->$"um1u.clo"$
um1u_lat_corners = a->$"um1u.cla"$
um1u_mask = b->$"um1u.msk"$
um1u_mask = 1 - um1u_mask

um1u_lat_c_reorder = um1u_lat_corners(ny2|:, nx2|:, nc2|:)
um1u_lon_c_reorder = um1u_lon_corners(ny2|:, nx2|:, nc2|:)

opt = True
opt@GridCornerLat = um1u_lat_c_reorder
opt@GridCornerLon = um1u_lon_c_reorder
opt@GridMask = um1u_mask
opt@ForceOverwrite = True

um1u_file = "um1u_scrip.nc"
curvilinear_to_SCRIP(um1u_file, um1u_lat, um1u_lon, opt)
delete(opt)

um1v_lat = a->$"um1v.lat"$
um1v_lon = a->$"um1v.lon"$
um1v_lon_corners = a->$"um1v.clo"$
um1v_lat_corners = a->$"um1v.cla"$
um1v_mask = b->$"um1v.msk"$
um1v_mask = 1 - um1v_mask

um1v_lat_c_reorder = um1v_lat_corners(ny3|:, nx3|:, nc3|:)
um1v_lon_c_reorder = um1v_lon_corners(ny3|:, nx3|:, nc3|:)

opt = True
opt@GridCornerLat = um1v_lat_c_reorder
opt@GridCornerLon = um1v_lon_c_reorder
opt@GridMask = um1v_mask
opt@ForceOverwrite = True

um1v_file = "um1v_scrip.nc"
curvilinear_to_SCRIP(um1v_file, um1v_lat, um1v_lon, opt)
delete(opt)


cice_lat = a->$"cice.lat"$
cice_lon = a->$"cice.lon"$
cice_lon_corners = a->$"cice.clo"$
cice_lat_corners = a->$"cice.cla"$
cice_mask = b->$"cice.msk"$
cice_mask = 1 - cice_mask

cice_lat_c_reorder = cice_lat_corners(nyo|:, nxo|:, nco|:)
cice_lon_c_reorder = cice_lon_corners(nyo|:, nxo|:, nco|:)

opt1 = True
opt1@GridCornerLat = cice_lat_c_reorder
opt1@GridCornerLon = cice_lon_c_reorder
opt1@GridMask = cice_mask
opt1@ForceOverwrite = True

cice_file = "cice_scrip.nc"
curvilinear_to_SCRIP(cice_file, cice_lat, cice_lon, opt1)

end