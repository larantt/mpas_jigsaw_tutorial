import numpy as np

# ---------------------------------------------------------------------------
# Region definition -- these are the knobs you tune.
# ---------------------------------------------------------------------------

hfun_min = 3.0     # Grid distance (km) in refinement region
hfun_max = 60.0    # Grid distance (km) outside refinement region

# Corners of the region of interest, (lon, lat) in degrees.
# The polygon through these points is inflated by halo_km, which also
# rounds off the corners. Add as many anchors as you like.
anchors = [
    (121.42, -17.1),   
    (151.2, -7.3),   
    (162.4, -27.9),
    (142.3, -39.45)   
]

halo_km = 100.0 #800.0    # The h_min zone extends this far beyond the polygon
trans_km = 1200.0  # Sets the OUTER slope: (hfun_max - hfun_min) / trans_km

# The transition is in two legs with different slopes, meeting at hfun_mid.
# Inner leg: hfun_min -> hfun_mid over inner_km, deliberately shallow so the
# finest resolutions change slowly. Outer leg: hfun_mid -> hfun_max at the
# original slope. The two legs meet at hfun_mid, so h(r) stays continuous.
hfun_mid = 12.0    # Grid distance (km) at which the slope changes
inner_km = 1000.0  # Width (km) of the shallow hfun_min -> hfun_mid leg

# Anywhere further than this from the region centre is forced to h_max.
# This guarantees a perfectly uniform background over the rest of the globe.
far_field_km = 15000.0

# TC track to overlay on the plot, for checking coverage. None to skip.
TC_DICT = {'name': 'frances', 'basin': 'all', 'year': 2004}

_field = None      # cached (interpolator) -- built once, on first use


# ---------------------------------------------------------------------------
# The mesh-size function
# ---------------------------------------------------------------------------

def get_hfun(longitude, latitude):
    # Distance (km) to the region: negative inside, positive outside
    r = region_distance(np.degrees(longitude), np.degrees(latitude))

    # Return grid distances
    return h(r)


def h(r):
   t_begin = 0.0            # Ramp begins at the region edge (r = 0)
   t_mid = inner_km         # Shallow inner leg ends here (h = h_mid)
   t_end = t_end_km()       # Outer leg ends here (h = h_max)
   h_min = hfun_min         # Grid distance (km) in refinement region
   h_mid = hfun_mid         # Grid distance (km) where the slope changes
   h_max = hfun_max         # Grid distance (km) outside region and transition

   hires_mask = r < t_begin
   inner_mask = np.logical_and(r >= t_begin, r < t_mid)
   transition_mask = np.logical_and(r >= t_mid, r < t_end)
   lowres_mask = r >= t_end

   ret = np.zeros_like(r)
   ret[hires_mask] = h_min
   ret[inner_mask] = h_min + (r[inner_mask] - t_begin) * (h_mid - h_min) / (t_mid - t_begin)
   ret[transition_mask] = h_mid + (r[transition_mask] - t_mid) * outer_slope()
   ret[lowres_mask] = h_max

   return ret


def outer_slope():
    """Slope (km per km) of the outer leg -- unchanged from the original."""
    return (hfun_max - hfun_min) / trans_km


def t_end_km():
    """Distance (km) at which the outer leg reaches h_max."""
    return inner_km + (hfun_max - hfun_mid) / outer_slope()


# ---------------------------------------------------------------------------
# Distance to the region
#
# Everything happens in a map projection centred on the region, where
# distances are in metres and shapely can be used directly. Because the
# projection is flat, there is no "other side of the globe" to worry about.
# ---------------------------------------------------------------------------

def region_distance(lon_deg, lat_deg):
    """Distance (km) to the region boundary. Negative inside, positive out."""
    global _field
    if _field is None:
        _field = build_field()

    # Look up the precomputed field. Longitudes must be in [-180, 180).
    lon = ((np.asarray(lon_deg) + 180.0) % 360.0) - 180.0
    lat = np.asarray(lat_deg)
    return _field((lat.flatten(), lon.flatten()))


def build_field():
    """Precompute distance-to-region on a lon/lat grid, once."""
    from shapely.geometry import Polygon, Point
    from scipy.interpolate import RegularGridInterpolator
    import pyproj

    # 1. Set up a projection centred on the middle of the anchors.
    #    Azimuthal equidistant: distances from the centre are true.
    lon0 = float(np.mean([a[0] for a in anchors]))
    lat0 = float(np.mean([a[1] for a in anchors]))
    proj = pyproj.Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m +datum=WGS84",
        always_xy=True,
    )

    # 2. Build the region: the anchor polygon, inflated by the halo.
    #    buffer() rounds the corners for us.
    corners_x, corners_y = proj.transform(
        [a[0] for a in anchors], [a[1] for a in anchors]
    )
    region = Polygon(zip(corners_x, corners_y)).buffer(halo_km * 1000.0)

    # 3. Make a lon/lat grid to evaluate the distance on.
    grid_res = 0.25
    lats = np.arange(-90.0, 90.0 + grid_res, grid_res)
    lons = np.arange(-180.0, 180.0 + grid_res, grid_res)
    grid_lon, grid_lat = np.meshgrid(lons, lats)
    grid_x, grid_y = proj.transform(grid_lon, grid_lat)

    # 4. Distance from each grid point to the region boundary,
    #    negated where the point is inside the region.
    dist = np.full(grid_lon.shape, far_field_km)

    for j in range(grid_lon.shape[0]):
        for i in range(grid_lon.shape[1]):
            x, y = grid_x[j, i], grid_y[j, i]

            # Points too far away (or that don't project) are background.
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            if np.hypot(x, y) > far_field_km * 1000.0:
                continue

            point = Point(x, y)
            d_m = region.exterior.distance(point)
            if region.contains(point):
                d_m = -d_m
            dist[j, i] = d_m / 1000.0          # metres -> km

    # 5. Wrap it up so we can look up any (lat, lon) later.
    return RegularGridInterpolator(
        (lats, lons), dist, bounds_error=False, fill_value=far_field_km
    )


# ---------------------------------------------------------------------------
# Kept from the original script
# ---------------------------------------------------------------------------

def geo_to_cartesian(lam, phi):
    x = np.cos(lam) * np.cos(phi)
    y = np.sin(lam) * np.cos(phi)
    z = np.sin(phi)

    return (x, y, z)


def unit_sphere_distance(p, q_arr):
    return np.arccos(q_arr @ p)


if __name__ == "__main__":
    from shapely.geometry import Polygon
    from matplotlib.backends.backend_pdf import PdfPages
    from tropycal import tracks, utils
    import cartopy
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    import matplotlib

    res = 0.5   # Resolution (degrees) of the plotting grid

    # Set to e.g. [90, 160, -5, 50] to zoom page 1 on the region;
    # None keeps the global view.
    ZOOM = None

    pdf = PdfPages('hfun.pdf')

    # ---------------- Page 1: map of h(x) ----------------------------------
    fig = plt.figure()
    proj = cartopy.crs.PlateCarree(0.0)
    ax = plt.axes(projection=proj)
    if ZOOM is None:
        ax.set_global()
    else:
        ax.set_extent(ZOOM, crs=cartopy.crs.PlateCarree())

    cmap = matplotlib.colormaps.get_cmap('Blues_r')

    cells = []
    for lat in np.arange(-90.0, 90.0, res):
        for lon in np.arange(-180.0, 180.0, res):
            lon_corners = [lon, lon + res, lon + res, lon]
            lat_corners = [lat, lat, lat + res, lat + res]
            cells.append(Polygon(np.array(list(zip(lon_corners,lat_corners)))))

    latc = np.arange(-90.0 + 0.5*res, 90.0, res)
    lonc = np.arange(-180.0 + 0.5*res, 180.0, res)
    lons, lats = np.meshgrid(lonc, latc)
    h_flat = get_hfun(np.deg2rad(lons), np.deg2rad(lats)).flatten()

    norm = plt.Normalize(vmin=np.min(h_flat), vmax=np.max(h_flat))

    ax.add_geometries(cells, cartopy.crs.PlateCarree(), facecolors=cmap(norm(h_flat)), edgecolors=cmap(norm(h_flat)), linewidth=0.03)

    gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
        linewidth=0.3, color='black', alpha=1.0, linestyle='--')

    gl.top_labels = False
    gl.left_labels = False

    xticks = np.arange(-180, 180, 15)
    yticks = np.arange(-90, 90, 15)

    gl.ylocator = mticker.FixedLocator(yticks)
    gl.xlocator = mticker.FixedLocator(xticks)

    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5, edgecolor='black')

    # Mark the anchor points defining the region
    ax.plot([a[0] for a in anchors], [a[1] for a in anchors],
            marker='x', color='k', markersize=6, mew=1.5, linestyle='none',
            transform=cartopy.crs.PlateCarree(), zorder=99)

    # Overlay the TC track to check the region covers it.
    # NB: add_tropycal() rebinds ax, so this must come after the cartopy calls.
    if TC_DICT is not None:
        print(f"getting track for {TC_DICT['name']}, {TC_DICT['year']}")
        basin = tracks.TrackDataset(basin=TC_DICT['basin'], include_btk=False,
                                    source='ibtracs')
        storm = basin.get_storm((TC_DICT['name'], TC_DICT['year']))
        print(storm)

        ax = utils.add_tropycal(ax)
        ax.plot_storm(storm, color='red', linewidth=1.5, zorder=100)

        ax.set_title(f"TC {TC_DICT['name'].title()}, {TC_DICT['year']} (IBTrACS)",
                     loc='right', fontweight='bold', fontsize=10)

    pdf.savefig(fig, bbox_inches='tight', transparent=True)
    plt.close(fig)

    # ---------------- Page 2: the transition function h(r) -----------------
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))

    t_end = t_end_km()
    r = np.linspace(-halo_km - 200.0, t_end + 500.0, 2000)
    hr = h(r)

    ax2.plot(r, hr, color='tab:blue', linewidth=2)

    ax2.axvline(0.0, color='red', linestyle='--', linewidth=1)
    ax2.axvline(inner_km, color='0.4', linestyle=':', linewidth=1)
    ax2.axvline(t_end, color='0.4', linestyle='--', linewidth=1)

    ax2.annotate('region edge\n(h_min zone ends)', xy=(0.0, hfun_min),
                 xytext=(-halo_km * 0.9, hfun_min + 0.45 * (hfun_max - hfun_min)),
                 fontsize=8, color='red',
                 arrowprops=dict(arrowstyle='->', color='red', lw=0.8))
    ax2.annotate(f'slope change\nh = {hfun_mid:g} km', xy=(inner_km, hfun_mid),
                 xytext=(inner_km * 0.55, hfun_mid + 0.30 * (hfun_max - hfun_min)),
                 fontsize=8, color='0.3',
                 arrowprops=dict(arrowstyle='->', color='0.3', lw=0.8))
    ax2.annotate('background h_max', xy=(t_end, hfun_max),
                 xytext=(t_end * 0.60, hfun_max - 0.18 * (hfun_max - hfun_min)),
                 fontsize=8, color='0.3',
                 arrowprops=dict(arrowstyle='->', color='0.3', lw=0.8))

    ax2.fill_between(r, hfun_min, hr, where=(r < 0.0), color='tab:blue', alpha=0.12)
    ax2.text(-halo_km * 0.5, hfun_min + 0.06 * (hfun_max - hfun_min),
             f'inside region\nh = {hfun_min:g} km', fontsize=8, ha='center', color='tab:blue')
    ax2.text(inner_km * 0.5, hfun_min + 0.16 * (hfun_max - hfun_min),
             f'inner leg\n{inner_km:g} km', fontsize=8, ha='center', color='tab:blue')
    ax2.text((inner_km + t_end) * 0.5, hfun_min + 0.72 * (hfun_max - hfun_min),
             f'outer leg\n{t_end - inner_km:.0f} km', fontsize=8, ha='center', color='tab:blue')

    ax2.set_xlabel('Distance to region boundary (km)   [negative = inside]')
    ax2.set_ylabel('Target grid distance h (km)')
    ax2.set_title(f'Mesh transition function: {hfun_min:g} km core -> {hfun_max:g} km background\n'
                  f'halo = {halo_km:g} km, inner leg {hfun_min:g}-{hfun_mid:g} km over {inner_km:g} km, '
                  f'then {outer_slope()*1000.0:.0f} m/km to {hfun_max:g} km', fontsize=9)
    ax2.grid(True, linewidth=0.3, alpha=0.5)
    ax2.set_ylim(0, hfun_max + 5.0)

    pdf.savefig(fig2, bbox_inches='tight')
    plt.close(fig2)

    pdf.close()

    print('Saving HFUN plot to hfun.pdf')