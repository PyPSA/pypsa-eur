import pandas as pd
import geopandas as gpd

def area(gdf):
    """Returns area of GeoDataFrame geometries in square kilometers."""
    return gdf.to_crs(epsg=3035).area.div(1e6)


def allocate_sequestration_potential(gdf, regions, attr='conservative estimate Mt', threshold=3):
    gdf = gdf.loc[gdf[attr] > threshold, [attr, "geometry"]]
    gdf["area_sqkm"] = area(gdf)
    overlay = gpd.overlay(regions, gdf, keep_geom_type=True)
    overlay["share"] = area(overlay) / overlay["area_sqkm"]
    adjust_cols = overlay.columns.difference({"name", "area_sqkm", "geometry", "share"})
    overlay[adjust_cols] = overlay[adjust_cols].multiply(overlay["share"], axis=0)
    gdf_regions = overlay.groupby("name").sum()
    gdf_regions.drop(["area_sqkm", "share"], axis=1, inplace=True)
    return gdf_regions.squeeze()


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_sequestration_potentials',
            simpl='',
            clusters="181"
        )

    # TODO move to config.yaml
    threshold = 3
    include_onshore = False

    gdf = gpd.read_file(snakemake.input.sequestration_potential)

    regions = gpd.read_file(snakemake.input.regions_offshore)
    if include_onshore:
        onregions = gpd.read_file(snakemake.input.regions_onshore)
        regions = pd.concat([regions, onregions]).dissolve(by='name').reset_index()

    attr = snakemake.config['sector']["sequestration_potential"]
    kwargs = dict(attr=attr, threshold=threshold) if isinstance(attr, str) else {}

    s = allocate_sequestration_potential(gdf, regions, **kwargs)

    s = s.where(s>threshold).dropna()

    s.to_csv(snakemake.output.sequestration_potential)
