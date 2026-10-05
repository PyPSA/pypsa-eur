# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


rule build_shapes_parquet:
    input:
        shapes=resources("onshore_regions.geojson"),
    output:
        shapes_parquet=resources("onshore_regions.parquet"),
    run:
        import geopandas as gpd
        import country_converter as coco

        cc = coco.CountryConverter()
        gdf = gpd.read_file(input.shapes)
        gdf = gdf.rename({"name": "shape_id"}, axis=1)
        iso2 = gdf["shape_id"].str[:2]
        gdf["country_id"] = cc.convert(names=iso2, src="ISO2", to="ISO3")
        gdf["shape_class"] = "land"
        gdf.to_parquet(output.shapes_parquet, index=False)


module module_euro_gas_grid:
    pathvars:
        logs="logs/module_euro_gas_grid",
        resources="resources/module_euro_gas_grid",
        results="resources/module_euro_gas_grid",
        user_shapes="resources/{shapes}/onshore_regions.parquet",
    snakefile:
        github(
            "modelblocks-org/module_euro_gas_grid",
            path="workflow/Snakefile",
            branch="v0.1.3",
        )
    config:
        config["euro_gas_grid"]


use rule * from module_euro_gas_grid as module_euro_gas_grid_*
