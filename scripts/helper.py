
import logging
logger = logging.getLogger(__name__)

#https://stackoverflow.com/questions/20833344/fix-invalid-polygon-in-shapely
#https://stackoverflow.com/questions/13062334/polygon-intersection-error-in-shapely-shapely-geos-topologicalerror-the-opera
#https://shapely.readthedocs.io/en/latest/manual.html#object.buffer
def clean_invalid_geometries(geometries):
    """Fix self-touching or self-crossing polygons; these seem to appear
due to numerical problems from writing and reading, since the geometries
are valid before being written in pypsa-eur/scripts/cluster_network.py"""
    for i,p in geometries.items():
        if not p.is_valid:
            logger.warning(f'Clustered region {i} had an invalid geometry, fixing using zero buffer.')
            geometries[i] = p.buffer(0)
