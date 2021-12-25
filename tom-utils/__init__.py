import intake
from ophyd import Signal

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

# Look up a driver class by its name in registry.
catalog_class = intake.registry["bluesky-mongo-normalized-catalog"]

sirepo_bluesky_catalog_instance = catalog_class(
        metadatastore_db="mongodb://localhost:27017/md",
        asset_registry_db="mongodb://localhost:27017/ar",
)



def rot_2d(a):
    return np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])

def xy_to_ae(dx,dy,ca=0,ce=0):

y,z = np.matmul()



return a,e