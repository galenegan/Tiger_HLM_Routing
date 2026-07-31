import xarray as xr
import pandas as pd
import numpy as np
ds = xr.open_dataset('/Users/ea-gegan/Documents/gitrepos/Tiger_HLM_Routing/data/runoff/total_runoff_2017.nc')


valid_time = pd.to_datetime(pd.date_range("2026-01-01", "2026-01-04", freq="H")).values
ro = np.zeros((len(valid_time),))
ro[12:25] = np.linspace(0, 10, 13)
ro[24:38] = np.linspace(10, 5, 14)
ro[38:] = 5
ro = ro.reshape(1, -1)

link_id = np.array([420558772])

ds_out = xr.Dataset(
    data_vars={
        "ro": (["LinkID", "valid_time"], ro)
    },
    coords={
        "LinkID": link_id,
        "valid_time": valid_time
    }
)
ds_out.to_netcdf('test_runoff.nc')