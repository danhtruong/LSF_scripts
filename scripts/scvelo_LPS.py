
#import library
import scvelo as scv

#set parameters
scv.set_figure_params()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

#load data
adata_merged = scv.read("/rsrch3/scratch/sarc_med_onco-rsch/dtruong4/LPS_scRNA/LPS_data_cell_velocyto.h5ad")

print('Running recover_dynamics')
scv.tl.recover_dynamics(adata_merged, n_jobs = 20)

print('Running velocity')
scv.tl.velocity(adata_merged, mode='dynamical')

print('Running velocity_graph')
scv.tl.velocity_graph(adata_merged, n_jobs = 20)

print('Writing data')
adata_merged.write("/rsrch3/scratch/sarc_med_onco-rsch/dtruong4/LPS_scRNA/LPS_data_cell_velocyto_calc_dynamic.h5ad")

print('Job done')