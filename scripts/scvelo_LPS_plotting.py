import os
cwd = os.getcwd()
print(cwd)
#import library
import scvelo as scv

#set parameters
scv.set_figure_params()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

#load data
adata_merged = scv.read("/rsrch3/scratch/sarc_med_onco-rsch/dtruong4/LPS_scRNA/LPS_data_cell_velocyto_calc_dynamic.h5ad")

print('Plotting velocity_embedding_stream')
#scv.pl.velocity_embedding_stream(adata_merged,
# basis='umap',
#  color = 'Annotations',
# save = 'lps_velocity_embedding_stream_dynamic.png')


#Identify important genes
#print('Identify important genes')
#scv.tl.rank_velocity_genes(adata_merged, groupby='Annotations', min_corr=.3)

#df = scv.DataFrame(adata_merged.uns['rank_velocity_genes']['names'])
#df.to_csv('LPS_rank_velocity_genes.csv')

scv.tl.latent_time(adata_merged)
scv.pl.scatter(adata_merged,
 color='latent_time', 
color_map='gnuplot', 
save = 'lps_velocity_embedding_stream_dynamic_latent_time.png')

print('Job done')