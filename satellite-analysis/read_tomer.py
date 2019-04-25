import pandas as pd

def read_tomer(galaxynumber):
    global tomer_list
    path = r'/Users/user1/Documents/rockstar_analysis/tomercatalogs/'

    cen_gal_cat = pd.read_pickle(path + 'cen_gal_cat.pkl')
    sim_cat     = pd.read_pickle(path + 'sim_table.pkl')
    sim_cat.set_index(['sid'], inplace=True)
    sgal_cat    = pd.read_pickle(path + 'sat_gal_table.pkl')
    sgal_cat.set_index(['sgal_id'], inplace=True)
    uniq_cat    = pd.read_pickle(path + 'tgal_tmp_thick_table.pkl')
    
    sgal_attributes_cat = pd.read_pickle(path + 'gal_R_attribute_cat.pkl')
    sgal_attributes_cat.set_index(['sgal_id'], inplace=True)
    
    important_col = ['center[0](code)', 'center[1](code)', 'center[2](code)', 'aexpn','galaxynumber','r_vir[kpc]', 'M_rvir[Msun]', 'Mdarkmatter(0.1rvir)', 'Mdarkmatter(0.2rvir)']
    cen_gal_cat_s = cen_gal_cat[cen_gal_cat['gen'] == 'VELA_v2'][cen_gal_cat['galaxynumber'] == galaxynumber][cen_gal_cat['type'] == 'Thick']

    tomer_list = cen_gal_cat_s[important_col]