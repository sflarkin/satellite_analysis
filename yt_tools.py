def filters(par_type):
    global ad, masses
    ad = ds.all_data()
    masses = yt.np.unique(ad[(par_type, 'particle_mass')])
    for index in range(len(masses)):
        global filter_name, index1
        index1 = index
        filter_name = par_type + str(index)
        def mass_filter(pfilter, data):
            filter = data[(pfilter.filtered_type, 'particle_mass')] == masses[index1]
            return filter
        yt.add_particle_filter(filter_name, function=mass_filter, filtered_type=par_type, requires=['particle_mass'])
        ds.add_particle_filter(filter_name) 

