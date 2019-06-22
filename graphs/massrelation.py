import numpy as np
import matplotlib.pyplot as plt
from satellite_analysis.massrelations import stellarmass_relations

def massrelation(stellarmass, darkmattermass, massratio, halomass_dir, a):
    x = np.linspace(10**8, 10**14, 3000)
    Behroozirelation = stellarmass_relations.BehrooziWechslerConroy2013(a/1000, (1/(a/1000) - 1), x)
    Pueblarelation = stellarmass_relations.PueblaPrimack2017(a/1000, (1/(a/1000) - 1), x)

    
    plt.figure(1, figsize=(20,10))
    plt.subplot(121)
    plt.xlabel('halomass')
    plt.ylabel('massratio')
    plt.xscale('log')
    #plt.yscale('log')
    plt.xlim(10**(7.5), 10**(14))
    plt.ylim(0, 0.10)
    #plt.errorbar(x, Behroozirelation[3], yerr=np.absolute(Behroozirelation[5]-Behroozirelation[3]))
    plt.plot(x, Behroozirelation[3], 'k', color='#CC4F1B')
    plt.fill_between(x, Behroozirelation[5], Behroozirelation[4], alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.plot(x, Pueblarelation[3], 'k', color='#1B2ACC')
    plt.fill_between(x, Pueblarelation[5], Pueblarelation[4], alpha=0.5, edgecolor='#1B2ACC', facecolor='#089FFF')

    plt.scatter(darkmattermass, massratio, s = 1)
    
    
    plt.subplot(122)
    plt.xlabel('halomass')
    plt.ylabel('stellarmass')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(10**(7.8), 10**(14))
    #plt.ylim(10**(2), 10**(10))
    plt.plot(x, 10**Behroozirelation[0], 'k', color='#CC4F1B')
    plt.fill_between(x, 10**Behroozirelation[2], 10**Behroozirelation[1], alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    plt.plot(x, 10**Pueblarelation[0], 'k', color='#1B2ACC')
    plt.fill_between(x, 10**Pueblarelation[2], 10**Pueblarelation[1], alpha=0.5, edgecolor='#1B2ACC', facecolor='#089FFF')
    plt.scatter(darkmattermass, stellarmass, s = 1)
    plt.savefig('%s/massrelation%s.png' % (halomass_dir, str(a)))
    plt.close()