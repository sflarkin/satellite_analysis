import numpy as np
import matplotlib.pyplot as plt


#need to use Mpeak for subhalos
def BehrooziWechslerConroy2013(a, z, halomass):
    epsilon1 = -1.777
    epsilon1e = [0.133, -0.146]
    epsilon2 = -0.006
    epsilon2e = [0.113, -0.361]
    epsilon3 = -0.000
    epsilon3e = [0.003, -0.104]
    epsilon4 = -0.119
    epsilon4e = [0.061, -0.012]
    
    M11 = 11.514
    M11e = [0.053, -0.009]
    M12 = -1.793
    M12e = [0.315, -0.330]
    M13 = -0.251
    M13e = [0.012, -0.125]
    
    alpha1 = -1.412 
    alpha1e = [0.020, -0.105]
    alpha2 = 0.731
    alpha2e = [0.344, -0.296]
    
    delta1 = 3.508
    delta1e = [0.087, -0.369]
    delta2 = 2.608
    delta2e = [2.446, -1.261]
    delta3 = -0.043
    delta3e = [0.958, -0.071]
    
    gamma1 = 0.316
    gamma1e = [0.076, -0.012]
    gamma2 = 1.319
    gamma2e = [0.584, -0.505]
    gamma3 = 0.279
    gamma3e = [0.256, -0.081]
    
    
    #these are the equations that generate the  intrinsic paramaters of https://arxiv.org/pdf/1207.6105.pdf
    #found on page 10 as well
    nu = np.exp(-4*a**2)
    
    def epsilon(a, z, nu, epsilon1, epsilon2, epsilon3, epsilon4, epsilon1e, epsilon2e, epsilon3e, epsilon4e):
        return((epsilon1 + epsilon1e) + ((epsilon2 + epsilon2e)*(a - 1) + (epsilon3 + epsilon3e)*z)*nu + (epsilon4 + epsilon4e)*(a - 1))
    
    def M1(a, z, nu, M11, M12, M13, M11e, M12e, M13e):
        return((M11 + M11e) + ((M12 + M12e)*(a - 1) + (M13 + M13e)*z)*nu)
        
    def alpha(a, z, nu, alpha1, alpha2, alpha1e, alpha2e):
        return((alpha1 + alpha1e) + ((alpha2 + alpha2e)*(a - 1))*nu)
    
    def delta(a, z, nu, delta1, delta2, delta3, delta1e, delta2e, delta3e):
        return((delta1 + delta1e) + ((delta2 + delta2e)*(a - 1) + (delta3 + delta3e)*z)*nu)
    
    def gamma(a, z, nu, gamma1, gamma2, gamma3, gamma1e, gamma2e, gamma3e):
        return((gamma1 + gamma1e) + ((gamma2 + gamma2e)*(a - 1) + (gamma3 + gamma3e)*z)*nu)
    
    #these are the acutal derived values of the above equations, to feed into f and log10stellarmass
    #_value is the value with no error
    #_error are the [max,min] errors after propogation, similar to gamma1e above, but instead of the +/-, the actual values
    #I used nested for loops since this was the fastest solution I thought of, def not the most efficent, prob save overall time due to execution
    epsilonvalue = epsilon(a, z, nu, epsilon1, epsilon2, epsilon3, epsilon4, 0, 0, 0, 0)
    epsilonerrorpropogation = []
    for epsilon1error in epsilon1e:
        for epsilon2error in epsilon2e:
            for epsilon3error in epsilon3e:
                for epsilon4error in epsilon4e:
                    epsilonerrorpropogation.append(epsilon(a, z, nu, epsilon1, epsilon2, epsilon3, epsilon4, epsilon1error, epsilon2error, epsilon3error, epsilon4error))
    epsilonerror = [max(epsilonerrorpropogation), min(epsilonerrorpropogation)]
    
    M1value = M1(a, z, nu, M11, M12, M13, 0, 0, 0)
    M1errorpropogation = []
    for M11error in M11e:
        for M12error in M12e:
            for M13error in M13e:
                M1errorpropogation.append(M1(a, z, nu, M11, M12, M13, M11error, M12error, M13error))
    M1error = [max(M1errorpropogation), min(M1errorpropogation)]
    
    alphavalue = alpha(a, z, nu, alpha1, alpha2, 0, 0)
    alphaerrorpropogation = []
    for alpha1error in alpha1e:
        for alpha2error in alpha2e:
            alphaerrorpropogation.append(alpha(a, z, nu, alpha1, alpha2, alpha1error, alpha2error))
    alphaerror = [max(alphaerrorpropogation), min(alphaerrorpropogation)]
    
    deltavalue = delta(a, z, nu, delta1, delta2, delta3, 0, 0, 0)
    deltaerrorpropogation = []
    for delta1error in delta1e:
        for delta2error in delta2e:
            for delta3error in delta3e:
                deltaerrorpropogation.append(delta(a, z, nu, delta1, delta2, delta3, delta1error, delta2error, delta3error))
    deltaerror = [max(deltaerrorpropogation), min(deltaerrorpropogation)]

    gammavalue = gamma(a, z, nu, gamma1, gamma2, gamma3, 0, 0, 0)
    gammaerrorpropogation = []
    for gamma1error in gamma1e:
        for gamma2error in gamma2e:
            for gamma3error in gamma3e:
                gammaerrorpropogation.append(gamma(a, z, nu, gamma1, gamma2, gamma3, gamma1error, gamma2error, gamma3error))
    gammaerror = [max(gammaerrorpropogation), min(gammaerrorpropogation)]
    
    #now that we have the max and min values from the error propogation,             
    def f(x, a, z, alpha, delta, gamma):
        return(-np.log10((10**(alpha*x) + 1)) + delta*(np.log10(1 + np.exp(x)))**gamma / (1 + np.exp(10**-x)))
    
    def log10stellarmass(halomass, a, z, epsilon, M1, alpha, delta, gamma):
        return((epsilon + M1) + f((np.log10(halomass) - M1), a, z, alpha, delta, gamma)  - f(0, a, z, alpha, delta, gamma))
    
    log10stellarmassvalue = log10stellarmass(halomass, a, z, epsilonvalue, M1value, alphavalue, deltavalue, gammavalue)
    log10stellarmasserrorpropogation = []
    for epsilonerrorvalue in epsilonerror:
        for M1errorvalue in M1error:
            for alphaerrorvalue in alphaerror:
                for deltaerrorvalue in deltaerror:
                    for gammaerrorvalue in gammaerror:
                        log10stellarmasserrorpropogation.append(log10stellarmass(halomass, a, z, epsilonerrorvalue, M1errorvalue, alphaerrorvalue, deltaerrorvalue, gammaerrorvalue))
        
    log10stellarmassplus = np.amax(log10stellarmasserrorpropogation, axis=0)
    log10stellarmassminus = np.amin(log10stellarmasserrorpropogation, axis=0)
    return([log10stellarmassvalue, log10stellarmassplus, log10stellarmassminus, (10**log10stellarmassvalue)/halomass, (10**log10stellarmassplus)/halomass, (10**log10stellarmassminus)/halomass])


def BehrooziWechslerConroy2013a(a, z, halomass):
    nu = np.exp(-4*a**2)
    
    log10epsilon = (-1.777 + 0) + ((-0.006 + 0)*(a - 1) + (-0.0000 + 0)*z)*nu + (-0.119 + 0)*(a - 1)
    log10epsilonplus = (-1.777 + 0.133) + ((-0.006 + 0.113)*(a - 1) + (-0.0000 + 0.003)*z)*nu + (-0.119 + 0.061)*(a - 1)
    log10epsilonminus = (-1.777 - 0.146) + ((-0.006 - 0.361)*(a - 1) + (-0.0000 - 0.104)*z)*nu + (-0.119 - 0.012)*(a - 1)
    
    log10M1 = (11.514 + 0) + ((-1.793 + 0)*(a - 1) + (-0.251 + 0)*z)*nu
    log10M1plus = (11.514 + 0.053) + ((-1.793 + 0.315)*(a - 1) + (-0.251 + 0.012)*z)*nu
    log10M1minus = (11.514 - 0.009) + ((-1.793 - 0.330)*(a - 1) + (-0.251 - 0.125)*z)*nu
    
    alpha = (-1.412 + 0) + ((-0.731 + 0)*(a - 1))*nu
    alphaplus = (-1.412 + 0.020) + ((-0.731 + 0.344)*(a - 1))*nu
    alphaminus = (-1.412 - 0.105) + ((-0.731 - 0.296)*(a - 1))*nu
    
    delta = (3.508 + 0) + ((2.608 + 0)*(a - 1) + (-0.043 + 0)*z)*nu
    deltaplus = (3.508 + 0.087) + ((2.608 + 2.446)*(a - 1) + (-0.043 + 0.958)*z)*nu
    deltaminus = (3.508 - 0.369) + ((2.608 - 1.261)*(a - 1) + (-0.043 - 0.071)*z)*nu
    
    gamma = (0.316 + 0) + ((1.319 + 0)*(a - 1) + (0.279 + 0)*z)*nu
    gammaplus = (0.316 + 0.076) + ((1.319 + 0.584)*(a - 1) + (0.279 + 0.256)*z)*nu
    gammaminus = (0.316 - 0.012) + ((1.319 - 0.505)*(a - 1) + (0.279 - 0.081)*z)*nu
    
    def f(x):
        return(-np.log10(10**(alpha*x) + 1) + delta*(np.log10(1 + np.exp(x)))**gamma / (1 + np.exp(10**-x)))
        
    def fplus(x):
        return(-np.log10(10**(alphaplus*x) + 1) + deltaplus*(np.log10(1 + np.exp(x)))**gammaplus / (1 + np.exp(10**-x)))
        
    def fminus(x):
        return(-np.log10(10**(alphaminus*x) + 1) + deltaminus*(np.log10(1 + np.exp(x)))**gammaminus / (1 + np.exp(10**-x)))
        
    log10stellarmass = (log10epsilon + log10M1) + f(np.log10(halomass) - log10M1) - f(0)
    log10stellarmassplus = (log10epsilonplus + log10M1plus) + f(np.log10(halomass) - log10M1plus) - fplus(0)
    log10stellarmassminus = (log10epsilonminus + log10M1minus) + f(np.log10(halomass) - log10M1minus) - fminus(0)
                                                                                                      
    #print(log10epsilon, log10epsilonplus, log10epsilonminus)
    #print(log10M1, log10M1plus, log10M1minus)
    #print(alpha, alphaplus, alphaminus)
    #print(delta, deltaplus, deltaminus)
    #print(gamma, gammaplus, gammaminus)
    return([log10stellarmass, log10stellarmassplus, log10stellarmassminus, (10**log10stellarmass)/halomass, (10**log10stellarmassplus)/halomass, (10**log10stellarmassminus)/halomass])
        
def PueblaPrimack2017(a, z, halomass):
    
    epsilon1 = -1.758
    epsilon1e = [0.040, -0.040]
    epsilon2 = 0.110
    epsilon2e = [0.166, -0.166]
    epsilon3 = -0.061
    epsilon3e = [0.029, -0.029]
    epsilon4 = -0.023
    epsilon4e = [0.009, -0.009]
    
    M11 = 11.548
    M11e = [0.049, -0.049]
    M12 = -1.297
    M12e = [0.225, -0.225]
    M13 = -0.026
    M13e = [0.043, -0.043]
    
    alpha1 = 1.975 
    alpha1e = [0.074, -0.074]
    alpha2 = 0.714
    alpha2e = [0.165, -0.165]
    alpha3 = 0.042
    alpha3e = [0.017, -0.017]
    
    delta1 = 3.390
    delta1e = [0.281, -0.281]
    delta2 = -0.472
    delta2e = [0.899, -0.899]
    delta3 = -0.931
    delta3e = [0.147, -0.147]
    
    gamma1 = 0.498
    gamma1e = [0.044, -0.044]
    gamma2 = -0.157
    gamma2e = [0.122, -0.122]
    
    nu = np.exp(-4*a**2)
    
    def p(x1,y1,z1):
        return y1*z1 - (x1*z1)/(1+z1)
    
    def epsilon(z, nu, epsilon1, epsilon2, epsilon3, epsilon4, epsilon1e, epsilon2e, epsilon3e, epsilon4e):
        return((epsilon1 + epsilon1e) + p((epsilon2 + epsilon2e), (epsilon3 + epsilon3e), z)*nu +
               p((epsilon4 + epsilon4e), 0, z))
    
    def M1(z, nu, M11, M12, M13, M11e, M12e, M13e):
        return((M11 + M11e) + p((M12 + M12e), (M13 + M13e), z)*nu)
    
    def alpha(z, nu, alpha1, alpha2, alpha3, alpha1e, alpha2e, alpha3e):
        return((alpha1 + alpha1e) + p((alpha2 + alpha2e), (alpha3 + alpha3e), z)*nu)
    
    def delta(z, nu, delta1, delta2, delta3, delta1e, delta2e, delta3e):
        return((delta1 + delta1e) + p((delta2 + delta2e), (delta3 + delta3e), z)*nu)
    
    def gamma(z, nu, gamma1, gamma2, gamma1e, gamma2e):
        return((gamma1 + gamma1e) + p((gamma2 + gamma2e), 0, z)*nu)
    
    epsilonvalue = epsilon(z, nu, epsilon1, epsilon2, epsilon3, epsilon4, 0, 0, 0, 0)
    epsilonerrorpropogation = []
    for epsilon1error in epsilon1e:
        for epsilon2error in epsilon2e:
            for epsilon3error in epsilon3e:
                for epsilon4error in epsilon4e:
                    epsilonerrorpropogation.append(epsilon(z, nu, epsilon1, epsilon2, epsilon3, epsilon4, epsilon1error, epsilon2error, epsilon3error, epsilon4error))
    epsilonerror = [max(epsilonerrorpropogation), min(epsilonerrorpropogation)]
    
    M1value = M1(z, nu, M11, M12, M13, 0, 0, 0)
    M1errorpropogation = []
    for M11error in M11e:
        for M12error in M12e:
            for M13error in M13e:
                M1errorpropogation.append(M1(z, nu, M11, M12, M13, M11error, M12error, M13error))
    M1error = [max(M1errorpropogation), min(M1errorpropogation)]
    
    alphavalue = alpha(z, nu, alpha1, alpha2, alpha3, 0, 0, 0)
    alphaerrorpropogation = []
    for alpha1error in alpha1e:
        for alpha2error in alpha2e:
            for alpha3error in alpha3e:
                alphaerrorpropogation.append(alpha(z, nu, alpha1, alpha2, alpha3, alpha1error, alpha2error, alpha3error))
    alphaerror = [max(alphaerrorpropogation), min(alphaerrorpropogation)]
    
    deltavalue = delta(z, nu, delta1, delta2, delta3, 0, 0, 0)
    deltaerrorpropogation = []
    for delta1error in delta1e:
        for delta2error in delta2e:
            for delta3error in delta3e:
                deltaerrorpropogation.append(delta(z, nu, delta1, delta2, delta3, delta1error, delta2error, delta3error))
    deltaerror = [max(deltaerrorpropogation), min(deltaerrorpropogation)]
    
    gammavalue = gamma(z, nu, gamma1, gamma2, 0, 0)
    gammaerrorpropogation = []
    for gamma1error in gamma1e:
        for gamma2error in gamma2e:
            gammaerrorpropogation.append(gamma(z, nu, gamma1, gamma2, gamma1error, gamma2error))
    gammaerror = [max(gammaerrorpropogation), min(gammaerrorpropogation)]
    
    #now that we have the max and min values from the error propogation,  
    
    def f(x, alpha, delta, gamma):
        return(-np.log10(10**(-alpha*x) + 1) + delta*(np.log10(1 + np.exp(x)))**gamma / (1 + np.exp(10**-x)))
    
    def log10stellarmass(halomass, epsilon, M1, alpha, delta, gamma):
        return((epsilon + M1) + f((np.log10(halomass) - M1), alpha, delta, gamma) - f(0, alpha, delta, gamma))
    
    log10stellarmassvalue = log10stellarmass(halomass, epsilonvalue, M1value, alphavalue, deltavalue, gammavalue)
    log10stellarmasserrorpropogation = []
    for epsilonerrorvalue in epsilonerror:
        for M1errorvalue in M1error:
            for alphaerrorvalue in alphaerror:
                for deltaerrorvalue in deltaerror:
                    for gammaerrorvalue in gammaerror:
                        log10stellarmasserrorpropogation.append(log10stellarmass(halomass, epsilonerrorvalue, M1errorvalue, alphaerrorvalue, deltaerrorvalue, gammaerrorvalue))                    
    log10stellarmassplus = np.amax(log10stellarmasserrorpropogation, axis=0)
    log10stellarmassminus = np.amin(log10stellarmasserrorpropogation, axis=0)
    return([log10stellarmassvalue, log10stellarmassplus, log10stellarmassminus, (10**log10stellarmassvalue)/halomass, (10**log10stellarmassplus)/halomass, (10**log10stellarmassminus)/halomass])

    
def PueblaPrimack2017a(a, z, halomass):
    nu = np.exp(-4*a**2)
    
    def p(x1,y1,z1):
        return y1*z1 - (x1*z1)/(1+z1)
    
    log10epsilon = (-1.758 + 0) + p((0.110 + 0), (-0.061 + 0), z)*nu  + p((-0.023 + 0), 0, z)
    log10epsilonplus = (-1.758 + 0.040) + p((0.110 + 0.166), (-0.061 + 0.029), z)*nu  + p((-0.023 + 0.009), 0, z)
    log10epsilonminus = (-1.758 - 0.040) + p((0.110 - 0.166), (-0.061 - 0.029), z)*nu  + p((-0.023 - 0.009), 0, z)
    
    log10M1 = (11.548 + 0) + p((-1.297 + 0), (-0.026 + 0), z)*nu
    log10M1plus = (11.548 + 0.049) + p((-1.297 + 0.225), (-0.026 + 0.043), z)*nu
    log10M1minus = (11.548 - 0.049) + p((-1.297 - 0.225), (-0.026 - 0.043), z)*nu
    
    alpha = (1.975 + 0) + p((0.714 + 0), (0.042 + 0), z)*nu
    alphaplus = (1.975 + 0.074) + p((0.714 + 0.165), (0.042 + 0.017), z)*nu
    alphaminus = (1.975 - 0.074) + p((0.714 - 0.165), (0.042 - 0.017), z)*nu
    
    delta = (3.390 + 0) + p((-0.472 + 0), (-0.931 + 0), z)*nu
    deltaplus = (3.390 + 0.281) + p((-0.472 + 0.899), (-0.931 + 0.147), z)*nu
    deltaminus = (3.390 - 0.281) + p((-0.472 - 0.899), (-0.931 - 0.147), z)*nu
    
    gamma = (0.498 + 0) + p((-0.157 + 0), 0, z)*nu
    gammaplus = (0.498 + 0.044) + p((-0.157 + 0.122), 0, z)*nu
    gammaminus = (0.498 - 0.044) + p((-0.157 - 0.122), 0, z)*nu
    
    #difference between these two is -alpha*x
    def f(x):
        return(-np.log10(10**(-alpha*x) + 1) + delta*(np.log10(1 + np.exp(x)))**gamma / (1 + np.exp(10**-x)))
        
    def fplus(x):
        return(-np.log10(10**(-alphaplus*x) + 1) + deltaplus*(np.log10(1 + np.exp(x)))**gammaplus / (1 + np.exp(10**-x)))
        
    def fminus(x):
        return(-np.log10(10**(-alphaminus*x) + 1) + deltaminus*(np.log10(1 + np.exp(x)))**gammaminus / (1 + np.exp(10**-x)))
    
    log10stellarmass = (log10epsilon + log10M1) + f(np.log10(halomass) - log10M1) - f(0)
    log10stellarmassplus = (log10epsilonplus + log10M1plus) + f(np.log10(halomass) - log10M1plus) - fplus(0)
    log10stellarmassminus = (log10epsilonminus + log10M1minus) + f(np.log10(halomass) - log10M1minus) - fminus(0)
    
    return([log10stellarmass, log10stellarmassplus, log10stellarmassminus, (10**log10stellarmass)/halomass, (10**log10stellarmassplus)/halomass, (10**log10stellarmassminus)/halomass])
    
    
    
    