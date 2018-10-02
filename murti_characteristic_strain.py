#Gravitational Observatory Advisory Team
import numpy as np
import matplotlib.pyplot as plt


def plot_characteristic_strain(mass1,mass2,distance):
    '''
    input: two masses in solar masses and luminosity distance in megaparsec

    output: plot of characteristic strain vs frequency
    return: chirp mass, array of dimensionless strain, array of characteristic strain
    '''

    #constants
    G = 6.67408e-11 # Gravitational constant
    c = 2.99e8      # speed of light

    #create frequency array
    freq = np.logspace(-5, 1, 1000)

    #convert distance from Mpc to meters
    distance = distance * 3.086e+22

    #convert solar masses to mass in kg
    mass1 *= 1.99e30
    mass2 *= 1.99e30

    #find chirp mass
    Mc = ((mass1 * mass2)**(3./5.)) / ((mass1 + mass2)**(1./5.))

    #find f dot
    f_dot = ((96. * c**3. * freq ) / (5. * G * Mc)) * ( (G/ c**3.) * np.pi * freq * Mc)**(8./3.)

    #find h0(dimensionless strain)
    dim_strain = ((4.* G * Mc )/ (c**2. * distance )) * ( G * np.pi * freq * Mc / c**3.)**(2./3.)

    #calculate characteristic strain
    char_strain = np.sqrt((freq**2 / f_dot)) * dim_strain

    #plot frequency vs characteristic strain
    plt.loglog(freq, char_strain, c='red')
    plt.title('characteristic strain vs frequency at redshift 3')
    plt.ylabel('characteristic strain')
    plt.yticks(np.arange(0, 1e-17, step=0.2))
    plt.xlabel('frequency (Hz)')
    plt.show()

    #return helpful values/arrays
    return Mc, dim_strain, char_strain

#from cosmological calculator
lum_distance = 25924.3 #mpc

m1 = ((1e7)/2) # solar masses
m2 = ((1e7)/2) # solar masses

#run function to produce plot
plot_characteristic_strain(m1,m2,lum_distance)
