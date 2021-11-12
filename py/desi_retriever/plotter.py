import matplotlib.pyplot as plt


def plot(D):
    plt.clf()
    plt.subplot(311)
    plt.plot(D['b_wavelength'], D['b_flux'], color='blue')
    plt.subplot(312)
    plt.plot(D['r_wavelength'], D['r_flux'], color='green')
    plt.subplot(313)
    plt.plot(D['z_wavelength'], D['z_flux'], color='red')
    plt.xlabel('Wavelength')
