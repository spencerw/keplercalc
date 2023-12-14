import numpy as np
from mathhelpers import cross, dot, nr, PQW

def hello(name):
    return 'Hello ' + name + '!'

def orb_params(snap, isHelio=False, mCentral=1.0):
    """
    Takes a Pynbody snapshot of particles and calculates gives them fields
    corresponding to Kepler orbital elements. Assumes that the central star
    is particle #0 and that the positions and velocities are in barycentric
    coordinates.

    Parameters
    ----------
    snap: SimArray
        A snapshot of the particles
    isHelio: boolean
        Skip frame transformation if the snap is already in heliocentric coordinates
    mCentral: float
        Mass of central star in simulation units. Need to provide if no star particle in snap

    Returns
    -------
    SimArray
        A copy of 'snap' with additional fields 'a', 'ecc', 'inc', 'asc_node',
        'omega' and 'M' added.
    """

    x = np.array(snap.d['pos'])
    x_h = x[1:] - x[0]
    v = np.array(snap.d['vel'])
    v_h = v[1:] - v[0]
    m1 = np.array(snap['mass'][0])
    pl = snap[1:]
    m2 = np.array(pl['mass'])

    if isHelio:
        x_h = x
        v_h = v
        pl = snap
        m1 = mCentral

    pl['a'], pl['e'], pl['inc'], pl['asc_node'], pl['omega'], pl['M'] \
        = cart2kep(x_h[:,0], x_h[:,1], x_h[:,2], v_h[:,0], v_h[:,1], v_h[:,2], m1, m2)
    
    return pl

def cart2kep(X, Y, Z, vx, vy, vz, m1, m2):
    """
    Convert an array of cartesian positions and velocities
    (in heliocentric coordinates) to an array of kepler orbital
    elements. Units are such that G=1. This function is
    fully vectorized.

    Parameters
    ----------
    X: numpy array of floats
        X position of body
    Y: numpy array of floats
        Y position of body
    Z: numpy array of floats
        Z position of body
    vx: numpy array of floats
        x velocity of body
    vy: numpy array of floats
        y velocity of body
    vz: numpy array of floats
        z velocity of body
    m1: float
        mass of central body
    m2: numpy array of floats
        mass of orbiting body

    Returns
    -------
    numpy array of floats
        a, e, inc, asc_node, omega, M - Semimajor axis, eccentricity
        inclination, longitude of ascending node, longitude of 
        perihelion and mean anomaly of orbiting body
    """

    mu = m1 + m2
    magr = np.sqrt(X ** 2. + Y ** 2. + Z ** 2.)

    hx, hy, hz = cross(X, Y, Z, vx, vy, vz)
    magh = np.sqrt(hx ** 2. + hy ** 2. + hz ** 2.)
    tmpx, tmpy, tmpz = cross(vx, vy, vz, hx, hy, hz)
    evecx = tmpx / mu - X / magr
    evecy = tmpy / mu - Y / magr
    evecz = tmpz / mu - Z / magr

    e = np.sqrt(evecx ** 2. + evecy ** 2. + evecz ** 2.)

    a = dot(hx, hy, hz, hx, hy, hz) / (mu * (1. - e ** 2.))

    ivec = [1., 0., 0.]
    jvec = [0., 1., 0.]
    kvec = [0., 0., 1.]

    inc = np.arccos(dot(kvec[0], kvec[1], kvec[2], hx, hy, hz) / magh)

    nx, ny, nz = cross(kvec[0], kvec[1], kvec[2], hx, hy, hz)
    nmag = np.sqrt(nx ** 2. + ny ** 2. + nz ** 2.)
    asc_node = np.where(inc == 0., 0., np.arccos(dot(ivec[0], ivec[1], ivec[2], nx, ny, nz) / nmag))
    asc_node[dot(nx, ny, nz, jvec[0], jvec[1], jvec[2]) < 0.] = 2. * np.pi - asc_node[
        dot(nx, ny, nz, jvec[0], jvec[1], jvec[2]) < 0.]

    omega = np.where(inc == 0., np.arctan2(evecy / e, evecx / e),
                     np.arccos(dot(nx, ny, nz, evecx, evecy, evecz) / (nmag * e)))
    omega[dot(evecx, evecy, evecz, kvec[0], kvec[1], kvec[2]) < 0.] = 2. * np.pi - omega[
        dot(evecx, evecy, evecz, kvec[0], kvec[1], kvec[2]) < 0.]

    theta = np.arccos(dot(evecx, evecy, evecz, X, Y, Z) / (e * magr))
    theta = np.where(dot(X, Y, Z, vx, vy, vz) < 0., 2 * np.pi - theta, theta)

    E = np.arccos((e + np.cos(theta)) / (1 + e * np.cos(theta)))
    E = np.where(np.logical_and(theta > np.pi, theta < 2 * np.pi), 2 * np.pi - E, theta)
    M = E - e * np.sin(E)

    return a, e, inc, asc_node % (2 * np.pi), omega, M

def kep2cart(sma, ecc, inc, Omega, omega, M, mass, m_central):
    """
    Convert a single set of kepler orbital elements into cartesian
    positions and velocities. Units are such that G=1.

    Parameters
    ----------
    sma: float
        semi-major axis body
    ecc: float
        eccentricity of body
    inc: float
        inclination of body
    Omega: float
        longitude of ascending node of body
    omega: float
        longitude of perihelion of body
    M: float
        Mean anomaly of body
    mass: float
        mass of body
    m_central: float
        mass of central body

    Returns
    -------
    float
        X, Y, Z, vx, vy, vz - Cartesian positions and velocities
        of the bodies
    """

    if inc == 0.:
        Omega = 0.
    if ecc == 0.:
        omega = 0.
    E = nr(M, ecc)
    
    X = sma * (np.cos(E) - ecc)
    Y = sma * np.sqrt(1. - ecc**2.) * np.sin(E)
    mu = m_central + mass
    n = np.sqrt(mu / sma**3.)
    Edot = n / (1. - ecc * np.cos(E))
    Vx = - sma * np.sin(E) * Edot
    Vy = sma * np.sqrt(1. - ecc**2.) * Edot * np.cos(E)
    
    Px, Py, Pz, Qx, Qy, Qz = PQW(Omega, omega, inc)

    # Rotate Positions
    x = X * Px + Y * Qx
    y = X * Py + Y * Qy
    z = X * Pz + Y * Qz

    # Rotate Velocities
    vx = Vx * Px + Vy * Qx
    vy = Vx * Py + Vy * Qy
    vz = Vx * Pz + Vy * Qz
    
    return x, y, z, vx, vy, vz

def kep2poinc(a, e, i, omega, Omega, M , m1, m2):
    """
    Convert kepler orbital elements into Poincaire variables.
    Can accept either single values or lists of coordinates

    Parameters
    ----------
    a: float
        semi-major axis body
    ecc: float
        eccentricity of body
    inc: float
        inclination of body
    Omega: float
        longitude of ascending node of body
    omega: float
        longitude of perihelion of body
    M: float
        Mean anomaly of body
    mass: float
        mass of body
    m_central: float
        mass of central body

    Returns
    -------
    float
        lam, gam, z, Lam, Gam, Z - Poincaire coordinates
    """

    mu = m1 + m2
    mustar = m1*m2/(m1 + m2)

    lam = M + omega + Omega
    gam = -omega - Omega
    z = -Omega
    Lam = mustar*np.sqrt(mu*a)
    Gam = mustar*np.sqrt(mu*a)*(1 - np.sqrt(1 - e**2))
    Z = mustar*np.sqrt(mu*a*np.sqrt(1 - e**2))*(1 - np.cos(i))

    return lam, gam, z, Lam, Gam, Z

def kep2del(a, e, i, omega, Omega, M , m1, m2):
    """
    Convert kepler orbital elements into Delunay variables.
    Can accept either single values or lists of coordinates

    Parameters
    ----------
    a: float
        semi-major axis body
    ecc: float
        eccentricity of body
    inc: float
        inclination of body
    Omega: float
        longitude of ascending node of body
    omega: float
        longitude of perihelion of body
    M: float
        Mean anomaly of body
    mass: float
        mass of body
    m_central: float
        mass of central body

    Returns
    -------
    float
        l, g, h, L, G, H - Delunay coordinates
    """

    L = np.sqrt((m1 + m2)*a)
    G = L*np.sqrt(1 - e**2)
    H = G*np.cos(i)
    lv = M
    g = omega
    h = Omega

    return lv, g, h, L, G, H

def kep2mdel(a, e, i, omega, Omega, M, m1, m2):
    """
    Convert kepler orbital elements into modified Delunay variables.
    Can accept either single values or lists of coordinates

    Parameters
    ----------
    a: float
        semi-major axis body
    ecc: float
        eccentricity of body
    inc: float
        inclination of body
    Omega: float
        longitude of ascending node of body
    omega: float
        longitude of perihelion of body
    M: float
        Mean anomaly of body
    mass: float
        mass of body
    m_central: float
        mass of central body

    Returns
    -------
    float
        Lam, P, Q, lam, p, q - Delunay coordinates
    """
    
    L = np.sqrt((m1 + m2)*a)
    G = L*np.sqrt(1 - e**2)
    H = G*np.cos(i)
    lv = M
    g = omega
    h = Omega

    Lam = L
    P = L - G
    Q = G - H
    lam = lv + g + h
    p = -g - h
    q = -h
    return Lam, P, Q, lam, p, q