import numpy
from . import quantities
from amuse.units.si import *
from amuse.units.derivedsi import *
from amuse.units import constants

# misc every day
minute = named('minute', 'min', 60.0 * s)
hour   = named('hour',   'hr',  60.0 * minute)
day    = named('day',    'day', 24.0 * hour)
yr     = named('year',   'yr',  365.242199 * day)
julianyr = named('julian yr','julianyr',365.25* day)
ms = named('meter per seconds', 'ms', m / s)
kms = named('kilometer per seconds', 'kms', km / s)

# units based on measured quantities
e = named('electron charge', 'e', constants.elementary_charge.as_unit())
eV=named('electron volt','eV', e*V)
MeV=named('mega electron volt','MeV', 1e6*eV)
GeV=named('giga electron volt','GeV', 1e9*eV)
E_h = named('hartree energy', 'E_h', constants.Hartree_energy.as_unit())
amu = named('atomic mass unit', 'amu', constants.u.as_unit())
Ry = named('rydberg unit', 'Ry', (constants.Rydberg_constant * constants.h * constants.c).as_quantity_in(eV).as_unit())

# astronomical units
angstrom = named('angstrom', 'angstrom', 1e-10*m)
AU =  named('astronomical unit', 'AU', 149597870691.0  * m)
au =  named('astronomical unit', 'au', 149597870691.0  * m)
AUd = named('AU per day','AUd', 149597870691.0  * m / day)
parsec=named('parsec','parsec', AU / numpy.tan(numpy.pi/(180*60*60)))
kpc=named('kilo parsec','kpc',10**3 * parsec)
Mpc=named('mega parsec','Mpc',10**6 * parsec)
Gpc=named('giga parsec','Gpc',10**9 * parsec)
lightyear = named('light year', 'ly', 9460730472580.8 * km)
#lightyear = named('light year', 'ly', c*julianyr)
LSun = named('solar luminosity', 'LSun', 3.839e26 * W) 
MSun = named('solar mass', 'MSun', 1.98892e30 * kg)
MJupiter = named('jupiter mass', 'MJupiter', 1.8987e27 * kg)
MEarth = named('earth mass', 'MEarth', 5.9722e24 * kg)
RSun = named('solar radius', 'RSun', 6.955e8 * m)
RJupiter = named('jupiter radius', 'RJupiter', 71492. * km)
REarth = named('earth radius', 'REarth',  6371.0088 * km) # IUGG mean radius
kyr = named('kilo year', 'kyr', 1000 * yr)
myr = named('million year', 'Myr', 1000000 * yr)
Myr = myr
gyr = named('giga (billion) year', 'Gyr', 1000000000 * yr)
Gyr = gyr
pc = parsec

# cgs units
g = named('gram','g', 1e-3 * kg)
cm = named('centimeter','cm',0.01*m)
erg = named('erg','erg', 1e-7 * J)
barye = named('barye', 'Ba', 0.1*Pa)

# imperial distance units
inch = named('inch', 'in', 0.0254 * m)
foot = named('foot', 'ft', 0.3048 * m)
mile = named('mile', 'mi', 1609.344 * m)

percent = named('percent', '%', 0.01 * none)
metallicity = core.none_unit('metallicity', 'metallicity')

string = core.string_unit('string', 'string')

stellar_type = core.enumeration_unit(
    'stellar_type',
    'stellar_type',
    None,
    [
        "deeply or fully convective low mass MS star",  # 0
        "Main Sequence star",  # 1
        "Hertzsprung Gap",  # 2
        "First Giant Branch",  # 3
        "Core Helium Burning",  # 4
        "First Asymptotic Giant Branch",  # 5
        "Second Asymptotic Giant Branch",  # 6
        "Main Sequence Naked Helium star",  # 7
        "Hertzsprung Gap Naked Helium star",  # 8
        "Giant Branch Naked Helium star",  # 9
        "Helium White Dwarf",  # 10
        "Carbon/Oxygen White Dwarf",  # 11
        "Oxygen/Neon White Dwarf",  # 12
        "Neutron Star",  # 13
        "Black Hole",  # 14
        "Massless Supernova",  # 15
        "Unknown stellar type",  # 16
        "Pre-main-sequence Star",  # 17
        "Planet",  # 18
        "Brown Dwarf",  # 19
    ]
)

#special unit for keys of particles
object_key = core.key_unit('object_key','key')

#angles
#rad=named('radian','rad',m/m) (defined in derivedsi.py)
pi=numpy.pi | rad
rev=named('revolutions','rev',(2*numpy.pi) * rad)
deg=named('degree','deg',(numpy.pi/180) *  rad)
arcmin=named('arcminutes', 'arcmin', (1./60) * deg)
arcsec=named('arcseconds', 'arcsec', (1./3600) * deg)
