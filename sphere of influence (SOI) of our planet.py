import spiceypy, os
# 1 AU in km

kerneldir = r"c:/Users/micha/OneDrive/Desktop/spaceproject/kernels/"
os.chdir(r"c:/Users/micha/OneDrive/Desktop/spaceproject/kernels/")
for i in os.listdir():
    filename = (os.path.join(kerneldir,i ))
    spiceypy.furnsh(filename)


ONE_AU = spiceypy.convrt(x=1, inunit='AU', outunit='km')  # a is the semimajor axis of the smaller object's (usually a planet's) orbit around the larger body (usually the Sun)

_, gm_sun_pre = spiceypy.bodvcd(bodyid=10, item='GM', maxn=1)

GM_SUN = gm_sun_pre[0]


# Set the G*M parameter of our planet
_, gm_earth_pre = spiceypy.bodvcd(bodyid=399, item='GM', maxn=1)
GM_EARTH = gm_earth_pre[0]

# Compute the SOI radius of the Earth
SOI_EARTH_R = ONE_AU * (GM_EARTH/GM_SUN) ** (2.0/5.0)

# Set one Lunar Distance (LD) in km (value from spaceweather.com)
ONE_LD = 384401.0 

print(f'SOI of the Earth in LD: {SOI_EARTH_R/ONE_LD}')

'''Earth + Moon	SOI radius 	0.929 (10^6 km) '''