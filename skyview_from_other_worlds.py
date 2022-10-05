

# track the position of the sun and the moon 
# (or any celestial object across the sky)
# during a period of time 
# from any location on earth or other solar system planets)
# to remember the how the code works see the solar system 3d animation files for some comments


import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import astropy.units as u
from astropy.time import Time
from astroquery.jplhorizons import Horizons


planets_id = [10,301]  # just the sun and the moon # sun is not a planet but I just don't want to change the vairable name
colors = ["lightyellow", "lavender","gold"]
mylocation_dicto = {}
mylocation_dicto["lat"], mylocation_dicto["lon"] ,mylocation_dicto["elevation"], mylocation_dicto["body"]= 41.2579028 , -74.73456666666667 , 0.34702294557097116, 399
# Change mylocation_dicto["body"] for different views of the sky from other planets around the solar system
time = Time("2022-10-03 6:55").jd  # 4 hours diff from NewYork


def elaz_to_xyz(postions):
    el, az = postions["EL"] , postions["AZ"]
    el   = (90 - el)* np.pi/180 # to radian
    az     = az* np.pi/180
    x = np.sin( el ) * np.cos( az )
    y = np.sin( el ) * np.sin( az )
    z = np.cos( el )
    return np.array([x,y,z])




fig = plt.figure()
ax = p3.Axes3D(fig)

# create data and lines: 
data = []
for P_id in planets_id:
    postions = Horizons(id=P_id, location=mylocation_dicto ,epochs = {'start': '2022-10-03 10:56 ', 'stop': '2022-10-03 22:35 ', 'step': "1m"}).ephemerides()  # from 6 am to 6 pm utc + 4) (NewYork time)
    data.append(elaz_to_xyz(postions ))
# constructors 
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]
timestamp = ax.text(.03, .94, 15,'Date: ', color='purple', transform=ax.transAxes, fontsize='x-large')


planet_index = 0  # see the solar system 3d code to remember why I created that variable
def update_lines(num, dataLines, lines):
    global time
    global planet_index
    for line, data in zip(lines, dataLines):
        planet_index += 1
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
        line.set_color(colors[planet_index-1])
        line.set_linewidth(2)
        timestamp.set_text('Date: ' + Time(time, format='jd', out_subfmt='str').iso)
        time += 1/(24*60)    ## min
    planet_index = 0
    return lines + [timestamp]

line_ani = animation.FuncAnimation(fig, update_lines,len(data[0][0]), fargs=(data, lines),
                                   interval=0.05, blit=False, repeat=False)


# create the earth and the sky dome: 

u, v = np.mgrid[0:360:20j, 0:90:10j]  # u is the direction (360 deg all around) and v is the elevation (90 deg is above)
u   = u* np.pi / 180
v     = v* np.pi / 180 
x_sphere = np.cos(u)*np.sin(v)
y_shere = np.sin(u)*np.sin(v)
z_sphere = np.cos(v)
ax.plot_wireframe(x_sphere, y_shere, z_sphere, cmap="twilight",  linewidth=0.5,zorder=0)
ax.plot_surface(x_sphere, y_shere, np.zeros((20,10)), color="black")

ax.set_frame_on(False)

# the limit is always -1 to 1 because the radius is 1 (of the unit sphere)
ax.set_xlim3d([-1, 1.0])
ax.set_xlabel('South-North')
ax.set_ylim3d([0-1, 1.0])
ax.set_ylabel('West-East')
ax.set_zlim3d([-1, 1.0])
ax.set_zlabel('Altitude')
ax.set_title('3D Sunset' , color="black",fontweight ="bold")
ax.set_facecolor("lightpink")
ax.set_navigate_mode("PAN")
ax.set_navigate_mode("ZOOM")

plt.style.use('dark_background')
plt.show()
