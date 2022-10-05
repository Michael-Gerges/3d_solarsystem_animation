import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = "2022-08-01"     
sim_end_date = "2023-08-06" 
sim_duration =  159                # orbit 53-day
day = 0
timestep, dt = "1d" , 1 # the dt is actually used to increment a list index so it must be integers


id_list = [-61]  # for juno and whatever other satellites 
colors = ['gray']
sizes = [6]


def get_data():
    obj1 = Horizons(id=id_list[0], location="@599", epochs = {'start': sim_start_date, 'stop': sim_end_date, 'step': timestep}).vectors()
    cor1y, cor1z  = obj1["y"].value, obj1["z"].value
    return cor1y, cor1z

cor1y, cor1z = get_data()  # get all the data at the same time so I don't use the API multiple times during the function calls

def getpositionofasatelliteonaday(day):
    global cor1y, cor1z
    return cor1y[day],  cor1z[day]


class Object:                  
    def __init__(self, name, rad, color, D_0 ):   # D_0 is initial position (x,y,z) 
        self.name = name
        self.D_0    = np.array(D_0, dtype=np.float)
        self.ys = []
        self.zs = []
        self.plot = ax.scatter(D_0[0], D_0[1], color=color, s=rad**2, edgecolors=None, zorder=10)
        self.line, = ax.plot([], [], color=color, linewidth=1.4)



class System:
    def __init__(self, thecenter):
        self.thecenter = thecenter
        self.satellites = []
        self.time = Time(sim_start_date).jd
        self.timestamp = ax.text(.03, .94, 'Date: ', color='w', transform=ax.transAxes, fontsize='x-large')
    def add_satellites(self, satellite):
        self.satellites.append(satellite)
    def evolve(self):           
        global day
        self.time += dt
        plots = []
        lines = []
        day += 1 
        for s in self.satellites:
            s.D_0[0] , s.D_0[1] = getpositionofasatelliteonaday(day)  
            s.ys.append(s.D_0[0])
            s.zs.append(s.D_0[1])
            s.plot.set_offsets( [s.D_0[0], s.D_0[1]] )
            s.line.set_xdata(s.ys)      # cast the y, and z dimensions on set_xdata and set_ydata because 'Line2D' object has no attribute 'set_zdata'
            s.line.set_ydata(s.zs)
            plots.append(s.plot)
            lines.append(s.line)
        self.timestamp.set_text('Date: ' + Time(self.time, format='jd', out_subfmt='str').iso)
        return plots + lines+ [self.timestamp]


plt.style.use('dark_background')
fig = plt.figure(figsize=[6, 6])
factor = 0.1 # AU # the actual area of space that the animation displays 
ax = plt.axes([0., 0., 1., 1.], xlim=(-factor, factor), ylim=(-factor,factor))
ax.set_aspect('equal')
ax.axis('off')

Jupiter_Juno_system = System(Object("Jupiter", 28, 'red', [0, 0, 0]))
for i, nasa_id in enumerate(id_list):
    y_0 , z_0 = getpositionofasatelliteonaday(0) # for initialization
    Jupiter_Juno_system.add_satellites(Object(nasa_id,  sizes[i], colors[i], [y_0, z_0] ))                            


def animate(_):
    return Jupiter_Juno_system.evolve()
ani = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=True, interval=20,)     
plt.show()





