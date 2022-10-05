##I don't have any idea what is going on with JWST

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from matplotlib import colors as mcolors

ISS_id =  -125544
JWST_id = -170 
moon_id = 301
earth_id = 399
l4_id = 3014 #EM-L4 (Lagrange point)
l5_id =  3015 #EM-L5 (Lagrange point)

idlist  = [3011, 3012,3013 , l4_id, l5_id, moon_id,  JWST_id] 
sizes =[9] * len(idlist)
colors = list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())

sim_start_date = "2021-12-26"     
sim_end_date = "2023-10-03" 
timestep, dt = "1d" , 1 
sim_interval = 5 # the less the faster
simulation_display_size =.03 # AU # the actual area of space that the animation displays 
rotating_around = "@399"  # EARTH

def get_data():
    global sim_start_date,  sim_end_date, timestep , idlist, rotating_around
    for index, object_id in enumerate(idlist):
        obj_vectors = Horizons(id=object_id, location=rotating_around, epochs = {'start': sim_start_date, 'stop': sim_end_date, 'step': timestep}).vectors()
        ranges = obj_vectors["range"].value
        globals()[f"cor_{index}_x"] = obj_vectors["x"].value
        globals()[f"cor_{index}_y"] = obj_vectors["y"].value
    variables_dict = {k: v for k, v in globals().items()  if  k.startswith("cor_")}
    return variables_dict , ranges  # ranges are not needed for the simulation but I got it just in case lol


variables_dict , ranges = get_data()


def getpositionofasatelliteonaday(object_index_in_idlist, timeindex):
    global variables_dict
    corx = variables_dict[f"cor_{object_index_in_idlist}_x"]
    cory = variables_dict[f"cor_{object_index_in_idlist}_y"]
    return corx[timeindex],  cory[timeindex]


class Object:                  
    def __init__(self, name, rad, color, D_0, range ):   # D_0 is initial position (x,y) 
        self.name = name                            # its actually the index in the id list
        self.range = range
        self.D_0    = np.array(D_0, dtype=np.float)
        self.xs = []
        self.ys = []
        self.plot = ax.scatter(D_0[0], D_0[1], cmap=color, s=rad, edgecolors=None, zorder=10)
        self.line, = ax.plot([], [], color=color, linewidth=1.4)



#def get_center_postion(time_index):
#    return 0,0



time_index = 0
class System:
    def __init__(self, center):
        self.center = center
        self.time = Time(sim_start_date).jd
        self.timestamp = ax.text(.03, .94, 'Date: ', color='w', transform=ax.transAxes, fontsize='x-large')
        self.satellites = []
    def add_satellites(self, satellite):
        self.satellites.append(satellite)
    def evolve(self):   
        global time_index
        mycenter = self.center
        self.time += dt
        plots , lines = [], [] 
        #a,b = get_center_postion(time_index)
        #self.center.plot.set_offsets([a,b])
        #plots.append(self.center.plot)
        for s in self.satellites:
            if time_index > len(variables_dict["cor_0_x"]) -1:
                break
            a,b = getpositionofasatelliteonaday(s.name, time_index )  
            s.xs.append(a)
            s.ys.append(b)
            s.plot.set_offsets( [a,b]  )
            s.line.set_xdata(s.xs)      
            s.line.set_ydata(s.ys)
            plots.append(s.plot)
            lines.append(s.line)
        
        
        self.timestamp.set_text('Date: ' + Time(self.time, format='jd', out_subfmt='str').iso )
        time_index += 1
        return plots + lines+ [self.timestamp]

plt.style.use('dark_background')
fig = plt.figure(figsize=[8, 8])
ax = plt.axes([0., 0., 1., 1.], xlim=(-simulation_display_size, simulation_display_size), ylim=(-simulation_display_size,simulation_display_size))
ax.set_aspect('equal')
ax.axis('off')

jwst_location_system = System(Object("bla bla bal", 225, 'blue', [0, 0], range=0)) 

for i, nasa_id in enumerate(idlist):
    x_0 , y_0 = getpositionofasatelliteonaday(i, 0) # for initialization
    jwst_location_system.add_satellites(Object(i,  sizes[i], colors[i], [x_0, y_0] , range=0))                            


def animate(_):
    return jwst_location_system.evolve()
ani = animation.FuncAnimation(fig, animate, repeat=False, frames=len(variables_dict["cor_0_x"]), blit=True, interval=sim_interval)     
plt.show()

