## please dont make fun of me I am trying to learn how to figure this out, I have no formal education of any of that, all self-taught, thanks


import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import astropy.units as u
plt.style.use('dark_background') # from astropy.visualization import astropy_mpl_style  or plt.style.use(astropy_mpl_style)    #


starttime = '2022-10-03 10:56'
endtime =  '2029-12-09 22:35'
timestep, dt = "1d" , 1 # if 1h or 1m divide by 24 or 24 and 60 
time = Time(starttime) + 4*u.hour  # for the timestamp: Eastern Daylight Time 

planets_id = [199,299, 399,499,599, 699, 799] 
colors = ["brown", "pink","blue", "red","orange", "yellow", "green"]

data = []
for P_id in planets_id:
    # location @10 is the center of the sun
    postions = Horizons(id=P_id, location="@10" ,epochs = {'start':starttime , 'stop': endtime, 'step': timestep}).vectors()  # from 6 am to 6 pm utc + 4)
    x = postions["x"]
    y = postions["y"]
    z = postions["z"]
    data.append(np.array([x,y,z] ))  
    
'''
    ### that will get 3 dimensional matrix the inner most layer is the time dimension from JPL positions["x"] , positions["y"] , positions["z"] that are basically lists of the variable x for a specific planet 
    the second layer (middle layer represented by the middle dimension of the 3D matrix) ... the second layer is the np array of xyz (which I constructed using: np.array([x,y,z]) for a planet at a point of time (which physically correlates to the spatial dimensions )
    the outer most layer is what the for loop take care of. It physically represents the dimension of the variety planets. Those planets each of which has its own time and space (dimensions)

'''

### after collecting the data, its time to build the figure

fig = plt.figure()
ax = p3.Axes3D(fig)

## use three costructors. those updates the points (scattered ponints of the moving planets) 
## and the line (that trace the movement of planets)
## and the timestamp that passes with each new update ... 

timestamp = ax.text(.03, .94, 15,'Date: ', color='purple', transform=ax.transAxes, fontsize='x-large')
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data] #  for planet amoung planets
plots = [ax.scatter(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1]) for dat in data]

'''
## I made the variable "planet index" to do line.set_color(colors[planet_index-1]) 
# 
# because the function update_lines() operates on the time dimension, with each call representing a jump in time frame
# ..... which means that, it iterates on the inner most dimension of the big matrix 'dataline'  
the for loop then takes it from here and iterate ("for each instance of time”) ... through

that’s because the for loop iterate on the outer most layer of the matrix dataLines the planets 
... (if there are 5 planets the loop will go 5 times for each planet when the function is called once ... (the function is called with progress of time)
... that means that the loop handles the outer most dimension of the datalines matrix

the plot and line constructors (scatter and plot functions) iterate through the spatial dimension xyz
... to give the information of each point to the graph to instruct the planet about the movement
... (actually, instructing the computer where to put the dots lol)
... this is the middle layer of the big matrix 

... consequently, through the iteration of the loop (much like the first for P_id in planets_id:), it just gets the data of planet after planet. 
 ..... slicing the datatline matrix to get 2 dimensional smaller matrices each with represents 4D of the planet (space and time dimensions)
 .... that slices that are handled inside doesn't contain the planet name ... 
 ... so I created a counter to see which planet is now its turn on the loop to get photogrpahed for the constructors 
 ... and used that counter to till the color list .... line.set_color(colors[planet_index-1]) ... which planet it is to color it accordingly



'''
def update_lines(num, dataLines, lines, plots):     # called with each instance of time
    global time     # sorry Einstein but in that case the time is universal not relative, lol  (defined outside all iterators)
    ## I made that variable to do line.set_color(colors[planet_index-1]) 
    planet_index = 0  # reset the planet index to get ready for a fresh instece of time through which we go through the planets again and again 
    for line,data, plot in zip(lines,dataLines, plots): # give me picture by picture  and out of those pictures I will get the information that pertain to that instance of time for which the function was called 
        planet_index += 1   # each planet gets its own for loop iteration
        # the plot scatter has 2 different construction methods                   
        line.set_data(data[0:2, :num])   # that is the plot taking the information of the first 2 spatial dimension        
        line.set_3d_properties(data[2, :num])       # and the third (he cannot take all three together lol (3 is a big amount tbh lol))
        
        plot.set_offsets(np.array(dataLines)[:,0:2,num])        # and the scatter constructor taking the first 2 spatial dims 
        plot.set_3d_properties(np.array(dataLines)[:,2,num],  zdir='z')     # and the final one lol 
        
        # for some reason this next two lines of codes are different: 
        line.set_color(colors[planet_index-1])      # apparently ax.plot jumps checks the points as it draws them while the ax.scatter swallow up the data till it get a full instance of time then draws that instance of time 
        plot.set_color(colors[planet_index-1]) 
        
        line.set_linewidth(2) # line.set_zorder(5) #line.set_solid_capstyle("round")

    ## after looping throgh planets we increase the time variable to make the time passes for the new call of the function (updaterf)

    timestamp.set_text('Date: ' + Time(time, format='jd', out_subfmt='str').iso)
    time += dt   # increased with call of the function that update it (not inside the loop that goes over planets)
    planet_index = 0        
    return lines, plots+[timestamp]



## The function animation.FuncAnimation askes the following questions:
#  What? 
# who's goanna till me how to draw, 
# How many times you want me to draw that (same as how many time data point we have), 
# (give me the data and the constructors for plot and scatter) 
# How fast?
line_ani = animation.FuncAnimation(fig, update_lines,len(data[0][0]), fargs=(data, lines,plots),
                                   interval=.0001, blit=False, repeat=False)

distance_drawn = 2 # Astronomical units (that’s the unit that JPL uses by default on Horizons)
ax.set_xlim3d([-distance_drawn, distance_drawn])
ax.set_ylim3d([-distance_drawn, distance_drawn])
ax.set_zlim3d([-distance_drawn, distance_drawn])
### in the dark mode: 
ax.set_axis_off()

ax.scatter(0, 0, 0,color="orange", s=200, zorder=10) ## don't forget to put a sun at the origin 

plt.show()


