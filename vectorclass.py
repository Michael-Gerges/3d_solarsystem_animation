## extremely inaccurate analytical calculation of the solar system  
# after one hour the difference in predicted and actual earth position is 185.29898071289062 61.61135292053223 -8.695259622531012
import math, itertools
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons
import numpy as np 

acceleration_step = 1  # second
number_of_steps = 3600 # 
epo =  {'start': '2022-09-27 00:00 ', 'stop': '2022-09-27 11:00', 'step': "1m"}
index_of_result = int(number_of_steps*acceleration_step / 60) # 60 is seconds in 1m

au = 1.496e11
gravitational_constant = 6.6743 * 10**-11  #* 2.98692e-34
display_log_base = 1.3#  1.3
sunmass = 1.989 * 10**30 
namelist = [ "MERCURY" ,"VENUS" ,"EARTH", "MOON", "MARS","JUPITER","SATURN","URANUS", "NEPTUNE","PLUTO" ]
masslist = [0.330,4.87,5.97,0.073,0.642,1898,568,86.8,102,0.0130] # 10**24
idlist = [199,299,399,499,599,699,799,899,999]


class Vector:
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Vector({self.x}, {self.y}, {self.z})"

    def __str__(self):
        return f"{self.x}i + {self.y}j + {self.z}k"

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        elif item == 2:
            return self.z
        else:
            raise IndexError("There are only three elements in the vector")

    def __add__(self, other):
        return Vector(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __sub__(self, other):
        return Vector(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __mul__(self, other):
        if isinstance(other, Vector):  # Vector dot product
            return (
                self.x * other.x
                + self.y * other.y
                + self.z * other.z
            )
        elif isinstance(other, (int, float)):  # Scalar multiplication
            return Vector(
                self.x * other,
                self.y * other,
                self.z * other,
            )
        else:
            raise TypeError("operand must be Vector, int, or float")

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Vector(
                self.x / other,
                self.y / other,
                self.z / other,
            )
        else:
            raise TypeError("operand must be int or float")

    def get_magnitude(self):
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        magnitude = self.get_magnitude()
        return Vector(
            self.x / magnitude,
            self.y / magnitude,
            self.z / magnitude,
        )
    def clone(self):
        return Vector(self.x, self.y, self.z)
# Scalar Product
            
    def dot(self,other):
        """Returns dot product of two vectors"""
        dotproduct = 0.0
        dotproduct += self.x*other.x
        dotproduct += self.y*other.y
        dotproduct += self.z*other.z
        return dotproduct
    
# Vector Product
    def cross(self,other):
        """Calculates cross product (self x other)"""

        cross = Vector(0.0,0.0,0.0)
        cross.x = self.y * other.z - self.z * other.y
        cross.y = self.z * other.x - self.x * other.z
        cross.z = self.x * other.y - self.y * other.x
        return cross

    
# Rotate around the Z axis

    def rotateX(self,angle):
        """Rotates the vector around the x axis"""
        oldvec = self.clone()
        self.x - oldvec.x
        self.y = oldvec.y*np.cos(angle) - oldvec.z*np.sin(angle)
        self.z = oldvec.y*np.sin(angle) + oldvec.z*np.cos(angle)

    def rotateY(self,angle):
        """Rotates the vector around the y axis"""
        oldvec = self.clone()
        self.x = oldvec.x*np.cos(angle) + oldvec.z*np.sin(angle)
        self.y = oldvec.y
        self.z = -oldvec.x*np.sin(angle) + oldvec.z*np.cos(angle)

    def rotateZ(self, angle):
        """Rotates the vector around the z axis"""
 
        oldvec = self.clone()

        self.x = oldvec.x*np.cos(angle) - oldvec.y*np.sin(angle)
        self.y = oldvec.x*np.sin(angle) + oldvec.y*np.cos(angle)
        self.z = oldvec.z

    


class SolarSystem:
    def __init__(self, size, projection_2d=False):
        self.size = size
        self.projection_2d = projection_2d
        self.bodies = []

        self.fig, self.ax = plt.subplots(1,1, subplot_kw={"projection": "3d"}, figsize=(self.size / 50, self.size / 50))
        self.fig.tight_layout()
        if self.projection_2d:
            self.ax.view_init(10, 0)
        else:
            self.ax.view_init(0, 0)

    def add_body(self, body):
        self.bodies.append(body)
    def update_all(self):
        self.bodies.sort(key=lambda item: item.position[0])
        for body in self.bodies:
            body.move()
            body.draw()
    def draw_all(self):
        if self.projection_2d:
            self.ax.xaxis.set_ticklabels([])
            self.ax.yaxis.set_ticklabels([])
            self.ax.zaxis.set_ticklabels([])
        else:
            self.ax.axis(False)
        plt.pause(0.001)
        self.ax.clear()
    def calculate_all_body_interactions(self):
        bodies_copy = self.bodies.copy()
        for idx, first in enumerate(bodies_copy):
            for second in bodies_copy[idx + 1:]:
                first.accelerate_due_to_gravity(second)


class SolarSystemBody:
    display_log_base = 10.3
    min_display_size = 0
    def __init__(
        self,
        solar_system,
        mass,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        self.solar_system = solar_system
        self.mass = mass
        self.position = position
        self.velocity = Vector(*velocity)
        self.display_size = max(
            math.log(self.mass, self.display_log_base),
            self.min_display_size,
        )
        self.color = "black"
        self.solar_system.add_body(self)
    def move(self):
        self.position = (
            self.position[0] + self.velocity[0],
            self.position[1] + self.velocity[1],
            self.position[2] + self.velocity[2],
        )
    def draw(self):
        postion_tuble = [self.position[0] / au , self.position[1] / au, self.position[2] / au] 
        self.solar_system.ax.plot(
            postion_tuble[0]  , postion_tuble[1] , postion_tuble[2] ,
            marker="o",
            color=self.color
        )
        if self.solar_system.projection_2d:
            self.solar_system.ax.plot(
                self.position[0],
                self.position[1],
                -self.solar_system.size / 2,
                marker="o",
                #markersize=self.display_size / 2,
                color=(.5, .5, .5),
            )
    def accelerate_due_to_gravity(self, other):
        global gravitational_constant 
        distance = Vector(*other.position) - Vector(*self.position)
        distance_mag = distance.get_magnitude()
        force_mag = gravitational_constant * (self.mass * other.mass / (distance_mag ** 2))
        force = distance.normalize() * force_mag 
        reverse = 1
        for body in self, other:
            acceleration = force / (body.mass/acceleration_step)          # accelaration in a day  # 86400
            body.velocity += acceleration * reverse
            reverse = -1


class Sun(SolarSystemBody):
    def __init__(
        self,
        solar_system,
        mass=10_000,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        super(Sun, self).__init__(solar_system, mass, position, velocity)
        self.color = "black"

class Planet(SolarSystemBody):
    colors = itertools.cycle([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
    def __init__(
        self,
        solar_system,
        mass=10,
        position=(0, 0, 0),
        velocity=(0, 0, 0),
    ):
        super(Planet, self).__init__(solar_system, mass, position, velocity)
        self.color = next(Planet.colors)



#### mock system (remove the gravitational_constant and au and set     display_log_base = 1.3 and)


#solar_system = SolarSystem(400 , projection_2d=True)
#suns = (
#    Sun(solar_system, position=(40, 40, 40), velocity=(6, 0, 6)),
#    Sun(solar_system, position=(-40, -40, 40), velocity=(-6, 0, -6)),
#)
#planets = (
#    Planet(
#        solar_system,
#        10,
#        position=(100, 100, 0),
#        velocity=(0, 5.5, 5.5),
#    ),
#    Planet(
#        solar_system,
#        20,
#        position=(0, 0, 0),
#        velocity=(-11, 11, 0),
#    ),
#)
#
#def go():
#    now = datetime.datetime.utcnow().minute
#    print(now)
#    while datetime.datetime.utcnow().minute < now +1 :
#        solar_system.calculate_all_body_interactions()
#        solar_system.update_all()
#        solar_system.draw_all()
#
#
#go()
#




#
#solar_system_2 = SolarSystem(50  )#, projection_2d=True)
#Sun(solar_system_2, mass=sunmass , position=(0, 0, 0), velocity=(0, 0, 0))



#MERCURY_vector = (Horizons(idlist[0], epochs=epo)).vectors()
#x = (MERCURY_vector["x"][0])
#y = (MERCURY_vector["y"][0])
#z = (MERCURY_vector["z"][0])
#vx = (MERCURY_vector["vx"][0])
#vy = (MERCURY_vector["vy"][0])
#vz = (MERCURY_vector["vz"][0])
#MERCURY = Planet(solar_system_2, masslist[0] * 10**24, position=(x * au, y* au, z* au), velocity=(vx* au/86400, vy *au/86400, vz* au/86400) )

#VENUS_vector = (Horizons(idlist[1], epochs=epo)).vectors()
#x = (VENUS_vector["x"][0])
#y = (VENUS_vector["y"][0])
#z = (VENUS_vector["z"][0])
#vx = (VENUS_vector["vx"][0])
#vy = (VENUS_vector["vy"][0])
#vz = (VENUS_vector["vz"][0])
#VENUS = Planet(solar_system_2, masslist[1] * 10**24, position=(x * au, y* au, z* au), velocity=(vx* au/86400, vy *au/86400, vz* au/86400) )


#EARTH_vector = (Horizons(idlist[2], epochs=epo)).vectors()
#x = (EARTH_vector["x"][0])
#y = (EARTH_vector["y"][0])
#z = (EARTH_vector["z"][0])
#vx = (EARTH_vector["vx"][0])
#vy = (EARTH_vector["vy"][0])
#vz = (EARTH_vector["vz"][0])
#end_vector = Vector(EARTH_vector["x"][index_of_result] *au   ,(EARTH_vector["y"][index_of_result]) * au,(EARTH_vector["z"][index_of_result]) * au)



#EARTH = Planet(solar_system_2, masslist[2] * 10**24, position=(x * au, y* au, z* au), velocity=(vx* au/86400, vy *au/86400, vz* au/86400) )

#MARS_vector = (Horizons(idlist[3], epochs=epo)).vectors()
#x = (MARS_vector["x"][0])
#y = (MARS_vector["y"][0])
#z = (MARS_vector["z"][0])
#vx = (MARS_vector["vx"][0])
#vy = (MARS_vector["vy"][0])
#vz = (MARS_vector["vz"][0])
#MARS = Planet(solar_system_2, masslist[3] * 10**24, position=(x * au, y* au, z* au), velocity=(vx* au/86400, vy *au/86400, vz* au/86400) )



#for i in range(number_of_steps):
#        solar_system_2.calculate_all_body_interactions()
#        solar_system_2.update_all()
#        solar_system_2.draw_all()
#



#print(MERCURY.position[0]- (MERCURY_vector["x"][index_of_result]) * au, MERCURY.position[1]  -(MERCURY_vector["y"][index_of_result]) *au, MERCURY.position[2] - (MERCURY_vector["z"][index_of_result]) * au)
#print(VENUS.position [0] - (VENUS_vector["x"][index_of_result]) *au   , VENUS.position[1]   - (VENUS_vector["y"][index_of_result]) * au, VENUS.position[2]    - (VENUS_vector["z"][index_of_result]) * au)
#print(EARTH.position[0]  - (EARTH_vector["x"][index_of_result]) *au   ,EARTH.position[1]   - (EARTH_vector["y"][index_of_result]) * au, EARTH.position[2]     - (EARTH_vector["z"][index_of_result]) * au)
#print(MARS.position[0]   - (MARS_vector["x"][index_of_result]) *au    ,MARS.position[1]     - (MARS_vector["y"][index_of_result]) * au, MARS.position[2]      - (MARS_vector["z"][index_of_result]) * au)

#print(MERCURY_vector["datetime_str"][index_of_result])



