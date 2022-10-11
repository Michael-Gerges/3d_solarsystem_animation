## extremely inaccurate analytical calculation of the solar system  
# after one hour the difference in predicted and actual earth position is 185.29898071289062 61.61135292053223 -8.695259622531012
import math, itertools
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons

acceleration_step = 1  # second
number_of_steps = 3600 # iterations where each iteration add dt of one second
epo =  {'start': '2022-09-27 00:00 ', 'stop': '2022-09-27 11:00', 'step': "1m"}
index_of_result = int(number_of_steps*acceleration_step / 60) # 60 is seconds in 1m

au = 1.496e11 # austronomical unit
gravitational_constant = 6.6743 * 10**-11  #* 2.98692e-34 (in au)
#display_log_base = 1.3#  1.3
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
            acceleration = force / (body.mass/acceleration_step)          
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

