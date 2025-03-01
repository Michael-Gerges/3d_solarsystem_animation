import numpy as np
from astroquery.jplhorizons import Horizons

# Function to fetch asteroid data for 2024 YR4
def get_trajectory(asteroid_id, start_date, end_date):
    # Use id_type 'smallbody' and set location to Earth's center ("@399")
    obj = Horizons(id=asteroid_id, id_type='smallbody', location='@399', 
                   epochs={'start': start_date, 'stop': end_date, 'step': '1h'})
    eph = obj.vectors()
    return eph['datetime_jd'], np.array(eph['x']), np.array(eph['y']), np.array(eph['z'])

# Fetch trajectory for 2024 YR4 over its close approach period
jd, x, y, z = get_trajectory("2024 YR4", "2032-12-01", "2032-12-31")

#import os
#
## Install astroquery
##os.system("pip install astroquery")
#
#import numpy as np
#import matplotlib.pyplot as plt
#from astroquery.jplhorizons import Horizons 
#
#astroid_dist = Horizons(id=99942, location="@399", epochs = {'start': '2029-04-13 00:00 ', 'stop': '2029-04-14 4:01 ', 'step': "1h"}).vectors()["range"] ## JPL units is in AU
##min(jupiter_dist)
#id_type='smallbody'
#
## Function to fetch asteroid data
#def get_trajectory(asteroid_id, start_date, end_date):
#    obj = Horizons(id=asteroid_id, id_type='smallbody' , location='@399', epochs={'start': start_date, 'stop': end_date, 'step': '1d'})
#    eph = obj.vectors()
#    return eph['datetime_jd'], eph['x'], eph['y'], eph['z']
#
## Define asteroid details
#asteroids = {
#    "99942 Apophis": (99942, "2029-04-01", "2029-04-30"),
#    "2024 YR4": ("2024 YR4", "2032-12-01", "2032-12-31")
#}
## Plotting
#fig = plt.figure(figsize=(8, 8))
#ax = fig.add_subplot(111, projection='3d')
#
#for name, (asteroid_id, start, end) in asteroids.items():
#    try:
#        jd, x, y, z = get_trajectory(asteroid_id, start, end)
#        ax.plot(x, y, z, label=name)
#    except Exception as e:
#        print(f"Error fetching data for {name}: {e}")
#
## Earth at (0,0,0)
#ax.scatter(0, 0, 0, color='blue', s=100, label="Earth")
#
#ax.set_xlabel("X (AU)")
#ax.set_ylabel("Y (AU)")
#ax.set_zlabel("Z (AU)")
#ax.legend()
#plt.title("Trajectories of 99942 Apophis & 2024 YR4")
#
## Show the plot
#plt.show()


#import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib.animation import FuncAnimation
#from astroquery.jplhorizons import Horizons
#
## Function to fetch asteroid data for 2024 YR4
#def get_trajectory(asteroid_id, start_date, end_date):
#    # Use id_type 'smallbody' and set location to Earth's center ("@399")
#    obj = Horizons(id=asteroid_id, id_type='smallbody', location='@399', 
#                   epochs={'start': start_date, 'stop': end_date, 'step': '1h'})
#    eph = obj.vectors()
#    return eph['datetime_jd'], np.array(eph['x']), np.array(eph['y']), np.array(eph['z'])
#
## Fetch trajectory for 2024 YR4 over its close approach period
#jd, x, y, z = get_trajectory("2024 YR4", "2032-12-01", "2032-12-31")
#
## Create a 3D plot for the animation
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111, projection='3d')
#
## Plot elements: full trajectory (as a red line) and current position (as a blue dot)
#line, = ax.plot([], [], [], 'r-', lw=2, label="2024 YR4 Trajectory")
#point, = ax.plot([], [], [], 'bo', markersize=8, label="Asteroid Position")
#
## Mark Earth at the origin
#ax.scatter(0, 0, 0, color='blue', s=100, label="Earth")
#
#ax.set_xlabel("X (AU)")
#ax.set_ylabel("Y (AU)")
#ax.set_zlabel("Z (AU)")
#ax.set_title("Animated Trajectory of 2024 YR4")
#ax.legend()
#
## Set fixed axis limits (adjust as needed for your data)
#ax.set_xlim(np.min(x)*1.1, np.max(x)*1.1)
#ax.set_ylim(np.min(y)*1.1, np.max(y)*1.1)
#ax.set_zlim(np.min(z)*1.1, np.max(z)*1.1)
#
#def init():
#    line.set_data([], [])
#    line.set_3d_properties([])
#    point.set_data([], [])
#    point.set_3d_properties([])
#    return line, point
#
#def update(frame):
#    # Update the trajectory line and current position marker
#    line.set_data(x[:frame], y[:frame])
#    line.set_3d_properties(z[:frame])
#    point.set_data(x[frame-1:frame], y[frame-1:frame])
#    point.set_3d_properties(z[frame-1:frame])
#    
#    # Optionally rotate the view for a dynamic perspective
#    ax.view_init(elev=30, azim=frame*0.5)
#    return line, point
#
## Create the animation: interval in milliseconds between frames
#ani = FuncAnimation(fig, update, frames=len(x), init_func=init, blit=True, interval=100)
#
#plt.show()
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from astroquery.jplhorizons import Horizons
from astropy.time import Time

# Function to fetch trajectory data (for objects like 2024 YR4)
def get_trajectory(asteroid_id, start_date, end_date, id_type='smallbody'):
    obj = Horizons(id=asteroid_id, id_type=id_type, location='@399',
                   epochs={'start': start_date, 'stop': end_date, 'step': '1h'})
    eph = obj.vectors()
    return np.array(eph['datetime_jd']), np.array(eph['x']), np.array(eph['y']), np.array(eph['z'])

# Define the time range for the animation (for example, December 2032)
start_date = "2032-12-01"
end_date   = "2032-12-31"

# Fetch trajectory for 2024 YR4 (the asteroid)
jd_ast, x_ast, y_ast, z_ast = get_trajectory("2024 YR4", start_date, end_date, id_type='smallbody')

# Fetch trajectory for the Sun (id '10') and Moon (id '301')
jd_sun, x_sun, y_sun, z_sun = get_trajectory("10", start_date, end_date, id_type=None)
jd_moon, x_moon, y_moon, z_moon = get_trajectory("301", start_date, end_date, id_type=None)

# Create a 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot full trajectory lines (asteroid in red, Sun in gold, Moon in cyan)
line_ast, = ax.plot(x_ast, y_ast, z_ast, 'r-', lw=1, label="2024 YR4 Trajectory")
line_sun, = ax.plot(x_sun, y_sun, z_sun, color='gold', lw=1, label="Sun Trajectory")
line_moon, = ax.plot(x_moon, y_moon, z_moon, color='c', lw=1, label="Moon Trajectory")

# Create markers for current positions
point_ast, = ax.plot([], [], [], marker='o', color='orange', markersize=8, linestyle='None', label="2024 YR4")
point_sun, = ax.plot([], [], [], marker='o', color='gold', markersize=10, linestyle='None', label="Sun")
point_moon, = ax.plot([], [], [], marker='o', color='c', markersize=8, linestyle='None', label="Moon")

# Plot Earth at the origin for reference
ax.scatter(0, 0, 0, color='blue', s=100, label="Earth")

# Add a timestamp text label in the upper left corner
timestamp_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes, color='black', fontsize=12)

# Set axis labels and title
ax.set_xlabel("X (AU)")
ax.set_ylabel("Y (AU)")
ax.set_zlabel("Z (AU)")
ax.set_title("Animated Trajectory of 2024 YR4 with Sun & Moon Reference")
ax.legend(loc="upper right")

# Set fixed axis limits (adjust these multipliers as needed)
ax.set_xlim(np.min(x_ast)*1.1, np.max(x_ast)*1.1)
ax.set_ylim(np.min(y_ast)*1.1, np.max(y_ast)*1.1)
ax.set_zlim(np.min(z_ast)*1.1, np.max(z_ast)*1.1)

def init():
    point_ast.set_data([], [])
    point_ast.set_3d_properties([])
    point_sun.set_data([], [])
    point_sun.set_3d_properties([])
    point_moon.set_data([], [])
    point_moon.set_3d_properties([])
    timestamp_text.set_text("")
    return point_ast, point_sun, point_moon, timestamp_text

def update(frame):
    # Update asteroid's current position
    point_ast.set_data(x_ast[frame:frame+1], y_ast[frame:frame+1])
    point_ast.set_3d_properties(z_ast[frame:frame+1])
    # Update Sun's current position
    point_sun.set_data(x_sun[frame:frame+1], y_sun[frame:frame+1])
    point_sun.set_3d_properties(z_sun[frame:frame+1])
    # Update Moon's current position
    point_moon.set_data(x_moon[frame:frame+1], y_moon[frame:frame+1])
    point_moon.set_3d_properties(z_moon[frame:frame+1])
    
    # Update timestamp label (convert current JD to ISO string)
    current_time = Time(jd_ast[frame], format='jd').iso
    timestamp_text.set_text(f"Time: {current_time}")
    
    # Rotate view slowly for dynamic perspective
    ax.view_init(elev=30, azim=frame*0.5)
    return point_ast, point_sun, point_moon, timestamp_text

# Create animation (interval is in milliseconds)
ani = FuncAnimation(fig, update, frames=len(x_ast), init_func=init, blit=True, interval=50)

plt.show()


import numpy as np
import plotly.graph_objects as go
from astroquery.jplhorizons import Horizons

def get_trajectory(obj_id, start_date, end_date, step='1h', id_type='smallbody', location='@sun'):
    """
    Fetch trajectory data from JPL Horizons.

    Parameters:
    -----------
    obj_id : str
        Object ID (e.g. '2024 YR4', '399' for Earth, '10' for Sun, etc.)
    start_date : str
        Start date in 'YYYY-MM-DD' format.
    end_date : str
        End date in 'YYYY-MM-DD' format.
    step : str
        Step size for ephemeris (e.g., '1h', '1d', etc.).
    id_type : str
        'smallbody', 'majorbody', or None. Usually 'smallbody' for asteroids, 'majorbody' for planets.
    location : str
        Reference location, e.g. '@sun' for heliocentric, '@0' for solar system barycenter, '@399' for geocentric.

    Returns:
    --------
    jd : np.array
        Julian dates for each ephemeris entry
    x, y, z : np.array
        Position coordinates (in AU, if using default Horizons settings)
    """
    obj = Horizons(
        id=obj_id,
        id_type=id_type,
        location=location,
        epochs={'start': start_date, 'stop': end_date, 'step': step}
    )
    vectors = obj.vectors()
    jd = np.array(vectors['datetime_jd'], dtype=float)
    x = np.array(vectors['x'], dtype=float)
    y = np.array(vectors['y'], dtype=float)
    z = np.array(vectors['z'], dtype=float)
    return jd, x, y, z

# ------------------------------------------------------------------------------
# 1) Define the time range for December 2032
# ------------------------------------------------------------------------------
start_date = "2032-12-01"
end_date   = "2032-12-31"

# ------------------------------------------------------------------------------
# 2) Get ephemeris data from Horizons
# ------------------------------------------------------------------------------
# 2.1) Asteroid 2024 YR4
jd_ast, x_ast, y_ast, z_ast = get_trajectory("2024 YR4", start_date, end_date, 
                                             step='1h', id_type='smallbody', 
                                             location='@sun')
# 2.2) Earth (planet center)
jd_earth, x_earth, y_earth, z_earth = get_trajectory("399", start_date, end_date, 
                                                     step='1h', id_type='majorbody', 
                                                     location='@sun')
# 2.3) Sun is effectively at origin in this frame
#     But let's also demonstrate how you'd fetch it explicitly (id=10, majorbody)
#     and confirm it is near (0,0,0).
jd_sun, x_sun, y_sun, z_sun = get_trajectory("10", start_date, end_date, 
                                            step='1h', id_type='majorbody',
                                            location='@sun')

# ------------------------------------------------------------------------------
# 3) Convert data to something manageable or keep in AU
#    (They should already be in AU from Horizons by default.)
# ------------------------------------------------------------------------------
# x_earth, y_earth, z_earth, etc. are already in AU.

# ------------------------------------------------------------------------------
# 4) Prepare the Plotly animation
# ------------------------------------------------------------------------------
num_frames = len(jd_earth)  # Earth & asteroid should have same length if steps match
# (If not, you'd need to handle mismatched arrays carefully.)

# Initialize figure
fig = go.Figure()

# Sun data: plot the Sun at each time step or just a big marker at the origin
# Because we’re in a sun-centered frame, (x_sun, y_sun, z_sun) should be near zero.
sun_marker = go.Scatter3d(
    x=[0], y=[0], z=[0],
    mode='markers',
    marker=dict(size=18, color='#FFA500', opacity=1),
    name="Sun"
)
fig.add_trace(sun_marker)

# Earth marker (initial position)
earth_marker = go.Scatter3d(
    x=[x_earth[0]], y=[y_earth[0]], z=[z_earth[0]],
    mode='markers',
    marker=dict(size=6, color='blue'),
    name="Earth"
)

# Asteroid marker (initial position)
ast_marker = go.Scatter3d(
    x=[x_ast[0]], y=[y_ast[0]], z=[z_ast[0]],
    mode='markers',
    marker=dict(size=4, color='red'),
    name="2024 YR4"
)

# Earth orbit trail (initially empty)
earth_trail = go.Scatter3d(
    x=[], y=[], z=[],
    mode='lines',
    line=dict(color='deepskyblue', width=2),
    name="Earth Orbit"
)

# Asteroid orbit trail (initially empty)
ast_trail = go.Scatter3d(
    x=[], y=[], z=[],
    mode='lines',
    line=dict(color='gray', width=2),
    name="Asteroid Orbit"
)

fig.add_trace(earth_trail)
fig.add_trace(ast_trail)
fig.add_trace(earth_marker)
fig.add_trace(ast_marker)

# ------------------------------------------------------------------------------
# 5) Build animation frames
# ------------------------------------------------------------------------------
frames = []
for k in range(1, num_frames):
    frames.append(go.Frame(
        data=[
            # Earth marker update
            go.Scatter3d(
                x=[x_earth[k]], 
                y=[y_earth[k]], 
                z=[z_earth[k]],
                mode='markers',
                marker=dict(size=6, color='blue')
            ),
            # Asteroid marker update
            go.Scatter3d(
                x=[x_ast[k]], 
                y=[y_ast[k]], 
                z=[z_ast[k]],
                mode='markers',
                marker=dict(size=4, color='red')
            ),
            # Earth trail update
            go.Scatter3d(
                x=x_earth[:k], 
                y=y_earth[:k], 
                z=z_earth[:k],
                mode='lines',
                line=dict(color='deepskyblue', width=2)
            ),
            # Asteroid trail update
            go.Scatter3d(
                x=x_ast[:k], 
                y=y_ast[:k], 
                z=z_ast[:k],
                mode='lines',
                line=dict(color='gray', width=2)
            ),
            # Sun marker stays the same (origin)
            go.Scatter3d(
                x=[0], 
                y=[0], 
                z=[0],
                mode='markers',
                marker=dict(size=18, color='#FFA500', opacity=1)
            )
        ],
        name=f"frame{k}"
    ))

fig.frames = frames

# ------------------------------------------------------------------------------
# 6) Layout and animation controls
# ------------------------------------------------------------------------------
fig.update_layout(
    title="2024 YR4 & Earth Orbit (Sun-Centered) - December 2032",
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False),
        bgcolor='black',
        # Adjust camera eye so you get a decent perspective
        camera=dict(
            eye=dict(x=1.8, y=1.8, z=0.9),
            up=dict(x=0, y=0, z=1)
        )
    ),
    showlegend=True,
    legend=dict(itemsizing='constant')
)

fig.update_layout(
    updatemenus=[
        dict(
            type="buttons",
            showactive=False,
            buttons=[
                dict(
                    label="Play",
                    method="animate",
                    args=[None, 
                          dict(frame=dict(duration=40, redraw=True), 
                               fromcurrent=True)]
                ),
                dict(
                    label="Pause",
                    method="animate",
                    args=[[None], 
                          dict(frame=dict(duration=0, redraw=False), 
                               mode="immediate")]
                )
            ]
        )
    ]
)

# Slider to scrub through frames
fig.update_layout(
    sliders=[dict(
        steps=[
            dict(
                method="animate",
                args=[
                    [f"frame{k}"],
                    dict(mode="immediate", frame=dict(duration=0, redraw=True))
                ],
                label=str(k)
            ) for k in range(num_frames)
        ],
        active=0,
        transition=dict(duration=0)
    )]
)

# ------------------------------------------------------------------------------
# 7) Show the figure!
# ------------------------------------------------------------------------------
fig.show()
