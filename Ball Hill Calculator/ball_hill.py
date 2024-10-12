import matplotlib.pyplot as plt
import numpy 
import math

class ramp:
    def __init__(self, hyp, theta, coeff):
        self.hyp = hyp
        self.height = hyp*math.degrees(math.sin(theta))
        self.base = hyp*math.degrees(math.cos(theta))
        self.theta = theta
        self.coeff = coeff
        if not 0 <= coeff <= 1:
            print(f'invalid friction coeff: {coeff} is greater than one or less than zero, defaulting to 0')
            coeff = 0

    def get_static_params(self):
        return self.hyp, self.height, self.base, self.theta, self.coeff
    
class ball:
    def __init__(self, mass, radius, g):
        self.mass = mass
        self.radius = radius 
        self.moment = 0.5*mass*radius*radius
        self.x_pos = 0
        self.vel = 0
        self.rot_vel= 0
        self.acc = 0
        self.rot_acc = 0

    def get_static_params(self):
        return self.mass, self.radius, self.moment
    def set_x_pos(self, new_x_pos):
        self.x_pos = new_x_pos
    def get_x_pos(self):
        return self.x_pos
    def set_vel(self, new_vel):
        self.vel = new_vel
    def get_vel(self):
        return self.vel
    def set_rot_vel(self, new_rot_vel):
        self.rot_vel = new_rot_vel
    def get_rot_vel(self):
        return self.rot_vel
    def set_acc(self, new_acc):
        self.acc = new_acc
    def get_acc(self):
        return self.acc
    def set_rot_acc(self, new_rot_acc):
        self.rot_acc = new_rot_acc
    def get_rot_acc(self):
        return self.rot_acc
    
def calculate_initial_forces(ball:ball, ramp:ramp, g):
    hyp, height, base, theta, coeff = ramp.get_static_params()
    mass, radius, moment = ball.get_static_params()
    #find forces
    f_g_x = mass*g*math.sin(math.radians(theta))
    f_g_y = mass*g*math.cos(math.radians(theta))
    f_s_max = coeff*f_g_y
    if f_g_x - f_s_max > 0:
        acc = (2/3) * g * math.sin(math.radians(theta))
        rot_acc = acc / radius
        ball.set_acc(acc)
        ball.set_rot_acc(rot_acc)
    else: 
        print("ball doesn't translate nor rotate since gravitational force is not enough to overcome friction force")

def calculate_energy(ball:ball, ramp:ramp):
    hyp, height, base, theta, coeff = ramp.get_static_params()
    mass, radius, moment = ball.get_static_params()
    x = ball.get_x_pos()
    vel = ball.get_vel()
    rot_vel = ball.get_rot_vel()
    if not base - x <= 0:
        gravitational_energy = mass*g*((radius) + (base - x)*math.tan(math.radians(theta))) 
    else:
        return 1
    kinetic_energy = 0.5*mass*vel*vel
    rotational_energy = 0.5*moment*rot_vel*rot_vel
    total_energy = kinetic_energy+rotational_energy+gravitational_energy
    return total_energy, kinetic_energy, rotational_energy, gravitational_energy
    
def update_params(ball: ball, ramp: ramp, dt):
    x, vel, acc = ball.get_x_pos(), ball.get_vel(), ball.get_acc()
    rot_vel, rot_acc = ball.get_rot_vel(), ball.get_rot_acc()
    vel = vel + acc*dt
    ball.set_vel(vel)
    x = x + vel*dt + 0.5*acc*dt*dt
    ball.set_x_pos(x)
    rot_vel = rot_vel + rot_acc*dt 
    ball.set_rot_vel(rot_vel)
    return x, vel, rot_vel

def plot(x_axis, y_axis, top_label, label_x, units_x, label_y, units_y, color="red"):
    plt.title(f'{top_label} Graph') 
    plt.xlabel(f'{label_x} {units_x}') 
    plt.ylabel(f'{label_y} {units_y}') 
    plt.plot(x_axis, y_axis, color) 
    plt.show()

hyp = 999 #length of ramp hyp in meters
theta = 30 #angle of ramp from the horizontal 
coeff = 0

mass = 5 #mass in kg 
radius = 0.15 #radius in meters 

g = 9.80665 #gravitational constant 
timesteps = 1000 #simulation timesteps 
dt = 0.01 #change in time in seconds 

ramp_instance = ramp(hyp, theta, coeff)
ball_instance = ball(mass, radius, g)
calculate_initial_forces(ball_instance, ramp_instance,g)

#visualization
total_energy_array = numpy.array([])
ke_array = numpy.array([])
rot_array = numpy.array([])
gpe_array = numpy.array([])
total_time_array = numpy.array([])
total_time=0
for i in range(timesteps):
    if not calculate_energy(ball_instance, ramp_instance) == 1:
        total_energy, kinetic_energy, rotational_energy, gravitational_energy = calculate_energy(ball_instance, ramp_instance)

        x, vel, rot_vel = update_params(ball_instance, ramp_instance,dt)

        print(f'total energy: {total_energy} k_e {kinetic_energy} r_e {rotational_energy} g_p_e {gravitational_energy}')
        print(f'x: {x} vel: {vel} rot_vel {rot_vel}')
        total_energy_array = numpy.append(total_energy_array, total_energy)
        ke_array = numpy.append(ke_array, kinetic_energy)
        rot_array = numpy.append(rot_array, rotational_energy)
        gpe_array = numpy.append(gpe_array, gravitational_energy)
        total_time = total_time+dt
        total_time_array = numpy.append(total_time_array, total_time)
    else: 
        print(f'end of ramp reached')
        break

plot(total_time_array, total_energy_array, 'Total Simulation Energy', 'Time', '(s)', 'Energy', '(J)')
plot(total_time_array, ke_array, 'Kinetic Simulation Energy', 'Time', '(s)', 'Energy', '(J)')
plot(total_time_array, rot_array,'Rotational Simulation Energy', 'Time', '(s)', 'Energy', '(J)')
plot(total_time_array, gpe_array,'Gravitational Potential Simulation Energy', 'Time', '(s)', 'Energy', '(J)')