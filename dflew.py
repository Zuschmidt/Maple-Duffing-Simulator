import numpy as np
import matplotlib.pyplot as plt

# Duffing equation parameters
def getTS(damping, stiffness, nonlinear, driving_force, freq):
    d = damping
    a = stiffness
    b = nonlinear
    df = driving_force
    w = freq

    # Initial conditions
    x0 = 0.0
    v0 = 0.0

    # Time points
    t_start = 0.0
    t_end = 100.0
    dt = (2*np.pi)/w/100
    num_steps = int((t_end - t_start) / dt)
    t = np.linspace(t_start, t_end, num_steps)

    # Duffing equation
    def duffing(t, x, v):
        dxdt = v
        dvdt = df * np.cos(w * t) - d * v - a * x - b * x**3
        return dxdt, dvdt

    # Runge-Kutta 4th order method
    def rk4_step(f, t, x, v, dt):
        k1_x, k1_v = f(t, x, v)
        k2_x, k2_v = f(t + 0.5 * dt, x + 0.5 * dt * k1_x, v + 0.5 * dt * k1_v)
        k3_x, k3_v = f(t + 0.5 * dt, x + 0.5 * dt * k2_x, v + 0.5 * dt * k2_v)
        k4_x, k4_v = f(t + dt, x + dt * k3_x, v + dt * k3_v)
        
        x_next = x + (dt / 6.0) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        v_next = v + (dt / 6.0) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)
        
        return x_next, v_next

    # Arrays to store the results
    x = np.zeros(num_steps)
    v = np.zeros(num_steps)

    # Initial conditions
    x[0] = x0
    v[0] = v0

    # Time-stepping using RK4
    for i in range(1, num_steps):
        x[i], v[i] = rk4_step(duffing, t[i-1], x[i-1], v[i-1], dt)
    return x

# Print time series
# for i in range(num_steps):
#     print(f"time: {t[i]}, x: {x[i]}, v: {v[i]}")
# print(x, v)

# Plot the results
# plt.plot(t, x, label='x(t)')
# plt.plot(t, v, label='v(t)')
# plt.xlabel('Time')
# plt.ylabel('Values')
# plt.legend()
# plt.title('Duffing Oscillator Time Series (Runge-Kutta)')
# plt.show()










import numpy.matlib as matlib
import numpy as np
import sys # to implement args later
import copy
import warnings

def LyE_W(x, Fs, tau, dim, evolve):
    """
    inputs  - x, time series
            - Fs, sampling frequency
            - tau, time lag
            - dim, embedding dimension
            - evolve, parameter of the same name from Wolf's 1985 paper. This
              code expects a number of frames as an input.
    outputs - out, matrix detailing variables at each iteration
            - LyE, largest lyapunov exponent
    """
    
    x = np.array([x])
    SCALEMX = (np.max(x)-np.min(x))/10
    ANGLMX = 30*np.pi/180
    ZMULT = 1

    DT = 1/Fs

    ITS = 0
    distSUM = 0

    if np.size(x, axis=0) == 1:
        m = dim
        N = np.size(x, axis=1)
        M = N-(m-1)*tau
        Y = np.zeros((M,m))
        
        for i in range(0,m):
            Y[:,i]=x[:,(0+i*tau):(M + i*tau)]

        NPT=np.size(x, axis=1)-(dim-1)*tau-evolve # Size of useable data
        Y=Y[0:NPT+evolve,:] 

    else:
        Y=np.array(x)
        NPT=np.size(Y, axis=0)-evolve

    out=np.zeros((int(np.floor(NPT/evolve)+1),9),dtype="object")
    thbest=0
    OUTMX=SCALEMX

    # Find first pair
            
    # Distance from current point to all other points
    current_point = 0

    Yinit = matlib.repmat(Y[current_point],NPT,1)
    Ydiff = (Yinit - Y[0:NPT,:])**2
    Ydisti = np.sqrt(np.sum(Ydiff,1))
        
    # Exclude points too close on path and close in distance
    range_exclude = np.arange(current_point-10,current_point+10+1)
    range_exclude = range_exclude[(range_exclude>=0) & (range_exclude < NPT)]
    Ydisti[Ydisti<=0] = np.nan
    Ydisti[range_exclude] = np.nan
        
    # find minimum distance point for first pair
    current_point_pair = np.argsort(Ydisti)[0]

    for i in range(0, NPT, evolve):
        current_point = i
        # calculate starting and evolved distance
        if current_point_pair + evolve < len(Y) and current_point + evolve < len(Y):
            start_dist = np.linalg.norm(Y[current_point,:] - Y[current_point_pair,:])
            end_dist = np.linalg.norm(Y[current_point+evolve,:] - Y[current_point_pair+evolve,:])
        else:
            start_dist = np.linalg.norm(Y[current_point,:] - Y[current_point_pair,:])
            end_dist = np.linalg.norm(Y[current_point+evolve,:] - Y[current_point_pair+evolve-1,:])
        
        # calculate total distance so far
        distSUM = distSUM + np.log2(end_dist/start_dist)/(evolve*DT) # DT is sampling rate?!
        ITS = ITS+1  # count iterations
        LyE=distSUM/ITS # max Lyapunov exponent
        
        # Store found pairs
        out[int(np.floor(i/evolve))] = [ITS,current_point,current_point_pair,start_dist,end_dist,LyE,OUTMX,(thbest*180/np.pi),(ANGLMX*180/np.pi)]
        
        ZMULT=1
        
        if end_dist < SCALEMX:
            current_point_pair = current_point_pair+evolve
            if current_point_pair > NPT:
                current_point_pair = current_point_pair - evolve
                flag = 1
                (current_point_pair,ZMULT,ANGLMX,thbest,OUTMX) = get_next_point(flag,Y,current_point,current_point_pair,NPT,evolve,SCALEMX,ZMULT,ANGLMX)
            continue
        # find point pairing for next iteration
        flag=0
        (current_point_pair,ZMULT,ANGLMX,thbest,OUTMX) = get_next_point(flag,Y,current_point,current_point_pair,NPT,evolve,SCALEMX,ZMULT,ANGLMX)
        
    return (out, LyE)

def get_next_point(flag, Y, current_point, current_point_pair, NPT, evolve,SCALEMX, ZMULT, ANGLMX):

    # Distance from evolved point to all other points
    Yinit = np.matlib.repmat(Y[current_point+evolve,:],NPT,1)
    Ydiff = (Yinit - Y[0:NPT,:])**2
    Ydisti = np.sqrt(np.sum(Ydiff,axis=1))

    # Exclude points too close on path and close in distance than noise
    range_exclude = np.arange(current_point+evolve-10,current_point+evolve+10+1)
    range_exclude = range_exclude[(range_exclude >= 0) & (range_exclude < NPT)]
    Ydisti[range_exclude] = np.nan

    if current_point_pair + evolve < len(Y) and current_point + evolve < len(Y):
        end_dist = np.linalg.norm(Y[current_point+evolve,:] - Y[current_point_pair+evolve,:])
    else:
        end_dist = np.linalg.norm(Y[current_point+evolve,:] - Y[current_point_pair+evolve-1,:])

    # Vector from evolved point to all other points
    Vnew = np.matlib.repmat(Y[current_point+evolve,:],NPT,1) - Y[:NPT,:]

    # Vector from evolved point to evolved point pair 
    if current_point_pair + evolve < len(Y) and current_point + evolve < len(Y):
        PT1 = Y[current_point+evolve,:]
        PT2 = Y[current_point_pair+evolve,:]
    else:
        PT1 = Y[current_point+evolve,:]
        PT2 = Y[current_point_pair+evolve-1,:]
    Vcurr = PT1-PT2

    # Angle between evolved pair vector and all other vectors
    # Clipping values of cosTheta to be within the range [-1, 1]
    cosTheta = np.clip(np.abs(np.divide(np.sum(Vcurr.T*Vnew, axis=1),(Ydisti*end_dist))), -1, 1)
    theta = np.arccos(cosTheta)

    # Search for next point
    # -1 Meaning point not found.
    next_point=-1
    while next_point == -1:
        (next_point,ZMULT,ANGLMX,thbest,SCALEMX)=find_next_point(flag,theta,Ydisti,SCALEMX,ZMULT,ANGLMX)
    
    return next_point, ZMULT, ANGLMX, thbest, SCALEMX

def find_next_point(flag,theta,Ydisti, SCALEMX, ZMULT, ANGLMX):

    # Restrict search based on distance and angle
    PotenDisti= np.copy(Ydisti)
    PotenDisti[(Ydisti<=0) | (theta>=ANGLMX)] = np.nan

    next_point=-1
    if flag==0:
        next_point = np.argsort(PotenDisti)[0]
        # if closest angle point is within angle range -> point found and reset
        # search space
        if PotenDisti[next_point] <= SCALEMX:
            ANGLMX = 30*np.pi/180
            thbest=np.abs(theta[next_point])
            return (next_point, ZMULT, ANGLMX, thbest, SCALEMX)
        else:
            next_point=-1
            flag=1
    if flag == 1:
        PotenDisti=np.copy(Ydisti)
        PotenDisti[Ydisti<=0] = np.nan
        next_point = np.argsort(PotenDisti)[0]
        thbest=ANGLMX
    
    return (next_point, ZMULT, ANGLMX, thbest, SCALEMX)




LE = []
coe = []
for i in range(1, 300):
    out, LyE = LyE_W(getTS(i/10, 1, 1, .855, 1), 1000, 1, 2, 400)
    LE.append(LyE)
    coe.append(i/10)

out, LyE = LyE_W(getTS(2.1, 1, 1, .855, 1), 1000, 1, 2, 400)
print(LyE)


# print(LE)
# plt.plot(coe, LE)
# plt.xlabel('Damping')
# plt.ylabel('LE')
# plt.show()