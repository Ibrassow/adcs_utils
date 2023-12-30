import matplotlib.pyplot as plt 
import numpy as np


def dark(enabled):
    if enabled:
        plt.style.use("dark_background")
    if not enabled:
        plt.style.use('default')


def plot3D(X, center=False, equal_axis=False):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z') 

    if center:
        ax.scatter(0,0,0,c='g',label="Center", marker='o', s=300)

    ax.plot3D(X[:,0], X[:,1], X[:,2], c='b')

    if equal_axis:
        max_val = np.max(np.abs(X))
        # Set the same scale for all axes
        ax.set_xlim([-max_val, max_val])
        ax.set_ylim([-max_val, max_val])
        ax.set_zlim([-max_val, max_val])

    plt.legend()
    plt.show()



def plotAttitude(X, dt):

    t = np.arange(0, X.shape[0]*dt, step=dt)

    plt.figure(1)
    plt.xlabel('Time [s]')
    plt.ylabel('Quaternion [/]')
    plt.plot(t, X[:,0:4])

    plt.figure(2)
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [rad/s]')
    plt.plot(t, X[:,4:7])

    
    plt.show()


def plot3axis(x, dt=None, legend=False, ylabel=None):

    if dt != None:
        t = np.arange(0, x.shape[0]*dt, step=dt)
        plt.plot(t, x[:,0], label="x")
        plt.plot(t, x[:,1], label="y")
        plt.plot(t, x[:,2], label="z")
        plt.xlabel('Time [s]')
    else:
        plt.plot(x[:,0], label="x")
        plt.plot(x[:,1], label="y")
        plt.plot(x[:,2], label="z")

    if ylabel != None:
        plt.ylabel(ylabel)
        
    if legend:
        plt.legend()

    
    plt.show()

def splot(x):
    plt.plot(x)
    plt.show()



def compare(w1, w2):

    plt.plot(w1, c='b', label="1")
    plt.plot(w2, c='r', label="2")
    plt.legend()
    plt.show()




def plot_scatter(x, t=None):
    plt.plot(x, ".")
    plt.show()


def scatter3D(x):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = x[:, 0]
    ys = x[:, 1]
    zs = x[:, 2]

    ax.scatter(xs, ys, zs)
    plt.show()


