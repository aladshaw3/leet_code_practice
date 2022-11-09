# Allowed imports
import random
import numpy as np

# Must have this installed
import matplotlib.pyplot as plt

# Helper Function: Computes `Error`
#   between target and state (of device)
#
#   Example:
#       current_error
#           = error_value(current_target, current_dev_state)
def error_value(target, state):
    return target - state

## Finish this class by adding to the
#    `` Functions...
class Controller(object):
    # Default constructor
    def __init__(self):
        self.current_action = 0
        self.old_action = 0
        self.older_action = 0

    ############ FINISH THIS FUNCTION #################
    def ReturnAction(self, older_time: float,
                        old_time: float,
                        current_time: float,

                        older_target: float,
                        old_target: float,
                        current_target: float,

                        older_dev_state: float,
                        old_dev_state: float,
                        current_dev_state: float) -> float:
        # Simplest
        #current_error = error_value(current_target, current_dev_state)
        #return current_error

        # PID
        e = error_value(current_target, current_dev_state)
        e_old = error_value(old_target, old_dev_state)
        e_older = error_value(older_target, older_dev_state)
        dt = current_time - old_time
        dt_old = old_time - older_time
        act_old = self.old_action

        Kp = 0.1
        TI = 5
        TD = 5

        KI = Kp / TI
        KD = Kp * TD
        # NOTE: Removed the old_action + ...
        #   This makes this return a differential (i.e., dAct = compute)
        if dt_old > 0:
            return Kp*(e-e_old) + KI*e*dt + KD*((e-e_old)/dt) + KD*((e_old-e_older)/dt_old)
        else:
            return Kp*(e-e_old) + KI*e*dt + 0.5*KD*((e-e_old)/dt)

    ## DO NOT MODIFY THIS SECTION
    #   This is automatically called in the simulation loop
    def setAction(self, action: float):
        self.older_action = self.old_action
        self.old_action = self.current_action
        self.current_action = action

## END Controller Class


## Target Signal Class
class Target(object):
    # Default Constructor
    def __init__(self, ival: float):
        self.current_target = ival
        self.old_target = ival
        self.older_target = ival

    # Set and update targets
    def setTarget(self, target: float):
        self.older_target = self.old_target
        self.old_target = self.current_target
        self.current_target = target


## Time information
class Time(object):
    # Default Constructor
    def __init__(self):
        self.current_time = 0
        self.old_time = 0
        self.older_time = 0

    # Set and update targets
    def setTime(self, time: float):
        self.older_time = self.old_time
        self.old_time = self.current_time
        self.current_time = time


## Device Object
class Device(object):
    # Default Constructor
    def __init__(self):
        self.current_state = 0
        self.old_state = 0
        self.older_state = 0

    # Device state based on information from controller
    def setDeviceState(self, cont_obj: Controller):
        temp_state = cont_obj.current_action + self.current_state
        noise = random.uniform(0.85, 1.15)
        delay = random.randint(0,1)

        delay=0
        noise=1

        if (delay==0):
            temp_state = temp_state*noise
        else:
            temp_state = self.current_state

        self.older_state = self.old_state
        self.old_state = self.current_state
        self.current_state = temp_state



if __name__=="__main__":
    # Test
    time_obj = Time()
    target_obj = Target(5)
    device_obj = Device()
    controller_obj = Controller()

    act_set = []
    dev_set = []
    x_set = []
    target_set = []

    act_set.append(controller_obj.current_action)
    dev_set.append(device_obj.current_state)
    x_set.append(time_obj.current_time)
    target_set.append(target_obj.current_target)

    # simulation loop
    for i in range(0,80):
        new_time = (i+1)/2
        time_obj.setTime(new_time)

        if new_time > 20:
            target_obj.setTarget(2)

        action = controller_obj.ReturnAction(time_obj.older_time,
                              time_obj.old_time,
                              time_obj.current_time,

                              target_obj.older_target,
                              target_obj.old_target,
                              target_obj.current_target,

                              device_obj.older_state,
                              device_obj.old_state,
                              device_obj.current_state)
        controller_obj.setAction(action)
        print(action)

        device_obj.setDeviceState(controller_obj)

        dev_set.append(device_obj.current_state)
        act_set.append(controller_obj.current_action)
        x_set.append(time_obj.current_time)
        target_set.append(target_obj.current_target)

    print("DONE")
    xvals = list(x_set)
    fig,ax = plt.subplots(figsize=(10,5))
    yvals = list(act_set)
    ax.plot(xvals,yvals,color='red',label='action')
    yvals = list(target_set)
    ax.plot(xvals,yvals,color='blue',label='target')
    yvals = list(dev_set)
    ax.plot(xvals,yvals,color='black',label='device')
    plt.tight_layout()
    plt.savefig("test_pid_simulator2.png")
    fig.show()
    print("\nDisplaying plot. Press enter to continue...(this closes the images)")
    input()
    plt.close()
