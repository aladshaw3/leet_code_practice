import numpy as np

# Must have this installed
import matplotlib.pyplot as plt

# Defines PID action var
def compute_action(Kp, TI, TD, e, e_old, e_older, action_old, dt, dt_old):
    KI = Kp / TI
    KD = Kp * TD
    if dt_old > 0:
        return action_old + Kp*(e-e_old) + KI*e*dt + KD*((e-e_old)/dt) + KD*((e_old-e_older)/dt_old)
    else:
        return action_old + Kp*(e-e_old) + KI*e*dt + 0.5*KD*((e-e_old)/dt)

# Defines error
def error_func(r, y):
    return r-y

# Assume action == measured value

# test
current_setpoint = 5
current_measured_value = 0
old_measured_value = 0
older_measured_value = 0

Kp = 0.1
TI = 5
TD = 5

time = 0
old_time = 0
old_action = 0

y_set = []
x_set = []
target_set = []

y_set.append(old_action)
x_set.append(old_time)
target_set.append(current_setpoint)

# test loop
for i in range(0,1000):
    new_time = (i+1)/2
    dt = new_time - time
    dt_old = time - old_time

    if new_time > 20:
        current_setpoint = 2

    new_error = error_func(current_setpoint, current_measured_value)
    if abs(new_error) < 1e-4:
        break

    error = error_func(current_setpoint, old_measured_value)
    old_error = error_func(current_setpoint, older_measured_value)

    new_action = compute_action(Kp, TI, TD, new_error, error,
                                old_error, old_action, dt, dt_old)

    older_measured_value = old_measured_value
    old_measured_value = current_measured_value
    current_measured_value = new_action
    old_action = new_action

    y_set.append(new_action)
    x_set.append(new_time)
    target_set.append(current_setpoint)

print("DONE")
xvals = list(x_set)
fig,ax = plt.subplots(figsize=(10,5))
yvals = list(y_set)
ax.plot(xvals,yvals)
yvals = list(target_set)
ax.plot(xvals,yvals)
plt.tight_layout()
plt.savefig("test_pid_simulator.png")
fig.show()
print("\nDisplaying plot. Press enter to continue...(this closes the images)")
input()
plt.close()
