import numpy as np
import math

def rotation(actions_taken, observation_0, observation_1, m, b):

    avoid_actions = np.argwhere(np.array(actions_taken) == 0).flatten()
    x_avoid = np.mean(observation_1[avoid_actions])
    y_avoid = np.mean(observation_0[avoid_actions])

    # find perpendicular line that passes through the mean of avoid points

    perpendicular_m = -1/m
    perpendicular_b = y_avoid - perpendicular_m * x_avoid

    # find intersection of lines
    x_line = ((b - perpendicular_b) / (perpendicular_m - m)).item()
    y_line = (m * x_line + b).item()

    avoid_vector = [x_avoid - x_line, y_avoid - y_line]
    #unit_vector = [1, 0]
    
    angle = math.degrees(math.acos(avoid_vector[0] / np.sqrt(avoid_vector[0]**2 + avoid_vector[1]**2)))

    return angle

def compute_reversal_time(total_rotation, idx_reversal):
    before = total_rotation[idx_reversal-45:idx_reversal]
    after = total_rotation[-45:]

    change = False
    t = 0

    for angle in total_rotation[idx_reversal:]:
        if not change:
            if angle > np.percentile(before, 95) or angle < np.percentile(before, 5):
                change = True
                t += 1
        else:
            if angle < np.percentile(after, 95) and angle > np.percentile(after, 5):
                break
            else:
                t += 1

    return t