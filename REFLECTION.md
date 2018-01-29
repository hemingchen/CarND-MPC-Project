# Reflection

## Selection of N

* Determines how far car trajectory is optimized into the future by MPC.
* If N is too big, then the computing demand will increase. However, since uncertainties always exist in the future, it would be a waste of computing resource to optimize the trajectory way into the future.
* If N is too small, then the MPC would end up with some very local optimal results, which might affect the overall performance of the car when following the reference trajectory.

Selected value: N = 15

In my test, the MPC worked fairly well for N between [10, 20]. When N was less than 10, the car didn't drive as smooth and needed frequent corrections to steering and throttle. When N was too big, the MPC spent more time computing the trajectory and eventually caused delay in control.

## Selection of dt

* Determines how frequently the optimal trajectory is computed.
* Small dt demands more frequent updates from MPC which might cause delay.
* Large dt leads to less frequent updates from MPC but can result in large error in following reference trajectory.

Selected value: dt = 0.1

In my test, dt around 0.1 worked fairly well. When dt was increased too much, e.g. 1.0, the car wouldn't follow the reference trajectory well. I didn't try any value less than 0.1 as that wouldn't be computational efficient since the delay in control is already 0.1.


## Preprocessing before launching MPC

The preprocessing contains the following steps:

1. The reference waypoints (ptsx, ptsy) are translated into car coordinate frame - lines 109-117 in main.cpp.

2. The translated reference waypoints are fit to a 3rd order polynomial - lines 120-121 in main.cpp.

3. The delta from simulator (j[1]["steering_angle"]) is multiplied by -1 in all equations before sent toMPC, since it is referred in the opposite direction than the dynamic model in the MPC.

## Handling of 100ms control delay

The 100ms control delay means the initial states sent to the MPC is not the current states but those 100ms later. Therefore, the car states are projected by another 100ms and the states after the projection are used as the actual initial states by the MPC. Code can be found at lines 126-147 in main.cpp.


