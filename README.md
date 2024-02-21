# Energy Saving in Mobile Networks using Dynamic Programming

Energy consumption within the telecommunications sector is drastically increasing mainly due to the deployment of a large number of lower-power base stations, called small cells, aiming to meet the growing demands of network users. Therefore, developing enery saving strategies are essential to address this issue.

This project proposes an energy saving scheme that allows for the deactivation and activation of mobile base stations based on traffic variability throughout the day, aiming to improve the overall energy efficiency of the mobile network. 

The proposed algorithm applies Dynamic Programming along with Certainty Equivalent Control to find an optimal control policy by turning on and off the deployed base stations. Additionally, the UCB1 algorithm is applied to find the optimal configuration of the network interference management mechanism, minimizing interference between base stations and thus improving channel capacity. This approach has low computational overhead and does not require a large amount of  network information.

Simulation results demonstrate the effectiveness in achieving energy savings. These are accentuated as the number of deployed base stations increases. Moreover, it is observed to be scalable with the number of deployed base stations, where the computation time remains almost constant and is entirely feasible within the temporal scale at which the mechanism must operate.

This is the code repository of the project. You can find more information about it here: https://repositorio.upct.es/handle/10317/6132
