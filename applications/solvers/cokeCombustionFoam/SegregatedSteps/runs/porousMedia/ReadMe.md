

### Differences with the LB case
0. The fluid flow is different since the SRT is not accuracy
1. The flow is solved every iteration when using DBS
2. The method to compute the reactive surface area is different
3. The reactive surface area is re-calculated every interaction. However, the LB only update the reactive surface area when one cell is burned out 
4. The coke volume is different. DBS:0.99, LB:1
5. The solid cp varies with time
6. The mass equation is different. DBS count the produced gas into the fluid, leading to that the increasing velocity when rho is constant, independent of the temperature

### To Develop
1. [Done] Switch to consider if considering the gas production when solving NS equation 
2. [Done] report the residual coke
3. [Done] complete loop when the residual coke reached a specified target 

### To Test
1. High temperature when the combustion in convection-limited regime with high Pe

## To investigate 
1. [Done] the surface area
2. the diffusivity